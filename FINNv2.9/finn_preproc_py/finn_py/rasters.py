"""Raster I/O and MODIS HDF handling — PyHDF version.

This module exists because many GDAL builds (conda-forge wheels,
some Homebrew/apt packages) ship without the HDF4 driver, so
``rasterio.open("MCD12Q1.A2020001.h08v05.061.*.hdf")`` raises
``RasterioIOError``.  We work around that by:

  1. Reading the raw arrays out of MODIS HDF4-EOS files with PyHDF.
  2. Parsing the ``StructMetadata.0`` global attribute to recover each
     tile's sinusoidal-meters transform.
  3. Writing each tile to a temporary single-tile GeoTIFF in
     sinusoidal CRS (rasterio writes GTiff without needing HDF4 support).
  4. Using rasterio.merge + warp to build the WGS84 mosaic / COG.

After step (3) we are out of HDF4 land and the rest of the pipeline
runs on plain GeoTIFFs.

Replaces the parts of ``code_anaconda/rst_import.py`` that built the
``raster.rst_<tag>`` PostGIS tables and the ``raster.wireframe`` tile
geometry table.
"""
from __future__ import annotations

import logging
import re
import tempfile
from pathlib import Path
from typing import Iterable, Sequence

import contextlib
import logging as _logging

import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.warp import (Resampling, calculate_default_transform, reproject)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Polar-warning suppression for sinu -> WGS84 warps
# ---------------------------------------------------------------------------
#
# MODIS's sinusoidal projection is only valid up to y = ±R·π/2 ≈ ±10,007,543 m
# (the geographic poles).  The MODIS tile grid is defined geometrically
# regardless of projection validity, so the polar tile rows (v=00, 01, 16, 17)
# have rows of pixels whose sinusoidal coordinates correspond to latitudes
# beyond ±90°.  When the warper inverse-transforms those pixels to WGS84,
# PROJ correctly reports "sinu: Invalid latitude" — once per invalid pixel,
# millions of times for a global mosaic.
#
# These are *warnings* (GDAL err_no=1 = CE_Warning), not errors.  GDAL writes
# nodata for the invalid pixels and continues.  The output is correct.
#
# rasterio routes GDAL warnings through a `rasterio._env` logger at INFO
# level, which makes them visible by default.  We install a temporary filter
# that drops just the "Invalid latitude" line during the warp and leaves
# every other rasterio/GDAL message alone.

class _ProjPolarFilter(_logging.Filter):
    """Suppress 'PROJ: sinu: Invalid latitude' warnings, nothing else."""

    def filter(self, record: _logging.LogRecord) -> bool:
        msg = record.getMessage()
        return "Invalid latitude" not in msg


@contextlib.contextmanager
def _quiet_proj_polar_warnings():
    """Context manager that filters out PROJ polar-latitude warnings
    while the warp runs.  Any other GDAL warnings/errors pass through.
    """
    target = _logging.getLogger("rasterio._env")
    f = _ProjPolarFilter()
    target.addFilter(f)
    try:
        yield
    finally:
        target.removeFilter(f)


def _import_pyhdf():
    """Lazy import — only needed for MODIS HDF reads, not for using
    pre-built mosaics, computing tile geometry, or running step2.
    """
    try:
        from pyhdf.SD import SD, SDC
    except ImportError as e:
        raise ImportError(
            "pyhdf is required to read MODIS HDF4 files but is not "
            "installed.  Install it with: "
            "`conda install -c conda-forge hdf4 pyhdf` (recommended) "
            "or `apt install libhdf4-dev libjpeg-dev && pip install pyhdf` "
            "on Debian/Ubuntu."
        ) from e
    return SD, SDC

# MODIS sinusoidal CRS as proj4.  Equivalent to EPSG-ish ID 9006842 used
# by the original FINN code — built explicitly so we don't depend on
# whether the local PROJ database has the entry.
MODIS_SINU_PROJ4 = (
    "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 "
    "+a=6371007.181 +b=6371007.181 +units=m +no_defs"
)
WGS84 = "EPSG:4326"

# MODIS HDF4-EOS subdataset names we use.
SDS_LCT_MCD12Q1 = "LC_Type1"
SDS_VCF_MOD44B = (
    "Percent_Tree_Cover",
    "Percent_NonTree_Vegetation",
    "Percent_NonVegetated",
)


# ---------------------------------------------------------------------------
# Low-level HDF4 reading via PyHDF
# ---------------------------------------------------------------------------

# ODL-ish key=value extractor for StructMetadata.0.
_RE_KV = re.compile(r"(\w+)\s*=\s*([^\n]+)")
_RE_TUPLE = re.compile(r"\(\s*([^)]+)\s*\)")


def _parse_struct_metadata(md: str) -> dict:
    """Pull the fields we need out of StructMetadata.0.

    StructMetadata.0 is a multi-page PVL/ODL-flavoured text blob attached
    as a global attribute on every HDF-EOS file.  We don't need to parse
    the whole thing — for our purposes the grid dims and the corner
    coordinates are enough to build a rasterio transform.
    """
    def grab(key: str, required: bool = True) -> str | None:
        m = re.search(rf"\b{key}\s*=\s*([^\n]+)", md)
        if m is None:
            if required:
                raise RuntimeError(f"StructMetadata.0 missing field {key!r}")
            return None
        return m.group(1).strip().rstrip(";").strip()

    def grab_tuple(key: str) -> tuple[float, float]:
        raw = grab(key)
        inner = _RE_TUPLE.search(raw)
        if inner is None:
            raise RuntimeError(f"StructMetadata.0: {key} is not a tuple: {raw!r}")
        parts = [p.strip() for p in inner.group(1).split(",")]
        return float(parts[0]), float(parts[1])

    xdim = int(grab("XDim"))
    ydim = int(grab("YDim"))
    ul_x, ul_y = grab_tuple("UpperLeftPointMtrs")
    lr_x, lr_y = grab_tuple("LowerRightMtrs")
    grid_name = grab("GridName").strip('"')
    projection = grab("Projection", required=False)

    if projection and "SNSOID" not in projection.upper():
        log.warning("expected sinusoidal grid, got Projection=%s", projection)

    return {
        "grid_name": grid_name,
        "xdim": xdim, "ydim": ydim,
        "ulx": ul_x, "uly": ul_y,
        "lrx": lr_x, "lry": lr_y,
    }


def read_modis_hdf(
    hdf_path: Path,
    sds_name: str,
) -> tuple[np.ndarray, "rasterio.transform.Affine", float | int | None]:
    """Read one named subdataset from a MODIS HDF4-EOS file.

    Returns (array, sinusoidal-meters transform, nodata).  The array is
    a NumPy view shaped (rows, cols); the transform places the upper-left
    corner of pixel (0, 0) at (ulx, uly) in MODIS sinusoidal meters and
    assumes square pixels of size ``(lrx - ulx) / xdim``.
    """
    SD, SDC = _import_pyhdf()
    hdf = SD(str(hdf_path), SDC.READ)
    try:
        try:
            struct_md_str = hdf.attributes()["StructMetadata.0"]
        except KeyError as e:
            raise RuntimeError(
                f"{hdf_path}: no StructMetadata.0 attribute — not an HDF-EOS file?"
            ) from e
        meta = _parse_struct_metadata(struct_md_str)

        # Subdataset
        try:
            sds = hdf.select(sds_name)
        except Exception as e:
            avail = list(hdf.datasets().keys())
            raise RuntimeError(
                f"{hdf_path}: subdataset {sds_name!r} not found; available: {avail}"
            ) from e
        try:
            data = sds.get()  # numpy array
            attrs = sds.attributes()
        finally:
            sds.endaccess()
    finally:
        hdf.end()

    # MODIS uses several spellings for fill / nodata.
    nodata = None
    for k in ("_FillValue", "FillValue", "fillvalue"):
        if k in attrs:
            nodata = attrs[k]
            break

    # Sanity check: array shape should match StructMetadata dims.  If it
    # does not, trust the array and stretch the transform proportionally.
    h, w = data.shape[-2:]
    if (h, w) != (meta["ydim"], meta["xdim"]):
        log.debug(
            "%s/%s: dims %s do not match StructMetadata (%s, %s) — using array dims",
            hdf_path.name, sds_name, (h, w), meta["ydim"], meta["xdim"],
        )

    transform = from_bounds(
        west=meta["ulx"], south=meta["lry"],
        east=meta["lrx"], north=meta["uly"],
        width=w, height=h,
    )
    return data, transform, nodata


# ---------------------------------------------------------------------------
# Per-tile HDF -> sinusoidal GeoTIFF (the bridge out of HDF4 land)
# ---------------------------------------------------------------------------

def hdf_to_sinu_tif(
    hdf_path: Path,
    sds_names: Sequence[str],
    out_tif: Path,
) -> Path:
    """Write a multi-band sinusoidal GeoTIFF for one MODIS HDF tile.

    All requested subdatasets must share dims and transform (they do for
    MOD44B's three VCF bands).  After this call the file is a plain
    GeoTIFF and the rest of the pipeline does not need HDF4 support.
    """
    arrays: list[np.ndarray] = []
    transform = None
    nodata = None
    dtype = None

    for sds in sds_names:
        arr, tf, nd = read_modis_hdf(hdf_path, sds)
        arrays.append(arr)
        if transform is None:
            transform = tf
            nodata = nd
            dtype = arr.dtype
        else:
            if arr.shape != arrays[0].shape:
                raise RuntimeError(
                    f"subdataset {sds!r} in {hdf_path} has shape {arr.shape}, "
                    f"first SDS had {arrays[0].shape}"
                )

    out_tif.parent.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff",
        "height": arrays[0].shape[0],
        "width": arrays[0].shape[1],
        "count": len(arrays),
        "dtype": str(dtype),
        "crs": MODIS_SINU_PROJ4,
        "transform": transform,
        "nodata": nodata,
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
        "compress": "deflate",
    }
    with rasterio.open(out_tif, "w", **profile) as dst:
        for i, (arr, sds) in enumerate(zip(arrays, sds_names), start=1):
            dst.write(arr, i)
            dst.set_band_description(i, sds)
    return out_tif


# ---------------------------------------------------------------------------
# Public mosaic API
# ---------------------------------------------------------------------------

def build_mosaic(
    hdf_paths: Iterable[Path],
    sds_names: str | Sequence[str],
    out_tif: Path,
    resolution_deg: float = 0.005,        # ~ 500 m near equator
    resampling: Resampling = Resampling.nearest,
    keep_intermediate: bool = False,
    num_threads: int = 0,                 # 0 = ALL_CPUS
    warp_mem_mb: int = 2048,
) -> Path:
    """Build a WGS84 GeoTIFF mosaic from MODIS HDF4 tiles.

    Memory-bounded: the destination is opened as a tiled GeoTIFF and
    written by GDAL's warper block by block, so global mosaics work the
    same as regional ones and peak RAM is bounded by ``warp_mem_mb``
    (default 2 GB).  Previously this function held the full sinusoidal
    mosaic *and* the full destination as numpy arrays, which OOM'd
    notebook kernels even for modest AOIs.

    Pipeline:
        HDF4 tile  ──PyHDF──>  numpy  ──rasterio.write──>  sinusoidal .tif
            (one per tile, in a temp dir)
        sinusoidal .tifs                ──GDAL VRT──>     virtual sinu mosaic
        virtual sinu mosaic             ──GDAL warp──>    tiled WGS84 GeoTIFF
            (streaming, block by block — no full-mosaic numpy array)

    For global mosaics, polar tile rows (v=00, 01, 16, 17) extend
    geometrically beyond y=±R·π/2 in the sinusoidal grid, where PROJ's
    inverse transform isn't defined (lat > 90°).  The warper writes
    nodata for those pixels and emits a "PROJ: sinu: Invalid latitude"
    warning at GDAL's CE_Warning level.  These are not errors — they
    fire once per invalid pixel, can total millions for a global build,
    and are noise — so we filter them out for the duration of the
    warp.  Any other GDAL warnings / errors still surface normally.
    """
    if isinstance(sds_names, str):
        sds_names = (sds_names,)
    sds_names = list(sds_names)
    out_tif = Path(out_tif)
    out_tif.parent.mkdir(parents=True, exist_ok=True)

    hdf_paths = [Path(p) for p in hdf_paths]
    if not hdf_paths:
        raise ValueError("no HDF inputs")

    workdir = Path(tempfile.mkdtemp(prefix="finn_mosaic_"))
    log.info("mosaicing %d HDF tiles -> %s (workdir=%s)",
             len(hdf_paths), out_tif, workdir)

    try:
        # Stage 1: HDF -> sinusoidal GeoTIFF per tile (already streaming —
        # each tile is small, ~5 MB for LCT, ~65 MB for VCF).
        sinu_tifs: list[Path] = []
        for p in hdf_paths:
            tif = workdir / f"{p.stem}.sinu.tif"
            hdf_to_sinu_tif(p, sds_names, tif)
            sinu_tifs.append(tif)

        # Stage 2: virtual mosaic over the sinu tiles.  The VRT is a
        # tiny XML file referencing the tiles; no pixel data is moved.
        vrt_path = workdir / "all.sinu.vrt"
        _build_vrt(sinu_tifs, vrt_path)

        # Stage 3: stream-reproject the VRT to a tiled WGS84 COG.
        # Both source and destination are open datasets, so GDAL's
        # warper iterates over destination blocks, reads only the
        # needed source windows, and writes block-by-block.
        with rasterio.open(vrt_path) as src:
            dst_transform, dst_w, dst_h = calculate_default_transform(
                src.crs, WGS84,
                src.width, src.height, *src.bounds,
                resolution=resolution_deg,
            )
            dtype = src.dtypes[0]
            src_nodata = src.nodata

            profile = {
                "driver": "GTiff",
                "height": dst_h,
                "width": dst_w,
                "count": src.count,
                "dtype": dtype,
                "crs": WGS84,
                "transform": dst_transform,
                "nodata": src_nodata,
                "tiled": True,
                "blockxsize": 512,
                "blockysize": 512,
                "compress": "deflate",
                "predictor": 2 if "int" in dtype else 1,
                "BIGTIFF": "IF_SAFER",
            }
            log.info("destination: %d x %d x %d band(s) %s, ~%.1f GB uncompressed",
                     dst_w, dst_h, src.count, dtype,
                     dst_w * dst_h * src.count * np.dtype(dtype).itemsize / 1e9)

            with rasterio.open(out_tif, "w", **profile) as dst:
                # Reproject each band separately — passing open Band
                # objects keeps the warp streaming.
                with _quiet_proj_polar_warnings():
                    for i in range(1, src.count + 1):
                        reproject(
                            source=rasterio.band(src, i),
                            destination=rasterio.band(dst, i),
                            src_transform=src.transform,
                            src_crs=src.crs,
                            dst_transform=dst_transform,
                            dst_crs=WGS84,
                            resampling=resampling,
                            src_nodata=src_nodata,
                            dst_nodata=src_nodata,
                            num_threads=num_threads or 0,
                            warp_mem_limit=warp_mem_mb,
                        )
                        # Preserve band description from the sinu mosaic.
                        desc = src.descriptions[i - 1] if src.descriptions else None
                        dst.set_band_description(i, desc or sds_names[i - 1])

    finally:
        if not keep_intermediate:
            import shutil
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            log.info("keeping intermediate dir %s", workdir)

    log.info("wrote %s", out_tif)
    return out_tif


def _build_vrt(sources: Sequence[Path], vrt_path: Path) -> Path:
    """Build a GDAL VRT over the given GeoTIFFs.

    Tries the in-process osgeo.gdal binding first (usually installed
    alongside rasterio in conda envs); falls back to invoking the
    gdalbuildvrt CLI tool if not.  Either way the result is a tiny
    .vrt XML file referencing the input tiles by path — rasterio.open
    then treats it as a single logical raster.
    """
    sources = [str(Path(p)) for p in sources]
    vrt_path = Path(vrt_path)

    # Try in-process osgeo.gdal first.
    try:
        from osgeo import gdal  # type: ignore
        gdal.UseExceptions()
        ds = gdal.BuildVRT(str(vrt_path), sources)
        if ds is None:
            raise RuntimeError("gdal.BuildVRT returned None")
        ds = None
        if vrt_path.exists():
            return vrt_path
    except ImportError:
        pass  # fall through to subprocess

    # Fall back to subprocess.
    import shutil as _shutil, subprocess
    if _shutil.which("gdalbuildvrt") is None:
        raise RuntimeError(
            "Cannot build a mosaic VRT: neither the osgeo.gdal Python "
            "binding nor the gdalbuildvrt CLI tool is available.  Install "
            "one of them, e.g. with `conda install -c conda-forge gdal`.  "
            "(On HPC systems this is usually already in the conda env "
            "alongside rasterio.)"
        )
    cmd = ["gdalbuildvrt", str(vrt_path), *sources]
    log.debug("falling back to subprocess: gdalbuildvrt %s ... (%d files)",
              vrt_path, len(sources))
    res = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            f"gdalbuildvrt failed (exit {res.returncode}):\n"
            f"  stderr: {res.stderr.strip()}"
        )
    if not vrt_path.exists():
        raise RuntimeError(f"gdalbuildvrt did not produce {vrt_path}")
    return vrt_path


# ---------------------------------------------------------------------------
# MODIS tile geometry (replaces raster.wireframe + modis_tile.py)
# ---------------------------------------------------------------------------

def modis_tile_polygons():
    """GeoDataFrame of MODIS sinusoidal tile boundaries in WGS84.

    No HDF reading involved — this is pure geometry.  Replaces the
    `raster.wireframe` table the SQL pipeline relied on.
    """
    import geopandas as gpd
    from pyproj import Transformer
    from shapely.geometry import Polygon
    from shapely.ops import transform as sh_transform

    TILE = 1111950.5196666666  # MODIS tile size in sinusoidal meters

    tx = Transformer.from_crs(MODIS_SINU_PROJ4, WGS84, always_xy=True)
    tiles = []
    for h in range(36):
        for v in range(18):
            x0 = (h - 18) * TILE
            y1 = (9 - v) * TILE
            x1 = x0 + TILE
            y0 = y1 - TILE
            sinu = Polygon([(x0, y0), (x0, y1), (x1, y1), (x1, y0), (x0, y0)])
            # densify so polar curvature is well represented in lon/lat
            sinu = sinu.segmentize(TILE / 20)
            wgs = sh_transform(tx.transform, sinu)
            if wgs.is_empty or not wgs.is_valid:
                continue
            tiles.append({
                "tilename": f"h{h:02d}v{v:02d}",
                "ih": h, "iv": v,
                "geometry": wgs,
            })
    return gpd.GeoDataFrame(tiles, geometry="geometry", crs=WGS84)


def tiles_needed(af_gdf, tile_gdf=None) -> list[str]:
    """Return the MODIS tile names that cover the given AF detections."""
    import geopandas as gpd
    if tile_gdf is None:
        tile_gdf = modis_tile_polygons()
    j = gpd.sjoin(af_gdf, tile_gdf, how="inner", predicate="within")
    return sorted(j["tilename"].unique().tolist())


# ---------------------------------------------------------------------------
# Convenience: extract h##v## from a MODIS filename
# ---------------------------------------------------------------------------

_RE_TILENAME = re.compile(r"\.h(\d{2})v(\d{2})\.", re.IGNORECASE)


def tilename_of(hdf_path: Path) -> str | None:
    """Return 'h08v05' (etc.) from a MODIS filename, or None if not found."""
    m = _RE_TILENAME.search(Path(hdf_path).name)
    return f"h{m.group(1)}v{m.group(2)}" if m else None
