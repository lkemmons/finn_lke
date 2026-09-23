"""AF input loading.

Replaces `code_anaconda/af_import.py` and the input-file logic of
`work_nrt.sec3_import_af`.  Reads MODIS/VIIRS active-fire products from
shapefile, zipped shapefile, CSV, or TXT and returns a single
GeoDataFrame ready to be turned into `work_pnt`.
"""
from __future__ import annotations

import datetime as dt
import logging
import re
import zipfile
from pathlib import Path
from typing import Iterable

import geopandas as gpd
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# FIRMS naming conventions we know how to recognise.
_RE_FIRE_ARCHIVE = re.compile(r"fire_archive_(.+)\.shp", re.IGNORECASE)
_RE_DL_FIRE = re.compile(r"DL_FIRE_(.+)\.zip", re.IGNORECASE)
_RE_NRT = re.compile(r"(MODIS_C6|VNP14IMGTDL_NRT)_(.+)\.shp", re.IGNORECASE)

# Columns we expect from FIRMS CSV/TXT inputs (and the dtypes ogr2ogr's
# .vrt template asked GDAL to enforce).  Missing columns are filled with
# defaults so the downstream code is uniform.
_AF_COLUMNS_DEFAULT = {
    "scan": "float64",
    "track": "float64",
    "acq_date": "object",
    "acq_time": "object",
    "satellite": "object",
    "instrument": "object",   # archive files supply this directly;
                              # NRT files don't, so step1.prep derives it
                              # from `satellite`.
    "confidence": "object",
    "version": "object",
    "brightness": "float64",
    "bright_t31": "float64",
    "bright_ti4": "float64",
    "bright_ti5": "float64",
    "frp": "float64",
    "daynight": "object",
    "type": "Int64",  # nullable int — archive only; mapped to anomtype
}


def _resolve_shp(path: Path) -> Path:
    """Mirror of `work_nrt.sec3_import_af`'s file-resolution logic.

    Returns a `.shp`/`.csv`/`.txt` path that actually exists, unzipping
    a sibling zip if needed.  Raises FileNotFoundError on failure.
    """
    if path.suffix.lower() in (".shp", ".csv", ".txt"):
        if path.exists():
            return path
        # Try to find a sibling .zip
        candidates: list[Path] = []
        m = _RE_FIRE_ARCHIVE.match(path.name)
        if m:
            candidates.append(path.parent / f"DL_FIRE_{m.group(1)}.zip")
        m = _RE_NRT.match(path.name)
        if m:
            candidates.append(path.with_suffix(".zip"))
        for z in candidates:
            if z.exists():
                log.info("Extracting %s -> %s", z, z.parent)
                with zipfile.ZipFile(z) as zf:
                    zf.extractall(z.parent)
                if path.exists():
                    return path
        raise FileNotFoundError(f"cannot find file: {path}")

    if path.suffix.lower() == ".zip":
        # Unzip into the same directory and look for a single .shp inside.
        with zipfile.ZipFile(path) as zf:
            shps = [n for n in zf.namelist() if n.lower().endswith(".shp")]
            if not shps:
                raise RuntimeError(f"no .shp inside {path}")
            zf.extractall(path.parent)
        return path.parent / shps[0]

    raise RuntimeError(f"unrecognized AF file extension: {path}")


def _read_one(path: Path) -> gpd.GeoDataFrame:
    """Read a single AF file as a WGS84 GeoDataFrame."""
    path = _resolve_shp(path)
    if path.suffix.lower() == ".shp":
        gdf = gpd.read_file(path)
    else:
        # CSV/TXT — replicate the .vrt that af_import.mk_vrt used to write.
        df = pd.read_csv(path)
        if "longitude" not in df.columns or "latitude" not in df.columns:
            raise RuntimeError(
                f"{path}: CSV must have 'longitude' and 'latitude' columns"
            )
        gdf = gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(df.longitude, df.latitude),
            crs="EPSG:4326",
        )

    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    elif gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs("EPSG:4326")

    # Lowercase columns for consistency.
    gdf.columns = [c.lower() if c != gdf.geometry.name else c for c in gdf.columns]

    # Backfill missing default columns.
    for col, dtype in _AF_COLUMNS_DEFAULT.items():
        if col not in gdf.columns:
            if dtype == "Int64":
                gdf[col] = pd.NA
            elif dtype == "float64":
                gdf[col] = np.nan
            else:
                gdf[col] = None
    return gdf


def load_af_files(
    paths: Iterable[Path],
    *,
    tropical_carryover_paths: Iterable[Path] | None = None,
    tropical_carryover_to_date: "dt.date | str | None" = None,
    tropical_lat_bounds: tuple[float, float] = (-30.0, 30.0),
) -> gpd.GeoDataFrame:
    """Load and concatenate multiple AF input files.

    Adds a `src_file` integer column (1-based) so the downstream code can
    keep track of which input each detection came from — equivalent to
    the `af_in_<i>` table layout.

    Tropical carryover
    ------------------
    If ``tropical_carryover_paths`` is supplied, those files are also
    loaded but only their detections within ``tropical_lat_bounds`` are
    kept, and their ``acq_date`` field is rewritten to
    ``tropical_carryover_to_date`` so the per-day pipeline groups them
    with the current day's detections.

    This replicates the original FINN's compensation for MODIS swath
    gaps in the tropics, **specifically for single-day workflows**.

    Note that ``step1.prep`` already duplicates each tropical MODIS
    detection into the *next* day within whatever data is loaded — but
    that mechanism is dormant in single-day runs because the previous
    day's file isn't in the input.  ``tropical_carryover_paths`` plugs
    that gap by loading yesterday's MODIS file explicitly.  The two
    mechanisms don't double-count: when ``date`` selects today,
    ``step1.prep``'s duplicates land on tomorrow and get dropped by the
    date filter.

    VIIRS doesn't need this because its wider swath leaves no
    comparable gap.

    Carryover detections are marked with `src_file=0` (rather than 1+)
    so they're identifiable later if you want to inspect / exclude them.
    """
    frames: list[gpd.GeoDataFrame] = []

    for i, p in enumerate(paths, start=1):
        gdf = _read_one(Path(p))
        gdf["src_file"] = i
        log.info("loaded %s: %d rows", p, len(gdf))
        frames.append(gdf)
    if not frames:
        raise ValueError("no AF files supplied")

    # ---- tropical carryover --------------------------------------------------
    if tropical_carryover_paths:
        if tropical_carryover_to_date is None:
            raise ValueError(
                "tropical_carryover_to_date is required when "
                "tropical_carryover_paths is non-empty"
            )
        # Normalise the target date to a 'YYYY-MM-DD' string that matches
        # what the FIRMS CSVs already store.
        if hasattr(tropical_carryover_to_date, "isoformat"):
            target_str = tropical_carryover_to_date.isoformat()
        else:
            target_str = str(tropical_carryover_to_date)

        lat_min, lat_max = tropical_lat_bounds
        carry_frames: list[gpd.GeoDataFrame] = []
        for p in tropical_carryover_paths:
            gdf = _read_one(Path(p))
            n_in = len(gdf)
            in_tropics = (
                (gdf.geometry.y >= lat_min) & (gdf.geometry.y <= lat_max)
            )
            gdf = gdf.loc[in_tropics].copy()
            gdf["acq_date"] = target_str
            gdf["src_file"] = 0
            log.info(
                "carryover %s: %d -> %d tropical rows, re-dated to %s",
                p, n_in, len(gdf), target_str,
            )
            if len(gdf):
                carry_frames.append(gdf)
        frames.extend(carry_frames)

    combined = pd.concat(frames, ignore_index=True)
    return gpd.GeoDataFrame(combined, geometry="geometry", crs="EPSG:4326")
