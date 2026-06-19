"""Raster preparation pipeline (local files only — no downloads).

Replaces ``preprocessor/code_bashinterface/work_raster.py``, but only
its raster-mosaicking responsibilities.  All inputs are expected to be
already on disk; the original's Earthdata download step is gone.

Inputs:
  - One or more **MCD12Q1** HDF tiles (MODIS land cover).
  - One or more **MOD44B** HDF tiles (MODIS VCF — three percentage
    bands: tree, herb, bare).

Outputs (placed in ``raster_dir``):
  - ``modlct_<year>.tif``  — single-band WGS84 mosaic of LCT.
  - ``modvcf_<year>.tif``  — three-band WGS84 mosaic of VCF.

The region polygon file (``regnum.gpkg``) is supplied separately and
must already exist in ``raster_dir`` before ``finn-nrt`` is run.  The
column carrying the region number defaults to ``Region_num``; override
via the ``variable_in`` field of the regnum RasterSpec in
:mod:`finn_py.config`.

These are the files ``finn-nrt`` expects in ``--raster_dir``:
``modlct_<year>.tif``, ``modvcf_<year>.tif``, ``regnum.gpkg``.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

from . import rasters
from .rasters import SDS_LCT_MCD12Q1, SDS_VCF_MOD44B

log = logging.getLogger(__name__)


def run_raster(
    year_rst: int,
    lct_files: Sequence[Path],
    vcf_files: Sequence[Path],
    raster_dir: Path,
    *,
    lct_resolution_deg: float = 0.005,      # ~ 500 m at equator (MCD12Q1)
    vcf_resolution_deg: float = 0.0025,     # ~ 250 m at equator (MOD44B)
) -> dict[str, Path]:
    """Mosaic local MODIS HDFs into the two GeoTIFFs ``finn-nrt`` consumes.

    Parameters
    ----------
    year_rst
        Used only as a label in the output filenames (``modlct_<year>.tif``
        etc.).  Must match the year of data the input HDFs actually cover.
    lct_files
        List of MCD12Q1 HDF paths.  All tiles will be mosaicked into one
        WGS84 GeoTIFF.
    vcf_files
        List of MOD44B HDF paths.  All three VCF bands are stacked into
        one 3-band GeoTIFF.
    raster_dir
        Output directory; created if missing.  The region polygon file
        (``regnum.gpkg``) is **not** produced here — drop it into this
        directory yourself before running the NRT pipeline.

    Returns
    -------
    dict
        ``{'lct': <path>, 'vcf': <path>}``.
    """
    raster_dir = Path(raster_dir)
    raster_dir.mkdir(parents=True, exist_ok=True)

    lct_files = [_check_exists(p, "lct") for p in lct_files]
    vcf_files = [_check_exists(p, "vcf") for p in vcf_files]

    if not lct_files:
        raise ValueError("no LCT (MCD12Q1) HDF files supplied")
    if not vcf_files:
        raise ValueError("no VCF (MOD44B) HDF files supplied")

    out: dict[str, Path] = {}

    # --- LCT (single band, thematic) ---
    lct_tif = raster_dir / f"modlct_{year_rst}.tif"
    log.info("building LCT mosaic from %d HDF tiles -> %s",
             len(lct_files), lct_tif)
    rasters.build_mosaic(
        lct_files, SDS_LCT_MCD12Q1, lct_tif,
        resolution_deg=lct_resolution_deg,
    )
    out["lct"] = lct_tif

    # --- VCF (three bands, continuous) ---
    vcf_tif = raster_dir / f"modvcf_{year_rst}.tif"
    log.info("building VCF mosaic from %d HDF tiles -> %s",
             len(vcf_files), vcf_tif)
    rasters.build_mosaic(
        vcf_files, SDS_VCF_MOD44B, vcf_tif,
        resolution_deg=vcf_resolution_deg,
    )
    out["vcf"] = vcf_tif

    return out


def _check_exists(p: Path, kind: str) -> Path:
    p = Path(p)
    if not p.exists():
        raise FileNotFoundError(f"{kind} input not found: {p}")
    return p
