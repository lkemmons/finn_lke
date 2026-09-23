"""CLI entry points.

Two commands matching the original FINN scripts:

  finn-nrt        — equivalent of work_nrt.py
  finn-raster     — equivalent of work_raster.py

Same flag names and defaults as the originals, so an existing wrapper
script (``finn-nrt -t mytag -y 2020 -o ./out myfires.shp``) continues to
work.
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
from pathlib import Path

from .config import FinnConfig
from .pipeline import run_nrt
from .raster_pipeline import run_raster


def _str2bool(v: str | bool) -> bool:
    if isinstance(v, bool):
        return v
    s = v.lower()
    if s in ("yes", "true", "t", "y", "1"):
        return True
    if s in ("no", "false", "f", "n", "0"):
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got {v!r}")


def _yyyyjjj(s: str | int | None) -> dt.date | None:
    if s is None:
        return None
    return dt.datetime.strptime(str(s), "%Y%j").date()


def _setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )


# ---------------------------------------------------------------------------
# finn-nrt   ( ~= work_nrt.py )
# ---------------------------------------------------------------------------

def nrt_main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="finn-nrt",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run the FINN preprocessor (pure-Python) on AF input(s).",
    )
    req = p.add_argument_group("required arguments")
    req.add_argument("-t", "--tag_af", required=True, type=str,
                     help="tag for AF processing")
    req.add_argument("-y", "--year_rst", required=True, type=int,
                     help="dataset year for raster")
    req.add_argument("-o", "--out_dir", required=True, type=Path,
                     help="output directory")

    req.add_argument("-d", "--date", required=True, type=int,
                     help="processing date as YYYYJJJ (e.g. 2026100)")
    p.add_argument("-s", "--summary_file", default=None, type=Path,
                   help="summary filename (appended to)")

    p.add_argument("--work_dir", default=Path("./work"), type=Path,
                   help="working directory for intermediate parquet files")
    p.add_argument("--raster_dir", default=Path("./rasters"), type=Path,
                   help="directory containing modlct_<year>.tif etc.")

    p.add_argument("--run_import", default=True, type=_str2bool,
                   nargs="?", const=True,
                   help="re-import AF and rebuild work_pnt [yes/no]")
    p.add_argument("--run_step1", default=True, type=_str2bool,
                   nargs="?", const=True,
                   help="run step1 (geometry processing) [yes/no]")
    p.add_argument("--run_step2", default=True, type=_str2bool,
                   nargs="?", const=True,
                   help="run step2 (lct/vcf identification) [yes/no]")
    p.add_argument("--date_definition", default="UTC", choices=("UTC", "LST"),
                   help="date definition for daily grouping")
    p.add_argument("--no_filter_persistent", dest="filter_persistent_sources",
                   action="store_false",
                   help="do NOT drop volcano/other persistent AF detections")
    p.add_argument("--workers", default=1, type=int,
                   help="parallel worker processes (1 = sequential)")
    p.add_argument("-v", "--verbose", action="store_true")

    p.add_argument("af_fnames", nargs="+", type=Path,
                   help="AF file name(s)")

    args = p.parse_args(argv)
    _setup_logging(args.verbose)

    cfg = FinnConfig(
        tag_af=args.tag_af,
        year_rst=args.year_rst,
        af_files=args.af_fnames,
        out_dir=args.out_dir,
        date=_yyyyjjj(args.date),
        work_dir=args.work_dir,
        raster_dir=args.raster_dir,
        filter_persistent_sources=args.filter_persistent_sources,
        date_definition=args.date_definition,
        run_import=args.run_import,
        run_step1=args.run_step1,
        run_step2=args.run_step2,
        workers=args.workers,
        summary_file=args.summary_file,
    )
    out = run_nrt(cfg)
    print(out, file=sys.stdout)
    return 0


# ---------------------------------------------------------------------------
# finn-raster   ( ~= work_raster.py )
# ---------------------------------------------------------------------------

def raster_main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="finn-raster",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Mosaic local MODIS HDF tiles into the GeoTIFF inputs "
                    "that finn-nrt expects in --raster_dir.  Does no "
                    "downloads — all inputs must already be on disk.",
    )
    req = p.add_argument_group("required arguments")
    req.add_argument("-y", "--year_rst", required=True, type=int,
                     help="year label for output filenames (modlct_<year>.tif)")
    req.add_argument("--lct", nargs="+", required=True, type=Path,
                     metavar="HDF",
                     help="MCD12Q1 land-cover HDF file(s); shell globs OK")
    req.add_argument("--vcf", nargs="+", required=True, type=Path,
                     metavar="HDF",
                     help="MOD44B VCF HDF file(s); shell globs OK")
    p.add_argument("--raster_dir", default=Path("./rasters"), type=Path,
                   help="output directory for the mosaics. The region "
                        "polygon file (regnum.gpkg) is NOT built here — "
                        "drop it into this directory yourself")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)
    _setup_logging(args.verbose)

    out = run_raster(
        year_rst=args.year_rst,
        lct_files=args.lct,
        vcf_files=args.vcf,
        raster_dir=args.raster_dir,
    )
    for k, v in out.items():
        print(f"{k}: {v}", file=sys.stdout)
    return 0
