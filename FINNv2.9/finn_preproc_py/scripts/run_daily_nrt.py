#!/usr/bin/env python3
"""Run the finn_py NRT pipeline for a single day.

Examples
--------
    # One day, locally
    ./run_daily_nrt.py --date 2026100 \\
        --af-root    /glade/derecho/scratch/emmons/finn3_inputs \\
        --raster-dir /glade/derecho/scratch/emmons/finn3_rasters \\
        --work-dir   /glade/derecho/scratch/emmons/finn3_work \\
        --out-dir    /glade/derecho/scratch/emmons/finn3_out

    # Fanned out via PBS array — see run_daily_nrt.pbs.template
    qsub run_daily_nrt.pbs.template

Output layout
-------------
For the requested date YYYYJJJ this produces:

    <work-dir>/af_<tag>_<YYYYJJJ>/   per-day working parquets + log
    <out-dir>/out_<tag>_<YYYYJJJ>_modlct_<year>_modvcf_<year>_regnum.csv
    <out-dir>/out_<tag>_<YYYYJJJ>_modlct_<year>_modvcf_<year>_regnum.gpkg
    <out-dir>/summary_<tag>_<YYYYJJJ>.txt

Exit code is 0 on success (or when output already exists and
``--skip-existing`` is on), 1 on failure.
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
import time
from pathlib import Path

# Allow `python scripts/run_daily_nrt.py ...` from a checkout without install.
_SCRIPT = Path(__file__).resolve()
_FINN_ROOT = _SCRIPT.parent.parent
if str(_FINN_ROOT) not in sys.path:
    sys.path.insert(0, str(_FINN_ROOT))

from finn_py.config import FinnConfig                              # noqa: E402
from finn_py.pipeline import run_nrt                               # noqa: E402

log = logging.getLogger("run_daily_nrt")


# ---------------------------------------------------------------------------
# Date conversion (YYYYJJJ <-> datetime.date)
# ---------------------------------------------------------------------------

def yyyyjjj_to_date(yyyyjjj: int) -> dt.date:
    s = f"{int(yyyyjjj):07d}"
    return dt.date(int(s[:4]), 1, 1) + dt.timedelta(days=int(s[4:]) - 1)


def date_to_yyyyjjj(d: dt.date) -> int:
    return d.year * 1000 + d.timetuple().tm_yday


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="run_daily_nrt.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    req = p.add_argument_group("required")
    req.add_argument("--date", required=True, type=int,
                     help="processing date as YYYYJJJ (e.g. 2026100)")
    req.add_argument("--af-root", required=True, type=Path,
                     help="directory containing FIRMS NRT input files")
    req.add_argument("--raster-dir", required=True, type=Path,
                     help="directory with pre-built modlct_<year>.tif, "
                          "modvcf_<year>.tif, regnum.gpkg")
    req.add_argument("--work-dir", required=True, type=Path,
                     help="root working directory (per-day subdir created here)")
    req.add_argument("--out-dir",  required=True, type=Path,
                     help="output directory")

    p.add_argument("--tag-prefix", default="nrt",
                   help="output tag = '<prefix>_<YYYYJJJ>'  (default: %(default)s)")
    p.add_argument("--year-rst", type=int, default=None,
                   help="year for raster file selection "
                        "(default: year derived from the date)")
    p.add_argument("--af-modis",
                   default="MODIS_C6_1_Global_MCD14DL_NRT_{yyyyjjj}.txt",
                   help="filename pattern under --af-root for MODIS AF inputs; "
                        "{yyyyjjj} is substituted.  Empty string = skip MODIS.")
    p.add_argument("--af-viirs",
                   default="SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_{yyyyjjj}.txt",
                   help="filename pattern for VIIRS AF inputs; "
                        "empty string = skip VIIRS.")
    p.add_argument("--date-definition", default="UTC", choices=["UTC", "LST"],
                   help="how to assign each detection to a calendar day "
                        "(default: %(default)s)")
    p.add_argument("--tropical-carryover", action="store_true",
                   help="include the previous day's MODIS detections "
                        "within --tropics-lat-bounds, re-labelled with "
                        "the current day's date.  Replicates the original "
                        "FINN's compensation for MODIS swath gaps; "
                        "without it, finn_py produces noticeably fewer "
                        "tropical polygons than FINN2.")
    p.add_argument("--tropics-lat-bounds", nargs=2, type=float,
                   default=[-23.5, 23.5],
                   metavar=("LAT_MIN", "LAT_MAX"),
                   help="latitude bounds for tropical carryover "
                        "(default: -23.5 23.5, matching the original "
                        "FINN's step1_prep_v7m.sql duplication threshold).")
    p.add_argument("--skip-existing", action="store_true", default=True,
                   help="skip if output file already exists (default)")
    p.add_argument("--no-skip-existing", action="store_false", dest="skip_existing")
    p.add_argument("-v", "--verbose", action="store_true")

    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# AF-file lookups
# ---------------------------------------------------------------------------

def find_af_files(yyyyjjj: int, af_root: Path,
                   modis_pat: str, viirs_pat: str) -> list[Path]:
    """Look up the AF input files that exist for a given date."""
    found = []
    for pat in (modis_pat, viirs_pat):
        if not pat:
            continue
        p = af_root / pat.format(yyyyjjj=f"{yyyyjjj:07d}")
        if p.exists():
            found.append(p)
        else:
            log.warning("AF file missing for %d: %s", yyyyjjj, p)
    return found


def find_carryover_modis(yyyyjjj: int, af_root: Path,
                          modis_pat: str) -> list[Path]:
    """Look up the *previous day's* MODIS file for tropical carryover."""
    if not modis_pat:
        return []
    prev = date_to_yyyyjjj(yyyyjjj_to_date(yyyyjjj) - dt.timedelta(days=1))
    p = af_root / modis_pat.format(yyyyjjj=f"{prev:07d}")
    if p.exists():
        log.info("tropical carryover from previous day: %s", p)
        return [p]
    log.warning("tropical carryover requested but previous-day MODIS "
                "file not found: %s", p)
    return []


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    yyyyjjj = args.date
    date = yyyyjjj_to_date(yyyyjjj)
    tag = f"{args.tag_prefix}_{yyyyjjj}"
    year_rst = args.year_rst or date.year

    log.info("=" * 70)
    log.info("DATE %d  (%s)", yyyyjjj, date.isoformat())
    log.info("=" * 70)

    # Skip if any output for this tag already exists.
    existing = list(args.out_dir.glob(f"out_{tag}_*.gpkg"))
    if args.skip_existing and existing:
        log.info("output exists; skipping: %s", existing[0].name)
        return 0

    af_files = find_af_files(yyyyjjj, args.af_root,
                              args.af_modis, args.af_viirs)
    if not af_files:
        log.error("no AF files found for %d", yyyyjjj)
        return 1
    log.info("AF inputs:")
    for f in af_files:
        log.info("  %s  (%.1f MB)", f, f.stat().st_size / 1e6)

    carryover_files = []
    if args.tropical_carryover:
        carryover_files = find_carryover_modis(
            yyyyjjj, args.af_root, args.af_modis,
        )

    cfg = FinnConfig(
        tag_af=tag,
        year_rst=year_rst,
        af_files=af_files,
        out_dir=args.out_dir,
        work_dir=args.work_dir,
        raster_dir=args.raster_dir,
        date=date,
        date_definition=args.date_definition,
        summary_file=args.out_dir / f"summary_{tag}.txt",
        tropical_carryover_files=carryover_files,
        tropical_carryover_to_date=date if carryover_files else None,
        tropical_lat_bounds=tuple(args.tropics_lat_bounds),
    )
    try:
        result = run_nrt(cfg)
        log.info("✓ produced %s", result)
        log.info("elapsed: %.1f min", (time.time() - t0) / 60)
        return 0
    except Exception as e:
        log.exception("date %d failed: %s", yyyyjjj, e)
        log.info("elapsed: %.1f min", (time.time() - t0) / 60)
        return 1


if __name__ == "__main__":
    sys.exit(main())
