#!/usr/bin/env python3
"""Run the finn_py pipeline against finalized FIRMS archive AF files.

Differences from ``run_daily_nrt.py``:

* Archive files contain many days (often a whole year) in one CSV/SHP,
  with an ``instrument`` column inside the data instead of encoded in the
  filename.  Pass the file paths verbatim with ``--af-files`` — no
  pattern substitution.
* Multiple input files can be combined in one run (e.g. MODIS C6 + VIIRS
  C2 archives for the same year).  The pipeline filters by ``--date``
  internally.
* The ``type`` column maps to ``anomtype``; ``--filter-persistent``
  (default on) excludes types 1/2/3 (volcanoes, static land sources,
  offshore) just as in the NRT pipeline.

Each run still processes one day — for a year of archive data, fan out
across days with a job array (one job per YYYYJJJ).  The driver loads
the archive files once per invocation; loading is cheap relative to step
1+2 processing.

Examples
--------
    # One day of 2020 archive data, both MODIS + VIIRS
    ./work_archive.py --date 2020196 \\
        --af-files /glade/.../fire_archive_M-C61_2020.csv \\
                   /glade/.../fire_archive_SV-C2_2020.csv \\
        --raster-dir /glade/.../rasters \\
        --work-dir   /glade/.../work \\
        --out-dir    /glade/.../out_archive

    # Submit as PBS array — see work_archive.pbs.template
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
import time
from pathlib import Path

_SCRIPT = Path(__file__).resolve()
_FINN_ROOT = _SCRIPT.parent.parent
if str(_FINN_ROOT) not in sys.path:
    sys.path.insert(0, str(_FINN_ROOT))

from finn_py.config import FinnConfig                              # noqa: E402
from finn_py.pipeline import run_nrt                               # noqa: E402

log = logging.getLogger("work_archive")


def yyyyjjj_to_date(yyyyjjj: int) -> dt.date:
    s = f"{int(yyyyjjj):07d}"
    return dt.date(int(s[:4]), 1, 1) + dt.timedelta(days=int(s[4:]) - 1)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="work_archive.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    req = p.add_argument_group("required")
    req.add_argument("--date", required=True, type=int,
                     help="processing date as YYYYJJJ (e.g. 2020196)")
    req.add_argument("--af-files", required=True, type=Path, nargs="+",
                     metavar="FILE",
                     help="one or more finalized FIRMS archive files "
                          "(.csv/.shp/.zip).  Pass the full path of each "
                          "file; no pattern substitution.")
    req.add_argument("--raster-dir", required=True, type=Path,
                     help="directory with pre-built modlct_<year>.tif, "
                          "modvcf_<year>.tif, regnum.gpkg")
    req.add_argument("--work-dir", required=True, type=Path,
                     help="root working directory")
    req.add_argument("--out-dir",  required=True, type=Path,
                     help="output directory")

    p.add_argument("--tag-prefix", default="arc",
                   help="output tag = '<prefix>_<YYYYJJJ>'  (default: %(default)s)")
    p.add_argument("--year-rst", type=int, default=None,
                   help="year for raster file selection "
                        "(default: year derived from the date)")
    p.add_argument("--date-definition", default="UTC", choices=["UTC", "LST"],
                   help="how to assign each detection to a calendar day "
                        "(default: %(default)s)")
    p.add_argument("--filter-persistent", action="store_true", default=True,
                   help="drop persistent / non-vegetation anomtypes "
                        "(volcano, static land source, offshore).  "
                        "Default on; archive's `type` column drives this.")
    p.add_argument("--no-filter-persistent", action="store_false",
                   dest="filter_persistent")
    p.add_argument("--tropical-carryover-file", type=Path, default=None,
                   metavar="FILE",
                   help="optional archive file for the previous day's "
                        "MODIS detections (used only if you're running "
                        "single-day from per-day archive slices; with one "
                        "big multi-month archive file this isn't needed "
                        "because the previous day's data is already in "
                        "the same file and step1.prep duplicates "
                        "tropical MODIS into the next day automatically).")
    p.add_argument("--tropics-lat-bounds", nargs=2, type=float,
                   default=[-23.5, 23.5],
                   metavar=("LAT_MIN", "LAT_MAX"),
                   help="latitude bounds for --tropical-carryover-file "
                        "(default: -23.5 23.5)")
    p.add_argument("--skip-existing", action="store_true", default=True,
                   help="skip if output file already exists (default)")
    p.add_argument("--no-skip-existing", action="store_false",
                   dest="skip_existing")
    p.add_argument("-v", "--verbose", action="store_true")

    return p.parse_args(argv)


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
    log.info("DATE %d  (%s)   archive run", yyyyjjj, date.isoformat())
    log.info("=" * 70)

    # Skip if any output for this tag already exists.
    existing = list(args.out_dir.glob(f"out_{tag}_*.gpkg"))
    if args.skip_existing and existing:
        log.info("output exists; skipping: %s", existing[0].name)
        return 0

    # Validate inputs upfront so we fail before the expensive load.
    missing = [p for p in args.af_files if not p.exists()]
    if missing:
        for p in missing:
            log.error("AF file not found: %s", p)
        return 1
    log.info("AF inputs:")
    for f in args.af_files:
        log.info("  %s  (%.1f MB)", f, f.stat().st_size / 1e6)

    carryover_files = []
    if args.tropical_carryover_file:
        if not args.tropical_carryover_file.exists():
            log.error("tropical-carryover file not found: %s",
                       args.tropical_carryover_file)
            return 1
        carryover_files = [args.tropical_carryover_file]
        log.info("tropical carryover: %s", args.tropical_carryover_file)

    cfg = FinnConfig(
        tag_af=tag,
        year_rst=year_rst,
        af_files=list(args.af_files),
        out_dir=args.out_dir,
        work_dir=args.work_dir,
        raster_dir=args.raster_dir,
        date=date,
        date_definition=args.date_definition,
        filter_persistent_sources=args.filter_persistent,
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
