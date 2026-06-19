#!/usr/bin/env python3
"""Build global MODIS LCT + VCF mosaics for one or more years.

Run as a single multi-year job, or as one job per year (job array) by
passing --year.  Years are independent so job-array parallelism is the
easy win: 24 simultaneous PBS jobs finish in roughly the wall-clock
time of one year.

Examples
--------
  # Single year, locally
  ./build_global_rasters.py --year 2020 \\
      --lct-root /glade/derecho/scratch/.../MCD12Q1 \\
      --vcf-root /glade/derecho/scratch/.../MOD44B \\
      --raster-dir /glade/.../rasters_global

  # All years 2002-2025 sequentially (one big job)
  ./build_global_rasters.py --year 2002-2025 [other args]

  # PBS job-array dispatch — see build_global.pbs.template
  ./build_global_rasters.py --year ${PBS_ARRAY_INDEX} [other args]

The region polygon file (``regnum.gpkg``, with a ``Region_num`` column)
must already be sitting in ``--raster-dir`` before any NRT run.

Expected directory layout for --lct-root / --vcf-root::

    <lct-root>/
      2002/  MCD12Q1.A2002001.*.hdf
      2003/  MCD12Q1.A2003001.*.hdf
      ...
    <vcf-root>/
      2002/  MOD44B.A2002065.*.hdf
      ...

If your layout is different, override --lct-glob / --vcf-glob.
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

# Allow `python build_global_rasters.py ...` from a checkout without install.
_SCRIPT = Path(__file__).resolve()
_FINN_ROOT = _SCRIPT.parent.parent           # finn_py/scripts/.. = finn_py/
if str(_FINN_ROOT) not in sys.path:
    sys.path.insert(0, str(_FINN_ROOT))

from finn_py.raster_pipeline import run_raster                       # noqa: E402

log = logging.getLogger("build_global_rasters")


# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------

def parse_year_spec(spec: str) -> list[int]:
    """Accept '2020', '2002-2025', or '2002,2005,2020'."""
    spec = spec.strip()
    if "-" in spec and "," not in spec:
        a, b = spec.split("-")
        return list(range(int(a), int(b) + 1))
    if "," in spec:
        return [int(x) for x in spec.split(",")]
    return [int(spec)]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="build_global_rasters.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    req = p.add_argument_group("required")
    req.add_argument("--year", required=True, type=parse_year_spec,
                     help="single year (2020), range (2002-2025), "
                          "or list (2002,2005,2020)")
    req.add_argument("--lct-root", required=True, type=Path,
                     help="root directory with one subdir per year of MCD12Q1 HDFs")
    req.add_argument("--vcf-root", required=True, type=Path,
                     help="root directory with one subdir per year of MOD44B HDFs")
    req.add_argument("--raster-dir", required=True, type=Path,
                     help="output dir for modlct_<year>.tif, modvcf_<year>.tif. "
                          "regnum.gpkg is NOT built here — drop it into this "
                          "directory yourself before running finn-nrt")

    p.add_argument("--lct-glob", default="{year}/MCD12Q1.A{year}*.hdf",
                   help="glob pattern under --lct-root; {year} is substituted")
    p.add_argument("--vcf-glob", default="{year}/MOD44B.A{year}*.hdf",
                   help="glob pattern under --vcf-root; {year} is substituted")
    p.add_argument("--lct-resolution-deg", type=float, default=0.005,
                   help="WGS84 pixel size for LCT mosaic (default ~500 m)")
    p.add_argument("--vcf-resolution-deg", type=float, default=0.005,
                   help="WGS84 pixel size for VCF mosaic.  Default 500 m even "
                        "though MOD44B is 250 m native — halves the global "
                        "VCF mosaic size and the zonal-mean over a 1 km "
                        "MODIS pixel is unaffected.")
    p.add_argument("--skip-existing", action="store_true", default=True,
                   help="don't rebuild a year whose outputs already exist (default)")
    p.add_argument("--no-skip-existing", action="store_false", dest="skip_existing")
    p.add_argument("--keep-intermediate", action="store_true",
                   help="keep the per-tile sinu GeoTIFFs (for debugging)")
    p.add_argument("-v", "--verbose", action="store_true")

    return p.parse_args(argv)


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

    args.raster_dir.mkdir(parents=True, exist_ok=True)

    overall_t0 = time.time()
    results: list[tuple[int, str, float]] = []  # (year, status, secs)

    for year in args.year:
        log.info("=" * 70)
        log.info("year %d  (output -> %s)", year, args.raster_dir)
        log.info("=" * 70)
        year_t0 = time.time()

        # Per-year output paths (two GeoTIFFs per year).
        out_lct = args.raster_dir / f"modlct_{year}.tif"
        out_vcf = args.raster_dir / f"modvcf_{year}.tif"

        if args.skip_existing and out_lct.exists() and out_vcf.exists():
            log.info("both outputs exist; skipping (use --no-skip-existing to override)")
            results.append((year, "skipped", 0.0))
            continue

        # Resolve input HDFs.
        lct_files = sorted(args.lct_root.glob(
            args.lct_glob.format(year=year)
        ))
        vcf_files = sorted(args.vcf_root.glob(
            args.vcf_glob.format(year=year)
        ))
        log.info("inputs: %d LCT HDFs, %d VCF HDFs", len(lct_files), len(vcf_files))

        if not lct_files:
            log.error("no LCT HDFs found under %s with pattern %s",
                      args.lct_root, args.lct_glob.format(year=year))
            results.append((year, "no_lct", 0.0))
            continue
        if not vcf_files:
            log.error("no VCF HDFs found under %s with pattern %s",
                      args.vcf_root, args.vcf_glob.format(year=year))
            results.append((year, "no_vcf", 0.0))
            continue

        try:
            run_raster(
                year_rst=year,
                lct_files=lct_files,
                vcf_files=vcf_files,
                raster_dir=args.raster_dir,
                lct_resolution_deg=args.lct_resolution_deg,
                vcf_resolution_deg=args.vcf_resolution_deg,
            )
            results.append((year, "ok", time.time() - year_t0))
        except Exception as e:
            log.exception("year %d failed: %s", year, e)
            results.append((year, f"error: {e!r}", time.time() - year_t0))

    # Summary
    log.info("")
    log.info("=" * 70)
    log.info("SUMMARY  (total wall time %.1f min)", (time.time() - overall_t0) / 60)
    log.info("=" * 70)
    log.info("year   status      elapsed")
    for year, status, secs in results:
        log.info("%4d   %-10s  %5.1f min", year, status, secs/60)

    n_bad = sum(1 for _, s, _ in results if s not in ("ok", "skipped"))
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
