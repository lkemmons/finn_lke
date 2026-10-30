#!/usr/bin/env python3
"""Append one daily-average day (24 replicated hourly rows) to the END
of each hourly-cadence FINN NetCDF.

Symmetric with ``concat_finn_daily_hourly.py`` — that tool prepends a
replicated daily day at the *start*; this one appends it at the *end*.

For each species with an hourly-cadence NetCDF in ``--hourly-dir``
(typically the output of ``concat_finn_daily_hourly.py``) and a
daily-average NetCDF in ``--daily-dir``:

  1. Determines the day to add — defaults to the day AFTER the hourly
     file's last date; override with ``--daily-date YYYYMMDD`` or
     ``--daily-index N``.
  2. Extracts that row from the daily file, rebases its time onto the
     hourly file's epoch/calendar, and replicates it across 24 hourly
     slots (``datesec = [0, 3600, ..., 82800]``).
  3. Concatenates those 24 rows to the end of the hourly file.
  4. Writes one output NetCDF per species.

Handles the same edge cases as the concat tool: different time-unit
epochs, generic ``fire``-style emission variable names, lat/lon
coord-vs-data-var inconsistencies.  All that logic is imported
directly from ``concat_finn_daily_hourly.py``.

Example
-------
    python append_daily_to_hourly.py \\
        --hourly-dir /path/to/concat_outputs   \\
        --daily-dir  /path/to/daily_avg        \\
        --out-dir    /path/to/extended

If the hourly file ends 2024-03-31 and the daily file contains
2024-04-01, this appends 24 hourly rows for Apr 1 to the end.
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
from pathlib import Path

import numpy as np
import xarray as xr

# Import shared helpers from the sister concat tool
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))
from concat_finn_daily_hourly import (
    build_species_map,
    prepare_daily, prepare_hourly,
    _parse_time_epoch, _time_to_dates,
)

log = logging.getLogger("append_daily_to_hourly")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="append_daily_to_hourly.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("--hourly-dir", type=Path, required=True,
                   help="directory of hourly-cadence NetCDFs to extend "
                        "(typically outputs of concat_finn_daily_hourly.py)")
    p.add_argument("--daily-dir", type=Path, required=True,
                   help="directory of daily-average NetCDFs — source of "
                        "the day being appended")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output directory (created if missing)")
    selector = p.add_mutually_exclusive_group()
    selector.add_argument("--daily-date", type=int, default=None,
                           metavar="YYYYMMDD",
                           help="pick this specific date from each daily "
                                "file (default: auto-detected as the day "
                                "AFTER the hourly file's last date)")
    selector.add_argument("--daily-index", type=int, default=None, metavar="N",
                           help="pick this 0-based time index from each "
                                "daily file (alternative to --daily-date)")
    p.add_argument("--species", nargs="+", default=None,
                   help="subset of species to process (default: all species "
                        "found in BOTH directories)")
    p.add_argument("--filename-template",
                   default="{species}_extended.nc",
                   help="output filename template using {species}, "
                        "{hourly_stem}, {daily_stem} (default: %(default)s)")
    p.add_argument("--overwrite", action="store_true",
                   help="overwrite existing output files (default: skip)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _next_yyyymmdd(date_int: int) -> int:
    """``20240331`` → ``20240401``.  Handles month/year rollover."""
    y, m, d = date_int // 10000, (date_int // 100) % 100, date_int % 100
    d_np = np.datetime64(f"{y:04d}-{m:02d}-{d:02d}", "D") + np.timedelta64(1, "D")
    return int(np.datetime_as_string(d_np, unit="D").replace("-", ""))


def _peek_hourly_end(path: Path) -> tuple[int, str, str]:
    """Return ``(last_date_int, time_units, calendar)`` from the hourly file."""
    with xr.open_dataset(path, decode_times=False) as ds:
        time_units = ds["time"].attrs.get(
            "units", "days since 1970-01-01 00:00:00")
        calendar = ds["time"].attrs.get("calendar", "gregorian")
        epoch = _parse_time_epoch(time_units)
        if "date" in ds.variables:
            last_date = int(ds["date"].values[-1])
        else:
            last_date = int(_time_to_dates(ds["time"].values[-1:], epoch)[0])
    return last_date, time_units, calendar


# ---------------------------------------------------------------------------
# Append + write
# ---------------------------------------------------------------------------

def append_and_write(hourly_ds: xr.Dataset, daily_ds: xr.Dataset, *,
                     out_path: Path, species: str,
                     hourly_src: Path, daily_src: Path,
                     seam_date: int) -> None:
    """Concat ``[hourly_ds, daily_ds]`` along time (daily at end); write NetCDF."""
    if hourly_ds.sizes["ncol"] != daily_ds.sizes["ncol"]:
        raise RuntimeError(
            f"ncol mismatch: hourly={hourly_ds.sizes['ncol']} "
            f"daily={daily_ds.sizes['ncol']}")

    h_units = hourly_ds["time"].attrs.get("units", "")
    d_units = daily_ds["time"].attrs.get("units", "")
    if d_units != h_units:
        raise RuntimeError(
            f"internal error: time units differ after rebase "
            f"(hourly={h_units!r}, daily={d_units!r})")

    h_last  = float(hourly_ds["time"].values[-1])
    d_first = float(daily_ds["time"].values[0])
    if d_first <= h_last:
        log.warning(
            "  appended day's first time (%.4f) <= hourly last time "
            "(%.4f); output may have overlapping/out-of-order timesteps",
            d_first, h_last)

    # HOURLY first, then DAILY — that's the whole difference from concat.
    combo = xr.concat(
        [hourly_ds, daily_ds],
        dim="time",
        data_vars="minimal",
        coords="minimal",
        compat="override",
        combine_attrs="drop_conflicts",
        join="exact",
    )

    combo["date"]    = combo["date"].astype("int32")
    combo["datesec"] = combo["datesec"].astype("int32")

    # Global attrs — carry hourly's forward, add append metadata
    combo.attrs = {**hourly_ds.attrs, **combo.attrs}
    combo.attrs["append_note"] = (
        f"appended 24 timesteps from daily-average day {seam_date} "
        f"(datesec 0..82800) to the end of the hourly file — uniform "
        f"hourly cadence maintained throughout."
    )
    combo.attrs["append_hourly_source"] = hourly_src.name
    combo.attrs["append_daily_source"]  = daily_src.name
    combo.attrs["append_seam_date"]     = np.int32(seam_date)
    combo.attrs["append_created"] = dt.datetime.utcnow().strftime(
        "%Y-%m-%dT%H:%M:%SZ")

    encoding = {species: {"zlib": True, "complevel": 4}}
    combo.to_netcdf(out_path, format="NETCDF4", encoding=encoding,
                    unlimited_dims=["time"])
    log.info("  wrote %s   (time=%d: %d existing + 24 appended daily)",
             out_path.name, combo.sizes["time"], hourly_ds.sizes["time"])


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

    for d in (args.hourly_dir, args.daily_dir):
        if not d.is_dir():
            log.error("not a directory: %s", d); return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log.info("scanning %s ...", args.hourly_dir)
    hourly_map = build_species_map(args.hourly_dir)
    log.info("  found %d species: %s", len(hourly_map), sorted(hourly_map))
    log.info("scanning %s ...", args.daily_dir)
    daily_map = build_species_map(args.daily_dir)
    log.info("  found %d species: %s", len(daily_map), sorted(daily_map))

    if args.species:
        target = list(args.species)
    else:
        target = sorted(set(daily_map) & set(hourly_map))
    log.info("processing %d species", len(target))

    ok, fail, skipped = 0, 0, 0
    for sp in target:
        log.info("--- %s ---", sp)
        if sp not in hourly_map:
            log.warning("  no hourly file for %s — skipping", sp); fail += 1; continue
        if sp not in daily_map:
            log.warning("  no daily file for %s — skipping", sp); fail += 1; continue

        out_name = args.filename_template.format(
            species=sp,
            hourly_stem=hourly_map[sp].stem,
            daily_stem=daily_map[sp].stem,
        )
        out_path = args.out_dir / out_name
        if out_path.exists() and not args.overwrite:
            log.info("  %s already exists — skipping (use --overwrite)",
                     out_path.name)
            skipped += 1
            continue

        try:
            hourly_last, hourly_units, hourly_cal = _peek_hourly_end(hourly_map[sp])
            log.info("  hourly last date: %d", hourly_last)

            # Determine the target date/index.  For append the *default*
            # is (hourly_last + 1 day) — we compute that here rather than
            # relying on prepare_daily's default (which is prepend-style
            # "day before the reference").
            if args.daily_date is not None:
                target_date_arg = args.daily_date
                target_index_arg = None
            elif args.daily_index is not None:
                target_date_arg = None
                target_index_arg = args.daily_index
            else:
                target_date_arg = _next_yyyymmdd(hourly_last)
                target_index_arg = None
                log.info("  auto-detected append date: %d", target_date_arg)

            daily_ds, chosen_date = prepare_daily(
                daily_map[sp],
                target_date=target_date_arg,
                target_index=target_index_arg,
                hourly_first_date=hourly_last,       # unused when target_date is set
                hourly_time_units=hourly_units,
                hourly_calendar=hourly_cal,
                species=sp,
            )
            hourly_ds = prepare_hourly(hourly_map[sp], species=sp)
            append_and_write(
                hourly_ds, daily_ds,
                out_path=out_path, species=sp,
                hourly_src=hourly_map[sp], daily_src=daily_map[sp],
                seam_date=chosen_date,
            )
            ok += 1
        except Exception as e:
            log.error("  failed: %s", e, exc_info=args.verbose)
            fail += 1

    log.info("done — %d ok, %d skipped (exists), %d failed", ok, skipped, fail)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
