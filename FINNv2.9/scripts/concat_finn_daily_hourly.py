#!/usr/bin/env python3
"""Concatenate a single daily-average day + hourly-average file, per species.

For each species that has both a daily-average NetCDF (in ``--daily-dir``)
and an hourly-average NetCDF (in ``--hourly-dir``), this tool:

  1. Picks ONE timestep from the daily file (auto-detected as the day
     before the hourly file's first date, or specified explicitly with
     ``--daily-date YYYYMMDD`` or ``--daily-index N``).
  2. Replicates that day's spatial pattern across 24 hourly rows with
     ``datesec = [0, 3600, …, 82800]`` and matching ``time`` values.
  3. Concatenates those 24 rows with the full hourly file along ``time``.
  4. Writes one output NetCDF per species — **uniform hourly cadence
     throughout**, so downstream tools that expect a consistent time
     step work unchanged.

Both inputs must be on the same unstructured (``ncol``) grid, with
emission variables of dims ``(time, ncol)``.  Files are paired by
looking at the emission variable's name inside each NetCDF — file
names don't need to match.

Example
-------
    python concat_finn_daily_hourly.py \\
        --daily-dir  /path/to/daily_avg    \\
        --hourly-dir /path/to/hourly_avg   \\
        --out-dir    /path/to/concat

With the hourly file starting 2024-02-02, this picks 2024-02-01 from
each daily file, replicates it to 24 hourly rows, and prepends to the
hourly data — giving a Feb 1 → Mar 31 output with 24 timesteps per day.

Output filenames default to ``<species>_concat.nc``; override with
``--filename-template`` — placeholders ``{species}``, ``{daily_stem}``,
``{hourly_stem}`` are substituted.
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import re
import sys
from pathlib import Path

import numpy as np
import xarray as xr

log = logging.getLogger("concat_finn_daily_hourly")

# Variable names that are NEVER the emission species itself.
_META_VARS = frozenset({
    "date", "datesec", "time", "hour",
    "lat", "lon", "latitude", "longitude",
    "area", "grid_area", "cell_area",
    "ncol", "rrfac",                        # non-emission ncol-only vars
})

# Variable names whose value we treat as generic — meaning the actual
# species has to come from the variable's attributes, not the name.
_GENERIC_EMISSION_NAMES = frozenset({
    "fire", "emis", "emission", "emissions", "flux", "emis_flux",
})

# Regex to pull the species out of a long_name like "CO fire emissions"
# or "PM2.5 emissions"; the species is the leading whitespace-terminated
# token.  Accepts letters, digits, underscores and dots (for e.g. PM2.5).
_LONGNAME_SPECIES_RE = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9_.]*)\s+(?:fire\s+)?emissions?\b",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="concat_finn_daily_hourly.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("--daily-dir", type=Path, required=True,
                   help="directory of daily-average FINN NetCDFs")
    p.add_argument("--hourly-dir", type=Path, required=True,
                   help="directory of hourly-average FINN NetCDFs")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output directory (created if missing)")
    selector = p.add_mutually_exclusive_group()
    selector.add_argument("--daily-date", type=int, default=None,
                           metavar="YYYYMMDD",
                           help="pick this specific date from each daily file "
                                "(default: auto-detected as the day before "
                                "the hourly file's first date)")
    selector.add_argument("--daily-index", type=int, default=None, metavar="N",
                           help="pick this 0-based time index from each daily "
                                "file (alternative to --daily-date)")
    p.add_argument("--species", nargs="+", default=None,
                   help="subset of species to process (default: all species "
                        "found in BOTH directories)")
    p.add_argument("--filename-template",
                   default="{species}_concat.nc",
                   help="output filename template using {species}, "
                        "{daily_stem}, {hourly_stem} (default: %(default)s)")
    p.add_argument("--overwrite", action="store_true",
                   help="overwrite existing output files (default: skip)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Species discovery
# ---------------------------------------------------------------------------

def identify_species(ds: xr.Dataset) -> str | None:
    """Return a species key from a dataset, or None if no emission variable.

    The emission variable is any data variable on ``(time, ncol)`` that
    isn't in ``_META_VARS``.  Its NAME is often the species (e.g. ``CO``,
    ``BC``, ``num_a1_so4``), in which case that's the key.

    But some FINN files use a **generic** name — usually ``fire`` — for
    the emission variable itself, and put the actual species in
    attributes.  When we detect a generic name we look at, in order:

      1. ``long_name`` — matches patterns like ``"CO fire emissions"``
         or ``"NOx_asNO emissions"`` (species is the leading token).
      2. ``map`` — matches patterns like ``"CO->CO"`` or
         ``"NO2->NOx_asNO"`` (species is the token after ``->``, i.e.
         the output/target species).

    If neither yields a species, we fall back to the variable name so
    the caller still gets *some* key — but pairing across a directory
    of generic-named files will then collapse to a single entry.
    """
    candidates = [
        v for v in ds.data_vars
        if v not in _META_VARS
        and "time" in ds[v].dims and "ncol" in ds[v].dims
    ]
    if not candidates:
        return None
    if len(candidates) > 1:
        log.warning("  multiple species-like variables found (%s); using %s",
                    candidates, candidates[0])
    var_name = candidates[0]

    if var_name.lower() not in _GENERIC_EMISSION_NAMES:
        return var_name

    var = ds[var_name]
    # (1) long_name — the most common source
    long_name = str(var.attrs.get("long_name", ""))
    m = _LONGNAME_SPECIES_RE.match(long_name)
    if m:
        return m.group(1)

    # (2) map — some FINN conventions use "input->output"
    map_attr = str(var.attrs.get("map", ""))
    if "->" in map_attr:
        _, _, right = map_attr.partition("->")
        right = right.strip()
        if right:
            return right

    # Fall back to the (generic) var name; user will see a "duplicate
    # species" warning if multiple files land on the same key.
    log.warning("  %s: variable name %r is generic and neither long_name "
                "nor map yielded a species; will use %r as the key",
                "<dataset>", var_name, var_name)
    return var_name


def build_species_map(directory: Path) -> dict[str, Path]:
    """Scan ``directory`` and return {species_name: filepath}."""
    out: dict[str, Path] = {}
    for path in sorted(directory.glob("*.nc")):
        try:
            with xr.open_dataset(path, decode_times=False) as ds:
                sp = identify_species(ds)
        except Exception as e:
            log.warning("cannot open %s: %s", path.name, e)
            continue
        if sp is None:
            log.debug("no species variable in %s; skipping", path.name)
            continue
        if sp in out:
            log.warning("duplicate species %r in %s and %s — keeping first (%s)",
                        sp, out[sp].name, path.name, out[sp].name)
            continue
        out[sp] = path
    return out


# ---------------------------------------------------------------------------
# Time / date / datesec helpers
# ---------------------------------------------------------------------------

def _parse_time_epoch(time_units_str: str) -> np.datetime64:
    """'days since 1970-01-01 00:00:00' → np.datetime64('1970-01-01')."""
    prefix = "days since "
    if not time_units_str.lower().startswith(prefix):
        raise RuntimeError(f"unsupported time units: {time_units_str!r}")
    rest = time_units_str[len(prefix):].strip()
    date_part = rest.split()[0]
    return np.datetime64(date_part)


def _yyyymmdd_int(dates_D: np.ndarray) -> np.ndarray:
    """np.datetime64[D] array → int32 array of YYYYMMDD."""
    return np.array(
        [int(np.datetime_as_string(d, unit="D").replace("-", ""))
         for d in dates_D],
        dtype="int32",
    )


def _time_to_dates(time_vals: np.ndarray, epoch: np.datetime64
                    ) -> np.ndarray:
    """Fractional days-since-epoch → int32 YYYYMMDD."""
    dates_dt = epoch + (time_vals * 86400).astype("timedelta64[s]")
    return _yyyymmdd_int(dates_dt.astype("datetime64[D]"))


def _dates_to_days_since_epoch(date_ints: np.ndarray,
                                 epoch: np.datetime64) -> np.ndarray:
    """int32 YYYYMMDD → float64 days-since-epoch (midnight of each date)."""
    epoch_D = epoch.astype("datetime64[D]")
    out = np.empty(len(date_ints), dtype="float64")
    for i, d in enumerate(date_ints.astype(int)):
        y, m, day = d // 10000, (d // 100) % 100, d % 100
        d_np = np.datetime64(f"{y:04d}-{m:02d}-{day:02d}", "D")
        out[i] = (d_np - epoch_D) / np.timedelta64(1, "D")
    return out


# ---------------------------------------------------------------------------
# Daily / hourly preparation
# ---------------------------------------------------------------------------

def _prev_yyyymmdd(date_int: int) -> int:
    """20240202 → 20240201.  Handles month/year rollover via datetime64."""
    y, m, d = date_int // 10000, (date_int // 100) % 100, date_int % 100
    d_np = np.datetime64(f"{y:04d}-{m:02d}-{d:02d}", "D") - np.timedelta64(1, "D")
    return int(np.datetime_as_string(d_np, unit="D").replace("-", ""))


def _peek_hourly(path: Path) -> tuple[int, str, str]:
    """Peek at an hourly file to get ``(first_date_int, time_units, calendar)``.

    The hourly file is treated as the reference: its time epoch and
    calendar are used for the concatenated output so the daily portion
    can be rebased to match.  This lets us mix files with different
    ``time.units`` strings (e.g. daily on 'days since 1950-01-01' and
    hourly on 'days since 1970-01-01') without touching the hourly data.
    """
    with xr.open_dataset(path, decode_times=False) as ds:
        time_units = ds["time"].attrs.get(
            "units", "days since 1970-01-01 00:00:00")
        calendar = ds["time"].attrs.get("calendar", "gregorian")
        epoch = _parse_time_epoch(time_units)
        if "date" in ds.variables:
            first_date = int(ds["date"].values[0])
        else:
            first_date = int(_time_to_dates(ds["time"].values[:1], epoch)[0])
    return first_date, time_units, calendar


# Calendars that our np.datetime64-based day arithmetic handles correctly
# (all equivalent for dates in the modern era).
_GREGORIAN_FAMILY = frozenset({"gregorian", "standard", "proleptic_gregorian"})


def _calendars_compatible(a: str, b: str) -> bool:
    return a.lower() in _GREGORIAN_FAMILY and b.lower() in _GREGORIAN_FAMILY


def _normalize_for_concat(ds: xr.Dataset, species: str) -> xr.Dataset:
    """Bring ``ds`` into a canonical shape so it concatenates cleanly.

    The two files we're merging can come from different pipelines with
    slightly different conventions.  Left alone, xr.concat trips on any
    inconsistency.  We fix the common ones here:

      * If the emission variable has a **generic** name (``fire`` etc.)
        that doesn't match the caller's ``species`` key, rename it in
        place so both datasets end up with the same variable name.
      * Promote ``lat`` and ``lon`` from data_vars to coordinates
        (some files use ``variable:coordinates = "lat lon"`` for this;
        others don't).  ``xr.concat`` errors out when a name is a coord
        in one dataset and a data var in the other.
      * Drop the ``ncol`` dim-coord variable if present — it's usually
        just an integer index and often only exists in one of the files.
      * Drop ``area`` / ``rrfac`` if present — same reasoning, they only
        appear on one side.  If you need cell areas for later diagnostics,
        get them from the SCRIP file, which is the authoritative source.
    """
    # 1. Rename generic emission variable to the actual species key
    for name in list(ds.data_vars):
        if name.lower() in _GENERIC_EMISSION_NAMES and name != species:
            log.debug("  renaming emission variable %r → %r", name, species)
            ds = ds.rename({name: species})
            break

    # 2. Promote lat/lon to coordinates
    to_promote = [v for v in ("lat", "lon", "latitude", "longitude")
                    if v in ds.data_vars]
    if to_promote:
        ds = ds.set_coords(to_promote)

    # 3. Drop nuisance variables that only appear on some files
    for name in ("ncol", "area", "rrfac"):
        if name in ds.variables:
            ds = ds.drop_vars(name)

    return ds


def prepare_daily(path: Path, *, target_date: int | None,
                  target_index: int | None, hourly_first_date: int,
                  hourly_time_units: str, hourly_calendar: str,
                  species: str
                  ) -> tuple[xr.Dataset, int]:
    """Load a daily file, pick a single row, and expand it into 24 hourly rows.

    Selection precedence:

      * ``target_index`` if not None → use that row directly.
      * else ``target_date`` if not None → look it up in ``date``.
      * else use ``hourly_first_date - 1`` (the day right before hourly starts).

    The reconstructed ``time`` coordinate uses **the hourly file's**
    ``units`` string and calendar, so the daily and hourly datasets are
    guaranteed to be compatible for ``xr.concat`` even if they were
    written against different epochs.

    Returns ``(expanded_dataset, chosen_date_int)``.  The expanded
    dataset has ``time`` of length 24, with:

      * ``time`` = day_offset (from hourly's epoch) + [0/24, 1/24, …, 23/24]
      * ``date`` = 24 copies of the chosen YYYYMMDD
      * ``datesec`` = 0, 3600, 7200, …, 82800
      * species variable = 24 identical copies of the chosen row
    """
    ds = xr.open_dataset(path, decode_times=False)
    nt = ds.sizes["time"]
    daily_units = ds["time"].attrs.get("units", "days since 1970-01-01 00:00:00")
    daily_calendar = ds["time"].attrs.get("calendar", "gregorian")
    daily_epoch = _parse_time_epoch(daily_units)

    # Calendar compatibility check — the day-arithmetic below is only
    # correct for the gregorian family.  Anything exotic (360_day,
    # noleap, julian, …) needs a proper cftime conversion; refuse rather
    # than silently produce wrong times.
    if not _calendars_compatible(daily_calendar, hourly_calendar):
        raise RuntimeError(
            f"{path.name}: calendar mismatch — daily={daily_calendar!r} "
            f"vs hourly={hourly_calendar!r}.  This tool assumes the "
            "gregorian family (gregorian/standard/proleptic_gregorian).  "
            "Convert one of the files to a gregorian-family calendar "
            "first (e.g. `cdo setcalendar,proleptic_gregorian in.nc out.nc`).")
    if daily_units != hourly_time_units:
        log.info("  rebasing daily time onto hourly epoch: "
                 "%r → %r", daily_units, hourly_time_units)

    # Get or derive the date array (used to locate the chosen row)
    if "date" in ds.variables:
        date_vals = ds["date"].values.astype(int)
    else:
        date_vals = _time_to_dates(ds["time"].values, daily_epoch).astype(int)

    # Decide which row to pick
    if target_index is not None:
        if not (0 <= target_index < nt):
            raise ValueError(f"{path.name}: --daily-index {target_index} out "
                             f"of range (file has {nt} timesteps)")
        idx = int(target_index)
    else:
        want = target_date if target_date is not None else _prev_yyyymmdd(hourly_first_date)
        matches = np.where(date_vals == want)[0]
        if matches.size == 0:
            raise ValueError(
                f"{path.name}: date {want} not found in daily file "
                f"(available: {int(date_vals[0])}..{int(date_vals[-1])})")
        idx = int(matches[0])

    chosen_date = int(date_vals[idx])
    log.info("  picked daily index %d (date=%d) → replicating to 24 hours",
             idx, chosen_date)

    # Extract single row (time dim size 1), then concat-replicate 24 times
    single = ds.isel(time=[idx])
    expanded = xr.concat(
        [single] * 24, dim="time",
        data_vars="minimal", coords="minimal",
        compat="override", join="exact",
    )

    # Compute the correct time / date / datesec for the 24 hourly rows.
    # Use HOURLY's epoch so the two datasets align.
    hourly_epoch = _parse_time_epoch(hourly_time_units)
    y, m, d = chosen_date // 10000, (chosen_date // 100) % 100, chosen_date % 100
    day_zero = np.datetime64(f"{y:04d}-{m:02d}-{d:02d}", "D")
    day_offset = float((day_zero - hourly_epoch.astype("datetime64[D]"))
                        / np.timedelta64(1, "D"))
    time_new    = (day_offset + np.arange(24) / 24.0).astype("float64")
    date_new    = np.full(24, chosen_date, dtype="int32")
    datesec_new = (np.arange(24) * 3600).astype("int32")

    # Attach the new time coord — with HOURLY's units + calendar
    time_attrs = {
        "units":     hourly_time_units,
        "calendar":  hourly_calendar,
        "long_name": "time",
        "standard_name": "time",
    }
    expanded = expanded.assign_coords(
        time=(("time",), time_new, time_attrs),
    )
    expanded = expanded.drop_vars("date", errors="ignore").assign(
        date=(("time",), date_new,
               {"long_name": "current date (YYYYMMDD)"}),
    )
    expanded = expanded.drop_vars("datesec", errors="ignore").assign(
        datesec=(("time",), datesec_new,
                  {"long_name": "current seconds of current date",
                   "units": "s"}),
    )

    # Final canonicalisation (rename generic emission var, promote
    # lat/lon to coords, drop nuisance vars).
    expanded = _normalize_for_concat(expanded, species)
    return expanded, chosen_date


def prepare_hourly(path: Path, *, species: str) -> xr.Dataset:
    """Open an hourly-average file; ensure ``date`` and ``datesec`` are
    both present (compute from ``time`` if not)."""
    ds = xr.open_dataset(path, decode_times=False)
    time_units = ds["time"].attrs.get("units", "days since 1970-01-01 00:00:00")
    epoch = _parse_time_epoch(time_units)

    if "date" not in ds.variables:
        date_vals = _time_to_dates(ds["time"].values, epoch)
        ds = ds.assign(date=(("time",), date_vals,
                              {"long_name": "current date (YYYYMMDD)"}))

    if "datesec" not in ds.variables:
        dates_dt = epoch + (ds["time"].values * 86400).astype("timedelta64[s]")
        datesec_vals = np.array(
            [int((d - d.astype("datetime64[D]"))
                    .astype("timedelta64[s]").astype(int))
             for d in dates_dt],
            dtype="int32",
        )
        ds = ds.assign(datesec=(("time",), datesec_vals,
                                  {"long_name": "current seconds of current date",
                                   "units": "s"}))
    return _normalize_for_concat(ds, species)


# ---------------------------------------------------------------------------
# Concatenation
# ---------------------------------------------------------------------------

def concat_and_write(daily_ds: xr.Dataset, hourly_ds: xr.Dataset, *,
                     out_path: Path, species: str,
                     daily_src: Path, hourly_src: Path,
                     seam_date: int) -> None:
    """Concat daily+hourly along time, write NetCDF."""
    # Sanity checks
    if daily_ds.sizes["ncol"] != hourly_ds.sizes["ncol"]:
        raise RuntimeError(
            f"ncol mismatch: daily={daily_ds.sizes['ncol']} "
            f"hourly={hourly_ds.sizes['ncol']}")

    d_units = daily_ds["time"].attrs.get("units", "")
    h_units = hourly_ds["time"].attrs.get("units", "")
    if d_units != h_units:
        # This should never happen — prepare_daily rebases to hourly's epoch.
        # If it does, something upstream is out of sync.
        raise RuntimeError(
            f"internal error: time units still differ after rebase "
            f"(daily={d_units!r}, hourly={h_units!r})")

    d_last  = float(daily_ds["time"].values[-1])
    h_first = float(hourly_ds["time"].values[0])
    if d_last >= h_first:
        log.warning(
            "  daily-last time (%.4f) >= hourly-first time (%.4f); output "
            "may have overlapping/out-of-order timesteps", d_last, h_first)

    # Concat.  ``coords='minimal'`` and ``data_vars='minimal'`` ensure only
    # variables that already contain the ``time`` dimension get concatenated;
    # the ``lat(ncol)`` and ``lon(ncol)`` coordinates stay attached to the
    # spatial dim (with ``compat='override'`` we skip an exhaustive bit-
    # equality check — the ncol-size match above and the exact-join below
    # are sufficient guardrails).
    combo = xr.concat(
        [daily_ds, hourly_ds],
        dim="time",
        data_vars="minimal",
        coords="minimal",
        compat="override",
        combine_attrs="drop_conflicts",
        join="exact",
    )

    # Force clean dtypes for the time-metadata columns
    combo["date"]    = combo["date"].astype("int32")
    combo["datesec"] = combo["datesec"].astype("int32")

    # Global attributes describing the concat
    combo.attrs = {**daily_ds.attrs, **combo.attrs}
    combo.attrs["concat_note"] = (
        f"first 24 timesteps replicate daily-average day {seam_date} "
        f"across hours 0..23 (datesec 0..82800); remaining "
        f"{hourly_ds.sizes['time']} timesteps are hourly averages from "
        f"the hourly file — uniform hourly cadence throughout."
    )
    combo.attrs["concat_daily_source"]  = daily_src.name
    combo.attrs["concat_hourly_source"] = hourly_src.name
    combo.attrs["concat_seam_date"]     = np.int32(seam_date)
    combo.attrs["concat_created"] = dt.datetime.utcnow().strftime(
        "%Y-%m-%dT%H:%M:%SZ")

    encoding = {species: {"zlib": True, "complevel": 4}}
    combo.to_netcdf(out_path, format="NETCDF4", encoding=encoding,
                    unlimited_dims=["time"])
    log.info("  wrote %s   (time=%d: 24 replicated daily + %d hourly)",
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

    for d in (args.daily_dir, args.hourly_dir):
        if not d.is_dir():
            log.error("not a directory: %s", d); return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log.info("scanning %s ...", args.daily_dir)
    daily_map = build_species_map(args.daily_dir)
    log.info("  found %d species: %s", len(daily_map), sorted(daily_map))
    log.info("scanning %s ...", args.hourly_dir)
    hourly_map = build_species_map(args.hourly_dir)
    log.info("  found %d species: %s", len(hourly_map), sorted(hourly_map))

    if args.species:
        target = list(args.species)
    else:
        target = sorted(set(daily_map) & set(hourly_map))
    log.info("processing %d species", len(target))

    ok, fail, skipped = 0, 0, 0
    for sp in target:
        log.info("--- %s ---", sp)
        if sp not in daily_map:
            log.warning("  no daily file for %s — skipping", sp); fail += 1; continue
        if sp not in hourly_map:
            log.warning("  no hourly file for %s — skipping", sp); fail += 1; continue

        out_name = args.filename_template.format(
            species=sp,
            daily_stem=daily_map[sp].stem,
            hourly_stem=hourly_map[sp].stem,
        )
        out_path = args.out_dir / out_name
        if out_path.exists() and not args.overwrite:
            log.info("  %s already exists — skipping (use --overwrite)",
                     out_path.name)
            skipped += 1
            continue

        try:
            hourly_first_dt, hourly_units, hourly_cal = _peek_hourly(hourly_map[sp])
            daily_ds, chosen_date = prepare_daily(
                daily_map[sp],
                target_date=args.daily_date,
                target_index=args.daily_index,
                hourly_first_date=hourly_first_dt,
                hourly_time_units=hourly_units,
                hourly_calendar=hourly_cal,
                species=sp,
            )
            hourly_ds = prepare_hourly(hourly_map[sp], species=sp)
            concat_and_write(
                daily_ds, hourly_ds,
                out_path=out_path, species=sp,
                daily_src=daily_map[sp], hourly_src=hourly_map[sp],
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
