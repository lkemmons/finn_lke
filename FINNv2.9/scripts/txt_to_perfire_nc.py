#!/usr/bin/env python3
"""Convert legacy FINN text CSV to per-fire NetCDF (unstructured by polyid).

Distinct from grid_txt_emissions.py, which grids each species onto a
regular or SCRIP grid.  This tool keeps the per-fire granularity:
every unique POLYID becomes one entry along a polyid dimension.

Input format
------------
Legacy FINN text output.  Columns:

  DAY, FIREID, POLYID, GENVEG, LATI, LONGI, AREA, BMASS, FRP,
  <species1>, <species2>, ..., [HOUR], datetimeUTC, datetimeLT,
  Country, TrendCountry

Country and TrendCountry are always dropped.  The HOUR column is
also dropped and its VALUES are IGNORED — HOUR is unreliable in some
FINN outputs, so the actual UTC hour is derived from datetimeUTC.

Mode is determined from HOUR column presence (a stable structural
property of the file, independent of HOUR's unreliable values):

  * HOUR column present -> hourly: output dims (polyid, hour), with
                                    hour indices derived from datetimeUTC.
  * HOUR column absent  -> daily: output dim (polyid,).

Override with --daily to force daily even when HOUR was present.

Species values are NOT unit-converted; they preserve whatever units the
source CSV uses (typically mol/hour for gases, kg/hour for aerosols,
number/hour for particle counts).  If an --mw-table is supplied, species
variables get units, species_type and molecular_weight_g_per_mol
attributes.
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

log = logging.getLogger("txt_to_perfire_nc")

# Columns that are never emission species and never kept as fire metadata.
_DROP_COLS = frozenset({
    "Country", "TrendCountry",
    "DAY", "datetimeUTC", "datetimeLT",   # redundant with date + hour
})

# Per-fire metadata columns (become (polyid,) vars in the output).
_FIRE_META_COLS = ("FIREID", "GENVEG", "LATI", "LONGI", "AREA", "BMASS", "FRP")

# Everything that's NOT a species column in the CSV.
_ALL_NON_SPECIES = frozenset({"POLYID", "HOUR"} | set(_DROP_COLS) | set(_FIRE_META_COLS))


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="txt_to_perfire_nc.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("inputs", type=Path, nargs="+",
                   help="one or more legacy FINN text CSV files")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output directory (created if missing)")
    p.add_argument("--mw-table", type=Path, default=None,
                   help="optional CSV with columns species,type,mw")
    p.add_argument("--daily", action="store_true",
                   help="force daily mode even when the input has multiple "
                        "rows per POLYID (default: auto-detected as hourly "
                        "if any POLYID has more than one row)")
    p.add_argument("--overwrite", action="store_true",
                   help="overwrite existing output files")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def read_csv(path: Path) -> pd.DataFrame:
    log.info("reading %s", path)
    df = pd.read_csv(path)
    for col in ("Country", "TrendCountry"):
        if col in df.columns:
            df = df.drop(columns=col)
    # HOUR column is deliberately NOT dropped here — the main() loop
    # uses its presence as the hourly-vs-daily mode signal (that's a
    # stable structural property of the file).  Its VALUES are
    # unreliable and never trusted; the actual UTC hour comes from
    # datetimeUTC.  HOUR is dropped in main() right after detection so
    # it can't leak into the output.
    for req in ("POLYID", "LATI", "LONGI"):
        if req not in df.columns:
            raise RuntimeError(f"{path}: missing required column {req!r}")
    log.info("  %d rows, %d columns", len(df), len(df.columns))
    return df


def read_mw_table(path: Path | None) -> dict:
    if path is None:
        return {}
    df = pd.read_csv(path, comment="#", skip_blank_lines=True)
    df.columns = [c.strip() for c in df.columns]
    if "species" not in df.columns or "type" not in df.columns:
        log.warning("MW table missing 'species' or 'type' column; ignoring")
        return {}
    if "mw" not in df.columns:
        df["mw"] = np.nan
    df = df.dropna(subset=["species", "type"], how="any")
    df["species"] = df["species"].astype(str).str.strip()
    df["type"]    = df["type"].astype(str).str.strip().str.lower()
    out = {}
    for _, r in df.iterrows():
        if not r["species"]:
            continue
        out[r["species"]] = {
            "type": r["type"],
            "mw":   float(r["mw"]) if pd.notna(r["mw"]) else float("nan"),
        }
    log.info("MW table: %d species (%s)", len(out), path.name)
    return out


def resolve_date(df: pd.DataFrame) -> tuple[int, str]:
    if "datetimeUTC" in df.columns and df["datetimeUTC"].notna().any():
        iso = str(df["datetimeUTC"].dropna().iloc[0])[:10]
        y, m, d = int(iso[:4]), int(iso[5:7]), int(iso[8:10])
        return y * 10000 + m * 100 + d, iso
    if "DAY" in df.columns:
        day = int(df["DAY"].iloc[0])
        year = dt.date.today().year
        d0 = dt.date(year, 1, 1) + dt.timedelta(days=day - 1)
        log.warning("using current year %d for Julian DAY=%d", year, day)
        return d0.year * 10000 + d0.month * 100 + d0.day, d0.isoformat()
    raise RuntimeError("no date info found (no datetimeUTC or DAY column)")


def _species_attrs(sp: str, mw_lookup: dict, per_time: str) -> dict:
    """per_time is 'hour' or 'day'."""
    if sp in mw_lookup:
        info = mw_lookup[sp]
        typ = info["type"]
        u_map = {"gas": f"mol/{per_time}", "aerosol": f"kg/{per_time}",
                 "number": f"number/{per_time}"}
        attrs = {"units": u_map.get(typ, f"per {per_time} (unspecified)"),
                 "species_type": typ,
                 "long_name": f"{sp} emission"}
        if np.isfinite(info.get("mw", np.nan)):
            attrs["molecular_weight_g_per_mol"] = float(info["mw"])
        return attrs
    return {
        "long_name": f"{sp} emission",
        "units": f"per {per_time} (source-file units; mol/{per_time} for gases, "
                 f"kg/{per_time} for aerosols, number/{per_time} for counts)",
    }


_FIRE_META_ATTRS = {
    "FIREID": {"long_name": "FINN fire ID"},
    "GENVEG": {"long_name": "generic vegetation type"},
    "LATI":   {"units": "degrees_north", "standard_name": "latitude"},
    "LONGI":  {"units": "degrees_east",  "standard_name": "longitude"},
    "AREA":   {"units": "m2", "long_name": "fire area"},
    "BMASS":  {"units": "kg m-2", "long_name": "biomass burned per unit area"},
    "FRP":    {"units": "MW", "long_name": "fire radiative power"},
}

_FIRE_META_RENAME = {
    "FIREID": "fireid", "GENVEG": "genveg",
    "LATI": "lat", "LONGI": "lon",
    "AREA": "area", "BMASS": "bmass", "FRP": "frp",
}


def _fire_meta_arrays(df: pd.DataFrame, first_idx: np.ndarray) -> dict:
    """Per-fire metadata (one value per polyid) from first-occurrence rows."""
    out = {}
    for col in _FIRE_META_COLS:
        if col in df.columns:
            out[col] = df[col].values[first_idx]
    return out


def build_hourly_ds(df: pd.DataFrame, mw_lookup: dict, source: str) -> xr.Dataset:
    # HOUR values in some FINN outputs are unreliable, so we derive the
    # UTC hour from datetimeUTC instead.  HOUR itself is ignored (and
    # never written to the output).
    if "datetimeUTC" not in df.columns:
        raise RuntimeError("hourly mode requires datetimeUTC column "
                            "(HOUR column is not trusted and cannot be used)")
    species_cols = [c for c in df.columns if c not in _ALL_NON_SPECIES]
    log.info("  hourly mode: %d species columns", len(species_cols))

    poly_vals = df["POLYID"].values.astype("int64")
    polyids, first_idx = np.unique(poly_vals, return_index=True)
    polyids = polyids.astype("int32")
    n_poly = len(polyids)
    log.info("  %d unique fires (POLYIDs)", n_poly)

    row_poly_idx = np.searchsorted(polyids, poly_vals)
    row_hour_idx = pd.to_datetime(df["datetimeUTC"]).dt.hour.values.astype("int64")

    valid = (row_hour_idx >= 0) & (row_hour_idx < 24)
    dropped = int((~valid).sum())
    if dropped:
        log.warning("  %d rows with UTC hour outside [0,23] dropped", dropped)
        row_poly_idx = row_poly_idx[valid]
        row_hour_idx = row_hour_idx[valid]

    combined = row_poly_idx.astype("int64") * 24 + row_hour_idx
    if len(np.unique(combined)) < len(combined):
        log.warning("  duplicate (POLYID, UTC-hour) rows found -- values will be summed")

    species_arrays = {}
    for sp in species_cols:
        arr = np.zeros((n_poly, 24), dtype="float32")
        vals = df[sp].values[valid].astype("float32") if dropped else \
               df[sp].values.astype("float32")
        np.add.at(arr, (row_poly_idx, row_hour_idx), vals)
        species_arrays[sp] = arr

    fire_meta = _fire_meta_arrays(df, first_idx)
    date_int, iso_date = resolve_date(df)
    log.info("  date: %d (%s)", date_int, iso_date)

    data_vars = {}
    for sp in species_cols:
        data_vars[sp] = (("polyid", "hour"), species_arrays[sp],
                          _species_attrs(sp, mw_lookup, "hour"))
    for col, vals in fire_meta.items():
        data_vars[_FIRE_META_RENAME[col]] = (
            ("polyid",), vals, _FIRE_META_ATTRS[col])
    data_vars["date"] = ((), np.int32(date_int),
                         {"long_name": "date as YYYYMMDD"})

    return xr.Dataset(
        data_vars=data_vars,
        coords={
            "polyid": (("polyid",), polyids, {"long_name": "FINN polygon ID"}),
            "hour":   (("hour",), np.arange(24, dtype="int32"),
                        {"long_name": "hour of day (UTC)", "units": "hour"}),
        },
        attrs={
            "title":          "FINN per-fire hourly emissions",
            "source":         source,
            "date_processed": iso_date,
            "Conventions":    "CF-1.8",
            "created":        dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "comment":        "One entry per (POLYID, UTC-hour). Species "
                              "values preserve source-file units. UTC "
                              "hour is derived from datetimeUTC; the "
                              "HOUR column in the CSV is ignored.",
        },
    )


def build_daily_ds(df: pd.DataFrame, mw_lookup: dict, source: str) -> xr.Dataset:
    if "HOUR" in df.columns:
        df = df.drop(columns=["HOUR"])
    species_cols = [c for c in df.columns if c not in _ALL_NON_SPECIES]
    log.info("  daily mode: %d species columns", len(species_cols))

    poly_vals = df["POLYID"].values.astype("int64")
    polyids, first_idx = np.unique(poly_vals, return_index=True)
    polyids = polyids.astype("int32")
    n_poly = len(polyids)
    log.info("  %d unique fires (POLYIDs)", n_poly)

    row_poly_idx = np.searchsorted(polyids, poly_vals)
    if n_poly < len(df):
        log.warning("  duplicate POLYID rows -- values will be summed")

    species_arrays = {}
    for sp in species_cols:
        arr = np.zeros(n_poly, dtype="float32")
        np.add.at(arr, row_poly_idx, df[sp].values.astype("float32"))
        species_arrays[sp] = arr

    fire_meta = _fire_meta_arrays(df, first_idx)
    date_int, iso_date = resolve_date(df)
    log.info("  date: %d (%s)", date_int, iso_date)

    data_vars = {}
    for sp in species_cols:
        data_vars[sp] = (("polyid",), species_arrays[sp],
                          _species_attrs(sp, mw_lookup, "day"))
    for col, vals in fire_meta.items():
        data_vars[_FIRE_META_RENAME[col]] = (
            ("polyid",), vals, _FIRE_META_ATTRS[col])
    data_vars["date"] = ((), np.int32(date_int),
                         {"long_name": "date as YYYYMMDD"})

    return xr.Dataset(
        data_vars=data_vars,
        coords={
            "polyid": (("polyid",), polyids, {"long_name": "FINN polygon ID"}),
        },
        attrs={
            "title":          "FINN per-fire daily emissions",
            "source":         source,
            "date_processed": iso_date,
            "Conventions":    "CF-1.8",
            "created":        dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "comment":        "One entry per POLYID (daily average).",
        },
    )


def write_dataset(ds: xr.Dataset, out_path: Path) -> None:
    encoding = {}
    exclude = {"date"} | set(_FIRE_META_RENAME.values())
    for v in ds.data_vars:
        if v not in exclude and ds[v].ndim >= 1:
            encoding[v] = {"zlib": True, "complevel": 4}
    ds.to_netcdf(out_path, format="NETCDF4", encoding=encoding)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    mw_lookup = read_mw_table(args.mw_table)

    ok, fail, skipped = 0, 0, 0
    for path in args.inputs:
        if not path.exists():
            log.error("input not found: %s", path); fail += 1; continue
        out_path = args.out_dir / f"{path.stem}.nc"
        if out_path.exists() and not args.overwrite:
            log.info("skip %s (exists -- use --overwrite)", out_path.name)
            skipped += 1
            continue
        try:
            df = read_csv(path)
            # Mode signal: presence of the HOUR column in the source CSV.
            # HOUR's *values* are unreliable and never trusted — the
            # actual UTC hour is derived from datetimeUTC in
            # build_hourly_ds.  We drop HOUR here after using its
            # presence as the signal, so it can't leak into any output.
            hourly_source = "HOUR" in df.columns
            if "HOUR" in df.columns:
                df = df.drop(columns=["HOUR"])
            hourly = hourly_source and not args.daily
            if hourly:
                ds = build_hourly_ds(df, mw_lookup, source=path.name)
            else:
                ds = build_daily_ds(df, mw_lookup, source=path.name)
            write_dataset(ds, out_path)
            log.info("  wrote %s  (dims: %s)", out_path.name, dict(ds.sizes))
            ok += 1
        except Exception as e:
            log.error("failed on %s: %s", path.name, e, exc_info=args.verbose)
            fail += 1

    log.info("done -- %d ok, %d skipped, %d failed", ok, skipped, fail)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
