#!/usr/bin/env python3
"""Grid FINN text-format emissions onto a regular lat-lon or SCRIP grid.

Legacy FINN pre-processor output is a per-day CSV with columns:

    DAY, FIREID, POLYID, GENVEG, LATI, LONGI, AREA, BMASS, FRP,
    <species1>, <species2>, ..., HOUR, datetimeUTC, datetimeLT,
    Country, TrendCountry

Each row is one (fire × hour) sample.  This script accumulates those
into gridded fields, converts to molecules cm⁻² s⁻¹ using a caller-
supplied molecular-weight table, and writes **one NetCDF per species**
with dimensions

    (time, lat, lon)               for --grid-resolution
    (time, ncol)                   for --scrip

where ``time`` has 24 entries per input file (one per hour of the day).
Companion variables ``date`` (YYYYMMDD, int) and ``datesec`` (seconds
of day, int) are repeated across every timestep, so multi-day
concatenations along ``time`` produce a proper CAM/CESM-style emission
input file.

Ignored columns
---------------
``Country`` and ``TrendCountry`` are always dropped.  Metadata columns
(``DAY``, ``FIREID``, ``POLYID``, ``GENVEG``, ``LATI``, ``LONGI``,
``AREA``, ``BMASS``, ``FRP``, ``HOUR``, ``datetimeUTC``, ``datetimeLT``)
are never gridded.

MW table
--------
The ``--mw-table`` CSV must have columns ``species,type,mw[,notes]``.
See ``scripts/finn_text_species_mw.csv`` for a starter file.  Species
that are in the input CSV but *not* in the MW table are skipped with
a warning.

Examples
--------
    # Global 0.1° × 0.1° grid, every species that has a MW entry
    python grid_txt_emissions.py \\
        FINNv2.9nrt_2024080_hourly.txt \\
        --mw-table scripts/finn_text_species_mw.csv \\
        --grid-resolution 0.1 0.1 \\
        --out-dir ./gridded

    # MPAS unstructured grid via SCRIP file, subset of species
    python grid_txt_emissions.py \\
        FINNv2.9nrt_2024080_hourly.txt \\
        --mw-table scripts/finn_text_species_mw.csv \\
        --scrip /glade/.../scrip_mxc.nc \\
        --grid-label mxc \\
        --species CO CO2 BC OC PM2.5 \\
        --out-dir ./gridded_mxc
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

# Re-use the geometry / SCRIP / binning primitives from grid_emissions.py.
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))
from grid_emissions import (
    R_EARTH_M, AVOGADRO,
    latlon_to_xyz, _read_scrip, _decimals_for_step,
)

log = logging.getLogger("grid_txt_emissions")

# Columns that must never be gridded even if they contain floats.
_METADATA_COLS = frozenset({
    "DAY", "FIREID", "POLYID", "GENVEG", "LATI", "LONGI",
    "AREA", "BMASS", "FRP", "HOUR",
    "datetimeUTC", "datetimeLT",
    "Country", "TrendCountry",
})


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="grid_txt_emissions.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("input", type=Path,
                   help="legacy FINN text emissions CSV (one file = one day)")
    p.add_argument("--mw-table", type=Path, required=True,
                   help="CSV with columns species,type,mw[,notes]. "
                        "type ∈ {gas, aerosol, number}.")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output directory (created if missing)")

    grid_grp = p.add_mutually_exclusive_group(required=True)
    grid_grp.add_argument("--grid-resolution", type=float, nargs=2,
                           metavar=("DLAT", "DLON"),
                           help="regular lat-lon grid resolution in degrees")
    grid_grp.add_argument("--scrip", type=Path,
                           help="path to a SCRIP-format grid description file")

    p.add_argument("--species", nargs="+", default=None,
                   help="species to grid (default: every non-metadata column "
                        "that appears in the MW table)")
    p.add_argument("--grid-label", default=None,
                   help="grid label embedded in output filenames "
                        "(default: '<DLAT>x<DLON>deg' or SCRIP stem)")
    p.add_argument("--filename-prefix", default=None,
                   help="output filename prefix (default: derived from input)")
    p.add_argument("--input-time-units", choices=["hour", "day"], default="hour",
                   help="whether each row's value represents per-hour "
                        "(default) or per-day emissions.  Only affects the "
                        "denominator when converting to molec/cm2/s.")
    p.add_argument("--units", choices=["flux", "total"], default="flux",
                   help="'flux' = molecules/cm2/s (default); "
                        "'total' = leave input rate untouched, sum per cell")
    p.add_argument("--time-epoch", default="1970-01-01",
                   help="CF time-units epoch (default: %(default)s)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# MW-table reader
# ---------------------------------------------------------------------------

def read_mw_table(path: Path) -> dict[str, dict]:
    """Parse the MW-table CSV; return {species: {type: str, mw: float|nan}}."""
    df = pd.read_csv(path, comment="#", skip_blank_lines=True)
    df.columns = [c.strip() for c in df.columns]
    for req in ("species", "type"):
        if req not in df.columns:
            raise RuntimeError(f"{path}: MW table missing required column '{req}'")
    if "mw" not in df.columns:
        df["mw"] = np.nan
    df = df.dropna(subset=["species", "type"], how="any")
    df["species"] = df["species"].astype(str).str.strip()
    df["type"]    = df["type"].astype(str).str.strip().str.lower()
    df = df[df["species"].astype(bool)]                     # drop empty rows

    valid = {"gas", "aerosol", "number"}
    bad_type = ~df["type"].isin(valid)
    if bad_type.any():
        log.warning("MW table: %d rows have invalid type; ignoring: %s",
                    int(bad_type.sum()), df.loc[bad_type, "species"].tolist())
        df = df[~bad_type]

    table = {}
    for _, r in df.iterrows():
        mw_val = float(r["mw"]) if pd.notna(r["mw"]) else float("nan")
        table[r["species"]] = {"type": r["type"], "mw": mw_val}
    log.info("MW table: %d species (gas=%d, aerosol=%d, number=%d)",
              len(table),
              sum(1 for v in table.values() if v["type"] == "gas"),
              sum(1 for v in table.values() if v["type"] == "aerosol"),
              sum(1 for v in table.values() if v["type"] == "number"))
    return table


# ---------------------------------------------------------------------------
# Text-CSV reader
# ---------------------------------------------------------------------------

def read_text_emissions(path: Path) -> pd.DataFrame:
    """Read a FINN legacy text-CSV.  Drops Country/TrendCountry and the
    unreliable HOUR column, sorts rows by datetimeUTC ascending."""
    log.info("reading %s", path)
    df = pd.read_csv(path)
    # Drop columns we're told to ignore
    for col in ("Country", "TrendCountry"):
        if col in df.columns:
            df = df.drop(columns=col)
    # HOUR values in some FINN outputs are unreliable — derive the UTC
    # hour from datetimeUTC instead, and drop HOUR outright so it
    # can't be used downstream by mistake or leaked into outputs.
    if "HOUR" in df.columns:
        df = df.drop(columns=["HOUR"])
    for req in ("datetimeUTC", "LATI", "LONGI", "DAY"):
        if req not in df.columns:
            raise RuntimeError(f"{path}: missing required column {req!r}")
    # Sort by datetimeUTC (ISO-formatted strings sort lexicographically
    # in chronological order for same-length values).
    df = df.sort_values("datetimeUTC", kind="mergesort").reset_index(drop=True)
    utc_hours = pd.to_datetime(df["datetimeUTC"]).dt.hour
    log.info("  %d rows, %d columns; UTC-hour range %d..%d "
             "(derived from datetimeUTC; HOUR column ignored)",
             len(df), len(df.columns), int(utc_hours.min()), int(utc_hours.max()))
    return df


def resolve_date(df: pd.DataFrame) -> tuple[int, str, int]:
    """Return (YYYYMMDD int, YYYY-MM-DD string, year int) for the file."""
    if "datetimeUTC" in df.columns and df["datetimeUTC"].notna().any():
        first = str(df["datetimeUTC"].dropna().iloc[0])[:10]        # YYYY-MM-DD
        y, m, d = int(first[:4]), int(first[5:7]), int(first[8:10])
    else:                                                             # fall back to DAY
        year = dt.date.today().year
        day = int(df["DAY"].iloc[0])
        d0 = dt.date(year, 1, 1) + dt.timedelta(days=day - 1)
        y, m, d = d0.year, d0.month, d0.day
        first = d0.isoformat()
    return y * 10000 + m * 100 + d, first, y


# ---------------------------------------------------------------------------
# Gridding — one hour at a time, one species at a time.
# ---------------------------------------------------------------------------

def _cell_areas_latlon(dlat: float, dlon: float
                        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (lat_edges, lat_centers, lon_centers, cell_area_m2[nlat])."""
    nlat = int(round(180.0 / dlat))
    nlon = int(round(360.0 / dlon))
    if abs(nlat * dlat - 180.0) > 1e-9 or abs(nlon * dlon - 360.0) > 1e-9:
        raise ValueError(
            f"resolution ({dlat}, {dlon}) doesn't divide (180, 360) cleanly")
    lat_edges = -90.0 + np.arange(nlat + 1) * dlat
    lon_edges = -180.0 + np.arange(nlon + 1) * dlon
    # Snap the output-coordinate arrays to the grid-step precision so
    # float64 accumulation fuzz doesn't leak into the NetCDF.  Edges
    # keep their raw values because they drive the binning searchsorted.
    _decs = _decimals_for_step(min(dlat, dlon))
    lat_centers = np.round((lat_edges[:-1] + lat_edges[1:]) / 2.0, _decs)
    lon_centers = np.round((lon_edges[:-1] + lon_edges[1:]) / 2.0, _decs)
    sin_n = np.sin(np.radians(lat_edges[1:]))
    sin_s = np.sin(np.radians(lat_edges[:-1]))
    area_strip = R_EARTH_M ** 2 * np.radians(dlon) * (sin_n - sin_s)   # (nlat,)
    return lat_edges, lat_centers, lon_edges, lon_centers, area_strip


def _bin_latlon(lat: np.ndarray, lon: np.ndarray, val: np.ndarray,
                lat_edges: np.ndarray, lon_edges: np.ndarray
                ) -> np.ndarray:
    """Sum ``val`` into (nlat, nlon) grid at the (lat, lon) bins."""
    nlat = len(lat_edges) - 1
    nlon = len(lon_edges) - 1
    grid = np.zeros((nlat, nlon), dtype="float64")
    if len(lat) == 0:
        return grid
    i_lat = np.clip(np.searchsorted(lat_edges, lat, side="right") - 1, 0, nlat - 1)
    i_lon = np.clip(np.searchsorted(lon_edges, lon, side="right") - 1, 0, nlon - 1)
    np.add.at(grid, (i_lat, i_lon), val)
    return grid


def _bin_scrip(lat: np.ndarray, lon: np.ndarray, val: np.ndarray,
                cell_xyz_tree, ncells: int) -> np.ndarray:
    """Sum ``val`` into ncells-sized grid via nearest-neighbor to cell centers."""
    grid = np.zeros(ncells, dtype="float64")
    if len(lat) == 0:
        return grid
    fire_xyz = latlon_to_xyz(lat, lon)
    _, i_cell = cell_xyz_tree.query(fire_xyz, k=1)
    np.add.at(grid, i_cell, val)
    return grid


def _scale_to_flux(per_time: np.ndarray, sp_info: dict,
                    cell_area_m2, seconds: float) -> np.ndarray:
    """Convert per-cell/per-time totals to molec/cm²/s (or num/cm²/s).

    ``cell_area_m2`` broadcasts with ``per_time`` (scalar / 1-D / 2-D).
    """
    typ = sp_info["type"]
    mw = sp_info["mw"]
    if typ == "gas":
        # mol/time → molec via Avogadro
        result = per_time * AVOGADRO
    elif typ == "aerosol":
        if not np.isfinite(mw) or mw <= 0:
            raise ValueError(f"aerosol species requires positive MW, got {mw}")
        # kg/time → molec via (1000 g/kg) / (MW g/mol) × N_A
        result = per_time * (1000.0 / mw) * AVOGADRO
    elif typ == "number":
        result = per_time                                    # already number/time
    else:
        raise ValueError(f"unknown type {typ!r}")
    cell_area_cm2 = cell_area_m2 * 1e4
    return result / cell_area_cm2 / seconds


def output_units_for(sp_info: dict, mode: str) -> str:
    typ = sp_info["type"]
    if mode == "total":
        return "kg/cell/hour" if typ == "aerosol" else (
                "number/cell/hour" if typ == "number" else "mol/cell/hour")
    return "molecules/cm2/s" if typ != "number" else "number/cm2/s"


# ---------------------------------------------------------------------------
# NetCDF writers
# ---------------------------------------------------------------------------

def _time_values(yyyymmdd: int, epoch_str: str
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Return per-hour arrays for the CAM/CESM emission-input convention.

    All three arrays have length 24 (one entry per hour of the file's
    day).  Time is fractional days since the CF epoch; ``date`` repeats
    the same YYYYMMDD int on every hour; ``datesec`` is 0, 3600,
    7200, …, 82800.
    """
    y, m, d = yyyymmdd // 10000, (yyyymmdd // 100) % 100, yyyymmdd % 100
    day_zero = dt.date(y, m, d)
    epoch = dt.datetime.strptime(epoch_str, "%Y-%m-%d").date()
    day_offset = (day_zero - epoch).days
    hours = np.arange(24)
    time_vals    = (day_offset + hours / 24.0).astype("float64")     # (24,)
    date_vals    = np.full(24, yyyymmdd, dtype="int32")              # (24,)
    datesec_vals = (hours * 3600).astype("int32")                    # (24,)
    return time_vals, date_vals, datesec_vals, \
           f"days since {epoch_str} 00:00:00"


def write_species_latlon(out_path: Path, species: str, sp_info: dict,
                          grid_hlL: np.ndarray, lat: np.ndarray,
                          lon: np.ndarray, *,
                          time_vals: np.ndarray, date_vals: np.ndarray,
                          datesec_vals: np.ndarray, time_units: str,
                          out_units: str, extra_attrs: dict) -> None:
    """Write a single-species NetCDF, dims (time, lat, lon).

    ``grid_hlL`` shape (24, nlat, nlon).
    """
    var_attrs = {"units": out_units,
                  "long_name": f"{species} emission",
                  "species_type": sp_info["type"]}
    if np.isfinite(sp_info["mw"]):
        var_attrs["molecular_weight_g_per_mol"] = float(sp_info["mw"])

    ds = xr.Dataset(
        data_vars={
            species:   (("time", "lat", "lon"),
                         grid_hlL.astype("float32"), var_attrs),
            "date":    (("time",), date_vals,
                         {"long_name": "current date (YYYYMMDD)"}),
            "datesec": (("time",), datesec_vals,
                         {"long_name": "current seconds of current date",
                          "units": "s"}),
        },
        coords={
            "time": (("time",), time_vals,
                     {"units": time_units, "calendar": "gregorian",
                      "long_name": "time", "standard_name": "time"}),
            "lat":  (("lat",), lat.astype("float64"),
                     {"units": "degrees_north", "standard_name": "latitude"}),
            "lon":  (("lon",), lon.astype("float64"),
                     {"units": "degrees_east", "standard_name": "longitude"}),
        },
        attrs=extra_attrs,
    )
    ds.to_netcdf(out_path, format="NETCDF4",
                  encoding={species: {"zlib": True, "complevel": 4}},
                  unlimited_dims=["time"])


def write_species_scrip(out_path: Path, species: str, sp_info: dict,
                         grid_hc: np.ndarray, center_lat: np.ndarray,
                         center_lon: np.ndarray, *,
                         time_vals: np.ndarray, date_vals: np.ndarray,
                         datesec_vals: np.ndarray, time_units: str,
                         out_units: str, extra_attrs: dict) -> None:
    """Write a single-species NetCDF, dims (time, ncol).

    ``grid_hc`` shape (24, ncol).
    """
    var_attrs = {"units": out_units,
                  "long_name": f"{species} emission",
                  "species_type": sp_info["type"]}
    if np.isfinite(sp_info["mw"]):
        var_attrs["molecular_weight_g_per_mol"] = float(sp_info["mw"])

    ds = xr.Dataset(
        data_vars={
            species:   (("time", "ncol"),
                         grid_hc.astype("float32"), var_attrs),
            "date":    (("time",), date_vals,
                         {"long_name": "current date (YYYYMMDD)"}),
            "datesec": (("time",), datesec_vals,
                         {"long_name": "current seconds of current date",
                          "units": "s"}),
        },
        coords={
            "time": (("time",), time_vals,
                     {"units": time_units, "calendar": "gregorian",
                      "long_name": "time", "standard_name": "time"}),
            "lat":  (("ncol",), center_lat.astype("float64"),
                     {"units": "degrees_north", "standard_name": "latitude"}),
            "lon":  (("ncol",), center_lon.astype("float64"),
                     {"units": "degrees_east", "standard_name": "longitude"}),
        },
        attrs=extra_attrs,
    )
    ds.to_netcdf(out_path, format="NETCDF4",
                  encoding={species: {"zlib": True, "complevel": 4}},
                  unlimited_dims=["time"])


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )

    if not args.input.exists():
        log.error("input file not found: %s", args.input); return 1
    if not args.mw_table.exists():
        log.error("MW-table file not found: %s", args.mw_table); return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Load MW table + input
    mw_table = read_mw_table(args.mw_table)
    df = read_text_emissions(args.input)
    date_int, date_str, year = resolve_date(df)
    log.info("processing date %d (%s)", date_int, date_str)

    # Which columns are candidate species?
    candidates = [c for c in df.columns if c not in _METADATA_COLS]
    if args.species:
        species_list = [s for s in args.species if s in candidates]
        missing = set(args.species) - set(candidates)
        if missing:
            log.warning("--species not present in input: %s", sorted(missing))
    else:
        species_list = candidates

    # Filter to those we have MW info for
    keep, skip = [], []
    for sp in species_list:
        if sp not in mw_table:
            skip.append((sp, "not in MW table"))
            continue
        info = mw_table[sp]
        if info["type"] == "aerosol" and not (np.isfinite(info["mw"]) and info["mw"] > 0):
            skip.append((sp, "aerosol without positive MW")); continue
        if info["type"] == "gas" and not np.isfinite(info["mw"]):
            skip.append((sp, "gas without MW")); continue
        keep.append(sp)
    log.info("gridding %d species: %s", len(keep), keep)
    if skip:
        for sp, reason in skip:
            log.warning("skip %s: %s", sp, reason)

    # Set up the grid geometry (once)
    if args.scrip:
        center_lat, center_lon, cell_area_m2 = _read_scrip(args.scrip)
        ncells = len(center_lat)
        from scipy.spatial import cKDTree
        cell_xyz_tree = cKDTree(latlon_to_xyz(center_lat, center_lon))
        grid_kind = "scrip"
        glabel = args.grid_label or args.scrip.stem
    else:
        dlat, dlon = args.grid_resolution
        lat_edges, lat_centers, lon_edges, lon_centers, area_strip_m2 = \
            _cell_areas_latlon(dlat, dlon)
        grid_kind = "latlon"
        glabel = args.grid_label or f"{dlat:g}x{dlon:g}deg"

    # Time / date / datesec — 1-D arrays of length 24
    time_vals, date_vals, datesec_vals, time_units = _time_values(
        date_int, args.time_epoch)

    # Pre-extract raw arrays (avoid repeated pandas indexing)
    hour_arr = pd.to_datetime(df["datetimeUTC"]).dt.hour.to_numpy(dtype=np.int32)
    lat_arr  = df["LATI"].to_numpy(dtype=np.float64)
    lon_arr  = df["LONGI"].to_numpy(dtype=np.float64)
    # Precompute row indices per hour
    per_hour_idx = [np.flatnonzero(hour_arr == h) for h in range(24)]

    seconds = 3600.0 if args.input_time_units == "hour" else 86400.0

    # Filename prefix
    prefix = args.filename_prefix or args.input.stem

    common_attrs = {
        "title":           "FINN text-format emissions gridded to lat-lon / SCRIP",
        "source":          f"gridded from {args.input.name}",
        "input_file":      str(args.input),
        "mw_table":        str(args.mw_table),
        "date_processed":  date_str,
        "input_time_units": args.input_time_units,
        "units_mode":      args.units,
        "grid_label":      glabel,
        "Conventions":     "CF-1.8",
        "created":         dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
    }

    # Iterate species
    for sp in keep:
        info = mw_table[sp]
        val_arr = df[sp].to_numpy(dtype=np.float64)
        if grid_kind == "latlon":
            hourly = np.zeros((24, len(lat_centers), len(lon_centers)), dtype="float64")
            for h in range(24):
                idx = per_hour_idx[h]
                if idx.size == 0: continue
                hourly[h] = _bin_latlon(
                    lat_arr[idx], lon_arr[idx], val_arr[idx],
                    lat_edges, lon_edges,
                )
            if args.units == "flux":
                hourly = _scale_to_flux(hourly, info,
                                         area_strip_m2[None, :, None], seconds)
            out_name = f"{prefix}_{sp}_{glabel}_{date_int}.nc"
            out_path = args.out_dir / out_name
            write_species_latlon(
                out_path, sp, info, hourly, lat_centers, lon_centers,
                time_vals=time_vals, date_vals=date_vals,
                datesec_vals=datesec_vals, time_units=time_units,
                out_units=output_units_for(info, args.units),
                extra_attrs={**common_attrs, "grid_resolution_deg": f"{dlat} x {dlon}"},
            )
        else:                                             # scrip
            hourly = np.zeros((24, ncells), dtype="float64")
            for h in range(24):
                idx = per_hour_idx[h]
                if idx.size == 0: continue
                hourly[h] = _bin_scrip(
                    lat_arr[idx], lon_arr[idx], val_arr[idx],
                    cell_xyz_tree, ncells,
                )
            if args.units == "flux":
                hourly = _scale_to_flux(hourly, info,
                                         cell_area_m2[None, :], seconds)
            out_name = f"{prefix}_{sp}_{glabel}_{date_int}.nc"
            out_path = args.out_dir / out_name
            write_species_scrip(
                out_path, sp, info, hourly, center_lat, center_lon,
                time_vals=time_vals, date_vals=date_vals,
                datesec_vals=datesec_vals, time_units=time_units,
                out_units=output_units_for(info, args.units),
                extra_attrs={**common_attrs, "scrip_file": str(args.scrip)},
            )

        nz = int(np.count_nonzero(hourly))
        tot = float(hourly.sum())
        log.info("  %s (%s): %d non-zero cell-hours, total %.3e %s → %s",
                  sp, info["type"], nz, tot,
                  output_units_for(info, args.units), out_path.name)

    log.info("done — %d species written to %s", len(keep), args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
