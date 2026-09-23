#!/usr/bin/env python3
"""Grid finn_py per-fire emissions onto a regular lat-lon grid or a
SCRIP-described unstructured grid.

Reads the NetCDF written by ``calc_emis_daily.py`` (one ``fire``
dimension, one variable per chemical species in mol/day or kg/day) and
accumulates each species into a 2-D (lat-lon) or 1-D (unstructured)
emission map.  Writes **one output file per species** so the caller can
concatenate days into time series later (the user's typical workflow).

Examples
--------
    # 0.1° × 0.1° regular grid
    python grid_emissions.py FINNv2.9.2_NRTmodvrs_MOZART_20260530.nc \\
        --out-dir ./gridded \\
        --grid-resolution 0.1 0.1

    # CESM-SE ne30 unstructured grid via SCRIP file
    python grid_emissions.py FINNv2.9.2_NRTmodvrs_MOZART_20260530.nc \\
        --out-dir ./gridded_ne30 \\
        --scrip /glade/.../ne30np4_pentagons.091226.nc \\
        --grid-label ne30

    # Just a subset of species; keep totals (mol or kg per cell per day)
    python grid_emissions.py …input.nc \\
        --species CO BC OC PM2.5 \\
        --grid-resolution 0.5 0.5 \\
        --units total

Output filename
---------------
``<prefix>_<species>_<grid-label>_<YYYYMMDD>.nc`` — one file per species,
each holding a single day.  The time coordinate is CF-compliant
("days since YEAR-01-01"), so ``ncrcat`` / ``cdo cat`` along ``time``
just works across days of the same year.

Units
-----
``--units flux`` (default): every species is converted to **molecules
cm⁻² s⁻¹** using the ``molecular_weight_g_per_mol`` attribute carried
on each species variable.  For variables already in ``mol/day``
(gases) the conversion is value × Avogadro / cell_area_cm² / 86400.
For variables in ``kg/day`` (aerosols, where the EF-table uses MW = 1
as a sentinel), mass is first turned into a mole-equivalent
(value × 1000 / MW), then converted the same way — so MW = 1 just
means "1 kg per (fake) mol", giving the per-cm² molecule count that
CAM-Chem fire input expects.

``--units total``: leaves the per-cell daily totals untouched (raw sum):
  * gases    → mol/cell/day
  * aerosols → kg/cell/day

Time coordinate
---------------
Per-day files carry a single ``time`` value as **days since 1970-01-01
00:00:00** on the Gregorian calendar, so daily files concatenate
cleanly across years and decode to the actual calendar date when read
with any CF-aware tool.  Override with ``--time-epoch YYYY-MM-DD`` if
your downstream workflow expects a different epoch.
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree

log = logging.getLogger("grid_emissions")

# Physical constants
R_EARTH_M = 6_371_008.8                    # mean radius of Earth (m)
SECONDS_PER_DAY = 86400.0
AVOGADRO = 6.02214076e23                   # molecules per mole (exact since SI 2019)


def _decimals_for_step(step: float) -> int:
    """Decimal places needed to represent ``step`` cleanly, plus a
    2-decimal safety margin.

    Used to snap lat/lon coordinate arrays to their nominal precision
    so float64 accumulation fuzz (e.g. ``-9.849999999999998`` where
    ``-9.85`` was intended) doesn't leak into the output NetCDF.
    """
    s = f"{step:.15g}"
    if "e" in s or "E" in s:                # scientific notation → default deep
        return 8
    frac_len = len(s.split(".")[-1]) if "." in s else 0
    return frac_len + 2


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="grid_emissions.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("input", type=Path,
                   help="per-fire NetCDF produced by calc_emis_daily.py")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output directory (created if missing)")

    grid_grp = p.add_mutually_exclusive_group(required=True)
    grid_grp.add_argument("--grid-resolution", type=float, nargs=2,
                           metavar=("DLAT", "DLON"),
                           help="regular lat-lon grid resolution in degrees")
    grid_grp.add_argument("--scrip", type=Path,
                           help="path to a SCRIP-format grid description file")

    p.add_argument("--grid-label", default=None,
                   help="grid label embedded in output filenames "
                        "(default: '<DLAT>x<DLON>deg' for lat-lon grids; "
                        "stem of the SCRIP filename for unstructured)")
    p.add_argument("--species", nargs="+", default=None,
                   help="species names to grid; default = every variable in "
                        "the input file whose units are mol/day or kg/day")
    p.add_argument("--units", choices=["flux", "total"], default="flux",
                   help="output units (default: %(default)s) — see header")
    p.add_argument("--filename-prefix", default=None,
                   help="output filename prefix (default: derived from input)")
    p.add_argument("--time-epoch", default=None, metavar="YYYY-MM-DD",
                   help="CF time-units epoch (default: <input-year>-01-01). "
                        "Use a fixed epoch like 2000-01-01 when concatenating "
                        "across multiple years.")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def latlon_to_xyz(lat_deg, lon_deg):
    """(lat, lon) in degrees → unit-sphere Cartesian coords for KDTree."""
    lat = np.radians(np.asarray(lat_deg))
    lon = np.radians(np.asarray(lon_deg))
    return np.column_stack([
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat),
    ])


def detect_species(ds: xr.Dataset) -> list[str]:
    """Return data_vars whose units are mol/day or kg/day."""
    return [
        v for v in ds.data_vars
        if str(ds[v].attrs.get("units", "")).lower() in ("mol/day", "kg/day")
    ]


# ---------------------------------------------------------------------------
# Lat-lon gridding
# ---------------------------------------------------------------------------

def _scale_to_molec_per_cm2_per_s(per_day_value, input_units, mw, cell_area_m2):
    """Convert per-cell per-day to molecules cm⁻² s⁻¹ in-place-ish.

    ``per_day_value`` is whatever units the calc_emis_daily.py NetCDF
    holds (``mol/day`` for gases, ``kg/day`` for aerosols with MW=1).
    ``cell_area_m2`` is either a scalar/1-D array broadcastable against
    the grid shape.

    Both paths use Avogadro; the difference is whether we already have
    moles (gases) or need to derive them from mass via MW (aerosols).
    """
    if input_units == "mol/day":
        # Gas: value is already mol/day, just × N_A for molec/day.  MW
        # is informational only here (already used by calc_emis_daily).
        molec_per_day = per_day_value * AVOGADRO
    elif input_units == "kg/day":
        # Aerosol: value is kg/day; convert to mol/day via MW (g/mol),
        # then × N_A.  MW=1 (the FINN sentinel) yields the CAM-Chem
        # "molecules" convention: 1 kg = 1000 fake-mol = 6.02e26 fake-molec.
        if mw is None or float(mw) <= 0:
            raise ValueError(
                f"can't convert {input_units} to molec/cm²/s without a "
                f"positive molecular_weight_g_per_mol attribute (got {mw!r})"
            )
        molec_per_day = per_day_value * (1000.0 / float(mw)) * AVOGADRO
    else:
        raise ValueError(f"unsupported input units {input_units!r}; "
                          f"expected 'mol/day' or 'kg/day'")

    cell_area_cm2 = cell_area_m2 * 1e4
    return molec_per_day / cell_area_cm2 / SECONDS_PER_DAY


def grid_to_latlon(fire_lat, fire_lon, fire_val, dlat, dlon,
                    *, units: str, input_units: str, mw):
    """Accumulate per-fire values onto a regular lat-lon grid.

    Cell convention: centers offset by half a step from the edges, so
    the grid covers [-90, 90] × [-180, 180) exactly.  For 0.1°: 1800
    latitudes × 3600 longitudes, with centers from -89.95° to 89.95° in
    lat and -179.95° to 179.95° in lon.

    Returns ``(grid[nlat,nlon], lat_centers, lon_centers)`` in float32.
    Output units depend on ``units``:
      * ``"total"`` → raw daily totals (mol/cell/day or kg/cell/day)
      * ``"flux"``  → molecules cm⁻² s⁻¹ (via :func:`_scale_to_molec_per_cm2_per_s`)
    """
    nlat = int(round(180.0 / dlat))
    nlon = int(round(360.0 / dlon))
    if abs(nlat * dlat - 180.0) > 1e-9 or abs(nlon * dlon - 360.0) > 1e-9:
        raise ValueError(
            f"resolution ({dlat}, {dlon}) doesn't divide (180, 360) cleanly")

    lat_edges = -90.0 + np.arange(nlat + 1) * dlat
    lon_edges = -180.0 + np.arange(nlon + 1) * dlon
    # Centers used for the output coordinate arrays.  Snap to the grid-
    # step precision so a call like `-90 + arange(1800)*0.1 + 0.05`
    # doesn't leak float64 fuzz (-9.849999999… instead of -9.85) into
    # the NetCDF.
    _decs = _decimals_for_step(min(dlat, dlon))
    lat_centers = np.round((lat_edges[:-1] + lat_edges[1:]) / 2.0, _decs)
    lon_centers = np.round((lon_edges[:-1] + lon_edges[1:]) / 2.0, _decs)

    grid = np.zeros((nlat, nlon), dtype="float64")
    if len(fire_lat):
        i_lat = np.clip(
            np.searchsorted(lat_edges, fire_lat, side="right") - 1, 0, nlat - 1)
        i_lon = np.clip(
            np.searchsorted(lon_edges, fire_lon, side="right") - 1, 0, nlon - 1)
        if (fire_lat.min() < -90 or fire_lat.max() > 90 or
                fire_lon.min() < -180 or fire_lon.max() > 180):
            log.warning("some fires have out-of-range lat/lon; clipped to grid edges")
        np.add.at(grid, (i_lat, i_lon), fire_val)

    if units == "flux":
        sin_n = np.sin(np.radians(lat_edges[1:]))
        sin_s = np.sin(np.radians(lat_edges[:-1]))
        dlon_rad = np.radians(dlon)
        area_strip = R_EARTH_M ** 2 * dlon_rad * (sin_n - sin_s)  # (nlat,) m²
        grid = _scale_to_molec_per_cm2_per_s(
            grid, input_units, mw, area_strip[:, None]
        )

    return grid.astype("float32"), lat_centers, lon_centers


# ---------------------------------------------------------------------------
# SCRIP unstructured gridding
# ---------------------------------------------------------------------------

def _read_scrip(scrip_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a SCRIP file and return (center_lat_deg, center_lon_deg,
    cell_area_m2).

    Handles common SCRIP variants: lat/lon in radians (the spec) or in
    degrees (commonly produced by CESM tools); ``grid_area`` always
    interpreted as steradians.  Variables can be lower- or upper-case.
    """
    ds = xr.open_dataset(scrip_path, decode_cf=False)
    # Look up variable names case-insensitively
    name_map = {v.lower(): v for v in ds.variables}
    for need in ("grid_center_lat", "grid_center_lon", "grid_area"):
        if need not in name_map:
            raise RuntimeError(
                f"{scrip_path}: missing required variable '{need}'; "
                f"found {list(ds.variables)}")

    clat = ds[name_map["grid_center_lat"]]
    clon = ds[name_map["grid_center_lon"]]
    area = ds[name_map["grid_area"]]

    # Detect degree vs radian for centers
    units_lat = str(clat.attrs.get("units", "radians")).lower()
    if "degree" in units_lat:
        center_lat = clat.values.astype("float64")
        center_lon = clon.values.astype("float64")
    else:                                       # radians (SCRIP default)
        center_lat = np.degrees(clat.values).astype("float64")
        center_lon = np.degrees(clon.values).astype("float64")

    # grid_area is in steradians (radians²) per the SCRIP spec
    cell_area_m2 = area.values.astype("float64") * R_EARTH_M ** 2

    ncol = len(center_lat)
    log.info("SCRIP %s: %d cells, total area = %.3f × Earth surface",
             scrip_path.name, ncol,
             cell_area_m2.sum() / (4 * np.pi * R_EARTH_M ** 2))
    ds.close()
    return center_lat, center_lon, cell_area_m2


def grid_to_scrip(fire_lat, fire_lon, fire_val,
                   center_lat, center_lon, cell_area_m2,
                   *, units: str, input_units: str, mw):
    """Assign each fire to its nearest SCRIP cell (great-circle) and sum.

    Returns ``grid[ncol]`` in float32.  Output units follow the same
    conventions as :func:`grid_to_latlon`.
    """
    ncol = len(center_lat)
    grid = np.zeros(ncol, dtype="float64")
    if len(fire_lat):
        cell_xyz = latlon_to_xyz(center_lat, center_lon)
        fire_xyz = latlon_to_xyz(fire_lat, fire_lon)
        tree = cKDTree(cell_xyz)
        _, i_cell = tree.query(fire_xyz, k=1)
        np.add.at(grid, i_cell, fire_val)

    if units == "flux":
        grid = _scale_to_molec_per_cm2_per_s(grid, input_units, mw, cell_area_m2)

    return grid.astype("float32")


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _output_units_for(input_units: str, mode: str) -> str:
    """Map (per-fire units, target mode) → output-variable units string."""
    if mode == "total":
        # mol/day stays mol/cell/day; kg/day stays kg/cell/day
        return input_units.replace("/day", "/cell/day")
    # flux mode: uniformly molecules cm⁻² s⁻¹ for everything
    return "molecules/cm2/s"


def _time_coord(date_str: str, epoch_str: str
                ) -> tuple[np.ndarray, str]:
    """Return (values, units) for a CF-compliant single-day time axis."""
    target = dt.datetime.strptime(date_str, "%Y-%m-%d").date()
    epoch  = dt.datetime.strptime(epoch_str, "%Y-%m-%d").date()
    days = (target - epoch).days
    return np.array([days], dtype="float64"), f"days since {epoch_str} 00:00:00"


def write_latlon_species(out_path: Path, species: str, grid_2d: np.ndarray,
                          lat_centers: np.ndarray, lon_centers: np.ndarray,
                          *, time_value: float, time_units: str,
                          date_yyyymmdd: int,
                          units_out: str, mw: float | None,
                          extra_attrs: dict) -> None:
    var_attrs = {"units": units_out,
                  "long_name": f"{species} emission"}
    if mw is not None:
        var_attrs["molecular_weight_g_per_mol"] = mw

    ds = xr.Dataset(
        data_vars={
            species: (("time", "lat", "lon"),
                       grid_2d[None, :, :], var_attrs),
            "date":  (("time",),
                       np.array([date_yyyymmdd], dtype="int32"),
                       {"long_name": "calendar date as YYYYMMDD"}),
        },
        coords={
            "time": (("time",), np.array([time_value], dtype="float64"),
                      {"units": time_units, "calendar": "gregorian",
                       "long_name": "time", "standard_name": "time"}),
            "lat":  (("lat",), lat_centers.astype("float64"),
                      {"units": "degrees_north", "standard_name": "latitude"}),
            "lon":  (("lon",), lon_centers.astype("float64"),
                      {"units": "degrees_east", "standard_name": "longitude"}),
        },
        attrs=extra_attrs,
    )
    ds.to_netcdf(out_path, format="NETCDF4",
                  encoding={species: {"zlib": True, "complevel": 4}},
                  unlimited_dims=["time"])


def write_unstructured_species(out_path: Path, species: str, grid_1d: np.ndarray,
                                center_lat: np.ndarray, center_lon: np.ndarray,
                                *, time_value: float, time_units: str,
                                date_yyyymmdd: int,
                                units_out: str, mw: float | None,
                                extra_attrs: dict) -> None:
    var_attrs = {"units": units_out,
                  "long_name": f"{species} emission"}
    if mw is not None:
        var_attrs["molecular_weight_g_per_mol"] = mw

    ds = xr.Dataset(
        data_vars={
            species: (("time", "ncol"),
                       grid_1d[None, :], var_attrs),
            "date":  (("time",),
                       np.array([date_yyyymmdd], dtype="int32"),
                       {"long_name": "calendar date as YYYYMMDD"}),
        },
        coords={
            "time": (("time",), np.array([time_value], dtype="float64"),
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
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )

    if not args.input.exists():
        log.error("input file not found: %s", args.input)
        return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Load per-fire NetCDF
    log.info("reading per-fire emissions from %s", args.input)
    fires = xr.open_dataset(args.input)
    n_fire = fires.sizes.get("fire", 0)
    log.info("  %d fire rows", n_fire)

    # Resolve the date of these fires (derived from the per-fire `date`
    # variable, which holds YYYYMMDD ints).  All fires in a calc_emis
    # output should share one date, but tolerate mixed dates by using the
    # earliest as the time stamp on the gridded output.
    if "date" in fires.variables and n_fire:
        date_ints = fires.date.values.astype("int64")
        if len(np.unique(date_ints)) > 1:
            log.warning("input contains multiple dates %s; using earliest",
                        sorted(set(date_ints.tolist()))[:5])
        d_int = int(date_ints.min())
        year, month, day = d_int // 10000, (d_int // 100) % 100, d_int % 100
        date_str = f"{year:04d}-{month:02d}-{day:02d}"
        ymd = f"{year:04d}{month:02d}{day:02d}"
    else:
        # Fall back to global attribute or today
        date_str = fires.attrs.get("date_processed", dt.date.today().isoformat())
        year = int(date_str[:4])
        ymd = date_str.replace("-", "")
        d_int = int(ymd)
    epoch_str = args.time_epoch or "1970-01-01"
    time_value, time_units = _time_coord(date_str, epoch_str)
    log.info("  time = %s  (= %.1f %s, gregorian calendar)",
             date_str, time_value[0], time_units)

    # Which species?
    species = args.species or detect_species(fires)
    if not species:
        log.error("no species found in %s — pass --species or check input file",
                  args.input)
        return 1
    log.info("species: %s", species)

    # Grid label
    if args.grid_label:
        glabel = args.grid_label
    elif args.scrip:
        glabel = args.scrip.stem
    else:
        dlat, dlon = args.grid_resolution
        # Format e.g. 0.1x0.1deg, 0.5x0.625deg
        glabel = (f"{dlat:g}x{dlon:g}deg"
                   .replace(" ", ""))

    # Output filename prefix
    if args.filename_prefix:
        prefix = args.filename_prefix
    else:
        # Use everything before the species label in the input filename,
        # e.g. FINNv2.9.2_NRTmodvrs_MOZART_20260530.nc
        # → FINNv2.9.2_NRTmodvrs_MOZART
        stem = args.input.stem
        # Strip trailing date if present
        bits = stem.split("_")
        if bits and len(bits[-1]) == 8 and bits[-1].isdigit():
            bits = bits[:-1]
        prefix = "_".join(bits) or stem

    # If SCRIP mode, pre-load the grid once
    grid_attrs_extra = {}
    if args.scrip:
        center_lat, center_lon, cell_area_m2 = _read_scrip(args.scrip)
        grid_attrs_extra["scrip_file"] = str(args.scrip)
        grid_attrs_extra["grid_label"] = glabel
    else:
        dlat, dlon = args.grid_resolution
        grid_attrs_extra["grid_resolution_deg"] = f"{dlat} x {dlon}"
        grid_attrs_extra["grid_label"] = glabel

    fire_lat = fires.lat.values.astype("float64") if n_fire else np.array([])
    fire_lon = fires.lon.values.astype("float64") if n_fire else np.array([])

    # Per-species loop
    for sp in species:
        if sp not in fires.variables:
            log.warning("species %s not in input; skipping", sp)
            continue
        in_units = str(fires[sp].attrs.get("units", "mol/day"))
        out_units = _output_units_for(in_units, args.units)
        mw = fires[sp].attrs.get("molecular_weight_g_per_mol")
        fire_val = (fires[sp].values.astype("float64")
                    if n_fire else np.array([]))

        out_name = f"{prefix}_{sp}_{glabel}_{ymd}.nc"
        out_path = args.out_dir / out_name

        common_attrs = {
            "title":         f"FINN gridded emissions of {sp}",
            "source":        f"gridded from {args.input.name}",
            "date_processed": date_str,
            "units_mode":     args.units,
            "Conventions":    "CF-1.8",
            **{k: v for k, v in fires.attrs.items()
                if k in ("finn_version", "sim_id", "species_set", "title")},
            **grid_attrs_extra,
        }

        if args.scrip:
            grid_1d = grid_to_scrip(
                fire_lat, fire_lon, fire_val,
                center_lat, center_lon, cell_area_m2,
                units=args.units, input_units=in_units, mw=mw,
            )
            write_unstructured_species(
                out_path, sp, grid_1d, center_lat, center_lon,
                time_value=float(time_value[0]), time_units=time_units,
                date_yyyymmdd=d_int,
                units_out=out_units, mw=mw, extra_attrs=common_attrs,
            )
            nz = int(np.count_nonzero(grid_1d))
            tot = float(grid_1d.sum())
        else:
            dlat, dlon = args.grid_resolution
            grid_2d, lat_c, lon_c = grid_to_latlon(
                fire_lat, fire_lon, fire_val, dlat, dlon,
                units=args.units, input_units=in_units, mw=mw,
            )
            write_latlon_species(
                out_path, sp, grid_2d, lat_c, lon_c,
                time_value=float(time_value[0]), time_units=time_units,
                date_yyyymmdd=d_int,
                units_out=out_units, mw=mw, extra_attrs=common_attrs,
            )
            nz = int(np.count_nonzero(grid_2d))
            tot = float(grid_2d.sum())

        log.info("  %s: %d non-zero cells, total %.3e %s → %s",
                  sp, nz, tot, out_units, out_path.name)

    fires.close()
    log.info("done — %d species written to %s", len(species), args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
