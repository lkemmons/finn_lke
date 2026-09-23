#!/usr/bin/env python3
"""Compute regional emissions totals from gridded FINN NetCDF and plot.

Reads gridded emissions produced by ``grid_emissions.py`` or
``grid_txt_emissions.py`` (or any similarly structured file) and
integrates each species over one or more geographic regions per
timestep.  Handles:

* Regular lat-lon grids (dims ``lat``, ``lon``) — cell areas computed
  from the lat coordinate spacing.
* Unstructured / SCRIP-described grids (dim ``ncol``) — cell areas read
  from a SCRIP file passed with ``--scrip``.
* Daily emissions (only ``date`` YYYYMMDD present) — each timestep is
  a per-day average, integrated over 86400 s.
* Hourly emissions (``date`` + ``datesec`` in seconds-of-day) — each
  timestep integrated over 3600 s.

Output
------
For each input NetCDF file:
* ``<input-stem>_totals.csv`` with columns ``time`` + one per region
  (values in kg where MW is available, else raw molecule/number counts).
* ``<input-stem>_totals.png`` — time series with one line per region.

Examples
--------
    # Global lat-lon file, default regions
    python compute_emis_totals.py \\
        FINN_v29_CO_0.1x0.1deg_20240320.nc \\
        --out-dir ./totals

    # SCRIP file, with a custom region added
    python compute_emis_totals.py \\
        FINN_v29_CO_mxc_20240320.nc \\
        --scrip /glade/.../scrip_mxc.nc \\
        --region conus -125 -66 24 50 \\
        --out-dir ./totals

    # Multiple species files at once (one plot each)
    python compute_emis_totals.py \\
        FINN_v29_CO_*.nc FINN_v29_BC_*.nc \\
        --out-dir ./totals --y-log
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
import matplotlib as mpl
import matplotlib.pyplot as plt

log = logging.getLogger("compute_emis_totals")

R_EARTH_M = 6_371_008.8
AVOGADRO  = 6.02214076e23

# Default regions covering the globe with some standard slices.
DEFAULT_REGIONS = {
    "global":           (-180.0, 180.0, -90.0,   90.0),
    "nh":               (-180.0, 180.0,   0.0,   90.0),
    "sh":               (-180.0, 180.0, -90.0,    0.0),
    "tropics":          (-180.0, 180.0, -23.5,  23.5),
    "nh_extratropics":  (-180.0, 180.0,  23.5,  90.0),
    "sh_extratropics":  (-180.0, 180.0, -90.0, -23.5),
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="compute_emis_totals.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("inputs", type=Path, nargs="+",
                   help="one or more gridded emissions NetCDFs")
    p.add_argument("--scrip", type=Path, default=None,
                   help="SCRIP file (required for unstructured/ncol grids)")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output directory")
    p.add_argument("--region", action="append", nargs=5, default=None,
                   metavar=("NAME", "W", "E", "S", "N"),
                   help="custom region (repeatable): name west east south "
                        "north (degrees).  If given, the default regions "
                        "are dropped unless --include-defaults is set.")
    p.add_argument("--include-defaults", action="store_true",
                   help="always include the default set of regions "
                        "(global/hemispheres/tropics) even when --region "
                        "is used")
    p.add_argument("--no-plot", action="store_true",
                   help="skip PNG output")
    p.add_argument("--y-log", action="store_true",
                   help="use log-scale y-axis on the plot")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# File I/O helpers
# ---------------------------------------------------------------------------

def _identify_species_var(ds: xr.Dataset) -> str:
    """Return the name of the flux variable (excluding known metadata)."""
    excluded = {"date", "datesec", "time", "hour",
                 "lat", "lon", "latitude", "longitude",
                 "area", "grid_area"}
    candidates = [
        v for v in ds.data_vars
        if v not in excluded and (ds[v].ndim >= 2 or "time" in ds[v].dims)
    ]
    if not candidates:
        # Fall back to any non-excluded data var
        candidates = [v for v in ds.data_vars if v not in excluded]
    if not candidates:
        raise RuntimeError(f"no species variable found; data_vars={list(ds.data_vars)}")
    if len(candidates) > 1:
        # Prefer one whose units look like a flux / mass rate.
        for c in candidates:
            u = str(ds[c].attrs.get("units", "")).lower()
            if any(t in u for t in ("molec", "molecules", "kg", "mol", "number/cm2")):
                return c
    return candidates[0]


def _timestamps_from_date_datesec(dates: np.ndarray, datesec: np.ndarray
                                    ) -> np.ndarray:
    """Combine YYYYMMDD int + seconds-of-day int → np.datetime64[s] array."""
    out = np.empty(len(dates), dtype="datetime64[s]")
    for i, (d, s) in enumerate(zip(dates.astype(int), datesec.astype(int))):
        y = d // 10000
        m = (d // 100) % 100
        day = d % 100
        base = np.datetime64(f"{y:04d}-{m:02d}-{day:02d}", "s")
        out[i] = base + np.timedelta64(int(s), "s")
    return out


def _timestamps_from_date(dates: np.ndarray) -> np.ndarray:
    """YYYYMMDD int → np.datetime64[D] array (midnight of each day)."""
    out = np.empty(len(dates), dtype="datetime64[s]")
    for i, d in enumerate(dates.astype(int)):
        y = d // 10000
        m = (d // 100) % 100
        day = d % 100
        out[i] = np.datetime64(f"{y:04d}-{m:02d}-{day:02d}", "s")
    return out


def infer_cadence(ds: xr.Dataset, timestamps: np.ndarray
                    ) -> tuple[int, bool]:
    """Return ``(duration_seconds_per_timestep, is_hourly)``.

    Detection precedence:
      1. If ``datesec`` variable is present, use consecutive datesec
         differences (fall back to 3600 s).
      2. Else if there are ≥2 timesteps, use the difference between the
         first two timestamps.
      3. Else assume daily (86400 s).
    """
    if "datesec" in ds.variables:
        ds_arr = ds["datesec"].values
        if len(ds_arr) >= 2:
            diffs = np.diff(ds_arr.astype(int))
            positive = diffs[diffs > 0]
            if positive.size > 0:
                return int(positive.min()), True
        return 3600, True

    if len(timestamps) >= 2:
        delta_s = int((timestamps[1] - timestamps[0]) / np.timedelta64(1, "s"))
        return delta_s, (delta_s < 86400)

    return 86400, False


def get_timestamps(ds: xr.Dataset) -> np.ndarray:
    """Return an array of np.datetime64[s], one per timestep."""
    if "date" in ds.variables and "datesec" in ds.variables:
        return _timestamps_from_date_datesec(ds["date"].values,
                                              ds["datesec"].values)
    if "date" in ds.variables:
        return _timestamps_from_date(ds["date"].values)
    if "time" in ds.coords:
        # xarray gives us datetime64 if time was CF-decoded, else raw floats
        tv = ds["time"].values
        if np.issubdtype(tv.dtype, np.datetime64):
            return tv.astype("datetime64[s]")
        raise RuntimeError("time coordinate is not CF-decoded and no "
                            "date/datesec variables present")
    raise RuntimeError("no time information found in file")


def detect_grid_type(ds: xr.Dataset) -> str:
    if "ncol" in ds.dims:
        return "unstructured"
    if "lat" in ds.dims and "lon" in ds.dims:
        return "latlon"
    raise RuntimeError(f"cannot detect grid type from dims {list(ds.dims)}")


# ---------------------------------------------------------------------------
# Cell-area computations
# ---------------------------------------------------------------------------

def compute_latlon_cell_areas(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Cell areas (m²) for a regular lat-lon grid, shape (nlat, nlon).

    Assumes uniform spacing.  If cells aren't uniform, use SCRIP instead.
    """
    nlat, nlon = len(lat), len(lon)
    dlat = float(np.abs(lat[1] - lat[0])) if nlat > 1 else 180.0 / nlat
    dlon = float(np.abs(lon[1] - lon[0])) if nlon > 1 else 360.0 / nlon
    lat_sorted = np.sort(lat)
    lat_edges = np.r_[lat_sorted - dlat / 2.0, lat_sorted[-1] + dlat / 2.0]
    sin_n = np.sin(np.radians(lat_edges[1:]))
    sin_s = np.sin(np.radians(lat_edges[:-1]))
    area_strip = R_EARTH_M ** 2 * np.radians(dlon) * (sin_n - sin_s)  # (nlat,)
    # If lat was descending, the strips are also computed in ascending
    # order — remap back to the original lat order.
    if lat[0] > lat[-1]:
        area_strip = area_strip[::-1]
    return np.repeat(area_strip[:, None], nlon, axis=1)


def read_scrip_areas(scrip_path: Path
                        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read ``(center_lat_deg, center_lon_deg, cell_area_m2)`` from SCRIP."""
    ds = xr.open_dataset(scrip_path, decode_cf=False)
    lat_units = str(ds["grid_center_lat"].attrs.get("units", "radians")).lower()
    if "rad" in lat_units:
        clat = np.degrees(ds["grid_center_lat"].values).astype("float64")
        clon = np.degrees(ds["grid_center_lon"].values).astype("float64")
    else:
        clat = ds["grid_center_lat"].values.astype("float64")
        clon = ds["grid_center_lon"].values.astype("float64")
    area = ds["grid_area"].values.astype("float64") * R_EARTH_M ** 2
    ds.close()
    # Normalise longitudes to [-180, 180] for consistent bbox masking
    clon = ((clon + 180.0) % 360.0) - 180.0
    return clat, clon, area


# ---------------------------------------------------------------------------
# Region masks + total integration
# ---------------------------------------------------------------------------

def bbox_mask_latlon(lat: np.ndarray, lon: np.ndarray,
                      bbox: tuple[float, float, float, float]) -> np.ndarray:
    """2-D bool mask (nlat, nlon) selecting cells inside ``bbox=(W,E,S,N)``.

    Half-open convention: ``south <= center_lat < north`` and analogous
    for lon.  This guarantees adjacent regions don't double-count cells
    that fall exactly on a shared boundary (e.g. a cell at lat=23.5°
    would otherwise be in both ``tropics`` and ``nh_extratropics``).
    """
    w, e, s, n = bbox
    lat_mask = (lat >= s) & (lat < n)
    # Handle wrap-around: if west > east, region crosses the dateline
    if w <= e:
        lon_mask = (lon >= w) & (lon < e)
    else:
        lon_mask = (lon >= w) | (lon < e)
    return lat_mask[:, None] & lon_mask[None, :]


def bbox_mask_ncol(center_lat: np.ndarray, center_lon: np.ndarray,
                    bbox: tuple[float, float, float, float]) -> np.ndarray:
    """1-D bool mask (ncol,), same half-open convention as bbox_mask_latlon."""
    w, e, s, n = bbox
    lat_mask = (center_lat >= s) & (center_lat < n)
    if w <= e:
        lon_mask = (center_lon >= w) & (center_lon < e)
    else:
        lon_mask = (center_lon >= w) | (center_lon < e)
    return lat_mask & lon_mask


def integrate_region_latlon(flux: np.ndarray, area_m2: np.ndarray,
                             mask: np.ndarray, duration_s: float
                             ) -> np.ndarray:
    """Integrate flux [molec/cm²/s] over cells where mask is True.

    ``flux`` shape (n_time, nlat, nlon); ``area_m2`` shape (nlat, nlon);
    ``mask`` shape (nlat, nlon).  Returns total molecules per timestep,
    shape (n_time,).
    """
    area_cm2 = area_m2 * 1e4                                   # m² → cm²
    weight = np.where(mask, area_cm2, 0.0)                     # (nlat, nlon)
    # sum over the spatial dims
    per_step = (flux * weight[None, :, :]).sum(axis=(1, 2))    # (n_time,)
    return per_step * duration_s


def integrate_region_ncol(flux: np.ndarray, area_m2: np.ndarray,
                           mask: np.ndarray, duration_s: float
                           ) -> np.ndarray:
    """Same, for unstructured (ncol) grids."""
    area_cm2 = area_m2 * 1e4
    weight = np.where(mask, area_cm2, 0.0)                     # (ncol,)
    per_step = (flux * weight[None, :]).sum(axis=1)            # (n_time,)
    return per_step * duration_s


# ---------------------------------------------------------------------------
# Unit conversion + plot
# ---------------------------------------------------------------------------

def totals_to_mass_kg(totals_molec: np.ndarray, mw_g_per_mol: float
                       ) -> np.ndarray:
    """Uniform gas & aerosol-with-MW=1 conversion:
    kg = molec × MW / (N_A × 1000).
    """
    return totals_molec * mw_g_per_mol / (AVOGADRO * 1000.0)


def pick_mass_unit(max_kg: float) -> tuple[str, float]:
    """Return (label, divisor_from_kg) matched to the magnitude."""
    if max_kg >= 1e12: return "Tg", 1e9
    if max_kg >= 1e9:  return "Gg", 1e6
    if max_kg >= 1e6:  return "Mg", 1e3
    return "kg", 1.0


def plot_time_series(timestamps: np.ndarray,
                      totals_by_region: dict[str, np.ndarray],
                      *, species: str, units_label: str,
                      out_path: Path, y_log: bool) -> None:
    """Line plot: one series per region."""
    fig, ax = plt.subplots(figsize=(11, 5.5))
    cmap = plt.get_cmap("tab10")
    ts = pd.to_datetime(timestamps)

    for i, (name, values) in enumerate(totals_by_region.items()):
        style = "-" if name == "global" else "--"
        lw = 2.0 if name == "global" else 1.4
        ax.plot(ts, values, style, linewidth=lw,
                color=cmap(i % 10), label=name, marker=".", markersize=4)

    ax.set_xlabel("time (UTC)")
    ax.set_ylabel(f"{species}  [{units_label}]")
    ax.set_title(f"{species} emissions per timestep, by region")
    ax.grid(True, alpha=0.3)
    if y_log:
        ax.set_yscale("log")
    ax.legend(loc="best", frameon=True, ncol=2 if len(totals_by_region) > 4 else 1)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def process_file(path: Path, *, regions: dict, scrip_data, args) -> None:
    log.info("--- %s ---", path.name)
    ds = xr.open_dataset(path, decode_times=False)

    species = _identify_species_var(ds)
    var = ds[species]
    units_in = str(var.attrs.get("units", "molecules/cm2/s"))
    mw = var.attrs.get("molecular_weight_g_per_mol")
    sp_type = var.attrs.get("species_type")
    log.info("  species=%s  units=%s  MW=%s  type=%s",
             species, units_in, mw, sp_type)

    timestamps = get_timestamps(ds)
    duration_s, is_hourly = infer_cadence(ds, timestamps)
    log.info("  %d timesteps, %s cadence (%d s per step): %s .. %s",
             len(timestamps), "hourly" if is_hourly else "daily",
             duration_s, timestamps[0], timestamps[-1])

    grid_type = detect_grid_type(ds)
    log.info("  grid type: %s", grid_type)

    # Cell areas + spatial coordinates
    if grid_type == "latlon":
        lat_vals = ds["lat"].values.astype("float64")
        lon_vals = ds["lon"].values.astype("float64")
        area_m2 = compute_latlon_cell_areas(lat_vals, lon_vals)
        # Normalise lon to [-180,180] for bbox masking
        lon_norm = ((lon_vals + 180.0) % 360.0) - 180.0
        flux = var.values.astype("float64")                  # (n_time, nlat, nlon)
    else:
        if scrip_data is None:
            raise RuntimeError(
                f"{path.name} is on an unstructured grid (ncol); pass "
                "--scrip <SCRIP file> so cell areas can be computed")
        clat, clon, area_m2 = scrip_data
        if len(clat) != ds.sizes["ncol"]:
            raise ValueError(
                f"SCRIP ncol ({len(clat)}) != emission file ncol "
                f"({ds.sizes['ncol']})")
        flux = var.values.astype("float64")                  # (n_time, ncol)

    # Compute totals per region
    totals_molec: dict[str, np.ndarray] = {}
    for name, bbox in regions.items():
        if grid_type == "latlon":
            mask = bbox_mask_latlon(lat_vals, lon_norm, bbox)
            per_step = integrate_region_latlon(flux, area_m2, mask, duration_s)
        else:
            mask = bbox_mask_ncol(clat, clon, bbox)
            per_step = integrate_region_ncol(flux, area_m2, mask, duration_s)
        totals_molec[name] = per_step

    # Convert to mass if we have a MW
    have_mw = (mw is not None) and float(mw) > 0
    is_number = (sp_type == "number") or (
        "number" in units_in.lower() and "molec" not in units_in.lower())

    if have_mw and not is_number:
        totals_mass_kg = {n: totals_to_mass_kg(v, float(mw))
                          for n, v in totals_molec.items()}
        max_kg = max((np.nanmax(np.abs(v)) for v in totals_mass_kg.values()),
                      default=0.0)
        unit_label, divisor = pick_mass_unit(max_kg)
        totals_plot = {n: v / divisor for n, v in totals_mass_kg.items()}
        y_label = f"{species}  [{unit_label} per timestep]"
    else:
        totals_mass_kg = None
        totals_plot = totals_molec
        unit_label = "number per timestep" if is_number else "molecules per timestep"
        y_label = f"{species}  [{unit_label}]"

    # --- Write CSV ---
    df = pd.DataFrame({"time": pd.to_datetime(timestamps)})
    # Also include YYYYMMDD and datesec for easy filtering
    df["date"] = np.array([int(f"{t.astype('datetime64[D]').astype(str).replace('-','')}")
                            for t in timestamps], dtype="int32")
    df["datesec"] = np.array(
        [int((t - t.astype("datetime64[D]")).astype("timedelta64[s]").astype(int))
         for t in timestamps], dtype="int32")
    # Totals in molecules (or number)
    for name, values in totals_molec.items():
        col = f"{name}_{'number' if is_number else 'molec'}"
        df[col] = values
    # Optional mass columns
    if totals_mass_kg is not None:
        for name, values in totals_mass_kg.items():
            df[f"{name}_kg"] = values

    csv_path = args.out_dir / f"{path.stem}_totals.csv"
    df.to_csv(csv_path, index=False)
    log.info("  wrote %s", csv_path.name)

    # Summary in log
    log.info("  totals over the file (all timesteps summed):")
    for name in regions:
        tot_m = float(totals_molec[name].sum())
        line = f"    {name:<18s}  {tot_m:.4e} {'number' if is_number else 'molec'}"
        if totals_mass_kg is not None:
            tot_kg = float(totals_mass_kg[name].sum())
            line += f"   {tot_kg / divisor:.4e} {unit_label.split()[0]}"
        log.info(line)

    # --- Plot ---
    if not args.no_plot:
        png_path = args.out_dir / f"{path.stem}_totals.png"
        plot_time_series(timestamps, totals_plot,
                          species=species, units_label=unit_label,
                          out_path=png_path, y_log=args.y_log)
        log.info("  wrote %s", png_path.name)

    ds.close()


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve regions
    regions: dict[str, tuple[float, float, float, float]] = {}
    if not args.region or args.include_defaults:
        regions.update(DEFAULT_REGIONS)
    if args.region:
        for name, w, e, s, n in args.region:
            regions[name] = (float(w), float(e), float(s), float(n))
    log.info("regions: %s", list(regions))

    # Pre-load SCRIP if needed (do it once)
    scrip_data = None
    if args.scrip:
        if not args.scrip.exists():
            log.error("SCRIP file not found: %s", args.scrip); return 1
        scrip_data = read_scrip_areas(args.scrip)
        log.info("SCRIP: %d cells, total area = %.4f × Earth surface",
                 len(scrip_data[0]),
                 scrip_data[2].sum() / (4 * np.pi * R_EARTH_M ** 2))

    # Process each input file
    for path in args.inputs:
        if not path.exists():
            log.error("input not found: %s", path); return 1
        try:
            process_file(path, regions=regions, scrip_data=scrip_data, args=args)
        except Exception as e:
            log.exception("failed to process %s: %s", path.name, e)
            return 1

    log.info("done — outputs in %s", args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
