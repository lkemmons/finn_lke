#!/usr/bin/env python3
"""Plot a map of unstructured-grid emissions using SCRIP cell polygons.

Reads a gridded emissions NetCDF (produced by ``grid_emissions.py
--scrip``) and the corresponding SCRIP file, then renders one
filled-polygon map per species — each cell coloured by its emission
value, **no interpolated contouring**.  Suitable for MPAS, FV3, CAM-SE
and any other SCRIP-described grid up to a few million cells.

Examples
--------
    # All species in the file, default global PlateCarree, log color
    python plot_emissions_scrip.py \\
        FINNv2.9.2_NRTmodvrs_MOZART_CO_mxc_20260530.nc \\
        /glade/.../scrip_mxc.nc \\
        --out-dir ./maps

    # Subset of species, Robinson projection, linear scale
    python plot_emissions_scrip.py emis.nc scrip.nc \\
        --out-dir ./maps \\
        --species CO BC OC PM2.5 \\
        --projection robinson \\
        --linear

    # Regional zoom
    python plot_emissions_scrip.py emis.nc scrip.nc \\
        --out-dir ./maps \\
        --extent -130 -60 20 55

Output
------
One PNG per (species, time index) at ``--out-dir/<species>_<YYYYMMDD>.png``.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature

log = logging.getLogger("plot_emissions_scrip")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="plot_emissions_scrip.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("emissions", type=Path,
                   help="gridded emissions NetCDF (one file per species, "
                        "with `time` and `ncol` dims) from "
                        "grid_emissions.py --scrip")
    p.add_argument("scrip", type=Path,
                   help="SCRIP grid file (provides cell corner geometry)")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output directory for PNGs (created if missing)")

    p.add_argument("--species", nargs="+", default=None,
                   help="species to plot (default: every data_var on `ncol` "
                        "whose units are molecules/cm2/s or similar)")
    p.add_argument("--time-index", type=int, default=0,
                   help="time index to plot when input has time > 1 "
                        "(default: 0)")

    p.add_argument("--projection", default="platecarree",
                   choices=["platecarree", "robinson", "mollweide",
                            "northpolar", "southpolar"],
                   help="map projection (default: %(default)s)")
    p.add_argument("--extent", nargs=4, type=float, default=None,
                   metavar=("WEST", "EAST", "SOUTH", "NORTH"),
                   help="lon/lat extent in degrees (default: global)")

    p.add_argument("--cmap", default="YlOrRd",
                   help="matplotlib colormap (default: %(default)s)")
    p.add_argument("--vmin", type=float, default=None,
                   help="color minimum (default: 1st percentile of nonzero "
                        "data, or 0 in linear mode)")
    p.add_argument("--vmax", type=float, default=None,
                   help="color maximum (default: 99.9th percentile)")
    p.add_argument("--log", action="store_true", default=True,
                   help="log-scale color (default)")
    p.add_argument("--linear", action="store_false", dest="log",
                   help="linear-scale color")
    p.add_argument("--mask-below", type=float, default=0.0,
                   help="hide cells whose value is <= this (default: 0)")

    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--width-inches", type=float, default=12.0)
    p.add_argument("--coastlines", action="store_true", default=True,
                   help="draw coastlines (default)")
    p.add_argument("--no-coastlines", action="store_false", dest="coastlines")
    p.add_argument("-v", "--verbose", action="store_true")

    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# SCRIP I/O + geometry
# ---------------------------------------------------------------------------

def load_scrip_polygons(scrip_path: Path):
    """Read corner geometry from a SCRIP file.

    Returns
    -------
    corner_lon_deg : ndarray, shape (ncells, ncorners)
        Polygon-vertex longitudes in degrees, **unwrapped relative to
        the cell center** — so every vertex is within 180° of its own
        center, even for cells that straddle the dateline.  Values can
        therefore fall outside [-180, 180]; that's intentional and
        cartopy handles it.
    corner_lat_deg : ndarray, shape (ncells, ncorners)
    center_lon_deg : ndarray, shape (ncells,)
    center_lat_deg : ndarray, shape (ncells,)
    """
    ds = xr.open_dataset(scrip_path, decode_cf=False)
    log.info("SCRIP %s: %d cells, %d corners",
             scrip_path.name,
             ds.sizes.get("grid_size", -1),
             ds.sizes.get("grid_corners", -1))

    # Detect units
    units = str(ds["grid_corner_lat"].attrs.get("units", "radians")).lower()
    if "rad" in units:
        clat = np.degrees(ds["grid_corner_lat"].values)
        clon = np.degrees(ds["grid_corner_lon"].values)
        center_lat = np.degrees(ds["grid_center_lat"].values)
        center_lon = np.degrees(ds["grid_center_lon"].values)
    else:                                            # already degrees
        clat = ds["grid_corner_lat"].values.astype("float64")
        clon = ds["grid_corner_lon"].values.astype("float64")
        center_lat = ds["grid_center_lat"].values.astype("float64")
        center_lon = ds["grid_center_lon"].values.astype("float64")

    # Normalize longitudes to [-180, 180] first
    clon = ((clon + 180.0) % 360.0) - 180.0
    center_lon = ((center_lon + 180.0) % 360.0) - 180.0

    # Unwrap each polygon's corners so all are within ±180° of the center
    diff = clon - center_lon[:, None]
    clon = np.where(diff >  180.0, clon - 360.0, clon)
    clon = np.where(diff < -180.0, clon + 360.0, clon)

    ds.close()
    return clon, clat, center_lon, center_lat


def build_polygon_array(clon: np.ndarray, clat: np.ndarray,
                         center_lon: np.ndarray, values: np.ndarray,
                         seam: float = 175.0
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Build a (n_polygons, ncorners, 2) array suitable for PolyCollection.

    Cells whose centre is within ``seam`` degrees of ±180° are duplicated:
    one copy with the unwrapped longitudes, one copy shifted by ∓360°.
    That way the dateline-spanning ring of cells renders on both sides of
    a PlateCarree map without a visible gap, and the colormap value is
    repeated for each copy.

    Returns ``(polygon_xy[N, ncorners, 2], polygon_value[N])`` where
    ``N == ncells + 2 * (cells near seam)``.
    """
    # Stack lon/lat into (ncells, ncorners, 2)
    poly_xy = np.stack([clon, clat], axis=-1)

    near_seam = np.abs(np.abs(center_lon) - 180.0) < (180.0 - seam)
    if not np.any(near_seam):
        return poly_xy, values

    # Duplicate near-seam cells with longitudes shifted to the other side.
    shift_xy = poly_xy[near_seam].copy()
    shift = np.where(center_lon[near_seam] > 0, -360.0, +360.0)[:, None, None]
    shift_xy[:, :, 0:1] += shift                       # only shift the lon column

    poly_xy_all = np.concatenate([poly_xy, shift_xy], axis=0)
    values_all = np.concatenate([values, values[near_seam]], axis=0)
    log.info("  duplicated %d dateline-spanning cells for full coverage",
             int(near_seam.sum()))
    return poly_xy_all, values_all


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def get_projection(name: str) -> ccrs.Projection:
    return {
        "platecarree":  ccrs.PlateCarree(),
        "robinson":     ccrs.Robinson(),
        "mollweide":    ccrs.Mollweide(),
        "northpolar":   ccrs.NorthPolarStereo(),
        "southpolar":   ccrs.SouthPolarStereo(),
    }[name]


def pick_color_limits(values: np.ndarray, args) -> tuple[float, float]:
    """Choose sensible vmin/vmax from data percentiles.

    For log scale: default to a 6-decade dynamic range below the 99.9th
    percentile.  This avoids the matplotlib LogNorm pathology where a
    single near-zero outlier (e.g. floating-point noise from a long
    Gaussian tail) makes the entire interesting range collapse onto one
    color.  Users can override either bound explicitly.
    """
    nz = values[values > args.mask_below]
    if len(nz) == 0:
        return (1.0, 10.0)                            # arbitrary; nothing to draw
    if args.log:
        vmax = args.vmax if args.vmax is not None else float(np.percentile(nz, 99.9))
        if args.vmin is not None:
            vmin = args.vmin
        else:
            # 1st percentile, but floored at vmax / 1e6 (6 decades) and
            # at the actual minimum, so the colorbar never silently
            # gets a 30-decade range from synthetic outliers.
            p1 = float(np.percentile(nz, 1))
            vmin = max(p1, vmax / 1e6)
        if vmin <= 0:
            vmin = max(nz.min(), vmax / 1e8)
    else:
        vmin = args.vmin if args.vmin is not None else 0.0
        vmax = args.vmax if args.vmax is not None else float(np.percentile(nz, 99.9))
    return vmin, vmax


def plot_one(species: str, values: np.ndarray, units: str,
              poly_xy: np.ndarray, poly_value: np.ndarray,
              args, *, title: str, out_path: Path) -> None:
    """Render a single map."""
    proj = get_projection(args.projection)
    aspect = 0.55 if args.projection == "platecarree" else 0.6
    fig = plt.figure(figsize=(args.width_inches, args.width_inches * aspect))
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    if args.extent is not None:
        west, east, south, north = args.extent
        ax.set_extent([west, east, south, north], crs=ccrs.PlateCarree())
    else:
        ax.set_global()

    # Mask cells we don't want to draw
    drawn_value = poly_value.astype("float64").copy()
    drawn_value[drawn_value <= args.mask_below] = np.nan

    vmin, vmax = pick_color_limits(values, args)
    log.info("  %s: vmin=%.3g, vmax=%.3g (%s)",
             species, vmin, vmax, "log" if args.log else "linear")

    norm = (mpl.colors.LogNorm(vmin=vmin, vmax=vmax, clip=False) if args.log
            else mpl.colors.Normalize(vmin=vmin, vmax=vmax))

    coll = PolyCollection(
        poly_xy,                                       # (npoly, ncorners, 2)
        array=drawn_value,
        cmap=args.cmap, norm=norm,
        edgecolors="none",
        antialiased=False,                             # faster + cleaner at scale
        transform=ccrs.PlateCarree(),
    )
    ax.add_collection(coll)

    if args.coastlines:
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4, color="0.2")
        ax.add_feature(cfeature.BORDERS, linewidth=0.2, color="0.4")

    ax.gridlines(draw_labels=(args.projection == "platecarree"),
                  linewidth=0.3, color="gray", alpha=0.5)

    # Horizontal colorbar below the map
    cbar = fig.colorbar(coll, ax=ax, orientation="horizontal",
                         shrink=0.8, pad=0.05, extend="both")
    cbar.set_label(f"{species}   [{units}]")

    ax.set_title(title, fontsize=11)
    fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight",
                  facecolor="white")
    plt.close(fig)
    log.info("  wrote %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def detect_species(ds: xr.Dataset) -> list[str]:
    """Return data_vars defined on the ncol dimension."""
    candidates = []
    for v in ds.data_vars:
        # Skip the `date` housekeeping variable
        if v == "date":
            continue
        if "ncol" in ds[v].dims:
            candidates.append(v)
    return candidates


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )

    if not args.emissions.exists():
        log.error("emissions file not found: %s", args.emissions); return 1
    if not args.scrip.exists():
        log.error("SCRIP file not found: %s", args.scrip); return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Load SCRIP polygons (once, shared across species)
    clon, clat, center_lon, center_lat = load_scrip_polygons(args.scrip)
    ncells = clon.shape[0]
    log.info("loaded %d cells, %d corners each", ncells, clon.shape[1])

    # Load emissions
    log.info("reading emissions from %s", args.emissions)
    emis = xr.open_dataset(args.emissions, decode_times=False)
    if "ncol" not in emis.dims:
        log.error("emissions file %s has no 'ncol' dimension — was it "
                  "produced with grid_emissions.py --scrip?", args.emissions)
        return 1
    if emis.sizes["ncol"] != ncells:
        log.error("ncol mismatch: emissions has %d cells, SCRIP has %d",
                  emis.sizes["ncol"], ncells); return 1

    nt = emis.sizes.get("time", 1)
    if args.time_index >= nt:
        log.error("--time-index %d out of range (file has %d times)",
                  args.time_index, nt); return 1
    ti = args.time_index

    # Date for filename + title
    if "date" in emis.variables:
        d_int = int(emis["date"].values[ti])
        date_str = f"{d_int:08d}"                     # YYYYMMDD
    else:
        date_str = "0"
    log.info("plotting time index %d (date %s)", ti, date_str)

    species_list = args.species or detect_species(emis)
    if not species_list:
        log.error("no species variables found in %s", args.emissions); return 1
    log.info("species to plot: %s", species_list)

    for sp in species_list:
        if sp not in emis.variables:
            log.warning("species %s not in input; skipping", sp); continue
        values = emis[sp].isel(time=ti).values.astype("float64")
        units = str(emis[sp].attrs.get("units", ""))
        log.info("=== %s (%s) — %d nonzero of %d cells, range %.3g .. %.3g",
                 sp, units,
                 int((values > 0).sum()), len(values),
                 float(values[values > 0].min()) if (values > 0).any() else 0.0,
                 float(values.max()) if values.size else 0.0)

        poly_xy, poly_value = build_polygon_array(clon, clat, center_lon, values)

        title = f"{sp}  {date_str}   ({units})"
        out_path = args.out_dir / f"{sp}_{date_str}.png"
        plot_one(sp, values, units, poly_xy, poly_value,
                  args, title=title, out_path=out_path)

    log.info("done — %d maps written to %s", len(species_list), args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
