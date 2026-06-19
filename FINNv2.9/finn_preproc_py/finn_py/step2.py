"""Step 2: raster zonal joins (LCT majority, VCF mean, regnum centroid).

Replaces:
  - code_anaconda/run_step2.py + the SQL it generates dynamically
  - the VCF-zonal-mean part of code_anaconda/run_vcf.py
"""
from __future__ import annotations

import logging
from collections.abc import Iterable
from pathlib import Path
from typing import Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask as rio_mask

from .config import RasterSpec

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Thematic (LCT) — majority class
# ---------------------------------------------------------------------------

def thematic_zonal_stats(
    polygons: gpd.GeoDataFrame, raster_path: Path, var: str,
) -> pd.DataFrame:
    """Per-polygon land-cover histogram, **long form**.

    Returns a DataFrame with one row per (polyid, LCT class) pair, i.e.
    a polygon containing two land-cover classes produces two rows.  This
    matches the original FINN's output schema (multiple rows per polygon
    when the polygon spans multiple LCT classes — see the user's diff:
    polyid 2664 had separate rows for ``v_lct=2, f_lct=0.867`` and
    ``v_lct=4, f_lct=0.133``).

    Columns
    -------
    polyid
        Index column on the input is propagated as a regular column.
    v_<var>
        The class value.
    f_<var>
        Fraction of the polygon covered by this class (0..1).
    cv_<var>
        Pixel count of this class.
    ct_<var>
        Total valid (non-nodata) pixel count in the polygon — same value
        on every row for a given polygon.
    r_<var>
        1-based rank of this class within the polygon (1 = majority).

    Polygons that don't overlap any valid raster pixel get a single
    null-filled row so the downstream join keeps them.

    Mirrors `mkcmd_insert_table_thematic` in run_step2.py.
    """
    rows: list[dict] = []
    with rasterio.open(raster_path) as src:
        nodata = src.nodata
        for poly_id, poly in zip(polygons.index, polygons.geometry):
            try:
                arr, _ = rio_mask(src, [poly.__geo_interface__],
                                  crop=True, all_touched=False, filled=True)
            except (ValueError, IndexError):
                # polygon does not overlap raster at all
                rows.append(_thematic_null(poly_id, var))
                continue
            data = arr[0]
            if nodata is not None:
                valid = data[data != nodata]
            else:
                valid = data.ravel()
            if valid.size == 0:
                rows.append(_thematic_null(poly_id, var))
                continue
            vals, counts = np.unique(valid, return_counts=True)
            order = np.argsort(-counts)  # descending by count
            vals, counts = vals[order], counts[order]
            tcnt = int(counts.sum())
            for rank, (v, c) in enumerate(zip(vals, counts), start=1):
                rows.append({
                    "polyid": poly_id,
                    f"v_{var}": int(v),
                    f"f_{var}": float(c) / tcnt,
                    f"cv_{var}": int(c),
                    f"ct_{var}": tcnt,
                    f"r_{var}": rank,
                })
    return pd.DataFrame(rows)


def _thematic_null(poly_id, var) -> dict:
    return {
        "polyid": poly_id,
        f"v_{var}": pd.NA, f"f_{var}": np.nan,
        f"cv_{var}": pd.NA, f"ct_{var}": pd.NA, f"r_{var}": pd.NA,
    }


# ---------------------------------------------------------------------------
# Continuous (VCF) — mean per band, no-data aware
# ---------------------------------------------------------------------------

def continuous_zonal_stats(
    polygons: gpd.GeoDataFrame,
    raster_path: Path,
    var_names: Sequence[str],
    bands: Sequence[int] | None = None,
    valid_range: tuple[float, float] | None = None,
) -> pd.DataFrame:
    """Per-polygon mean of each band.

    Returns a DataFrame indexed by polygons.index with columns
    `v_<name>` for each name in `var_names`.  `bands` defaults to
    1..len(var_names).

    The mean is taken over pixels that are (a) not equal to the
    GeoTIFF nodata tag and (b) inside ``valid_range`` if supplied.
    This is how MOD44B's special codes (200=water, 251/252=cloud/snow,
    253=fill) get excluded — without that masking they'd pollute the
    means and bias `v_tree` etc. upward for coastal / cloudy fires.

    Mirrors `mkcmd_insert_table_continuous` in run_step2.py (which uses
    `ST_SummaryStatsAgg(clp, band, true)`), but with the extra
    valid-range filter the PostGIS code didn't have.
    """
    if bands is None:
        bands = list(range(1, len(var_names) + 1))
    rows: list[dict] = []
    with rasterio.open(raster_path) as src:
        nodata = src.nodata
        for poly_id, poly in zip(polygons.index, polygons.geometry):
            try:
                arr, _ = rio_mask(
                    src, [poly.__geo_interface__],
                    crop=True, all_touched=False, filled=True,
                    indexes=bands,
                )
            except (ValueError, IndexError):
                rows.append({"polyid": poly_id,
                             **{f"v_{n}": np.nan for n in var_names}})
                continue
            row = {"polyid": poly_id}
            for name, band in zip(var_names, range(arr.shape[0])):
                data = arr[band]
                # Build the mask of pixels to *include* in the mean.
                mask = np.ones(data.shape, dtype=bool)
                if nodata is not None:
                    mask &= (data != nodata)
                if valid_range is not None:
                    lo, hi = valid_range
                    mask &= (data >= lo) & (data <= hi)
                if not mask.any():
                    row[f"v_{name}"] = np.nan
                    continue
                row[f"v_{name}"] = float(
                    np.nanmean(data[mask].astype("float64"))
                )
            rows.append(row)
    return pd.DataFrame(rows).set_index("polyid")


# ---------------------------------------------------------------------------
# Polygon (regnum) — region of polygon's centroid
# ---------------------------------------------------------------------------

def polygons_zonal_join(
    polygons: gpd.GeoDataFrame,
    region_layer_path: Path,
    var: str,
    variable_in: str,
) -> pd.DataFrame:
    """Per-polygon region number by spatial join on centroid.

    Mirrors `mkcmd_insert_table_polygons` in run_step2.py.
    """
    regions = gpd.read_file(region_layer_path)
    if variable_in not in regions.columns:
        raise RuntimeError(
            f"region layer is missing '{variable_in}' column; got {list(regions.columns)}"
        )
    cent = polygons.copy()
    cent["__cent__"] = cent.geometry.centroid
    cent = cent.set_geometry("__cent__")
    j = gpd.sjoin(cent, regions[[variable_in, "geometry"]],
                  how="left", predicate="within")
    # If a centroid hits multiple region polygons (e.g. boundaries),
    # keep the first match deterministically.
    j = j[~j.index.duplicated(keep="first")]
    out = pd.DataFrame({f"v_{var}": j[variable_in].astype("Int64")},
                       index=polygons.index)
    out.index.name = "polyid"
    return out


# ---------------------------------------------------------------------------
# Driver — assemble output table for one day
# ---------------------------------------------------------------------------

def zonal_join_one_day(
    work_div_day: gpd.GeoDataFrame,
    rasters: Iterable[RasterSpec],
    *,
    work_pnt: gpd.GeoDataFrame | None = None,
) -> gpd.GeoDataFrame:
    """For one day's `work_div` rows, compute every raster join and
    return the wide output GeoDataFrame.

    Polygons that span multiple LCT classes produce multiple output rows
    (one per class), with all other columns (geometry, VCF means, regnum,
    FRP, area) repeated.  This matches the original FINN's CSV schema.

    Equivalent of running the per-day `step2_work_*.sql` script.
    """
    if not len(work_div_day):
        return work_div_day.iloc[0:0].copy()

    polygons = work_div_day.set_index("polyid")
    # Per-polygon parts (indexed by polyid).  We may have at most one
    # long-form (multi-row-per-polyid) part: the thematic raster.
    thematic_long: pd.DataFrame | None = None
    per_polygon_parts: list[pd.DataFrame] = []

    for r in rasters:
        if r.kind == "thematic":
            if thematic_long is not None:
                # We assume at most one thematic raster, like the
                # original FINN's standard recipe.  If you ever add a
                # second, you'd need a many-to-many merge here.
                raise NotImplementedError(
                    "multiple thematic rasters not supported; the "
                    "original FINN only used one"
                )
            thematic_long = thematic_zonal_stats(polygons, r.path, r.variable)
        elif r.kind == "continuous":
            per_polygon_parts.append(continuous_zonal_stats(
                polygons, r.path, r.variables, valid_range=r.valid_range,
            ))
        elif r.kind == "polygons":
            per_polygon_parts.append(
                polygons_zonal_join(polygons, r.path, r.variable, r.variable_in)
            )
        elif r.kind == "input":
            if work_pnt is None:
                raise RuntimeError("kind='input' needs work_pnt")
            sub = work_pnt[work_pnt["polyid"].notna()][["polyid", r.variable_in]]
            agg = sub.groupby("polyid")[r.variable_in].mean()
            per_polygon_parts.append(pd.DataFrame({f"v_{r.variable}": agg}))
        else:
            raise RuntimeError(f"unknown raster kind: {r.kind}")

    # Start with the polygons; centroid coords are constant per polyid.
    polygons = polygons.copy()
    cents = polygons.geometry.centroid
    polygons["cen_lon"] = cents.x
    polygons["cen_lat"] = cents.y

    # Join every per-polygon part (1:1) onto the polygon table.
    for part in per_polygon_parts:
        polygons = polygons.join(part, how="left")

    # Now expand by the thematic part (1:many).  If thematic_long is
    # None (no thematic raster in the recipe), just emit polygons as-is.
    polygons = polygons.reset_index()
    if thematic_long is None:
        return polygons

    merged = polygons.merge(thematic_long, on="polyid", how="left")
    # Preserve geometry column as GeoSeries with CRS after the merge.
    merged = gpd.GeoDataFrame(
        merged, geometry="geom",
        crs=work_div_day.crs if hasattr(work_div_day, "crs") else "EPSG:4326",
    )
    return merged
