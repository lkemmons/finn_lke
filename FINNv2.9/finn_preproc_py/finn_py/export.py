"""Final output writers.

Replaces ``code_anaconda/export_shp.py``.  Writes the final polygon
table out as both CSV (no geometry, for tabular consumers) and
GeoPackage / Shapefile (with geometry, for GIS consumers).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import geopandas as gpd
import pandas as pd

log = logging.getLogger(__name__)


def write_output(
    gdf: gpd.GeoDataFrame,
    out_dir: Path,
    out_basename: str,
    *,
    fields: Sequence[str] = (
        "v_lct", "f_lct", "v_tree", "v_herb", "v_bare", "v_regnum", "v_frp",
    ),
    date_definition: str = "UTC",
    write_shapefile: bool = False,
    write_geopackage: bool = True,
) -> dict[str, Path]:
    """Write the final output as CSV (+ GeoPackage and/or Shapefile).

    The CSV mirrors the column order produced by
    ``code_anaconda/export_shp.main`` (polyid, fireid, cen_lon, cen_lat,
    acq_date_<lst|utc>, area_sqkm, alg_agg, <fields>); the geospatial
    file additionally carries the polygon geometry and the cleanids
    array (the latter only in GeoPackage — shapefiles can't hold lists).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    acq_col = f"acq_date_{date_definition.lower()}"
    paths: dict[str, Path] = {}

    csv_path = out_dir / f"{out_basename}.csv"
    csv_cols = ["polyid", "fireid", "cen_lon", "cen_lat",
                acq_col, "area_sqkm", "alg_agg", *fields]
    df = gdf.copy()
    df[acq_col] = df["acq_date_use"]
    # Only the columns we have — be permissive if some join failed.
    csv_cols = [c for c in csv_cols if c in df.columns]
    df[csv_cols].to_csv(csv_path, index=False)
    paths["csv"] = csv_path
    log.info("wrote %s (%d rows)", csv_path, len(df))

    # Geometry-bearing output(s).
    geo_cols = csv_cols + ([] if "cleanids" not in df.columns else ["cleanids"])
    # cleanids is a list[int]; shapefiles can't store that, GPKG can.
    if write_geopackage:
        gpkg_path = out_dir / f"{out_basename}.gpkg"
        out_gdf = df.set_geometry("geom") if "geom" in df.columns else df.set_geometry(df.geometry.name)
        # cleanids as JSON string is safest across drivers
        if "cleanids" in out_gdf.columns:
            out_gdf = out_gdf.copy()
            out_gdf["cleanids"] = out_gdf["cleanids"].map(
                lambda v: ",".join(map(str, v)) if isinstance(v, list) else None
            )
        out_gdf[geo_cols + [out_gdf.geometry.name]].to_file(gpkg_path, driver="GPKG")
        paths["gpkg"] = gpkg_path
        log.info("wrote %s", gpkg_path)

    if write_shapefile:
        shp_path = out_dir / f"{out_basename}.shp"
        out_gdf = df.set_geometry("geom") if "geom" in df.columns else df.set_geometry(df.geometry.name)
        # Drop cleanids (can't store list); shapefile field-name truncation also applies.
        shp_cols = [c for c in csv_cols if c != "cleanids"]
        out_gdf[shp_cols + [out_gdf.geometry.name]].to_file(
            shp_path, driver="ESRI Shapefile"
        )
        paths["shp"] = shp_path
        log.info("wrote %s", shp_path)

    return paths
