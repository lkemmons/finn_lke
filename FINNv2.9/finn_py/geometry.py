"""Geometry helpers.

Most of these wrap shapely/pyproj idioms that show up repeatedly in
step1/step2.  Kept here to avoid scattering coordinate-system trivia.
"""
from __future__ import annotations

import numpy as np
import shapely
from pyproj import Geod
from shapely.geometry import Polygon, box

_GEOD = Geod(ellps="WGS84")


def spheroidal_area_km2(geom) -> float:
    """Equivalent to PostGIS `ST_Area(geom, true) / 1e6`.

    Uses pyproj's WGS84 geodesic.  Accepts a single geometry; returns the
    absolute value (geometry_area_perimeter returns a signed area whose
    sign depends on ring orientation).
    """
    if geom is None or geom.is_empty:
        return 0.0
    area_m2, _ = _GEOD.geometry_area_perimeter(geom)
    return abs(area_m2) / 1.0e6


def spheroidal_area_km2_vec(geoms) -> np.ndarray:
    """Vectorised version returning a numpy array."""
    out = np.empty(len(geoms), dtype="float64")
    for i, g in enumerate(geoms):
        out[i] = spheroidal_area_km2(g)
    return out


def lat_lon_box(lon: float, lat: float, dx: float, dy: float) -> Polygon:
    """A lon/lat rectangle.  `dx`, `dy` are half-widths in degrees."""
    return box(lon - dx, lat - dy, lon + dx, lat + dy)


def make_pixel_dxdy(
    lat: np.ndarray,
    fire_size_km: np.ndarray,
    pix_scan_km: np.ndarray | None = None,
    pix_track_km: np.ndarray | None = None,
    pixfac: float = 1.1,
    great_circle_km: float = 2.0 * np.pi * 6370.997,
) -> dict[str, np.ndarray]:
    """Half-widths (degrees) for the small fire square and the pixel
    envelope, using the same latitude correction as step1a_work_v7m.sql.

    If `pix_scan_km` / `pix_track_km` are omitted, the pixel box collapses
    onto the fire box (the conservative definition used by step 1b).
    """
    cos_lat = np.cos(np.radians(lat))
    fire_dx = 0.5 * fire_size_km * 360.0 / great_circle_km / cos_lat
    fire_dy = 0.5 * fire_size_km * 360.0 / great_circle_km
    if pix_scan_km is None or pix_track_km is None:
        pix_dx, pix_dy = fire_dx.copy(), fire_dy.copy()
    else:
        pix_dx = pixfac * 0.5 * pix_scan_km * 360.0 / great_circle_km / cos_lat
        pix_dy = pixfac * 0.5 * pix_track_km * 360.0 / great_circle_km
    return {"fire_dx": fire_dx, "fire_dy": fire_dy,
            "pix_dx": pix_dx, "pix_dy": pix_dy}


def fill_small_holes(
    polygon, area_threshold_deg2: float
) -> Polygon:
    """Drop interior rings whose envelope area or area is below threshold.

    Matches step 2.4 of step1a_work_v7m.sql:
        - first fill all holes (take exterior ring)
        - then re-punch holes whose envelope area > threshold AND area > threshold.
    `area_threshold_deg2` is in degree^2 (the SQL uses (1/240)^2).
    """
    if polygon.is_empty or polygon.geom_type == "Point":
        return polygon
    if polygon.geom_type == "MultiPolygon":
        return shapely.ops.unary_union(
            [fill_small_holes(p, area_threshold_deg2) for p in polygon.geoms]
        )
    # Polygon
    if not polygon.interiors:
        return polygon
    keep_rings = []
    for ring in polygon.interiors:
        hole = Polygon(ring)
        if (hole.envelope.area > area_threshold_deg2
                and hole.area > area_threshold_deg2):
            keep_rings.append(ring)
    return Polygon(polygon.exterior, holes=keep_rings)


def polsby_popper(geom) -> float:
    """Polsby-Popper compactness.  1 = circle, 0 = degenerate line.

    Equivalent to step1_prep_v7m.sql's st_polsbypopper(geom).
    """
    if geom is None or geom.is_empty:
        return 0.0
    a = geom.area
    p = geom.length
    if a == 0 or p == 0:
        return 0.0
    return 4.0 * np.pi * a / (p * p)
