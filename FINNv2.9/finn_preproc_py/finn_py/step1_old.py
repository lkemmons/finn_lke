"""Step 1: build work_pnt, cluster into work_lrg, subdivide into work_div.

Replaces:
  - code_anaconda/step1_prep_v7m.sql        (prep)
  - code_anaconda/step1a_work_v7m.sql       (aggressive aggregation, per day)
  - code_anaconda/step1b_prep_v7m.sql       (no-op once prep is in-mem)
  - code_anaconda/step1b_work_v7m.sql       (conservative + Voronoi subdiv)
  - code_anaconda/run_step1a.py / run_step1b.py / run_vcf.py drivers
"""
from __future__ import annotations

import datetime as dt
import logging
from typing import Sequence

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import shapely
from scipy.spatial import KDTree, Voronoi
from shapely.geometry import MultiPoint, Point, Polygon, box
from shapely.ops import unary_union

from .config import (
    FIRE_SIZE_MODIS_KM, FIRE_SIZE_VIIRS_KM, MODIS_CONFIDENCE_MIN,
    PERSISTENT_ANOMTYPES, PIXFAC, SKIM_RADIUS_DEG, SMALL_HOLE_AREA_DEG2,
    TROPICS_LAT_DEG, VIIRS_CONFIDENCE_DROP,
)
from .geometry import (
    fill_small_holes, make_pixel_dxdy, polsby_popper, spheroidal_area_km2_vec,
)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# PREP — replaces step1_prep_v7m.sql Parts 1, 3
# ---------------------------------------------------------------------------

def _instrument_from_satellite(sat) -> str | None:
    """Mirror of step1_prep's get_instrument(satellite character)."""
    if sat is None or (isinstance(sat, float) and np.isnan(sat)):
        return None
    if isinstance(sat, bool):  # GDAL sometimes coerces 'N' to False
        return None if sat else "VIIRS"
    s = str(sat).strip().upper()
    if not s:
        return None
    first = s[0]
    return {"T": "MODIS", "A": "MODIS", "N": "VIIRS"}.get(first)


def _is_confident(instrument: str | None, confidence) -> bool:
    """Instrument-specific confidence rule (cf. tbl_flddefs in prep)."""
    if instrument is None or confidence is None:
        return False
    if instrument == "MODIS":
        try:
            return int(confidence) >= MODIS_CONFIDENCE_MIN
        except (TypeError, ValueError):
            return False
    if instrument == "VIIRS":
        return str(confidence).strip().lower()[:1] != VIIRS_CONFIDENCE_DROP
    return False


def _acq_datetime_lst(date, time_hhmm: str, lon: float) -> pd.Timestamp:
    """Mirror of get_acq_datetime_lst(date, time, lon)."""
    h = int(time_hhmm[:2])
    m = int(time_hhmm[2:4])
    base = pd.Timestamp(date) + pd.Timedelta(hours=h, minutes=m)
    offset_hours = int(round(lon / 15.0))
    return base + pd.Timedelta(hours=offset_hours)


def _normalize_time(t) -> str:
    """Normalize 'HH:MM' or 'HHMM' or int 1234 to 'HHMM'."""
    if t is None or (isinstance(t, float) and np.isnan(t)):
        return "0000"
    s = str(t).replace(":", "").strip()
    return s.zfill(4)[:4]


def prep(
    af: gpd.GeoDataFrame,
    *,
    date_definition: str = "UTC",
    date_range: tuple[dt.date, dt.date] | None = None,
    filter_persistent_sources: bool = True,
) -> gpd.GeoDataFrame:
    """Build `work_pnt` from raw AF detections.

    Equivalent of `step1_prep_v7m.sql` Part 3, minus everything about
    table creation and SQL function definitions.
    """
    df = af.copy()

    # Geometry, lon/lat (denormalised).
    df["lon"] = df.geometry.x.astype("float64")
    df["lat"] = df.geometry.y.astype("float64")
    df = df.rename(columns={df.geometry.name: "geom_pnt"})
    df = df.set_geometry("geom_pnt")

    # Numeric scan/track (some FIRMS products store as strings).
    df["scan"] = pd.to_numeric(df["scan"], errors="coerce")
    df["track"] = pd.to_numeric(df["track"], errors="coerce")

    # Instrument: archive files supply this column directly ('MODIS' /
    # 'VIIRS'); NRT files don't (each file is one instrument), so derive
    # from the satellite character.  Where both are present, the explicit
    # column wins; we only fill in nulls from the derivation.
    derived = df["satellite"].map(_instrument_from_satellite)
    if "instrument" in df.columns:
        inst = df["instrument"].astype("string").str.upper().str.strip()
        # Normalise common variants: 'Modis', 'V', 'VIIRS' etc.
        inst = inst.replace({"V": "VIIRS", "M": "MODIS"})
        df["instrument"] = inst.where(inst.isin(["MODIS", "VIIRS"]), derived)
    else:
        df["instrument"] = derived
    assert df["instrument"].notna().all(), "could not derive instrument for some rows"
    df["confident"] = [
        _is_confident(i, c) for i, c in zip(df["instrument"], df["confidence"])
    ]

    # Dates.
    df["acq_date_utc"] = pd.to_datetime(df["acq_date"]).dt.date
    df["acq_time_utc"] = df["acq_time"].map(_normalize_time)
    df["acq_datetime_lst"] = [
        _acq_datetime_lst(d, t, lon)
        for d, t, lon in zip(df["acq_date_utc"], df["acq_time_utc"], df["lon"])
    ]
    df["acq_date_lst"] = df["acq_datetime_lst"].dt.date

    if date_definition == "UTC":
        df["acq_date_use"] = df["acq_date_utc"]
    elif date_definition == "LST":
        df["acq_date_use"] = df["acq_date_lst"]
    else:
        raise ValueError(f"date_definition must be 'UTC' or 'LST': {date_definition!r}")

    # Anomaly type (FIRMS 'Type' field — 0 if absent).
    df["anomtype"] = (
        df["type"].astype("Int64") if "type" in df.columns else 0
    )
    has_type_col = df["anomtype"].notna().any()

    df["frp"] = pd.to_numeric(df.get("frp"), errors="coerce")
    df["rawid"] = np.arange(len(df), dtype="int64")

    # --- Tropics duplication (step1_prep "dup tropics") ---
    trop_mask = (df["instrument"] == "MODIS") & (df["lat"].abs() <= TROPICS_LAT_DEG)
    if trop_mask.any():
        dup = df.loc[trop_mask].copy()
        oneday = pd.Timedelta(days=1)
        dup["acq_date_utc"] = dup["acq_date_utc"] + oneday
        dup["acq_date_lst"] = dup["acq_date_lst"] + oneday
        dup["acq_datetime_lst"] = dup["acq_datetime_lst"] + oneday
        dup["acq_date_use"] = dup["acq_date_use"] + oneday
        df = pd.concat([df, dup], ignore_index=True)
        log.info("duplicated %d tropical MODIS rows to next day", trop_mask.sum())

    # --- Date filter ---
    if date_range is not None:
        first, last = date_range
        if first is not None:
            df = df[df["acq_date_use"] >= first]
        if last is not None:
            df = df[df["acq_date_use"] <= last]

    # --- Confidence filter ---
    n0 = len(df)
    df = df[df["confident"]]
    log.info("dropped low confidence: %d -> %d", n0, len(df))

    # --- Persistent-source filter (only if Type column is present) ---
    if filter_persistent_sources and has_type_col:
        n0 = len(df)
        df = df[~df["anomtype"].isin(PERSISTENT_ANOMTYPES)]
        log.info("dropped persistent (anomtype in %s): %d -> %d",
                 PERSISTENT_ANOMTYPES, n0, len(df))

    # Final PK + downstream-needed columns.
    df = df.reset_index(drop=True)
    df["cleanid"] = np.arange(len(df), dtype="int64")
    for col in ("fireid1", "fireid2", "ndetect1", "ndetect2",
                "alg_agg", "polyid"):
        if col not in df.columns:
            df[col] = pd.NA

    keep = [
        "cleanid", "rawid", "src_file", "geom_pnt", "lon", "lat",
        "scan", "track", "acq_date_utc", "acq_time_utc",
        "acq_date_lst", "acq_datetime_lst", "acq_date_use",
        "instrument", "confident", "anomtype", "frp",
        "fireid1", "ndetect1", "fireid2", "ndetect2", "alg_agg", "polyid",
    ]
    keep = [c for c in keep if c in df.columns]
    out = gpd.GeoDataFrame(df[keep], geometry="geom_pnt", crs="EPSG:4326")
    return out


# ---------------------------------------------------------------------------
# STEP 1A — aggressive aggregation, per day
# ---------------------------------------------------------------------------

def _fire_size_km(instrument: pd.Series) -> np.ndarray:
    return instrument.map(
        {"MODIS": FIRE_SIZE_MODIS_KM, "VIIRS": FIRE_SIZE_VIIRS_KM}
    ).to_numpy(dtype="float64")


def _adjacency_pairs(
    lon: np.ndarray, lat: np.ndarray,
    pix_dx: np.ndarray, pix_dy: np.ndarray,
) -> np.ndarray:
    """Pairs (i, j) with i<j whose pixel envelopes overlap in lon/lat.

    Equivalent to step1a's `tbl_adj_det`:
        a.geom_pix && b.geom_pix
        AND abs(a.lon-b.lon) < a.pix_dx+b.pix_dx
        AND abs(a.lat-b.lat) < a.pix_dy+b.pix_dy

    We use a KDTree on (lon, lat) with a conservative search radius
    (the largest possible euclidean distance at which any two rectangles
    in the input set could still touch) and then filter with the exact
    rectangle test in Python.  Much faster than O(n²) and correct for
    any input.
    """
    if len(lon) == 0:
        return np.zeros((0, 2), dtype=np.int64)
    pts = np.column_stack([lon, lat])
    # Upper bound on |a-b| for any pair to still have overlapping
    # envelopes: 2 * (max(pix_dx) + max(pix_dy)).  Loose, but cheap to
    # filter exactly afterwards and avoids missing edge cases.
    radius = 2.0 * (float(np.max(pix_dx)) + float(np.max(pix_dy)))
    tree = KDTree(pts)
    pairs = tree.query_pairs(r=radius, output_type="ndarray")
    if len(pairs) == 0:
        return pairs.astype(np.int64)
    # exact rectangular overlap test
    i, j = pairs[:, 0], pairs[:, 1]
    keep = (np.abs(lon[i] - lon[j]) < (pix_dx[i] + pix_dx[j])) & \
           (np.abs(lat[i] - lat[j]) < (pix_dy[i] + pix_dy[j]))
    return pairs[keep].astype(np.int64)


def _connected_components(
    n_nodes: int, edges: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Returns (fireid, ndetect) arrays of length n_nodes.

    fireid = smallest local index within the node's component.
    Lone nodes get their own index as fireid and ndetect=1.

    This is the Python equivalent of `pnt2grp`.
    """
    g = nx.Graph()
    g.add_nodes_from(range(n_nodes))
    if len(edges):
        g.add_edges_from(edges.tolist())
    fireid = np.empty(n_nodes, dtype=np.int64)
    ndetect = np.empty(n_nodes, dtype=np.int64)
    for cc in nx.connected_components(g):
        m = min(cc)
        n = len(cc)
        for v in cc:
            fireid[v] = m
            ndetect[v] = n
    return fireid, ndetect


def step1a_one_day(
    work_pnt: gpd.GeoDataFrame, oned: dt.date
) -> gpd.GeoDataFrame:
    """Aggressive clustering for one day.  Returns a `work_lrg1`-shape
    GeoDataFrame for that day (one row per fireid).

    Also writes `fireid1` and `ndetect1` back into `work_pnt` *in place*
    on the subset rows.

    Equivalent of step1a_work_v7m.sql, except that it returns the day's
    new `work_lrg1` rows instead of `INSERT INTO work_lrg1`.
    """
    mask = work_pnt["acq_date_use"] == oned
    if not mask.any():
        return _empty_work_lrg()

    day = work_pnt.loc[mask].copy().reset_index(drop=False)  # 'index' = work_pnt row id
    n = len(day)

    fire_size = _fire_size_km(day["instrument"])
    dxdy = make_pixel_dxdy(
        lat=day["lat"].to_numpy(),
        fire_size_km=fire_size,
        pix_scan_km=day["scan"].to_numpy(),
        pix_track_km=day["track"].to_numpy(),
        pixfac=PIXFAC,
    )

    # Small-fire rectangles (output geometry) and pixel rectangles (adjacency).
    lon = day["lon"].to_numpy()
    lat = day["lat"].to_numpy()
    geom_sml = shapely.box(
        lon - dxdy["fire_dx"], lat - dxdy["fire_dy"],
        lon + dxdy["fire_dx"], lat + dxdy["fire_dy"],
    )

    pairs = _adjacency_pairs(lon, lat, dxdy["pix_dx"], dxdy["pix_dy"])

    # local fireid (within day): smallest local index per component
    local_fireid, ndetect = _connected_components(n, pairs)

    # Convert local fireid back to cleanid-based fireid for compatibility.
    cleanid = day["cleanid"].to_numpy()
    # fireid in step1 is "smallest cleanid in the group"
    fireid = np.empty(n, dtype=np.int64)
    for cc_min_local in np.unique(local_fireid):
        members = local_fireid == cc_min_local
        fireid[members] = cleanid[members].min()

    # Push fireid1 / ndetect1 back into work_pnt.
    work_pnt.loc[mask, "fireid1"] = fireid
    work_pnt.loc[mask, "ndetect1"] = ndetect

    # Build large polygons.
    # For ndetect>1 components: union of pairwise convex hulls (matches
    # step1a's `tbl_adj_det.geom_pair = st_convexhull(st_collect(...))`
    # then `st_union(geom_pair) group by fireid`).
    # For lone detections: just geom_sml.
    large_rows: list[dict] = []
    if len(pairs) > 0:
        # Build pair convex hulls in bulk via shapely.
        pair_hulls = []
        pair_fireids = []
        # group pairs by fireid (use local_fireid of one endpoint)
        local_fid_for_pair = local_fireid[pairs[:, 0]]
        for fid_local in np.unique(local_fid_for_pair):
            edge_mask = local_fid_for_pair == fid_local
            edges_for_fid = pairs[edge_mask]
            hulls = []
            for ii, jj in edges_for_fid:
                hulls.append(
                    shapely.convex_hull(
                        shapely.union(geom_sml[ii], geom_sml[jj])
                    )
                )
            polygon = unary_union(hulls)
            # Fill small holes (step 2.4 of original)
            polygon = fill_small_holes(polygon, SMALL_HOLE_AREA_DEG2)
            members = local_fireid == fid_local
            members_cid = cleanid[members]
            large_rows.append({
                "fireid": int(members_cid.min()),
                "geom_lrg": polygon,
                "acq_date_use": oned,
                "ndetect": int(members.sum()),
                "alg_agg": pd.NA,
            })
    # Lone detections
    lone_mask = ndetect == 1
    for k in np.flatnonzero(lone_mask):
        large_rows.append({
            "fireid": int(cleanid[k]),
            "geom_lrg": geom_sml[k],
            "acq_date_use": oned,
            "ndetect": 1,
            "alg_agg": pd.NA,
        })

    out = gpd.GeoDataFrame(large_rows, geometry="geom_lrg", crs="EPSG:4326")
    if len(out):
        out["area_sqkm"] = spheroidal_area_km2_vec(out["geom_lrg"].values)
    else:
        out["area_sqkm"] = pd.Series(dtype="float64")
    return out


def _empty_work_lrg() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "fireid": pd.Series(dtype="int64"),
            "geom_lrg": gpd.GeoSeries([], crs="EPSG:4326"),
            "acq_date_use": pd.Series(dtype="object"),
            "ndetect": pd.Series(dtype="int64"),
            "area_sqkm": pd.Series(dtype="float64"),
            "alg_agg": pd.Series(dtype="Int64"),
        },
        geometry="geom_lrg",
        crs="EPSG:4326",
    )


# ---------------------------------------------------------------------------
# Tree-cover gating between 1a and 1b (replaces run_vcf.py)
# ---------------------------------------------------------------------------

def set_alg_agg_from_tree(
    work_pnt: gpd.GeoDataFrame,
    work_lrg1: gpd.GeoDataFrame,
    v_tree: pd.Series,           # indexed by fireid
    threshold_pct: float = 50.0,
) -> None:
    """In-place update of `alg_agg` on work_pnt and work_lrg1.

    alg_agg = 1 (aggressive) if v_tree >= threshold else 2 (conservative).

    Polygons with v_tree = NaN (no valid VCF pixels — e.g. the polygon
    is entirely covered by MOD44B's water/cloud/snow/fill codes, or
    falls outside the VCF mosaic's coverage) default to alg_agg = 2.
    This matches the original FINN's behavior; downstream emissions
    calculations use LCT to handle the no-VCF case.

    `v_tree` is computed by `step2.continuous_zonal_stats` against the
    aggressive work_lrg1 polygons.
    """
    # NaN >= threshold is False, so np.where assigns 2 to NaN entries.
    # We then fillna(2) defensively to catch any fireids in work_lrg1
    # that didn't appear in v_tree's index at all (would otherwise map
    # to NaN through the .map call below).
    alg = pd.Series(
        np.where(v_tree.values >= threshold_pct, 1, 2),
        index=v_tree.index, dtype="Int64",
    )
    mapped = work_lrg1["fireid"].map(alg)
    n_unknown = int(mapped.isna().sum())
    work_lrg1["alg_agg"] = mapped.fillna(2).astype("Int64")
    if n_unknown:
        log.info(
            "alg_agg: %d polygon(s) had no v_tree match; defaulted to 2 (conservative)",
            n_unknown,
        )

    # propagate to work_pnt via fireid1, again with fillna(2) as a safety net
    fid_to_alg = dict(zip(work_lrg1["fireid"], work_lrg1["alg_agg"]))
    work_pnt["alg_agg"] = (
        work_pnt["fireid1"].map(fid_to_alg).fillna(2).astype("Int64")
    )


# ---------------------------------------------------------------------------
# STEP 1B — conservative aggregation + Voronoi subdivision, per day
# ---------------------------------------------------------------------------

def _pnt2drop(
    edges: np.ndarray, invdist: np.ndarray, n_nodes: int
) -> tuple[list[int], list[list[int]]]:
    """Iteratively eliminate the highest-density point from each
    connected component until no edges remain.  Ties get merged to a
    centroid placeholder.

    This is a faithful Python re-port of the PL/Python `pnt2drop`
    function in step1_prep_v7m.sql.

    Returns:
        todrop:  list of node indices that were dropped
        torepl:  parallel list, each item is a list of *other* nodes
                 that were merged with the dropped one (centroid replacement)
    """
    g = nx.Graph()
    for (l, r), inv in zip(edges, invdist):
        g.add_edge(int(l), int(r), invdist=float(inv))

    todrop: list[int] = []
    torepl_pairs: list[list[int]] = []
    added: dict[int, list[int]] = {}  # added node id -> dropped originals
    iadd = 0

    for cc_nodes in nx.connected_components(g):
        cc = g.subgraph(cc_nodes).copy()
        # Sanity bound on iterations.
        for _ in range(2 * (cc.size() + 1)):
            if cc.size() == 0:
                break
            # node with max sum of invdist among neighbors
            m = 0.0
            mydrop: list[int] = []
            for node, nbrs in cc.adjacency():
                s = sum(d["invdist"] for d in nbrs.values())
                if s > m:
                    mydrop = [node]
                    m = s
                elif s == m:
                    mydrop.append(node)
            if not mydrop:
                break

            # Ties (groups of 2+ tied nodes) get merged to a placeholder.
            sg = cc.subgraph(mydrop)
            myrepl = [list(c) for c in nx.connected_components(sg) if len(c) > 1]
            for repl in myrepl:
                newnbr: dict[int, float] = {}
                for nd in repl:
                    for nbr in cc.neighbors(nd):
                        if nbr in repl:
                            continue
                        w = cc.get_edge_data(nd, nbr)["invdist"]
                        newnbr[nbr] = max(newnbr.get(nbr, 0.0), w)
                if newnbr:
                    iadd -= 1
                    for nbr, w in newnbr.items():
                        cc.add_edge(iadd, nbr, invdist=w)
                    added[iadd] = repl

            todrop.extend(mydrop)
            torepl_pairs.extend(myrepl)
            cc.remove_nodes_from(mydrop)

    # If a placeholder we added was later dropped, undo its replacement.
    for idrp in todrop:
        if idrp >= 0:
            continue
        if added.get(idrp) in torepl_pairs:
            torepl_pairs.remove(added[idrp])

    # Recursively resolve placeholders to originals.
    def getorig(ns):
        out = []
        for n_ in ns:
            if n_ >= 0:
                out.append(n_)
            else:
                out.extend(getorig(added[n_]))
        return out

    torepl2: dict[int, list[int]] = {}
    for repl in torepl_pairs:
        orig = getorig(repl)
        mn = min(orig)
        torepl2[mn] = [x for x in orig if x != mn]

    others = [torepl2.get(d, []) for d in todrop]
    return todrop, others


def _voronoi_finite_polygons_2d(vor: Voronoi, radius: float | None = None):
    """SciPy Voronoi gives infinite ridges; reconstruct finite polygons
    by closing them off at a far radius.  Adapted from the recipe used
    in step1_prep_v7m.sql's st_voronoi_python.
    """
    new_regions = []
    new_vertices = vor.vertices.tolist()
    center = vor.points.mean(axis=0)
    if radius is None:
        radius = np.ptp(vor.points, axis=0).max() * 3

    all_ridges: dict[int, list[tuple[int, int, int]]] = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]
        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]
        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                continue
            t = vor.points[p2] - vor.points[p1]
            t = t / np.linalg.norm(t)
            n = np.array([-t[1], t[0]])
            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius
            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())
        # sort vertices counter-clockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = [v for _, v in sorted(zip(angles, new_region))]
        new_regions.append(new_region)
    return new_regions, np.asarray(new_vertices)


def _voronoi_cells(seeds_xy: np.ndarray) -> list[Polygon]:
    """One finite Voronoi polygon per seed point.  Requires >=4 points.

    Deprecated — kept only for reference.  Use `_voronoi_cells_robust`
    below: the radius computed by `_voronoi_finite_polygons_2d` from the
    seed extent is too small whenever seeds are tightly clustered
    relative to the parent polygon.
    """
    vor = Voronoi(seeds_xy)
    regions, vertices = _voronoi_finite_polygons_2d(vor)
    return [Polygon(vertices[r]) for r in regions]


def _voronoi_cells_robust(
    seeds_xy: np.ndarray,
    parent_bounds: tuple[float, float, float, float] | None = None,
) -> list[Polygon]:
    """One Voronoi cell per *unique* seed, robust to small n and tight clustering.

    Uses four "ghost" seed points placed well outside the parent's
    bounding box, so scipy.Voronoi yields *finite* cells for the real
    seeds without us guessing a clip radius from the seed extent (which
    collapses when seeds are tightly clustered relative to the parent —
    the bug that caused fires with 100% area shortfall).  Handles n=1,
    2, 3, 4+ uniformly.

    Duplicate seeds (exactly co-located) collapse to a single cell;
    any "duplicate" detections are reattached to the surviving polygon
    by the downstream sjoin (which uses ``predicate='intersects'``).
    """
    n = len(seeds_xy)
    if n == 0:
        return []

    # Dedupe — scipy.Voronoi can't handle duplicate input points.
    unique = np.unique(seeds_xy, axis=0)
    nu = len(unique)

    # Sizing for the ghost box: must be wider than the parent polygon AND
    # wider than the seed spread, so all real cells stay finite.
    if parent_bounds is not None:
        px = parent_bounds[2] - parent_bounds[0]
        py = parent_bounds[3] - parent_bounds[1]
    else:
        px = py = 0.0
    seed_x = float(np.ptp(unique[:, 0]))
    seed_y = float(np.ptp(unique[:, 1]))
    pad = max(px, py, seed_x, seed_y, 1e-3) * 100.0 + 1.0

    if nu == 1:
        # All seeds at the same location → one cell (caller clips to parent).
        c = unique[0]
        return [box(c[0] - pad, c[1] - pad, c[0] + pad, c[1] + pad)]

    minx, miny = unique.min(axis=0)
    maxx, maxy = unique.max(axis=0)
    ghost = np.array([
        [minx - pad, miny - pad],
        [maxx + pad, miny - pad],
        [maxx + pad, maxy + pad],
        [minx - pad, maxy + pad],
    ])
    pts = np.vstack([unique, ghost])

    try:
        vor = Voronoi(pts)
    except Exception as e:
        log.warning("Voronoi failed for %d unique seeds + ghosts (%s); "
                    "falling back to single cell covering all", nu, e)
        return [box(minx - pad, miny - pad, maxx + pad, maxy + pad)]

    cells: list[Polygon] = []
    for i in range(nu):
        region_idx = vor.point_region[i]
        region = vor.regions[region_idx]
        if -1 in region or len(region) < 3:
            cells.append(box(minx - pad, miny - pad, maxx + pad, maxy + pad))
            continue
        verts = vor.vertices[region]
        cells.append(Polygon(verts))
    return cells


def _custom_cutter(seeds_xy: np.ndarray) -> list[Polygon]:
    """Custom partition for 2- or 3-point cases (Voronoi needs >=4).

    Replaces `st_cutter_py` in step1_prep_v7m.sql.  Returns one polygon
    per input seed, each large enough to clip against the parent's
    large polygon.
    """
    n = len(seeds_xy)
    if n == 1:
        # Single point — return a large bounding box; caller clips.
        r = 1.0
        c = seeds_xy[0]
        return [box(c[0] - r, c[1] - r, c[0] + r, c[1] + r)]
    if n == 2:
        # Perpendicular bisector splits the world into two half-planes.
        a, b = seeds_xy
        mid = (a + b) / 2.0
        d = b - a
        # bisector direction
        n_vec = np.array([-d[1], d[0]])
        n_vec /= np.linalg.norm(n_vec)
        radius = max(1.0, np.linalg.norm(d) * 10)
        far1 = mid + n_vec * radius
        far2 = mid - n_vec * radius
        # build two big triangles from each seed
        return [
            Polygon([a, far1, far2]),
            Polygon([b, far2, far1]),
        ]
    if n == 3:
        # Use the circumcenter as the shared vertex; each cell is the
        # union of two infinite half-strips, approximated by a large
        # triangle as in `st_cutter_py`.
        pts = seeds_xy
        cc = _circumcenter(pts[0], pts[1], pts[2])
        # midpoints of each opposing edge
        midpts = [(pts[(i + 1) % 3] + pts[(i + 2) % 3]) / 2.0 for i in range(3)]
        # far points from cc through each midpoint
        radius = np.ptp(pts, axis=0).max() * 10
        far = [cc + (m - cc) / np.linalg.norm(m - cc) * radius for m in midpts]
        # cell for seed i is bounded by far[(i+1)%3], cc, far[(i+2)%3]
        polys = []
        for i in range(3):
            polys.append(Polygon([cc, far[(i + 1) % 3], pts[i], far[(i + 2) % 3]]))
        return polys
    raise ValueError(f"_custom_cutter only handles n<=3, got n={n}")


def _circumcenter(a, b, c):
    ax, ay = a; bx, by = b; cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if d == 0:
        # Collinear — fall back to centroid
        return (a + b + c) / 3.0
    ux = ((ax ** 2 + ay ** 2) * (by - cy)
          + (bx ** 2 + by ** 2) * (cy - ay)
          + (cx ** 2 + cy ** 2) * (ay - by)) / d
    uy = ((ax ** 2 + ay ** 2) * (cx - bx)
          + (bx ** 2 + by ** 2) * (ax - cx)
          + (cx ** 2 + cy ** 2) * (bx - ax)) / d
    return np.array([ux, uy])


def step1b_one_day(
    work_pnt: gpd.GeoDataFrame,
    work_lrg2: gpd.GeoDataFrame,        # conservative-aggregation large polys
    oned: dt.date,
) -> gpd.GeoDataFrame:
    """Build sub-divided polygons for one day.

    Equivalent of step1b_work_v7m.sql STEP 3.  `work_lrg2` is expected
    to contain that day's conservative `work_lrg2` polygons (built the
    same way as work_lrg1 but with pix_dx/dy = fire_dx/dy).

    Returns a work_div-shape GeoDataFrame for the day.
    """
    pnt_mask = work_pnt["acq_date_use"] == oned
    if not pnt_mask.any() or not len(work_lrg2):
        return _empty_work_div()

    day = work_pnt.loc[pnt_mask, [
        "cleanid", "lon", "lat", "fireid2", "alg_agg",
    ]].copy().reset_index(drop=True)

    # ---- 3.1: tbl_close (near-table at SKIM_RADIUS_DEG) ----
    pts = day[["lon", "lat"]].to_numpy()
    tree = KDTree(pts)
    pairs = tree.query_pairs(r=SKIM_RADIUS_DEG, output_type="ndarray")
    # Only pairs sharing the same fireid get to participate in skim.
    if len(pairs):
        fids = day["fireid2"].to_numpy()
        same_fid = fids[pairs[:, 0]] == fids[pairs[:, 1]]
        pairs = pairs[same_fid]
    if len(pairs):
        # invdist
        d = np.linalg.norm(pts[pairs[:, 0]] - pts[pairs[:, 1]], axis=1)
        invdist = np.where(d == 0, 1.0 / (1.0 / 60.0 / 60.0), 1.0 / d)
    else:
        invdist = np.zeros(0)

    # ---- 3.2 + 3.3: pnt2drop, per fireid ----
    todrop_idx: set[int] = set()
    fillers: list[tuple[int, np.ndarray]] = []  # (id_kept, centroid_xy)

    if len(pairs):
        df_p = pd.DataFrame({"l": pairs[:, 0], "r": pairs[:, 1],
                              "inv": invdist, "fid": day["fireid2"].to_numpy()[pairs[:, 0]]})
        for fid, sub in df_p.groupby("fid"):
            todrop, others = _pnt2drop(
                sub[["l", "r"]].to_numpy(), sub["inv"].to_numpy(),
                n_nodes=len(day),
            )
            todrop_idx.update(todrop)
            for d_idx, others_list in zip(todrop, others):
                if others_list:
                    group = [d_idx] + list(others_list)
                    centroid = pts[group].mean(axis=0)
                    fillers.append((min(group), centroid))

    # ---- 3.4: assemble skim point set ----
    keep_mask = ~np.isin(np.arange(len(day)), list(todrop_idx))
    seeds_xy = pts[keep_mask].copy()
    seeds_fid = day["fireid2"].to_numpy()[keep_mask].copy()
    # add fillers
    for kept_min, c in fillers:
        seeds_xy = np.vstack([seeds_xy, c])
        # filler's fireid = kept_min's fireid
        seeds_fid = np.append(seeds_fid, day["fireid2"].iloc[kept_min])

    # ---- 3.5 / 3.6 / 3.8: build a polygon per seed, grouped by fireid ----
    out_rows: list[dict] = []
    lrg2_by_fid = dict(zip(work_lrg2["fireid"], work_lrg2["geom_lrg"]))
    lrg2_aux = work_lrg2.set_index("fireid")[["alg_agg", "acq_date_use"]].to_dict("index")

    for fid in np.unique(seeds_fid):
        sel = seeds_fid == fid
        seeds = seeds_xy[sel]
        parent = lrg2_by_fid.get(int(fid))
        if parent is None or parent.is_empty:
            continue

        # Robust Voronoi handles n=1/2/3/4+ uniformly and uses ghost
        # points sized against the parent's bounding box (not the seed
        # extent), so tightly clustered seeds inside a larger parent
        # still produce cells that fully tile the parent on intersection.
        cells = _voronoi_cells_robust(seeds, parent_bounds=parent.bounds)

        for cell in cells:
            try:
                sub = parent.intersection(cell)
            except Exception:
                continue
            if sub.is_empty:
                continue
            if sub.geom_type == "GeometryCollection":
                sub = unary_union(
                    [g for g in sub.geoms if g.geom_type in ("Polygon", "MultiPolygon")]
                )
            if sub.is_empty:
                continue
            aux = lrg2_aux.get(int(fid), {})
            out_rows.append({
                "fireid": int(fid),
                "geom": sub,
                "acq_date_use": aux.get("acq_date_use", oned),
                "alg_agg": aux.get("alg_agg", pd.NA),
                "cleanids": None,    # filled in below
            })

    if not out_rows:
        return _empty_work_div()

    work_div = gpd.GeoDataFrame(out_rows, geometry="geom", crs="EPSG:4326")
    work_div = work_div.reset_index(drop=True)
    work_div["polyid"] = work_div.index.astype("int64")

    # Drop badly-degenerate geometries (step 3.7 'funky geometries').
    pp = work_div["geom"].map(polsby_popper)
    work_div = work_div[pp > 0]

    # 3.9: back-link cleanids using sjoin.  We use 'intersects' rather
    # than 'within' so points that land exactly on a Voronoi cell edge
    # are still attributed to a polygon.  When a point lies on a shared
    # boundary it may match more than one polygon — we deterministically
    # pick the smallest polyid.
    pnt_day = work_pnt.loc[pnt_mask, ["cleanid", "geom_pnt"]].set_geometry("geom_pnt")
    joined = gpd.sjoin(
        pnt_day, work_div[["polyid", "geom"]].set_geometry("geom"),
        how="inner", predicate="intersects",
    )
    # Drop duplicate (cleanid, polyid) pairs and resolve duplicate cleanids
    # by keeping the smallest polyid.
    joined = (
        joined.sort_values(["cleanid", "polyid"])
              .drop_duplicates(subset=["cleanid"], keep="first")
    )
    cleanids_by_poly = joined.groupby("polyid")["cleanid"].apply(list)
    work_div["cleanids"] = work_div["polyid"].map(cleanids_by_poly)

    # area
    work_div["area_sqkm"] = spheroidal_area_km2_vec(work_div["geom"].values)

    # 3.11: push polyid back to work_pnt
    polyid_by_clean = {}
    for poly_id, lst in zip(work_div["polyid"], work_div["cleanids"]):
        if isinstance(lst, list):
            for c in lst:
                polyid_by_clean[c] = poly_id
    work_pnt["polyid"] = work_pnt["cleanid"].map(polyid_by_clean).astype("Int64")

    # ----------------------------------------------------------------------
    # Diagnostics — surface common polygon-construction problems.  These
    # log warnings with enough context to localize a bad fire; they don't
    # raise.
    # ----------------------------------------------------------------------
    _log_step1b_diagnostics(work_pnt, work_lrg2, work_div, oned)

    return work_div[
        ["polyid", "fireid", "cleanids", "geom", "acq_date_use",
         "area_sqkm", "alg_agg"]
    ].reset_index(drop=True)


def _log_step1b_diagnostics(
    work_pnt: gpd.GeoDataFrame,
    work_lrg2: gpd.GeoDataFrame,
    work_div: gpd.GeoDataFrame,
    oned: dt.date,
) -> None:
    """Log warnings about polygon-coverage and area-conservation problems.

    Three checks:
      1. Detections on `oned` whose work_pnt.polyid is NaN — they fell
         outside every sub-polygon (shouldn't happen if step1a/1b are
         consistent).
      2. fireids in work_lrg2 that produced zero sub-polygons in work_div.
      3. Per-fireid: sum(sub-polygon areas) vs parent (work_lrg2) area.
         A big shortfall means cells were dropped or sub-polygons were
         clipped too aggressively.
    """
    day_mask = work_pnt["acq_date_use"] == oned
    day = work_pnt.loc[day_mask]
    n_day = len(day)

    # (1) detections with no polyid
    n_orphan = int(day["polyid"].isna().sum())
    if n_orphan:
        sample = (day[day["polyid"].isna()]
                  .head(5)[["cleanid", "lon", "lat", "fireid2"]]
                  .to_dict("records"))
        log.warning(
            "step1b %s: %d/%d detections (%.1f%%) not inside any sub-polygon. "
            "First few: %s",
            oned, n_orphan, n_day, 100.0*n_orphan/max(1, n_day), sample,
        )

    # (2) fireids that should have produced polygons but didn't
    fids_in   = set(int(f) for f in work_lrg2["fireid"].unique())
    fids_out  = set(int(f) for f in work_div["fireid"].unique())
    missing   = fids_in - fids_out
    if missing:
        log.warning(
            "step1b %s: %d fireid(s) in work_lrg2 produced zero sub-polygons "
            "(first few: %s)",
            oned, len(missing), list(sorted(missing))[:5],
        )

    # (3) area conservation per fire
    if len(work_div):
        sum_by_fid = work_div.groupby("fireid")["area_sqkm"].sum()
        parent_by_fid = work_lrg2.set_index("fireid")["area_sqkm"]
        joined = pd.concat([parent_by_fid.rename("parent"),
                            sum_by_fid.rename("subs")], axis=1).dropna()
        if len(joined):
            joined["rel_shortfall"] = (joined["parent"] - joined["subs"]) / joined["parent"]
            # >5% shortfall flagged as suspicious; >50% as serious.
            bad = joined[joined["rel_shortfall"] > 0.05]
            if len(bad):
                worst = bad.sort_values("rel_shortfall", ascending=False).head(5)
                log.warning(
                    "step1b %s: %d fire(s) lost >5%% of parent area in "
                    "subdivision (worst: %s)",
                    oned, len(bad),
                    [(int(fid), f"parent={r.parent:.3f} km², "
                                f"subs={r.subs:.3f} km², "
                                f"shortfall={100*r.rel_shortfall:.1f}%")
                     for fid, r in worst.iterrows()],
                )


def _empty_work_div() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "polyid": pd.Series(dtype="int64"),
            "fireid": pd.Series(dtype="int64"),
            "cleanids": pd.Series(dtype="object"),
            "geom": gpd.GeoSeries([], crs="EPSG:4326"),
            "acq_date_use": pd.Series(dtype="object"),
            "area_sqkm": pd.Series(dtype="float64"),
            "alg_agg": pd.Series(dtype="Int64"),
        },
        geometry="geom", crs="EPSG:4326",
    )


# ---------------------------------------------------------------------------
# Conservative aggregation (work_lrg2) — same shape as step1a but with
# pix_dx/dy = fire_dx/dy (the only difference between 1a and 1b's STEP 2).
# ---------------------------------------------------------------------------

def step1b_aggregate_one_day(
    work_pnt: gpd.GeoDataFrame, oned: dt.date
) -> gpd.GeoDataFrame:
    """Conservative version of `step1a_one_day` — produces work_lrg2."""
    mask = work_pnt["acq_date_use"] == oned
    if not mask.any():
        return _empty_work_lrg()
    day = work_pnt.loc[mask].copy().reset_index(drop=False)
    n = len(day)

    fire_size = _fire_size_km(day["instrument"])
    # No scan/track here — pix_dx/dy = fire_dx/dy (conservative)
    dxdy = make_pixel_dxdy(
        lat=day["lat"].to_numpy(),
        fire_size_km=fire_size,
        pix_scan_km=None, pix_track_km=None,
    )
    lon = day["lon"].to_numpy()
    lat = day["lat"].to_numpy()
    geom_sml = shapely.box(
        lon - dxdy["fire_dx"], lat - dxdy["fire_dy"],
        lon + dxdy["fire_dx"], lat + dxdy["fire_dy"],
    )
    pairs = _adjacency_pairs(lon, lat, dxdy["pix_dx"], dxdy["pix_dy"])
    local_fireid, ndetect = _connected_components(n, pairs)
    cleanid = day["cleanid"].to_numpy()
    fireid = np.empty(n, dtype=np.int64)
    for fl in np.unique(local_fireid):
        members = local_fireid == fl
        fireid[members] = cleanid[members].min()
    work_pnt.loc[mask, "fireid2"] = fireid
    work_pnt.loc[mask, "ndetect2"] = ndetect

    large_rows: list[dict] = []
    if len(pairs) > 0:
        local_fid_for_pair = local_fireid[pairs[:, 0]]
        for fl in np.unique(local_fid_for_pair):
            edge_mask = local_fid_for_pair == fl
            edges_for_fid = pairs[edge_mask]
            hulls = []
            for ii, jj in edges_for_fid:
                hulls.append(shapely.convex_hull(shapely.union(geom_sml[ii], geom_sml[jj])))
            polygon = unary_union(hulls)
            polygon = fill_small_holes(polygon, SMALL_HOLE_AREA_DEG2)
            members = local_fireid == fl
            # alg_agg was set on work_pnt by set_alg_agg_from_tree before
            # this function ran; every point in a fireid2 component
            # belongs to one fireid1 component (conservative ⊆ aggressive),
            # so they share the same alg_agg.  Just take the first.
            member_idx = int(np.flatnonzero(members)[0])
            alg = day["alg_agg"].iloc[member_idx] if "alg_agg" in day.columns else pd.NA
            large_rows.append({
                "fireid": int(cleanid[members].min()),
                "geom_lrg": polygon,
                "acq_date_use": oned,
                "ndetect": int(members.sum()),
                "alg_agg": alg,
            })
    lone_mask = ndetect == 1
    for k in np.flatnonzero(lone_mask):
        alg = day["alg_agg"].iloc[int(k)] if "alg_agg" in day.columns else pd.NA
        large_rows.append({
            "fireid": int(cleanid[k]),
            "geom_lrg": geom_sml[k],
            "acq_date_use": oned,
            "ndetect": 1,
            "alg_agg": alg,
        })
    out = gpd.GeoDataFrame(large_rows, geometry="geom_lrg", crs="EPSG:4326")
    if len(out):
        out["area_sqkm"] = spheroidal_area_km2_vec(out["geom_lrg"].values)
    else:
        out["area_sqkm"] = pd.Series(dtype="float64")
    return out
