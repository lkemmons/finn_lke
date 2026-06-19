"""End-to-end NRT pipeline.

Public entry point ``run_nrt(config)`` does the same job as
``work_nrt.main(...)`` in the original code:

  1. load AF input files
  2. build work_pnt
  3. step 1a per day  -> work_lrg1
  4. tree-cover join  -> alg_agg flag on work_lrg1 / work_pnt
  5. step 1b per day  -> work_lrg2, work_div
  6. step 2 per day   -> wide output table
  7. export to CSV + GeoPackage and summarise
"""
from __future__ import annotations

import datetime as dt
import json
import logging
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from . import af_io, export, step1, step2
from .config import FinnConfig, RasterSpec

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Simple run logger — replaces the database's tbl_log + summarize_log.sql
# ---------------------------------------------------------------------------

class _Logger:
    """Append-only event log with row-count deltas, JSON-serialisable."""

    def __init__(self, path: Path | None):
        self.events: list[dict] = []
        self._stack: list[tuple[str, str, int, float]] = []
        self.path = path

    def checkin(self, event: str, table: str, count_before: int, message: str | None = None) -> None:
        self._stack.append((event, table, count_before, time.time()))
        log.info("%s: start (%s n=%d)%s", event, table, count_before,
                 f" :: {message}" if message else "")

    def checkout(self, count_after: int) -> None:
        event, table, before, t0 = self._stack.pop()
        elapsed = time.time() - t0
        self.events.append({
            "event": event, "table": table,
            "n_before": before, "n_after": count_after,
            "n_delta": count_after - before,
            "elapsed_s": round(elapsed, 3),
            "t_finish": dt.datetime.now().isoformat(timespec="seconds"),
        })
        log.info("%s: done (%s n=%d, %+d, %.2fs)",
                 event, table, count_after, count_after - before, elapsed)

    def flush(self) -> None:
        if self.path is None:
            return
        self.path.write_text(json.dumps(self.events, indent=2))


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_nrt(cfg: FinnConfig) -> Path:
    """Run the full pipeline.  Returns the path to the output GeoPackage."""
    cfg.schema_dir.mkdir(parents=True, exist_ok=True)
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    runlog = _Logger(cfg.schema_dir / "tbl_log.json")

    # -----------------------------------------------------------------------
    # 1. Load AF + build work_pnt
    # -----------------------------------------------------------------------
    work_pnt_path = cfg.schema_dir / "work_pnt.parquet"
    if cfg.run_import or not work_pnt_path.exists():
        runlog.checkin("import_af", "af_in", 0,
                       message=",".join(str(p) for p in cfg.af_files))
        af = af_io.load_af_files(
            cfg.af_files,
            tropical_carryover_paths=cfg.tropical_carryover_files or None,
            tropical_carryover_to_date=cfg.tropical_carryover_to_date,
            tropical_lat_bounds=cfg.tropical_lat_bounds,
        )
        runlog.checkout(len(af))

        runlog.checkin("prep_work_pnt", "work_pnt", 0)
        date_range = (cfg.date, cfg.date) if cfg.date is not None else None
        work_pnt = step1.prep(
            af,
            date_definition=cfg.date_definition,
            date_range=date_range,
            filter_persistent_sources=cfg.filter_persistent_sources,
        )
        work_pnt.to_parquet(work_pnt_path)
        runlog.checkout(len(work_pnt))
    else:
        log.info("re-using existing work_pnt at %s (run_import=False)", work_pnt_path)
        work_pnt = gpd.read_parquet(work_pnt_path)

    # -----------------------------------------------------------------------
    # 2. Step 1a per day -> work_lrg1
    # -----------------------------------------------------------------------
    work_lrg1_path = cfg.schema_dir / "work_lrg1.parquet"
    work_div_path = cfg.schema_dir / "work_div.parquet"

    dates = sorted(work_pnt["acq_date_use"].dropna().unique())
    log.info("unique dates in work_pnt: %d (%s..%s)",
             len(dates), dates[0] if dates else "-", dates[-1] if dates else "-")

    if cfg.run_step1:
        # Step 1a (aggressive)
        runlog.checkin("step1a", "work_lrg1", 0)
        lrg1_frames: list[gpd.GeoDataFrame] = []
        for d in dates:
            lrg1_frames.append(step1.step1a_one_day(work_pnt, d))
        work_lrg1 = (
            pd.concat(lrg1_frames, ignore_index=True)
            if lrg1_frames else step1._empty_work_lrg()
        )
        work_lrg1 = gpd.GeoDataFrame(work_lrg1, geometry="geom_lrg",
                                     crs="EPSG:4326")
        runlog.checkout(len(work_lrg1))

        # Tree cover join (step 2-style zonal mean against aggressive lrg)
        runlog.checkin("vcf_treecover", "work_lrg1.alg_agg", len(work_lrg1))
        v_tree = _aggressive_tree_cover(work_lrg1, cfg.rasters)
        step1.set_alg_agg_from_tree(
            work_pnt, work_lrg1, v_tree,
            threshold_pct=cfg.tree_cover_threshold_pct,
        )
        runlog.checkout(int(work_lrg1["alg_agg"].notna().sum()))

        work_lrg1.to_parquet(work_lrg1_path)
        work_pnt.to_parquet(work_pnt_path)  # alg_agg / polyid updated in place

        # Step 1b (conservative aggregation + Voronoi subdivision)
        runlog.checkin("step1b", "work_div", 0)
        div_frames: list[gpd.GeoDataFrame] = []
        lrg2_frames: list[gpd.GeoDataFrame] = []
        for d in dates:
            lrg2 = step1.step1b_aggregate_one_day(work_pnt, d)
            lrg2_frames.append(lrg2)
            div = step1.step1b_one_day(work_pnt, lrg2, d)
            div_frames.append(div)
        work_div = (
            pd.concat(div_frames, ignore_index=True)
            if div_frames else step1._empty_work_div()
        )
        # Re-assign global polyid (concatenation produced per-day local IDs).
        work_div = work_div.reset_index(drop=True)
        work_div["polyid"] = work_div.index.astype("int64")
        work_div = gpd.GeoDataFrame(work_div, geometry="geom", crs="EPSG:4326")
        work_div.to_parquet(work_div_path)

        # Persist work_lrg2 too (the conservative parent polygons step1b
        # subdivides into work_div) so it's available for debugging — same
        # pattern as work_lrg1 above.
        work_lrg2_path = cfg.schema_dir / "work_lrg2.parquet"
        work_lrg2 = (
            pd.concat(lrg2_frames, ignore_index=True)
            if lrg2_frames else step1._empty_work_lrg()
        )
        work_lrg2 = gpd.GeoDataFrame(work_lrg2, geometry="geom_lrg",
                                       crs="EPSG:4326")
        work_lrg2.to_parquet(work_lrg2_path)

        runlog.checkout(len(work_div))
    else:
        log.info("skipping step 1 — reading cached parquet files")
        work_lrg1 = gpd.read_parquet(work_lrg1_path)
        work_div = gpd.read_parquet(work_div_path)

    # -----------------------------------------------------------------------
    # 3. Step 2 per day -> wide output
    # -----------------------------------------------------------------------
    out_path = cfg.schema_dir / "out.parquet"
    if cfg.run_step2:
        runlog.checkin("step2", "out", 0)
        out_frames: list[gpd.GeoDataFrame] = []
        for d in dates:
            day = work_div[work_div["acq_date_use"] == d]
            if not len(day):
                continue
            out_frames.append(step2.zonal_join_one_day(
                day, cfg.rasters, work_pnt=work_pnt,
            ))
        out_gdf = (
            pd.concat(out_frames, ignore_index=True)
            if out_frames else work_div.iloc[0:0].copy()
        )
        out_gdf = gpd.GeoDataFrame(out_gdf, geometry="geom", crs="EPSG:4326")
        out_gdf.to_parquet(out_path)
        runlog.checkout(len(out_gdf))
    else:
        log.info("skipping step 2 — reading cached out.parquet")
        out_gdf = gpd.read_parquet(out_path)

    # -----------------------------------------------------------------------
    # 4. Export
    # -----------------------------------------------------------------------
    basename = f"out_{cfg.tag_af}_modlct_{cfg.year_rst}_modvcf_{cfg.year_rst}_regnum"
    paths = export.write_output(
        out_gdf, cfg.out_dir, basename,
        date_definition=cfg.date_definition,
        write_geopackage=True, write_shapefile=False,
    )

    runlog.flush()
    _write_summary(cfg, runlog, paths)
    return paths.get("gpkg", paths["csv"])


def _aggressive_tree_cover(
    work_lrg1: gpd.GeoDataFrame, rasters: list[RasterSpec]
) -> pd.Series:
    """Per-fireid mean of MOD44B 'tree' band, used to gate alg_agg."""
    vcf = next((r for r in rasters if r.tag.startswith("modvcf_")), None)
    if vcf is None:
        log.warning("no modvcf raster configured; defaulting v_tree=NaN (all fires get alg_agg=2)")
        return pd.Series(np.nan, index=work_lrg1["fireid"], name="v_tree")
    polygons = work_lrg1.set_index("fireid").set_geometry("geom_lrg")
    df = step2.continuous_zonal_stats(
        polygons, vcf.path, var_names=("tree",),
        valid_range=vcf.valid_range,
    )
    df.index.name = "fireid"
    return df["v_tree"]


def _write_summary(cfg: FinnConfig, runlog: _Logger, paths: dict[str, Path]) -> None:
    if cfg.summary_file is None:
        return
    cfg.summary_file.parent.mkdir(parents=True, exist_ok=True)
    with cfg.summary_file.open("a") as fp:
        fp.write(f"\n# FINN-py run {dt.datetime.now().isoformat(timespec='seconds')}\n")
        fp.write(f"tag_af={cfg.tag_af}  year_rst={cfg.year_rst}\n")
        fp.write(f"af_files={','.join(str(p) for p in cfg.af_files)}\n")
        fp.write(f"date={cfg.date}  date_def={cfg.date_definition}\n")
        for ev in runlog.events:
            fp.write(
                f"  {ev['event']:<24s} {ev['table']:<16s} "
                f"{ev['n_before']:>10d} -> {ev['n_after']:>10d}  "
                f"({ev['n_delta']:+d})  {ev['elapsed_s']:>8.2f}s\n"
            )
        fp.write(f"outputs: {paths}\n")
