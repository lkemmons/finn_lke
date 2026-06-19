"""Run configuration and hard-wired constants.

All the magic numbers from `step1_prep_v7m.sql` and friends live here so
the algorithm behavior is self-documented and easy to override.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Sequence


# ---------------------------------------------------------------------------
# Hard-wired algorithm constants (matching step1*_v7m.sql / run_step2 v8b)
# ---------------------------------------------------------------------------

# Earth radius used in degree<->km conversion in step1a/b prep.  Matches the
# value used inside step1a_work_v7m.sql.  Note this is *not* the WGS84
# semi-major axis; the original chose a slightly smaller value.
EARTH_RADIUS_KM = 6370.997
GREAT_CIRCLE_KM = 2.0 * 3.14159265358979323846 * EARTH_RADIUS_KM

# Fire-pixel sizes (km).  Hard-wired in step1a_work_v7m.sql.
FIRE_SIZE_MODIS_KM = 1.0
FIRE_SIZE_VIIRS_KM = 0.375

# "pixfac" — scaling applied to scan/track when building the pixel envelope
# used for the adjacency test (only in step 1a; step 1b uses fire_size for
# both fire_dx and pix_dx).  Hard-wired in step1a_work_v7m.sql.
PIXFAC = 1.1

# Tropics latitude band where MODIS detections are duplicated to the next
# UTC day.  Hard-wired in step1_prep_v7m.sql.
TROPICS_LAT_DEG = 23.5

# Skim radius for step 1b's pnt2drop (deg ~ 0.5 arcmin ~ 0.93 km).
SKIM_RADIUS_DEG = 0.5 / 60.0

# Small-hole threshold for step 1a hole-filling (deg^2).  Equivalent to
# roughly a 15-arcsec square (~0.5 km on a side at equator).
SMALL_HOLE_AREA_DEG2 = (1.0 / 240.0) ** 2

# Tree-cover threshold (%) above which a fire is aggregated aggressively.
# Set by run_vcf.py's algorithm_merge_aggressive_threshold_tree_cover kwarg.
TREE_COVER_AGGRESSIVE_THRESHOLD_PCT = 50.0

# Confidence rules per instrument.  Originally encoded in tbl_flddefs.
MODIS_CONFIDENCE_MIN = 20      # confidence >= 20
VIIRS_CONFIDENCE_DROP = "l"    # drop 'l' (low) confidence

# Anomaly types that 'filter_persistent_sources' drops (volcano / other).
PERSISTENT_ANOMTYPES = (1, 2)


# ---------------------------------------------------------------------------
# Run configuration
# ---------------------------------------------------------------------------

DateDefinition = Literal["UTC", "LST"]


@dataclass
class RasterSpec:
    """Defines one raster join target for step 2.

    Mirrors the dicts in `work_common.sec1_user_config`'s `rasters` list.
    """
    tag: str
    kind: Literal["thematic", "continuous", "polygons", "input"]
    # For 'thematic' and 'polygons': single attribute name.
    # For 'continuous': list of band names (MODIS VCF = ['tree','herb','bare']).
    # For 'polygons': also need `variable_in` (source column on the polygon).
    variable: str | None = None
    variables: Sequence[str] | None = None
    variable_in: str | None = None
    path: Path | None = None  # path to the COG/shapefile on disk
    # Pixels outside [valid_min, valid_max] are excluded from zonal stats
    # *in addition to* the GeoTIFF nodata tag.  Used for raster products
    # that encode special states with in-range-uint8 sentinels — e.g.
    # MOD44B VCF: 200 = water, 251 = cloud, 252 = snow/ice, 253 = fill.
    # None = trust the GeoTIFF nodata tag and use everything else.
    valid_range: tuple[float, float] | None = None


def default_rasters(year_rst: int, raster_dir: Path) -> list[RasterSpec]:
    """The standard 3-raster recipe used by work_nrt.py."""
    tag_lct = f"modlct_{year_rst}"
    tag_vcf = f"modvcf_{year_rst}"
    return [
        RasterSpec(
            tag=tag_lct, kind="thematic", variable="lct",
            path=raster_dir / f"{tag_lct}.tif",
        ),
        RasterSpec(
            tag=tag_vcf, kind="continuous", variables=("tree", "herb", "bare"),
            path=raster_dir / f"{tag_vcf}.tif",
            # MOD44B encodes water=200, cloud=251, snow=252, fill=253 in the
            # same uint8 as the 0-100 percent values.  Cap at 100 so those
            # don't pollute the zonal means.
            valid_range=(0, 100),
        ),
        RasterSpec(
            tag="regnum", kind="polygons", variable="regnum",
            variable_in="Region_num",
            path=raster_dir / "regnum.gpkg",
        ),
        # FRP isn't a raster — it's averaged from work_pnt.frp into each
        # output polygon.  Produces the v_frp column to match the
        # original FINN's output CSV.
        RasterSpec(
            tag="frp", kind="input", variable="frp", variable_in="frp",
        ),
    ]


@dataclass
class FinnConfig:
    """Run-level configuration.  Equivalent of `sec1_user_config`.

    The pipeline runs for a single calendar day.  ``date`` selects which
    day's AF detections to process; the working dir, output filenames,
    and the date filter all key off it.
    """

    tag_af: str
    year_rst: int
    af_files: list[Path]
    out_dir: Path

    # Single processing date (Python datetime.date).  Required.
    date: object | None = None  # datetime.date

    # Working / output directories.
    work_dir: Path = field(default_factory=lambda: Path("./work"))
    raster_dir: Path = field(default_factory=lambda: Path("./rasters"))

    # Algorithm options.
    filter_persistent_sources: bool = True
    date_definition: DateDefinition = "UTC"
    tree_cover_threshold_pct: float = TREE_COVER_AGGRESSIVE_THRESHOLD_PCT

    # Tropical carryover (replicates the original FINN's previous-day
    # MODIS duplication in the tropics that compensates for MODIS swath
    # gaps).  Set ``tropical_carryover_files`` to the list of previous-
    # day MODIS files and ``tropical_carryover_to_date`` to the date
    # they should be re-labelled as.  When ``run_daily_nrt.py`` is
    # invoked with ``--tropical-carryover``, it auto-fills these from
    # yesterday's filename.  Defaults to ±23.5° to match the original
    # FINN's step1_prep_v7m.sql duplication threshold exactly.
    tropical_carryover_files: list[Path] = field(default_factory=list)
    tropical_carryover_to_date: object | None = None  # datetime.date | None
    tropical_lat_bounds: tuple[float, float] = (-23.5, 23.5)

    # Pipeline toggles.
    run_import: bool = True
    run_step1: bool = True
    run_step2: bool = True

    # Parallelism.  0 or 1 = sequential.
    workers: int = 1

    # Summary file (appended to).
    summary_file: Path | None = None

    @property
    def schema_dir(self) -> Path:
        """Per-run working directory; replacement for PG schema `af_<tag>`."""
        return self.work_dir / f"af_{self.tag_af}"

    @property
    def rasters(self) -> list[RasterSpec]:
        return default_rasters(self.year_rst, self.raster_dir)
