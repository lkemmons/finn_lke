# finn_py — pure-Python FINN toolkit

A Python reimplementation of the FINN2 wildfire-emissions preprocessor
(no PostgreSQL / PostGIS / psql / IDL dependency), plus companion tools
for calculating emissions, gridding onto regular or unstructured grids,
concatenating files across time, computing regional totals and
plotting maps.

Everything lives on disk: MODIS land-cover and VCF rasters are
Cloud-Optimized GeoTIFFs, per-run "schemas" are GeoParquet, per-fire
and gridded emissions are NetCDF.  Nothing writes to a database.

See `DESIGN.md` for the design decisions and the table-by-table
mapping of what was replaced with what.

## Install

```bash
# system libs (HDF4 — GDAL's HDF4 driver is usually missing)
sudo apt install libhdf4-dev libjpeg-dev       # Debian/Ubuntu
brew install hdf4                              # macOS
conda install -c conda-forge hdf4 pyhdf        # conda

# Python package (editable)
pip install -e .
```

Extra dependencies used only by the emissions / plot / totals scripts:
`xarray`, `netCDF4`, `matplotlib`, `cartopy` — install these too if
you use the corresponding tool.

On casper HPC at NCAR, just use `conda activate npl`.

## What's in the package

### Core preprocessor pipeline (AF → burned-area polygons)

| Script                     | Purpose                                                               |
|----------------------------|-----------------------------------------------------------------------|
| `build_global_rasters.py`  | Build the annual LCT + VCF GeoTIFFs from MODIS HDF tiles              |
| `work_archive.py`          | For multi-day FIRMS archive files of Active Fires                     |
| `run_daily_nrt.py`         | Daily-driver for the near-real-time AF → polygons pipeline            |

### Emissions calculation and gridding

| Script                     | Input                                          | Output                                             |
|----------------------------|------------------------------------------------|----------------------------------------------------|
| `calc_emis_daily.py`       | polygons (`out_*.csv`) + fuel-load, EF (.csv)  | Per-fire NetCDF, `(fire,)` dim                     |
| `grid_emissions.py`        | Per-fire NetCDF (from `calc_emis_daily.py`)    | Gridded NetCDF, `(time, lat, lon)` or `(time, ncol)` |

### Analysis and visualization

| Script                      | What it does                                                   |
|-----------------------------|----------------------------------------------------------------|
| `preview_rasters.py`        | Map annual LCT and VCF GeoTIFF rasters                         |
| `compute_emis_totals.py`    | Sum a gridded file over regions per timestep → CSV + line plot |
| `plot_emissions_scrip.py`   | Polygon-fill map of an unstructured SCRIP-grid emissions file  |

## Usage by task

### 1. Build the annual rasters (do for each year prior to year of fires)

Build:
```bash
python scripts/build_global_rasters.py \
    --lct-root /data/MCD12Q1 --vcf-root /data/MOD44B \
    --regions  /data/regions/All_Countries.shp --region-column "Region_num" \
    --raster-dir ./rasters --year 2002
```

On NCAR Derecho, `build_global.pbs.template` fans this out as a 24-job
array covering 2002–2025. (Not tested)

Expected input layout:
```
<lct-root>/2002/MCD12Q1.A2002001.*.hdf
<lct-root>/2003/MCD12Q1.A2003001.*.hdf
...
<vcf-root>/2002/MOD44B.A2002065.*.hdf
```

The pipeline picks up the year-matching TIFs automatically (via
`modlct_<year>.tif` / `modvcf_<year>.tif`).

### 2. Run AF → burned-area polygons

Single date, NRT (each day is a separate FIRMS file):
```bash
python scripts/run_daily_nrt.py \
    --date 2024061 \
    --af-root  /glade/derecho/scratch/... \
    --raster-dir ./rasters \
    --work-dir ./work --out-dir ./out
```

Single date, from a multi-day FIRMS archive file:
```bash
python scripts/work_archive.py \
    --date 2025061 \
    --af-files /glade/.../fire_archive_M-C61_20241231_20260101.csv \
               /glade/.../fire_archive_SV-C2_20241231_20260101.csv \
    --raster-dir ./rasters --work-dir ./work --out-dir ./out
```

Both scripts write `out_<tag>_<yyyyjjj>_*.gpkg` and `.csv` plus a
`summary_<tag>.txt`.

#### Tropical carryover

The original FINN includes every tropical MODIS detection (|lat| ≤
23.5°) from the previous calendar day, compensating for MODIS swath gaps.
Both (nrt and archive) scripts run this "dup tropics" step by default.

Turn it off with `--no-tropical-carryover`.  In NRT this also skips
the optional yesterday-file lookup that `--tropical-carryover` opts
into.  The two flags are mutually exclusive.

### 3. Calculate per-fire emissions

`calc_emis_daily.py` is a direct port of the IDL
`finn2_9_2_calc_emis_nrt_daily.pro`, with NetCDF output instead of CSV.
All algorithmic logic (genveg assignment, VCF cleanup, fuel-load
lookups, biomass-burned formulas, gas-vs-aerosol unit conversion) is
preserved exactly.

```bash
python scripts/calc_emis_daily.py 2024061 2024 \
    --path-inputs ../finn_inputs \
    --path-in     /glade/.../finn2.9nrt_output \
    --path-out    /glade/.../finn2.9.2nrt_emis \
    --sim-id      NRTmodvrs
```

Output: `FINNv2.9.2_NRTmodvrs_MOZART_YYYYMMDD.nc` with one variable per
species, each shape `(fire,)`.  Gas species carry units `mol/day`,
aerosols `kg/day`, number species `number/day`.  Every species variable
gets a `molecular_weight_g_per_mol` attribute.

### Additional tools

#### Grid emissions onto lat-lon or SCRIP

For per-fire NetCDF (from `calc_emis_daily.py`):
```bash
# 0.1° x 0.1° regular grid, all species
python scripts/grid_emissions.py \
    FINNv2.9.2_NRTmodvrs_MOZART_20240301.nc \
    --grid-resolution 0.1 0.1 --out-dir ./gridded

# Unstructured grid defined by SCRIP file
python scripts/grid_emissions.py \
    FINNv2.9.2_NRTmodvrs_MOZART_20240301.nc \
    --scrip /glade/.../ne30pg3_scrip_170611.nc \
    --grid-label ne30pg3 --out-dir ./gridded_ne30pg3
```

#### Compute regional totals and plot time series

```bash
python scripts/compute_emis_totals.py \
    FINNv2.9.2_NRTmodvrs_MOZART_CO_0.1x0.1deg_20240301.nc \
    --out-dir ./totals
```

Default regions: global, NH, SH, tropics, NH-extratropics,
SH-extratropics (half-open cell-mask convention, so regions don't
double-count cells on their shared boundary).  Add custom regions with
`--region name WEST EAST SOUTH NORTH` (repeatable).  Auto-detects
daily vs hourly cadence via `datesec` presence.  Converts to mass in
kg using the `molecular_weight_g_per_mol` attribute (uniform formula
handles both gases and MW=1 aerosols).

Output: `<input-stem>_totals.csv` (columns per region, values in kg
auto-scaled to Mg/Gg/Tg based on magnitude) and
`<input-stem>_totals.png` (line plot per region).

For SCRIP-grid input pass `--scrip <scrip.nc>` so cell areas can be
computed correctly.

#### Map an emissions file on unstructured grid

```bash
python scripts/plot_emissions_scrip.py \
    FINNv2.9.2_NRTmodvrs_MOZART_CO_ne30_20240301.nc \
    /glade/.../scrip_ne30.nc \
    --out-dir ./maps
```

Renders each SCRIP cell as its own filled polygon (no interpolation).
Handles polygon dateline crossings by duplicating near-seam cells
shifted ±360°.  Log color scale by default, bounded at 6 decades below
the 99.9th-percentile value.  Projections: PlateCarree (default),
Robinson, Mollweide, North/South polar.  Regional zoom via `--extent
WEST EAST SOUTH NORTH`.

## Project layout

```
finn_py/
├── DESIGN.md
├── README.md
├── pyproject.toml
├── finn_py/
│   ├── __init__.py
│   ├── config.py            FinnConfig + hard-wired algorithm constants
│   ├── geometry.py          spheroidal area, hole-fill, Polsby-Popper
│   ├── af_io.py             load shp / csv / zipped FIRMS inputs
│   ├── rasters.py           PyHDF read, sinu→WGS84 mosaicking, COG output
│   ├── step1.py             prep work_pnt + step 1a + step 1b
│   ├── step2.py             zonal joins for LCT / VCF / regnum
│   ├── export.py            CSV + GeoPackage writers
│   ├── pipeline.py          run_nrt(): end-to-end orchestration
│   ├── raster_pipeline.py   run_raster(): mosaic local HDFs (no downloads)
│   └── cli.py               argparse-driven entry points
└── scripts/
    ├── build_global_rasters.py           batch: build per-year global mosaics
    ├── run_daily_nrt.py                  daily NRT driver
    ├── work_archive.py                   multi-day-archive driver
    ├── calc_emis_daily.py                IDL→Python emissions calc
    ├── grid_emissions.py                 per-fire NetCDF → gridded NetCDF
    ├── compute_emis_totals.py            regional totals CSV + plot
    ├── plot_emissions_scrip.py           SCRIP polygon map
    └── preview_rasters.py                plot LCT & VCF GeoTIFFs
```

## Function-by-function map to the original (preprocessor)

| Original SQL/Python                                  | New Python                                                |
|------------------------------------------------------|-----------------------------------------------------------|
| `af_import.main`                                     | `af_io.load_af_files`                                     |
| `step1_prep_v7m.sql` Part 3 (load work_pnt, filters) | `step1.prep`                                              |
| `step1a_work_v7m.sql`                                | `step1.step1a_one_day`                                    |
| `run_vcf.py` tree-cover join                         | `step1.set_alg_agg_from_tree`                             |
| `step1b_work_v7m.sql` STEP 2                         | `step1.step1b_aggregate_one_day`                          |
| `step1b_work_v7m.sql` STEP 3                         | `step1.step1b_one_day`                                    |
| PL/Python `pnt2grp`                                  | `step1._connected_components` (networkx)                  |
| PL/Python `pnt2drop`                                 | `step1._pnt2drop` (networkx)                              |
| PL/Python `st_voronoi_python`                        | `step1._voronoi_cells_robust` (scipy.spatial.Voronoi)     |
| `st_cutter_py` (2/3 point cutter)                    | `step1._custom_cutter`                                    |
| `mkcmd_insert_table_thematic` (LCT majority)         | `step2.thematic_zonal_stats`                              |
| `mkcmd_insert_table_continuous` (VCF mean)           | `step2.continuous_zonal_stats`                            |
| `mkcmd_insert_table_polygons` (regnum centroid)      | `step2.polygons_zonal_join`                               |
| `export_shp.main`                                    | `export.write_output`                                     |
| `rst_import.Importer` (HDF → PostGIS raster)         | `rasters.build_mosaic` (HDF → COG, PyHDF-based)           |
| `modis_tile.py` + `raster.wireframe`                 | `rasters.modis_tile_polygons` / `rasters.tiles_needed`    |
| `tbl_log` + `summarize_log.sql`                      | `pipeline._Logger` (JSON event log + summary text)        |
| `finn2_9_2_calc_emis_nrt_daily.pro` (IDL)            | `scripts/calc_emis_daily.py`                              |
