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

Extra deps used only by the emissions / plot / totals scripts:
`xarray`, `netCDF4`, `matplotlib`, `cartopy` — install these too if
you use the corresponding tool.

On casper HPC at NCAR, just use `conda activate npl`.

## What's in the package

### Core preprocessor pipeline (AF → burned-area polygons)

| Script                     | Purpose                                                               |
|----------------------------|-----------------------------------------------------------------------|
| `work_raster.py`           | Build the annual LCT + VCF GeoTIFFs from MODIS HDF tiles (regional)   |
| `build_global_rasters.py`  | Same but global (streaming warp, bounded RAM); PBS job-array template |
| `work_nrt.py`              | Original single-shot AF → polygons CLI shim                           |
| `run_daily_nrt.py`         | Daily-driver for the AF → polygons pipeline (NRT)                     |
| `work_archive.py`          | Same, for multi-day FIRMS archive files                               |

### Emissions calculation and gridding

| Script                     | Input                                          | Output                                             |
|----------------------------|------------------------------------------------|----------------------------------------------------|
| `calc_emis_daily.py`       | polygons (`out_*.gpkg`) + fuel-load / EF CSVs  | Per-fire NetCDF, `(fire,)` dim                     |
| `grid_emissions.py`        | Per-fire NetCDF (from `calc_emis_daily.py`)    | Gridded NetCDF, `(time, lat, lon)` or `(time, ncol)` |
| `grid_txt_emissions.py`    | Legacy FINN text CSV (for FINNv2.7-hourly)     | Gridded NetCDF, `(time, lat, lon)` or `(time, ncol)` |
| `txt_to_perfire_nc.py`     | Legacy FINN text CSV (hourly or daily)         | Per-fire NetCDF, `(polyid, hour)` or `(polyid,)`   |

### File combination

| Script                        | What it does                                                        |
|-------------------------------|---------------------------------------------------------------------|
| `concat_finn_daily_hourly.py` | Prepend 24 hourly copies of one daily-average day to an hourly file |
| `append_daily_to_hourly.py`   | Append 24 hourly copies of one daily-average day to an hourly file  |

Together they let you extend an hourly-cadence file by one day at
either end, using data from a daily-average source file — while
keeping uniform hourly cadence throughout (developed for FINNv2.7-hourly).

### Analysis and visualization

| Script                      | What it does                                                   |
|-----------------------------|----------------------------------------------------------------|
| `compute_emis_totals.py`    | Sum a gridded file over regions per timestep → CSV + line plot |
| `plot_emissions_scrip.py`   | Polygon-fill map of an unstructured SCRIP-grid emissions file  |

### Supporting data

| File                              | Purpose                                              |
|-----------------------------------|------------------------------------------------------|
| `finn_text_species_mw.csv`        | Species → type (gas/aerosol/number) + MW lookup      |

## Usage by task

### 1. Build the annual rasters (do for each year of fires)

Regional build:
```bash
finn-raster -y 2020 \
    --lct      /data/MCD12Q1/2020/MCD12Q1.A2020001.*.hdf \
    --vcf      /data/MOD44B/2020/MOD44B.A2020065.*.hdf \
    --regions  /data/regions/All_Countries.shp \
    --raster_dir ./rasters
# produces:  ./rasters/modlct_2020.tif
#            ./rasters/modvcf_2020.tif
#            ./rasters/regnum.gpkg
```

Global multi-year build (streaming warp; peak RAM ~ `--warp-mem-mb`,
default 2 GB):
```bash
python scripts/build_global_rasters.py \
    --lct-root /data/MCD12Q1 --vcf-root /data/MOD44B \
    --regions  /data/regions/All_Countries.shp \
    --raster-dir ./rasters --years 2002 2025
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

### 4. Grid emissions onto lat-lon or SCRIP

For per-fire NetCDF (from `calc_emis_daily.py`):
```bash
# 0.1° x 0.1° regular grid, all species
python scripts/grid_emissions.py \
    FINNv2.9.2_NRTmodvrs_MOZART_20240301.nc \
    --grid-resolution 0.1 0.1 --out-dir ./gridded

# MPAS SCRIP grid
python scripts/grid_emissions.py \
    FINNv2.9.2_NRTmodvrs_MOZART_20240301.nc \
    --scrip /glade/.../scrip_mxc.nc \
    --grid-label mxc --out-dir ./gridded_mxc
```

For legacy FINN text CSVs (needs `--mw-table` for unit annotation):
```bash
python scripts/grid_txt_emissions.py \
    FINN_v2.7_MOZART_hourly_2024033.txt \
    --mw-table scripts/finn_text_species_mw.csv \
    --grid-resolution 0.1 0.1 --out-dir ./gridded
```

Both produce one NetCDF per species with **CAM/CESM emission-input
conventions**: single unlimited `time` dim, plus `date(time)` (int32
YYYYMMDD) and `datesec(time)` (int32 seconds-of-day) companion
variables.  Units default to `molecules/cm2/s` for gases and aerosols,
`number/cm2/s` for particle-count species.  Multi-day concatenation is
a plain `xr.concat` along `time`.

### 5. Extend an hourly file with one more day

For FINNv2.7-hourly: 
Both tools pick one day from a daily-average source file, replicate
its spatial pattern across 24 hourly slots, rebase onto the hourly
file's epoch, then concat.  Result: uniform hourly cadence throughout.

Prepend the day before the hourly file's first date:
```bash
python scripts/concat_finn_daily_hourly.py \
    --daily-dir  /glade/.../daily_avg    \
    --hourly-dir /glade/.../hourly_avg   \
    --out-dir    /glade/.../concat
# default seam date = day BEFORE hourly's first date
```

Append the day after the hourly file's last date:
```bash
python scripts/append_daily_to_hourly.py \
    --hourly-dir /glade/.../concat       \
    --daily-dir  /glade/.../daily_avg    \
    --out-dir    /glade/.../extended
# default seam date = day AFTER hourly's last date
```

Both handle:
- **Cross-epoch time units** (e.g. daily `days since 1950-01-01` +
  hourly `days since 1970-01-01`): daily is rebased onto hourly's
  epoch automatically.
- **Cross-calendar** (gregorian / standard / proleptic_gregorian
  treated as equivalent; other calendars refused with a clear error).
- **Generic emission variable names**: some FINN files call the
  emission variable `fire` with the actual species in a `long_name`
  or `map` attribute; species pairing across a directory of `fire`
  variables works via those attributes.
- **Coord/data-var inconsistencies**: `lat` and `lon` are promoted to
  coordinates before concat, so files that declare
  `fire:coordinates = "lat lon"` mix cleanly with files that don't.

### 6. Convert legacy FINN text CSVs to per-fire NetCDF

Some workflows need to keep the per-fire structure rather than grid
onto a raster.  `txt_to_perfire_nc.py` writes one entry per unique
POLYID.

```bash
# Hourly file (auto-detected from HOUR-column presence) → (polyid, hour)
python scripts/txt_to_perfire_nc.py \
    FINN_v2.7_hourly_2024033.txt \
    --mw-table scripts/finn_text_species_mw.csv \
    --out-dir ./perfire

# Daily file (no HOUR column) → (polyid,)
python scripts/txt_to_perfire_nc.py \
    FINN_v2.7_daily_2024033.txt \
    --out-dir ./perfire

# Force daily even if HOUR present (sums across the 24 hours per fire)
python scripts/txt_to_perfire_nc.py hourly_file.txt --daily --out-dir ./perfire
```

The **HOUR column's values are ignored** (unreliable in some FINN
outputs) and the column is never written to output.  The actual UTC
hour is derived from `datetimeUTC`.  The presence of the HOUR column
is used only as a hint that the file is hourly.

### 7. Compute regional totals and plot time series

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

### 8. Map an unstructured (SCRIP) grid file

```bash
python scripts/plot_emissions_scrip.py \
    FINNv2.9.2_NRTmodvrs_MOZART_CO_mxc_20240301.nc \
    /glade/.../scrip_mxc.nc \
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
    ├── work_nrt.py                       drop-in CLI shim (single run)
    ├── work_raster.py                    drop-in CLI shim
    ├── build_global_rasters.py           batch: build per-year global mosaics
    ├── build_global.pbs.template         Derecho PBS 24-job array template
    ├── run_daily_nrt.py                  daily NRT driver
    ├── run_daily_nrt.pbs.template        Derecho PBS template
    ├── work_archive.py                   multi-day-archive driver
    ├── work_archive.pbs.template         Derecho PBS template
    ├── calc_emis_daily.py                IDL→Python emissions calc
    ├── grid_emissions.py                 per-fire NetCDF → gridded NetCDF
    ├── grid_txt_emissions.py             legacy text CSV → gridded NetCDF
    ├── txt_to_perfire_nc.py              legacy text CSV → per-fire NetCDF
    ├── concat_finn_daily_hourly.py       prepend one daily day to hourly file
    ├── append_daily_to_hourly.py         append one daily day to hourly file
    ├── compute_emis_totals.py            regional totals CSV + plot
    ├── plot_emissions_scrip.py           SCRIP polygon map
    └── finn_text_species_mw.csv          species → type + MW lookup
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
