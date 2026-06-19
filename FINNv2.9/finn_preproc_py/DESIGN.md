# FINN preprocessor — pure-Python rewrite

A reimplementation of `preprocessor/code_bashinterface/work_nrt.py` and
`work_raster.py` (and the `code_anaconda/*.sql` files they depend on)
in pure Python — no PostgreSQL, no PostGIS, no `psql`, no `ogr2ogr`.

## What's mapped to what

| Original                                                  | Pure-Python replacement                                  |
|-----------------------------------------------------------|----------------------------------------------------------|
| PostgreSQL schema `af_<tag>`                              | Working dir `work/<tag>/` with GeoParquet files          |
| PostGIS `raster` schema (`rst_modlct_*`, `rst_modvcf_*`)  | `rasters/<tag>/*.tif` Cloud-Optimized GeoTIFFs           |
| `ogr2ogr ... PG:dbname=finn`                              | `geopandas.read_file` / `to_parquet`                     |
| `psql -f step1_prep_v7m.sql`                              | `finn_py.step1.prep`                                     |
| `psql -f step1a_work_v7m.sql -v oned=...`                 | `finn_py.step1.step1a_one_day`                           |
| `psql -f step1b_work_v7m.sql -v oned=...`                 | `finn_py.step1.step1b_one_day`                           |
| `run_vcf.py`'s on-the-fly SQL                             | `finn_py.step1.join_tree_cover` (between 1a and 1b)      |
| `psql -f step2_work_*_v8b.sql` (generated)                | `finn_py.step2.zonal_join_one_day`                       |
| `ST_Union`, `ST_Intersection`, `ST_Difference`, …         | `shapely.ops.unary_union`, `.intersection`, …            |
| `ST_DWithin` adjacency                                    | `scipy.spatial.cKDTree.query_pairs`                      |
| `ST_Voronoi` / `st_voronoi_python`                        | `scipy.spatial.Voronoi` (kept verbatim from PL/Python)   |
| `pnt2grp` (connected components)                          | `networkx.connected_components` (verbatim)               |
| `pnt2drop` (point skimming)                               | direct port (verbatim)                                   |
| `raster2pgsql ... | psql` (MODIS HDF -> PostGIS)          | `rasters.build_mosaic` (PyHDF -> sinu GeoTIFF -> warp -> COG) |
| `ST_Clip` + `ST_ValueCount` (LCT majority class)          | `rasterio.mask.mask` + `numpy.unique(return_counts=True)`|
| `ST_Clip` + `ST_SummaryStatsAgg` (VCF mean)               | `rasterio.mask.mask` + `numpy.nanmean`                   |
| Region polygon centroid join                              | `geopandas.sjoin` on centroids                           |
| `tbl_log` + `summarize_log.sql`                           | structured-log JSON file per run                         |
| `export_shp.py`                                           | `GeoDataFrame.to_file` + `.to_csv`                       |

## Working-table schemas (GeoParquet)

`work_pnt.parquet` — one row per AF detection that survives initial filtering:

| column            | dtype                | notes                                          |
|-------------------|----------------------|------------------------------------------------|
| cleanid           | int64 (index)        | PK; assigned after all filtering               |
| rawid             | int64                | row-number across input files                  |
| src_file          | int                  | which input file (1, 2, …)                     |
| geom_pnt          | geometry (Point)     | WGS84                                          |
| lon, lat          | float64              | from geometry, denormalized for speed          |
| scan, track       | float64              | from MODIS/VIIRS product                       |
| acq_date_utc      | date                 |                                                |
| acq_time_utc      | string (HHMM)        |                                                |
| acq_date_lst      | date                 | approximate local solar time                   |
| acq_datetime_lst  | timestamp            |                                                |
| acq_date_use      | date                 | whichever the user selected; the only one used |
| instrument        | string               | 'MODIS' or 'VIIRS'                             |
| confident         | bool                 | instrument-specific rule                       |
| anomtype          | int                  | 0–3 from FIRMS 'Type' field (0 if absent)      |
| frp               | float64              |                                                |
| fireid1           | int64                | aggressive cluster id                          |
| ndetect1          | int                  |                                                |
| fireid2           | int64                | conservative cluster id                        |
| ndetect2          | int                  |                                                |
| alg_agg           | int                  | 1 = aggressive, 2 = conservative               |
| polyid            | int64                | back-link to work_div                          |

`work_lrg1.parquet` / `work_lrg2.parquet` / `work_lrg.parquet` — one row per
aggregated burned area (large polygon):

| column        | dtype                  |
|---------------|------------------------|
| fireid        | int64 (index)          |
| geom_lrg      | geometry (Polygon)     |
| acq_date_use  | date                   |
| ndetect       | int                    |
| area_sqkm     | float64                |
| alg_agg       | int                    |
| v_tree        | float64 (work_lrg only)|

`work_div.parquet` — one row per Voronoi sub-divided sub-polygon
(the unit on which step 2 rasters are joined):

| column        | dtype                |
|---------------|----------------------|
| polyid        | int64 (index)        |
| fireid        | int64                |
| cleanids      | list[int64]          | AF detections inside this sub-polygon |
| geom          | geometry (Polygon)   |
| acq_date_use  | date                 |
| area_sqkm     | float64              |
| alg_agg       | int                  |

`out.gpkg` — final deliverable, mirrors `work_div` plus joined raster fields
`v_lct`, `f_lct`, `v_tree`, `v_herb`, `v_bare`, `v_regnum`, `cen_lon`,
`cen_lat`.

## Pipeline order (matches original)

1. **Load AF**           `af_io.load_af_files(fnames)` → DataFrame
2. **Build work_pnt**    `step1.prep(...)` → `work_pnt.parquet`
3. **Step 1a per day**   `step1.step1a_one_day(...)` → appends to `work_lrg1.parquet`
4. **Tree-cover join**   `step1.join_tree_cover(...)` → updates `work_lrg1.alg_agg`
5. **Step 1b per day**   `step1.step1b_one_day(...)` → appends to `work_lrg2`/`work_div`
6. **Step 2 per day**    `step2.zonal_join_one_day(...)` → appends to `out`
7. **Export**            `export.write(...)` → `.gpkg` + `.csv`

## Things that change behaviorally

- **MODIS HDFs are read with PyHDF, not GDAL.** Many GDAL builds
  (conda-forge wheels, recent Homebrew bottles) ship without the HDF4
  driver, so `rasterio.open(hdf_path)` raises. We use PyHDF for raw
  array access and parse `StructMetadata.0` in Python to recover each
  tile's sinusoidal transform, then convert each tile to a plain
  intermediate GeoTIFF before mosaicking. The downstream pipeline never
  touches HDF4 again. System lib: `libhdf4-dev` (Linux) /
  `brew install hdf4` (macOS) / `conda install -c conda-forge hdf4
  pyhdf`.
- **`raster.wireframe` is dropped.** The MODIS tile geometry is computed
  in-memory from the sinusoidal projection definition; no on-disk table.
- **`testpy()` is dropped.** The Python environment is the runtime, not
  PL/Python inside Postgres, so there's nothing to test.
- **`work_div_oned` is no longer temporary** — there is no transaction
  scope. Per-day work is held in-memory and appended to the working
  parquet file at the end of each day.
- **No `tbl_log`** — replaced by a `logging` stream and a JSON summary file
  written by `pipeline._Logger`.
- **Concurrency model.** The original is sequential per-day with the DB
  doing the heavy work. The Python version can use `concurrent.futures`
  to process days in parallel since per-day inputs are independent after
  prep; this is opt-in via `--workers N`.

## Things deliberately preserved

- **Hard-wired constants** (`pixfac = 1.1`, fire size 1.0 km MODIS /
  0.375 km VIIRS, small-hole threshold `(1/240)² deg²`, 0.5 arcmin skim
  distance, 23.5° tropics latitude, 50% tree-cover threshold for
  `alg_agg`) — same numbers, same `_v7m` semantics.
- **Tropics duplication** for MODIS and the consequent first-day drop.
- **Two-algorithm structure**: aggressive then conservative aggregation
  with a tree-cover-gated `alg_agg` flag in between.
- **Per-day independence** of step 1a/1b/2 work.
- **CLI flag set** of `work_nrt.py` — same names, same defaults, so an
  existing wrapper script continues to work.
