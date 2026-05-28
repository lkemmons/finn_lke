# FINN Preprocessor — Detailed Description of `work_nrt.py` and the `code_anaconda` Pipeline

## 1. What the preprocessor does

The Fire INventory from NCAR (FINN) preprocessor turns satellite active-fire (AF) point detections from MODIS and VIIRS into a daily, georeferenced layer of burned-area polygons, each tagged with land-cover, vegetation-cover-fraction (VCF), and region attributes. Those polygons are then handed downstream to the FINN emissions code, which multiplies area by land-cover- and region-specific emission factors to produce trace-gas and aerosol emissions for atmospheric chemistry models.

`preprocessor/code_bashinterface/work_nrt.py` is the top-level driver script. It is the bash-line equivalent of the legacy Jupyter notebook workflow, intended for daily near-real-time (NRT) operation. It does not do any geometry itself; it imports a handful of Python modules from `preprocessor/code_anaconda/`, and each of those modules in turn launches `psql` against a PostgreSQL+PostGIS database where almost all of the heavy lifting happens. The Python layer is essentially a thin orchestrator that:

1. parameterizes SQL files with the user's tag and date range,
2. shells out to `psql -f <script>.sql -v tag=... -v oned=...`,
3. writes per-day log files, and
4. at the end, exports the result out of the database as CSV and shapefile.

The database has two main groups of objects: a `raster` schema (shared) that holds the imported MODIS land-cover and VCF rasters plus a MODIS tile wireframe, and a per-run `af_<tag>` schema that holds the working tables built from a particular set of AF inputs.

## 2. Top-level structure of `work_nrt.py`

The file is short (~275 lines) and contains four functions plus an `argparse` block.

### 2.1 Imports and module wiring

```python
sys.path = sys.path + ['../code_anaconda']
import af_import
import run_step1a as runner_step1a
import run_step1b as runner_step1b
import run_vcf   as runner_vcf
import run_step2 as runner_step2
import export_shp
import run_extra
import work_common as common
```

`work_common` is sibling code in `code_bashinterface/`; it owns the database connection environment variables and configuration defaults. Every other import is from `code_anaconda/`. Each "runner" module is a thin shell around one or more `*.sql` files that lives next to it.

### 2.2 Command-line interface

`argparse` exposes the following:

| Flag | Meaning |
|---|---|
| `-t / --tag_af` (required) | A short string identifying this run. It becomes the suffix of the PostgreSQL schema `af_<tag>` and the prefix of every output filename. |
| `-y / --year_rst` (required) | Calendar year of the MODIS raster collections to use (MCD12Q1 land-cover and MOD44B VCF). Translated to `tag_lct = 'modlct_<year>'` and `tag_vcf = 'modvcf_<year>'`. |
| `-o / --out_dir` (required) | Directory for final CSV/SHP and summary file. |
| `-fd / --first_day`, `-ld / --last_day` | Date range as `YYYYJJJ` (year + day-of-year). Optional; if absent, all days present in the AF input are processed. |
| `-s / --summary_file` | A text file the script appends per-run diagnostics to. |
| `--run_import`, `--run_step1`, `--run_step2` | Boolean stage toggles (default true) so a user can re-run only part of the pipeline. |
| `af_fnames` (positional, 1+) | One or more AF input files — typically MODIS or VIIRS shapefiles delivered by FIRMS as `fire_archive_*.shp` or NRT files like `MODIS_C6_*.shp` / `VNP14IMGTDL_NRT_*.shp`. Zip files and CSV/TXT are also accepted. |

The custom `str2bool` helper lets each toggle accept the usual `yes/true/y/1` / `no/false/n/0` spellings.

### 2.3 `main(...)`

`main` is the entry point both for the CLI and for callers that import `work_nrt` directly. It:

1. Creates `out_dir` if needed.
2. Converts `first_day`/`last_day` from `YYYYJJJ` integers into `datetime.date` objects via `datetime.datetime.strptime(..., '%Y%j').date()`.
3. Calls `common.sec1_user_config(tag_af, af_fnames, year_rst)` which returns a dict of "hard-wired" defaults (raster definitions, `filter_persistent_sources = True`, the date-definition mode, etc.) and **dumps it into module globals** via `globals().update(user_config)`. After this line, names like `tag_lct`, `tag_vcf`, `tag_regnum`, `rasters`, `filter_persistent_sources`, and `date_definition` become module-level globals visible to `sec3_import_af`, `sec6_process_activefire`, and `sec7_export_output`. (This is unusual Python style — those functions take their parameters via globals rather than arguments — but it is how the original notebook code factored cleanly into helper sections.)
4. Calls `common.sec2_check_environment(out)`, which prints versions of PostgreSQL/PostGIS, runs `testpy.sql` to confirm the PL/Python extension can import NumPy + networkx, and calls `rst_import.prep_modis_tile()` to make sure the `raster.wireframe` table (MODIS sinusoidal tile boundaries) is populated.
5. If `run_import` is true, calls `sec3_import_af(out)`.
6. Calls `sec6_process_activefire(first_day, last_day, run_step1, run_step2)`.
7. Calls `sec7_export_output(out_dir, summary_file)`.

### 2.4 `sec3_import_af(out)` — getting AF data into the DB

This function takes the user's `af_fnames` list and turns each entry into an actual `.shp` file that can be loaded. The logic is mostly pre-flight bookkeeping because FIRMS data arrives in different shapes:

- It walks each filename and inspects the extension.
- If the file already exists and is `.shp`/`.csv`/`.txt`, it is used directly.
- If it is `.zip`, it is unzipped into the same directory and the corresponding `.shp` is selected.
- If it is `.shp` but the file is **missing**, three recovery strategies are tried:
  - If it is the named sample `fire_archive_M6_28864.shp`, the corresponding global archive is `wget`'d from S3.
  - Otherwise the function tries to match `fire_archive_<arcname>.shp` and look for a sibling `DL_FIRE_<arcname>.zip` to expand.
  - Otherwise it tries the NRT pattern `(MODIS_C6|VNP14IMGTDL_NRT)_*.shp` and looks for the corresponding `.zip`.
- Anything that does not match any of those raises a `RuntimeError`.

Once every entry has been resolved to a real `.shp`, the function calls `af_import.main(tag_af, af_fnames)` which is what actually loads the data into the database, then prints a row count for each loaded table by shelling out `psql -c 'select count(*) from "af_<tag>".af_in_<i>;'`.

There is a TODO comment noting that this stage is destructive — `af_import.main` unconditionally drops and recreates the schema — and that a future version should list existing tables and ask the user before clobbering.

### 2.5 `sec6_process_activefire(...)` — the heart of the pipeline

This is the orchestration of the daily geometry build. It performs three or four ordered stages, each delegated to a `run_*` module:

```
runner_step1a.main(tag_af, ..., first_day, last_day, date_definition)
runner_vcf.main(tag_af, rasters)
runner_step1b.main(tag_af, ..., first_day, last_day, date_definition)
runner_step2.main(tag_af, rasters, first_day, last_day)
```

The pre-amble validates `first_day` and `last_day` against the AF data range by calling `af_import.get_dates('af_' + tag_af, combined=True)`, raising if the requested range is fully outside what the input files contain. It deliberately keeps one extra day before `first_day` because of MODIS-tropics carry-over (see step 1 prep below).

`run_step1` and `run_step2` are honored as Boolean toggles. There is no `--run_vcf` toggle — vcf is bundled into the `run_step1` toggle in this script, because `run_vcf.main` is called between the step-1a and step-1b passes.

### 2.6 `sec7_export_output(out_dir, summary_file)`

After all processing, this function calls `export_shp.main(...)` to dump the final output table out to both CSV and `.shp`, and `run_extra.summarize_log(...)` + `run_extra.db_use_af(...)` to write a log summary and a per-table disk-usage table to the summary file. Output base name is `out_<tag_af>_<tag_lct>_<tag_vcf>_<tag_regnum>.shp` — i.e., the filename carries the run tag and the year of the raster datasets and the region tag.

## 3. `work_common.py` — environment configuration

Although `work_common.py` is not in `code_anaconda/`, every `code_anaconda` module reads the environment variables it sets, so it is functionally part of the pipeline.

Driven by the `FINN_DRIVER` env var (`use_native`, `use_docker`, or `from_inside_docker`), `work_common`:

- Sets `PGDATABASE`, `PGUSER`, `PGPASSWORD`, `PGPORT`, and `PGHOST` so every later `psql` / `psycopg2` call connects to the right server.
- Sets `raster_download_rootdir` for the MODIS HDF cache.
- Prompts (or pulls from env) `EARTHDATAUSER` and `EARTHDATAPW` credentials (NASA Earthdata, required to fetch MODIS).
- Appends `/usr/pgsql-11/bin` to `PATH`.
- Adds `../code_anaconda` to `sys.path` and imports `rst_import`.

`sec1_user_config(tag_af, af_fnames, year_rst)` returns a dict-of-locals that includes:

- `filter_persistent_sources = True` — drops AF "Type 1" (volcano) and "Type 2" (other static thermal anomalies) when the Type field is present. Has no effect on NRT data because that flag is only in the archive product.
- `date_definition` — `'LST'` or `'UTC'`. Decides whether daily slicing happens by approximate Local Solar Time (UTC + round(longitude/15) hours) or by UTC.
- `tag_lct`, `tag_vcf`, `tag_regnum` — the raster-dataset identifiers used to name PostgreSQL tables in the `raster` schema.
- `rasters` — a list of three dicts describing the rasters to be joined against polygons: a `thematic` raster (MODIS LCT), a `continuous` raster with three bands `tree/herb/bare` (MODIS VCF), and a `polygons` layer (`regnum`, the FINN region-number polygon).
- Several `wipe_*` flags that default to False (the NRT script does not clean up between runs).

## 4. `code_anaconda` — the modules `work_nrt.py` calls

What follows is what each imported module does, ordered the way `work_nrt.py` calls them.

### 4.1 `af_import.py` — load AF points into PostGIS

`main(tag, fnames)`:

1. Drops and recreates schema `af_<tag>` with `psql -c "DROP SCHEMA IF EXISTS ... CASCADE; CREATE SCHEMA ...;"`.
2. For each filename, builds an `ogr2ogr` command that loads the data straight into PostgreSQL as table `af_<tag>.af_in_<i>` (i = 1, 2, ...). Loader options used:
   - `-f PostgreSQL -overwrite`
   - `-lco SPATIAL_INDEX=GIST` (or `=YES` for GDAL < 2.4)
   - `-lco GEOMETRY_NAME=geom`, `-lco FID=gid`, `-lco SCHEMA=af_<tag>`
3. For `.csv`/`.txt` inputs there are no Point features natively, so `mk_vrt(fname)` first writes an `OGRVRTDataSource` wrapper (`fname.vrt`) that defines a `wkbPoint` layer using `longitude` and `latitude` as `PointFromColumns`, declares the WGS84 SRS, and types all the FIRMS columns (`scan`, `track`, `acq_date`, `acq_time`, `satellite`, `confidence`, `version`, `brightness`, `bright_t31`, `bright_ti4`, `bright_ti5`, `frp`, `daynight`). `ogr2ogr` then reads the VRT as if it were a shapefile.

The module also has utility query helpers used elsewhere:

- `gdal_vernum_sys()` parses `gdal-config --version` to decide which `-lco SPATIAL_INDEX` syntax to emit.
- `check_raster_contains_fire(rst, fire)` returns a dict counting how many fire points fall inside a raster skeleton's coverage.
- `get_tiles_needed(schema, combined)` SELECTs from each `af_in_<i>` and `raster.wireframe` to find which MODIS sinusoidal tiles the fires fall into (used by `rst_import` to decide which HDF tiles to download).
- `get_lnglat(schema, combined)` returns `(lon, lat)` arrays from each `af_in_<i>`.
- `get_dates(schema, combined)` returns the distinct `acq_date` values per AF table — this is what `work_nrt.sec6_process_activefire` uses to validate `first_day` / `last_day`.

### 4.2 `run_step1a.py` — aggressive aggregation pass

`run_step1a.main(tag, first_day, last_day, ver='v7m', ...)`:

1. Sets `PGCLIENTENCODING='utf-8'` (mac/docker workaround).
2. **Prep**: runs `psql -f step1_prep_<ver>.sql -v tag=... -v filter_persistent_sources=... -v date_range=... -v date_definition=...`. Re-tries up to 3 times with back-off if it fails. Output to `log.step1a.prep`.
3. If `first_day`/`last_day` are not supplied, `get_first_last_day(tag)` queries `af_<tag>.work_pnt` for the min/max `acq_date_lst`, with a special case: if any MODIS detections exist between ±30° latitude, the *first* day is incremented by one because of the "tropics duplication" performed in prep (every tropical MODIS detection is duplicated to the next UTC day, so the first day in the table is incomplete and must be discarded).
4. **Work**: iterates over each date in `[first_day, last_day]` and runs `psql -f step1a_work_<ver>.sql -v tag=... -v oned='YYYY-MM-DD'`, writing `log.step1a.o<YYYYMMDD>` for each day.

#### What `step1_prep_v7m.sql` does (1153 lines)

This is the big preparatory script — it both builds the schema **and** loads all the points into the master `work_pnt` working table.

**Part 1: Tables.** Drops/creates:
- `work_pnt` — one row per AF detection, with raw + cleaned ids, geometry (point, "small" fire square, "pixel" actual MODIS/VIIRS footprint), date in UTC and LST, instrument, confidence boolean, anomaly type, FRP, and two pairs of `(fireid<N>, ndetect<N>)` for the two clustering algorithms.
- `work_lrg1` — one row per "large" aggregated fire polygon (aggressive grouping).
- `tbl_flddefs` — instrument-specific column definitions (the rule that MODIS confidence is "> 20" but VIIRS confidence is "!= 'l'").
- `tbl_options` — the run-time options passed in from Python (`filter_persistent_sources`, `date_range`, `date_definition`).
- `tbl_log` — a structured log table populated throughout by `log_checkin(event, table, nrec)` / `log_checkout(id, nrec)` functions that wrap each major mutation and record how many rows changed and how long it took.

**Part 2: Functions and Types.** Several PL/pgSQL and PL/Python helper objects:
- `log_checkin`, `log_checkout`, `log_purge` for the logging pattern above.
- `testpy()` (PL/Python) that imports `sys`, `numpy`, and `networkx` and returns their versions — `sec2_check_environment` runs it to make sure the database has the right Python environment.
- Type `p2grp(fireid int, lhs int, rhs int, ndetect int)` and function `pnt2grp(lhs int[], rhs int[]) RETURNS SETOF p2grp` (PL/Python). Given an edge list between adjacent detections, it builds a `networkx.Graph`, extracts connected components, and emits one row per (edge, connected-component-min-id, component-size). This is how multiple near-by points get a shared `fireid`.
- Type `p2drp(id int, others int[])` and function `pnt2drop(lhs, rhs, invdist)` (PL/Python). Used by step 1b. Iteratively eliminates the highest-density point from each connected component until no edges remain, recording which points get dropped and which clusters of tied-score points get replaced by their centroid. This is the "skim" used to thin clusters before Voronoi tessellation.
- Type `list_of_polygon_coords(x float[], y float[], pos int)` and function `st_voronoi_python(x[], y[])` (PL/Python via `scipy.spatial.Voronoi`). Computes Voronoi cells in Python and returns finite polygon coordinates back to SQL — used because PostGIS's own `ST_VoronoiPolygons` was unsatisfactory for this use case at the time.
- `st_cutter_python` / `st_cutter_py` — special-case "cutter" geometries when there are only 2 or 3 points (Voronoi requires ≥ 4 in 2D), implementing circumcircle-based hand-built triangulations and obtuse-triangle handling.
- `st_polsbypopper(geom)` and `st_polsbypopper(geom, use_spheroid)` — the Polsby-Popper compactness measure (1 = circle, 0.78 = square, →0 = elongated). Used later to filter "funky" geometries.
- `get_acq_datetime_lst(date, time, lon)` — converts UTC date+time to approximate local-solar-time by adding `round(lon/15)` hours.
- `time_to_char(time)` — normalizes acquisition time strings (strips colons from `HH:MM`).
- `get_instrument(satellite)` — maps the FIRMS satellite letter to `MODIS` or `VIIRS` (`T`/`A` → MODIS, `N` → VIIRS).

**Part 3: Processing.**
- Builds `af_ins`, the list of `af_in_<i>` tables that were just loaded, recording whether each has the "Type" column. If none do, `filter_persistent_sources` is forced to false (you cannot filter what is not there).
- A `do $$ ... $$` block loops over each `af_in_<i>` table and `EXECUTE`s a dynamic `INSERT INTO work_pnt SELECT ... FROM af_in_<i>`. Each row gets converted to common columns: `rawid` (a monotonically-numbered id across all input tables), the geometry, lon/lat, scan/track, both UTC and LST dates, instrument, a `confident` boolean from the instrument-specific rule, the anomaly `type` (or 0 if absent), and `frp`. `acq_date_use` is either UTC or LST depending on the run option.
- "Dup tropics": every MODIS detection with `abs(lat) <= 23.5` is duplicated with all dates shifted forward by one day. This makes day-N processing also see tropical day-(N-1) detections, supporting MODIS's two daily overpasses' incomplete day coverage. This is precisely why `run_step1a.get_first_last_day` skips the first day.
- "Drop by date": deletes rows whose `acq_date_use` falls outside the `date_range` (as a `daterange`), if a range was supplied.
- "Drop low confidence": deletes rows where `confident = false`.
- "Drop persistent": if the option is on, deletes rows where `anomtype = 1` (volcano) or `2` (other persistent).
- Adds a serial `cleanid` primary key.

#### What `step1a_work_v7m.sql` does (~357 lines)

This script is the per-day worker. The Python wrapper invokes it once for each date with `-v oned='YYYY-MM-DD'`.

- **Step 1.** Copies the day's `work_pnt` rows into a temporary `work_pnt_oned`. Computes per-row `fire_size` (1.0 km for MODIS, 0.375 km for VIIRS), then half-side spans `fire_dx, fire_dy` (in degrees, latitude-corrected) and `pix_dx, pix_dy` (scaled by a hard-wired `pixfac = 1.1` and by the per-row `scan`/`track`). Builds two square geometries per detection: `geom_sml` (the "fire" square of fire_size km) and `geom_pix` (the larger pixel footprint used for the adjacency test). Adds a GIST index on `geom_pix`.
- **Step 2.1.** Builds the near-table `tbl_adj_det` by self-joining `work_pnt_oned`: two detections are neighbors when their `geom_pix` bounding boxes overlap (`&&`, uses GIST), they are within combined pixel distance (`ST_DWithin`), and they share `acq_date_use`. For each such pair it records `lhs`/`rhs` (the two cleanids) and the convex hull of the two `geom_sml` rectangles as `geom_pair`.
- **Step 2.2.** Wraps `tbl_adj_det` into a single array per day, calls `pnt2grp(lhs, rhs)` which goes into PL/Python + networkx to find connected components, and writes `(fireid, lhs, rhs, ndetect)` back into `tbl_togrp`. The smallest cleanid in each component becomes that group's `fireid`. Those fields are then copied back to `tbl_adj_det`, to `work_pnt_oned` (as `fireid1` / `ndetect1`), and finally to the persistent `work_pnt`. Detections with no neighbor get `fireid1 = cleanid` and `ndetect1 = 1` (lone detections).
- **Step 2.3.** Builds `work_lrg_oned`: for each `fireid` of two or more detections, the polygon is `ST_Union(geom_pair)` across the edges; lone detections get inserted as their `geom_sml` square.
- **Step 2.4.** Fills small holes inside the aggregated polygons: any interior ring with area < `(1/240)² square degrees` (~ 15 arcsec, roughly 0.5 km × 0.5 km) is dropped and the outer ring is rebuilt; larger holes are preserved by `ST_Difference`.
- **Step 2.5.** Updates `work_lrg_oned` with `acq_date_use`, `ndetect` (from `work_pnt_oned`), and `area_sqkm = st_area(geom_lrg, true) / 1e6` (spheroidal area in km²). Appends to the persistent `work_lrg1`.

After step 1a finishes the run on every day, `work_lrg1` is fully populated with aggressive-aggregation polygons.

### 4.3 `run_vcf.py` — early VCF/tree-cover join (used between 1a and 1b)

`run_vcf.main(tag_af, rasters, first_day, last_day, ..., algorithm_merge_aggressive_threshold_tree_cover=50)`:

This module looks for `rasters[*]['tag'].startswith('modvcf_')` and, for that one, performs an early raster-zonal mean of the `tree` band against the aggressive `work_lrg1` polygons. That tree-cover number then drives the choice between aggressive and conservative aggregation in step 1b:

- `alg_agg = 1` ("aggressive") if `v_tree >= 50`,
- `alg_agg = 2` ("conservative") otherwise.

The intuition is that in densely forested fires you want the aggressive algorithm to glom adjacent detections together (forest fires are large, contiguous), while in grass/savanna fires you want each detection treated more conservatively.

Unlike step 1a, `run_vcf.py` doesn't ship the SQL on disk — it **generates** the SQL with three helper functions and writes it to disk for `psql -f`:

- `mkcmd_create_table_oned()` — temporary table holding `:oned` date and pulling that day's slice from `work_lrg1` into `work_lrg_oned`.
- `mkcmd_create_table_continuous(tag_tbl, ['tree'], schema)` — creates `tbl_modtree_<year>` with columns `fireid`, `v_tree`, `acq_date_use`.
- `mkcmd_insert_table_continuous(...)` — the actual raster-vs-polygon zonal stats. Uses `ST_Clip(rast, geom_lrg)` to cut the raster to each polygon, then `ST_SummaryStatsAgg(clp, band, exclude_nodata=true)` to compute the mean. The CTE shape is `piece(fireid, clp, acq_date_use)` → grouped `st_summarystatsagg` → insert into `tbl_modtree_<year>`. Comments warn about multi-band-vs-touching-polygon edge cases (PostGIS tickets 3725/3730).
- `mkcmd_create_table_output(...)` / `mkcmd_insert_table_output(...)` — assembles a `work_tree` table mirroring `work_lrg_oned` plus the new `v_tree` field.
- The post-script (`mkcmd_post`, executed once after all days) updates `work_pnt.alg_agg` and `work_lrg1.alg_agg` based on `v_tree`'s comparison to the threshold.

Three .sql files are written to the current working directory: `stepvcf_prep_<tag_lct>_<tag_vcf>_<tag_regnum>_dev1.sql`, `stepvcf_work_<...>_dev1.sql`, `stepvcf_post_<...>_dev1.sql`. The Python then `psql -f`'s prep, then work (one call per day), then post.

### 4.4 `run_step1b.py` — conservative aggregation + Voronoi subdivision

`run_step1b.main(...)` mirrors the structure of `run_step1a.main` almost line for line: same `get_first_last_day` helper (literally duplicated), same retry/log machinery, but it runs `step1b_prep_<ver>.sql` and `step1b_work_<ver>.sql`.

`step1b_prep_v7m.sql` (110 lines) is mostly comments — most of its table creations were already done by `step1_prep`. It creates the conservative-aggregation table `work_lrg2`, a combined `work_lrg`, and the final `work_div` (sub-divided polygon table). `work_div` carries `polyid serial primary key`, `fireid`, `cleanids integer[]` (the AF detections inside each sub-polygon), `geom`, `acq_date_use`, `area_sqkm`, `alg_agg`.

`step1b_work_v7m.sql` (790 lines) does the per-day Voronoi sub-division and is structured as:

- **Step 1** — Same as 1a's step 1: load that day's `work_pnt` rows, compute `pix_dx`, `pix_dy`, `geom_sml`, `geom_pix`. Note that in 1b, `pix_dx`/`pix_dy` are computed using `fire_size`, not `scan`/`track`. This is the "conservative" definition (smaller adjacency radius).
- **Step 2** — Same pattern as 1a (adjacency table → `pnt2grp` → component fireids → union polygons → fill small holes). This rebuilds `work_lrg2` with conservative grouping.
- **Step 3** — Sub-division into individual fire-detection-sized polygons. This is the algorithmic heart of FINN:
  - **3.1** `tbl_close` — near-table at the much tighter threshold of 0.5 arcmin (about 0.9 km). Records `(lhs, rhs, invdist, fireid)`. Used to spot redundant near-coincident detections.
  - **3.2** `tbl_toskim` — calls `pnt2drop(lhs, rhs, invdist)` (the PL/Python iterative "kill the highest-degree node" algorithm) to identify points to drop or merge before tessellation, so the Voronoi seeds are well-spaced.
  - **3.3** `tmp_fillers` — for ties identified by `pnt2drop`, computes the centroid of each tied cluster as the replacement seed.
  - **3.4** `tmp_skmpnt` — the final seed set: original points minus dropped points, plus filler centroids.
  - **3.5** Voronoi tessellation for `fireid` groups with > 3 seeds. Calls `st_voronoi_python(x, y)`, then intersects each Voronoi cell with the corresponding `work_lrg2.geom_lrg` to get the per-detection sub-polygon. Inserts into `work_div_oned`.
  - **3.6** Custom cutter (`st_cutter_py`) for the 2- and 3-seed cases that Voronoi cannot handle in 2D.
  - **3.7** Drops "funky" geometries (using `st_polsbypopper` and other checks; some logic is commented out — the v7m algorithm is more permissive here than earlier versions).
  - **3.8** For singletons (`npnts = 1`), the original `geom_lrg` becomes the sub-polygon as-is.
  - **3.9** Back-link: for each sub-polygon, finds which `cleanid` points fall inside it (`ST_Within`) and stores them in the `cleanids integer[]` column. This is how the original AF detections survive into the final output.
  - **3.10** Pushes `work_div_oned` into the persistent `work_div`.
  - **3.11** (experimental) Updates `work_pnt.polyid` to point each detection back at its sub-polygon.

At the end of step 1b, `work_div` holds one row per *final* burned-area polygon: a fireid (the cluster), a polyid (the sub-cell), the input cleanids, the geometry, the date, the area, and the algorithm flag.

### 4.5 `run_step2.py` — raster joins for LCT, VCF, region number

`run_step2.main(tag_af, rasters, first_day, last_day, ...)` builds, on the fly, a pair of SQL scripts (`step2_prep_<rasters>_v8b.sql`, `step2_work_<rasters>_v8b.sql`) and runs them. It walks the `rasters` list — three entries for typical NRT runs (LCT, VCF, regnum) — and, depending on `kind`, calls one of three builder pairs:

- **`thematic` (LCT)** — `mkcmd_create_table_thematic` / `mkcmd_insert_table_thematic`. Creates `tbl_modlct_<year>(polyid, v_lct, f_lct, cv_lct, ct_lct, r_lct, acq_date_use)`. The insert does `ST_Clip(rast, geom)` to clip the LCT raster to each `work_div_oned` polygon, then `ST_ValueCount(clp)` to compute the histogram of land-cover classes inside the polygon, and finally selects the *majority* class (`(pvc).value` with max `count`), recording the majority value (`v_lct`), the fraction of the polygon covered by that class (`f_lct`), the pixel count of the majority class (`cv_lct`), the total pixel count (`ct_lct`), and the rank (`r_lct`).
- **`continuous` (VCF)** — `mkcmd_create_table_continuous` / `mkcmd_insert_table_continuous`. Like the early-VCF table in `run_vcf.py` but for three bands (`tree`, `herb`, `bare`): `ST_SummaryStatsAgg` per band, inserting `v_tree`, `v_herb`, `v_bare`.
- **`polygons` (regnum)** — `mkcmd_create_table_polygons` / `mkcmd_insert_table_polygons`. Joins each `work_div_oned` polygon to the region polygon layer by `ST_Intersects(r.geom, ST_Centroid(d.geom))` and copies `region_num` as `v_regnum`. Notice it uses the *centroid* of the polygon as the test point — small polygons get a single deterministic region label.
- **`input`** — `mkcmd_create_table_input` / `mkcmd_insert_table_input`. Not used by the standard NRT recipe but supported: averages a numeric attribute from `work_pnt` (joined via `polyid`) over each sub-polygon. Allows piping per-detection input like FRP into the polygon table.

After each raster has its own `tbl_*` table, `mkcmd_create_table_output` and `mkcmd_insert_table_output` create and populate the final wide table `out_<tag_lct>_<tag_vcf>_<tag_regnum>` with one row per polygon containing: `polyid`, `fireid`, `cleanids`, `geom`, `cen_lon`/`cen_lat` (centroid), `acq_date_use`, `area_sqkm`, `alg_agg`, plus all the joined fields `v_lct`, `f_lct`, `r_lct`, `v_tree`, `v_herb`, `v_bare`, `v_regnum`.

### 4.6 `export_shp.py` — write the final outputs

`export_shp.main(odir, schema, tblname, flds, shpname, date_definition)`:

- Writes `out_<...>.csv` via `psql -c "\COPY (SELECT polyid, fireid, cen_lon, cen_lat, acq_date_use AS acq_date_<lst|utc>, area_sqkm, alg_agg, <flds> FROM "<schema>"."<tblname>") TO '<csvname>' DELIMITER ',' CSV HEADER"`.
- Writes `out_<...>.shp` via `ogr2ogr -f "ESRI Shapefile" -overwrite ... -sql 'select * from <schema>.<tblname>'`. Geometry is included in the shapefile, which is the deliverable handed to downstream emissions code.

### 4.7 `run_extra.py` — diagnostics and cleanup

This is the "post-processing and cleanup" toolbox. `work_nrt.py` only calls two of its functions:

- `summarize_log(tag, out_file)` — runs `summarize_log.sql` against the `af_<tag>.tbl_log` table and writes the result to the summary file. The SQL builds three CTEs that group `tbl_log` by event and reconstruct row-count deltas across step transitions (e.g., the "agg to large" event reports how many work_pnt rows went into how many work_lrg rows). It produces a clean per-step audit trail.
- `db_use_af(tag_af, outfile)` — runs `pg_total_relation_size` queries on every table in the `af_<tag>` schema and writes both the per-table size and the schema total to the summary file.

The other functions (`clean_db_af`, `purge_db_af`, `disk_use_raster`, `db_use_raster`, `purge_hdf_raster`, `purge_tif_raster`, `purge_db_reaster`) are tools for clearing intermediate state — used by other driver scripts in `code_bashinterface/` (e.g., reset/wipe utilities) but not by `work_nrt.py`.

### 4.8 Modules touched indirectly

`work_nrt.py` does not import these directly, but they are part of the system:

- `rst_import.py` — the long (929-line) MODIS raster downloader/importer. `work_common.sec2_check_environment` calls `rst_import.prep_modis_tile()` to make sure `raster.wireframe` is populated. The full raster import (HDF download from Earthdata, mosaic to GeoTIFF, `raster2pgsql` load into the `raster` schema) is handled separately by `do_everything.py`-style driver scripts, not by `work_nrt.py`. NRT operation assumes the rasters are already loaded.
- `modis_tile.py` — generates the MODIS sinusoidal tile geometry that `prep_modis_tile` loads into `raster.wireframe`.
- `polygon_import.py` — utility for loading a polygon shapefile (such as `regnum`) into PostgreSQL via `ogr2ogr`. The region polygon also must be loaded before `work_nrt.py` is run.
- `downloader.py` — Earthdata HTTP downloader used by `rst_import`.

## 5. Pipeline summary — what happens end-to-end on one call

A complete invocation like

```bash
work_nrt.py -t mytag -y 2020 -o ./out -fd 2024100 -ld 2024110 myfires.shp
```

unfolds as follows:

1. **Environment setup** (`work_common`): PG* environment variables are set, Earthdata creds prompted, `sys.path` extended.
2. **Configuration** (`sec1_user_config`): the raster tag triplet (`modlct_2020`, `modvcf_2020`, `regnum`) and the rasters list are baked into module globals.
3. **Environment check** (`sec2_check_environment`): PostgreSQL/PostGIS versions printed, `testpy()` runs to confirm numpy + networkx in PL/Python, `prep_modis_tile()` makes sure `raster.wireframe` exists.
4. **AF import** (`sec3_import_af` → `af_import.main`): `myfires.shp` is unzipped if needed, then loaded as `af_mytag.af_in_1` via `ogr2ogr`. Schema `af_mytag` is dropped and recreated.
5. **Step 1 prep** (`run_step1a.main` first half → `step1_prep_v7m.sql`): the per-run schema is fully built (functions, types, log table). All AF tables are imported into `work_pnt`, tropics are duplicated forward one day, the date window is enforced, low-confidence and (optionally) persistent-source detections are dropped.
6. **Step 1a work** (`run_step1a.main` second half → `step1a_work_v7m.sql` once per day): aggressive clustering, large-polygon construction, small-hole filling. `work_lrg1` gets populated.
7. **VCF tree-cover join** (`run_vcf.main` → generated `stepvcf_*.sql`): tree cover percentage is averaged inside each `work_lrg1` polygon; if it exceeds 50%, the polygon is flagged `alg_agg = 1` (aggressive). Otherwise it is flagged `alg_agg = 2` (conservative).
8. **Step 1b** (`run_step1b.main` → `step1b_prep_v7m.sql`, then `step1b_work_v7m.sql` once per day): conservative clustering builds `work_lrg2`; then for each large polygon, points are skimmed via `pnt2drop`, Voronoi or custom-cutter tessellation is performed via `st_voronoi_python` / `st_cutter_py`, the sub-polygons are intersected with the parent, "funky" geometries are filtered, and sub-polygons are written to `work_div` with back-links to the originating AF detections.
9. **Step 2** (`run_step2.main` → generated `step2_*.sql`): for each polygon in `work_div`, raster zonal majority (LCT), zonal means (VCF), and centroid-region (regnum) are computed and joined into a wide output table `out_modlct_2020_modvcf_2020_regnum`.
10. **Export** (`sec7_export_output` → `export_shp.main`): the output table is dumped to `./out/out_mytag_modlct_2020_modvcf_2020_regnum.csv` and the matching `.shp`.
11. **Summary** (`run_extra.summarize_log` + `db_use_af`): the `tbl_log` audit and a per-table disk-usage report are appended to the summary file.

At that point the user has, in `./out`, a CSV plus shapefile of one polygon per FINN-sub-divided burned area, with date, fire ID, area in km², majority land-cover class and its fraction, mean tree/herb/bare percent, region number, and the algorithm flag — exactly the inputs the FINN emissions code consumes downstream.

## 6. Design notes and idioms worth knowing

A few things become clearer when you read the whole pipeline rather than any one script:

- **The DB is the workflow engine.** Python is glue. Each step writes its results back into the same `af_<tag>` schema so the next step can read them. There is no in-process Python state that crosses stages; if you stop and re-run with `--run_step2=true --run_import=false --run_step1=false`, the previous step's state is just there in the database.
- **Per-day temporary tables.** Both `step1a_work` and `step1b_work` create `work_pnt_oned`, `work_lrg_oned`, `work_div_oned` as TEMPORARY tables, do all the heavy SQL on them, then `INSERT INTO` the persistent counterpart at the end. This keeps memory and lock scope bounded per day and lets the work scripts be invoked once per date by Python in a `for date in dates:` loop — a natural unit of parallelism even though `run_step1a` doesn't actually parallelize them.
- **PL/Python for graph and Voronoi work.** Connected-component analysis and Voronoi tessellation are done in PL/Python (`networkx`, `scipy.spatial.Voronoi`) because SQL/PostGIS doesn't naturally handle either at the scale and quality needed. The `testpy()` function is the canary that catches a Python-binding mismatch *before* the long run.
- **Two clustering algorithms in parallel.** `work_pnt` carries `fireid1`/`ndetect1` (aggressive) *and* `fireid2`/`ndetect2` (conservative), and `alg_agg` is per-polygon. The choice between them is made by VCF tree-cover (the call between 1a and 1b). This is why step 1 is split into two halves with VCF in the middle.
- **Time vs space convention.** `acq_date_use` is whichever the user picked at run time (LST or UTC) and is the *only* date used downstream. `acq_date_lst`/`acq_date_utc` are kept around for audit but never used for grouping after `step1_prep`.
- **Tropics duplication.** MODIS's twice-daily polar overpasses leave tropical equator-region detections potentially split across the UTC-day boundary. Duplicating all `abs(lat) <= 23.5` MODIS rows forward one day, then discarding the first day of the output, sidesteps this gracefully.
- **Logging is structured.** Almost every major mutation is wrapped in `log_checkin(event, table, count_before, oned)` ... `log_checkout(id, count_after)`. `summarize_log.sql` then reports row-count deltas for events like `'agg to large'`, `'subdiv'`, `'join modlct_<year>'`, and `'merge all'`. This is how a daily-run operator audits what changed without re-running the pipeline.
- **The version suffix.** `step1_prep_v7m.sql`, `step1a_work_v7m.sql`, `step1b_*_v7m.sql` are version `v7m` of FINN's geometry algorithm. `run_step1a` and `run_step1b` have a long block of `if vorimp == 'scipy_fixcutter_...': ver='v7?'` branches kept for backward compatibility, but the script defaults to `ver='v7m'` and never reaches them in the NRT flow. `run_step2` runs `v8b`.
- **NRT vs archive.** "NRT" means "near-real-time": the daily FIRMS NRT product lacks the `Type` field, so `filter_persistent_sources` quietly becomes a no-op when the prep script discovers that no `af_in_*` table has a `type` column. Running the same code against archived `fire_archive_*.shp` (which does have `type`) will actually drop volcanoes.
