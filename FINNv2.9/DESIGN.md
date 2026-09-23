# finn_py — design notes

This document explains the design choices behind the pure-Python
FINN rewrite and the surrounding toolkit.  It's meant for anyone
picking up the code who wants to know *why* the code is shaped the
way it is — not just how to use it (that's `README.md`).

## Scope

The project started as a rewrite of the FINN preprocessor
(`preprocessor/code_bashinterface/work_nrt.py`, `work_raster.py`, and
the `code_anaconda/*.sql` files they call) with no PostgreSQL, PostGIS,
psql or ogr2ogr dependency.  It has grown to cover the emissions
calculation (a port of the IDL `finn2_9_2_calc_emis_nrt_daily.pro`) and
downstream tools for gridding, extending, verifying and visualizing
FINN NetCDF outputs.

All the parts are independent scripts that talk through files.  There
is no daemon, no database and no in-memory state that spans script
invocations.

## Design axes

### 1. Storage is on disk, in open formats

- **PostgreSQL schema `af_<tag>`** → working dir `work/<tag>/` with
  GeoParquet files.  One parquet per stage (`work_pnt.parquet`,
  `work_lrg1.parquet`, `work_lrg2.parquet`, `work_div.parquet`).
- **PostGIS raster** (`rst_modlct_*`, `rst_modvcf_*`) → Cloud-Optimized
  GeoTIFF in a `rasters/` directory.  One TIF per year per product.
- **Emissions**: NetCDF with `xarray`/`netCDF4`.  Either per-fire
  (`fire` or `polyid` dim) or gridded (`lat`/`lon` or `ncol`).
- **Regional / country polygons**: GeoPackage.

Everything is inspectable with standard tools (`geopandas`,
`ogrinfo`, `ncdump`, `gdalinfo`, `xarray`).  No proprietary or
process-lifetime state.

### 2. One script = one job

Each script in `scripts/` does one thing and communicates via files.
No script imports from another script.  Shared logic lives in the
`finn_py/` package.  This means:

- You can rerun any stage in isolation.
- Every stage's output is human-readable (or at least tool-readable
  without needing the whole toolkit installed).
- Adding a new stage is just adding a script; no rewiring.

### 3. HDF4 via PyHDF, not GDAL

Many GDAL builds (conda-forge wheels, recent Homebrew bottles) ship
without the HDF4 driver.  `rasterio.open(hdf_path)` then raises.  We
use PyHDF for raw array access and parse `StructMetadata.0` in Python
to recover each tile's sinusoidal transform, then convert each tile
to a plain intermediate GeoTIFF before mosaicking.  The downstream
pipeline never touches HDF4 again.

## Preprocessor (AF → burned-area polygons)

### What's mapped to what

| Original                                                  | Pure-Python replacement                                       |
|-----------------------------------------------------------|---------------------------------------------------------------|
| PostgreSQL schema `af_<tag>`                              | Working dir `work/<tag>/` with GeoParquet files               |
| PostGIS `raster` schema (`rst_modlct_*`, `rst_modvcf_*`)  | `rasters/<tag>/*.tif` Cloud-Optimized GeoTIFFs                |
| `ogr2ogr ... PG:dbname=finn`                              | `geopandas.read_file` / `to_parquet`                          |
| `psql -f step1_prep_v7m.sql`                              | `finn_py.step1.prep`                                          |
| `psql -f step1a_work_v7m.sql -v oned=...`                 | `finn_py.step1.step1a_one_day`                                |
| `psql -f step1b_work_v7m.sql -v oned=...`                 | `finn_py.step1.step1b_one_day`                                |
| `run_vcf.py`'s on-the-fly SQL                             | `finn_py.step1.set_alg_agg_from_tree`                         |
| `psql -f step2_work_*_v8b.sql`                            | `finn_py.step2.zonal_join_one_day`                            |
| `ST_Union`, `ST_Intersection`, `ST_Difference`, …         | `shapely.ops.unary_union`, `.intersection`, …                 |
| `ST_DWithin` adjacency                                    | `scipy.spatial.cKDTree.query_pairs`                           |
| `ST_Voronoi` / `st_voronoi_python`                        | `scipy.spatial.Voronoi` with 4-corner ghost points            |
| `pnt2grp` (connected components)                          | `networkx.connected_components`                               |
| `pnt2drop` (point skimming)                               | `step1._pnt2drop` (networkx)                                  |
| `raster2pgsql ... | psql`                                 | `rasters.build_mosaic` (PyHDF → sinu GeoTIFF → warp → COG)    |
| `ST_Clip` + `ST_ValueCount` (LCT majority)                | `rasterio.mask.mask` + `numpy.unique(return_counts=True)`     |
| `ST_Clip` + `ST_SummaryStatsAgg` (VCF mean)               | `rasterio.mask.mask` + `numpy.nanmean`                        |
| Region polygon centroid join                              | `geopandas.sjoin` on centroids                                |
| `tbl_log` + `summarize_log.sql`                           | structured JSON event log + human-readable summary text       |
| `export_shp.py`                                           | `GeoDataFrame.to_file` + `.to_csv`                            |

### Working-table schemas (GeoParquet)

`work_pnt.parquet` — one row per AF detection that survives initial
filtering:

| column            | dtype               | notes                                                                  |
|-------------------|---------------------|------------------------------------------------------------------------|
| cleanid           | int64 (index)       | PK; assigned after all filtering                                       |
| rawid             | int64               | row number across input files                                          |
| src_file          | int                 | which input file (1, 2, …)                                             |
| geom_pnt          | geometry (Point)    | WGS84                                                                  |
| lon, lat          | float64             | denormalized from geometry                                             |
| scan, track       | float64             | from MODIS/VIIRS product                                               |
| acq_date_utc      | date                |                                                                        |
| acq_time_utc      | string (HHMM)       |                                                                        |
| acq_date_lst      | date                | approximate local solar time                                           |
| acq_datetime_lst  | timestamp           |                                                                        |
| acq_date_use      | date                | whichever the user selected; the only one used downstream              |
| instrument        | string              | 'MODIS' or 'VIIRS'                                                     |
| confident         | bool                | instrument-specific rule                                               |
| anomtype          | int                 | 0–3 from FIRMS 'Type' field (0 if absent)                              |
| frp               | float64             |                                                                        |
| fireid1, ndetect1 | int64, int          | aggressive cluster id and its detection count                          |
| fireid2, ndetect2 | int64, int          | conservative cluster id and its detection count                        |
| alg_agg           | int                 | 1 = aggressive, 2 = conservative (tree-cover gated)                    |
| polyid            | int64               | back-link to work_div                                                  |

`work_lrg1.parquet` / `work_lrg2.parquet` / `work_lrg.parquet` — one
row per aggregated burned area:

| column        | dtype                   |
|---------------|-------------------------|
| fireid        | int64 (index)           |
| geom_lrg      | geometry (Polygon)      |
| acq_date_use  | date                    |
| ndetect       | int                     |
| area_sqkm     | float64                 |
| alg_agg       | int                     |
| v_tree        | float64 (work_lrg only) |

`work_div.parquet` — one row per Voronoi sub-polygon (the unit on
which step 2 rasters are joined):

| column        | dtype                | notes                              |
|---------------|----------------------|------------------------------------|
| polyid        | int64 (index)        |                                    |
| fireid        | int64                |                                    |
| cleanids      | list[int64]          | AF detections inside this sub-polygon |
| geom          | geometry (Polygon)   |                                    |
| acq_date_use  | date                 |                                    |
| area_sqkm     | float64              |                                    |
| alg_agg       | int                  |                                    |

`out.gpkg` — final deliverable, mirrors `work_div` plus joined raster
fields `v_lct`, `f_lct`, `v_tree`, `v_herb`, `v_bare`, `v_regnum`,
`cen_lon`, `cen_lat`.  For LCT multi-row output, `f_lct` is the
fraction of that polygon covered by that LCT class (rows within a
polyid sum to 1.0).

### Pipeline order (matches original)

1. **Load AF**           `af_io.load_af_files(fnames)` → DataFrame
2. **Build work_pnt**    `step1.prep(...)` → `work_pnt.parquet`
3. **Step 1a per day**   `step1.step1a_one_day(...)` → appends to `work_lrg1`
4. **Tree-cover join**   `step1.set_alg_agg_from_tree(...)` → sets `alg_agg`
5. **Step 1b per day**   `step1.step1b_one_day(...)` → appends to `work_lrg2` / `work_div`
6. **Step 2 per day**    `step2.zonal_join_one_day(...)` → appends to `out`
7. **Export**            `export.write_output(...)` → `.gpkg` + `.csv`

### Tropical carryover

MODIS has swath gaps in the tropics that make fires appear only every
1–2 days.  The original FINN's `step1_prep_v7m.sql` compensates by
duplicating every tropical MODIS detection (|lat| ≤ 23.5°) into the
next calendar day.  Without this, finn_py produces noticeably fewer
tropical polygons than FINN2.

This behavior is preserved in `step1.prep`, gated on the
`duplicate_tropical_modis` kwarg (default True).  The `FinnConfig`
field of the same name is threaded through by `pipeline.run_nrt`.
Both `run_daily_nrt.py` and `work_archive.py` expose a
`--no-tropical-carryover` flag that sets it to False.

Two pieces of context matter for understanding the design:

1. **Automatic (step1.prep) vs. explicit (carryover_files).** step1.prep
   duplicates all tropical MODIS forward one day, using rows *already
   in the loaded DataFrame*.  For a multi-day archive that's enough:
   yesterday's rows are already in the file.  For NRT (one file per
   day) that's not enough — yesterday isn't loaded — so the driver
   optionally loads yesterday's file separately and passes it as
   `tropical_carryover_files`.  `af_io.load_af_files` then filters
   those rows to the tropics and re-dates them to today.
2. **Don't double-count.** If you pass a multi-day archive as *both*
   `af_files` and `tropical_carryover_files`, the pipeline would (a)
   read yesterday's rows once via `af_files`, (b) read them again via
   `tropical_carryover_files` and re-date to today, then (c)
   step1.prep would duplicate the re-dated rows forward one more day.
   `work_archive.py` therefore never passes archive files as
   `tropical_carryover_files`.

### Things that change behaviorally vs. the original

- **MODIS HDFs are read with PyHDF, not GDAL** (see §"HDF4 via PyHDF").
- **`raster.wireframe` is dropped.**  Tile geometry is computed
  in-memory from the sinusoidal projection definition; no on-disk table.
- **`testpy()` is dropped.**  Python is the runtime, not PL/Python
  inside Postgres.
- **`work_div_oned` is no longer temporary.**  Per-day work is held
  in-memory and appended to the working parquet at the end of each day.
- **No `tbl_log`** — replaced by a `logging` stream and a JSON summary
  file written by `pipeline._Logger`.
- **Concurrency model.**  The original is sequential per-day with the
  DB doing the heavy work.  The Python version's per-day loop is
  sequential too, but per-day inputs are independent after prep and can
  be dispatched to a `ProcessPoolExecutor` — the `--workers` CLI flag
  is plumbed for this though the parallel executor is not yet wired in.

### Things deliberately preserved

- **Hard-wired constants** (`pixfac = 1.1`, fire size 1.0 km MODIS /
  0.375 km VIIRS, small-hole threshold `(1/240)² deg²`, 0.5 arcmin
  skim distance, 23.5° tropics latitude, 50% tree-cover threshold for
  `alg_agg`) — same numbers, same `_v7m` semantics.
- **Tropics duplication** for MODIS.
- **Two-algorithm structure**: aggressive then conservative
  aggregation, with a tree-cover-gated `alg_agg` flag in between.
- **Per-day independence** of step 1a / 1b / 2 work.

## Emissions calculation

`scripts/calc_emis_daily.py` is a direct port of
`finn2_9_2_calc_emis_nrt_daily.pro` (IDL).  Every algorithmic detail
(genveg assignment, VCF cleanup that maps MOD44B fill codes 200/251/
252/253 to sensible fallbacks, fuel-load lookups, the three-branch
biomass formula, urban-lct-13 mutation to 10/8/1/5, the gas-vs-aerosol
unit conversion) is preserved with the same numeric behavior.

The **only** intentional change vs. the IDL is the output format:
NetCDF instead of a CSV.  This gives us free metadata (units, MW,
species type) and free tooling downstream (`xarray`, `netCDF4`).

Per-fire NetCDF layout:

```
dimensions:
    fire = <N fires>
variables:
    float32 <species>(fire)    // one variable per emission species
        <species>:units = "mol/day"      (gas)
                    | "kg/day"           (aerosol)
                    | "number/day"       (number)
        <species>:molecular_weight_g_per_mol = ...   // 1 for aerosols
```

The MW = 1 sentinel for aerosols is the CAM-Chem convention: it makes
the mass-to-"molecule" formula uniform across gases and aerosols
(`mass_kg = molec × MW / (Avogadro × 1000)` recovers the original kg
for MW=1 and the true mass for real MWs).

## Gridding conventions (CAM/CESM)

Both `grid_emissions.py` (per-fire NetCDF input) and
`grid_txt_emissions.py` (legacy text CSV input) produce output with
the **CAM/CESM emission-input schema**:

```
dimensions:
    time = UNLIMITED     // one entry per hour (24 per day) or per day
    lat  = ...           // regular grid, OR
    ncol = ...           // unstructured (SCRIP)
variables:
    double  time(time)       // fractional days since epoch
        time:units = "days since 1970-01-01 00:00:00"
        time:calendar = "gregorian"
    int32   date(time)       // YYYYMMDD, repeated across the 24 hours of a day
    int32   datesec(time)    // seconds of day: 0, 3600, ..., 82800
        datesec:units = "s"
    double  lat(lat), lon(lon)    // regular grid
        // OR
    double  lat(ncol), lon(ncol)  // unstructured (cell centers from SCRIP)
    float32 <species>(time, lat, lon)   // or (time, ncol)
        <species>:units = "molecules/cm2/s"   // (or "number/cm2/s" for particle counts)
        <species>:molecular_weight_g_per_mol = ...
```

Key properties:

- **Single `time` dim** — no separate `date` and `hour` dims.
  Concatenating across days is a plain `xr.concat` along `time`.
  `date` and `datesec` are companion variables that CAM/MOZART reads
  directly.
- **`date`/`datesec` are consistent with `time`.**
  `time == date_as_days_since_epoch + datesec / 86400` holds exactly
  (the daily portion of a concat file is normalized to midnight of
  each date so this stays true).
- **Uniform flux units** across gas / aerosol / number species so
  downstream tools don't need to branch.  For MW=1 aerosols the
  "molecules/cm2/s" number is 1000/MW × Avogadro × (kg/hour) /
  (cell_cm² × 3600 s); the inverse recovers kg.

### Cell-center precision fix

Coordinate values like `-90 + 0.1 × 800 + 0.05` accumulate float64
fuzz (`-9.849999999999998` instead of `-9.85`).  `_decimals_for_step`
in `grid_emissions.py` computes the decimal precision the grid step
needs plus a two-decimal safety margin, and rounds the cell-center
arrays to that precision before writing them to NetCDF.  Values in
the file are then the "expected" clean decimals with no fuzz.

Only the *center* arrays are rounded — edges keep their raw values so
the binning `searchsorted` stays exact.

## Legacy text-CSV conventions

Legacy FINN text files carry two conventions that need special
handling:

### 1. HOUR values are unreliable (FINNv2.7-hourly)

Some FINN text outputs have the `HOUR` column but its values don't
match the actual UTC hour of each row.  `datetimeUTC` is trustworthy
(ISO timestamp).  So the design is:

- **HOUR column presence** = signal that the file is *intended* to be
  hourly (stable structural fact).
- **HOUR values** = ignored.  Actual UTC hour comes from
  `pd.to_datetime(df["datetimeUTC"]).dt.hour`.
- **HOUR column** = dropped before writing, so it can't leak into
  outputs.

Both `grid_txt_emissions.py` and `txt_to_perfire_nc.py` follow this
rule.  For `txt_to_perfire_nc.py`, mode detection uses HOUR-column
presence (before it's dropped); hour binning uses `datetimeUTC`.

### 2. Species aren't in variable names

Some FINN files (both daily and hourly) call the emission variable
`fire` and put the actual species in a variable attribute:

- `fire:map = "CO->CO"` (or `"NO2->NOx_asNO"` — output species on the right)
- `fire:long_name = "CO fire emissions"` (species is the leading token)

`concat_finn_daily_hourly.py` and `append_daily_to_hourly.py` extract
the species from these attributes when the variable name is generic
(`fire`, `emis`, `emission`, `emissions`, `flux`).  Species pairing
across two directories of `fire`-named files then just works — the
key isn't the variable name, it's the extracted species.

## Concat / append tools

Both tools (`concat_finn_daily_hourly.py` and
`append_daily_to_hourly.py`) share the same building blocks and
differ only in *where* the replicated day goes:

- Concat: prepend 24 hourly rows to the start of the hourly file.
- Append: append 24 hourly rows to the end.

Common behavior:

- **Rebase daily onto hourly's epoch.**  Real-world files have
  different time units — daily on `days since 1950-01-01`, hourly on
  `days since 1970-01-01`.  Because the tools rebuild the daily's
  `time` coord from scratch from the chosen date and the hourly
  file's units string, the daily is silently rebased.
- **Calendar check.**  Gregorian family (`gregorian`, `standard`,
  `proleptic_gregorian`) treated as equivalent since they differ only
  pre-1583.  Other calendars (`360_day`, `noleap`, `julian`) refused
  with a clear error pointing at `cdo setcalendar`.
- **Normalize coord/data-var status.**  If `lat`/`lon` are data
  variables in one file and coordinates in the other (which happens
  when `fire:coordinates = "lat lon"` appears in one but not the
  other), promote them to coordinates before concat.  Drop `ncol`,
  `area` and `rrfac` which appear inconsistently across FINN file
  conventions.
- **Never double-count.**  The concat tool never passes archive files
  as `tropical_carryover_files`; the pipeline handles duplication
  internally via `step1.prep`.

## Regional totals

`compute_emis_totals.py` sums a gridded file over one or more regions
per timestep.  Two design choices worth mentioning:

### Half-open region masks

Regions are defined as bounding boxes `(W, E, S, N)`.  The mask uses
half-open intervals: `south ≤ lat < north` and `west ≤ lon < east`.

Without this, a cell at exactly lat=23.5° in a 1° grid would be in
both `tropics` (`lat ≤ 23.5`) and `nh_extratropics` (`lat ≥ 23.5`) —
double-counted.  With half-open, the cell falls into exactly one
region (in this case `nh_extratropics`).  The invariants
`global = NH + SH` and `global = tropics + nh_extratropics +
sh_extratropics` then hold to machine precision.

### Uniform mass conversion

For every species with an MW attribute, mass in kg is computed as
`molec × MW / (Avogadro × 1000)`.  This works uniformly:

- For gases (real MW), `molec / N_A = mol`, and `mol × MW / 1000` = kg.
- For aerosols (MW=1 sentinel), `molec × 1 / (N_A × 1000)` recovers
  the original kg because the flux conversion used
  `(kg/hour) × 1000/MW × N_A`.

Number species (no MW) stay in raw count units.

## Plot conventions

`plot_emissions_scrip.py` renders each SCRIP cell as its own filled
polygon (`matplotlib.collections.PolyCollection`) — no interpolation.
Design details:

- **7 corners with padding.**  MPAS SCRIP files use up to 7 corners
  per cell; pentagons fill the 7th slot by repeating the last vertex.
  matplotlib draws zero-length edges silently, so no explicit handling
  is needed.
- **Dateline unwrap.**  Corner longitudes are unwrapped relative to
  each cell's center so a cell that spans the dateline has corners
  within ±180° of its center.  Cells whose centers fall within a
  configurable seam distance (default 5°) of ±180° are drawn twice —
  once with the unwrapped longitudes, once shifted by ±360° — so the
  seam renders cleanly on both edges of a PlateCarree map.
- **Log color scale bounded at 6 decades.**  Default `vmin` is the
  greater of the 1st-percentile value and `vmax / 1e6`, so a single
  outlier can't crush the whole visible range onto one color.

## Publications and audit trail

Every long-running script emits a structured log line for each
significant step (via Python `logging`), plus a human-readable
`summary_<tag>.txt` per run.  The intent is that a reproduction check
one year later should need only the input files and the log — no
notebook or oral history required.

## Things still to wire up

- **Parallel per-day dispatch** in `pipeline.run_nrt`.  Per-day inputs
  are independent after `prep`; a `ProcessPoolExecutor` on the parquet
  files would use available cores.  The `--workers` CLI flag is
  plumbed but the executor isn't hooked up yet.
- **`tropical_lat_bounds` plumbing.**  The `FinnConfig` field is
  respected by `af_io.load_af_files` (for explicit carryover), but
  `step1.prep` still uses the module-level `TROPICS_LAT_DEG = 23.5`
  constant for its automatic duplication.  Fine as long as you want
  the default; would need a small change to make the config field
  authoritative.
- **CRS-mismatch retries** at raster edges.  The original SQL caught
  those by polyid; here we return NaN for that polygon and log it.
- **Non-gregorian calendars in concat/append.**  Currently refused
  with a clear error.  cftime-based handling is tractable but a
  bigger change.

