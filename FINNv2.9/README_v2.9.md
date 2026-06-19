# finn_py — pure-Python FINN preprocessor

A reimplementation of the FINN preprocessor (`work_nrt.py`, `work_raster.py`
and the `code_anaconda/*.sql` they call) with no PostgreSQL / PostGIS / psql
dependency. The MODIS land-cover and VCF rasters live as Cloud-Optimized
GeoTIFFs on disk; the per-run "schema" lives as GeoParquet files in a
working directory.

See `DESIGN.md` for the table-by-table mapping of what was replaced with
what.


## Building global mosaics for many years

Some workflows need a global LCT + VCF mosaic per year (e.g. building a
2002–2025 archive once and reusing it across daily NRT runs).  Two
things make the global case different from regional NRT:

1. **`build_mosaic` streams.**  The destination is opened as a tiled
   GeoTIFF and GDAL's warper writes it block-by-block — peak RAM is
   bounded by the `warp_mem_mb` argument (default 2 GB) regardless of
   how big the output is.  A global LCT mosaic at 500 m is fine in 2 GB;
   a global VCF mosaic at 250 m is also fine.

2. **Years are independent.**  `scripts/build_global_rasters.py` is a
   batch driver that loops over years (or runs one year at a time for
   job-array dispatch).  See `scripts/build_global.pbs.template` for
   a Casper PBS template that fans out 2002–2025 as a 24-job array.

Expected input layout::

    <lct-root>/2002/MCD12Q1.A2002001.*.hdf
    <lct-root>/2003/MCD12Q1.A2003001.*.hdf
    ...
    <vcf-root>/2002/MOD44B.A2002065.*.hdf
    ...

(override with `--lct-glob` / `--vcf-glob` if your layout differs).

Output layout (one set of three files in the same dir, one LCT + one
VCF tif per year, plus a single `regnum.gpkg` that all years share)::

    <raster-dir>/modlct_2002.tif
    <raster-dir>/modvcf_2002.tif
    <raster-dir>/modlct_2003.tif
    ...
    <raster-dir>/regnum.gpkg

Then `finn-nrt -y 2003 --raster_dir <raster-dir> ...` picks up the
matching year automatically (the run uses `modlct_<year>.tif` /
`modvcf_<year>.tif`).

### A note on GDAL availability

The streaming mosaic step builds a tiny VRT file referencing the
per-tile sinusoidal GeoTIFFs and then warps it.  That requires either
the `osgeo.gdal` Python bindings *or* the `gdalbuildvrt` CLI tool.  On
HPC conda environments both are typically present; on pip-only
installs you may have neither and will get a clear error pointing at
`conda install -c conda-forge gdal`.

## Usage (regional / single year)

This package does not download anything.  You are responsible for
fetching the MODIS HDF tiles and the region polygon file yourself (e.g.
with `earthaccess`, `pyMODIS`, or NASA's wget recipes) and pointing
`finn-raster` at the local paths.

```bash
# Step A — mosaic local MODIS HDF tiles into the inputs finn-nrt expects
finn-raster -y 2020 \
    --lct      /data/MCD12Q1/2020/MCD12Q1.A2020001.*.hdf \
    --vcf      /data/MOD44B/2020/MOD44B.A2020065.*.hdf \
    --regions  /data/regions/All_Countries.shp \
    --raster_dir ./rasters
# produces:  ./rasters/modlct_2020.tif
#            ./rasters/modvcf_2020.tif
#            ./rasters/regnum.gpkg

# Step B — run the AF -> burned-area-polygon pipeline
finn-nrt -t mytag -y 2020 -o ./out \
    --work_dir   ./work \
    --raster_dir ./rasters \
    -fd 2020001 -ld 2020031 \
    myfires.shp

# outputs in ./out/:
#   out_mytag_modlct_2020_modvcf_2020_regnum.csv
#   out_mytag_modlct_2020_modvcf_2020_regnum.gpkg
```

## Project layout

```
finn_preproc_py/
├── DESIGN.md
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
    ├── work_nrt.py                drop-in CLI shim
    ├── work_raster.py             drop-in CLI shim
    ├── build_global_rasters.py    batch driver: build per-year global mosaics
    └── build_global.pbs.template  Casper PBS job-array template (2002-2025)
```

## Function-by-function map to the original

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
| PL/Python `st_voronoi_python`                        | `step1._voronoi_cells` (scipy.spatial.Voronoi)            |
| `st_cutter_py` (2/3 point cutter)                    | `step1._custom_cutter`                                    |
| `mkcmd_insert_table_thematic` (LCT majority)         | `step2.thematic_zonal_stats`                              |
| `mkcmd_insert_table_continuous` (VCF mean)           | `step2.continuous_zonal_stats`                            |
| `mkcmd_insert_table_polygons` (regnum centroid)      | `step2.polygons_zonal_join`                               |
| `export_shp.main`                                    | `export.write_output`                                     |
| `rst_import.Importer` (HDF → PostGIS raster)         | `rasters.build_mosaic` (HDF → COG, PyHDF-based)           |
| `modis_tile.py` + `raster.wireframe`                 | `rasters.modis_tile_polygons` / `rasters.tiles_needed`    |
| `tbl_log` + `summarize_log.sql`                      | `pipeline._Logger` (JSON event log + summary text)        |

## Things still to wire up

A few things are intentionally minimal in this first cut:

- **Parallelism.** `FinnConfig.workers` and the `--workers` CLI flag are
  plumbed but the per-day loop in `pipeline.run_nrt` is sequential; the
  per-day calls (`step1a_one_day`, `step1b_*_one_day`, `step2.zonal_join_one_day`)
  are pure functions that can be dispatched to a `ProcessPoolExecutor`
  once you've decided how to share `work_pnt` between workers (a thin
  read-only mmap on the parquet, most likely).
- **Validation of MOD44B band order.** I assume the three bands in the
  output mosaic land in the order
  `(Percent_Tree_Cover, Percent_NonTree_Vegetation, Percent_NonVegetated)`
  matching `SDS_VCF_MOD44B`.  Worth double-checking against a real file
  before trusting `v_tree`/`v_herb`/`v_bare`.
- **No CRS-mismatch retries** when a polygon barely misses a raster
  edge.  The original SQL caught those by polyid; here we just return
  NaN for that polygon and log it.
- **Input acquisition is your problem.** This package will not download
  MODIS HDFs or the region polygons.  Use whatever you already use —
  `earthaccess`, `pyMODIS`, `wget`, the LP DAAC website — and hand
  `finn-raster` the local paths.
