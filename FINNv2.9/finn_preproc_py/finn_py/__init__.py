"""Pure-Python FINN preprocessor.

A reimplementation of the FINN `code_bashinterface` + `code_anaconda`
pipeline that does not require PostgreSQL / PostGIS / psql / ogr2ogr.

Public entry points:
    finn_py.pipeline.run_nrt(...)        — full AF → polygons pipeline
    finn_py.raster_pipeline.run_raster() — mosaic local MODIS HDFs

`run_raster` takes paths to MODIS HDFs already on disk and does not
download anything — you are expected to fetch the inputs yourself
(e.g. with `earthaccess`, `pyMODIS`, NASA's `wget` recipes, or your
existing tooling).

See DESIGN.md for the table-to-file mapping.
"""

from .config import FinnConfig, default_rasters

__all__ = ["FinnConfig", "default_rasters"]
