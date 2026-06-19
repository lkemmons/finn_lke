#!/usr/bin/env python3
"""work_raster.py — pure-Python replacement.

Drop-in for ``preprocessor/code_bashinterface/work_raster.py``.
"""
import sys
from finn_py.cli import raster_main

if __name__ == "__main__":
    sys.exit(raster_main())
