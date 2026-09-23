#!/usr/bin/env python3
"""work_nrt.py — pure-Python replacement.

Drop-in for ``preprocessor/code_bashinterface/work_nrt.py``.  Uses the
same flag set so any existing wrapper script keeps working.
"""
import sys
from finn_py.cli import nrt_main

if __name__ == "__main__":
    sys.exit(nrt_main())
