#!/usr/bin/env python3
"""Quick look at a work_*.parquet file from a finn_py run.

Usage:
    python peek_parquet.py <file.parquet>             # summary
    python peek_parquet.py <file.parquet> --head 10   # first 10 rows
    python peek_parquet.py <file.parquet> --where 'fireid == 5453'
    python peek_parquet.py <file.parquet> --cols polyid fireid area_sqkm
"""
import argparse, sys
from pathlib import Path
import geopandas as gpd
import pandas as pd

p = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('path', type=Path)
p.add_argument('--head', type=int, default=0,
               help='show first N rows (default: just summary)')
p.add_argument('--where', default=None,
               help='pandas .query() expression to filter rows')
p.add_argument('--cols', nargs='+', default=None,
               help='restrict columns to this list')
args = p.parse_args()

g = gpd.read_parquet(args.path)
if args.where:
    g = g.query(args.where)

print(f'file:    {args.path}')
print(f'rows:    {len(g):,}')
print(f'crs:     {g.crs}')
print(f'columns: {list(g.columns)}')
print(f'dtypes:\n{g.dtypes.to_string()}')

if 'area_sqkm' in g.columns:
    a = g['area_sqkm']
    print(f'\narea_sqkm: total={a.sum():.1f}  '
          f'median={a.median():.3f}  max={a.max():.3f}')
if 'fireid' in g.columns:
    print(f'distinct fireids: {g["fireid"].nunique():,}')
if 'polyid' in g.columns:
    print(f'distinct polyids: {g["polyid"].nunique():,}')
if 'acq_date_use' in g.columns:
    print(f'date range: {g.acq_date_use.min()} .. {g.acq_date_use.max()}')

if args.head:
    cols = args.cols or [c for c in g.columns if c not in ('geom', 'geom_pnt', 'geom_sml')]
    print()
    print(g[cols].head(args.head).to_string(index=False))
