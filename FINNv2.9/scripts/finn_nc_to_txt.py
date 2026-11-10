#!/usr/bin/env python3
"""Write FINNv2.9.2 per-fire netCDF files as an old-format FINN text file.

The text file can be read by the Fortran ``fire_emis`` (grid_finn_fire_emis) program.
Input netCDF files are those written by ``calc_emis_daily.py``: one ``fire`` dimension,
variables date, lat, lon, area [m2], bmass [kg/m2], genveg, polyid, fireid, and one
variable per species. Here the species (gases AND aerosols) are in mol/day.

Output header (comma separated, one fire per line, in the order of the netCDF ``fire``
dimension - nothing is sorted):

    DAY,POLYID,FIREID,GENVEG,LATI,LONGI,AREA,BMASS,<species ...>

  * DAY is the day of the year. The Fortran code requires DAY to be non-decreasing down the
    file and takes the year from its start_date, so one text file holds one year. Each
    netCDF file is expected to hold a single day; if several files are given they are
    appended in chronological order (by their dates), and the script stops with an error
    if the rows would not be in non-decreasing day order.
  * AREA is in m2. BMASS (biomass burned, kg/m2) is written as an ordinary column; the
    Fortran code treats any column after AREA as a species, so BMASS is simply never
    mapped to an output variable.
  * Gases stay in mol/day. Species named with --aerosols are converted to kg/day:
        kg/day = mol/day * MW[g/mol] / 1000
    MW comes from each variable's ``molecular_weight_g_per_mol`` attribute, or from
    --mw NAME=VALUE. (A MW of 1 is the old "this is an aerosol" flag, not a real
    weight, so it is rejected.)

If the input spans more than one year, one file per year is written, named
<output stem>_<YYYY><suffix>.

Example
-------
python finn_nc_to_txt.py '/data/FINNv2.9.2_modvrs_MOZART-NEIVA_2018*.nc' \\
    -o FINNv2.9.2_2018.txt --aerosols BC OC
"""
from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

log = logging.getLogger("finn_nc_to_txt")

NON_SPECIES = {"day", "date", "fireid", "polyid", "genveg", "lat", "lon", "area", 
               "bmass", "frp"}
FIXED_COLUMNS = ["DAY", "FIREID", "POLYID", "GENVEG", "LATI", "LONGI", "AREA", "BMASS", "FRP"]
MAX_NAME_LEN = 16            # length of the species-name variable in the Fortran code


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+", help="FINNv2.9.2 netCDF file(s) or glob(s)")
    p.add_argument("-o", "--output", required=True, help="output text file")
    p.add_argument("--aerosols", nargs="+", default=[], metavar="SPECIES",
                   help="species in mol/day that are aerosols and must be written in kg/day, "
                        "e.g. --aerosols BC OC")
    p.add_argument("--mw", action="append", default=[], metavar="NAME=VALUE",
                   help="molecular weight (g/mol) for an aerosol, overriding the file's "
                        "molecular_weight_g_per_mol attribute (repeatable)")
    p.add_argument("--species", nargs="+", default=None, metavar="SPECIES",
                   help="write only these species, in this order (default: all, in file order)")
    p.add_argument("--exclude", nargs="+", default=[], metavar="SPECIES",
                   help="species to leave out")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def expand(patterns):
    files = []
    for pat in patterns:
        hits = sorted(glob.glob(pat))
        if not hits:
            raise SystemExit(f"no files match {pat!r}")
        files += hits
    return sorted(set(files))


def species_in_file(nc) -> list[str]:
    """Per-fire variables other than the fixed columns, in file order."""
    return [n for n, v in nc.variables.items()
            if v.dimensions == ("fire",) and n.lower() not in NON_SPECIES]


def pick(names, wanted, what):
    """Resolve user-given names against the file's, exact first then case-insensitive."""
    out = []
    for w in wanted:
        if w in names:
            out.append(w)
            continue
        ci = [n for n in names if n.lower() == w.lower()]
        if len(ci) != 1:
            raise SystemExit(f"{what} {w!r} not found in the netCDF files; "
                             f"species available: {names}")
        out.append(ci[0])
    return out


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(levelname)s :: %(message)s")

    files = expand(args.inputs)

    # ---- species list, units and molecular weights from the first file ----------
    with netCDF4.Dataset(files[0]) as nc:
        for need in ("date", "lat", "lon", "area", "genveg", "polyid", "fireid", "bmass"):
            if need not in nc.variables:
                raise SystemExit(f"{files[0]}: variable {need!r} not found - is this a "
                                 "FINNv2.9.2 per-fire file?")
        all_species = species_in_file(nc)
        units = {s: str(getattr(nc[s], "units", "")) for s in all_species}
        file_mw = {s: getattr(nc[s], "molecular_weight_g_per_mol", None) for s in all_species}

    species = all_species if args.species is None else pick(all_species, args.species, "species")
    exclude = set(pick(all_species, args.exclude, "--exclude species"))
    species = [s for s in species if s not in exclude]
    if not species:
        raise SystemExit("no species left to write")
    too_long = [s for s in species if len(s) > MAX_NAME_LEN]
    if too_long:
        log.warning("species names longer than %d characters are truncated by the Fortran "
                    "program and will not match its maps: %s", MAX_NAME_LEN, too_long)

    aerosols = set(pick(all_species, args.aerosols, "--aerosols species"))
    not_written = aerosols - set(species)
    if not_written:
        log.warning("--aerosols species not being written: %s", sorted(not_written))
    mw_over = {}
    for item in args.mw:
        name, sep, val = item.partition("=")
        if not sep:
            raise SystemExit(f"--mw expects NAME=VALUE, got {item!r}")
        mw_over[pick(all_species, [name], "--mw species")[0]] = float(val)

    factor = np.ones(len(species))                 # multiplier applied to each species column
    for k, s in enumerate(species):
        if s in aerosols:
            if units[s].lower().startswith("kg"):
                log.warning("%s: units are already %r; written without conversion", s, units[s])
                continue
            mw = mw_over.get(s, file_mw.get(s))
            if mw is None:
                raise SystemExit(f"aerosol {s}: no molecular_weight_g_per_mol attribute in the "
                                 f"file; give --mw {s}=<g/mol>")
            mw = float(mw)
            if mw == 1.0:
                raise SystemExit(f"aerosol {s}: molecular weight is 1.0 (the old aerosol flag); "
                                 f"give the real weight with --mw {s}=<g/mol>")
            factor[k] = mw / 1000.0
            log.info("aerosol %-10s mol/day -> kg/day with MW = %g g/mol", s, mw)
        elif units[s] and not units[s].lower().startswith("mol"):
            log.warning("%s: units are %r, not mol/day; written unchanged", s, units[s])

    # ---- order the files by their first date, then stream them to the output ----------
    info = []
    for f in files:
        with netCDF4.Dataset(f) as nc:
            missing = set(species) - set(species_in_file(nc))
            if missing:
                raise SystemExit(f"{f}: species {sorted(missing)} found in {files[0]} are "
                                 "missing from this file")
            d = np.asarray(nc["date"][:]).astype(np.int64)
        if d.size:
            info.append((int(d.min()), f))
    info.sort()

    out = Path(args.output)
    header = ",".join(FIXED_COLUMNS + species) + "\n"
    paths: dict[int, Path] = {}
    last_day: dict[int, int] = {}
    nrows: dict[int, int] = {}

    def path_for(year: int, multi: bool) -> Path:
        if year not in paths:
            paths[year] = (out.with_name(f"{out.stem}_{year}{out.suffix}") if multi else out)
            paths[year].parent.mkdir(parents=True, exist_ok=True)
            paths[year].write_text(header)
            last_day[year] = 0
            nrows[year] = 0
        return paths[year]

    # find all years first (cheap) so the naming rule is known before writing
    all_years: set[int] = set()
    for _, f in info:
        with netCDF4.Dataset(f) as nc:
            all_years |= set(np.unique(np.asarray(nc["date"][:]).astype(np.int64) // 10000).tolist())
    multi = len(all_years) > 1
    if multi:
        log.info("input spans years %s: writing one file per year", sorted(all_years))

    for _, f in info:
        with netCDF4.Dataset(f) as nc:
            date = np.asarray(nc["date"][:]).astype(np.int64)
            if date.size == 0:
                continue
            def col(name, dtype):
                return np.ma.filled(nc[name][:], 0).astype(dtype)
            data = {
                "date": date,
                "POLYID": col("polyid", np.int64),
                "FIREID": col("fireid", np.int64),
                "GENVEG": col("genveg", np.int64),
                "LATI": col("lat", np.float64),
                "LONGI": col("lon", np.float64),
                "AREA": col("area", np.float64),
                "BMASS": col("bmass", np.float64),
                "FRP": col("frp", np.float64),
            }
            for k, s in enumerate(species):
                data[s] = np.nan_to_num(np.ma.filled(nc[s][:].astype(np.float64), 0.0)) * factor[k]
        df = pd.DataFrame(data)
        # day of year from the calendar date
        dts = pd.to_datetime(df["date"].astype(str), format="%Y%m%d")
        df.insert(0, "DAY", dts.dt.dayofyear.to_numpy(dtype=np.int64))
        df["year"] = (df["date"] // 10000).to_numpy()
        # rows stay in fire-dimension order; groupby only splits by year (sort=False keeps order)
        for year, chunk in df.groupby("year", sort=False):
            p = path_for(int(year), multi)
            days = chunk["DAY"].to_numpy()
            if days[0] < last_day[year] or (np.diff(days) < 0).any():
                raise SystemExit(f"{f}: day of year is not non-decreasing in the order written "
                                 f"(year {year}, after day {last_day[year]}); the Fortran program "
                                 "needs non-decreasing days. Files are expected to be single-day.")
            last_day[year] = int(days[-1])
            chunk[FIXED_COLUMNS + species].to_csv(
                p, mode="a", header=False, index=False, float_format="%.9g")
            nrows[year] += len(chunk)
        log.info("%s: %d fires", Path(f).name, len(df))

    for year, p in sorted(paths.items()):
        log.info("wrote %s (%d fires, year %d)", p, nrows[year], year)
    if not paths:
        raise SystemExit("no fires found in the input files")
    return 0


if __name__ == "__main__":
    sys.exit(main())

