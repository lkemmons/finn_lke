#!/usr/bin/env python3
"""Calculate FINNv2.9.2 NRT daily emissions from finn_py fire output.

Direct port of ``finn2_9_2_calc_emis_nrt_daily.pro`` (IDL) to Python.
The single intentional change vs. the IDL is the output format: NetCDF
instead of a CSV text file.  All algorithmic logic — genveg
assignment, VCF cleanup, fuel-load lookups, biomass-burned formulas,
the gas-vs-aerosol unit conversion — is preserved exactly.

Usage
-----
    python calc_emis_daily.py YYYYJJJ YEAR_RST [--options]

Examples
--------
    python calc_emis_daily.py 2026150 2026

    python calc_emis_daily.py 2026150 2026 \\
        --path-inputs   /glade/work/emmons/FINN_python/FINNv2.9nrt/emis_calc/finn_inputs \\
        --path-in       /glade/derecho/scratch/emmons/finn2.9nrt_output \\
        --path-out      /glade/derecho/scratch/emmons/finn2.9.2nrt_emis \\
        --tag-fire      v2.9nrt \\
        --finn-version  v2.9.2 \\
        --sim-id        NRTmodvrs

Inputs read
-----------
``<path-inputs>/Fuel_LOADS_NEW_022019.csv``
    Per-global-region fuel loads (g/m²) for 5 veg classes.
``<path-inputs>/LCTFuelLoad_fuel4_revisit20190521.csv``
    US-specific per-LCT tree + herb fuel loads, used when
    ``v_regnum == 1`` (North America).
``<path-inputs>/EFs_byGenVeg_NEIVA_MOZART_c20260209.csv``
    Emission factors (g of species per kg dry biomass) for genveg 1-6,9.
``<path-in>/out_<tag-fire>_<date_lab>_modlct_<year_rst>_modvcf_<year_rst>_regnum.csv``
    finn_py polygon output for the requested day.

Output written
--------------
``<path-out>/FINN<finnver>_<simid>_MOZART_<YYYYMMDD>.nc``
    NetCDF file with one ``fire`` dimension; per-fire columns (day,
    fireid, polyid, genveg, lat, lon, area, bmass, frp) plus one
    variable per chemical species.
``<path-out>/logs/LOG_calcemis_FINN<finnver>_<simid>_<YYYYMMDD>_<creation>.txt``
    Log of skipped fires and processing summary, matching the IDL log
    format.
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

log = logging.getLogger("calc_emis_daily")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="calc_emis_daily.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("date_lab", type=str,
                   help="emissions date as YYYYJJJ (e.g. 2026150)")
    p.add_argument("year_rst", type=str,
                   help="raster year used in the input filename "
                        "(e.g. 2026)")
    p.add_argument("--path-inputs", type=Path,
                   default=Path("/glade/work/emmons/FINN_python/"
                                "FINNv2.9nrt/emis_calc/finn_inputs"),
                   help="directory with fuel loads + EFs CSV files")
    p.add_argument("--path-in", type=Path,
                   default=Path("/glade/derecho/scratch/emmons/"
                                "finn2.9nrt_output"),
                   help="directory with finn_py polygon CSV outputs")
    p.add_argument("--path-out", type=Path,
                   default=Path("/glade/derecho/scratch/emmons/"
                                "finn2.9.2nrt_emis"),
                   help="output directory (NetCDF + logs)")
    p.add_argument("--tag-fire", default="v2.9nrt",
                   help="tag in the finn_py input filename "
                        "(default: %(default)s)")
    p.add_argument("--finn-version", default="v2.9.2",
                   help="FINN version label for the output filename")
    p.add_argument("--sim-id", default="NRTmodvrs",
                   help="simulation-id label for the output filename")
    p.add_argument("--file-label", default="MOZART",
                   help="species-set label for the output filename "
                        "(matches the EF table)")
    p.add_argument("--file-fuelloads", default="Fuel_LOADS_NEW_022019.csv")
    p.add_argument("--file-usfuel",    default="LCTFuelLoad_fuel4_revisit20190521.csv")
    p.add_argument("--file-efs",       default="EFs_byGenVeg_NEIVA_MOZART_c20260209.csv")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Input-table readers
# ---------------------------------------------------------------------------

def read_fuel_loads(path: Path) -> dict[str, np.ndarray]:
    """Read per-region fuel loads (g/m²).  Mirrors the IDL `read_csv` of
    Fuel_LOADS_NEW_022019.csv.

    Returns a dict with keys: ``region``, ``tf``, ``te``, ``bf``,
    ``ws``, ``gr`` — each a 1-D array indexed in the same order the
    file lists the regions (1..N).  The IDL code uses ``ireg = globreg - 1``
    to index into these; we preserve that convention.
    """
    df = pd.read_csv(path)
    log.info("  %s contains: %s", path.name, list(df.columns))
    if df.columns[0] != "GlobalRegion" or df.columns[5] != "SavannaGrasslands":
        raise RuntimeError(f"{path}: unexpected column layout {list(df.columns)}")
    return {
        "region": df.iloc[:, 0].to_numpy(),
        "tf":     df.iloc[:, 1].to_numpy(dtype="float64"),  # tropical forest
        "te":     df.iloc[:, 2].to_numpy(dtype="float64"),  # temperate forest
        "bf":     df.iloc[:, 3].to_numpy(dtype="float64"),  # boreal forest
        "ws":     df.iloc[:, 4].to_numpy(dtype="float64"),  # woody savanna
        "gr":     df.iloc[:, 5].to_numpy(dtype="float64"),  # grass/savanna
    }


def read_us_fuel_loads(path: Path) -> dict[int, tuple[float, float]]:
    """Read per-LCT US fuel loads.  Returns ``{lct: (tree, herb)}``."""
    df = pd.read_csv(path)
    log.info("  %s contains: %s", path.name, list(df.columns))
    if df.columns[2].upper() != "HERB":
        raise RuntimeError(f"{path}: expected HERB as 3rd column, got {df.columns[2]}")
    return {
        int(row.iloc[0]): (float(row.iloc[1]), float(row.iloc[2]))
        for _, row in df.iterrows()
    }


def read_emission_factors(path: Path) -> dict:
    """Read the EFs_byGenVeg CSV.

    Layout (matches IDL parser):
      row 1: header / title (ignored)
      row 2: column names — first 2 cols are GenVegType + GenVegDescript,
             remaining cols are species
      row 3: molecular weights (g/mol) on the species columns;
             MW=1 is a flag for "aerosol — keep in kg, no mol conversion"
      rows 4+: one row per genveg type, with the EF values in g/kg-dm

    Returns a dict with: ``species`` (np.ndarray[str]), ``mws``
    (np.ndarray[float]), ``genveg`` (np.ndarray[int]), ``ef``
    (2-D array shape (ntype, nspec)).
    """
    with path.open() as f:
        lines = [ln.rstrip("\n") for ln in f]
    header = lines[1].split(",")           # column names
    species = np.array(header[2:])
    mws = np.array(lines[2].split(",")[2:], dtype="float64")
    genveg, ef_rows = [], []
    for ln in lines[3:]:
        if not ln.strip():
            continue
        parts = ln.split(",")
        genveg.append(int(parts[0]))
        ef_rows.append([float(x) for x in parts[2:]])
    ef = np.array(ef_rows, dtype="float64")
    log.info("  EF species: %s", list(species))
    log.info("  EF genveg types: %s", genveg)
    return {
        "species": species,
        "mws":     mws,
        "genveg":  np.array(genveg, dtype="int32"),
        "ef":      ef,
    }


# ---------------------------------------------------------------------------
# Per-fire processing (genveg, VCF adjust, biomass burned)
# ---------------------------------------------------------------------------

# LCT codes that are non-burnable; remove fires assigned to them.
# The IDL filters: lct >= 17, lct <= 0, lct == 15.
# (15 = permanent snow/ice; 17 = water; 0/255 = unclassified.)
_INVALID_LCT = lambda lct: (lct >= 17) or (lct <= 0) or (lct == 15)


def _assign_genveg_and_lct(lct: int, lat: float, tree: float) -> tuple[int, int]:
    """Return ``(genveg, lct)`` for one fire.

    Faithful port of the IDL ``case lct of …`` block.  Note that the
    Urban branch (lct=13) mutates ``lct`` in addition to setting
    ``genveg``, which matters because the US fuel-load lookup below
    indexes into ``lcttree[lct]`` / ``lctherb[lct]``.
    """
    if lct == 1:
        return (5 if lat > 50 else 6, lct)
    if lct == 2:
        return (3 if -23.5 <= lat <= 23.5 else 4, lct)
    if lct == 3:
        return (5 if lat > 50 else 4, lct)
    if lct == 4:
        return (4, lct)
    if lct == 5:
        if lat > 50:
            return (5, lct)
        return (3 if -23.5 <= lat <= 23.5 else 4, lct)
    if lct in (6, 7, 8):
        return (2, lct)
    if lct in (9, 10, 11):
        return (1, lct)
    if lct == 12:
        return (9, lct)
    if lct == 13:                                     # Urban — reclassify
        if tree < 40.0:
            return (1, 10)                             # grassland
        if tree < 60.0:
            return (2, 8)                              # woody savanna
        # tree >= 60 → forest, lct → mixed forest
        if lat > 50:
            return (5, 1)                              # boreal, evergreen needleleaf
        return ((3 if -30 <= lat <= 30 else 4), 5)     # tropical or temperate forest
    if lct == 14:
        return (1, lct)
    if lct == 16:
        return (1, lct)
    return (-1, lct)


def _adjust_vcf(lct: int, tree: float, herb: float, bare: float
                ) -> tuple[float, float, float, float, str | None]:
    """Apply the IDL's VCF cleanup steps.

    Returns ``(tree, herb, bare, totcov, log_msg_or_None)``.
    Two adjustments happen here:

    1. Scale tree/herb/bare to sum to 100 if they don't already (within
       a small tolerance).
    2. Reassign cover values when ``bare >= 99.9`` (no VCF information
       — guess from LCT).
    """
    # Clip negatives (legacy -9999 fill values)
    if tree < 0.: tree = 0.
    if herb < 0.: herb = 0.
    if bare < 0.: bare = 0.

    totcov = tree + herb + bare
    msg = None

    # Scale to 100
    if totcov > 101. or totcov < 99.:
        if totcov > 0:
            totcov_orig = totcov
            tree = tree * 100. / totcov
            herb = herb * 100. / totcov
            bare = bare * 100. / totcov
            totcov = tree + herb + bare
            msg = f"totcov adjusted from {totcov_orig:.1f} to {totcov:.1f}"

    # 100% bare → guess from LCT
    if bare >= 99.9:
        if lct <= 5:
            tree, herb, bare = 60., 40., 0.
        elif (6 <= lct <= 8) or (lct == 11) or (lct == 14):
            tree, herb, bare = 50., 50., 0.
        elif lct in (9, 10, 12, 13, 16):
            tree, herb, bare = 20., 80., 0.
        msg = (msg + "; " if msg else "") + (
            f"100% bare reassigned by LCT={lct} to T/H/B={tree:.0f}/"
            f"{herb:.0f}/{bare:.0f}")

    return tree, herb, bare, totcov, msg


def _compute_bmass(tree: float, herb: float, bare: float,
                    coarsebm: float, herbbm: float) -> float:
    """Biomass burned (g-dm/m²).  Follows IDL's three-branch formula."""
    if tree > 60:                              # FOREST
        CF1, CF3 = 0.30, 0.90
        return ((herb / 100.) * herbbm * CF3
                + (tree / 100.) * (herbbm * CF3 + coarsebm * CF1))
    if tree > 40:                              # WOODLAND
        CF3 = np.exp(-0.013 * tree)
        CF1 = 0.30
        return ((herb / 100.) * herbbm * CF3
                + (tree / 100.) * (herbbm * CF3 + coarsebm * CF1))
    # GRASSLAND (tree <= 40)
    CF3 = 0.98
    return ((herb / 100.) * herbbm * CF3
            + (tree / 100.) * herbbm * CF3)


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------

def process_fires(
    df: pd.DataFrame,
    fuel: dict, us_fuels: dict, year_emis: int,
    log_fp,
) -> pd.DataFrame:
    """Process the finn_py polygon CSV into per-emission rows.

    Returns a DataFrame with columns: day, date, polyid, fireid, lat,
    lon, area, bmass, genveg, frp — one row per surviving fire (i.e.,
    one per LCT-class within each polygon that passes all filters).
    """
    n_in = len(df)
    log.info("# input rows: %d", n_in)
    print(f"{n_in} input fires", file=log_fp)

    out_rows: list[dict] = []
    iskip_yr  = 0
    iskip_reg = 0

    # Pull the columns we need
    have_regnum = "v_regnum" in df.columns
    for i, row in enumerate(df.itertuples(index=False)):
        # --- date ---
        try:
            yy, mm, dd = (int(s) for s in str(row.acq_date_utc).split("-"))
        except Exception as e:
            print(f"input row {i} unparseable date {row.acq_date_utc!r}: {e}",
                  file=log_fp)
            continue
        if yy != year_emis:
            iskip_yr += 1
            continue
        date_int = yy * 10000 + mm * 100 + dd
        jday = (dt.date(yy, mm, dd) - dt.date(yy, 1, 1)).days + 1

        # --- core polygon attrs ---
        polyid   = int(row.polyid)
        fireid   = int(row.fireid)
        lon      = float(row.cen_lon)
        lat      = float(row.cen_lat)
        area     = float(row.area_sqkm)
        # alg_agg available as row.alg_agg if needed
        lct      = int(row.v_lct)  if not pd.isna(row.v_lct)  else 0
        flct     = float(row.f_lct) if not pd.isna(row.f_lct) else 0.
        tree     = float(row.v_tree) if not pd.isna(row.v_tree) else 0.
        herb     = float(row.v_herb) if not pd.isna(row.v_herb) else 0.
        bare     = float(row.v_bare) if not pd.isna(row.v_bare) else 0.
        globreg  = int(row.v_regnum) if (have_regnum and not pd.isna(row.v_regnum)) else 0
        frp      = float(row.v_frp)  if not pd.isna(row.v_frp)  else 0.

        # --- region filter ---
        if globreg < 1 or globreg > 12:
            print(f"Fire {i} removed. global region: {globreg} "
                  f"lon, lat: {lon:.1f} {lat:.1f}", file=log_fp)
            iskip_reg += 1
            continue

        # --- LCT filter ---
        if _INVALID_LCT(lct):
            print(f"Fire {i} removed: lct={lct}", file=log_fp)
            continue

        # --- VCF cleanup ---
        tree, herb, bare, totcov, vcf_msg = _adjust_vcf(lct, tree, herb, bare)
        if totcov >= 240. or totcov < 1.:
            print(f"Fire {i} removed. totcov={totcov:.0f}", file=log_fp)
            continue
        if vcf_msg:
            print(f"Fire {i}: {vcf_msg}", file=log_fp)

        # --- genveg + (maybe) lct reassignment ---
        genveg, lct = _assign_genveg_and_lct(lct, lat, tree)
        if genveg <= 0:
            print(f"Fire {i} no genveg set. lat,lon,LCT = "
                  f"{lat:.3f},{lon:.3f},{lct}", file=log_fp)
            continue

        # --- fuel loads ---
        ireg = globreg - 1
        bmass1 = -1.
        if genveg == 1: bmass1 = fuel["gr"][ireg]
        elif genveg == 2: bmass1 = fuel["ws"][ireg]
        elif genveg == 3: bmass1 = fuel["tf"][ireg]
        elif genveg == 4: bmass1 = fuel["te"][ireg]
        elif genveg == 5: bmass1 = fuel["bf"][ireg]
        elif genveg == 6: bmass1 = fuel["te"][ireg]
        elif genveg == 9: bmass1 = 902.
        # Boreal in Southern Asia (region 11) → use temperate
        if genveg == 5 and globreg == 11:
            bmass1 = fuel["te"][ireg]

        if bmass1 < 0:
            print(f"BMASS1<0: Fire {i} removed. genveg={genveg} "
                  f"globreg={globreg} ireg={ireg}", file=log_fp)
            continue

        coarsebm = bmass1
        herbbm   = fuel["gr"][ireg]
        # North America: override with LCT-specific fuels
        if globreg == 1:
            tree_us, herb_us = us_fuels.get(lct, (coarsebm, herbbm))
            coarsebm, herbbm = tree_us, herb_us

        # --- biomass burned ---
        bmass = _compute_bmass(tree, herb, bare, coarsebm, herbbm)
        bmass = bmass / 1000.                            # g/m² → kg/m²

        # --- effective burned area (m²) ---
        areanow = area * flct * 1e6                      # km² → m²
        area_bare = areanow * (bare / 100.)
        areanow = areanow - area_bare
        if areanow < 1.:
            print(f"area=0. area,flct,bare: {areanow:.3f} {flct} {bare}",
                  file=log_fp)
            continue

        out_rows.append({
            "day":    jday,
            "date":   date_int,
            "polyid": polyid,
            "fireid": fireid,
            "lat":    lat,
            "lon":    lon,
            "area":   areanow,
            "bmass":  bmass,
            "genveg": genveg,
            "frp":    frp,
        })

    log.info("kept %d / %d fire-rows", len(out_rows), n_in)
    print(f"# fires skipped because wrong year: {iskip_yr}", file=log_fp)
    print(f"# fires skipped because no region assigned: {iskip_reg}", file=log_fp)
    print(f"# fire-rows with emissions: {len(out_rows)}", file=log_fp)
    if n_in:
        print(f"% of total fire-rows saved: {100.*len(out_rows)/n_in:.1f}",
              file=log_fp)

    out = pd.DataFrame(out_rows)
    if len(out):
        out = out.sort_values("day").reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# Vectorised emission calculation + NetCDF writer
# ---------------------------------------------------------------------------

def compute_emissions(processed: pd.DataFrame, ef: dict, log_fp
                       ) -> np.ndarray:
    """Compute per-fire per-species emissions.

    Returns shape ``(n_fires, n_species)``.  Units: ``mol/fire/day`` for
    gas-phase species (MW != 1), ``kg/fire/day`` for aerosols
    (MW == 1, the placeholder used in the EF table).
    """
    n = len(processed)
    nsp = len(ef["species"])
    if n == 0:
        return np.zeros((0, nsp), dtype="float32")

    # Map each fire's genveg to a row index into ef["ef"].
    genveg_to_row = {int(g): i for i, g in enumerate(ef["genveg"])}
    rows = []
    drop_mask = np.zeros(n, dtype=bool)
    for k, g in enumerate(processed["genveg"].astype(int).to_numpy()):
        idx = genveg_to_row.get(int(g))
        if idx is None:
            print(f"no EF for genveg={g}", file=log_fp)
            drop_mask[k] = True
            rows.append(0)
        else:
            rows.append(idx)
    rows = np.array(rows, dtype="int64")
    ef_per_fire = ef["ef"][rows, :]                       # (n, nsp)

    area  = processed["area"].to_numpy(dtype="float64")[:, None]    # (n,1)
    bmass = processed["bmass"].to_numpy(dtype="float64")[:, None]   # (n,1)

    # For gases (MW != 1):
    #   emis = EF * 1e-3 * area * bmass / (MW * 1e-3)
    # For aerosols (MW == 1):
    #   emis = EF * 1e-3 * area * bmass            (no MW conversion)
    denom = np.where(ef["mws"] != 1.0, ef["mws"] * 1e-3, 1.0)        # (nsp,)
    emis = ef_per_fire * 1e-3 * area * bmass / denom[None, :]        # (n, nsp)

    # Zero-out rows we couldn't find an EF for; the variable values
    # become 0, but we keep the row so the per-fire indices align with
    # the polygon/fireid columns.
    if drop_mask.any():
        emis[drop_mask, :] = 0.
        log.warning("zeroed %d / %d fires that had no matching EF", drop_mask.sum(), n)

    return emis.astype("float32")


def write_netcdf(processed: pd.DataFrame, emis: np.ndarray, ef: dict,
                  out_path: Path, *,
                  date_processed: str, input_file: Path,
                  finn_version: str, sim_id: str, file_label: str) -> None:
    """Write the output NetCDF with one variable per species."""
    ds = xr.Dataset(
        data_vars={
            "day":    (("fire",), processed["day"].to_numpy(dtype="int16"),
                        {"long_name": "julian day of year (1..366)"}),
            "date":   (("fire",), processed["date"].to_numpy(dtype="int32"),
                        {"long_name": "calendar date as YYYYMMDD"}),
            "polyid": (("fire",), processed["polyid"].to_numpy(dtype="int64"),
                        {"long_name": "polygon id from finn_py output"}),
            "fireid": (("fire",), processed["fireid"].to_numpy(dtype="int64"),
                        {"long_name": "fire cluster id from finn_py"}),
            "lat":    (("fire",), processed["lat"].to_numpy(dtype="float32"),
                        {"units": "degrees_north"}),
            "lon":    (("fire",), processed["lon"].to_numpy(dtype="float32"),
                        {"units": "degrees_east"}),
            "area":   (("fire",), processed["area"].to_numpy(dtype="float32"),
                        {"units": "m2",
                         "long_name": "burned area (vegetation only, "
                                       "bare fraction removed)"}),
            "bmass":  (("fire",), processed["bmass"].to_numpy(dtype="float32"),
                        {"units": "kg/m2",
                         "long_name": "biomass burned per unit area "
                                       "(dry matter)"}),
            "genveg": (("fire",), processed["genveg"].to_numpy(dtype="int16"),
                        {"long_name": "generic vegetation class (1..9)",
                         "flag_values": "1 2 3 4 5 6 9",
                         "flag_meanings": "grassland shrub tropical_forest "
                                          "temperate_forest boreal_forest "
                                          "temperate_evergreen crop_generic"}),
            "frp":    (("fire",), processed["frp"].to_numpy(dtype="float32"),
                        {"units": "MW",
                         "long_name": "fire radiative power (mean over "
                                       "polygon's source detections)"}),
        },
        attrs={
            "title":           f"FINN{finn_version} {sim_id} daily emissions",
            "source":          "computed from finn_py polygon output",
            "input_file":      str(input_file),
            "finn_version":    finn_version,
            "sim_id":          sim_id,
            "species_set":     file_label,
            "date_processed":  date_processed,
            "created":         dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "Conventions":     "CF-1.8",
        },
    )

    # One variable per species.  Distinguish gas vs aerosol by MW.
    for i, sp in enumerate(ef["species"]):
        mw = float(ef["mws"][i])
        if mw == 1.0:
            units = "kg/day"
            longname = f"emitted mass of {sp}"
        else:
            units = "mol/day"
            longname = f"emitted moles of {sp}"
        ds[sp] = (("fire",), emis[:, i],
                   {"units": units, "long_name": longname,
                    "molecular_weight_g_per_mol": mw})

    # Light compression on all numeric variables.
    encoding = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
    ds.to_netcdf(out_path, encoding=encoding, format="NETCDF4")
    log.info("wrote %s  (%d fire-rows × %d species)",
             out_path, len(processed), len(ef["species"]))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )

    # Date parsing  (YYYYJJJ → YYYY-MM-DD + YYYYMMDD label)
    year = int(args.date_lab[:4])
    doy = int(args.date_lab[4:7])
    target_date = dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)
    sdate_emis = target_date.strftime("%Y%m%d")
    log.info("=" * 70)
    log.info("FINN%s emissions for %s (%s)  sim=%s",
             args.finn_version, args.date_lab, target_date.isoformat(), args.sim_id)
    log.info("=" * 70)
    t0 = time.time()

    # Set up outputs
    args.path_out.mkdir(parents=True, exist_ok=True)
    (args.path_out / "logs").mkdir(parents=True, exist_ok=True)
    today = dt.date.today().strftime("%Y%m%d")
    logfile = args.path_out / "logs" / (
        f"LOG_calcemis_FINN{args.finn_version}_{args.sim_id}_"
        f"{sdate_emis}_{today}.txt"
    )
    log.info("writing log file: %s", logfile)

    # Input filename
    file_in = args.path_in / (
        f"out_{args.tag_fire}_{args.date_lab}_modlct_{args.year_rst}_"
        f"modvcf_{args.year_rst}_regnum.csv"
    )
    log.info("reading: %s", file_in)
    if not file_in.exists():
        log.error("input file not found: %s", file_in)
        return 1

    # Output filename
    out_nc = args.path_out / (
        f"FINN{args.finn_version}_{args.sim_id}_{args.file_label}_"
        f"{sdate_emis}.nc"
    )

    with logfile.open("w") as log_fp:
        print(f"FINN{args.finn_version} emissions calc, started "
              f"{dt.datetime.utcnow().isoformat()}Z", file=log_fp)
        print(f"date: {target_date.isoformat()} ({args.date_lab})", file=log_fp)
        print(f"input: {file_in}", file=log_fp)
        print(f"output: {out_nc}", file=log_fp)
        print("", file=log_fp)

        # Read input tables
        log.info("reading fuel/EF tables from %s", args.path_inputs)
        fuel     = read_fuel_loads(args.path_inputs / args.file_fuelloads)
        us_fuels = read_us_fuel_loads(args.path_inputs / args.file_usfuel)
        ef       = read_emission_factors(args.path_inputs / args.file_efs)

        # Log EF summary as IDL did
        for i, sp in enumerate(ef["species"]):
            print(f"  {sp:>12s}  MW={ef['mws'][i]:>5.1f}  "
                  f"EFs={ef['ef'][:, i]}", file=log_fp)
        print("", file=log_fp)

        # Read polygon CSV
        df_in = pd.read_csv(file_in)
        log.info("read %d rows from %s", len(df_in), file_in.name)
        print(f"input file columns: {list(df_in.columns)}", file=log_fp)

        processed = process_fires(df_in, fuel, us_fuels, year, log_fp)

        # Emissions
        log.info("computing emissions for %d fire-rows × %d species",
                 len(processed), len(ef["species"]))
        emis = compute_emissions(processed, ef, log_fp)

        # NetCDF output
        write_netcdf(
            processed, emis, ef, out_nc,
            date_processed=target_date.isoformat(),
            input_file=file_in,
            finn_version=args.finn_version,
            sim_id=args.sim_id,
            file_label=args.file_label,
        )

        elapsed = time.time() - t0
        msg = f"Running time: {int(elapsed//60)}:{int(elapsed%60):02d}"
        log.info(msg)
        print(msg, file=log_fp)
        print(f"Completed at: {dt.datetime.utcnow().isoformat()}Z", file=log_fp)

    log.info("done — %s", out_nc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
