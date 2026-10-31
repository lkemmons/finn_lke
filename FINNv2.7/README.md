# FINNv2.7

FINNv2.7 was created specifically for analysis of the ASIA-AQ field experiment
over Southeast Asia in Feb-Mar 2024 by Nattamon Maneenoi, UCLA.

Hourly emissions were created based on geostationary satellite fire detections.

A few tools were created to handle these unique files.

### Gridding and concatenating files

 
| Script                     | Input                                          | Output                                             |
|----------------------------|------------------------------------------------|----------------------------------------------------|
| `grid_txt_emissions.py`    | Legacy FINN text CSV (for FINNv2.7-hourly)     | Gridded NetCDF, `(time, lat, lon)` or `(time, ncol)` |
| `concat_finn_daily_hourly.py` | Prepend 24 hourly copies of one daily-average day to an hourly file |
| `append_daily_to_hourly.py`   | Append 24 hourly copies of one daily-average day to an hourly file  |

Together they let you extend an hourly-cadence file by one day at
either end, using data from a daily-average source file — while
keeping uniform hourly cadence throughout (developed for FINNv2.7-hourly).


### Gridding v2.7-hourly txt files

For legacy FINN text CSVs (needs `--mw-table` for unit annotation):
```bash
python scripts/grid_txt_emissions.py \
    FINN_v2.7_MOZART_hourly_2024033.txt \
    --mw-table scripts/finn_text_species_mw.csv \
    --grid-resolution 0.1 0.1 --out-dir ./gridded
```

Produces one NetCDF per species with **CAM/CESM emission-input
conventions**: single unlimited `time` dim, plus `date(time)` (int32
YYYYMMDD) and `datesec(time)` (int32 seconds-of-day) companion
variables.  Units default to `molecules/cm2/s` for gases and aerosols,
`number/cm2/s` for particle-count species.  Multi-day concatenation is
a plain `xr.concat` along `time`.


### Extend an hourly file with one more day

For FINNv2.7-hourly: 
Both tools pick one day from a daily-average source file, replicate
its spatial pattern across 24 hourly slots, rebase onto the hourly
file's epoch, then concat.  Result: uniform hourly cadence throughout.

Prepend the day before the hourly file's first date:
```bash
python scripts/concat_finn_daily_hourly.py \
    --daily-dir  /glade/.../daily_avg    \
    --hourly-dir /glade/.../hourly_avg   \
    --out-dir    /glade/.../concat
# default seam date = day BEFORE hourly's first date
```

Append the day after the hourly file's last date:
```bash
python scripts/append_daily_to_hourly.py \
    --hourly-dir /glade/.../concat       \
    --daily-dir  /glade/.../daily_avg    \
    --out-dir    /glade/.../extended
# default seam date = day AFTER hourly's last date
```

Both handle:
- **Cross-epoch time units** (e.g. daily `days since 1950-01-01` +
  hourly `days since 1970-01-01`): daily is rebased onto hourly's
  epoch automatically.
- **Cross-calendar** (gregorian / standard / proleptic_gregorian
  treated as equivalent; other calendars refused with a clear error).
- **Generic emission variable names**: some FINN files call the
  emission variable `fire` with the actual species in a `long_name`
  or `map` attribute; species pairing across a directory of `fire`
  variables works via those attributes.
- **Coord/data-var inconsistencies**: `lat` and `lon` are promoted to
  coordinates before concat, so files that declare
  `fire:coordinates = "lat lon"` mix cleanly with files that don't.

## Legacy text-CSV conventions

Legacy FINN text files carry two conventions that need special
handling:

### 1. HOUR values are unreliable (FINNv2.7-hourly)

The FINNv2.7 hourly text outputs have the `HOUR` column
but its values don't match the actual UTC hour of each row.
`datetimeUTC` is trustworthy (ISO timestamp).
So the design for `grid_txt_emissions.py` is:

- **HOUR column presence** = signal that the file is *intended* to be
  hourly (stable structural fact).
- **HOUR values** = ignored.  Actual UTC hour comes from
  `pd.to_datetime(df["datetimeUTC"]).dt.hour`.
- **HOUR column** = dropped before writing, so it can't leak into
  outputs.


### 2. Species aren't in variable names

The gridded (from txt) FINNv2.5 or v2.7 files (both daily and hourly)
call the emission variable `fire` and put the actual species in a
variable attribute:

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

