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
