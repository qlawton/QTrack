#!/usr/bin/env python

# # Example: tracking easterly waves globally
#
# The same four steps as `example_tracking.py`, but on a global grid and with wave
# genesis allowed in every basin instead of only over Africa. Longitude is treated
# as periodic throughout, so a wave is followed continuously across the
# international dateline rather than being lost there.
#
# Requirements are the same as the standard workflow -- 700 hPa u and v wind,
# 6-hourly, 1x1 degree -- except that the grid must span all 360 degrees of
# longitude. `prep_data` warns if the longitude axis it produces is not contiguous.

import time

import qtrack
from qtrack.curvvort import compute_curvvort
from qtrack.tracking import run_postprocessing, run_tracking

tstart = time.time()

# ### NAME FILES
prepped_data_save = "adjusted_data_global.nc"
curv_file_out = "curv_vort_global.nc"
AEW_raw_save_file = "AEW_tracks_raw_global.nc"
AEW_final_nc_file = "AEW_tracks_post_processed_global.nc"
AEW_final_obj_file = "AEW_tracks_post_processed_global.pkl"
year_in = 2010

# ### DOWNLOAD EXAMPLE DATA (OPTIONAL)
# A 10-day global ERA5 period in September 2010. `era5_2010` (July-September) works
# too and gives waves long enough to travel much further.
qtrack.download_examples("era5_2010_10day", "")

# ### Prep data
qtrack.prep_data(data_in="era5_700_wind_global_2010_10day.nc", data_out=prepped_data_save)

# ### Curvature vorticity calculation
# On a global grid the winds are wrapped across the seam before the zonal gradient
# and the radial average, so there is no hole in the field at 180.
compute_curvvort(prepped_data_save, curv_file_out, njobs_in=-1)

# ### Tracking step
# initiation_bounds=None allows a wave to be initiated anywhere. The default,
# (-35, 40), is the classic African easterly wave genesis region -- keep it if you
# only want AEWs, and they will still be followed westward around the globe.
#
# left_right_bounds defaults to None (the whole globe). regions defaults to
# qtrack.regions.global_regions(), which keeps the Africa-through-Caribbean sector
# exactly as it has always been and applies general ocean settings elsewhere.
run_tracking(
    input_file=curv_file_out,
    save_file=AEW_raw_save_file,
    initiation_bounds=None,
)

# ### Postprocessing step
# Pass the same `regions` value used above. hov_lon_limits=None (the default) makes
# the Hovmoller span the data; set it to (-100, 40) for the classic Atlantic view.
run_postprocessing(
    input_file=AEW_raw_save_file,
    real_year_used=year_in,
    curv_data_file=curv_file_out,
    save_obj_file=AEW_final_obj_file,
    save_nc_file=AEW_final_nc_file,
    hov_save_file="final_hovmoller_global.png",
)

# ### WHICH TRACKS CROSSED THE DATELINE?
# AEW_lon_unwrapped is AEW_lon with the 360 degree jumps removed, so a crossing
# track is a straight line rather than a discontinuity. Its values may lie outside
# -180 to 180.
import numpy as np  # noqa: E402
import xarray as xr  # noqa: E402

tracks = xr.open_dataset(AEW_final_nc_file)
lon = tracks["AEW_lon"].values
unwrapped = tracks["AEW_lon_unwrapped"].values
for system in range(lon.shape[0]):
    valid = ~np.isnan(lon[system])
    if valid.sum() < 2:
        continue
    if np.abs(np.diff(lon[system][valid])).max() > 180:
        first, last = unwrapped[system][valid][[0, -1]]
        print(f"System {system + 1} crossed 180: {first:.1f} to {last:.1f} degrees unwrapped ({first - last:.0f} degrees of westward travel)")

# ### OUTPUT TIME
tend = time.time()
print("Time to run computation: " + str(round(tend - tstart, 2)) + " seconds | " + str(round((tend - tstart) / 60, 2)) + " minutes")
print()
