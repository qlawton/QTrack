#!/usr/bin/env python

# # Example: Using QTrack to tracking AEWs in reanalysis data
#
# This example will demonstrate the basic functions of AEW tracking from the Qtrack module. In order to track AEWs, we must proceed through a few steps:
#
# 1. Download or location data for tracking. Data requirments: 700hPa wind (u and v component), 6 hourly timesteps, 1 degree by 1 degree resolution.
# 2. `prep_data` function to prep the input data for tracking.
# 4. `curvvort` Compute the curvature vorticity from the wind fields, and then take radial averages at each gridpoint to smooth it.
# 5. `tracking` Run the AEW tracking code.
# 6. `postproc` Run postprocessing on AEW data, including the creation of a netCDF4 file.


import time

import qtrack
from qtrack.curvvort import compute_curvvort
from qtrack.tracking import run_postprocessing, run_tracking

tstart = time.time()


# ### NAME FILES
prepped_data_save = "adjusted_data.nc"
curv_file_out = "curv_vort_era5_test.nc"
AEW_raw_save_file = "AEW_tracks_raw.nc"
AEW_final_nc_file = "AEW_tracks_post_processed.nc"
AEW_final_obj_file = "AEW_tracks_post_processed.pkl"
year_in = 2021

# ### DOWNLOAD EXAMPLE DATA (OPTIONAL)
qtrack.download_examples("mpas_2021092400", "")

# ### Prep data

qtrack.prep_data(data_in="mpas_30km_run_2021092400.nc", data_out=prepped_data_save, cut_lev_val=70000)

# ### Curvature vorticity calculation
data_file_in = prepped_data_save  # "prepped_data_for_tracking.nc"#"analysis_and_forecast_GFS_2024062612.nc"
compute_curvvort(data_file_in, curv_file_out, njobs_in=-1)

# ### AEW Tracking step
run_tracking(input_file=curv_file_out, save_file=AEW_raw_save_file)


# ### AEW Postprocessing step
run_postprocessing(input_file=AEW_raw_save_file, real_year_used=year_in, curv_data_file=curv_file_out, save_obj_file=AEW_final_obj_file, save_nc_file=AEW_final_nc_file)

# ### OUTPUT TIME
tend = time.time()
elapsed_time = tstart - tend
print("Time to run computation: " + str(round(tend - tstart, 2)) + " seconds | " + str(round((tend - tstart) / 60, 2)) + " minutes")
print()
