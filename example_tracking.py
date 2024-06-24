#!/usr/bin/env python
# coding: utf-8

# # Example: Using QTrack to tracking AEWs in reanalysis data
# 
# This example will demonstrate the basic functions of AEW tracking from the Qtrack module. In order to track AEWs, we must proceed through a few steps:
# 
# 1. Download or location data for tracking. Data requirments: 700hPa wind (u and v component), 6 hourly timesteps, 1 degree by 1 degree resolution. 
# 2. `prep_data` function to prep the input data for tracking. 
# 3. `nondivwind` (optional) if global data selected, user can run an optional divergent/non-divergent wind decomposition.
# 4. `curvvort` Compute the curvature vorticity from the wind fields, and then take radial averages at each gridpoint to smooth it. 
# 5. `tracking` Run the AEW tracking code. 
# 6. `postproc` Run postprocessing on AEW data, including the creation of a netCDF4 file. 


import qtrack
from qtrack.nondivwind import compute_nondiv_wind
from qtrack.curvvort import compute_curvvort
from qtrack.tracking import run_tracking, run_postprocessing
import time
tstart = time.time()


# ### Download Example ERA5 data from 2010
# The following helper script will obtain example data to test the tracker on. Available datasets include:
# - "era5_2010" ERA5 wind data from the 2010 AEW season. 

qtrack.download_examples("era5_2010_10day", "")

# ### Prep data (not completed yet)


# ### Non-divergent wind calculation (not working)
# #### **WARNING: ONLY RUN THIS STEP ON FULLY GLOBAL DATA**
# nondiv_data_file_in = "era5_700_wind_global_2010.nc"
# nondiv_data_file_out = "era5_nondiv_700_global_2010.nc"
# compute_nondiv_wind(nondiv_data_file_in, nondiv_data_file_out)


# ### Curvature vorticity calculation
data_file_in = "era5_700_wind_global_2010_10day.nc"
curv_file_out = "curv_vort_era5_test.nc"
compute_curvvort(data_file_in, curv_file_out, njobs_in = 6)


# ### AEW Tracking step
AEW_raw_save_file = 'AEW_tracks_raw.nc'
run_tracking(input_file = curv_file_out, save_file = AEW_raw_save_file)


# ### AEW Postprocessing step
AEW_final_nc_file = 'AEW_tracks_post_processed.nc'
AEW_final_obj_file = 'AEW_tracks_post_processed.pkl'
year_in = 2010
run_postprocessing(input_file = AEW_raw_save_file, real_year_used = year_in, curv_data_file = curv_file_out, save_obj_file = AEW_final_obj_file, save_nc_file = AEW_final_nc_file)

# ### OUTPUT TIME
tend = time.time()
elapsed_time = tstart - tend
print('Time to run computation: '+ str(round(tend - tstart, 2))+ ' seconds | '+str(round((tend - tstart)/60,2))+ ' minutes')
print()



