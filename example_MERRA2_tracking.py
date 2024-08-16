#!/usr/bin/env python

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


import time

import numpy as np

from qtrack.curvvort import compute_curvvort
from qtrack.tracking import run_postprocessing, run_tracking

tstart = time.time()

year_list = np.arange(1981, 2022, 1)
gen_dir = '/data/qlawton/MERRA2/'
for year_i in range(len(year_list)):
    year_in = year_list[year_i]
    print('Running on Year: '+str(year_in))
    # ### NAME FILES
    file_in = gen_dir+'NONDIV/nondiv_wind_MERRA2_700hPa_wind_year_'+str(year_in)+'.nc'
    prepped_data_save = file_in #gen_dir+'adjusted_nondiv_wind_MERRA2_700hPa_wind_year_'+str(year_in)+'.nc'
    curv_file_out = gen_dir+"CURV_VORT/curv_vort_MERRA2_700hPa_year_"+str(year_in)+'.nc'
    AEW_raw_save_file = gen_dir+'TRACKING/AEW_tracks_raw_year_'+str(year_in)+'.nc'
    AEW_final_nc_file = gen_dir+'TRACKING/AEW_tracks_post_processed_year_'+str(year_in)+'.nc'
    AEW_final_obj_file = gen_dir+'TRACKING/obj_AEW_tracks_post_processed_year_'+str(year_in)+'.pkl'
    hov_save_out = gen_dir+'TRACKING/HOV/hovmoller_diagram_year_'+str(year_in)+'.png'

    # ### Prep data

    #qtrack.prep_data(data_in = file_in,
    #                data_out = prepped_data_save, cut_lev_val = 70000)

    # ### Curvature vorticity calculation
    data_file_in = prepped_data_save #"prepped_data_for_tracking.nc"#"analysis_and_forecast_GFS_2024062612.nc"
    compute_curvvort(data_file_in, curv_file_out, njobs_in = 40, nondiv_wind=True)

    # ### AEW Tracking step
    run_tracking(input_file = curv_file_out, save_file = AEW_raw_save_file)


    # ### AEW Postprocessing step
    run_postprocessing(input_file = AEW_raw_save_file, real_year_used = year_in, curv_data_file = curv_file_out, save_obj_file = AEW_final_obj_file, save_nc_file = AEW_final_nc_file, TC_merge_dist = 500, hov_save_file = hov_save_out, TC_pairing = True)

    # ### OUTPUT TIME
    tend = time.time()
    elapsed_time = tstart - tend
    print('Time to run computation: '+ str(round(tend - tstart, 2))+ ' seconds | '+str(round((tend - tstart)/60,2))+ ' minutes')
    print()
