#!/bin/bash
# This shell script will run all of the python scripts for AEW tracking in sequence. Note that it is required that you
# specify an initial input file, for convenience, this should be placed in the same directory as these scripts.
# There are also several assorted subdirectoires that are included for automatic organization. Feel free to look at the
# scripts and change the output directories/names to something that is more convenient for you. Several scripts also have
# the ability to change the output location or input file name: please see below for more details on that.

# ---- PYTHON ENVIRONMENT ----
# First, this is specific code to run on the NCAR Cheyenne system so it should be modified for your use.
# There are specific module requirements, so you may have to check the included scripts to see which you need
# to include in your envionment.

# ------ FIRST, COMPUTE THE NON-DIVERGENT WINDS. -----
# You need u and v winds on a 1x1 degree grid. It is highly recommended that at least 10 days of data are included.
# The tracking needs a bit of spinup, so if you are running this on model output, I recommend you append at least
# 1 week of analysis/reanalysis data prior to the first model timestep.

# NOTE: The following non-divergent wind script REQUIRES a global domain (uses spherical harmonics). Otherwise, you
# should skip this step. This step likely has minor ramifications for the final output.
echo "RUNNING: NON-DIVERGENT WIND COMPUTATION"
python non_divergent_wind.py EXAMPLE_DATA/wind_season_lowres_700_2010_global.nc #non_divergent_wind.py [INPUT FILE NAME (required)] [OUTPUT NAME (optional)]I
#Default output path is "CURV_VORT/HELMHOLTZ/wind_700_helmholtz.nc"
echo "RUNNING: CURVATURE VORTICITY RADIAL AVERAGING"
python SAVE_CURV_VORT_SHELL_NON_DIV.py #SAVE_CURV_VORT_SHELL_NON_DIV.py [INPUT FILE NAME (optional)] [OUTPUT FILE NAME (optional)]
#Default input path is "CURV_VORT/HELMHOLTZ/wind_700_helmholtz.nc", default output path is "CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc"
echo "RUNNING: ACTUAL TRACKER"
python AEW_TRACKING_CODE.py #AEW_TRACKING_CODE.py [INPUT FILE NAME (optional)] [WIND INPUT FILE NAME (optional)]
#Default input path is "CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc", default wind input path is "wind_for_tracking.nc"
echo "RUNNING: Post-Processing of AEW Tracks"
python AEW_postprocessing.py 2010 #AEW_postprocessing.py [YEAR OF DATA (required)]

echo "TRACKING COMPLETED"
#This is the end of the shell script.
