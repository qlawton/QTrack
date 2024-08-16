# General Information
This repo contains code to run the AEW Tracker described by Lawton et al. (2022), which can be found here: https://doi.org/10.1175/MWR-D-21-0321.1.

Existing AEW tracks for ERA5 data from 1979 - 2022 are already produced and available here: https://osf.io/j4hpq/

## How the tracker works
In principle the AEW tracker uses an input file containing 700hPa zonal (u) and meridional (v) wind and computes a curvature vorticity (CV) term that is smoothed using a 600km average around each gridpoint. Local maxima in this CV are then tracked in two different ways. Over the African continent, meridional averages of CV are used to find longitudinal local maxima, which are used as a first guess for a centroid location. Set thresholds of CV value define when a new AEW is initialized and whether an existing one continues being tracked. Over the Atlantic, the first guess for an existing AEW's location is estimated by the AEW's propagation speed, and it is later found again using a centroid. A post processing script cleans up some of the tracks.

More detailed information on the tracking, not provided in the Lawton et al. (2022) paper or here, can be found in an online technical guide here: https://osf.io/6hqy5

## Input files
You need u and v winds on a 1x1 degree grid with 6-hourly temporal outputs. This is the optimal resolution tested in Lawton et al. (2022), but it can be adjusted manually in the code if necessary to change. Furthermore, it is highly recommended that at least 10 days of data are included. The tracking needs a bit of spinup, so if you are running this on model output, I recommend you append at least 1 week of analysis/reanalysis data prior to the first model timestep.

A sample model data/reanalysis combination script is provided under `TOOLS/wind_combine.py`

## Output AEW Track files
Two different kinds of output files are produced, a netCDF file and a python pickle file (with AEWs saved as objects). It is recommended that you use netCDF if possible. However, if you need to use the pickle files for whatever reason, you will need to have the `AEW_module.py` script available locally when you open them.

Details on the AEW objects and AEW_module.py are available here: https://osf.io/jnv5u

## Python environment
There are specific module requirements that are  listed in the "REQUIRED_MODULES" file.

# Details on script contents

## Main shell scripts (`run_tracker.sh`)
This shell script will run all of the python scripts for AEW tracking in sequence. Note that it is required that you specify an initial input file, for convenience, this should be placed in the same directory as these scripts. There are also several assorted subdirectoires that are included for automatic organization. Feel free to look at the scripts and change the output directories/names to something that is more convenient for you. Several scripts also have the ability to change the output location or input file name: please see below for more details on that.

## Non-Divergent Wind Calculation ( `non_divergent_wind.py`)
**This step is optional and can be skipped.** It should be emphasized that the following non-divergent wind script REQUIRES a global domain (as it uses spherical harmonics). Otherwise, you should skip this step. This step likely has minor ramifications for the final output. The necessary windspharm module is also tricky to install, and so if you have trouble I recommend skilling this step.

`python non_divergent_wind.py [INPUT FILE NAME (required)] [OUTPUT NAME (optional)]`

Default output path is "CURV_VORT/HELMHOLTZ/wind_700_helmholtz.nc"

## Computing Curvature Vorticity and Averaging (`SAVE_CURVE_VORT_SHELL_NON_DIV.py`)
This step computes the CV and takes the gridpoint averages. It is slow. However, there is an option within the script to use multiprocessing, which is highly recommended.

`python SAVE_CURV_VORT_SHELL_NON_DIV.py [INPUT FILE NAME (optional)] [OUTPUT FILE NAME (optional)]`

Default input path is "CURV_VORT/HELMHOLTZ/wind_700_helmholtz.nc", default output path is "CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc"

## Actual AEW Tracking Step (`AEW_TRACKING_CODE.py`)

This step runs the AEW tracker on the computed CV output from the previous steps. It is fairly quick to run.

`python AEW_TRACKING_CODE.py [INPUT FILE NAME (optional)] [WIND INPUT FILE NAME (optional)]`

Default input path is "CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc", default wind input path is "wind_for_tracking.nc"

## Post-Processing of AEW Tracks (`AEW_postprocessing.py`)
This step computes the netCDF4 files and saves the data there. It also tries to clean up the tracked AEW data. Importantly, there is a setting that eliminates AEWs that are not at least 2 days long. This can be adjusted within the script if necessary. There is also a feature to identify developing AEWs using HURDAT data (ONLY use this for reanalysis inputs), but this is buggy at the moment. This is why it is required that a year be input into this call, as this allows the code to select the correct year for comparision.

`python AEW_postprocessing.py [YEAR OF DATA (required)]`
