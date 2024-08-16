#Quinton Lawton, University of Miami, 2022
# ## Combine 700hPa wind files between ERA5 and MPAS, ultimate goal to create one large CV file for tracking


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pandas as pd

### Settings
case_in = 'VICTOR'
run_in = 'CPEX_RUN'
year_in = '2021'
cut_start = '2021-09-24-00' #The start time of the MPAS data
cut_end = '2021-10-03-19' #The end of the MPAS data
cut_start_early = '2021-09-23-18' #6 hours before the MPAS start -- basically, where the ERA5 data will end when merged
era5_start = '2021-09-01-00' #Where we want to start the ERA5 data from when merging
save_data = True # Save output of data

out_name = 'wind_for_tracking' #What you want to call
hr_delta = 1 #Time delta of the MPAS output (hours)
hr_plot_delta = 6 #Time delta of DESIRED output (hours, 6 hours for AEW tracking is good)

### *************** WARNING: HARDWIRED PATHS, CHANGE  **************
mpas_in = '/glade/u/home/qlawton/scratch/mpas_runs/'+case_in+'/'+run_in+'/latlon_iso_for_VP.nc' #MPAS data in (1x1, converted)
era5_in = '/glade/u/home/qlawton/DATA/AEW_Tracks/wind_season_lowres_700_'+year_in+'_B1-6hr.nc' #ERA5 data in
mpas_orig = '/glade/scratch/rberrios/cpex-aw/2021092400/diag.2021-09-24_00.00.00.nc' #Random non-converted diag file, for reference

### *************** WHAT YOU WANT TO SAVE THE FILE AS **************
out_file = '/glade/u/home/qlawton/scratch/mpas_runs/'+case_in+'/'+run_in+'/'+out_name+'.nc'
### ----------------------------------------------------------------

### LOAD IN THE DATA
mpas_xr = xr.open_dataset(mpas_in, chunks = 'auto')
era5_xr = xr.open_dataset(era5_in, chunks = 'auto').sel(time = slice(era5_start, cut_start_early))
orig_xr = xr.open_dataset(mpas_orig)
u_iso_levels_real = orig_xr['u_iso_levels']

## Start adjusting the mpas data (NOTE: These could changed based on the
mpas_xr.coords['u_iso_levels'] = u_iso_levels_real.values #Assign coordinate since this isn't also output
mpas_xr = mpas_xr.sel(u_iso_levels = 700*100) #We want the 700hPa level

## Next, we want to round the coordinates and rename the variables in the MPAS files
lat_round = np.round(mpas_xr.coords['latitude'].values)
lon_round = np.round(mpas_xr.coords['longitude'].values)
#(Rounding is necessary because sometimes convert_mpas outputs lons like 19.9999999999 instead of 20)

# Replace with the real value
mpas_xr.coords['longitude'] = lon_round
mpas_xr.coords['latitude'] = lat_round

# Rename the variables (this depends on your data files
mpas_xr = mpas_xr.rename({'umeridional_isobaric':'v', 'uzonal_isobaric':'u'})
mpas_xr = mpas_xr.rename_dims({'Time':'time'})

#For later slicing
lat_min_mpas = mpas_xr.latitude.values[0]
lat_max_mpas = mpas_xr.latitude.values[-1]

#Note that the latitudes are reversed in ERA5 hence the denotation, but double check your data
lat_min = era5_xr.latitude.values[-1]
lat_max = era5_xr.latitude.values[0]
lon_min = era5_xr.longitude.values[0]
lon_max = era5_xr.longitude.values[-1]

#Flip and cut down the era5 data
era5_xr = era5_xr.reindex(latitude=era5_xr.latitude[::-1]).sel(latitude = slice(lat_min_mpas, lat_max_mpas))

#Finally, cut down the longitude of the era5_data
mpas_xr = mpas_xr.sel(longitude = slice(lon_min, lon_max), latitude = slice(lat_min, lat_max))


##### --- UPDATING TIME FOR MPAS ----- #####
### We also want to make a list of times to represented the start and end dates contained in the files
cut_start_dt = datetime.strptime(cut_start, '%Y-%m-%d-%H')
cut_end_dt = datetime.strptime(cut_end, '%Y-%m-%d-%H')
length_dt = int((cut_end_dt - cut_start_dt).total_seconds()/3600/hr_delta)+1 #Not sure about the plus one here but it works....
timedelta_out = timedelta(hours=hr_delta)

date_list_dt = [cut_start_dt + timedelta_out*hour for hour in range(length_dt)]

date_list = [datetime.strftime(in_dt, '%Y-%m-%d-%H') for in_dt in date_list_dt]
date_list_6hr = date_list[::6]
date_list_dt_6hr = date_list_dt[::6]

### Now put it in the data and cut it down
mpas_xr.coords['time'] = date_list_dt
mpas_xr = mpas_xr.sel(time = date_list_dt_6hr)


##### --- MERGE THE TWO ARRAYS AND SAVE --- #####
if save_data == True:
    merged = xr.concat([era5_xr, mpas_xr], dim = 'time')
    merged = merged.fillna(0)
    merged.to_netcdf(out_file)
