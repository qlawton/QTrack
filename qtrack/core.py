# -*- coding: utf-8 -*-
"""
Core Functions for the QTrack algorithm and associated functions. 

Written by Quinton Lawton. 
"""

# IMPORT STATEMENTS
import numpy as np
import datetime
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import time as tm
import imageio
from numpy import dtype
import os
import dill as pickle
import xarray as xr
from joblib.externals.loky import set_loky_pickler
from joblib import Parallel, delayed
set_loky_pickler("dill")

class season:
    def __init__(self, year, AEW_group):
        self.year = year
        self.AEW_group = AEW_group
    def number_of_waves(self):
        return len(self.AEW_group)
    def get_wave(self, wavenumber):
        return self.AEW_group[wavenumber-1]
    def waves_with_TC(self):
        waves_TC = []
        for i in range(len(self.AEW_group)):
            TC_ans = self.AEW_group[i].connected_TC
            if TC_ans == True:
                waves_TC.append((i+1))
        return waves_TC
            
        
class AEW:
    def __init__(self, year, number, time, lon, lat, smooth_lon, smooth_lat, strength, over_africa, connected_TC, connected_TC_name = 'N/A', genesis_time = 'N/A'):
        self.year = year
        self.number = number
        self.time = time
        self.lon = lon
        self.lat = lat
        self.strength = strength
        self.over_africa = over_africa
        self.connected_TC = connected_TC
        self.connected_TC_name = connected_TC_name
        self.TC_genesis_time = genesis_time
        self.smooth_lon = smooth_lon
        self.smooth_lat = smooth_lat
        
    def get_data(self):
        return self.time, self.lon, self.lat
    
class AEW_CCKW:
    def __init__(self, year, number, time, lon, lat, smooth_lon, smooth_lat, 
                 strength, over_africa, connected_TC, TB_lat_range, VP_lat_range,
                 CCKW_TB, CCKW_VP, sig_TB_peak, sig_VP_peak, sig_TB_valley, sig_VP_valley,
                 connected_TC_name = 'N/A', genesis_time = 'N/A'):
        self.year = year
        self.number = number
        self.time = time
        self.lon = lon
        self.lat = lat
        self.strength = strength
        self.over_africa = over_africa
        self.connected_TC = connected_TC
        self.connected_TC_name = connected_TC_name
        self.TC_genesis_time = genesis_time
        self.smooth_lon = smooth_lon
        self.smooth_lat = smooth_lat
        
        #Also include the new stuff we wnat
        self.CCKW_TB = CCKW_TB
        self.CCKW_VP = CCKW_VP
        self.sig_TB_peak = sig_TB_peak
        self.sig_VP_peak = sig_VP_peak
        self.sig_TB_valley = sig_TB_valley
        self.sig_VP_valley = sig_VP_valley
        self.TB_lat_range = TB_lat_range
        self.VP_lat_range = VP_lat_range
def download_examples(input_str, output_dir=""):
    """
    This script will download example files. Available files are:
    - "era5_2010"

    """
    import gdown
        
    if input_str == "era5_2010":
        out_str = "era5_700_wind_global_2010.nc"
        print('Downloading ERA5 Test case to: '+output_dir+out_str)
        g_id = "https://drive.google.com/uc?id=1aLJpb1UAFTdmZRX56bHVn25EzuTcL8d_" 
        gdown.download(g_id, output_dir+out_str)  
        print('Downloaded Successfully.')
    else:
        raise Exception("Invalid example data name provided. Please check that you are inputting a valid file name.")
        
def COMPUTE_CURV_VORT_NON_DIV_UPDATE(data_in, data_out, res, radius, njobs, nondiv = True, SAVE_IMAGE = False, SAVE_OUTPUT = True, gif_dir_in = ''):
    import numpy as np
    import time
    dir_ani_frame = 'frames_R'+str(radius)
    
    # ### IMPORTANT DIRECTORIES AND CUSTOMIZATIONS 
    gif_dir = gif_dir_in
    data_dir = ''
    data_in = data_in
    data_out = data_out
    
    if SAVE_IMAGE == True:
        print('Setting Up | Output to: '+data_out)
        try:
            os.mkdir(gif_dir+dir_ani_frame)
        except:
            print('Directory could not be created -- may already exist')
        
    save_dir = gif_dir+dir_ani_frame+'/'
    
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        #print(idx)
        return idx, array[idx]
    
    def get_dist_meters(lon, lat):
        earth_circ = 6371*2*np.pi*1000 #earth's circumference in meters
        lat_met = earth_circ/360 #get the number of meters in a degree latitude (ignoring "bulge")
        lat_dist = np.gradient(lat, axis=0)*lat_met
        lon_dist = np.gradient(lon, axis=1)*np.cos(np.deg2rad(lat))*lat_met
        return lon_dist, lat_dist 
    
    def curv_vort(u, v, dx, dy):
        V_2 = (u**2+v**2)
        curv_vort_raw = (1/V_2)*(u*u*np.gradient(v, axis=1)/dx - v*v*np.gradient(u, axis=0)/dy  
                                 - v*u*np.gradient(u, axis=1)/dx + u*v*np.gradient(v, axis=0)/dy)
    
        return curv_vort_raw
    
    def static_background(zoom = False):
        #PROJECTION SETTINGS
        dataproj = ccrs.PlateCarree()
        
        #Define the extent of our static atlantic plot
        if zoom == True:
            lon1 = -70
            lon2 = 30
            lat1 = 5
            lat2 = 30
        else: #Modify to region of interest
            lon1 = -60
            lon2 = 0
            lat1 = 0
            lat2 = 40
        
        lon_delta = (lon2-lon1)/5
        lat_delta = (lat2-lat1)/5
        
        #Actually plot this
        fig=plt.figure(figsize=(15, 15))
        ax=plt.subplot(111, projection=dataproj)
        ax.set_extent([lon1, lon2, lat1, lat2],ccrs.PlateCarree())
        ax.coastlines('50m', linewidth=1.5)
        ax.add_feature(cfeature.STATES, linewidth=1.0)
        ax.add_feature(cfeature.BORDERS, linewidth=1.0)
        gl = ax.gridlines(color='gray',alpha=0.5,draw_labels=True) 
        gl.xlabels_top, gl.ylabels_right = False, False
        gl.xlabel_style, gl.ylabel_style = {'fontsize': 16}, {'fontsize': 16}
        gl.xlocator = mticker.FixedLocator([lon1,lon1+lon_delta, lon1+2*lon_delta, lon1+3*lon_delta, lon1+4*lon_delta, lon1+5*lon_delta])
        gl.ylocator = mticker.FixedLocator([lat1,lat1+lat_delta, lat1+2*lat_delta, lat1+3*lat_delta, lat1+4*lat_delta, lat1+5*lat_delta])
        gl.xformatter = LongitudeFormatter(zero_direction_label=True)
        gl.yformatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(gl.xformatter)
        ax.yaxis.set_major_formatter(gl.yformatter)
        ax.stock_img()
    
        return fig, ax
        
    def GetBG(lon, lat, data_file, res, radius):
        #from joblib import Parallel, delayed

        start= tm.time()
    
        def get_dist_meters(lon, lat):
            earth_circ = 6371*2*np.pi*1000 #earth's circumference in meters
            lat_met = earth_circ/360 #get the number of meters in a degree latitude (ignoring "bulge")
            lat_dist = np.gradient(lat, axis=0)*lat_met
            lon_dist = np.gradient(lon, axis=1)*np.cos(np.deg2rad(lat))*lat_met
            return lon_dist, lat_dist 
        
        #Get the distance array for later use
        dx, dy = get_dist_meters(lon, lat)
        
        #BUFFER -- let's give this a radius*lat buffer for longitude as well
        earth_circ = 6371*2*np.pi*1000 #earth's circumference in meters
        lat_met = earth_circ/360 #get the number of meters in a degree latitude (ignoring "bulge")
        lat_km_lat = lat_met/1000/1
        lat_km_lon = lat_met/1000/2
            
        buffer = int(np.ceil(radius/lat_km_lat/res))
        buffer_lon = int(np.ceil(radius/lat_km_lon/res))
        
        i_bounds = np.arange((buffer+1), np.shape(data_file)[0]-buffer)
        j_bounds = np.arange((buffer_lon+1), np.shape(data_file)[1]-buffer_lon)
        
    
        #RADIUS IS IN KILOMETERS
        data_shape= np.shape(data_file)
        #DEFINE ACTUAL BOUNDS ------------------------------------------------------
            #This is done because averaging must allow for a complete radius to
            #accurately consider values. Thus, a buffer is put on the bounds considered
            #equal to the radius of the area considered 
        avg_grid= np.zeros(data_shape)
    
        def rad_mask(i, j, dx, dy, radius):
                start = tm.time()
    
                earth_circ = 6371*2*np.pi*1000 #earth's circumference in meters
                lat_met = earth_circ/360 #get the number of meters in a degree latitude (ignoring "bulge")
                lat_km_lat = lat_met/1000/1 #Updated to fix high latitude errors...
                lat_km_lon = lat_met/1000/2 #Updated to work up to 60N...
                res = 1
                buffer = int(np.ceil(radius/lat_km_lat/res))
                buffer_j = int(np.ceil(radius/lat_km_lon/res))
                boolean_array = np.zeros(np.shape(dx), dtype=bool)
                #print(buffer)
                #i_array = np.zeros(np.shape(dy))
                #j_array = np.zeros(np.shape(dx))
                dy = -dy
    
                # Before doing this, we recognize something -- there is a max radius in the x and y direction that will
                # be computed. To save on computational time, we need to slice out the extraneous points that we already
                # know from this NOT to be in the circle. How? We know the absolute highest distance between gridspaces
                # will be the maximum value of latitude on the elliptical earth... ~ 111 km. So slice out roughly a 112km "box" 
                # plus one gridbox for buffer. 
                i_st = (i-(buffer+1))
                i_end = (i+(buffer+1))
                j_st = (j-(buffer_j+1))
                j_end = (j+(buffer_j+1))
                #print(i_st, i_end)
                #print(j_st, j_end)
    
                new_i = i-i_st
                new_j = j-j_st
    
                dy_slc = dy[i_st:i_end,j_st:j_end]
                dx_slc = dx[i_st:i_end,j_st:j_end]
    
                #print(np.shape(dy_slc), np.shape(dy))
    
    
                i_array_sub = np.zeros(np.shape(dy_slc))
                j_array_sub = np.zeros(np.shape(dx_slc))
    
                i_array_sub[new_i,:] = 0
                i_array_sub[(new_i+1):,:] = np.add.accumulate(dy_slc[(new_i+1):,:])
                i_array_sub[(new_i-1)::-1,:] = np.add.accumulate(dy_slc[(new_i-1)::-1,:])
    
                j_array_sub[:,new_j] = 0
                j_array_sub[:,(new_j+1):] = np.add.accumulate(dx_slc[:,(new_j+1):], axis=1)
                j_array_sub[:,(new_j-1)::-1] = np.add.accumulate(dx_slc[:,(new_j-1)::-1], axis=1)
    
                radial_array = (np.sqrt(np.square(i_array_sub)+np.square(j_array_sub))/1000) < radius # in km
                boolean_array[i_st:i_end,j_st:j_end] = radial_array
                end = tm.time()
                #print(end - start)
                return boolean_array
        
        for y in i_bounds: #This will iterate over ALL given points and calculate a background average
            for x in j_bounds:
                outMask = rad_mask(y,x, dx, dy, radius)
                avg_grid[y,x]= np.mean(data_file[outMask == True])
            #Parallel(n_jobs=-1)(delayed(run_code)(y, x) for x in j_bounds)
                    
        end =tm.time()
        #print('Time elapsed:'+str(end-start)+' s')
        
        return avg_grid
    
    ### ----------------------- ACTUAL COMPUTATIONAL BLOCK --------------------------- ###
    # First, import the data
    print('Starting Computation of Radial Averaged CV...')
    nc_file = xr.open_dataset(data_in)
    nc_file_alt = Dataset(data_in, 'r')
    
    time = nc_file_alt.variables['time'][:]
    time_units = nc_file_alt.variables['time'].units
    nclat = nc_file['latitude'].values
    nclon = nc_file['longitude'].values
    
    #Find, slice out only the area we are interested in (to reduce file size and prevent memory overuse/dumps!)
    lat = nclat
    lon = nclon
    
    if nondiv == False:
        u_wnd = nc_file['u'].values
        v_wnd = nc_file['v'].values
    else:
        u_wnd = nc_file['upsi'].values
        v_wnd = nc_file['vpsi'].values
    
    LON, LAT = np.meshgrid(lon, lat)
    dx, dy = get_dist_meters(LON, LAT)
    
    out_array = np.zeros((np.shape(time)[0], np.shape(LAT)[0], np.shape(LAT)[1]))
        
    def run_loop(slc_num, radius):
        #file_list = []
        print('Timestep number: '+str(slc_num))
        out_name = data_dir + 'curv_temp_data_'+str(slc_num)+'.npy'
        curv_vort_data = curv_vort(u_wnd[slc_num,:,:], v_wnd[slc_num,:,:], dx, dy)
        set_radius = radius
        curv_array = GetBG(LON, LAT, curv_vort_data,res, set_radius)
                
        if SAVE_IMAGE == True:
            curv_cont = [-12,-11, -10,-9, -8,-7, -6, -5, -4, -3, -2,-1, 0,1, 2,3, 4,5, 6,7, 8,9, 10,11, 12]
            fig, ax = static_background(zoom = True)
            valid_time_dt = (num2date(time[slc_num], time_units, only_use_cftime_datetimes = False))
            valid_time = datetime.datetime.strftime(valid_time_dt, '%Y%m%d%H')
            ax.set_title('DATA_OUT:'+str(set_radius)+'km Radially-Averaged Curvature Vorticity (10e-6 s-1)', {"fontsize": 16}, loc='left')
            ax.set_title(f'VALID: {valid_time}', {"fontsize": 16}, loc='right')
            layer1 = ax.contourf(lon, lat, curv_array*10**6, curv_cont, cmap = matplotlib.cm.get_cmap('RdGy_r'))
            cbar = fig.colorbar(layer1, orientation='horizontal', pad=0.05, shrink=1, aspect=30,extendrect=True, ticks = curv_cont)
    
            temp_file = save_dir+'curv_vort_helmholz'+'_R'+str(set_radius)+'_'+str(slc_num)+'.png'
            #file_list.append(temp_file)
            fig.savefig(temp_file, bbox_inches='tight')
            plt.close()
        
        if SAVE_OUTPUT == True:
            np.save(out_name, curv_array)
            
            
    start = tm.time()
    Parallel(n_jobs=njobs)(delayed(run_loop)(i, radius) for i in np.arange(len(time)))
    end = tm.time()
    print('Time to run computation: '+ str(end - start))
    
    if SAVE_IMAGE == True:
        file_list = []
        for i in np.arange(len(time)):
            temp_file = save_dir+'curv_vort_helmholz_R'+str(radius)+'_'+str(i)+'.png'
            file_list.append(temp_file)
        
        with imageio.get_writer(gif_dir+'curv_vort_helmholtz_ani_R'+str(radius)+'.gif', mode='I') as writer:
            for filename in file_list:
                image = imageio.imread(filename)
                writer.append_data(image)
            
        
    if SAVE_OUTPUT == True:
        for i in np.arange(len(time)):
            in_file = data_dir + 'curv_temp_data_'+str(i)+'.npy'
            temp_file = np.load(in_file)
            out_array[i, :, :] = temp_file[:,:]
            os.remove(in_file)
        #Finally, save the np file
        nlat = np.size(lat)
        nlon = np.size(lon)
        # open a netCDF file to write
        ncout = Dataset(data_out, 'w', format='NETCDF4')
        
        # define axis size
        ncout.createDimension('time', None)  # unlimited
        ncout.createDimension('latitude', nlat)
        ncout.createDimension('longitude', nlon)
        
        # create time axis
        time_data = ncout.createVariable('time', dtype('float64').char, ('time',))
        time_data.long_name = 'time'
        time_data.units = time_units
        time_data.calendar = 'gregorian'
        time_data.axis = 'T'
        
        # create latitude axis
        latitude = ncout.createVariable('latitude', dtype('float64').char, ('latitude'))
        latitude.standard_name = 'latitude'
        latitude.long_name = 'latitude'
        latitude.units = 'degrees_north'
        latitude.axis = 'Y'
        
        # create longitude axis
        longitude = ncout.createVariable('longitude', dtype('float64').char, ('longitude'))
        longitude.standard_name = 'longitude'
        longitude.long_name = 'longitude'
        longitude.units = 'degrees_east'
        longitude.axis = 'X'
        
        # create variable array
        curvout = ncout.createVariable('curv_vort', dtype('float64').char, ('time', 'latitude', 'longitude'))
        curvout.long_name = 'Radially Averaged Curvature Vorticity'
        curvout.units = 's**-1'
        
        
        # copy axis from original dataset
        time_data[:] = time[:]
        longitude[:] = lon[:]
        latitude[:] = lat[:]
        curvout[:,:,:] = out_array[:,:,:]
        ncout.close()