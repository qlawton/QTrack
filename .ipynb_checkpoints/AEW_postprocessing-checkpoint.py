# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 14:08:46 2020

@author: Quin7
"""

import numpy as np
import datetime
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import time as tm
from netCDF4 import Dataset, num2date, date2num
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.signal import savgol_filter
import sys
smooth_len = 7

out_attr = 'preliminary_tracks'
#year_used = 2006
year_used = sys.argv[1]
TC_distance = 500 #km
lat_cut = 25
rad_used = 600
merge_distance = 500 #Km
connect_distance = 700 #km #Distance potentially "broken" waves can be connected if their end points imply westward propagation
connect_distance_back = 200 #km #Same, but for waves that are near stationary or move in the wrong direction
connect_step = 2 #km
days_remove = 2
duplicate_removal = True
save_hovmoller = True
save_data = True
save_data_nc = True
run_animation = True
animation_smooth = True

### Animation settings ### ----------------------------------------------------
raw_ani = True #True will plot "raw" centers, false will plot smoothed. (Default: True)
wide_ani = True #Include up to gulf with true, otherwise stick over Atlantic with false (Default: True)
superwide_setting = False
wind_overlay = True #IF ANIMATION CREATED, decide if wind vectors overlain or not (Default: True)
wnd_skip = 2 #Aesthetically, controls how often wind vectors are plotted 
### ---------------------------------------------------------------------------


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    #print(idx)
    return idx, array[idx]

def haversine(lon1, lat1, lon2, lat2):
    from math import radians, cos, sin, asin, sqrt
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km

class season:
    def __init__(self, year, AEW_group):
        self.year = year
        self.AEW_group = AEW_group
    def number_of_waves(self):
        return len(AEW_group)
    def get_wave(self, wavenumber):
        return self.AEW_group[wavenumber-1]
    def waves_with_TC(self):
        waves_TC = []
        for i in range(len(self.AEW_group)):
            TC_ans = self.AEW_group[i].connected_TC
            if TC_ans == True:
                waves_TC.append((i+1))
        return waves_TC
            
def static_background(zoom = False, wider = False, superwide = False):
    #PROJECTION SETTINGS
    dataproj = ccrs.PlateCarree()
    
    #Define the extent of our static atlantic plot
    if zoom == True:
        lon1 = -70
        lon2 = 30
        lat1 = 5
        lat2 = 30
    if wider == True:
        lon1 = -110
        lon2 = 40
        lat1 = 0
        lat2 = 25
    elif superwide == True: #Modify to region of interest
        lon1 = -110
        lon2 = 40
        lat1 = 5
        lat2 = 35
    
    lon_delta = (lon2-lon1)/5
    lat_delta = (lat2-lat1)/5
    
    #Actually plot this
    fig=plt.figure(figsize=(15, 5))
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

def AEW_duplicate(lon_data, num_cutoff, close_cutoff, perc):
    lon_in = lon_data.copy()
    rm_list = []
    merge_list = []
    for row in range(np.shape(lon_in)[0]):
        #print(row)
        if row in rm_list:
            continue #This means this row has already been flagged for removal so don't bother
        else:
            sample_data = lon_in[row,:] #Sample data for row currently in
            for new_row in range(row+1,np.shape(lon_in)[0]): #Look at the data for upcoming rows
                new_data = lon_in[new_row,:]
                equal_num = np.abs(new_data-sample_data)<=close_cutoff
                equal_num = len(equal_num[equal_num == True])
                #print(old_equal_num, equal_num)
                if equal_num !=0:
                    first_i = np.where(np.abs(new_data-sample_data)<=close_cutoff)[0][0]
                    first_i2 = np.argwhere(new_data-sample_data<=close_cutoff)[0][0]
                    #print(first_i, first_i2)
                    beg_new = new_data[:(first_i+1)]
                    beg_old = sample_data[:(first_i+1)]
                    beg_new = beg_new[~np.isnan(beg_new)]
                    beg_old = beg_old[~np.isnan(beg_old)]
                    prior_num_new = len(beg_new[~np.isnan(beg_new)])
                    prior_num_old = len(beg_old[~np.isnan(beg_old)])
                    
                    # NEW EDITED SECTION.... TEST OUT
                    if prior_num_new != 0 and prior_num_old != 0:
                        if prior_num_new<=prior_num_old and (np.abs(beg_new[-1] - beg_new[0])<=5): #If the "new" data is greater or equal to 20, continue
                            #print(beg_new, beg_old)
                            rm_list.append(new_row)
                        elif prior_num_old<=prior_num_new and (np.abs(beg_old[-1] - beg_old[0])<=5): #If the old data also fits these requirements along with new data...
                            #print(beg_new, beg_old)
                            rm_list.append(row)

    return rm_list, merge_list
        
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
        #self.strength = strength
        
    def get_data(self):
        return self.time, self.lon, self.lat

datafile = 'TRACKING/prelim_aew_tracks_out.nc'
save_hov = 'TRACKING/final_hovmoller.png'

ncfile = Dataset(datafile, 'r')

# load in data
time_data = ncfile.variables['time'][:]
time_units = ncfile.variables['time'].units
lon = ncfile.variables['longitude'][:]
AEW_lon = ncfile.variables['AEW_lon'][:]
AEW_lat = ncfile.variables['AEW_lat'][:]
AEW_lon_smooth = ncfile.variables['AEW_lon_smooth'][:]
AEW_lat_smooth = ncfile.variables['AEW_lat_smooth'][:]
curv_array = ncfile.variables['curv_data_mean'][:]

reg_list = [num2date(i, time_units, only_use_cftime_datetimes = False) for i in time_data]

## WAVE MERGING

final_merge_list = np.array([])
for existing in range(np.shape(AEW_lon)[0]): #First iterate over all the waves
        AEW_lon_slc = AEW_lon[existing,:]
        AEW_lat_slc = AEW_lat[existing,:]

        for new_existing in range(existing+1,np.shape(AEW_lon)[0]):
            existing_pair = False
            new_AEW_lon_slc = AEW_lon[new_existing, :]
            new_AEW_lat_slc = AEW_lat[new_existing, :]

            for row in range(np.shape(final_merge_list)[0]):
                test = final_merge_list[row,:]
                if existing in test and new_existing in test:
                    existing_pair = True
                    break

            if existing_pair == True:
                #print('did it')
                continue
            for slc_num in range(len(new_AEW_lon_slc)):

                lon1 = AEW_lon_slc[slc_num]
                lat1 = AEW_lat_slc[slc_num]
                lon2 = new_AEW_lon_slc[slc_num]
                lat2 = new_AEW_lat_slc[slc_num]

                if np.isnan(lon1) or np.isnan(lon2) or np.isnan(lat1) or np.isnan(lat1): #If a NaN time for any wave, skip
                    continue

                two_wave_distance = haversine(lon1, lat1, lon2, lat2) #Check distance
                new_wave_data = new_AEW_lon_slc[~np.isnan(new_AEW_lon_slc)] #Pull out non-nan data
                old_wave_data = AEW_lon_slc[~np.isnan(AEW_lon_slc)] #Pull out non-nan data

                #Find lengths of each
                new_wave_last = new_wave_data[-1]
                old_wave_last = old_wave_data[-1]

                if ~np.isnan(two_wave_distance) and two_wave_distance < merge_distance: #If wave is found to be merging

                    if new_wave_last >= old_wave_last:
                        temp_merge = [existing, new_existing]
                        AEW_lon[new_existing, slc_num:] = AEW_lon[existing, slc_num:]
                        AEW_lat[new_existing, slc_num:] = AEW_lat[existing, slc_num:]
                    elif old_wave_last>= new_wave_last:
                        temp_merge = [new_existing, existing]
                        AEW_lon[existing, slc_num:] = AEW_lon[new_existing, slc_num:]
                        AEW_lat[existing, slc_num:] = AEW_lat[new_existing, slc_num:]           

                    #Finally, concatinate 
                    if not final_merge_list.any():
                        final_merge_list = np.array(temp_merge).reshape(1,2)
                    else:
                        final_merge_list = np.vstack((final_merge_list, temp_merge))
                    break

# WAVE RECONNECT
final_reconnect_list = np.array([])
for existing in range(np.shape(AEW_lon)[0]): #First iterate over all the waves
        AEW_lon_slc = AEW_lon[existing,:]
        AEW_lat_slc = AEW_lat[existing,:]

        #Pull out last position on wave track
        first_pos = np.argwhere(~np.isnan(AEW_lon_slc))[0]
        last_pos = np.argwhere(~np.isnan(AEW_lon_slc))[-1]

        lon_last = AEW_lon_slc[last_pos]
        lat_last = AEW_lat_slc[last_pos]

        #Pull out first position on wave track
        lon_first = AEW_lon_slc[first_pos]
        lat_first = AEW_lat_slc[first_pos]

        for new_existing in range(existing+1,np.shape(AEW_lon)[0]):
            new_AEW_lon_slc = AEW_lon[new_existing,:]
            new_AEW_lat_slc = AEW_lat[new_existing,:]
            new_first_pos = np.argwhere(~np.isnan(new_AEW_lon_slc))[0]
            new_last_pos = np.argwhere(~np.isnan(new_AEW_lon_slc))[-1]
            for tm_step in range(4):
                #First test to see if their last position matches up with our first position

                if first_pos <= 0+tm_step or new_last_pos!=(first_pos-1-tm_step): #If out of range, or if the timestep in question is not a wave's first or last position
                    pass
                elif lon_first>=-17 and tm_step>connect_step: #If over land, we only want to check 18 hours 
                    pass
                else:
                    new_lon_last = new_AEW_lon_slc[first_pos-1-tm_step]
                    new_lat_last = new_AEW_lat_slc[first_pos-1-tm_step]
                    test_distance1 = haversine(new_lon_last, new_lat_last, lon_first, lat_first)
                    if test_distance1<=connect_distance and lon_first<=new_lon_last: #IF the "first" position being connected is further west... bigger range
                        #print('BOOM!')
                        temp_arr = [existing, new_existing]
                        #print(existing, new_existing)
                        if not final_reconnect_list.any():
                            final_reconnect_list = np.array(temp_arr).reshape(1,2)
                        else:
                            final_reconnect_list = np.vstack((final_reconnect_list, temp_arr))
                    #WHOLE NEW BLOCK, ADDED LOGISTICS 
                    elif test_distance1<=connect_distance_back and lon_first>new_lon_last: #Same but for "backwards" wave
                        #print('BOOM!')
                        temp_arr = [existing, new_existing]
                        #print(existing, new_existing)
                        if not final_reconnect_list.any():
                            final_reconnect_list = np.array(temp_arr).reshape(1,2)
                        else:
                            final_reconnect_list = np.vstack((final_reconnect_list, temp_arr))

                if last_pos >= (np.shape(AEW_lon)[1]-1-tm_step) or new_first_pos!=(last_pos+1+tm_step): #If out of range, or if the timestep in question is not a wave's first or last position
                    pass
                elif lon_last>=-17 and tm_step>connect_step: #If over land, we only want to check 18 hours 
                    pass
                else:
                    new_lon_first = new_AEW_lon_slc[last_pos+1+tm_step]
                    new_lat_first = new_AEW_lat_slc[last_pos+1+tm_step]
                    test_distance2 = haversine(new_lon_first, new_lat_first, lon_last, lat_last)
                    #print(test_distance2)
                    if test_distance2<=connect_distance and new_lon_first<=lon_last:
                        #print(test_distance2)
                        #print('BOOM!')
                        temp_arr = [existing, new_existing]
                        #print(existing, new_existing)
                        #Finally, concatinate 
                        if not final_reconnect_list.any():
                            final_reconnect_list = np.array(temp_arr).reshape(1,2)
                        else:
                            final_reconnect_list = np.vstack((final_reconnect_list, temp_arr))
                    elif test_distance2<=connect_distance_back and new_lon_first>lon_last:
                        #print(test_distance2)
                        #print('BOOM!')
                        temp_arr = [existing, new_existing]
                        #print(existing, new_existing)
                        #Finally, concatinate 
                        if not final_reconnect_list.any():
                            final_reconnect_list = np.array(temp_arr).reshape(1,2)
                        else:
                            final_reconnect_list = np.vstack((final_reconnect_list, temp_arr))


#print(first_pos, last_pos)
#print(AEW_lon_slc[170])
#print(np.shape(AEW_lon_slc))
#print(final_reconnect_list)
connect_del_list = []
for row in range(np.shape(final_reconnect_list)[0]): #For each pair of waves to be reconnected
    wave1 = final_reconnect_list[row,0]
    wave2 = final_reconnect_list[row,1]

    #Pull out the actual waves
    wave1_lon = AEW_lon[wave1, :]
    wave1_lat = AEW_lat[wave1, :]
    wave2_lon = AEW_lon[wave2, :]
    wave2_lat = AEW_lat[wave2, :]

    for element in range(len(wave1_lon)):
        if np.isnan(wave1_lon[element]): #If an element in wave1 is a NaN, replace with that of wave2
            wave1_lon[element] = wave2_lon[element]
            wave1_lat[element] = wave2_lat[element]
    AEW_lon[wave1, :] = wave1_lon
    AEW_lat[wave1, :] = wave1_lat
    if wave2 not in connect_del_list:
        connect_del_list.append(wave2)

   #print(merge_list)

def cleanup_AEW(AEW_lon_in, AEW_lat_in, temporal_res, days_remove):
    def nan_helper(y):
        """Helper to handle indices and logical indices of NaNs. CREDIT TO: user "eat" on StackOverflow

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]

    del_list = []
    '''Oftentimes, a wave is "lost" to the tracker for 3-9 hours due to not meeting a threshold requirement, but then
    re-acheiving this requirement later on. This leaves annoying "nan" values in the missing places. Thus, this script
    will fill these holes with linear interpolated values.'''
    cutoff_len = int(np.round(days_remove*24/temporal_res))
    for row in range(np.shape(AEW_lon_in)[0]):
        AEW_lon_slc = AEW_lon_in[row,:]
        AEW_lat_slc = AEW_lat_in[row,:]
        # First, find the indices of the first 'non-nan' value and the last real value
        real_val = np.argwhere(~np.isnan(AEW_lon_slc))
        first = real_val[0][0]
        last = real_val[-1][0]

        #Linear interpolation of "intermediate" points
        nans, x = nan_helper(AEW_lon_slc)
        AEW_lon_slc[nans]= np.interp(x(nans), x(~nans), AEW_lon_slc[~nans])

        nans, x = nan_helper(AEW_lat_slc)
        AEW_lat_slc[nans]= np.interp(x(nans), x(~nans), AEW_lat_slc[~nans])

        AEW_lon_in[row,:] = AEW_lon_slc
        AEW_lat_in[row,:] = AEW_lat_slc

        #But this resulted in interpolations before and after the data -- thus, need to use "first, last" to restore to Nans
        AEW_lon_in[row, 0:first] = np.NaN
        AEW_lon_in[row, (last+1):] = np.NaN
        AEW_lat_in[row, 0:first] = np.NaN
        AEW_lat_in[row, (last+1):] = np.NaN
    return AEW_lon_in, AEW_lat_in
new_AEW_lon, new_AEW_lat = cleanup_AEW(AEW_lon, AEW_lat, 6, days_remove)
new_AEW_lon = np.delete(new_AEW_lon, connect_del_list, axis = 0)
new_AEW_lat = np.delete(new_AEW_lat, connect_del_list, axis = 0)
AEW_lat = new_AEW_lat
AEW_lon = new_AEW_lon

if duplicate_removal == True: 
    rm_dup_list, merge_list = AEW_duplicate(AEW_lon,1, 1, 0.7)

    AEW_lon= np.delete(AEW_lon, rm_dup_list, axis = 0)
    AEW_lat = np.delete(AEW_lat, rm_dup_list, axis = 0)

    #Adjust merge list values to be accurate
    for i in range(len(merge_list)):
        reduce_by = len(np.array(rm_dup_list)<merge_list[i])
        merge_list[i] = merge_list[i] - reduce_by


##TC Wave Pairing
linked_TC_name = []
linked_TC_wave = []
linked_TC_wave_time = []

#     season_pull = hurdat_atl.get_season(int(year_used))
#     TC_id_list = season_pull.summary()['id']
#     #print(TC_id_list)
#     # TC_id_list = []
#     for TC_num in range(len(TC_id_list)):


#         current_storm = season_pull.get_storm(TC_id_list[TC_num])
#         stm_type = current_storm.type
#         current_name = current_storm.name
#         current_date = current_storm.date[(stm_type!='DB') & (stm_type!='LO') & (stm_type!='WV')]
#         current_lat = current_storm.lat[(stm_type!='DB') & (stm_type!='LO') & (stm_type!='WV')]
#         current_lon = current_storm.lon[(stm_type!='DB') & (stm_type!='LO') & (stm_type!='WV')]

#         if (current_date[0] < datetime.datetime(int(year_used), 7, 1)) or np.min(current_lat)>=lat_cut:
#             continue

#         genesis_i = np.argmin(np.abs(np.array(reg_list) - current_date[0]))
#         TC_gen_lon = current_lon[0]
#         TC_gen_lat = current_lat[0]
#         TC_gen_time = current_date[0]
#         #print(current_name, TC_gen_lon, TC_gen_lat, genesis_i)

#         for wave_num in range(np.shape(AEW_lon)[0]):
#             wave_gen_lon = AEW_lon[wave_num, genesis_i]
#             wave_gen_lat = AEW_lat[wave_num, genesis_i]
#             if ~np.isnan(wave_gen_lon):
#                 gen_wave_distance = haversine(TC_gen_lon, TC_gen_lat, wave_gen_lon, wave_gen_lat)
#                 if gen_wave_distance <= TC_distance:
#                     linked_TC_wave_time.append(TC_gen_time)
#                     linked_TC_name.append(current_name)
#                     linked_TC_wave.append(wave_num)
#                     continue
#                 else:
#                     pass
#     TC_wave_frame = np.array([linked_TC_name, linked_TC_wave])   

## Find strength
curvdata = 'CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc'
curv_file = Dataset(curvdata, 'r')

time = curv_file.variables['time'][:]
time_units = curv_file.variables['time'].units
nclat = curv_file.variables['latitude'][:]
nclon = curv_file.variables['longitude'][:]

#Find, slice out only the area we are interested in (to reduce file size 
#and prevent memory overuse/dumps!)
# lon_st, n = find_nearest(nclon, -65)
# lon_end, n = find_nearest(nclon, 60)
# lat_end, n = find_nearest(nclat, 50)
# lat_st, n = find_nearest(nclat, -20)
lat = nclat#[lat_st:lat_end]
lon = nclon#[lon_st:lon_end]

curv_vort = curv_file.variables['curv_vort']#[:,lat_st:lat_end,lon_st:lon_end]

AEW_strength = np.ones(np.shape(AEW_lon))*np.NaN
for tm_i in range(np.shape(AEW_lon)[1]):
    curv_vort_slice = curv_vort[tm_i,:,:]
    for storm in range(np.shape(AEW_lon)[0]):
        lon_out = AEW_lon[storm, tm_i]
        lat_out = AEW_lat[storm, tm_i]
        if ~np.isnan(lon_out) and ~np.isnan(lat_out):
            lon_i,n = find_nearest(lon, lon_out) 
            lat_i,n = find_nearest(lat, lat_out)
            curv_value = curv_vort_slice[lat_i, lon_i]

            AEW_strength[storm, tm_i] = curv_value

#Finally, get objects
AEW_final_list = []
reg_list_edit = np.array(reg_list)

#print(TC_linked_list, TC_linked_num)

#del season
for slc_num in range(np.shape(AEW_lon)[0]):
    year_in = int(year_used)
    connected_TC = False #Hard wired, need to fix
    if np.nanmax(AEW_lon[slc_num,:])>-17:
        over_africa = True
    else:
        over_africa = False

    # if slc_num in TC_linked_num: #If the wave number is included in the wave list for linked TCs
    #     connected_TC = True
    #     name_i = np.where(TC_linked_num == slc_num)
    #     TC_connect_name = TC_linked_list[name_i][0]
    #     TC_genesis_time = linked_TC_wave_time[name_i]
    # else:
    
    ### NOTE: FOR NOW, REMOVING THE TC CONNECTION FUNCTION TO MAKE THIS TRACKER MORE UNIVERSAL. 
    # IF YOU WANT THIS FEATURE, PLEASE REACH OUT TO QUINTON AT QUINTON.LAWTON@RSMAS.MIAMI.EDU
    connected_TC = False
    TC_connect_name = 'N/A'
    TC_genesis_time = 'N/A'
    lon_in = AEW_lon[slc_num,:][~np.isnan(AEW_lon[slc_num, :])]
    lat_in = AEW_lat[slc_num,:][~np.isnan(AEW_lon[slc_num, :])]
    smooth_lon_in = savgol_filter(lon_in, smooth_len, 2, mode = 'nearest')
    smooth_lat_in = savgol_filter(lat_in, smooth_len, 2, mode = 'nearest')
    time_in = reg_list_edit[:][~np.isnan(AEW_lon[slc_num, :])]
    strength_in = AEW_strength[slc_num,:][~np.isnan(AEW_lon[slc_num, :])]

    AEW_object = AEW(year_in,slc_num+1,time_in, lon_in, lat_in, smooth_lon_in, smooth_lat_in, strength_in, over_africa, connected_TC, TC_connect_name, TC_genesis_time)  
    AEW_final_list.append(AEW_object)

season_object = season(int(year_used), AEW_final_list)
#season_object.waves_with_TC()

AEW_lon_filter = np.ones(np.shape(AEW_lon))*np.NaN
AEW_lat_filter = np.ones(np.shape(AEW_lat))*np.NaN
for row in range(np.shape(AEW_lon)[0]):
    AEW_lon_i = np.argwhere(~np.isnan(AEW_lon[row,:]))
    AEW_lon_st = AEW_lon_i[0][0]
    AEW_lon_end = AEW_lon_i[-1][0]+1
    AEW_lat_i = np.argwhere(~np.isnan(AEW_lat[row,:]))
    AEW_lat_st = AEW_lat_i[0][0]
    AEW_lat_end = AEW_lat_i[-1][0]+1
    AEW_lon_filter_pull = savgol_filter(AEW_lon[row,AEW_lon_st:AEW_lon_end], smooth_len, 2, mode = 'nearest')
    AEW_lat_filter_pull = savgol_filter(AEW_lat[row, AEW_lat_st:AEW_lat_end], smooth_len, 2, mode = 'nearest')   
    AEW_lon_filter[row,AEW_lon_st:AEW_lon_end] = AEW_lon_filter_pull
    AEW_lat_filter[row,AEW_lon_st:AEW_lon_end] = AEW_lat_filter_pull

if save_data == True:
    import pickle
    data_save_out = 'TRACKING/AEW_tracks_post_processed.pkl'
    pickle.dump(season_object, open(data_save_out, 'wb'))
    print('Saved')

if save_data_nc == True: ## Also want to save the raw output as NC file for other scripts to use
    outfile = 'TRACKING/AEW_tracks_post_processed.nc'

    nlat = np.size(lat)
    nlon = np.size(lon)
    nsys = np.shape(AEW_lon)[0]
    system_in = np.arange(nsys)+1

    # open a netCDF file to write
    ncout = Dataset(outfile, 'w', format='NETCDF4')

    # define axis size
    ncout.createDimension('time', None)  # unlimited
    ncout.createDimension('system', nsys)
    ncout.createDimension('latitude', nlat)
    ncout.createDimension('longitude', nlon)


    # create latitude axis
    latitude = ncout.createVariable('latitude', np.dtype('float64').char, ('latitude'))
    latitude.standard_name = 'latitude'
    latitude.long_name = 'latitude'
    latitude.units = 'degrees_north'
    latitude.axis = 'Y'

    # create longitude axis
    longitude = ncout.createVariable('longitude', np.dtype('float64').char, ('longitude'))
    longitude.standard_name = 'longitude'
    longitude.long_name = 'longitude'
    longitude.units = 'degrees_east'
    longitude.axis = 'X'

    # create time axis
    time_data = ncout.createVariable('time', np.dtype('float64').char, ('time'))
    time_data.long_name = 'time'
    time_data.units = time_units
    time_data.calendar = 'gregorian'
    time_data.axis = 'T'

    # create system
    system_data = ncout.createVariable('system', np.dtype('float64').char, ('system'))
    system_data.standard_name = 'system'
    system_data.long_name = 'system'
    system_data.units = 'number'
    system_data.axis = 'Y'


    # create variable array
    AEW_longitude = ncout.createVariable('AEW_lon', np.dtype('float64').char, ('system', 'time'))
    AEW_longitude.long_name = 'Longitude of AEW Tracks'
    AEW_longitude.units = 'degrees'

    AEW_latitude = ncout.createVariable('AEW_lat', np.dtype('float64').char, ('system', 'time'))
    AEW_latitude.long_name = 'Latitude of AEW Tracks'
    AEW_latitude.units = 'degrees'

    AEW_longitude_smooth = ncout.createVariable('AEW_lon_smooth', np.dtype('float64').char, ('system', 'time'))
    AEW_longitude_smooth.long_name = 'Smoothed Longitude of AEW Tracks'
    AEW_longitude_smooth.units = 'degrees'

    AEW_latitude_smooth = ncout.createVariable('AEW_lat_smooth', np.dtype('float64').char, ('system', 'time'))
    AEW_latitude_smooth.long_name = 'Smoothed Latitude of AEW Tracks'
    AEW_latitude_smooth.units = 'degrees'

    curv_data_mean = ncout.createVariable('curv_data_mean', np.dtype('float64').char, ('time','longitude'))
    curv_data_mean.long_name = 'Averaged Non-Divergent Curvature Vorticity (5-20N)'
    curv_data_mean.units = 's**-1' 

    # copy axis from original dataset
    time_data[:] = time[:]
    system_data[:] = system_in[:]
    longitude[:] = lon[:]
    latitude[:] = lat[:]
    AEW_longitude[:,:] = AEW_lon[:,:]
    AEW_latitude[:,:] = AEW_lat[:,:]
    AEW_longitude_smooth[:,:] = AEW_lon_filter[:,:]
    AEW_latitude_smooth[:,:] = AEW_lat_filter[:,:]
    curv_data_mean[:,:] = curv_array[:,:]
    ncout.close()

# -------------------------------------------------------------------------
# END DATA SAVE

#     if run_animation == True:
#         wind_dir = '/glade/scratch/qlawton/mpas_runs/'+system+'/'+year+'/wind_for_tracking_full.nc'
#         wind_file = Dataset(wind_dir, 'r')

#         wnd_lon = wind_file.variables['longitude'][:]
#         wnd_lat = wind_file.variables['latitude'][:]
#         wlon_st, n = find_nearest(wnd_lon, -60)
#         wlon_end, n = find_nearest(wnd_lon, 60)
#         wlat_end, n = find_nearest(wnd_lat, 50)
#         wlat_st, n = find_nearest(wnd_lat, -20)
#         wlat = wnd_lat[wlat_st:wlat_end]
#         wlon = wnd_lon[wlon_st:wlon_end]

#         uwnd = wind_file.variables['u'][:,wlat_st:wlat_end, wlon_st:wlon_end]
#         vwnd = wind_file.variables['v'][:,wlat_st:wlat_end, wlon_st:wlon_end]


#         import imageio, os
#         save_dir = '/glade/work/qlawton/DATA/AEW/PROCESSED/Animation/'
#         #att = info_name
#         dir_ani_frame = 'frames_'+system+'_'+year+'_R'+str(rad_used)+'edit_curvfix/'
#         GIF_dir = save_dir
#         file_list = []
#         time_list = []
#         lon_list = []
#         #mean_array = np.zeros((np.shape(time)[0],np.shape(lon)[0]) )

#         try:
#             os.mkdir(save_dir+dir_ani_frame)
#         except:
#             print('Directory could not be created -- may already exist')




#         for slc_num in range(len(time)):
#             curv_vort_data = curv_vort[slc_num, :, :]
#             uwnd_slice = uwnd[slc_num,:,:]
#             vwnd_slice = vwnd[slc_num,:,:]
#             if animation_smooth == True:
#                 AEW_lon_slc = AEW_lon_filter[:,slc_num]
#                 AEW_lat_slc = AEW_lat_filter[:,slc_num]
#             else:
#                 AEW_lon_slc = AEW_lon[:,slc_num]
#                 AEW_lat_slc = AEW_lat[:,slc_num]                
#             curv_cont = [-12,-11, -10,-9, -8,-7, -6, -5, -4, -3, -2,-1, 0,1, 2,3, 4,5, 6,7, 8,9, 10,11, 12]
#             fig, ax = static_background(zoom = True, wider = wide_ani, superwide = superwide_setting)
#             valid_time_dt = (num2date(time[slc_num], time_units, only_use_cftime_datetimes = False))
#             valid_time = datetime.datetime.strftime(valid_time_dt, '%Y%m%d%H')
#             ax.set_title('ERA-5:'+str(rad_used)+'km Non-Divergent Curvature Vorticity (10e-6 s-1) with AEW Tracks', {"fontsize": 16}, loc='left')
#             ax.set_title('VALID: {valid_time}', {"fontsize": 16}, loc='right')
#             layer1 = ax.contourf(lon, lat, curv_vort_data*10**6, curv_cont, cmap = matplotlib.cm.get_cmap('RdGy_r'))
#             aews = ax.scatter(AEW_lon_slc, AEW_lat_slc, marker = '*', s = 150, color = 'b')
#             cbar = fig.colorbar(layer1, orientation='horizontal', pad=0.1, shrink=1, aspect=30,extendrect=True, ticks = curv_cont)

#             if wind_overlay == True:
#                 X, Y = np.meshgrid(wlon, wlat)
#                 skip = (slice(None,None,wnd_skip), slice(None,None,wnd_skip))
#                 ax.quiver(X[skip], Y[skip], uwnd_slice[skip], vwnd_slice[skip], scale = 600, color = 'k')

#             temp_file = save_dir+dir_ani_frame+'smoothed_animation_banding_raw_'+str(rad_used)+'_'+str(slc_num)+'.png'
#             file_list.append(temp_file)
#             fig.savefig(temp_file)
#             plt.close()
#         with imageio.get_writer(save_dir+'smoothed_AEW_'+year+'.gif', mode='I', duration = 0.2) as writer:
#             for filename in file_list:
#                 image = imageio.imread(filename)
#                 writer.append_data(image)


if save_hovmoller == True:
    curv_cont_hov = [-12,-11, -10,-9, -8,-7, -6, -5, -4, -3, -2,-1, 0,1, 2,3, 4,5, 6,7, 8,9, 10,11, 12]

    fig = plt.figure(figsize = (10,8))
    ax = fig.add_subplot(111)
    bg = ax.contourf(lon, reg_list, curv_array*10**6, curv_cont_hov, cmap = matplotlib.cm.get_cmap('RdGy_r'),
                    extend = 'both')

    for wave_num in range(np.shape(AEW_lon)[0]):
        AEW_lon_plot = np.ma.masked_array(AEW_lon[wave_num,:], AEW_lat[wave_num,:]>25)
        if np.nanmax(AEW_lon[wave_num,:])<-17:
            ind_col = 'dimgrey'
            z_d = 4
        else:
            ind_col = 'k'
            z_d = 5
        plt.plot(AEW_lon_plot, reg_list, color = ind_col, linewidth = 3, zorder = z_d)
#         for TC_num in range(len(TC_id_list)):
#             #Storm type array

#             #current_lat = current_storm.lat
#             #current_lon = current_storm.lon
#             #current_name = current_storm.name

#             current_storm = season_pull.get_storm(TC_id_list[TC_num])
#             stm_type = current_storm.type
#             current_name = current_storm.name
#             current_date = current_storm.date[(stm_type!='DB') & (stm_type!='LO') & (stm_type!='WV')]
#             current_lat = current_storm.lat[(stm_type!='DB') & (stm_type!='LO') & (stm_type!='WV')]
#             current_lon = current_storm.lon[(stm_type!='DB') & (stm_type!='LO') & (stm_type!='WV')]

#             invest_date = current_storm.date[(stm_type=='DB') | (stm_type=='LO') | (stm_type=='WV')]
#             invest_lat = current_storm.lat[(stm_type=='DB') | (stm_type=='LO') | (stm_type=='WV')]
#             invest_lon = current_storm.lon[(stm_type=='DB') | (stm_type=='LO') | (stm_type=='WV')]

#             if (current_date[0] > datetime.datetime(int(year_used), 7, 1)) and np.min(current_lat)<=lat_cut:
#                 plt.text(current_lon[0]-15, current_date[0], current_name, color = 'b', fontsize = 10, weight = 'bold')
#             try:
#                 ax.scatter(current_lon[current_lat<=lat_cut], current_date[current_lat<=lat_cut], s=30, marker = '+', color='b', zorder = 7)
#                 ax.scatter(invest_lon[invest_lat<=lat_cut], invest_date[invest_lat<=lat_cut], s=30, marker = '+', color='g', zorder = 7)
#             except:
#                 print('fail')
#                 pass
    #ax.set_ylim([datetime.date(int(year), 7, 1), datetime.date(int(year), 9, 30)])
    plt.gca().invert_yaxis()
    ax.set_title("AEW Tracks with Curvature Vorticity", fontsize = 20)
    ax.set_xlabel('Longitude', fontsize = 20)
    ax.set_xlim(-100, 40)
    cbar = plt.colorbar(bg)
    cbar.set_label('Curv. Vort. (Avg. 5-20N, 1e-6)', fontsize = 12)
    plt.savefig(save_hov)

