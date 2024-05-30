# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 09:23:42 2020

@author: Quin7
"""

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