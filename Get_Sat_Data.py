# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 19:08:36 2022

@author: Ida Olsen

Getting AMSR2 and ECMWF data at co-located locations
"""

import numpy as np
import scipy.spatial as ss 
import pyproj
import netCDF4
import os
import datetime as dt
import matplotlib.pyplot as plt

## home made modules
from Read_ECMWF import ECMWF_data
from read_AMSR2 import AMSR2_data
from Functions_temp import fill_data

comb_info = fill_data()

def plot_t2m(self, obj, all_sites=False):
    # hours since 1900-01-01 00:00:00.0
    date_start = dt.datetime(1900, 1, 1, 0, 0, 0)
    times = [date_start + dt.timedelta(hours = h.astype(list)) for h in self.time]
    date_to_str = [t.strftime("%m-%d-%H") for t in times]
    plt.figure(figsize=(10,6))
    
    plt.grid()
    plt.scatter(date_to_str, self.t2m_inv)
    # tsur = [obj.Tsur for t in times]
    # plt.plot(date_to_str, tsur)
    # print(obj.Tsur)
    plt.xticks(date_to_str[::6], rotation=45)
    plt.xlabel('Date (month-day-hour)', fontsize=12)
    plt.ylabel('T2m [k]', fontsize=12)
    plt.title('Temperatures at 2m height during CB campaign at site: ' + obj.site, fontsize=14)
    # plt.ylim([242, 255])
    plt.savefig(os.path.join("C:/Users/Ida Olsen/Documents/Speciale_2022/figures/CB_figures",
                             "t2m_change_site_" + obj.site + ".png"), bbox_inches='tight')
    plt.show()

def Get_AMSR_ECMWF(landfree=True, original_AMSR2=False, plots=False):   
    ## path to data
    Base_path = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/"
    
    # define projection
    # Define standad projection for the Lambert Azimuthal Equeal Area projection
    Pollat0  = 90.
    Pollon0  = 0.
    ell = 'WGS84'
    PolProj=pyproj.Proj(proj='laea', ellps=ell, datum='WGS84', lat_0=Pollat0, lon_0=Pollon0, units='m')
    
    ## Get ECMW data
    ECMWF_s = ECMWF_data()
    ## project ECMWF data 
    ECx,ECy = PolProj(ECMWF_s.longrid.flatten(),ECMWF_s.latgrid.flatten())
    # print('projection of ECMWF succesfull')
    
    # Make list with AMSR2 and ECMWF classes for each site
    AMSR2 = []
    AMSR2_orig = []
    ECMWF = []
    for s in comb_info:
        if s.site == 'AK1':
            filename = Base_path + "5000003659957/AMSR_U2_L3_SeaIce25km_B04_20220420.he5"
        elif s.site == 'AK2':
            filename = Base_path + "5000003659957/AMSR_U2_L3_SeaIce25km_B04_20220421.he5"
        elif s.site == 'AK3':
            filename = Base_path + "5000003659957/AMSR_U2_L3_SeaIce25km_B04_20220423.he5"
        elif s.site == 'AK4':
            filename = Base_path + "5000003659957/AMSR_U2_L3_SeaIce25km_B04_20220424.he5"
    
        ## Get AMSR2 data A factor 10 too high TBs?
        AMSR2_s = AMSR2_data(filename)
        AMSR2_s.lat = AMSR2_s.lat.flatten()
        AMSR2_s.lon = AMSR2_s.lon.flatten()
        
        ## Extract subset of AMSR2 data based on site data
        latmin = np.amin(s.lat)- 1.0
        latmax = np.amax(s.lat)+ 1.0
        lonmin = np.amin(s.lon)- 5.0
        lonmax = np.amax(s.lon)+ 5.0
    
        ## make co-located AMSR class
        AMSR2_col = AMSR2_data(filename)
    
        ii = np.where([AMSR2_s.lat[ii]>latmin and AMSR2_s.lat[ii]<latmax and AMSR2_s.lon[ii]>lonmin and AMSR2_s.lon[ii]<lonmax for ii in range(AMSR2_s.lat.size)])[0]
        AMSR2_col.lat = AMSR2_s.lat[ii]
        AMSR2_col.lon = AMSR2_s.lon[ii]
        AMSR2_col.DataSetH     = np.concatenate((AMSR2_s.DataSetH),axis=0)[ii]
        AMSR2_col.DataSetV     = np.concatenate((AMSR2_s.DataSetV),axis=0)[ii]
        AMSR2_col.DataSetH_asc = np.concatenate((AMSR2_s.DataSetH_asc),axis=0)[ii] 
        AMSR2_col.DataSetV_asc = np.concatenate((AMSR2_s.DataSetV_asc),axis=0)[ii] 
        AMSR2_col.DataSetH_dec = np.concatenate((AMSR2_s.DataSetH_dec),axis=0)[ii] 
        AMSR2_col.DataSetV_dec = np.concatenate((AMSR2_s.DataSetV_dec),axis=0)[ii] 
        
        ## Append AMSR2 subset before co-location
        AMSR2_orig.append(AMSR2_s)
        
        ## project AMSR, obs and ECMWF data and find nearest neighbour
        AMSR2_col.projx,AMSR2_col.projy = PolProj(AMSR2_col.lon,AMSR2_col.lat)
        #print('projection of AMSR2 succesfull')
    
        obsx,obsy = PolProj(s.lon,s.lat)
        #print('projection of observation data succesfull')
    
        treeAMSR2 = ss.cKDTree(list(zip(AMSR2_col.projx,AMSR2_col.projy)))
        treeECMWF = ss.cKDTree(list(zip(ECx,ECy)))
        
        if landfree==True: ## landfree AMSR2 value
            lat_AMSR = 68.5
            lon_AMSR = -103.5
            obsx_AMSR,obsy_AMSR = PolProj(lon_AMSR,lat_AMSR)
        else:
            obsx_AMSR,obsy_AMSR = obsx, obsy 

        KDstructAMSR2 = treeAMSR2.query(list(zip([obsx_AMSR], [obsy_AMSR])))[1]
        KDstructECMWF = treeECMWF.query(list(zip([obsx], [obsy])))[1]
        #print('KDstructs made')
        
        ## get co-located AMSR data
        AMSR2_col.lat = AMSR2_s.lat[KDstructAMSR2]
        AMSR2_col.lon = AMSR2_s.lon[KDstructAMSR2]
        AMSR2_col.DataSetH     = AMSR2_col.DataSetH[KDstructAMSR2]
        AMSR2_col.DataSetV     = AMSR2_col.DataSetV[KDstructAMSR2]
        AMSR2_col.DataSetH_asc = AMSR2_col.DataSetH_asc[KDstructAMSR2]
        AMSR2_col.DataSetV_asc = AMSR2_col.DataSetV_asc[KDstructAMSR2]
        AMSR2_col.DataSetH_dec = AMSR2_col.DataSetH_dec[KDstructAMSR2]
        AMSR2_col.DataSetV_dec = AMSR2_col.DataSetV_dec[KDstructAMSR2]
        
        ## Get co-located ECMWF data
        ECMWF_col = ECMWF_data()
        ECMWF_col.sst  = ECMWF_s.sst.flatten()[KDstructECMWF]
        ECMWF_col.skt  = ECMWF_s.skt.flatten()[KDstructECMWF]
        ECMWF_col.t2m  = ECMWF_s.t2m.flatten()[KDstructECMWF]
        ECMWF_col.tcc  = ECMWF_s.tcc.flatten()[KDstructECMWF]
        ECMWF_col.tcw  = ECMWF_s.tcw.flatten()[KDstructECMWF]
        ECMWF_col.tcwv = ECMWF_s.tcwv.flatten()[KDstructECMWF]
        ECMWF_col.u10  = ECMWF_s.u10.flatten()[KDstructECMWF]
        ECMWF_col.v10  = ECMWF_s.v10.flatten()[KDstructECMWF]
        ECMWF_col.latgrid  = ECMWF_s.latgrid.flatten()[KDstructECMWF]
        ECMWF_col.longrid  = ECMWF_s.longrid.flatten()[KDstructECMWF]
        ## temp data
        ECMWF_t2m = ECMWF_data(t2m=True)
        ECMWF_col.time = ECMWF_t2m.time
        ECMWF_col.t2m_inv = np.array(ECMWF_t2m.t2m_inv).reshape(120, 25*9)[:,KDstructECMWF].flatten()
        
        ## Append classes to list
        AMSR2.append(AMSR2_col)
        ECMWF.append(ECMWF_col)
        if plots==True:
            plot_t2m(ECMWF_col, s)
            # AMSR2 and site plots
            ## reference points (e.g. location of sites)
            points = [np.array((s.lat, s.lon)) for s in comb_info]
            AMSR2_s.plot(ii, 1, 'V', plot_point=points, obs=True)
    if original_AMSR2 == True:
        return [AMSR2_orig, ii, AMSR2, ECMWF]
    else:
        return [AMSR2, ECMWF]
#%% MAKE temperature development plots for ECMWF data
# [AMSR2, ECMWF] = Get_AMSR_ECMWF(plots=True)
# for (E, s) in zip(ECMWF, comb_info):
#     plot_t2m(E, s.site)
# #%% Define point class to use for AMSR2 plots
# class point():
#     def __init__(self, index, sitenum, ref):
        
#         self.x = AMSR2_orig[sitenum].lat[index]
#         self.y = AMSR2_orig[sitenum].lon[index]

#         A = [(self.x[i],self.y[i]) for i in range(len(self.x))]

#         distances = np.linalg.norm(A-ref, axis=1)
#         self.index = min_index = np.argmin(distances)

#         self.point = [self.x[min_index], self.y[min_index]]
#%% AMSR2 and site plots 

# points = [np.array((s.lat, s.lon)) for s in comb_info]
# AMSR2_orig[0].plot(ii, 1, 'V', plot_point=points, obs=True)

# ## define reference position of landfree point
# point_plot = point(ii, 0, np.array(( 68.5,-103.5)))
# AMSR2_orig[0].plot(ii, 0, 'V', plot_point=point_plot)
# AMSR2_orig[0].plot(ii, 0, 'V', plot_point=point(ii, 0, np.array((comb_info[0].lat,comb_info[0].lon))))