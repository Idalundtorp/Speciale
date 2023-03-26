# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:25:17 2022

@author: Ida Olsen
"""


import numpy as np
import datetime as dt
import sys
import scipy.spatial as ss 
import pyproj
import os
import datetime as dt

from polar_plots_pred import plot
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#%% Data from CB
''' 
Data is from Cambridge Bay and consist of 4 snowpits where SMP, 
salinity and roughness measurements are made. 
'''
## AK1 info
AK1_lat = 69.08897	
AK1_lon = -105.50422						
AK1_time = '20/04/2022'	

## AK2 info
AK2_lat = 68.99985
AK2_lon = -105.47468		
AK2_time = '21/04/2022'


## AK3 info
AK3_lat = 68.91067
AK3_lon = -105.44359					
AK3_time = '23/04/2022'	

## AK4 info
AK4_lat = 68.82158
AK4_lon = -105.41325		
AK4_time = '24/04/2022'	

#%% Plot changes in temperature from ECMWF at 2m level 
def plot_t2m(self, site, all_sites=False):
    # hours since 1900-01-01 00:00:00.0
    date_start = dt.datetime(1900, 1, 1, 0, 0, 0)
    times = [date_start + dt.timedelta(hours = h.astype(list)) for h in self.time]
    date_to_str = [t.strftime("%m-%d-%H") for t in times]
    t2m_time = [np.concatenate((t2m))[KDstructECMWF[0]] for t2m in self.t2m_inv]
    plt.figure(figsize=(10,6))
    
    plt.grid()
    plt.scatter(date_to_str, t2m_time)
    plt.xticks(date_to_str[::2], rotation=45)
    plt.xlabel('Date (month-day-hour)', fontsize=12)
    plt.ylabel('T2m [k]', fontsize=12)
    plt.title('Temperatures at 2m height during CB campaign at site: ' + pit, fontsize=14)
    plt.ylim([242, 255])
    plt.savefig(os.path.join("C:/Users/Ida Olsen/Documents/Speciale_2022/figures/CB_figures",
                             "t2m_change_site_" + site + ".png"), bbox_inches='tight')
    plt.show()
#%% Get satelitte data

for pit in ['AK1', 'AK2', 'AK3', 'AK4']:
    # from read_AMSR2 import Get_AMSR2_data
    path = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/5000003659957/"
    if pit == 'AK1':
        filename = path + "AMSR_U2_L3_SeaIce25km_B04_20220420.he5"
        obslat = AK1_lat
        obslon = AK1_lon
    elif pit == 'AK2':
        filename = path + "AMSR_U2_L3_SeaIce25km_B04_20220421.he5"
        obslat = AK2_lat
        obslon = AK2_lon
    elif pit == 'AK3':
        filename = path + "AMSR_U2_L3_SeaIce25km_B04_20220423.he5"
        obslat = AK3_lat
        obslon = AK3_lon
    elif pit == 'AK4':
        filename = path + "AMSR_U2_L3_SeaIce25km_B04_20220424.he5"
        obslat = AK4_lat
        obslon = AK4_lon
    
    ## Get AMSR2 data A factor 10 too high TBs?
    from read_AMSR2 import AMSR2_data
    AMSR2 = AMSR2_data(filename)
    
    
    ## Get ECMW data
    from Read_ECMWF import ECMWF_data
    import netCDF4
    ECMWF = ECMWF_data()
    
    latmin = np.amin(obslat)- 1.0 #0.5
    latmax = np.amax(obslat)+ 1.0 # 0.5
    lonmin = np.amin(obslon)- 5.0 #2.0
    lonmax = np.amax(obslon)+ 5.0 # 2.0
    
    index = []
    for ii in range(AMSR2.lat.size):
        if AMSR2.lat.flatten()[ii]>latmin and AMSR2.lat.flatten()[ii]<latmax and AMSR2.lon.flatten()[ii]>lonmin and AMSR2.lon.flatten()[ii]<lonmax: 
            index = np.append(index, ii)     
    
    # Define standad projection for the Lambert Azimuthal Equeal Area projection
    Pollat0  = 90.
    Pollon0  = 0.
    #Pollatts = 70.
    ell = 'WGS84'
    PolProj=pyproj.Proj(proj='laea', ellps=ell, datum='WGS84', lat_0=Pollat0, lon_0=Pollon0, units='m')  #
                                                                                    # k_0=Scale factor at natural origin (normally 1.0)
    
    CSx,CSy = PolProj(AMSR2.lon.flatten()[index.astype(int)],AMSR2.lat.flatten()[index.astype(int)])
    print('projection of AMSR2 succesfull')
    
    ECx,ECy = PolProj(ECMWF.longrid.flatten(),ECMWF.latgrid.flatten())
    print('projection of ECMWF succesfull')
        
    obsx,obsy = PolProj(obslon,obslat)
    print('projection of observation data succesfull')
    

    treeAMSR2 = ss.cKDTree(list(zip(CSx, CSy)))
    
    treeECMWF = ss.cKDTree(list(zip(ECx,ECy)))
    
    KDstructAMSR2 = treeAMSR2.query(list(zip([obsx], [obsy])))[1]
    
    KDstructECMWF = treeECMWF.query(list(zip([obsx], [obsy])))[1]
    
    print('KDstructs made')
    
    #%% Get co-located TBs
      
    ## Make collocated AMSR2 data class
    AMSR2_col = AMSR2_data(filename)
    
    AMSR2_col.lat = np.zeros((5))
    AMSR2_col.lon = np.zeros((5))
    
    AMSR2_col.DataSetH = np.zeros((5,1))
    AMSR2_col.DataSetV = np.zeros((5,1))
    AMSR2_col.DataSetH_asc = np.zeros((5,1))
    AMSR2_col.DataSetV_asc = np.zeros((5,1))
    AMSR2_col.DataSetH_dec = np.zeros((5,1))
    AMSR2_col.DataSetV_dec = np.zeros((5,1))
    
    AMSR2_col.lat = [np.mean(AMSR2.lat.flatten()[index.astype(int)][KDstructAMSR2]) for i in range(5)]  
    AMSR2_col.lon = [np.mean(AMSR2.lon.flatten()[index.astype(int)][KDstructAMSR2]) for i in range(5)]  
    AMSR2_col.DataSetH = [np.mean(AMSR2.DataSetH[:,:,i].flatten()[index.astype(int)][KDstructAMSR2]) for i in range(5)]  
    AMSR2_col.DataSetV = [np.mean(AMSR2.DataSetV[:,:,i].flatten()[index.astype(int)][KDstructAMSR2]) for i in range(5)]  
    AMSR2_col.DataSetH_asc = [np.mean(AMSR2.DataSetH_asc[:,:,i].flatten()[index.astype(int)][KDstructAMSR2]) for i in range(5)]  
    AMSR2_col.DataSetV_asc = [np.mean(AMSR2.DataSetV_asc[:,:,i].flatten()[index.astype(int)][KDstructAMSR2]) for i in range(5)]  
    AMSR2_col.DataSetH_dec = [np.mean(AMSR2.DataSetH_dec[:,:,i].flatten()[index.astype(int)][KDstructAMSR2]) for i in range(5)]  
    AMSR2_col.DataSetV_dec = [np.mean(AMSR2.DataSetV_dec[:,:,i].flatten()[index.astype(int)][KDstructAMSR2]) for i in range(5)]  

    parentDirectory = "C:/Users/Ida Olsen/Documents/Speciale_2022/scripts/npy_files/"
    np.save(os.path.join(parentDirectory,pit + '_TBV'), AMSR2_col.DataSetV)
    np.save(os.path.join(parentDirectory,pit + '_TBV_err'), abs(np.array(AMSR2_col.DataSetV_asc)-np.array(AMSR2_col.DataSetV_dec)))
    #np.save(pit + '_TBV_dec', AMSR2_col.DataSetV_dec)
    
    np.save(os.path.join(parentDirectory,pit + '_TBH'), AMSR2_col.DataSetH)
    np.save(os.path.join(parentDirectory,pit + '_TBH_err'), abs(np.array(AMSR2_col.DataSetH_asc)-np.array(AMSR2_col.DataSetH_dec)))
    # np.save(pit + '_TBHasc', AMSR2_col.DataSetH_asc)
    # np.save(pit + '_TBHdec', AMSR2_col.DataSetH_dec)
    
        
    ## Get co-located ECMWF data

    ECMWF_col = ECMWF_data()
    for ii in range(len(KDstructECMWF)):
        ECMWF_col.sst  = np.mean(ECMWF.sst.flatten()[KDstructECMWF[ii]])
        ECMWF_col.skt  = np.mean(ECMWF.skt.flatten()[KDstructECMWF[ii]])
        ECMWF_col.t2m  = np.mean(ECMWF.t2m.flatten()[KDstructECMWF[ii]])
        ECMWF_col.tcc  = np.mean(ECMWF.tcc.flatten()[KDstructECMWF[ii]])
        ECMWF_col.tcw  = np.mean(ECMWF.tcw.flatten()[KDstructECMWF[ii]])
        ECMWF_col.tcwv = np.mean(ECMWF.tcwv.flatten()[KDstructECMWF[ii]])
        ECMWF_col.u10  = np.mean(ECMWF.u10.flatten()[KDstructECMWF[ii]])
        ECMWF_col.v10  = np.mean(ECMWF.v10.flatten()[KDstructECMWF[ii]])
        
        ECMWF_col.latgrid  = ECMWF.latgrid.flatten()[KDstructECMWF[ii]]
        ECMWF_col.longrid  = ECMWF.longrid.flatten()[KDstructECMWF[ii]]
    
    np.save(os.path.join(parentDirectory,pit + '_ECMWF_tcwv'), np.array(ECMWF_col.tcwv))
    # np.save(os.path.join(parentDirectory,pit + '_ECMWF_skt'), np.array(ECMWF_col.skt))
    np.save(os.path.join(parentDirectory,pit + '_ECMWF_t2m'), np.array(ECMWF_col.t2m))
    
    plot_t2m(ECMWF_col, pit)
    

    #%% Find AMSR2 point with no land influence
    class point():
        def __init__(self, ref):
            
            self.x = AMSR2.lat.flatten()[index.astype(int)]
            self.y = AMSR2.lon.flatten()[index.astype(int)]
    
            A = [(self.x[i],self.y[i]) for i in range(len(self.x))]
    
            distances = np.linalg.norm(A-ref, axis=1)
            self.index = min_index = np.argmin(distances)
    
            self.point = [self.x[min_index], self.y[min_index]]
    #%%  
    # pointsx = [AK1_lat, AK2_lat, AK3_lat, AK4_lat]
    # pointsy = [AK1_lon, AK2_lon, AK3_lon, AK4_lon]

    
    # ## Find closest AMSR2 points
    # points = [np.array(( pointsx[i],pointsy[i])) for i in range(len(pointsx))]
    
    
    
    # AMSR2.plot(index.astype(int)[0], 0, 'V', plot_point=points, obs=True)
    
    
    point_plot = point(np.array(( 68.5,-103.5)))
    
    # AMSR2.plot(index.astype(int), 0, 'V', plot_point=point_plot)
    
    # AMSR2.plot(index.astype(int), 0, 'V', plot_point=point(np.array((obslat,obslon))))
    
    # AMSR2.plot_asc(index.astype(int), 0, 'V')
    
    # AMSR2.plot_dec(index.astype(int), 0, 'V')
    
    AMSR2_landfree = AMSR2_data(filename)
    
    AMSR2_landfree.DataSetH = np.zeros((5,1))
    AMSR2_landfree.DataSetV = np.zeros((5,1))
    AMSR2_landfree.DataSetH_asc = np.zeros((5,1))
    AMSR2_landfree.DataSetV_asc = np.zeros((5,1))
    AMSR2_landfree.DataSetH_dec = np.zeros((5,1))
    AMSR2_landfree.DataSetV_dec = np.zeros((5,1))
    
    for i in range(AMSR2_landfree.DataSetH.size):
        AMSR2_landfree.DataSetH[i] = AMSR2.DataSetH[:,:,i].flatten()[index.astype(int)][point_plot.index]
        AMSR2_landfree.DataSetV[i] = AMSR2.DataSetV[:,:,i].flatten()[index.astype(int)][point_plot.index]
        AMSR2_landfree.DataSetH_asc[i] = AMSR2.DataSetH_asc[:,:,i].flatten()[index.astype(int)][point_plot.index]
        AMSR2_landfree.DataSetV_asc[i] = AMSR2.DataSetV_asc[:,:,i].flatten()[index.astype(int)][point_plot.index]
        AMSR2_landfree.DataSetH_dec[i] = AMSR2.DataSetH_dec[:,:,i].flatten()[index.astype(int)][point_plot.index]
        AMSR2_landfree.DataSetV_dec[i] = AMSR2.DataSetV_dec[:,:,i].flatten()[index.astype(int)][point_plot.index]
    
    np.save(os.path.join(parentDirectory,pit + '_TBV_Landfree'), AMSR2_landfree.DataSetV)
    np.save(os.path.join(parentDirectory,pit + '_TBV_Landfree_err'), abs(np.array(AMSR2_landfree.DataSetV_asc)-np.array(AMSR2_landfree.DataSetV_dec)))
    
    # np.save(pit + '_TBV_Landfree_asc', AMSR2_landfree.DataSetV_asc)
    # np.save(pit + '_TBV_Landfree_dec', AMSR2_landfree.DataSetV_dec)
    
    np.save(os.path.join(parentDirectory,pit + '_TBH_Landfree'), AMSR2_landfree.DataSetH)
    np.save(os.path.join(parentDirectory,pit + '_TBH_Landfree_err'), abs(np.array(AMSR2_landfree.DataSetH_asc)-np.array(AMSR2_landfree.DataSetH_dec)))
    
    # np.save(pit + '_TBH_Landfree_asc', AMSR2_landfree.DataSetH_asc)
    # np.save(pit + '_TBH_Landfree_dec', AMSR2_landfree.DataSetH_dec)
    
    
    # frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
    
    # plt.figure()
    # plot4,=plt.plot(frequencies, AMSR2_landfree.DataSetV,'.-',label='AMSR2 V')
    # plot4,=plt.plot(frequencies, AMSR2_landfree.DataSetV_asc,'.-',label='AMSR2 V asc')
    # plot4,=plt.plot(frequencies, AMSR2_landfree.DataSetV_dec,'.-',label='AMSR2 V dec')
    # plot3,=plt.plot(frequencies, AMSR2_landfree.DataSetH,'.-',label='AMSR2 H')
    # plot3,=plt.plot(frequencies, AMSR2_landfree.DataSetH_asc,'.-',label='AMSR2 H asc')
    # plot3,=plt.plot(frequencies, AMSR2_landfree.DataSetH_dec,'.-',label='AMSR2 H dec')
    
    # plot4,=plt.plot(frequencies, AMSR2_col.DataSetV,'.-',label=' col. AMSR2 V')
    # plot4,=plt.plot(frequencies, AMSR2_col.DataSetH,'.-',label=' col AMSR2 H')
    # plot4,=plt.plot(frequencies, AMSR2_col.DataSetV_asc,'.-',label=' col. AMSR2 V asc')
    # plot4,=plt.plot(frequencies, AMSR2_col.DataSetH_asc,'.-',label=' col AMSR2 H asc')
    # plot4,=plt.plot(frequencies, AMSR2_col.DataSetV_dec,'.-',label=' col. AMSR2 V dec')
    # plot4,=plt.plot(frequencies, AMSR2_col.DataSetH_dec,'.-',label=' col AMSR2 H dec')
    
    # plt.legend(bbox_to_anchor=[1,1.05])
    # plt.show()
