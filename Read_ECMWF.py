# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 09:28:02 2022

@author: Ida Olsen
"""

import netCDF4
import matplotlib.pyplot as plt
import numpy as np 



class ECMWF_data():
     
    # init method or constructor
    def __init__(self, t2m=False):
        ## path to data
        Base_path = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/"
        if t2m==False:
            f = netCDF4.Dataset(Base_path + 'ECMWF_data2.nc')
        else:
            f = netCDF4.Dataset(Base_path + 'ECMWF_t2m.nc')
            
        self.lat   = np.array(f.variables['latitude'])
        self.lon   = np.array(f.variables['longitude'])
        self.time  = np.array(f.variables['time'])        # "hours since 1900-01-01 00:00:00.0"
        try:
            self.sst   = np.mean(f.variables['sst'], axis=0)  # sea surface temperature
            self.t2m   = np.mean(f.variables['t2m'], axis=0)  # temperature 2m above surface
            self.tcc   = np.mean(f.variables['tcc'], axis=0)  # total cloud cover
            self.tcw   = np.mean(f.variables['tcw'], axis=0)  # total cloud liquid water
            self.tcwv  = np.mean(f.variables['tcwv'], axis=0) # total columnar water vapour
            self.skt   = np.mean(f.variables['skt'], axis=2)  # skin temperature
            self.u10   = np.mean(f.variables['u10'], axis=0)  # 10 metre U wind component
            self.v10   = np.mean(f.variables['v10'], axis=0)  # 10 metre V wind component
            self.KD_index = []
            
            self.latgrid, self.longrid = np.meshgrid(self.lat, self.lon)
        except:
            pass
        
        self.t2m_inv = f.variables['t2m']
    def plot(self, var):
        # grid_lat = np.zeros((self.lat.size, self.lon.size))
        # grid_lon = np.zeros((self.lat.size, self.lon.size))
        # for i in range(len(self.lon)):
        #     grid_lat[:, i] = self.lat
        # for i in range(len(self.lat)):
        #     grid_lon[i, :] = self.lon
        
        plt.imshow(var[0,:,:])
        plt.title('ECMW variable')
        plt.colorbar()
    
    def Make_atmosphere_correction(self):
        x = self.tcwv.flatten()
        self.Tdn = np.zeros((5, self.tcwv.size))
        ## linear regression for Tup=Tdn
        self.Tdn[0,:] = 0.0297*x+4.0639
        self.Tdn[1,:] = 0.0865*x+4.3262
        self.Tdn[2,:] = 0.7729*x+5.6202
        self.Tdn[3,:] = 2.2411*x+7.2606
        self.Tdn[4,:] = 0.8701*x+18.208

        ## linear regressions for transmissivity
        tau69  = np.mean(-1e-5*x+0.983)
        tau107 = np.mean(-0.0002*x+0.9814)
        tau187 = np.mean(-0.0028*x+0.9755)
        tau238 = np.mean(-0.0082*x+0.9683)
        tau365 = np.mean(-0.0028*x+0.9211)
        self.tau    = np.array([tau69,tau107,tau187,tau238,tau365])
        Tsp    = 2.7 #K - background radiation

# data = ECMWF_data(f)

# data.plot(data.tcwv)

