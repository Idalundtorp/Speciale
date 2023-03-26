# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 13:13:39 2022

@author: Ida Olsen

Make atmosphere correction
"""

## Get ECMW data
from Read_ECMWF import ECMWF_data
import netCDF4
import matplotlib.pyplot as plt
f = netCDF4.Dataset('ECMWF_data.nc')
ECMWF = ECMWF_data(f)

ECMWF.Make_atmosphere_correction()

for i in range(len(ECMWF.Tdn[:,0])):
    plt.scatter(ECMWF.tcwv, ECMWF.Tdn[i,:], label='')
plt.show()