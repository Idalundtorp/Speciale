# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 12:04:21 2023

@author: Ida Olsen
"""

import numpy as np
import matplotlib.pyplot as plt

## Ocean temperature
T_w = (273.15-1.8)  # K
ks = 0.25 # Thermal conductivity of snow [w/(m*K)]
k = 2.10
SD = 0.13
SIT_array = [0.5, 1.0, 1.5, 2.0]
colors = ['C' + str(i) for i in range(len(SIT_array))]
Tsur = 253
plt.figure(figsize=(6,6))
for SIT,c in zip(SIT_array, colors):
    T_sii = (k * SD * T_w + ks * SIT * Tsur) / (k * SD + ks * SIT)
    
      
    from sklearn.linear_model import LinearRegression
     
    x = [0, SD]
    y = [Tsur, T_sii]
    #print(x)
    # print(y)
    model = LinearRegression().fit(np.array(x).reshape(-1,1), np.array(y))
    
     
    ice_thick = np.array([0, -SIT])
    t_ice = np.array([T_sii, T_w])
    snow_thick = [SD, 0]
    t_snow = [Tsur, T_sii]
    
    # print(t_snow)
    #print(snow_thick)
    # print(t_ice)
    
    
    plt.plot(t_snow, snow_thick,'.-', color=c, label='SIT: ' + str(SIT)+ ' [m]')
    plt.plot(t_ice, ice_thick,'.-', color=c)
# plt.ylim([-0.1,0.1])
plt.grid()
plt.legend(loc='lower left', fontsize=12)
plt.xlabel('Temperature [K]', fontsize=12)
plt.ylabel('Depth [m]', fontsize=12)
plt.title('Calculated temperature profile of Snow/ice system')


## Ocean temperature
T_w = (273.15-1.8)  # K
ks = 0.25 # Thermal conductivity of snow [w/(m*K)]
k = 2.10
SIT = 1.0
SD_array = [0.05, 0.10, 0.15, 0.20, 0.25]
colors = ['C' + str(i) for i in range(len(SIT_array))]
Tsur = 260
plt.figure(figsize=(6,6))
for SD,c in zip(SD_array, colors):
    T_sii = (k * SD * T_w + ks * SIT * Tsur) / (k * SD + ks * SIT)
    
      
    from sklearn.linear_model import LinearRegression
     
    x = [0, SD]
    y = [Tsur, T_sii]
    #print(x)
    # print(y)
    model = LinearRegression().fit(np.array(x).reshape(-1,1), np.array(y))
    
     
    ice_thick = np.array([0, -SIT])
    t_ice = np.array([T_sii, T_w])
    snow_thick = [SD, 0]
    t_snow = [Tsur, T_sii]
    
    # print(t_snow)
    #print(snow_thick)
    # print(t_ice)
    
    
    plt.plot(t_snow, snow_thick,'.-', color=c, label='SD: ' + str(SD) + ' [m]')
    plt.plot(t_ice, ice_thick,'.-', color=c)
# plt.ylim([-0.1,0.1])
plt.grid()
plt.legend(loc='lower left', fontsize=12)
plt.xlabel('Temperature [K]', fontsize=12)
plt.ylabel('Depth [m]', fontsize=12)
plt.title('Calculated temperature profile of Snow/ice system')

## Ocean temperature
T_w = (273.15-1.8)  # K
ks = 0.25 # Thermal conductivity of snow [w/(m*K)]
k = 2.10
SIT = 1.0
SD = 0.13
ks_array = [0.10, 0.20, 0.30, 0.40]
colors = ['C' + str(i) for i in range(len(ks_array))]
Tsur = 260
plt.figure(figsize=(6,6))
for ks,c in zip(ks_array, colors):
    T_sii = (k * SD * T_w + ks * SIT * Tsur) / (k * SD + ks * SIT)
    
      
    from sklearn.linear_model import LinearRegression
     
    x = [0, SD]
    y = [Tsur, T_sii]
    #print(x)
    # print(y)
    model = LinearRegression().fit(np.array(x).reshape(-1,1), np.array(y))
    
     
    ice_thick = np.array([0, -SIT])
    t_ice = np.array([T_sii, T_w])
    snow_thick = [SD, 0]
    t_snow = [Tsur, T_sii]
    
    # print(t_snow)
    #print(snow_thick)
    # print(t_ice)
    
    
    plt.plot(t_snow, snow_thick,'.-', color=c, label='ks: ' + str(ks) + ' [w/Km]')
    plt.plot(t_ice, ice_thick,'.-', color=c)
# plt.ylim([-0.1,0.1])
plt.grid()
plt.legend(loc='lower left', fontsize=12)
plt.xlabel('Temperature [K]', fontsize=12)
plt.ylabel('Depth [m]', fontsize=12)
plt.title('Calculated temperature profile of Snow/ice system')

    
