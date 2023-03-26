# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:05:13 2022

@author: Ida Olsen

Get AMSR variation around the line in the fjord
"""
import numpy as np
from read_AMSR2 import AMSR2_data
import pyproj
import matplotlib.pyplot as plt
from matplotlib import cm, colors

pit = 'AK4'
#%% Get satelitte data

## AK4 info
AK4_lat = 68.82158
AK4_lon = -105.41325		
AK4_time = '24/04/2022'

# from read_AMSR2 import Get_AMSR2_data
path = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/5000003659957/"

filename = path + "AMSR_U2_L3_SeaIce25km_B04_20220424.he5"
obslat = AK4_lat
obslon = AK4_lon

#%% 
AMSR2f = AMSR2_data(filename)

latmin = np.amin(obslat)- 1.0 # 0.5
latmax = np.amax(obslat)+ 1.0 # 0.5
lonmin = np.amin(obslon)- 5.0 # 2.0
lonmax = np.amax(obslon)+ 5.0 # 2.0

index = []
for ii in range(AMSR2f.lat.size):
    if AMSR2f.lat.flatten()[ii]>latmin and AMSR2f.lat.flatten()[ii]<latmax and AMSR2f.lon.flatten()[ii]>lonmin and AMSR2f.lon.flatten()[ii]<lonmax: 
        index = np.append(index, ii)

#%% Find AMSR2 points along line
class point():
    def __init__(self, ref):
        
        self.x = AMSR2f.lat.flatten()[index.astype(int)]
        self.y = AMSR2f.lon.flatten()[index.astype(int)]

        A = [(self.x[i],self.y[i]) for i in range(len(self.x))]

        distances = np.linalg.norm(A-ref, axis=1)
        self.index = min_index = np.argmin(distances)

        self.point = [self.x[min_index], self.y[min_index]]
#%% Define search points
# start points (left corner) 69.05, -106.51
stepx = 0.11
stepy = 1.5/3

pointsx = [69.05 - stepx*i for i in range(11) ]
pointsy = [-106.51 + stepy*i for i in range(11) ]

colors = ['C' + str(int(num)) for num in range(11)]

## Find closest AMSR2 points
points = [point(np.array(( pointsx[i],pointsy[i]))) for i in range(len(pointsx))]

AMSR2f.plot(0, 0, 'V', plot_point=points, c=colors)

#%% 
AMSR2_landfree = AMSR2_data(filename)

AMSR2_landfree.DataSetH = np.zeros((5,11))
AMSR2_landfree.DataSetV = np.zeros((5,11))

point_index = [points[i].index for i in range(len(points))]

color = np.zeros((11))
for i in range(len(points)):
    color[i] = i

for i in range(AMSR2_landfree.DataSetH.shape[0]):
    AMSR2_landfree.DataSetH[i] = AMSR2f.DataSetH[:,:,i].flatten()[index.astype(int)][point_index]
    AMSR2_landfree.DataSetV[i] = AMSR2f.DataSetV[:,:,i].flatten()[index.astype(int)][point_index]

np.save(pit + '_TBV_Landfree', AMSR2_landfree.DataSetV)
np.save(pit + '_TBH_Landfree', AMSR2_landfree.DataSetH)


frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9])/1e9;  # [Hz]

freq_array = [f for f in frequencies for i in range(11)]

plt.figure()

colors = ['C' + str(int(num)) for num in color]
for i in range(AMSR2_landfree.DataSetV.shape[1]):
    plt.plot(frequencies, AMSR2_landfree.DataSetV[:,i], c=colors[i], label='Point: ' + str(i))
    plt.plot(frequencies, AMSR2_landfree.DataSetH[:,i], c=colors[i])

plt.ylim([210,260])
plt.legend(loc='center left',bbox_to_anchor=([1.0,0.6]))
plt.title('Fjord examination',fontsize=16)
plt.xlabel('freq GHz',fontsize=14)
plt.ylabel('TB',fontsize=14)
plt.grid()
plt.show()
