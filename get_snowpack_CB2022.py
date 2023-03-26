# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 21:47:19 2022

@author: Ida Olsen
"""
#%%
import os
parrent_dir = "C:\\Users\\Ida Olsen\\Documents\\Speciale_2022"
import sys
sys.path.append(parrent_dir + '\\smrt')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
# Import saline snow permittivity model
# from smrt.permittivity.saline_snow import saline_snow_permittivity_geldsetzer09 as ssp
from smrt.permittivity.saline_snow import saline_snow_permittivity_scharien_with_stogryn95 as ssp
from smrt import make_ice_column, make_snowpack, sensor_list
from smrt.core.globalconstants import PSU
from smrt import make_model
#%% Make snowpit class
from Snowpit_class import Snowpit_info
from Get_Sat_Data import Get_AMSR_ECMWF
from roughness import roughness
from test_roughness import test_rough
from Make_plots import plot_freq
from Functions_temp import load_ref_data, compute_bias_std
pd.set_option('mode.chained_assignment',None)
#Get Data info
Base_path = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/"
with open(Base_path + 'Code/smrt_in-out/smp_profile_CB', 'rb') as pickle_file:
    smp_profile = pickle.load(pickle_file)
#%%
# Make new SMRT configuration with returns
m = make_model("iba", "dort")

frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]

beamorientation = 55. # [degrees]
sensor = sensor_list.passive(frequencies, beamorientation)

CB_sp = []
count2 = 0
# SD_mean_site
labels_comb = []
colors_comb = []
## Plot all together
data_comb = np.zeros((4*4,5))
dataErr_comb = np.zeros((4*4,5))

name_prev = False
Calc_temp = False
from Functions_temp import Calc_temp_profiles
comb_info = Calc_temp_profiles(plot=True)
[AMSR, ECMWF] = Get_AMSR_ECMWF()

for (s, A, E) in zip(comb_info, AMSR, ECMWF):
    
    ## plot results
    profile_count = -1
    for name, profile in smp_profile[s.site].items():
        
        profile_count +=1
        
        ## If we use measured temperatures
        if Calc_temp==False:
            temp_ice = s.Tcal_ice[profile_count] # 250 K -very wrong :) 
            # temp_ice = profile.temperature.tail(1)
            temp_snow = profile.temperature
        ## Test output with calculated temp profiles
        else: ## If we use calculated temperature profile
            #print('temp array', s.temp[profile_count])
            #print('Tsur',s.Tsur)
            temp_ice = s.Tcal_ice[profile_count]
            temp_snow = s.Tcal_snow[profile_count]
            #print('Tsnow', temp_snow)

        ## Make icepack        
        sea_ice_density = 900
        icepack = make_ice_column(ice_type='firstyear',
                thickness=[s.SIT_m], temperature=temp_ice, 
                microstructure_model='exponential',
                corr_length=0.15e-3,
                brine_inclusion_shape='spheres',
                density=sea_ice_density,
                salinity=4*PSU,
                add_water_substrate="ocean")

        # Convert to exp. corr length
        debye = 1.5
        lex_array = debye * 4 * (1 - profile.density / 917) / (profile.ssa * 917)
        
        ## Make snowpack
        # sal = np.array(profile.salinity)
        # sal[-1] += 8
        # print(temp_snow)
        # print(profile.thick)
        snowpack = make_snowpack(profile.thick, microstructure_model='exponential',
                        ice_permittivity_model=ssp,
                        density=profile.density , 
                        corr_length=lex_array,
                        temperature = temp_snow,
                        salinity=profile.salinity*PSU)

        medium = snowpack + icepack
        CB_sp.append(medium)                     
    
    #% Run simulations
    # run simulations
    res = m.run(sensor, CB_sp)
    
    # apply roughness and atmosphere correction
    #from test_roughness import test_rough
    # test_rough(s, A, E, res)
    [TB_SMRT_V,TB_SMRT_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(E.tcwv,res,np.mean(s.Tcal_ice),site=s.site,several=True)
    
    #% make plots
    [data, dataErr] = load_ref_data(s.site)
    ## mean data
    data[4,:] = TB_SMRT_V
    data[5,:] = TB_SMRT_H
    data_comb[count2:count2+4,:] = data[2:6,:]
    
    #data[0,:] = np.mean(res.TbV(), axis=1)
    #data[1,:] = np.mean(res.TbH(), axis=1)
    
    
    ## Errors
    dataErr[4,:] = np.std(res.TbV(), axis=1)
    dataErr[5,:] = np.std(res.TbH(), axis=1)
    dataErr_comb[count2:count2+4,:] = dataErr[2:6,:]
    
    colors = ['g', 'r', 'purple', 'brown', 'b', 'orange']
    labels = ['', '', 'AMSR2 V landfree', 'AMSR2 H landfree', 'TBV SMRT', 'TBH SMRT']
    plot_freq(s, data, labels, cs=colors, err=dataErr, title='TBs simulated at CB ' + s.site + ',' + s.date)
    
    count2 += 4
    
    ## print bias
    compute_bias_std(s,data)

## two colors
#colors_comb = np.append(colors_comb,[['b', 'brown'] for i in range(int(data_comb.shape[0]/2))])
## 4 colors
colors_comb = np.append(colors_comb,[[colors[2:6]] for i in range(int(data_comb.shape[0]/4))])
plot_freq(s, data_comb, cs=colors_comb, labels=labels[2:6])

d = [np.mean(data_comb[i::4,:], axis=0) for i in range(4)]
d_err = [np.mean(dataErr_comb[i::4,:], axis=0) for i in range(4)]
plot_freq(s, np.array(d), labels[2:6], cs=colors[2:6], err=np.array(d_err), title='Average simulated and measured TBs at CB')

compute_bias_std(s,d,nullindex=False)
#%% Make plots of variables
SD_comb = []
rho_comb = []
# for c in comb_info:
    # c.plot_var(c.SD_tot, ' SD ', 'm')
    # c.plot_var(c.SD_tot, ' SD ', 'm')
    # c.plot_var(c.corr_m*10**3, ' Correlation length ', 'mm')
    # c.plot_var(c.SIT, ' SIT ', 'kg/m3')
    # c.plot_var(c.sal_m*1e3, ' Salinity ', 'kg/m3')
    # c.plot_indv(c.sal, ' Salinity PSU', plot_all=True)
    # c.plot_indv(c.rho, ' density [kg/m^3]', plot_all=True)
    # c.plot_var(c.SIT, ' SIT ', 'm')
    # c.plot_var(c.SD, ' SIT ', 'm')
#     var = c.corr
#     varname = 'corr'
#     for i in range(len(c.SD)):
#         SD_cum = np.array(np.cumsum(c.SD[i]))
#         varss = np.array(var[i])
#         SD_comb.append((SD_cum))
#         rho_comb.append((varss))
    
# plt.figure()
# plt.gca().invert_yaxis()

# # plt.scatter( np.concatenate(rho_comb), np.concatenate(SD_comb))
# z = np.polyfit(np.concatenate(rho_comb), np.concatenate(SD_comb), 1)
# p = np.poly1d(z)
# xx = np.linspace(np.min(np.concatenate(rho_comb)), np.max(np.concatenate(rho_comb)), 100)
# plt.plot( np.concatenate(rho_comb)*1e3, np.concatenate(SD_comb), '.', xx*1e3, p(xx), '-')
# plt.title('profiles of ' + varname)
# plt.ylabel('SD [m]')
# plt.xlabel(varname + '[mm]')
# plt.grid()
# plt.show()

