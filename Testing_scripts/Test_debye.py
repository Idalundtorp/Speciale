# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 17:01:57 2022

@author: Ida Olsen
"""
import numpy as np
parrent_dir = "C:\\Users\\Ida Olsen\\Documents\\Speciale_2022"
import sys
sys.path.append(parrent_dir + '\\smrt')
import pandas as pd
import pickle

from smrt.permittivity.saline_snow import saline_snow_permittivity_geldsetzer09 as ssp
from smrt import make_ice_column, make_interface, make_snowpack, sensor_list
from smrt.core.globalconstants import PSU
from smrt import make_model
from smrt.interface.iem_fung92_brogioni10 import IEM_Fung92_Briogoni10


# Make new SMRT configuration with returns
m = make_model("iba", "dort")

frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]

beamorientation = 55. # [degrees]
sensor = sensor_list.passive(frequencies, beamorientation)

pd.set_option('mode.chained_assignment',None)


with open('smrt_in-out/smp_profile_CB', 'rb') as pickle_file:
    smp_profile = pickle.load(pickle_file)


# def Test_debye_value(smp_profile):
CB_sp = []
temp = []

debyeV_out = np.zeros(5)
debyeH_out = np.zeros(5)

ErrV_ref = np.ones(5)*100
ErrH_ref = np.ones(5)*100

# sea_ice_density = 910
# icepack = make_ice_column(ice_type='firstyear',
#         thickness=[2], temperature=250, 
#         microstructure_model='independent_sphere',
#         radius=1e-3,
#         brine_inclusion_shape='spheres',
#         density=sea_ice_density,
#         salinity=32*PSU,
#         add_water_substrate=True)

name_prev = False
Calc_temp = True
from Functions_temp import Calc_temp_profiles
comb_info = Calc_temp_profiles(plot=False)
[AMSR, ECMWF] = Get_AMSR_ECMWF()

for (s, E) in zip(comb_info, ECMWF):
    ## plot results
    plt.figure()
    profile_count = -1
    for name, profile in smp_profile[s.site].items():
        profile_count +=1
        
        ## If we use measured temperatures
        if Calc_temp==False:
            temp_ice = s.Tcal_ice[profile_count] # 250 K -very wrong :) 
            temp_snow = profile.temperature
        ## Test output with calculated temp profiles
        else: ## If we use calculated temperature profile
            temp_ice = s.Tcal_ice[profile_count]
            temp_snow = s.Tcal_snow[profile_count]

        ## find optimal debye value
        debye_array = np.linspace(0.8, 4, 10)
        for debye in debye_array:
        
        # for name, profile in smp_profile[site].items():
            # Convert to exp. corr length
            lex_array = debye * 4 * (1 - profile.density / 917) / (profile.ssa * 917)
    
            ## Make icepack        
            sea_ice_density = 900
            icepack = make_ice_column(ice_type='firstyear',
                    thickness=[s.SIT_m], temperature=temp_ice, 
                    microstructure_model='exponential',
                    corr_length=0.15e-3,
                    brine_inclusion_shape='spheres',
                    density=sea_ice_density,
                    salinity=7*PSU,
                    add_water_substrate="ocean")
            
            snowpack = make_snowpack(profile.thick, microstructure_model='exponential',
                            ice_permittivity_model=ssp, density=profile.density , 
                            corr_length=lex_array, temperature = profile.temperature,
                            salinity=profile.salinity *PSU)
    
            medium = snowpack + icepack
            CB_sp.append(medium)
            temp = np.append(temp, np.mean(profile.temperature))
            
            
    
        #run simulations
        res = m.run(sensor, CB_sp)
        from roughness import roughness
        [TB_SMRT_V,TB_SMRT_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(E.tcwv,res,np.mean(s.Tcal_ice),site=s.site,several=True)

        ErrV = TB_SMRT_V-satdataV_landfree.flatten()
        print(ErrV)
        ErrH = TB_SMRT_H-satdataH_landfree.flatten()
        print(ErrH)
        for i in range(ErrV.size):
            if abs(ErrV[i]) < abs(ErrV_ref[i]):
                ErrV_ref[i] = ErrV[i]
                debyeV_out[i] = debye
        
            if abs(ErrH[i]) < abs(ErrH_ref[i]):
                ErrH_ref[i] = ErrH[i]
                debyeH_out[i] = debye

        np.save(s.site + 'Optimal_freq_debyeV', debyeV_out)
        np.save(s.site + 'Optimal_freq_debyeH', debyeH_out)
        

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
Calc_temp = True
from Functions_temp import Calc_temp_profiles
comb_info = Calc_temp_profiles(plot=False)
[AMSR, ECMWF] = Get_AMSR_ECMWF()

for (s, E) in zip(comb_info, ECMWF):
    ## find optimal debye value
    debye_array = np.linspace(0.8, 4, 10)
    for debye in debye_array:
        ## plot results
        plt.figure()
        profile_count = -1
        for name, profile in smp_profile[s.site].items():
            profile_count +=1
            
            ## If we use measured temperatures
            if Calc_temp==False:
                temp_ice = s.Tcal_ice[profile_count] # 250 K -very wrong :) 
                temp_snow = profile.temperature
            ## Test output with calculated temp profiles
            else: ## If we use calculated temperature profile
                temp_ice = s.Tcal_ice[profile_count]
                temp_snow = s.Tcal_snow[profile_count]
    
            ## Make icepack        
            sea_ice_density = 900
            icepack = make_ice_column(ice_type='firstyear',
                    thickness=[s.SIT_m], temperature=temp_ice, 
                    microstructure_model='exponential',
                    corr_length=0.15e-3,
                    brine_inclusion_shape='spheres',
                    density=sea_ice_density,
                    salinity=7*PSU,
                    add_water_substrate="ocean")
    
            # Convert to exp. corr length
            debye = 1.75
            lex_array = debye * 4 * (1 - profile.density / 917) / (profile.ssa * 917)
            
            ## Make snowpack
            snowpack = make_snowpack(profile.thick, microstructure_model='exponential',
                            ice_permittivity_model=ssp, density=profile.density , 
                            corr_length=lex_array, temperature = temp_snow,
                            salinity=profile.salinity *PSU)
    
            medium = snowpack + icepack
            CB_sp.append(medium)                     
        
        #% Run simulations
        # run simulations
        res = m.run(sensor, CB_sp)
        
        # apply roughness and atmosphere correction
        [TB_SMRT_V,TB_SMRT_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(E.tcwv,res,np.mean(s.Tcal_ice),site=s.site,several=True)
        
        #% make plots
        [data, dataErr] = load_ref_data(s.site)
        ## mean data
        data[4,:] =TB_SMRT_V
        data[5,:] = TB_SMRT_H
        data_comb[count2:count2+4,:] = data[2:6,:]
        
        
        ## Errors
        dataErr[4,:] = np.std(res.TbV(), axis=1)
        dataErr[5,:] = np.std(res.TbH(), axis=1)
        dataErr_comb[count2:count2+4,:] = dataErr[2:6,:]
        
        colors = ['g', 'r', 'purple', 'brown', 'b', 'orange']
        labels = ['', '', 'AMSR2 V landfree', 'AMSR2 H landfree', 'TBV SMRT', 'TBH SMRT']
        plot_freq(s, data, labels, cs=colors, err=dataErr, title='Tbs simulated at CB at ' + s.site + s.date)
        
        count2 += 4
    ## two colors
    #colors_comb = np.append(colors_comb,[['b', 'brown'] for i in range(int(data_comb.shape[0]/2))])
    ## 4 colors
    colors_comb = np.append(colors_comb,[[colors[2:6]] for i in range(int(data_comb.shape[0]/4))])
    plot_freq(s, data_comb, cs=colors_comb, labels=labels[2:6])
    #Test_debye_value(smp_profile)