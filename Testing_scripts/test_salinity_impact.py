# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 22:59:06 2022

@author: Ida Olsen

Test salinity
"""
import numpy as np
parrent_dir = "C:\\Users\\Ida Olsen\\Documents\\Speciale_2022"
import sys
import os
# sys.path.append(parrent_dir + '\\smrt')
import pickle
import matplotlib.pyplot as plt

from smrt.permittivity.saline_snow import  saline_snow_permittivity_scharien_with_stogryn95 as ssp
from smrt import make_ice_column, make_snowpack, sensor_list
from smrt.core.globalconstants import PSU
from smrt import make_model


#%% Make snowpit class
from Get_Sat_Data import Get_AMSR_ECMWF
from roughness import roughness


# pd.set_option('mode.chained_assignment',None)
#Get Data info
Base_path = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/"
with open(Base_path + 'Code/smrt_in-out/smp_profile_CB', 'rb') as pickle_file:
    smp_profile = pickle.load(pickle_file)

# Make new SMRT configuration with returns
m = make_model("iba", "dort")
frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
beamorientation = 55. # [degrees]
sensor = sensor_list.passive(frequencies, beamorientation)


CB_sp = []
count2 = 0
Calc_temp = True

from Functions_temp import Calc_temp_profiles
comb_info = Calc_temp_profiles(plot=False)
[AMSR, ECMWF] = Get_AMSR_ECMWF()

for (s, E, A) in zip(comb_info, ECMWF, AMSR):
    ## plot results
    plt.figure(figsize=(7,6))
    frequencies = np.array([6.925,10.65,18.7,23.8,36.5]);
    plot3,=plt.plot(frequencies, A.DataSetV.flatten(),'.-', linewidth = 2, color='k', label='AMSR2 TB')
    plot4,=plt.plot(frequencies, A.DataSetH.flatten(),'.-', linewidth = 2, color='k')
    #plot3,=plt.plot(frequencies, A.DataSetV.flatten(),'.-', color='b', label = 'AMSR2 V')            
    #plot4,=plt.plot(frequencies, A.DataSetH.flatten(),'.-', color='orange', label = 'AMSR2 H')
    
    ## find optimal debye value
    salinities = np.linspace(3.0, 15.0, 3)
    colors = ['C' + str(int(num)) for num in range(len(salinities))]
    for sal, c in zip(salinities, colors):
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
    
    
            # Convert to exp. corr length
            debye = 1.5
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
            
            ## Make snowpack
            salinity = np.array(profile.salinity)
            salinity[-1] = sal
            snowpack = make_snowpack(profile.thick, microstructure_model='exponential',
                            ice_permittivity_model=ssp, density=profile.density ,
                            corr_length=lex_array, temperature = temp_snow,
                            salinity=salinity *PSU)
    
            medium = snowpack + icepack
            CB_sp.append(medium)
    
        #run simulations
        res = m.run(sensor, CB_sp)
        [TB_SMRT_V,TB_SMRT_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(E.tcwv,res,np.mean(s.Tcal_ice),site=s.site,several=True)
         
        plot1,=plt.plot(frequencies, TB_SMRT_V,'.-', linewidth = debye, color=c, label = 'Bottom salinity [PSU] : ' + str(np.round(sal,2)))
        plot2,=plt.plot(frequencies, TB_SMRT_H,'.-', linewidth = debye, color=c)
    plt.ylim([200,260])
    plt.grid()
    plt.legend(loc='lower right',fontsize=12)
    plt.title('Salinity sensitivity study at site ' + s.site,fontsize=16)
    plt.xlabel('Frequency [GHz]',fontsize=14)
    plt.ylabel('TB [K]',fontsize=14)
    plt.savefig( os.getcwd() + '/figures/' + s.site + '_sal_comparisons.PNG', bbox_inches='tight')
    plt.show()