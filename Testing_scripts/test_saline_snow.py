# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:19:18 2022

@author: Ida Olsen
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import cmath
import numpy as np
from numpy import sqrt, cos, sin, exp
# from mpmath import sec
from scipy.constants import c

import os
parrent_dir = "C:\\Users\\Ida Olsen\\Documents\\Speciale_2022"
import sys
sys.path.append(parrent_dir + '\\smrt')
from smrt.core.globalconstants import PSU
from smrt import make_ice_column, make_interface, make_snowpack, sensor_list
from smrt import make_model

from smrt.permittivity.saline_snow import saline_snow_permittivity_geldsetzer09
from smrt.permittivity.saline_snow import saline_snow_permittivity_scharien_with_stogryn71
from smrt.permittivity.saline_snow import saline_snow_permittivity_scharien_with_stogryn95

def test_saline_snow(frequency,density, temperature, salinity):
    
    eps_g = saline_snow_permittivity_geldsetzer09(frequency, density, temperature, salinity)
 
    print("eps geldsetzer:", np.round(eps_g,3))
    eps_s71 = saline_snow_permittivity_scharien_with_stogryn71(frequency, density, temperature, salinity)
    print("eps scharien 71:",  np.round(eps_s71,3))
    eps_s95 = saline_snow_permittivity_scharien_with_stogryn95(frequency, density, temperature, salinity) 
    print("eps scharien 95:",  np.round(eps_s95,3))
    return eps_s95
def test_saline_snow_TB():    
    Base_path = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/"
    with open(Base_path + 'Code/smrt_in-out/smp_profile_CB', 'rb') as pickle_file:
        smp_profile = pickle.load(pickle_file)
    
    # Make new SMRT configuration with returns
    m = make_model("iba", "dort")
    frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
    beamorientation = 55. # [degrees]
    sensor = sensor_list.passive(frequencies, beamorientation)
    
    
    Calc_temp = True
    from Functions_temp import Calc_temp_profiles
    comb_info = Calc_temp_profiles()
    
    count_sal_model = -1
    
    CB_sp = []
    
    names = ['Geldsetzer', 'sharien71', 'sharien95']
    colors = ['b', 'g', 'orange']
    sea_ice_density = 900
    for ssp in [saline_snow_permittivity_geldsetzer09, saline_snow_permittivity_scharien_with_stogryn71, saline_snow_permittivity_scharien_with_stogryn95]:
        count_sal_model += 1
        c = colors[count_sal_model]
        sites = ['AK1', 'AK2', 'AK3', 'AK4']
        for site in sites:
            for i in range(len(sites)):
                if comb_info[i].site==site:
                    s = comb_info[i]
            ## plot results
            # plt.figure()
            profile_count = -1
            for name, profile in smp_profile[site].items():
                profile_count +=1
                
                ## If we use measured temperatures
                if Calc_temp==False:
                    temp_ice = s.Tcal_ice[profile_count] # 250 K -very wrong :) 
                    temp_snow = profile.temperature
                ## Test output with calculated temp profiles
                else: ## If we use calculated temperature profile
                    temp_ice = s.Tcal_ice[profile_count]
                    temp_snow = s.Tcal_snow[profile_count]
    
                
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
        
                snowpack = make_snowpack(profile.thick, microstructure_model='exponential',
                                ice_permittivity_model=ssp, density=profile.density , 
                                corr_length=lex_array, temperature = temp_snow,
                                salinity=profile.salinity*PSU)
        
                medium = snowpack + icepack
                CB_sp.append(medium)
    
    
        # run simulations
        res = m.run(sensor, CB_sp)
        
        ## Plot results
        frequencies =np.array([6.925,10.65,18.7,23.8,36.5]);
        plt.plot(frequencies, np.mean(res.TbV(), axis=1),'.-', color=c, label = 'TBV SMRT' + names[count_sal_model])
        plt.plot(frequencies, np.mean(res.TbH(), axis=1),'.-', color=c, label = 'TBH SMRT'+ names[count_sal_model])
        plt.ylim([200,260])
        plt.grid()
        plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
        plt.title('Salinity models test',fontsize=16)
        plt.xlabel('freq GHz',fontsize=14)
        plt.ylabel('TB',fontsize=14)
    plt.savefig(os.getcwd() + '/figures/' + s.site + '_saline_snow_model_test.PNG', bbox_inches='tight')
    plt.show()

# test_saline_snow_TB()
# test_saline_snow()

## Testing microstructure values
def test_debye_vals():
    from smrt.permittivity.saline_snow import saline_snow_permittivity_scharien_with_stogryn95 as ssp
    Base_path = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/"
    with open(Base_path + 'Code/smrt_in-out/smp_profile_CB', 'rb') as pickle_file:
        smp_profile = pickle.load(pickle_file)
    
    # Make new SMRT configuration with returns
    m = make_model("iba", "dort")
    frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
    beamorientation = 55. # [degrees]
    sensor = sensor_list.passive(frequencies, beamorientation)
    
    
    Calc_temp = True
    from Functions_temp import Calc_temp_profiles
    comb_info = Calc_temp_profiles()
    
    count_sal_model = -1
    
    CB_sp = []
    
    debye_vals = np.linspace(0.5,3.0,5)
    colors = ['b', 'g', 'orange', 'red', 'k']
    sea_ice_density = 900
    for debye in debye_vals:
        count_sal_model += 1
        c = colors[count_sal_model]
        sites = ['AK1', 'AK2', 'AK3', 'AK4']
        
        for site in sites:
            plt.figure()
            for i in range(len(sites)):
                if comb_info[i].site==site:
                    s = comb_info[i]
            profile_count = -1
            for name, profile in smp_profile[site].items():
                profile_count +=1
                
                ## If we use measured temperatures
                if Calc_temp==False:
                    temp_ice = s.Tcal_ice[profile_count] # 250 K -very wrong :) 
                    temp_snow = profile.temperature
                ## Test output with calculated temp profiles
                else: ## If we use calculated temperature profile
                    temp_ice = s.Tcal_ice[profile_count]
                    temp_snow = s.Tcal_snow[profile_count]
    
                icepack = make_ice_column(ice_type='firstyear',
                        thickness=[s.SIT_m], temperature=temp_ice, 
                        microstructure_model='exponential',
                        corr_length=0.15e-3,
                        brine_inclusion_shape='spheres',
                        density=sea_ice_density,
                        salinity=7*PSU,
                        add_water_substrate="ocean")
        
                # Convert to exp. corr length
               #  debye = 1.75
                lex_array = debye * 4 * (1 - profile.density / 917) / (profile.ssa * 917)
        
                snowpack = make_snowpack(profile.thick, microstructure_model='exponential',
                                ice_permittivity_model=ssp, density=profile.density , 
                                corr_length=lex_array, temperature = temp_snow,
                                salinity=profile.salinity*PSU)
        
                medium = snowpack + icepack
                CB_sp.append(medium)
    
    
            # run simulations
            res = m.run(sensor, CB_sp)
            
            ## Plot results
            frequencies =np.array([6.925,10.65,18.7,23.8,36.5]);
            plt.plot(frequencies, np.mean(res.TbV(), axis=1),'.-', color=c, label = 'TBV SMRT debye val: ' + str(debye))
            plt.plot(frequencies, np.mean(res.TbH(), axis=1),'.-', color=c, label = 'TBH SMRT debye val:'+ str(debye))
            plt.ylim([200,260])
            plt.grid()
            plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
            plt.title('Debye Values Test',fontsize=16)
            plt.xlabel('freq GHz',fontsize=14)
            plt.ylabel('TB',fontsize=14)
            plt.savefig(os.getcwd() + '/figures/' + s.site + '_debye_vals_test.PNG', bbox_inches='tight')
    plt.show()

def snells_law(eps1, eps2, theta1):
    # print('theta to layer ' + str(n) + ':', np.rad2deg(np.arcsin(sin(theta1)*(sqrt(eps0)/sqrt(eps1)))))
    return np.arcsin(sin(theta1)*(sqrt(eps2)/sqrt(eps1)))

def reflectivityH(eps1, eps2, theta1):
    theta2 = snells_law(eps1, eps2, theta1)
    nom = sqrt(eps1)*cos(theta1)-sqrt(eps2)*cos(theta2)
    denom = sqrt(eps1)*cos(theta1)+sqrt(eps2)*cos(theta2)
    R = (nom/denom).real**2+(nom/denom).imag**2
    return R, theta2


def reflectivityV(eps1, eps2, theta1):
    theta2 = snells_law(eps1, eps2, theta1)
    nom = sqrt(eps2)*cos(theta1)-sqrt(eps1)*cos(theta2)
    denom = sqrt(eps2)*cos(theta1)+sqrt(eps1)*cos(theta2)
    R = (nom/denom).real**2+(nom/denom).imag**2
    return R, theta2
    
def emission_specular_surface(Reflectivity):
    return 1-Reflectivity

f1 = 7*10e9 # Hz 
d1 = 200 # kg/m^3
t1 = 260 # K
s1  = 13*PSU # PSU

f0 = 7*10e9 # Hz 
d0 = 350 # kg/m^3
t0 = 260 # K
s0  = 0.5*PSU # PSU

fice = 7*10e9 # Hz 
dice = 900 # kg/m^3
tice = 265 # K
sice  = 7*PSU # PSU
    
# test_debye_vals()
# epsice = test_saline_snow(fice,dice,tice,sice)
epsice = 3.485902483144466+0.048482768955471656j
print(epsice)
eps1 = test_saline_snow(f1,d1,t1,s1)
print(eps1)
# eps1 = 1.2558465956229954+0.01724856436252578j
# eps1 = 1.01633713281697062+0.01001195747407692j #270
#eps0 = test_saline_snow(f0,d0,t0,s0)

# theta0 = snells_law(eps1.real, eps0.real, theta1, 0)

refH, to = reflectivityH(epsice, eps1, np.deg2rad(44))
print('Ref V:', refH)
print('TV', 1-refH)
refV, to = reflectivityV(epsice, eps1, np.deg2rad(44))
print('Ref H:', refV)
print('TH:', 1-refV)

import numpy as np
import matplotlib.pyplot as plt


plt.figure(figsize=(7,6))
# for theta in np.linspace(40, 60, 2):
theta=50
rH=[]
rV=[]
tV=[]
tH=[]
tVV=[]
tHH=[]
tVVV=[]
tHHH=[]
for dens in np.linspace(150, 400, 10):
    f1 = 7*10e9 # Hz 
    d1 = dens # kg/m^3
    t1 = 260 # K
    s1  = 8.6*PSU # PSU

    eps0 = test_saline_snow(f1,349,256,0.9*PSU)
    eps1 = test_saline_snow(f1,d1,t1,s1)
    thetaice = np.deg2rad(theta)
    refH, theta_out = reflectivityH(epsice, eps1, thetaice)
    refHH, theta_out2 = reflectivityH(eps1, eps0, theta_out)
    refHHH, to = reflectivityH(eps0, complex(1,0), theta_out2)
    TH = 1-refH
    THHH = 1-refHHH
    THH = 1-refHH
    #print('Ref V:', refH)
    #print('TV', 1-refH)
    refV, to = reflectivityV(epsice, eps1, thetaice)
    refVV, theta_out2 = reflectivityV(eps1, eps0, theta_out)
    refVVV, to = reflectivityV(eps0, complex(1,0), theta_out2)
    TV = 1-refV
    TVV = 1-refVV
    TVVV = 1-refVVV
    rH.append(refH)
    rV.append(refV)
    tV.append(TV)
    tH.append(TH)
    tVV.append(TVV)
    tVVV.append(TVVV)
    tHH.append(THH)
    tHHH.append(THHH)

Density = np.linspace(150, 400, 10)
# plt.ylim(0.6,1.0)
plt.figure(figsize=(7,6))
#plt.plot(Density,tHHH, '-o', label='V air/snow')
plt.plot(Density,tHH, '-o', label='V snow/snow')
plt.plot(Density,tH, '-o', label='V snow/ice')
plt.legend(loc='lower right', fontsize=12)
plt.ylim([0.98, 1.004])
plt.xlabel('Density [kg/m3]', fontsize=12)
plt.ylabel('Transmissivity', fontsize=12)
plt.title('Transmissivity across snow/ice, snow/snow and air/snow interfaces', fontsize=13)
plt.grid()
plt.show()

plt.figure(figsize=(7,6))
# plt.plot(Density,tVVV, '-o', label='H air/snow')
plt.plot(Density,tVV, '-o', label='H snow/snow')
plt.plot(Density,tV, '-o', label='H snow/ice')
plt.legend(loc='best', fontsize=12)
plt.ylim([0.85, 1.005])
plt.xlabel('Density [kg/m3]', fontsize=12)
plt.ylabel('Transmissivity', fontsize=12)
plt.title('Transmissivity across snow/ice, snow/snow and air/snow interfaces', fontsize=13)
plt.grid()
plt.show()



plt.figure(figsize=(7,6))
# for theta in np.linspace(40, 60, 2):
theta=50
rH=[]
rV=[]
tV=[]
tH=[]
tVV=[]
tHH=[]
tVVV=[]
tHHH=[]
for dens in np.linspace(150, 400, 10):
    f1 = 7*10e9 # Hz 
    d1 = dens # kg/m^3
    t1 = 260 # K
    s1  = 8.6*PSU # PSU

    eps0 = test_saline_snow(f1,dens,256,0.9*PSU)
    eps1 = test_saline_snow(f1,338,t1,s1)
    thetaice = np.deg2rad(theta)
    refH, theta_out = reflectivityH(epsice, eps1, thetaice)
    refHH, theta_out2 = reflectivityH(eps1, eps0, theta_out)
    refHHH, to = reflectivityH(eps0, complex(1,0), theta_out2)
    TH = 1-refH
    THHH = 1-refHHH
    THH = 1-refHH
    #print('Ref V:', refH)
    #print('TV', 1-refH)
    refV, to = reflectivityV(epsice, eps1, thetaice)
    refVV, theta_out2 = reflectivityV(eps1, eps0, theta_out)
    refVVV, to = reflectivityV(eps0, complex(1,0), theta_out2)
    TV = 1-refV
    TVV = 1-refVV
    TVVV = 1-refVVV
    rH.append(refH)
    rV.append(refV)
    tV.append(TV)
    tH.append(TH)
    tVV.append(TVV)
    tVVV.append(TVVV)
    tHH.append(THH)
    tHHH.append(THHH)

Density = np.linspace(150, 400, 10)
#plt.scatter(Density,rH)
#plt.scatter(Density,rV)
# plt.ylim(0.6,1.0)
plt.figure(figsize=(7,6))
plt.plot(Density,tHHH, '-o', label='V air/snow')
plt.plot(Density,tHH, '-o', label='V snow/snow')
#plt.plot(Density,tH, '-o', label='V snow/ice')
plt.legend(loc='lower right', fontsize=12)
plt.ylim([0.98, 1.004])
plt.xlabel('Density [kg/m3]', fontsize=12)
plt.ylabel('Transmissivity', fontsize=12)
plt.title('Transmissivity across snow/ice, snow/snow and air/snow interfaces', fontsize=13)
plt.grid()
plt.show()




plt.figure(figsize=(7,6))
plt.plot(Density,tVVV, '-o', label='H air/snow')
plt.plot(Density,tVV, '-o', label='H snow/snow')
#plt.plot(Density,tV, '-o', label='H snow/ice')
plt.legend(loc='best', fontsize=12)
plt.ylim([0.85, 1.005])
plt.xlabel('Density [kg/m3]', fontsize=12)
plt.ylabel('Transmissivity', fontsize=12)
plt.title('Transmissivity across snow/ice, snow/snow and air/snow interfaces', fontsize=13)
plt.grid()
plt.show()

#%%
plt.figure(figsize=(7,6))

rH=[]
rV=[]
tV=[]
tH=[]
tVV=[]
tHH=[]
thetaice = np.deg2rad(50)
sal = np.array([0.5, 0.7, 1, 2, 3, 5, 7, 10, 15, 20])*PSU
for s in sal:
    f1 = 7*10e9 # Hz 
    d1 = 338 # kg/m^3
    d0 = 348
    t1 = 260 # K
    t0 = 257
    s1  = s # PSU

    eps0 = test_saline_snow(f1,d0,t0,s1)
    eps1 = test_saline_snow(f1,d1,t1,s1)
    thetaice = np.deg2rad(theta)
    refH, theta_out = reflectivityH(epsice, eps1, thetaice)
    refHH, theta_out2 = reflectivityH(eps1, eps0, theta_out)
    TH = 1-refH
    THH = 1-refHH
    #print('Ref V:', refH)
    #print('TV', 1-refH)
    refV, theta_out = reflectivityV(epsice, eps1, thetaice)
    refVV, theta_out2 = reflectivityV(eps1, eps0, theta_out)
    TV = 1-refV
    TVV = 1-refVV
    rH.append(refH)
    rV.append(refV)
    tV.append(TV)
    tH.append(TH)
    tVV.append(TVV)
    tHH.append(THH)

#plt.scatter(Density,rH)
#plt.scatter(Density,rV)
# plt.ylim(0.6,1.0)
plt.plot(sal*1000,tH, '-o', label='V')
plt.plot(sal*1000,tV, '-o', label='H')
plt.plot(sal*1000,tHH, '-o', label='V')
plt.plot(sal*1000,tVV, '-o', label='H')
plt.legend(loc='best', fontsize=12)
plt.xlabel('Salinity [PSU]', fontsize=12)
plt.ylabel('Transmissivity', fontsize=12)
plt.title('Transmissivity across snow/ice interface', fontsize=13)
# plt.ylim([0.82,1.01])
plt.grid()
plt.show()


f = 18e9
plt.figure(figsize=(6,6))
sal = np.array([0.5, 0.7, 1, 2, 3, 5, 7, 10, 15, 20])*PSU
D = []
for s1 in sal:
    eps = test_saline_snow(f,338,260,s1)
    D.append((3e8*sqrt(eps.real))/(f*2*np.pi*eps.imag))
plt.plot(sal*1000, np.array(D)*100,'o-', label = 'D:' + str(np.round(D*100,1)))
plt.grid()
plt.ylim([0,80])
#plt.legend()
plt.xlabel('salinity [PSU]', fontsize=12)
plt.ylabel('Penetration depth [cm]', fontsize=12)
plt.title('Penetration depth at 18GHz into saline bottom layer', fontsize=13)



plt.figure(figsize=(6,6))
freq = np.array([6.9, 10.7, 18.7, 23.8, 36.5])*1e9
D=[]
for f in freq:
    # eps = test_saline_snow(f,270,264,13*PSU)
    eps = test_saline_snow(f,338,260,9*PSU)
    D.append((3e8*np.sqrt(eps.real))/(f*2*np.pi*eps.imag))
plt.plot(freq*1e-9, np.array(D)*100, 'o-', label = 'D:' + str(np.round(D*100,1)))
plt.grid()
plt.ylim([0,20])
#plt.legend()
plt.xlabel('Frequency [GHz]')
plt.ylabel('Penetration depth [cm]')
plt.title('Penetration depth into saline bottom snow layer')


plt.figure(figsize=(6,6))
freq = np.array([6.9, 10.7, 18.7, 23.8, 36.5])*1e9
D=[]
for f in freq:
    eps = test_saline_snow(f,350,261,0.5*PSU)
    D.append((3e8*np.sqrt(eps.real))/(f*2*np.pi*eps.imag))
plt.plot(freq*1e-9, np.array(D)*100, 'o-', label = 'D:' + str(np.round(D*100,1)))
plt.grid()
#plt.legend()
plt.xlabel('Frequency [GHz]')
plt.ylabel('Penetration depth [cm]')
plt.title('Penetration depth into center slightly saline snow layer')
# rH=[]
# rV=[]
# tV=[]
# tH=[]
# for sal in np.linspace(5, 20, 10):
#     f1 = 7*10e9 # Hz 
#     d1 = 270 # kg/m^3
#     t1 = 260 # K
#     s1  = sal*PSU # PSU

#     eps1 = test_saline_snow(f1,d1,t1,s1)
#     thetaice = np.deg2rad(44)
#     refH = reflectivityH(epsice, eps1, thetaice)
#     TH = 1-refH
#     #print('Ref V:', refH)
#     #print('TV', 1-refH)
#     refV = reflectivityV(epsice, eps1, thetaice)
#     TV = 1-refV
#     rH.append(refH)
#     rV.append(refV)
#     tV.append(TV)
#     tH.append(TH)

# Salinity = np.linspace(1, 20, 10)
# #plt.scatter(Density,rH)
# #plt.scatter(Density,rV)
# # plt.ylim(0.6,1.0)
# plt.scatter(Salinity,tH, label='TV, ')
# plt.scatter(Salinity,tV, label='TH, ')
# plt.legend(loc='best')
# plt.xlabel('Salinity [PSU]')
# plt.ylabel('Transmissivity')
# plt.grid()
# plt.show()
