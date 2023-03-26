# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:28:34 2022

@author: Ida Olsen


The equations assume the interface between the media is flat and that the media are homogeneous and isotropic.
The incident light is assumed to be a plane wave.
source: wikipedia
"""
import numpy as np
## Reflectance p-polarized=vertical, s-polarized=horizontal

def R(theta_inc, eps_m1, eps_m2):
    """

    Parameters
    ----------
    theta_inc : FLOAT
        INCIDENCE ANGLE
    eps_m1 : float
        permitivity origin medium
    eps_m2 : float
        permitivity entering medium

    Returns
    -------
    None.

    """
    n1 = np.sqrt(eps_m1)
    n2 = np.sqrt(eps_m2)
    
    RV = ((n1*np.cos(theta_inc)-n2*np.sqrt(1-n1/n2*np.sin(theta_inc)**2))/(
            n1*np.cos(theta_inc)+n2*np.sqrt(1-n1/n2*np.sin(theta_inc)**2)))**2
    
    RH = ((n1*np.sqrt(1-n1/n2*np.sin(theta_inc)**2)-n2*np.cos(theta_inc))/(
            n1*np.sqrt(1-n1/n2*np.sin(theta_inc)**2)+n2*np.cos(theta_inc)))**2
    
    [TV,TH] = T(RV, RH)
    return [RV, RH, TV, TH]
    
def T(RV, RH):
    TV = 1 - RV
    TH = 1 - RH
    
    return [TV, TH]

def Snell(theta_inc, eps_m1, eps_m2):
    theta_r = theta_inc
    n1 = np.sqrt(eps_m1)
    n2 = np.sqrt(eps_m2)
    theta_t = np.arcsin(n1/n2*np.sin(theta_inc))
    
    return [theta_r*180/PI, theta_t*180/PI]

def Absorption(eps, k0):
    import cmath
    
    n_mm = -cmath.sqrt(eps).imag
    alpha = k0*n_mm
    return alpha

def Print_info(n, RV, RH, TV, TH, theta_r, theta_t, a):
    print('-----------Interface number',n,'------------')
    print('Reflectance H pol', RH)
    print('Reflectance V pol', RV)
    print('Transmittance H pol', TH)
    print('Transmittance V pol', TV)
    print('Angle of reflection', theta_r)
    print('Angle of refraction', theta_t)
    print('-----------------------------------------')
    
    print('Information for frequency 18.7 GHz')
    print('Absorption coefficient of layer Ka [Np/m]', 2*a)
    print('Eqvivalent penetration depth dp [cm]', 1/(2*a)*100, '\n')

## constants
PI = 3.14
eps_air_re = 1

## ice values
eps_i_re = 3 # sea ice permitivity FYI for these freq. around 3.0 - 4.0 (this is old FYI so probably closer to 3)
eps_i_im = -0.25 # range approx 0.1 -0.3
eps_i = complex(eps_i_re,eps_i_im)


## snow values
#eps geldsetzer: (1.858184392216508+0.015434684176853451j)
#eps scharien 71: (1.5818353514194177+0.006379949159705199j)
#eps scharien 95: (1.5707024528917344+0.013437192783908376j)

rho_s = 350  # kg/m^3
vi = rho_s/916.7  # volume factor (ice to air)
eps_s_re = (1+0.84*vi)/(1-0.42*vi) # permitivity of snow

## the penetration depth and attenuation in the snow is highly impacted by the salinity
eps_i_im_s = 0.07 # closer to salt free ice e.g. pure ice
eps_s_im = -(0.34*vi*eps_i_im_s)/(1-0.42*vi)**2

eps_s = complex(eps_s_re,eps_s_im)

## Ocean values
eps_w = complex(60,-60)

## penetration depth ice
Freqs = [6.9, 10.7, 18.7, 23.8, 36.5]
y0 = [(3*10**8)/(f*10**9) for f in Freqs]
k0 = [(2*PI)/y for y in y0]
dp_ice = [(np.sqrt(eps_i.real)/(k*abs(eps_i.imag)))*100 for k in k0] #convertion to cm

# print('Penetration depth into ice [cm]:', dp_ice[3])

## penetration depth snow
dp_s = [(np.sqrt(eps_s.real)/(k*abs(eps_s.imag)))*100 for k in k0] #convertion to cm

# print('Penetration depth into snow [cm]:', dp_s[3])

## --------------------------------------------------------------
# Layer 1

[RV, RH, TV, TH] = R(55/180*PI, eps_air_re, eps_s_re)
[theta_r, theta_t] = Snell(55/180*PI, eps_air_re, eps_s_re)

a = Absorption(eps_s, k0[3])
Print_info(str(1)+' Air/Snow',RV, RH, TV, TH, theta_r, theta_t, a)
# print('dp [cm]:', 1/(2*a)*100)
## --------------------------------------------------------------
# Layer 2

[RV, RH, TV, TH] = R(theta_t/180*PI, eps_s_re, eps_i_re)
[theta_r, theta_t] = Snell(theta_t/180*PI, eps_s_re, eps_i_re)
a = Absorption(eps_i, k0[3])

Print_info(str(2)+' Snow/Ice',RV, RH, TV, TH, theta_r, theta_t, a)

## --------------------------------------------------------------
# Layer 3

[RV, RH, TV, TH] = R(theta_t/180*PI, eps_i.real, eps_w.real)
[theta_r, theta_t] = Snell(theta_t/180*PI, eps_i.real, eps_w.real)
a = Absorption(eps_w, k0[3])

Print_info(str(3)+ ' Ice/Ocean',RV, RH, TV, TH, theta_r, theta_t, a)

## --------------------------------------------------------------
# Layer 2 return

[RV, RH, TV, TH] = R(theta_r/180*PI, eps_i.real, eps_s.real)
[theta_r, theta_t] = Snell(theta_r/180*PI, eps_i.real, eps_s.real)
a = Absorption(eps_w, k0[3])

Print_info(str(2)+ ' Ice/Snow return',RV, RH, TV, TH, theta_r, theta_t, a)

## --------------------------------------------------------------
# Layer 1 return

[RV, RH, TV, TH] = R(theta_t/180*PI, eps_s.real, eps_air_re)
[theta_r, theta_t] = Snell(theta_t/180*PI, eps_s.real, eps_air_re)
a = Absorption(eps_s, k0[3])

Print_info(str(1)+ ' Snow/Air return',RV, RH, TV, TH, theta_r, theta_t, a)