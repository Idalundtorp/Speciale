# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:01:51 2023

@author: Ida Olsen
"""
import numpy as np
import re
import matplotlib.pyplot as plt

## Import SMRT packages
from smrt.permittivity.saline_snow import saline_snow_permittivity_scharien_with_stogryn95 as ssp
from smrt import make_ice_column, make_snowpack, sensor_list
from smrt.core.globalconstants import PSU
from smrt import make_model
from smrt import SMRTError

## self made modules
from roughness import roughness

def temp_profile(self, plot=False):
    '''
    This function computes a two piece linear temperature profile using ...
    The snow/ice interface temperature is computed based on ...

    Parameters
    ----------
    k : float
        Thermal conductivity of ice (MYI or FYI)
    SD : float
        total snow depth
    ice : float
        sea ice thickness
    tsk : float
        surface temprature.
    nl_snow : int
        Number of layers in snow
    nl_ice : int
        Number of layers in ice

    Returns
    -------
    list
        Snow/ice interface temperature, temperature profile of snow and of ice.
    '''

    ## Ocean temperature
    T_w = (273.15-1.8)  # K
    ks = 0.25 # Thermal conductivity of snow [w/(m*K)]
    
    ## snow ice interface temperature
    self.T_sii = (self.k * self.sd_tot * T_w + ks * self.sit *
                self.Tsur) / (self.k * self.sd_tot + ks * self.sit)
     
    from sklearn.linear_model import LinearRegression
    
    # print(self.sd_tot)
    # print('Tsii')
    # print(self.T_sii)
    
    x = [0, np.array(self.sd_tot)]
    y = [np.array(self.Tsur), np.array(self.T_sii)]
    # print('x and y')
    # print(x)
    # print(y)
    model = LinearRegression().fit(np.array(x).reshape(-1,1), np.array(y))
    
    # if self.nl_snow == 2:
        # xx = np.array([self.sd[0]/2, self.sd[0] + self.sd[1]/2]).reshape(-1,1)
    if self.nl_snow == 2:
        xx = np.array([self.sd[0]/2, self.sd[0] + self.sd[1]/2]).reshape(-1,1)
    else:
        xx = np.array([self.sd[0]/2, self.sd[0] + self.sd[1]/2, self.sd[0] + self.sd[1] + self.sd[2]/2]).reshape(-1,1)
    # print('Snow array', xx)
    # else: # snow thinner than 3cm one layer
    # xx = np.array([self.sd/2]).reshape(-1,1)
    self.Tsnow = model.predict(xx)
    

    if plot==True:
        # ice/snow thickness arrays for showing temperature profiles    
        ice_thick = np.array([0, -self.sit])
        t_ice = np.array([self.T_sii, T_w])
        snow_thick = [np.sum(self.sd), self.sd[0] + self.sd[1]/2, self.sd[0]/2, 0]
        t_snow = [self.Tsur] + [Ts for Ts in self.Tsnow] + [self.T_sii]

        plt.figure()
        plt.plot(t_snow, snow_thick,'.-')
        plt.plot(t_ice, ice_thick,'.-')
        plt.grid()
        plt.xlabel('Temperature [K]')
        plt.ylabel('Depth [m]')

class OE_data():
    def __init__(self, icetype):
        # super().__init__(icetype)
        self.sd = []
        self.sit = []
        self.T_sii = []
        self.Tsur = []
        self.Tsnow = []
        self.k = 2.10
        self.nl_snow = 2
        self.nl_ice = 1
    def Calc_temp(self):
        temp_profile(self, plot=False)



def forward_OIB_AMSR2(X, 
                      
                      self,
                      variables,
                      nl_snow,
                      n, pol,
                      logarithmic=True,
                      avg=False,
                      sim=False,
                      know_SD=False,
                      know_SIT=False):
    """
    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    variables : list
        A list of variables which are being optimized
    
    Returns
    -------
    TYPE
        DESCRIPTION.
    
    """

    ## OIB data - prepared inputs for OE (based on average and std of OIB data)
    FYI_data = self 
    # OE class to compute temperature profile
    OE = OE_data('FYI')
    
    # string of variables
    variables = [re.sub(r'\d+','',v) for v in variables]
    
    
    ## Snow variables
    ind_SD = [i for i, j in enumerate(variables) if j == 'SD']
    ind_sal = [i for i, j in enumerate(variables) if j == 'sal']
    ind_rho = [i for i, j in enumerate(variables) if j == 'rho']
    ind_lex = [i for i, j in enumerate(variables) if j == 'lex']
    
    
    # Snow Depth
    SD = np.zeros(nl_snow)
    # Correlation length
    lex = np.zeros(nl_snow)
    # salinity
    sal = np.zeros(nl_snow)
    # Density
    rho = np.zeros(nl_snow)
    
    if logarithmic==True:
        SD_tot = 10**X[ind_SD] if 'SD' in variables else np.mean(np.sum(FYI_data, axis=1))
        lex[0] = 10**X[ind_lex] if 'lex' in variables else 0.24*1e-3
        # sal[2] = 10**X[ind_sal] if 'sal' in variables else 13*PSU
    else:
        # Snow Depth
        SD_tot = X[ind_SD] if 'SD' in variables else np.mean(np.sum(FYI_data, axis=1))
        lex[0] = X[ind_lex] if 'lex' in variables else 0.24*1e-3
            
    ### the case where I link variables
    if nl_snow==3:             
        # Make 3 snow layers with 25% DH
        SD[0] = SD_tot[0]-SD_tot[0]*0.25
        SD[1] = SD_tot[0]*0.25 - 0.01
        SD[2] = 0.01
    
        # Correlation length with DH (0.1 mm higher in bottom)
        lex[1] = lex[0] + 0.13*1e-3
        lex[2] = lex[1] + 0.13*1e-3
        
        # salinity low salinity in top layers
        sal[2] = X[ind_sal] if 'sal' in variables else 13*PSU
        sal[0] = 0.001*PSU
        sal[1] = 0.5*PSU
        
        # density - WS followed by 2 layers NS-DH
        # rho = X[ind_rho] if 'rho' in variables else [350, 270, 270]
        rho[1] = X[ind_rho] if 'rho' in variables else 270
        rho[0] = 350
        rho[2] = rho[1]

    ### the case where i link variables
    if nl_snow==2:
             
        # Make 2 snow layers with 25% DH
        SD[0] = SD_tot[0]-SD_tot[0]*0.25
        SD[1] = SD_tot[0]*0.25
    
        # Correlation length with DH (0.1 mm higher in bottom)
        lex[1] = lex[0] + 0.15*1e-3
        # lex[1] = 0.204e-3
        
        ##
        #self.sd_rho_sim = [335, 270]
        #self.sd_lex_sim = [0.223e-3, 0.204e-3]
        #self.sd_sal_sim = [0.01*PSU, 10*PSU]
        
        # salinity low salinity in top layers
        sal[1] = X[ind_sal] if 'sal' in variables else 6*PSU
        # sal[1] = X[ind_sal] if 'sal' in variables else 10*PSU
        sal[0] = 0.02*PSU
        
        # density - WS followed by 2 layers NS-DH
        # rho = X[ind_rho] if 'rho' in variables else [350, 270, 270]
        rho[1] = X[ind_rho] if 'rho' in variables else 270
        rho[0] = 350
        # rho[0] = 335



    # Temperature of snow surface (this variable is assumed to be known e.g. from ECMWF)
    if avg==True:
        Tsur = X[variables.index('Tsur')] if 'Tsur' in variables else np.mean(FYI_data.skt)
    else:
        Tsur = X[variables.index('Tsur')] if 'Tsur' in variables else FYI_data.skt[n]

    ## Ice variables
    
    # Sea ice thickness
    if avg==True:
        SIT = X[variables.index('SIT')] if 'SIT' in variables else np.mean(FYI_data.sit)
    else:
        SIT = X[variables.index('SIT')] if 'SIT' in variables else FYI_data.sit[n]
        
    #print(SD)
    OE.nl_snow=nl_snow
    OE.sd = SD
    OE.sd_tot = np.sum(SD)
    OE.Tsur = Tsur # air temperature
    OE.sit = SIT
    OE.Calc_temp()
    # print('---------------------------------')
    # print(OE.Tsnow)
    # print(OE.T_sii)
    # print('---------------------------------')
    
    
    # print('---------------------------------')   
    # print(SD)
    # print(lex)
    # print(sal)
    # print(rho)
    # print('---------------------------------')
    
    ## SMRT SETUP
    
    ## Make icepack           
    sea_ice_density = 900
    icepack = make_ice_column(ice_type='firstyear',
            thickness=[SIT], 
            temperature=OE.T_sii, 
            microstructure_model='exponential',
            corr_length=0.15e-3,
            brine_inclusion_shape='spheres',
            density=sea_ice_density,
            salinity=7*PSU,
            add_water_substrate="ocean")
    
    ## Make snowpack
    snowpack = make_snowpack(SD, 
                    microstructure_model='exponential',
                    ice_permittivity_model=ssp, 
                    density=rho, 
                    corr_length=lex, 
                    temperature = OE.Tsnow,
                    salinity=sal)
    
    medium = snowpack + icepack
    
    
    m = make_model("iba", "dort")
    
    frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
    
    beamorientation = 55. # [degrees]
    sensor = sensor_list.passive(frequencies, beamorientation)
    
    #% Run simulations
    res = m.run(sensor, medium)
    
    # #% adjust for atmospheric influence

    m = X[variables.index('slope')] if 'slope' in variables else 0.19
    # print(m)
    
    if avg==True:
        [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(np.mean(FYI_data.tcwv),res,OE.T_sii, m=m)
    else:
        [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(FYI_data.tcwv[n],res,OE.T_sii, m=m)
    
    if sim==False:
        if nl_snow==2:
            Bias_TBV = [0.99187069, 2.02032188, 0.80001409, -0.23933172, -0.13909897]
            Bias_TBH = [2.71943228, 3.55059095, -0.3034888,  -1.31527995, -3.17481278]
        elif nl_snow==3:
            Bias_TBV = [0.97801158, 2.11916614,  1.05152169,  0.00541775, -0.18454027]
            Bias_TBH = [1.95766929, 3.41963567,  0.34257621, -0.53548227, -2.51303051]
            Bias_TBV =np.array(Bias_TBV) + np.array([0.07013689, 0.08566306, -0.05492397, -0.26809489, -0.96714545])
            Bias_TBH =np.array(Bias_TBH) + np.array([0.42418694, 0.58594604, 0.70326948, 0.58685952, 0.04449245])
            
            #Bias_TBV = [0.546, 1.753,  0.749, -0.020, 0.088]
            #Bias_TBH = [1.565, 3.240, 0.454, -0.097, -1.803]
            # Bias_TBV = [ 0, 0, 0, 0, 0]
        # Bias_TBV = [ 0.43621623,1.67760632,0.63745003,-0.46549513,-1.40957493]
        # Bias_TBV = [ 0.79518149 , 1.95518354,  0.78687942, -0.44976116, -1.73411494  ]
        # Bias_TBH = [ 2.12858105, 3.72325291,  0.72247579, -0.28446336, -3.19752185]
        # Bias_TBV = [ 0, 0, 0, 0, 0]
        # Bias_TBH = [ 0, 0, 0, 0, 0]  

    elif sim==True:
        Bias_TBV = [ 0, 0, 0, 0, 0]
        Bias_TBH = [ 0, 0, 0, 0, 0]   

    if pol=='V':
        return TB_AMSR_V + Bias_TBV
    elif pol=='VH0':
        y_out = np.array([TB_AMSR_V + Bias_TBV])
        return np.append(y_out, TB_AMSR_H[0]+Bias_TBH[0])
    else:
        return np.concatenate(([TB_AMSR_V + Bias_TBV,TB_AMSR_H + Bias_TBH])) # np.concatenate((res.TbV(),res.TbH()))
    
def forward_OIB_AMSR2_new(X, 
                      
                      self,
                      variables,
                      nl_snow,
                      n, pol,
                      logarithmic=False,
                      avg=False):
    """
    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    variables : list
        A list of variables which are being optimized
    
    Returns
    -------
    TYPE
        DESCRIPTION.
    
    """


    # OE class to compute temperature profile
    OE = OE_data('FYI')
    
    # string of variables
    variables = [re.sub(r'\d+','',v) for v in variables]
    
    
    ## Snow variables
    ind_SD = [i for i, j in enumerate(variables) if j == 'SD']
    ind_sal = [i for i, j in enumerate(variables) if j == 'sal']
    ind_rho = [i for i, j in enumerate(variables) if j == 'rho']
    ind_lex = [i for i, j in enumerate(variables) if j == 'lex']
    
    
    # Snow Depth
    SD = np.zeros(nl_snow)
    # Correlation length
    lex = np.zeros(nl_snow)
    # salinity
    sal = np.zeros(nl_snow)
    # Density
    rho = np.zeros(nl_snow)
    
    if logarithmic==True:
        SD_tot = 10**X[ind_SD] if 'SD' in variables else np.mean(self.sd)
        lex[0] = 10**X[ind_lex] if 'lex' in variables else 0.24*1e-3
        # sal[2] = 10**X[ind_sal] if 'sal' in variables else 13*PSU
    else:
        # Snow Depth
        SD_tot = X[ind_SD] if 'SD' in variables else np.mean(self.sd)
        lex[0] = X[ind_lex] if 'lex' in variables else 0.24*1e-3
            
    ### the case where I link variables
    if nl_snow==3:             
        # Make 3 snow layers with 25% DH
        SD[0] = SD_tot[0]-SD_tot[0]*0.25
        SD[1] = SD_tot[0]*0.25 - 0.01
        SD[2] = 0.01
    
        # Correlation length with DH (0.1 mm higher in bottom)
        lex[1] = lex[0] + 0.13*1e-3
        lex[2] = lex[0] + 0.13*1e-3
        
        # salinity low salinity in top layers
        sal[2] = X[ind_sal] if 'sal' in variables else 13*PSU
        sal[0] = 0.001*PSU
        sal[1] = 0.5*PSU
        
        # density - WS followed by 2 layers NS-DH
        # rho = X[ind_rho] if 'rho' in variables else [350, 270, 270]
        rho[1] = X[ind_rho] if 'rho' in variables else 270
        rho[0] = 350
        rho[2] = rho[1]

    # Temperature of snow surface (this variable is assumed to be known e.g. from ECMWF)
    if avg==True:
        Tsur = X[variables.index('Tsur')] if 'Tsur' in variables else np.mean(self.skt)
    else:
        Tsur = X[variables.index('Tsur')] if 'Tsur' in variables else self.skt[n]

    ## Ice variables
    
    # Sea ice thickness
    if avg==True:
        SIT = X[variables.index('SIT')] if 'SIT' in variables else np.mean(self.sit)
    else:
        SIT = X[variables.index('SIT')] if 'SIT' in variables else self.sit[n]
        
    #print(SD)
    OE.nl_snow=nl_snow
    OE.sd = SD
    OE.sd_tot = np.sum(SD)
    OE.Tsur = Tsur # air temperature
    OE.sit = SIT
    OE.Calc_temp()
    # print('---------------------------------')
    # print(OE.Tsnow)
    # print(OE.T_sii)
    # print('---------------------------------')
    
    
    # print('---------------------------------')   
    # print(SD)
    # print(lex)
    # print(sal)
    # print(rho)
    # print('---------------------------------')
    
    ## SMRT SETUP
    
    ## Make icepack           
    sea_ice_density = 900
    icepack = make_ice_column(ice_type='firstyear',
            thickness=[SIT], 
            temperature=OE.T_sii, 
            microstructure_model='exponential',
            corr_length=0.15e-3,
            brine_inclusion_shape='spheres',
            density=sea_ice_density,
            salinity=6*PSU,
            add_water_substrate="ocean")
    
    ## Make snowpack
    snowpack = make_snowpack(SD, 
                    microstructure_model='exponential',
                    ice_permittivity_model=ssp, 
                    density=rho, 
                    corr_length=lex, 
                    temperature = OE.Tsnow,
                    salinity=sal)
    
    medium = snowpack + icepack
    
    
    m = make_model("iba", "dort")
    
    frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
    
    beamorientation = 55. # [degrees]
    sensor = sensor_list.passive(frequencies, beamorientation)
    
    #% Run simulations
    res = m.run(sensor, medium)
    
    # #% adjust for atmospheric influence

    m = X[variables.index('slope')] if 'slope' in variables else 0.19
    
    if avg==True:
        [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(np.mean(self.tcwv),res,OE.T_sii, m=m)
    else:
        [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(self.tcwv[n],res,OE.T_sii, m=m)
    

    Bias_TBV = np.array([0.546, 1.753, 0.749, -0.020, 0.088])
    Bias_TBH = np.array([1.565, 3.240, 0.454, -0.097, -1.803])
 

    if pol=='V':
        return TB_AMSR_V + Bias_TBV
    elif pol=='VH0':
        y_out = np.array([TB_AMSR_V + Bias_TBV])
        return np.append(y_out, TB_AMSR_H[0]+Bias_TBH[0])
    else:
        return np.concatenate(([TB_AMSR_V + Bias_TBV,TB_AMSR_H + Bias_TBH])) # np.concatenate((res.TbV(),res.TbH()))