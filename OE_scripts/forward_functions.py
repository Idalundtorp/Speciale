# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:37:20 2022

@author: Ida Olsen
"""
import numpy as np
import sys
import re
import matplotlib.pyplot as plt

## append paths
sys.path.append('C:/Users/Ida Olsen/Documents/Speciale_2022/scripts')
parrent_dir = "C:\\Users\\Ida Olsen\\Documents\\Speciale_2022"
sys.path.append(parrent_dir + '\\smrt')
sys.path.append(parrent_dir + '\\Ida_passive_SMRT\\Ida_passive_SMRT\\Code')

## Import SMRT packages
from smrt.permittivity.saline_snow import saline_snow_permittivity_scharien_with_stogryn95 as ssp
# from smrt.permittivity.saline_snow import saline_snow_permittivity_geldsetzer09 as ssp
from smrt import make_ice_column, make_snowpack, sensor_list
from smrt.core.globalconstants import PSU
from smrt import make_model
from smrt import SMRTError
# from smrt.interface.iem_fung92_brogioni10 import IEM_Fung92_Briogoni10

## self made modules
from Functions_temp import Calc_temp_profiles
from roughness import roughness
from Functions_temp import temp_profile
import prepare_inputs_OE as piO
from roughness import Get_slope
## Fill data into class
comb_info = Calc_temp_profiles(plot=False)

sites = ['AK1', 'AK2', 'AK3', 'AK4']
for site in sites:
    for i in range(len(sites)):
        if comb_info[i].site==site:
            s = comb_info[i]

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
    #print(self.sd_tot)
    #print(self.sit)
    self.T_sii = (self.k * self.sd_tot * T_w + ks * self.sit *
                self.Tsur) / (self.k * self.sd_tot + ks * self.sit)
    
    # print(self.T_sii)
   # print(self.Tsur)
     
    from sklearn.linear_model import LinearRegression
    
    x = [0, self.sd_tot]
    y = [self.Tsur, self.T_sii]
    #print(x)
    # print(y)
    model = LinearRegression().fit(np.array(x).reshape(-1,1), np.array(y))
    
    if self.nl_snow == 2:
        # print(self.sd)
        xx = np.array([self.sd[0]/2, self.sd[0] + self.sd[1]/2]).reshape(-1,1)
    else: # snow thinner than 3cm one layer
        xx = np.array([self.sd/2]).reshape(-1,1)
    self.Tsnow = model.predict(xx)
    
    #print(self.Tsnow)
    #print(self.T_sii)
    
    

    if plot==True:
        # ice/snow thickness arrays for showing temperature profiles    
        ice_thick = np.array([0, -self.sit])
        t_ice = np.array([self.T_sii, T_w])
        snow_thick = [np.sum(self.sd), self.sd[0] + self.sd[1]/2, self.sd[0]/2, 0]
        t_snow = [self.Tsur] + [Ts for Ts in self.Tsnow] + [self.T_sii]
        
        # print(t_snow)
        #print(snow_thick)
        # print(t_ice)
        
        plt.figure()
        plt.plot(t_snow, snow_thick,'.-')
        plt.plot(t_ice, ice_thick,'.-')
        # plt.ylim([-0.1,0.1])
        plt.grid()
        plt.xlabel('Temperature [K]')
        plt.ylabel('Depth [m]')

def forward_mean(X, variables=[], sitenum=[]):
    """
    

    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    variables : list
        A list of variables which are being optimized
    sitenum : integer
        integer of the site number to use
    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # print(variables)
    # print(sitenum)
    # print(X)
    s = comb_info[sitenum]
    SD = X[variables.index('SD')] if 'SD' in variables else np.mean(s.SD_tot)
    SIT = X[variables.index('SIT')] if 'SIT' in variables else s.SIT
    lex = X[variables.index('lex')] if 'lex' in variables else np.mean(s.corr_m)
    sal = X[variables.index('sal')] if 'sal' in variables else np.mean(s.sal_m)
    rho = X[variables.index('rho')] if 'rho' in variables else np.mean(s.rho_m)

    
    ## temperatures calculated
    temp_ice = np.mean(np.array(s.Tcal_ice, dtype=object).flatten()) ## Mean snow/ice interface temperature
    temp_snow = np.mean(np.concatenate((np.array(s.Tcal_snow, dtype=object))))

    ## Make icepack        
    sea_ice_density = 900
    icepack = make_ice_column(ice_type='firstyear',
            thickness=[SIT], temperature=temp_ice, 
            microstructure_model='independent_sphere',
            radius=1e-3,
            brine_inclusion_shape='spheres',
            density=sea_ice_density,
            salinity=7*PSU,
            add_water_substrate=True)

    ## Make snowpack
    snowpack = make_snowpack([SD], microstructure_model='exponential',
                    ice_permittivity_model=ssp, density=rho, 
                    corr_length=lex, temperature = temp_snow,
                    salinity=sal)

    medium = snowpack + icepack

    
    m = make_model("iba", "dort")

    frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]

    beamorientation = 55. # [degrees]
    sensor = sensor_list.passive(frequencies, beamorientation)

    #% Run simulations
    res = m.run(sensor, medium)
    
    #% adjust for atmospheric influence
    tcwv = np.load(s.site + '_ECMWF_tcwv.npy', allow_pickle=True)
    [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(tcwv,res,temp_ice,site=s.site, several=False)
    return np.concatenate((TB_AMSR_V,TB_AMSR_H))


def forward(X, variables, nl_snow, sitenum=1):
    """
    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    variables : list
        A list of variables which are being optimized
    sitenum : integer
        integer of the site number to use
    Returns
    -------
    TYPE
        DESCRIPTION.
    
    """
    # print(X[0])
    variables = [re.sub(r'\d+','',v) for v in variables]
    
    #print(nl_snow)
    ## define inputs
    s = comb_info[sitenum]
    ## total snowdepth is split into nl_snow layers
    ind_SD = [i for i, j in enumerate(variables) if j == 'SD']
    SD = X[ind_SD] if 'SD' in variables else np.ones(nl_snow)*np.mean(s.SD_tot)/nl_snow
    ## salinity is high at the bottom of the snowpack
    ind_sal = [i for i, j in enumerate(variables) if j == 'sal']
    sal = X[ind_sal] if 'sal' in variables else np.linspace(np.min(np.concatenate((s.sal))),
                                                                            np.max(np.concatenate((s.sal))),
                                                                            nl_snow)
    ## density is close to constant, correlation length also show no depth dependence
    ind_rho = [i for i, j in enumerate(variables) if j == 'rho']
    rho = X[ind_rho] if 'rho' in variables else np.ones(nl_snow)*np.mean(s.rho_m)
    
    ind_lex = [i for i, j in enumerate(variables) if j == 'lex']
    lex = X[ind_lex] if 'lex' in variables else np.ones(nl_snow)*np.mean(s.corr_m)
    
    ## temperatures calculated
    SIT = X[variables.index('SIT')] if 'SIT' in variables else s.SIT_m
    [T_sii,temp_profile_snow,temp_profile_ice] = temp_profile(2.10, np.mean(s.SD_tot), s.SIT_m, s.Tsur, nl_snow)
    
    
    ## Make icepack        
    sea_ice_density = 900
    icepack = make_ice_column(ice_type='firstyear',
            thickness=[SIT], temperature=T_sii, 
            microstructure_model='independent_sphere',
            radius=1e-3,
            brine_inclusion_shape='spheres',
            density=sea_ice_density,
            salinity=7*PSU,
            add_water_substrate=True)
    
    # print('SD is:', SD)
    # print('lex is:', lex)
    # print('sal is:', sal)

    ## Make snowpack
    snowpack = make_snowpack(SD, microstructure_model='exponential',
                    ice_permittivity_model=ssp, density=rho, 
                    corr_length=lex, temperature = temp_profile_snow,
                    salinity=sal)
    
    medium = snowpack + icepack
    
    
    m = make_model("iba", "dort")
    
    frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
    
    beamorientation = 55. # [degrees]
    sensor = sensor_list.passive(frequencies, beamorientation)
    
    #% Run simulations
    res = m.run(sensor, medium)
    
    #% adjust for atmospheric influence
    tcwv = np.load(s.site + '_ECMWF_tcwv.npy', allow_pickle=True)
    
    m = X[variables.index('slope')] if 'slope' in variables else Get_slope(s.site)
    [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(tcwv,res,T_sii, m=m, site=s.site, several=False)
    return np.concatenate((TB_AMSR_V,TB_AMSR_H))

from prepare_inputs_OE import OIB_data
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
     
def forward_OIB(X, self, variables, nl_snow, n, nls, pol, logarithmic=True, avg=False, sim=False):

    
    # X = np.array([0.20, 0.02, 250, 260])
    # variables = ['SD','SD1', 'Tsnow', 'Tsnow']
    # nl_snow=2
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
    # print(X)
    # print('------------')
    if np.isnan(X).any():
        print(X)
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
    
    if np.isnan(np.array(X[ind_lex])):
        X[ind_lex] = np.log10(np.mean(FYI_data.sd_lex, axis=0))[0]
        print('nan')
        
    # print('rho index:', ind_rho)
    # rho = X[ind_rho] if 'rho' in variables else np.mean(FYI_data.sd_rho, axis=0)
    rho = X[ind_rho] if 'rho' in variables else [335, 270]
    if nls==2:
        SD = X[ind_SD] if 'SD' in variables else np.mean(FYI_data.sd, axis=0)
        #lex = X[ind_lex] if 'lex' in variables else np.mean(FYI_data.sd_lex, axis=0)
        #sal = X[ind_sal] if 'sal' in variables else np.mean(FYI_data.sd_sal, axis=0)
        
        lex = X[ind_lex] if 'lex' in variables else [0.223e-3, 0.204e-3]
        sal = X[ind_sal] if 'sal' in variables else [0.01*PSU, 4*PSU]
        
        
        if logarithmic==True:
            # if n<1:
            #     print('Assuming logarithmic input!')
            
            SD = 10**X[ind_SD] if 'SD' in variables else np.mean(FYI_data.sd, axis=0)
            lex = 10**X[ind_lex] if 'lex' in variables else np.mean(FYI_data.sd_lex, axis=0)
            sal = 10**X[ind_sal] if 'sal' in variables else np.mean(FYI_data.sd_sal, axis=0)

    elif nls==1:
        rho = np.mean(FYI_data.sd_rho, axis=0)
        # print(rho)
        rho[1] = X[ind_rho] if 'rho' in variables else rho[1]
        SD  = np.mean(FYI_data.sd, axis=0)
        lex = np.mean(FYI_data.sd_lex, axis=0)
        sal = np.mean(FYI_data.sd_sal, axis=0)

        
        # print(lex)
        # print(X[ind_lex])

        
        ## convert OE value to linear units
        if logarithmic==True:
            # if n<1:
            #     print('Assuming logarithmic input!')
            SD[0] = 10**X[ind_SD] if 'SD' in variables else SD[0]
            lex[0] = 10**X[ind_lex] if 'lex' in variables else lex[0]
            sal[1] = 10**X[ind_sal] if 'sal' in variables else sal[1]
        # elif avg==False:
        #     print('using real values for unknowns')
        #     rho = X[ind_rho] if 'rho' in variables else FYI_data.sd_rho[n]
        #     lex = X[ind_lex] if 'lex' in variables else FYI_data.sd_lex[n]
        #     sal = X[ind_sal] if 'sal' in variables else FYI_data.sd_sal[n]
        #     SD[0]  = X[ind_SD] if 'SD' in variables else SD[0]
            # print(lex)
        else:
            lex[0] = X[ind_lex] if 'lex' in variables else lex[0]
            sal[1] = X[ind_sal] if 'sal' in variables else sal[1]           
            SD[0]  = X[ind_SD] if 'SD' in variables else SD[0]


    #print('correlation length is:', lex*1e3)
    #print('salinity is:', sal*1e3)
        
        # print(lex)
        # print(sal)
        #print(SD)

    
        
    # Temperature of snow surface
    Tsur = X[variables.index('Tsur')] if 'Tsur' in variables else FYI_data.skt[n]
    # print(Tsur)
    if avg==True:
        Tsur = X[variables.index('Tsur')] if 'Tsur' in variables else np.mean(FYI_data.skt)
        Tsur = Tsur
        # print(Tsur)
        # if avg==True:
        #     print('T surface is:', Tsur)
        # Tsur=np.mean(FYI_data.skt)
        # print(Tsur)
    # Tsur =  np.mean(FYI_data.skt[n])
    ## Ice variables
    # Sea ice thickness
    SIT = X[variables.index('SIT')] if 'SIT' in variables else FYI_data.sit[n] # np.mean(FYI_data.sit)
    # if logarithmic==True:
    #     SIT = 10**X[variables.index('SIT')] if 'SIT' in variables else np.mean(FYI_data.sit, axis=0)
    
    # Temperature of snow/ice interface
    Tice = X[variables.index('Tice')] if 'Tice' in variables else  np.mean(FYI_data.Tice, axis=0)
    
    # Used to calculate the temperature profile
    # print(SD)
    # if SD[0] < 0:
    #     #print(SD)
    #     SD[0] = abs(SD[0])
    # elif SD[1] <0:
    #     SD[1] = abs(SD[1])
        
    #print(SD)
    OE.sd = SD
    OE.sd_tot = np.sum(SD)
    OE.Tsur = Tsur # air temperature
    OE.sit = SIT
    OE.Calc_temp()
    # print('---------------------------------')
    # print(SIT)
    # print(OE.T_sii)
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
    # print(SD)
    # print(rho)
    # print(lex)
    # print(sal)
    # print(Tsur)
    snowpack = make_snowpack(SD, 
                    microstructure_model='exponential',
                    ice_permittivity_model=ssp, 
                    density=rho, 
                    corr_length=lex, 
                    temperature = OE.Tsnow,
                    salinity=sal)
    # print(SD)
    # print(rho)
    # print(lex)
    # print(sal)
    # print(OE.Tsnow)
    
    medium = snowpack + icepack
    
    
    m = make_model("iba", "dort")
    
    frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
    
    beamorientation = 55. # [degrees]
    sensor = sensor_list.passive(frequencies, beamorientation)
    
    #% Run simulations
    res = m.run(sensor, medium)
    
    # #% adjust for atmospheric influence

    m = X[variables.index('slope')] if 'slope' in variables else 0.08 # 0.19 # 0.08
    # print(m)
    
    if avg==False:
        [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(FYI_data.tcwv[n],res,OE.T_sii, m=m)
    else:
        [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(np.mean(FYI_data.tcwv),res,OE.T_sii, m=m)
    
    if sim==False:
        # print('right bias')
        #Bias_TBV = [ 0.43621623,1.67760632,0.63745003,-0.46549513,-1.40957493]
        #Bias_TBH = [ 2.22099808, 4.14636338, 1.67438677, 0.88539753, -1.76596956]
        # Bias_TBV = [1.07411439, 2.31576401, 1.15553819, -0.17905842, -3.30731626]
        #Bias_TBV = [0, 1.124714, -0.346821, -1.997182, -6.157140]
        Bias_TBV =[1.08011672,2.24896883,1.00845678,-0.34775236,-2.6658168]
        Bias_TBH =[3.64055518, 5.14219181, 2.0641708, 0.96405109, -2.87242774]

        #Bias_TBV = [ 0, 0, 0, 0, 0]
        #Bias_TBH = [ 0, 0, 0, 0, 0]
        #Bias_TBV = [ -2.08, -0.46, -0.41, -0.55, 0.3]
        #Bias_TBH = [ 0.22, 0.65, -2.48, -2.6, -3.32]
    
    ## Test bias!!!
    if avg==True and sim==True:
        # print('yes')
        #Bias_TBV = [-0.58746677,-0.53409001,-0.45314495,-0.46169971,-0.44553977]
        #Bias_TBH = [-0.82481766,-0.99011067,-1.19102826,-1.24453696,-1.19335456]
        Bias_TBV = [ 0, 0, 0, 0, 0]
        Bias_TBH = [ 0, 0, 0, 0, 0]
        
        

    elif sim==True and avg==False:
        Bias_TBV = [ 0, 0, 0, 0, 0]
        Bias_TBH = [ 0, 0, 0, 0, 0]   
    # print(res.TbV())
    if pol=='V':
        return TB_AMSR_V + Bias_TBV
    else:
        return np.concatenate(([TB_AMSR_V + Bias_TBV,TB_AMSR_H + Bias_TBH])) # np.concatenate((res.TbV(),res.TbH()))


def forward_lex_sal(X, variables, nl_snow, n, nls, avg=False, sim=False):
    
    # X = np.array([0.20, 0.02, 250, 260])
    # variables = ['SD','SD1', 'Tsnow', 'Tsnow']
    # nl_snow=2
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
    if np.isnan(X).any():
        print(X)
    ## OIB data - prepared inputs for OE (based on average and std of OIB data)
    FYI_data = piO.OIB_data('FYI')
    # OE class to compute temperature profile
    OE = OE_data('FYI')
    
    # string of variables
    variables = [re.sub(r'\d+','',v) for v in variables]
    
    
    ## Snow variables
    ind_sal = [i for i, j in enumerate(variables) if j == 'sal']
    ind_lex = [i for i, j in enumerate(variables) if j == 'lex']
    
    if nls==2:
        lex = X[ind_lex] if 'lex' in variables else np.mean(FYI_data.sd_lex, axis=0)
        sal = X[ind_sal] if 'sal' in variables else np.mean(FYI_data.sd_sal, axis=0)
        

    elif nls==1:
        lex = np.mean(FYI_data.sd_lex, axis=0)
        sal = np.mean(FYI_data.sd_sal, axis=0)
        
        lex[0] = X[ind_lex] if 'lex' in variables else lex[0]
        sal[0] = X[ind_sal] if 'sal' in variables else sal[0]

    
        
    # Temperature of snow surface
    Tsur = np.mean(FYI_data.skt[n])
    SIT = FYI_data.sit[n]
    SD = FYI_data.sd[n]
    rho = np.mean(FYI_data.sd_rho)

    OE.sd = SD
    OE.sd_tot = np.sum(SD)
    OE.Tsur = Tsur # air temperature
    OE.sit = SIT
    OE.Calc_temp()
 
    ## Make icepack           
    sea_ice_density = 900
    icepack = make_ice_column(ice_type='firstyear',
            thickness=[SIT], 
            temperature=OE.T_sii, 
            microstructure_model='exponential',
            corr_length=0.3e-3,
            brine_inclusion_shape='spheres',
            density=sea_ice_density,
            salinity=7*PSU,
            add_water_substrate="ocean")
    
    # if we have negative salinity: 
    if (np.array(sal)<0).any():
        print(sal)
    
        ## Make snowpack
        snowpack = make_snowpack(SD, 
                        microstructure_model='exponential',
                        ice_permittivity_model=ssp,
                        density=rho, 
                        corr_length=lex, 
                        temperature = OE.Tsnow,
                        salinity=sal)
    else:
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
    m = X[variables.index('slope')] if 'slope' in variables else 0.19 # 0.08
    [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(np.mean(FYI_data.tcwv),res,OE.T_sii, m=m)
    
    if sim==False:
        # print('right bias')
        Bias_TBV = [ -2.08, -0.46, -0.41, -0.55, 0.3]
        Bias_TBH = [ 0.22, 0.65, -2.48, -2.6, -3.32] 
      
    return np.concatenate(([TB_AMSR_V + Bias_TBV,TB_AMSR_H + Bias_TBH])) # np.concatenate((res.TbV(),res.TbH()))
