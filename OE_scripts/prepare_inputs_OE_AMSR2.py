# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 14:45:20 2023

@author: Ida Olsen
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from smrt.core.globalconstants import PSU

## use data from CB to get varians of rho and lex
from Functions_temp import fill_data

#%% OIB data
class OIB_data():
    def __init__(self, icetype, new=False):
        if icetype=='FYI':
            datapath = "C:/Users/Ida Olsen/Documents/Speciale_2022/Optimal_Estimation/OE_ida/Simulated_data/TB_simulated.npy"
            if new==True:
                datapath = "C:/Users/Ida Olsen/Documents/Speciale_2022/Optimal_Estimation/OE_ida/Simulated_data/AMSR2_TB_simulated.npy"
        # plot path
        self.pp = "C:/Users/Ida Olsen/Documents/Speciale_2022/figures"
        
        #if new==False: 
        with open(datapath, 'rb') as f:
            self.TbV_sim = np.load(f)
            self.TbH_sim = np.load(f)
            self.sd = np.load(f)
            self.sit = np.load(f)
            self.sd_rho = np.load(f)
            self.sd_lex = np.load(f)
            self.sd_sal = np.load(f)
            self.Tsnow = np.load(f)
            self.t2m = np.load(f)
            self.skt = np.load(f)
            self.tcwv = np.load(f)
            self.Tice = np.load(f)
            self.OIB_surf_rough = np.load(f)
            self.TbV_AMSR_all = np.load(f)
            self.TbH_AMSR_all = np.load(f)
            self.lat = np.load(f)
            self.lon = np.load(f)  
            # print(self.lon)
            self.date = np.load(f, allow_pickle=True)



#%%
def prepare_inputs_FYI_AMSR2(self, variables, logarithmic, n=-1, nl_snow=3):
    """
    Parameters
    ----------
    Nd : TYPE
        Number of datapoints to use
    variables : TYPE
        String of variables 
    sim: BOOL
        True is simulated data
        False is real data
    
    Returns
    -------
    None.
    
    """
    
    ## First guesses

    
    ## construction 2-layer snowpack with salt in bottom of snowpacl
    # snow

    if nl_snow==3:
        
        ## START GUESSES
    

        # comes in in log units
        SD_ap = np.mean(self.sd)
        lex_ap = 0.24*1e-3
        sal_ap = 13*PSU
        rho_ap = 270
        ## other start guesses 
        Tsur_ap = np.mean(self.skt, axis=0)       # K
        # ice
        SIT_ap = np.mean(self.sit, axis=0)         # m

    elif nl_snow==2:
        SD_ap = np.mean(self.sd)
        lex_ap = 0.25*1e-3
        # sal_ap = 4*PSU
        sal_ap = 6*PSU
        rho_ap = 270
        
        ## other start guesses 
        Tsur_ap = np.mean(self.skt, axis=0)       # K
        # ice
        SIT_ap = np.mean(self.sit, axis=0)         # m
        # # Make 3 snow layers with 25% DH
        # SD_tot  = np.mean(np.sum(self.sd))
        # sd0 = SD_tot-SD_tot*0.25
        # sd1 = SD_tot*0.25
        # SD_ap = [sd0, sd1]
        # # winslab layer with two DH layers below
        # lex_ap = [0.23*1e-3, 0.33*1e-3]
        # sal_ap = [0.001*PSU, 5*PSU]
        # rho_ap = [350, 270]

    if n!=-1:
        ## SD initial guess based on GR  
        ## make SD guess from cavalleri based on GR
        # TB6V = self.TbV_AMSR_all[n][0]
        # TB18V = self.TbV_AMSR_all[n][2]
        # TB36V = self.TbV_AMSR_all[n][4]
        # GR_69_187_V  = (TB18V-TB6V)/(TB18V+TB6V)
        # GR_187_37_V  = (TB36V-TB18V)/(TB36V+TB18V)
        # print(GR_69_187_V)
        # # SD_ap = 0.1926 - 5.53*GR_69_187_V
        # # print(SD_ap)
        # SD_ap = 0.029 - 7.824*GR_187_37_V
        # print(SD_ap)
        # print('gradient ratio used')
        # SD_ap = np.log10(SD_ap)
        # Tsur_ap = self.skt[n]
        # SIT_ap = self.sit[n]
        print('know value guesses - not used')
    # print(SD_ap)  
    ## VARIANSES

    # varians from CB
    snowpit_data = fill_data()
    rho_CB = np.concatenate(np.concatenate([pit.rho for pit in snowpit_data]))
    lex_CB = np.concatenate(np.concatenate([pit.corr for pit in snowpit_data]))
    sal_CB = np.concatenate(np.concatenate([pit.sal for pit in snowpit_data]))
    
    if logarithmic==True:
        lex_CB = np.log10(lex_CB)
        # sal_CB = np.log10(sal_CB)
        lex_ap = np.log10(0.24*1e-3)
        # sal_ap = np.log10(13*PSU)
    
    var_lex = np.std(lex_CB)**2
    # var_sal = np.std(sal_CB)**2
    var_rho = np.std(rho_CB)**2
    var_SD = np.std(self.sd)**2
    # print(var_SD)
    var_sal = np.std(sal_CB)**2 # (2*PSU)**2
    
    var_SIT = np.std(self.sit)**2
    var_Tsur = np.std(self.skt)**2
    var_slope = (0.03)**2
    

    m_ap = 0.19

    var = variables # [var for var in variables
    varnames = variables
    
    # start guesses
    x_ap = np.zeros(len(var)) # starting guess array
    variances = np.zeros(len(var)) # variances of variables
    varian = []
    
    # i for s in used to flatten the list
    x_lowerlimit = dict.fromkeys(varnames)
    x_upperlimit = dict.fromkeys(varnames)
    x_info = {}
    
    ## SNOW VARIABLES
    if nl_snow>1:
        
        # Snow depth
        if 'SD' in variables:
            ind = [i for i, j in enumerate(var) if j == 'SD']
            x_ap[ind] = SD_ap
            variances[ind] = var_SD
            for layer in list(range(nl_snow)):
                if logarithmic==True:
                    x_lowerlimit['SD'] = np.log10(0.0)
                    x_upperlimit['SD'] = np.log10(1.0)
                else:
                    x_lowerlimit['SD'] = 0.0
                    x_upperlimit['SD'] = 1.0                 
                #x_info['SD' + str(layer) + ' cm'] = np.round(SD_ap[layer] * 1e2, 2)
        
        # correlation length of snow
        if 'lex' in variables:
            ind = [i for i, j in enumerate(var) if j == 'lex']
            x_ap[ind] = lex_ap
            variances[ind] = var_lex
            for layer in list(range(nl_snow)):
                if logarithmic==True:
                    x_lowerlimit['lex'] = np.log10(0.01e-3 )
                    x_upperlimit['lex'] = np.log10(1.0e-3)
                else:
                    x_lowerlimit['lex'] = 0.01e-3
                    x_upperlimit['lex'] = 1.0e-3
                #x_info['lex' + str(layer) + ' mm'] = np.round(lex_ap[layer] * 1e3, 2)
    
        # salinity of snow
        if 'sal' in variables:
            ind = [i for i, j in enumerate(var) if j == 'sal']
            x_ap[ind] = sal_ap
            variances[ind] = var_sal
            for layer in list(range(nl_snow)):
                #if logarithmic==True:
                #    x_lowerlimit['sal'] = np.log10(0.001 * PSU)
                #    x_upperlimit['sal' ] = np.log10(20.0 * PSU)
                #else:
                x_lowerlimit['sal'] = 0.001 * PSU
                x_upperlimit['sal' ] = 30.0 * PSU
                #x_info['sal' + str(layer) + ' PSU'] = np.round(sal_ap[layer]*1e3, 2)
        
        # Density of snow
        if 'rho' in variables:
            ind = [i for i, j in enumerate(var) if j == 'rho']
            x_ap[ind] = rho_ap
            variances[ind] = var_rho
            for layer in list(range(nl_snow)):
                x_lowerlimit['rho'] = 150 # kg/m3
                x_upperlimit['rho'] = 450 # kg/m3
                #x_info['sal' + str(layer) + ' kg/m3'] = np.round(rho_ap[layer], 2)
    

    if 'Tsur' in variables:
        ind = [i for i, j in enumerate(var) if j == 'Tsur']
        x_ap[ind] = Tsur_ap
        variances[ind] = var_Tsur
        x_lowerlimit['Tsur'] = 273.15 - 30  # K
        x_upperlimit['Tsur'] = 273.15 - 5  # K
        x_info['Tsur K'] = np.round(Tsur_ap, 2)
        
    ## ICE VARIABLES
    if 'SIT' in variables:
        ind = [i for i, j in enumerate(var) if j == 'SIT']
        x_ap[ind] = SIT_ap
        variances[ind] = var_SIT
        x_lowerlimit['SIT'] = 0.0
        x_upperlimit['SIT'] = 5.0
        x_info['SIT cm'] = SIT_ap * 1e2
    
    # roughness of snow/ice interface 
    if 'slope' in variables:
        ind = [i for i, j in enumerate(var) if j == 'slope']
        x_ap[ind] = m_ap
        variances[ind] = var_slope
        x_lowerlimit['slope'] = 0.01     # m/m
        x_upperlimit['slope'] = 0.4   # m/m  
        x_info['slope'] = 0.19 # m_ap

    # if logarithmic==True:
    #     for key, count in zip(x_lowerlimit.keys(), range(len(x_lowerlimit.keys()))):
                
    #         if 'SD' in key or 'lex' in key  or'sal' in key:
    #             print('Convertin '+key+' variances to log10')
    #             x_lowerlimit[key]=np.log10(x_lowerlimit[key])
    #             x_upperlimit[key]=np.log10(x_upperlimit[key])
    #             #print(varian[count])
    #             varian[count] = np.log10(varian[count])
    #             #print(10**varian[count])
                
    ## fill covariance matrix
    # x_cov = np.cov(np.array(varian))
    x_cov = np.diag(variances)

    return [x_ap, x_cov, x_lowerlimit, x_upperlimit, x_info]

def prepare_inputs_FYI_AMSR2_new(self, variables, logarithmic, n=-1, nl_snow=3):
    """
    Parameters
    ----------
    Nd : TYPE
        Number of datapoints to use
    variables : TYPE
        String of variables 
    sim: BOOL
        True is simulated data
        False is real data
    
    Returns
    -------
    None.
    
    """
    
    ## First guesses

    
    ## construction 2-layer snowpack with salt in bottom of snowpacl
    # snow

    if nl_snow==3:
        
        ## START GUESSES
    

        # comes in in log units
        SD_ap = np.mean(self.sd)
        lex_ap = 0.24*1e-3
        sal_ap = 13*PSU
        rho_ap = 270
 
    ## VARIANSES

    # varians from CB
    snowpit_data = fill_data()
    rho_CB = np.concatenate(np.concatenate([pit.rho for pit in snowpit_data]))
    lex_CB = np.concatenate(np.concatenate([pit.corr for pit in snowpit_data]))
    sal_CB = np.concatenate(np.concatenate([pit.sal for pit in snowpit_data]))
    
    # if logarithmic==True:
    #     lex_CB = np.log10(lex_CB)
    #     # sal_CB = np.log10(sal_CB)
    #     lex_ap = np.log10(0.24*1e-3)
    #     # sal_ap = np.log10(13*PSU)
    
    var_lex = np.std(lex_CB)**2
    var_rho = np.std(rho_CB)**2
    var_SD = np.std(self.sd)**2
    var_sal = (2*PSU)**2
    
    var_SIT = np.std(self.sit)**2
    var_Tsur = np.std(self.skt)**2
    var_slope = (0.03)**2
    
    ## other start guesses 
    Tsur_ap = np.mean(self.skt)       # K
    # ice
    SIT_ap = np.mean(self.sit)         # m
    m_ap = 0.19

    var = variables # [var for var in variables
    varnames = variables
    
    # start guesses
    x_ap = np.zeros(len(var)) # starting guess array
    variances = np.zeros(len(var)) # variances of variables
    varian = []
    
    # i for s in used to flatten the list
    x_lowerlimit = dict.fromkeys(varnames)
    x_upperlimit = dict.fromkeys(varnames)
    x_info = {}
    
    ## SNOW VARIABLES
    if nl_snow>1:
        
        # Snow depth
        if 'SD' in variables:
            ind = [i for i, j in enumerate(var) if j == 'SD']
            x_ap[ind] = SD_ap
            variances[ind] = var_SD
            for layer in list(range(nl_snow)):
                if logarithmic==True:
                    x_lowerlimit['SD'] = np.log10(0.0)
                    x_upperlimit['SD'] = np.log10(1.0)
                else:
                    x_lowerlimit['SD'] = 0.0
                    x_upperlimit['SD'] = 1.0                 
                #x_info['SD' + str(layer) + ' cm'] = np.round(SD_ap[layer] * 1e2, 2)
        
        # correlation length of snow
        if 'lex' in variables:
            ind = [i for i, j in enumerate(var) if j == 'lex']
            x_ap[ind] = lex_ap
            variances[ind] = var_lex
            for layer in list(range(nl_snow)):
                if logarithmic==True:
                    x_lowerlimit['lex'] = np.log10(0.01e-3 )
                    x_upperlimit['lex'] = np.log10(1.0e-3)
                else:
                    x_lowerlimit['lex'] = 0.01e-3
                    x_upperlimit['lex'] = 1.0e-3
                #x_info['lex' + str(layer) + ' mm'] = np.round(lex_ap[layer] * 1e3, 2)
    
        # salinity of snow
        if 'sal' in variables:
            ind = [i for i, j in enumerate(var) if j == 'sal']
            x_ap[ind] = sal_ap
            variances[ind] = var_sal
            for layer in list(range(nl_snow)):
                #if logarithmic==True:
                #    x_lowerlimit['sal'] = np.log10(0.001 * PSU)
                #    x_upperlimit['sal' ] = np.log10(20.0 * PSU)
                #else:
                x_lowerlimit['sal'] = 0.001 * PSU
                x_upperlimit['sal' ] = 30.0 * PSU
                #x_info['sal' + str(layer) + ' PSU'] = np.round(sal_ap[layer]*1e3, 2)
        
        # Density of snow
        if 'rho' in variables:
            ind = [i for i, j in enumerate(var) if j == 'rho']
            x_ap[ind] = rho_ap
            variances[ind] = var_rho
            for layer in list(range(nl_snow)):
                x_lowerlimit['rho'] = 150 # kg/m3
                x_upperlimit['rho'] = 450 # kg/m3
                #x_info['sal' + str(layer) + ' kg/m3'] = np.round(rho_ap[layer], 2)
    

    if 'Tsur' in variables:
        ind = [i for i, j in enumerate(var) if j == 'Tsur']
        x_ap[ind] = Tsur_ap
        variances[ind] = var_Tsur
        x_lowerlimit['Tsur'] = 273.15 - 50  # K
        x_upperlimit['Tsur'] = 273.15 - 2  # K
        x_info['Tsur K'] = np.round(Tsur_ap, 2)
        
    ## ICE VARIABLES
    if 'SIT' in variables:
        ind = [i for i, j in enumerate(var) if j == 'SIT']
        x_ap[ind] = SIT_ap
        variances[ind] = var_SIT
        x_lowerlimit['SIT'] = 0.0
        x_upperlimit['SIT'] = 5.0
        x_info['SIT cm'] = SIT_ap * 1e2
    
    # roughness of snow/ice interface 
    if 'slope' in variables:
        ind = [i for i, j in enumerate(var) if j == 'slope']
        x_ap[ind] = m_ap
        variances[ind] = var_slope
        x_lowerlimit['slope'] = 0.001     # m/m
        x_upperlimit['slope'] = 0.45   # m/m  
        x_info['slope'] = 0.19 # m_ap

    x_cov = np.diag(variances)

    return [x_ap, x_cov, x_lowerlimit, x_upperlimit, x_info]

