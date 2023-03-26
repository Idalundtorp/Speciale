# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 13:43:49 2022

@author: Ida Olsen
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from smrt.core.globalconstants import PSU

from roughness import Get_slope

# site, variance=False, plot=False

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

            
    def plot(self, plot=True, plot_idv=False):
        if plot==True:
            plt.figure(figsize=(8,6))
            frequencies=[6.925e9,10.65e9,18.7e9,23.8e9,36.5e9];
            freqs = [f*1e-9 for f in frequencies]
            y = np.mean(self.TbV_sim, axis=0)
            err = np.std(self.TbV_sim, axis=0)
            plt.plot(freqs, y, c='r',label = 'TBV')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='r', facecolor='r', alpha=0.3)
            
            y = np.mean(self.TbH_sim, axis=0)
            err = np.std(self.TbH_sim, axis=0)
            plt.plot(freqs,y, c='b', label = 'TBH')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='b', facecolor='b', alpha=0.3)
    
            
            y = np.mean(self.TbV_AMSR_all, axis=0)
            err = np.std(self.TbV_AMSR_all, axis=0)
            plt.plot(freqs, y, c='k', label = 'TBV AMSR')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='k', facecolor='k', alpha=0.3)
            
            y = np.mean(self.TbH_AMSR_all, axis=0)
            err = np.std(self.TbH_AMSR_all, axis=0)
            
            plt.plot(freqs, y, c='green', label = 'TBH AMSR')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='g', facecolor='g', alpha=0.3)
            plt.grid()
            plt.ylim([200,260])
            # plt.legend(['TBV', 'TBH', 'TBV AMSR', 'TBH AMSR'])
            plt.legend()
            plt.xlabel('Freq GHz')
            plt.ylabel('TB [K]')
            plt.title('Simulated FYI TBs from OIB data 2013-2019')
            plt.savefig(os.path.join(self.pp, 'OE_figures/Simulated_data.png'), bbox_inches='tight')
            plt.show()
            
            if plot_idv==True:
                plt.figure(figsize=(8,6))
                for TBV, TBH in zip(self.TbV_sim, self.TbH_sim):
                    plt.scatter(freqs, TBV, c='r')
                    plt.scatter(freqs, TBH, c='b')
                plt.legend(['TBV', 'TBH'])
                plt.grid()
                plt.ylim([220,270])
                plt.xlabel('Freq GHz')
                plt.ylabel('TB [K]')
                plt.title('Simulated FYI TBs from OIB data 2013-2019, nr data: ' + str(self.sd.shape[0]))
                plt.savefig(os.path.join(self.pp, 'OE_figures/Simulated_data.png'), bbox_inches='tight')
                plt.show()


def prepare_inputs_FYI(self, variables, logarithmic, Nd=100, nl_snow=2, sim=False, nls=2):
    # nl_snow = 2
    # Nd = 100
    # variables= ["SD", "lex", "sal", "Tsur", "SIT"]
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
    #s = OIB_data('FYI')
    #self.plot(plot=False, plot_idv=False)
    
    ## First guesses
    # snow
    SD_ap  = np.mean(self.sd, axis=0)         # m
    lex_ap = np.mean(self.sd_lex, axis=0)     # m
    sal_ap = np.mean(self.sd_sal, axis=0)
    rho_ap = np.mean(self.sd_rho, axis=0)
    ## make log values
    if logarithmic==True:
        print('Logarithmic convertion made')
        SD_ap  = np.log10(SD_ap)
        sal_ap = np.log10(sal_ap)
        lex_ap = np.log10(lex_ap)
    
    rho_ap = np.mean(self.sd_rho, axis=0)     # kg/m^3
    Tsur_ap = np.mean(self.skt, axis=0)       # K
    # ice
    SIT_ap = np.mean(self.sit, axis=0)         # m
    Tice_ap = np.mean(self.Tice, axis=0)      # K
    m = []
    for site in ['AK1', 'AK2', 'AK3', 'AK4']:
        m.append(Get_slope(site, variance=False)[0])
    m_ap = np.mean(m)
    
    # SD_ap = [0.10, 0.03]
    # lex_ap = [0.223e-3, 0.204e-3]
    # sal_ap = [0.01*PSU, 4*PSU]
    # rho_ap = [335, 270]
    # m_ap = 0.19

    
    if nls==2:
        var = [i for s in [nl_snow*[var] if var != 'slope' and var != 'Tsur' and var != 'SIT' and var !='Tice' else [var] for var in variables ] for i in s]
        varnames = [v + str(i) if v != 'slope' and v != 'Tsur' and v != 'SIT' and v !='Tice' else v for v,i in zip(var, list(range(nl_snow))*len(var))]
    elif nls==1:
        var = np.array([[var] for var in variables]).flatten()
        varnames = var

        ## make array with correct number of snow elements
        
    # var = [i for s in [nl_snow*[var] if  var != 'lex' and var != 'sal' and var != 'slope' and var != 'Tsur' and var != 'SIT' and var !='Tice' else [var] for var in variables ] for i in s]
    ## make array of unique varnames
    # varnames = [v + str(i) if v != 'lex' and v != 'sal' and v != 'slope' and v != 'Tsur' and v != 'SIT' and v !='Tice' else v for v,i in zip(var, list(range(nl_snow))*len(var))]
    x_ap = np.zeros(len(var)) # starting guess array
    variances = np.zeros(len(var)) # variances of variables
    varian = []
    
    # i for s in used to flatten the list
    x_lowerlimit = dict.fromkeys(varnames)
    x_upperlimit = dict.fromkeys(varnames)
    x_info = {}
    
    ## SNOW VARIABLES
    if nls==2:
        
        # Snow depth
        if 'SD' in variables:
            varian.append((np.sum(self.sd, axis=1)))
            varian.append((np.sum(self.sd, axis=1)))
            
            ind = [i for i, j in enumerate(var) if j == 'SD']
            x_ap[ind] = SD_ap
            # variances[ind] = SD_var
            for layer in list(range(nl_snow)):
                x_lowerlimit['SD' + str(layer)] = 0.0
                x_upperlimit['SD' + str(layer)] = 1.0
                x_info['SD' + str(layer) + ' cm'] = np.round(SD_ap[layer] * 1e2, 2)
        
        # correlation length of snow
        if 'lex' in variables:
            varian.append([lex[0] for lex in self.sd_lex])
            varian.append([lex[1] for lex in self.sd_lex])
            ind = [i for i, j in enumerate(var) if j == 'lex']
            x_ap[ind] = lex_ap
            #variances[ind] = lex_var
            for layer in list(range(nl_snow)):
                x_lowerlimit['lex' + str(layer)] = 0.01e-3 
                x_upperlimit['lex' + str(layer)] = 1.5e-3
                x_info['lex' + str(layer) + ' mm'] = np.round(lex_ap[layer] * 1e3, 2)
    
        # salinity of snow
        if 'sal' in variables:
            varian.append([sal[0] for sal in self.sd_sal])
            varian.append([sal[1] for sal in self.sd_sal])
            ind = [i for i, j in enumerate(var) if j == 'sal']
            x_ap[ind] = sal_ap
            #variances[ind] = sal_var
            for layer in list(range(nl_snow)):
                x_lowerlimit['sal' + str(layer)] = 0.001 * PSU
                x_upperlimit['sal' + str(layer)] = 20.0 * PSU
                x_info['sal' + str(layer) + ' PSU'] = np.round(sal_ap[layer]*1e3, 2)
    elif nls==1:
        # Snow depth
        if 'SD' in variables:
            varian.append((np.sum(self.sd, axis=1)))
            ind = [i for i, j in enumerate(var) if j == 'SD']
            x_ap[ind] = SD_ap[0]
            #variances[ind] = SD_var[0]
            x_lowerlimit['SD'] = 0.0
            x_upperlimit['SD' ] = 1.0
            x_info['SD cm'] = np.round(SD_ap * 1e2, 2)
        
        # correlation length of snow
        if 'lex' in variables:
            varian.append([lex[1]+lex[0] for lex in self.sd_lex])
            ind = [i for i, j in enumerate(var) if j == 'lex']
            x_ap[ind] = lex_ap[0]
            print(x_ap[ind])
            #variances[ind] = lex_var[0]
            x_lowerlimit['lex'] = 0.01e-3 
            x_upperlimit['lex'] = 0.9e-3
            x_info['lex mm'] = np.round(lex_ap * 1e3, 2)
        
    
        
        # salinity of snow
        if 'sal' in variables:
            #varian.append([sal[0] for sal in self.sd_sal])
            varian.append([sal[1] for sal in self.sd_sal])
            ind = [i for i, j in enumerate(var) if j == 'sal']
            x_ap[ind] = sal_ap[1]
            #variances[ind] = sal_var[0]
            x_lowerlimit['sal'] = 0.00 * PSU
            x_upperlimit['sal'] = 20.0 * PSU
            x_info['sal PSU'] = np.round(sal_ap*1e3, 2)



        # Density of snow
        if 'rho' in variables and 'rho1' not in variables:
            varian.append([rho[1] for rho in self.sd_rho])
            #varian.append([rho[1] for rho in self.sd_rho])
            ind = [i for i, j in enumerate(var) if j == 'rho']
            x_ap[ind] = rho_ap[1]
            #variances[ind] = rho_var
            x_lowerlimit['rho' ] = 150  # kg/m^3
            x_upperlimit['rho' ] = 450  # kg/m^3
            x_info['rho kg/m^3'] = np.round(rho_ap, 2)
                
                # Density of snow
        if 'rho1' in variables:
            varian.append([rho[1] for rho in self.sd_rho])
            #varian.append([rho[1] for rho in self.sd_rho])
            ind = [i for i, j in enumerate(var) if j == 'rho1']
            x_ap[ind] = rho_ap[1]
            print(x_ap)
            x_lowerlimit['rho1' ] = 150  # kg/m^3
            x_upperlimit['rho1' ] = 450  # kg/m^3
            x_info['rho kg/m^3'] = np.round(rho_ap, 2)
    

    if 'Tsur' in variables:
        varian.append(self.skt)
        ind = [i for i, j in enumerate(var) if j == 'Tsur']
        x_ap[ind] = Tsur_ap
        #variances[ind] = Tsur_var
        x_lowerlimit['Tsur'] = 273.15 - 30  # K
        x_upperlimit['Tsur'] = 273.15 - 5  # K
        x_info['Tsur K'] = np.round(Tsur_ap, 2)    
    ## ICE VARIABLES
    if 'SIT' in variables:
        varian.append(self.sit)
        ind = [i for i, j in enumerate(var) if j == 'SIT']
        x_ap[ind] = SIT_ap
        #variances[ind] = SIT_var
        x_lowerlimit['SIT'] = 0.0
        x_upperlimit['SIT'] = 5.0
        x_info['SIT cm'] = SIT_ap * 1e2
    
    # temperature at snow/ice interface
    if 'Tice' in variables:
        # varian.append(self.Tice)
        varian.append(self.Tice)
        ind = [i for i, j in enumerate(var) if j == 'Tice']
        x_ap[ind] = Tice_ap
        #variances[ind] = Tice_var
        x_lowerlimit['Tice'] = 273.15 - 20 # K
        x_upperlimit['Tice'] = 273.15  # K
        x_info['Tice K'] = Tice_ap
    
    # print(np.cov(np.array(varian)))
    #print(np.cov(np.array(varian))[0,0]*1e2)
    #print(np.cov(np.array(varian))[1,1]*1e3)
    # print(np.sqrt(np.cov(np.array(varian))[0,0])*1e2)
    ## Check if logaritmic is true
    if logarithmic==True:
        for key, count in zip(x_lowerlimit.keys(), range(len(x_lowerlimit.keys()))):
                
            if 'SD' in key or 'lex' in key  or'sal' in key:
                print('Convertin '+key+' variances to log10')
                x_lowerlimit[key]=np.log10(x_lowerlimit[key])
                x_upperlimit[key]=np.log10(x_upperlimit[key])
                #print(varian[count])
                varian[count] = np.log10(varian[count])
                #print(10**varian[count])
                
    ## fill covariance matrix
    x_cov = np.cov(np.array(varian))

    
    # roughness of snow/ice interface 
    if 'slope' in variables:
        ind = [i for i, j in enumerate(var) if j == 'slope']
        x_ap[ind] = 0.19  #m_ap
        #variances[ind] = m_var
        x_lowerlimit['slope'] = 0.01     # m/m
        x_upperlimit['slope'] = 0.3   # m/m  
        x_info['slope'] = 0.19 #m_ap
        x_covariance = np.zeros((x_cov.shape[0]+1,x_cov.shape[1]+1))
        # x_covariance[:-1, :-1] = x_cov
        x_covariance[-1, -1] = 0.05**2
        x_cov = x_covariance
    
    # print(x_ap.shape)
    # print(x_cov.shape)
    
    #print(x_lowerlimit)
    #print(x_upperlimit)
    return [x_ap, x_cov, x_lowerlimit, x_upperlimit, x_info]


def prepare_inputs(s, variables):

    
    SIT_ap = 2.0     # m
    SD_ap  = 0.10    # m
    lex_ap = 0.17e-3 # m
    sal_ap = 3 * PSU
    rho_ap = 300     # kg/m^3
    Tice_ap = 260    # K
    Tsnow_ap = 253   # K
    m_ap = Get_slope(s.site, variance=False)[0]  # roughness measured at site
    
    SD_var    = np.std(s.SD_tot)**2
    lex_var   = np.std(s.corr_m)**2
    sal_var   = np.std(s.sal_m)**2
    rho_var   = np.std(s.rho_m)**2
    Tice_var  =  np.std(s.Tcal_ice)**2
    Tsnow_var =  np.std(np.concatenate((s.Tcal_snow)))**2
    SIT_var   =  0.10
    m_var = Get_slope(s.site, variance=True)[1]  # variability of roughness values
    
    x_ap = np.zeros(len(variables)) # starting guess for SD and SIT
    variances = np.zeros(len(variables)) # variances of variables
    
    ## Lower and upper limits of variables
    x_lowerlimit = dict.fromkeys(variables)
    x_upperlimit = dict.fromkeys(variables)
    x_info = {}
    
    if 'SD' in variables:
        variances[variables.index('SD')] = SD_var
        x_ap[variables.index('SD')] = SD_ap
        x_lowerlimit['SD'] = 0.0  # m
        x_upperlimit['SD'] = 1.0  # m
        x_info['SD cm'] = SD_ap*1e2  # cm
        # x_info['SD'] = x_info.pop('SD cm')
    if 'lex' in variables:
        variances[variables.index('lex')] = lex_var
        x_ap[variables.index('lex')] = lex_ap
        x_lowerlimit['lex'] = 0.05e-3  # m
        x_upperlimit['lex'] = 0.7e-3   # m
        x_info['lex mm'] = lex_ap*1e3     # mm
        # x_info['lex'] = x_info.pop('lex mm')
    if 'sal' in variables:
        variances[variables.index('sal')] = sal_var
        x_ap[variables.index('sal')] = sal_ap
        x_lowerlimit['sal'] = 0.0 * PSU  
        x_upperlimit['sal'] = 20.0 * PSU 
        x_info['sal PSU'] = sal_ap*1e3       # PSU
        # x_info['sal'] = x_info.pop('sal PSU')
    if 'rho' in variables:
        variances[variables.index('rho')] = rho_var
        x_ap[variables.index('rho')] = rho_ap
        x_lowerlimit['rho'] = 150  # kg/m^3
        x_upperlimit['rho'] = 450  # kg/m^3
        x_info['rho kg/m3'] = rho_ap     # kg/m^3
        # x_info['rho'] = x_info.pop('rho kg/m3')
    ## temperature at snow/ice interface
    if 'Tice' in variables:
        variances[variables.index('Tice')] = Tice_var
        x_ap[variables.index('Tice')] = Tice_ap
        x_lowerlimit['Tice'] = 273.15 - 20 # K
        x_upperlimit['Tice'] = 273.15  # K
        x_info['Tice K'] = Tice_ap     # K
    if 'Tsnow' in variables:
        variances[variables.index('Tsnow')] = Tsnow_var
        x_ap[variables.index('Tsnow')] = Tsnow_ap
        x_lowerlimit['Tsnow'] = 243.15 - 30  # K
        x_upperlimit['Tsnow'] = 273.15 - 5  # K
        x_info['Tsnow K'] = Tsnow_ap     # K
    if 'SIT' in variables:
        variances[variables.index('SIT')] = SIT_var
        x_ap[variables.index('SIT')] = SIT_ap
        x_lowerlimit['SIT'] = 0.1   # m
        x_upperlimit['SIT'] = 6.0   # m    
        x_info['SIT cm'] = SIT_ap *1e2 # cm
    if 'slope' in variables:
        variances[variables.index('slope')] = m_var
        x_ap[variables.index('slope')] = m_ap
        x_lowerlimit['slope'] = 0     # m/m
        x_upperlimit['slope'] = 1.0   # m/m  
        x_info['slope'] = m_ap

    ## fill covariance matrix
    x_cov = np.diag(variances)
    
    return [x_ap, x_cov, x_lowerlimit, x_upperlimit, x_info]


def prepare_inputs_FYI_OIB2(self, variables, Nd=100, nl_snow=2, sim=False, nls=2):
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
    # snow
    SD_ap  = np.mean(self.sd, axis=0)         # m

    lex_ap = np.mean(self.sd_lex, axis=0)     # m
    sal_ap = np.mean(self.sd_sal, axis=0)
    ## convert to logarithmic units
    sal_ap = np.log10(sal_ap)
    lex_ap = np.log10(lex_ap)
    
    rho_ap = np.mean(self.sd_rho, axis=0)     # kg/m^3
    Tsur_ap = np.mean(self.skt, axis=0)       # K
    # ice
    SIT_ap = np.mean(self.sit, axis=0)         # m
    Tice_ap = np.mean(self.Tice, axis=0)      # K
    m = []
    for site in ['AK1', 'AK2', 'AK3', 'AK4']:
        m.append(Get_slope(site, variance=False)[0])
    m_ap = 0.08 #vnp.mean(m)
    
    ## expected variances
    # snow
    SD_var    = np.array([np.std(self.sd)**2] * nl_snow)
    lex_var   = np.std(self.sd_lex, axis=0)**2
    sal_var   = np.std(self.sd_sal, axis=0)**2
    rho_var   = np.std(self.sd_rho, axis=0)**2
    Tsur_var = np.std(self.skt)**2
    # ice
    SIT_var   = np.std(self.sit, axis=0)**2
    m_var = np.std(m)
    
    if nls==2:
        var = [i for s in [nl_snow*[var] if var !='SD' and var != 'slope' and var != 'Tsur' and var != 'SIT' and var !='Tice' else [var] for var in variables ] for i in s]
        varnames = [v + str(i) if v != 'SD' and v != 'slope' and v != 'Tsur' and v != 'SIT' and v !='Tice' else v for v,i in zip(var, list(range(nl_snow))*len(var))]
    elif nls==1:
        var = np.array([[var] for var in variables]).flatten()
        varnames = var

        ## make array with correct number of snow elements
        
    # var = [i for s in [nl_snow*[var] if  var != 'lex' and var != 'sal' and var != 'slope' and var != 'Tsur' and var != 'SIT' and var !='Tice' else [var] for var in variables ] for i in s]
    ## make array of unique varnames
    # varnames = [v + str(i) if v != 'lex' and v != 'sal' and v != 'slope' and v != 'Tsur' and v != 'SIT' and v !='Tice' else v for v,i in zip(var, list(range(nl_snow))*len(var))]
    x_ap = np.zeros(len(var)) # starting guess array
    variances = np.zeros(len(var)) # variances of variables
    varian = []
    
    # i for s in used to flatten the list
    x_lowerlimit = dict.fromkeys(varnames)
    x_upperlimit = dict.fromkeys(varnames)
    x_info = {}
    
    ## SNOW VARIABLES
    if nls==2:
        # Snow depth
        if 'SD' in variables:
            varian.append((np.sum(self.sd, axis=1)))
            ind = [i for i, j in enumerate(var) if j == 'SD']
            x_ap[ind] = SD_ap[0]
            variances[ind] = SD_var[0]
            x_lowerlimit['SD'] = 0.0
            x_upperlimit['SD' ] = 1.0
            x_info['SD cm'] = np.round(SD_ap * 1e2, 2)
        
        # correlation length of snow
        if 'lex' in variables:
            varian.append([lex[0] for lex in self.sd_lex])
            varian.append([lex[1] for lex in self.sd_lex])
            ind = [i for i, j in enumerate(var) if j == 'lex']
            x_ap[ind] = lex_ap
            variances[ind] = lex_var
            for layer in list(range(nl_snow)):
                x_lowerlimit['lex' + str(layer)] = 0.01e-3 
                x_upperlimit['lex' + str(layer)] = 1.5e-3
                x_info['lex' + str(layer) + ' mm'] = np.round(lex_ap[layer] * 1e3, 2)
    
        # salinity of snow
        if 'sal' in variables:
            varian.append([sal[0] for sal in self.sd_sal])
            varian.append([sal[1] for sal in self.sd_sal])
            ind = [i for i, j in enumerate(var) if j == 'sal']
            x_ap[ind] = sal_ap
            variances[ind] = sal_var
            for layer in list(range(nl_snow)):
                x_lowerlimit['sal' + str(layer)] = 0.001 * PSU
                x_upperlimit['sal' + str(layer)] = 20.0 * PSU
                x_info['sal' + str(layer) + ' PSU'] = np.round(sal_ap[layer]*1e3, 2)
        # Density of snow
        if 'rho' in variables:
            varian.append([rho[0] for rho in self.sd_rho])
            varian.append([rho[1] for rho in self.sd_rho])
            ind = [i for i, j in enumerate(var) if j == 'rho']
            x_ap[ind] = rho_ap
            variances[ind] = rho_var
            for layer in list(range(nl_snow)):
                x_lowerlimit['rho' + str(layer)] = 150  # kg/m^3
                x_upperlimit['rho' + str(layer)] = 450  # kg/m^3
                x_info['rho' + str(layer) + ' kg/m^3'] = np.round(rho_ap[layer], 2)
    elif nls==1:
        # Snow depth
        if 'SD' in variables:
            varian.append([sd[0] for sd in self.sd])
            ind = [i for i, j in enumerate(var) if j == 'SD']
            x_ap[ind] = SD_ap[0]
            variances[ind] = SD_var[0]
            x_lowerlimit['SD'] = 0.0
            x_upperlimit['SD' ] = 1.0
            x_info['SD cm'] = np.round(SD_ap * 1e2, 2)
        
        # correlation length of snow
        if 'lex' in variables:
            varian.append([np.mean(lex) for lex in np.log10(self.sd_lex)])
            ind = [i for i, j in enumerate(var) if j == 'lex']
            x_ap[ind] = lex_ap[0]
            variances[ind] = lex_var[0]
            x_lowerlimit['lex'] = np.log10(0.01e-3) 
            x_upperlimit['lex'] = np.log10(0.9e-3)
            x_info['lex mm'] = np.round(lex_ap * 1e3, 2)
        
    
        
        # salinity of snow
        if 'sal' in variables:
            #varian.append([sal[0] for sal in self.sd_sal])
            # varian.append([np.mean(sal) for sal in self.sd_sal])
            varian.append([np.mean(sal) for sal in np.log10(self.sd_sal)])
            ind = [i for i, j in enumerate(var) if j == 'sal']
            x_ap[ind] = sal_ap[0]
            variances[ind] = sal_var[0]
            #x_lowerlimit['sal'] = 0.00 * PSU
            #x_upperlimit['sal'] = 20.0 * PSU
            x_lowerlimit['sal'] = np.log10(0.00 * PSU)
            x_upperlimit['sal'] = np.log10(20.0 * PSU)
            x_info['sal PSU'] = np.round(sal_ap*1e3, 2)


    
        # Density of snow
        if 'rho' in variables:
            # varian.append([rho[0] for rho in self.sd_rho])
            varian.append([np.mean(rho) for rho in self.sd_rho])
            ind = [i for i, j in enumerate(var) if j == 'rho']
            x_ap[ind] = rho_ap[1]
            variances[ind] = rho_var[1]
            for layer in list(range(nl_snow)):
                x_lowerlimit['rho' + str(layer)] = 150  # kg/m^3
                x_upperlimit['rho' + str(layer)] = 450  # kg/m^3
                x_info['rho' + str(layer) + ' kg/m^3'] = np.round(rho_ap[layer], 2)
    

    if 'Tsur' in variables:
        varian.append(self.skt)
        ind = [i for i, j in enumerate(var) if j == 'Tsur']
        x_ap[ind] = Tsur_ap
        variances[ind] = Tsur_var
        x_lowerlimit['Tsur'] = 273.15 - 30  # K
        x_upperlimit['Tsur'] = 273.15 - 5  # K
        x_info['Tsur K'] = np.round(Tsur_ap, 2)    
    ## ICE VARIABLES
    if 'SIT' in variables:
        varian.append(self.sit)
        ind = [i for i, j in enumerate(var) if j == 'SIT']
        x_ap[ind] = SIT_ap
        variances[ind] = SIT_var
        x_lowerlimit['SIT'] = 0.0
        x_upperlimit['SIT'] = 5.0
        x_info['SIT cm'] = SIT_ap * 1e2
    
    ## fill covariance matrix
    x_cov = np.cov(np.array(varian))
    # roughness of snow/ice interface 
    if 'slope' in variables:
        ind = [i for i, j in enumerate(var) if j == 'slope']
        x_ap[ind] = m_ap
        variances[ind] = m_var
        x_lowerlimit['slope'] = 0     # m/m
        x_upperlimit['slope'] = 1.0   # m/m  
        x_info['slope'] = m_ap
        x_covariance = np.zeros((x_cov.shape[0]+1,x_cov.shape[1]+1))
        x_covariance[:-1, :-1] = x_cov
        x_covariance[-1, -1] = 0.005
        x_cov = x_covariance

    return [x_ap, x_cov, x_lowerlimit, x_upperlimit, x_info]