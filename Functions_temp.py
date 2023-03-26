# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 13:18:25 2022

@author: Ida Olsen

Compute "ideal" temperature profile
"""
import numpy as np
import pickle
import matplotlib.pyplot as plt
import os
import sys
import warnings
warnings.filterwarnings("ignore")
from smrt.core.globalconstants import PSU

# home made
from Snowpit_class import Snowpit_info
from polar_plots_pred import plot as pp
#append path
sys.path.append('C:/Users/Ida Olsen/Documents/Speciale_2022/scripts')
data_path = 'C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/Code/'

def haversine(file, Olat,Olon, Dlat,Dlon):

    radius = 6371e3  # km

    d_lat = np.radians(Dlat - Olat)
    d_lon = np.radians(Dlon - Olon)
    a = (np.sin(d_lat / 2.) * np.sin(d_lat / 2.) +
         np.cos(np.radians(Olat)) * np.cos(np.radians(Dlat)) *
         np.sin(d_lon / 2.) * np.sin(d_lon / 2.))

    c = 2. * np.arctan2(np.sqrt(a), np.sqrt(1. - a))
    d = radius * c
    # if 'AlongTrack' in file:
    #     if d_lon > 0: # if flight is west of pit
    #         d = -d
    # elif 'AcrossTrack' in file:
    #     if d_lat > 0: # if flight is south
    #         d = -d
    return d

def Get_SIT(self, site, plot=False):
    " Get SIT values taken at CB"
    
    directory = 'C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/Ice_thickness_Cambay/EM31 Data - comma delimited'
    files = os.listdir(directory)
    ice_comb = []
    data_loc = np.where([site in file for file in files])[0]
    for loc in data_loc:
        file = os.path.join(directory, files[loc])
        data = np.genfromtxt(file, delimiter=',', names=True, encoding='utf8')
        ice_comb.append(data['hi_m_'])

        self.SIT = np.concatenate((ice_comb))
        if plot==True: ## plot distribution of SIT
            sit = data['hi_m_']
            ## find count closest to location
            coords_flight = np.array([(lat,lon) for lat,lon in zip(data['lat'], data['lon'])])
            site_coord = np.array((self.lat, self.lon))
            distances = np.linalg.norm(coords_flight-site_coord, axis=1)
            dist = [haversine(file, self.lat, self.lon, Dlat,Dlon)for Dlat,Dlon in zip(data['lat'], data['lon'])]
            min_index = np.argmin(abs(np.array(dist)))
            #print(site_coord)
            #print(coords_flight[min_index])
            
            fig,ax = plt.subplots(figsize=(7,4))
            plt.scatter(dist, sit, s=5)
            plt.scatter(dist[min_index], np.mean(self.SIT), s=50, color='r', label='Location'+ self.site)
            plt.xlabel('distance [m] from pit')
            plt.ylabel('SIT [m]')
            plt.ylim(1,6)
            plt.grid()
            plt.legend()
            
            # plt.xlim(0,np.max(list(range(len(data['lon'])))))
            
            ## polar plot
            extent = [np.min(data['lon']), np.max(data['lon']),np.min(data['lat']), np.max(data['lat'])]
            if "AcrossTrack" in file:
                plt.title('Distribution of SIT at site: '+ site + ', AcrossTrack flight')
                title = 'AcrossTrack Flight path at site: '+ self.site
                plt.savefig(os.getcwd() + '/Testing_scripts/figures/' + self.site + 'SIT_ AcrossTrack.PNG')
            else:
                plt.title('Distribution of SIT at site: '+ site + ', AlongTrack flight')
                title = 'AlongTrack Flight path at site: '+ self.site
                plt.savefig(os.getcwd() + '/Testing_scripts/figures/' + self.site + 'SIT_ AlongTrack.PNG')
            pp(self, data['lat'], data['lon'], dist,ylabel='distance [m] from pit', title=title, s=10)
            plt.show()

# Load refernece data from AMSR with unc.
def load_ref_data(site, land_data=False):
    " load satelitte data (ECMWF and AMSR2)"
    Base_path = "C:/Users/Ida Olsen/Documents/Speciale_2022/scripts/npy_files/"
    satdataH = np.load(Base_path + site + '_TBH.npy')
    satdataV = np.load(Base_path + site + '_TBV.npy')
    satdataHErr = np.load(Base_path + site + '_TBH_err.npy')
    satdataVErr = np.load(Base_path + site + '_TBV_err.npy')
    
    satdataH_landfree = np.load(Base_path + site + '_TBH_Landfree.npy')
    satdataV_landfree = np.load(Base_path + site + '_TBV_Landfree.npy')
    satdataH_landfreeErr = np.load(Base_path + site + '_TBH_Landfree_Err.npy')
    satdataV_landfreeErr = np.load(Base_path + site + '_TBV_Landfree_Err.npy')
    
    data = np.zeros((6,len(satdataV)))
    dataErr = np.zeros((6,len(satdataV)))
    ## mean data
    if land_data==True:
        data[0,:] = satdataV
        data[1,:] = satdataH
    data[2,:] = satdataV_landfree.flatten()
    data[3,:] = satdataH_landfree.flatten()
    
    ## Errors
    if land_data==True:
        dataErr[0,:] = satdataVErr
        dataErr[1,:] = satdataHErr
    dataErr[2,:] = satdataV_landfreeErr.flatten()
    dataErr[3,:] = satdataH_landfreeErr.flatten()
    
    return [data, dataErr]

# Fill snowpit class with data
def fill_data():
    with open(data_path + 'smrt_in-out/smp_profile_CB', 'rb') as pickle_file:
        smp_profile = pickle.load(pickle_file)
    
    lats = [69.08897, 68.99985, 68.91067,68.82158]			
    lons = [-105.50422,-105.47468,-105.44359,-105.41325]
    dates = [' 20/04/2022',' 21/04/2022',' 23/04/2022',' 24/04/2022']
    Tsurfaces = [-22+273.15, -20+273.15, -18+273.15, -20+273.15]
    sites = ['AK1', 'AK2', 'AK3', 'AK4']
    
    comb_info =[];
    name_prev = False
    count = -1
    for site, lat, lon, date, Tsurface in zip(sites, lats, lons, dates, Tsurfaces):
        count +=1
        s = Snowpit_info(site)
        s.Tsur = Tsurface
        s.lat = lat
        s.lon = lon
        s.date = date
        s.Tsur_t2m = np.load(data_path + site + '_ECMWF_t2m.npy', allow_pickle=True)
        s.Tsur_skt = np.load(data_path + site + '_ECMWF_skt.npy', allow_pickle=True)
        Get_SIT(s, site, plot=False) ## Get SOT from data
        # print(s.Tsur_ECMWF)
        for name, profile in smp_profile[site].items():
            # Convert to exp. corr length
            debye = 1.5
            lex_array = debye * 4 * (1 - profile.density / 917) / (profile.ssa * 917)
            ## Get mean values
            if name == name_prev:
                s.SD   = np.append(s.SD, profile.thick)
                s.rho  = np.append(s.rho,profile.density)
                s.corr = np.append(s.corr,lex_array)
                s.temp = np.append(s.temp,profile.temperature)
                s.sal  = np.append(s.sal,profile.salinity *PSU)
            else:
                s.SD.append(profile.thick)
                s.rho.append(profile.density)
                s.corr.append(lex_array)
                s.temp.append(profile.temperature)
                s.sal.append(profile.salinity *PSU)
            
            s.get_mean()
            name_prev = name                         
        
        comb_info.append(s)
        ## Calculate SIT without snow
        s.SIT_m = s.SIT_m - np.mean(s.SD_tot)
        # print(site + ' SIT:', s.SIT_m)
    return comb_info


# Fill snowpit class with data
def fill_data_ALERT_EUREKA():
    with open(data_path + 'ALERT_EUREKA/smp_process_ida.txt', 'rb') as pickle_file:
        smp_profile = pickle.load(pickle_file)   
    comb_info_ALERT =[];
    # comb_info_EUREKA =[];
   
    for count, site in zip(range(len(smp_profile.keys())), smp_profile.keys()):
        # if 'Eureka' in site:
        #     s = Snowpit_info(site)
        #     s.name = []
        #     # source - https://tc.copernicus.org/articles/14/4323/2020/
        #     if '2' or '3' in smp_profile.keys():
        #         s.date = ' 08/04/2020'
        #     elif '5' or '6' in smp_profile.keys():
        #         s.date = ' 13/04/2020'
        #     elif '7' or '8' in smp_profile.keys():
        #         s.date = ' 12/04/2020'
        if 'Alert' in site:
            s = Snowpit_info(site)
            s.name = []
            # source - https://tc.copernicus.org/articles/14/4323/2020/
            if '2' or '3' in site:
                s.date = ' 11/04/2020'
            elif '5' or '6' in site:
                s.date = ' 13/04/2020'
            elif '7' or '8' in site:
                s.date = ' 12/04/2020'
            for p in site:    
                
                # s.Tsur = Tsurface (har vi temperatur somewhere?)
                # s.date = date og dato??
                s.name.append(p['name'])
                s.lat.append(p['lat'])
                s.lon.append(p['lon'])
                
                s.distance.append(p['smp_data']['distance'])
                s.dens.append(p['smp_data']['density'])
                
                # Convert to exp. corr length
                debye = 1.75
                lex_array = debye * 4 * (1 - p['smp_data']['density'] / 917) / (p['smp_data']['ssa'] * 917)
                s.corr.append(lex_array)
                s.SD.append(p['smp_data']['thick'])
                
                s.get_mean()                         
            
            comb_info_ALERT.append(s)
    return comb_info_ALERT


# Make temperature profile
def temp_profile(self, k, SD, ice, tsk, nl_snow, nl_ice=1, plot=False):
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
    # ks = 0.4  # Thermal conductivity of snow [w/(m*K)]
    # print(self.rho)
    rho = np.concatenate(self.rho)
    # print(rho)
    ks_array = 2.5*1e-6*rho**2-1.23*1e-4*rho+0.024
    # print('ks', np.mean(ks))
    ks = np.mean(ks_array)
    # ks = 0.4
    print(ks)
    
    
    ## snow ice interface temperature
    T_sii = (k * SD * T_w + ks * ice *
                tsk) / (k * SD + ks * ice)
                                                                
    
    tsk = tsk  # Surface temperature
    
    ## make temperature profile of snow and ice and concatenate
    temp_profile_ice = np.linspace(T_sii, T_w, nl_ice)
    temp_profile_snow = np.linspace(tsk, T_sii, nl_snow+2)[1:-1]
    

    
    if plot==True:
        # ice/snow thickness arrays for showing temperature profiles    
        ice_thick = np.linspace(-ice, 0, nl_ice)
        snow_thick = np.linspace(0,  nl_snow)
        
        plt.figure()
        plt.plot(temp_profile_snow, snow_thick,'.-')
        plt.plot(np.array([T_sii,temp_profile_ice]), np.array([0, ice_thick]),'.-')
        plt.grid()
        plt.xlabel('Temperature [K]')
        plt.ylabel('Depth [m]')
        
    
    return [T_sii,temp_profile_snow,temp_profile_ice]

def Calc_temp_profiles(plot=False):
    comb_info = fill_data()
    for c in comb_info:
        k_MYI = 1.88 # MYI
        k_FYI = 2.10 # FYI
        
        for profile in c.SD:
            # print(c.Tsur)
            [T_sii,temp_profile_snow,temp_profile_ice] = temp_profile(c, k_FYI, np.sum(profile), c.SIT_m, c.Tsur, len(profile))
            
            c.Tcal_snow.append(temp_profile_snow)
            c.Tcal_ice.append(T_sii)
            
        import matplotlib.colors as mcolors
        
        ## plot temperature profiles for each site
        # plots are saved to figures/profiles folder
        if plot==True:
            plt.figure()
            ice_thick = np.array([0, -c.SIT_m])
            T_w = (273.15-1.8)  # K
            
            colors = np.array([c for c in mcolors.CSS4_COLORS])
            for i in range(len(c.SD)):
                snow_thick = np.linspace(np.sum(c.SD[i]),0, len(c.SD[i])+1)
                plt.plot(np.append(c.Tcal_snow[i], c.Tcal_ice[i]), snow_thick,'.-', color=colors[i])
                plt.plot(np.array([c.Tcal_ice[i],T_w]) , ice_thick,'.-', color=colors[i])
            plt.grid()
            plt.xlabel('Temperature [K]')
            plt.ylabel('Depth [m]')
            plt.title('Temperature profiles for site: ' + c.site)
            plt.xlim(250, 274)
            savepath = "C:/Users/Ida Olsen/Documents/Speciale_2022/figures/profiles"
            savefile = os.path.join(savepath, c.site + '_Temperature_profiles.png')
            plt.savefig(savefile)

    return comb_info

def compute_bias_std(self, data, nullindex=True):
    if nullindex==False:
        biasV = np.round(data[0]-data[2], 2)
        biasH = np.round(data[1]-data[3], 2)
    else:
        biasV = np.round(data[2]-data[4], 2)
        biasH = np.round(data[3]-data[5], 2)
    
    import tabulate
    # header = ['6.9', '10.7', '18.7', '23.8', '36.5']
    data = [[self.site + ' Bias, TBV', biasV[0],biasV[1],biasV[2],biasV[3],biasV[4]],
            [self.site + ' Bias, TBH', biasH[0],biasH[1],biasH[2],biasH[3],biasH[4]]]
    print(tabulate.tabulate(data, tablefmt='presto'))
    
    if self.site=='AK1':
        file1 = open('bias_Tmeasured.txt', 'w')
        l = ['6.9  & 10.7 & 18.7 & 23.8 & 36.5 \n']
        file1.writelines(l)
    else:
        file1 = open('bias_test.txt', 'a')
    ll = [str(bV) + ' &' for bV in biasV] + ['\n']
    lll = [str(bH) + ' &' for bH in biasH] + ['\n']
    
    file1.writelines(ll)
    file1.writelines(lll)
    file1.close()
  

    # return biasV, biasH
    
    

# snowpit_data = fill_data()
# from Simulated_data import SMRT_simulations
# SIM_FYI = SMRT_simulations("FYI", snowpit_data, os.path.dirname(os.getcwd()))
# SIM_FYI.set_vars(plot=False)
# # SIM_FYI.simulated_data(num_simulations=100, save=True)



# ## Get CB values
# SD_CB = np.concatenate([pit.SD for pit in snowpit_data])

# lex_CB = np.concatenate([pit.corr for pit in snowpit_data])
# rho_CB = np.concatenate([pit.rho for pit in snowpit_data])
# sal_CB = np.concatenate([pit.sal for pit in snowpit_data])

# SD_CB_tot = np.concatenate([pit.SD_tot for pit in snowpit_data])

# ## extract values to be a part of top and bottom indexes
# index_top = []
# index_bottom = []
# count = 0
# for sd in SIM_FYI.sd:
#     count +=1
#     # which profiles to include - those within a 10 cm interval of the true depth
#     index1 = [SD_CB_tot <= sd+0.05]
#     index2 = [SD_CB_tot >= sd-0.05]
#     index = np.logical_and(index1, index2)[0]
#     # index = (np.array(index1) & np.array(index2)).flatten()
    
#     ## calculate cummulated sum of each profile
#     cumsum = np.concatenate([np.cumsum(sd_cb) for sd_cb in SD_CB[index]])
    
#     ## include corr lengths of all values above 3cm of sd
#     index_t = [cs<=(sd-0.03) for cs in cumsum]
#     ## incldue corr lengths below 3cm of snowdepth e.g.
#     index_b = [cs>=(sd-0.03) for cs in cumsum]
    
#     # get correlation lengths
#     lex_top = np.concatenate(lex_CB[index])[index_t]
#     lex_bottom = np.concatenate(lex_CB[index])[index_b]
    
#     # get densities
#     rho_top = np.concatenate(rho_CB[index])[index_t]
#     rho_bottom = np.concatenate(rho_CB[index])[index_b]
    
#     # get densities
#     # index_tsal = [sal<=(sd-0.03) for sal in sal_CB[index]]
#     # index_bsal = [sal<=(sd-0.03) for sal in sal_CB[index]]
#     sal_top = np.concatenate([sal[:-1] for sal in sal_CB[index]])
#     sal_bottom = np.concatenate([sal[-1:] for sal in sal_CB[index]])
#     cumsum_sal_t = np.concatenate([np.cumsum(sd_cb)[:-1] for sd_cb in SD_CB[index]])
#     cumsum_sal_b = np.concatenate([np.cumsum(sd_cb)[-1:] for sd_cb in SD_CB[index]])
    
    
    
#     if count % 30 == 0: 
#         cumsum_ref = np.concatenate([np.cumsum(sd_cb) for sd_cb in SD_CB])
#         plt.figure()
#         plt.gca().invert_yaxis()
#         plt.scatter(np.concatenate(lex_CB)*1e3, cumsum_ref, s=10, c='blue', label='lex all')
#         x = np.linspace(np.min(np.concatenate(lex_CB)*1e3),np.max(np.concatenate(lex_CB)*1e3), len(np.concatenate(lex_CB)))
#         y = [sd for lex in np.concatenate(lex_CB)]
#         yb = [sd-0.1 for lex in np.concatenate(lex_CB)]
#         yt = [sd+0.1 for lex in np.concatenate(lex_CB)]
#         plt.title('Correlation lengths for profile with SD, ' + str(sd) + ' [m]' )
#         plt.plot(x, y, label='SD OIB')
#         plt.plot(x, yb, label='lower lim CB SD')
#         plt.plot(x, yt, label='upper lim CB SD')
#         plt.scatter(lex_top*1e3, cumsum[index_t], s=10, c='red', label='lex top layer')
#         plt.scatter(lex_bottom*1e3, cumsum[index_b], s=10, c='green', label='lex bottom layer')
#         plt.xlim([np.min(np.concatenate(lex_CB)*1e3),np.max(np.concatenate(lex_CB)*1e3)])
#         plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
#         plt.xlabel('lex - unit mm')
#         plt.ylabel('SD [m]')
#         plt.show()
        
#         plt.figure()
#         plt.gca().invert_yaxis()
#         plt.scatter(np.concatenate(rho_CB), cumsum_ref, s=10, c='blue', label='rho all')
#         x = np.linspace(np.min(np.concatenate(rho_CB)),np.max(np.concatenate(rho_CB)), len(np.concatenate(rho_CB)))
#         y = [sd for rho in np.concatenate(rho_CB)]
#         yb = [sd-0.1 for rho in np.concatenate(rho_CB)]
#         yt = [sd+0.1 for rho in np.concatenate(rho_CB)]
#         plt.title('Densities for profile with SD, ' + str(sd) + ' [m]' )
#         plt.plot(x, y, label='SD OIB')
#         plt.plot(x, yb, label='lower lim CB SD')
#         plt.plot(x, yt, label='upper lim CB SD')
#         plt.scatter(rho_top, cumsum[index_t], s=10, c='red', label='rho top layer')
#         plt.scatter(rho_bottom, cumsum[index_b], s=10, c='green', label='rho bottom layer')
#         plt.xlim([np.min(np.concatenate(rho_CB)),np.max(np.concatenate(rho_CB))])
#         plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
#         plt.xlabel('rho - unit mm')
#         plt.ylabel('SD [m]')
#         plt.show()
        
#         plt.figure()
#         plt.gca().invert_yaxis()
#         plt.scatter(np.concatenate(sal_CB), cumsum_ref, s=10, c='blue', label='sal all')
#         x = np.linspace(np.min(np.concatenate(sal_CB)),np.max(np.concatenate(sal_CB)), len(np.concatenate(sal_CB)))
#         y = [sd for sal in np.concatenate(sal_CB)]
#         yb = [sd-0.1 for sal in np.concatenate(sal_CB)]
#         yt = [sd+0.1 for sal in np.concatenate(sal_CB)]
#         plt.title('salinities for profile with SD, ' + str(sd) + ' [m]' )
#         plt.plot(x, y, label='SD OIB')
#         plt.plot(x, yb, label='lower lim CB SD')
#         plt.plot(x, yt, label='upper lim CB SD')
#         plt.scatter(sal_top, cumsum_sal_t, s=10, c='red', label='sal top layer')
#         plt.scatter(sal_bottom, cumsum_sal_b, s=10, c='green', label='sal bottom layer')
#         plt.xlim([np.min(np.concatenate(sal_CB)),np.max(np.concatenate(sal_CB))])
#         plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
#         plt.xlabel('sal - unit mm')
#         plt.ylabel('SD [m]')
#         plt.show()
        