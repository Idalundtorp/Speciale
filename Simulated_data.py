# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 11:48:11 2022

@author: Ida Olsen
"""
import netCDF4 as nc
import os
import numpy as np
import matplotlib.pyplot as plt

## SMRT modules
from smrt import make_snowpack, make_ice_column, make_model, sensor_list
# import datetime as dt
#import smrt
#from smrt.substrate.reflector import make_reflector
from smrt.core.globalconstants import  PSU
# from smrt.permittivity.saline_snow import saline_snow_permittivity_geldsetzer09 as ssp
from smrt.permittivity.saline_snow import saline_snow_permittivity_scharien_with_stogryn95 as ssp
# from smrt.core.sensitivity_study import SensitivityStudy

## Home made modules     
from polar_plots_pred import plot_one as po
from polar_plots_pred import polar_plot as pp
from roughness import roughness
from Make_plots import plot_simulation_parameters_from_CB
from Functions_temp import fill_data
# plot path
pp = "C:/Users/Ida Olsen/Documents/Speciale_2022/figures"
  
class data():
    
    # init method or constructor
    def __init__(self, filename, dtype, ice_type=-1, empty=False):
        if empty == False:
            names = ["latitude","longitude","time","identification","sd_mean","sd_std","SDln","sit_mean","sit_std","SITln","FRB","FRBstd","Tsur","Tair","wsd","wdensity","lat1","lon1","time1","reference-id","69GHzH","69GHzV","73GHzH","73GHzV","107GHzH","107GHzV","187GHzH","187GHzV","238GHzH","238GHzV","365GHzH","365GHzV","890GHzH","890GHzV","Earth Incidence","Earth Azimuth","scanpos","upstreamfile","timediff","latitude2","longitude2","time2","reference-id2","upstreamfile","msl","u10","v10","ws","t2m","skt","istl1","istl2","istl3","istl4","sst","d2m","tcwv","tclw","tciw","ssrd","strd","e","tp","sf","fal","siconc","latitude3","longitude3","time3","reference-id3","upstreamfile","sigma_40","sigma_40_mask","nb_samples","warning","std","latitude4","longitude4","time4","reference-id4","upstreamfile","SMOS_Tbv","SMOS_Tbh","SMOS_RMSE_v","SMOS_RMSE_h","SMOS_nmp","SMOS_data1loss"]
            #print(filename)
            if "-NERSC-" in filename:
                data1 = np.genfromtxt(filename, names=True, dtype=dtype,
                                     skip_header=1, delimiter=',', encoding='utf8')
            else:
                data1 = np.genfromtxt(filename, names=names, dtype=dtype,
                                     skip_header=1, delimiter=',', encoding='utf8')
    
            # self.subset   = []
            self.sd       = data1['sd_mean']
            self.sit      = data1['sit_mean']
    
            index_SD  = ~np.isnan(self.sd)
            index_SIT = ~np.isnan(self.sit)
            index = [(isd and isit) for isd,isit in zip(index_SD, index_SIT)]
            self.index = index
            if ice_type!=-1:
                index = np.where(self.index)[0][ice_type]
                self.index = index
                # print(index)
    
            self.sd       = data1['sd_mean'][index]
            self.sit      = data1['sit_mean'][index]
            self.lat      = data1['latitude'][index]
            self.lon      = data1['longitude'][index]
            self.date     = data1['time'][index]
            self.id       = data1['identification'][index]
            self.sd_std   = data1['sd_std'][index]
            self.sit_std  = data1['sit_std'][index]
            try:
                self.slope    = data1['surface_roughness_mean'][index]
                self.slope_std= data1['surface_roughness_std'][index]
            except:
                self.slope    = [np.nan for el in range(len(self.lat))]
                self.slope_std= [np.nan for el in range(len(self.lat))]
            self.t2m      = data1['t2m'][index]
            self.skt      = data1['skt'][index]
            self.tcwv     = data1['tcwv'][index]
            
    
            ## AMSR2 V
            self.tb6V     = data1['69GHzV'][index]
            self.tb7V     = data1['73GHzV'][index]
            self.tb10V    = data1['107GHzV'][index]
            self.tb18V    = data1['187GHzV'][index]
            self.tb23V    = data1['238GHzV'][index]
            self.tb36V    = data1['365GHzV'][index]
            self.tb89V    = data1['890GHzV'][index]
            
            ## AMSR2 H
            self.tb6H     = data1['69GHzH'][index]
            self.tb7H     = data1['73GHzH'][index]
            self.tb10H    = data1['107GHzH'][index]
            self.tb18H    = data1['187GHzH'][index]
            self.tb23H    = data1['238GHzH'][index]
            self.tb36H    = data1['365GHzH'][index]
            self.tb89H    = data1['890GHzH'][index]
            
            ## Gradient ratio
            self.GR_69_187_V  = (self.tb18V-self.tb6V)/(self.tb18V+self.tb6V)
            self.GR_187_365_V = (self.tb36V-self.tb18V)/(self.tb36V+self.tb18V)
            self.sigma40      = data1['sigma_40'][index]
            
            self.TbV_AMSR = [[tb6V, tb10V, tb18V, tb23V, tb36V] for tb6V, tb10V, tb18V, tb23V, tb36V
                             in zip(self.tb6V, self.tb10V, self.tb18V, self.tb23V, self.tb36V)]
            self.TbH_AMSR = [[tb6H, tb10H, tb18H, tb23H, tb36H] for tb6H, tb10H, tb18H, tb23H, tb36H
                             in zip(self.tb6H, self.tb10H, self.tb18H, self.tb23H, self.tb36H)]
        elif empty == True:

            # self.subset   = []
            self.sd       = []
            self.sit      = []
            self.index    = []
            self.lat      = []
            self.lon      = []
            self.date     = []
            self.id       = []
            self.sd_std   = []
            self.sit_std  = []
            self.slope     = []
            self.slope_std = []
            self.t2m      = []
            self.skt      = []
            self.tcwv     = []
            ## AMSR2 V
            self.tb6V     = []
            self.tb7V     = []
            self.tb10V    = []
            self.tb18V    = []
            self.tb23V    = []
            self.tb36V    = []
            self.tb89V    = []
            
            ## AMSR2 H
            self.tb6H     = []
            self.tb7H     = []
            self.tb10H    = []
            self.tb18H    = []
            self.tb23H    = []
            self.tb36H    = []
            self.tb89H    = []
            
            ## Gradient ratio
            self.GR_69_187_V  = []
            self.GR_187_365_V = []
            self.sigma40      = []
            
            self.TbV_AMSR = []
            self.TbH_AMSR = []
            
    def concatenate(self, obj2):
        " Concatenate two class objects e.g. combine to one"
        # update directory of obj1 with informaiton from obj2
        # self.__dict__.update(obj2.__dict__)
        
        for dict1, dict2 in zip(self.__dict__.items(), obj2.__dict__.items()):
            # dict[0] = keys e.g. class names
            # dict[1] = values in keys e.g. values of a given key
            setattr(self, dict1[0], np.append(dict1[1],dict2[1]))

    def plot(self, varname):

        if varname in self.__dict__.keys():
            var=self.__dict__.get(varname)
            title = ('OIB distribution of ' + varname.upper())
            if self.subset == "FYI":
                # clim = [1, 2]
                clim = [240,260]
            elif self.subset == "MYI":
                clim = [2, 6]
            else: # snow values
                clim = [0, 1]
            po(self.lat, self.lon, var, title + " " + self.subset, ylabel=varname.upper() + ' [m]', clim=clim, s=10)
        else:
            print('ERROR variable in in class')
            print('Class attributes are:', self.__dict__.keys())

def Get_OIB_data(directory_path, plot=False):

    # directory_path = os.path.dirname(os.path.dirname(os.getcwd()))
    sourcedir = directory_path + "/RRDP_v2.2/NERSC_OIB/2013-2019/*.text"
    file = "QSCAT-vs-SMAP-vs-SMOS-vs-ASCAT-vs-AMSR2-vs-ERA-vs-NERSC-OIB-20150403.TEXT"
    ifile = os.path.join(directory_path, "RRDP_v2.2/NERSC_OIB/2013-2019/" + file)
    dtype = float
    files = os.listdir(directory_path +  "/RRDP_v2.2/NERSC_OIB/2013-2019/")
    ## make empty class object
    d_final_FYI = data(ifile, dtype, [])
    d_final_MYI = data(ifile, dtype, [])
    
    ## loop throguh files
    for file in files:
        filename = os.path.join(directory_path, "RRDP_v2.2/NERSC_OIB/2013-2019/" + file)
        data_np = np.genfromtxt(filename, skip_header=1, delimiter=',', encoding='utf8')
        if "-NERSC-" in filename:
            dtype = [float if (e!=2 and e!=3 and e!=21 and e!=22 and e!=23) else object for e in list(range(data_np.shape[1]))]
        else:
            dtype = [float if (e!=2 and e!=3) else object for e in list(range(data_np.shape[1]))]
    
        d = data(filename, dtype)            
        #% FYI/MYI based on scatter  
        ## exclude the 10% of measurements with the heighest std
        exclude_points_SD = int(len(d.sd_std)/10)
        exclude_points_SIT = int(len(d.sit_std)/10)
        
        ## argpartition (orders based on size)
        ind_SD = np.argpartition(d.sd_std, -exclude_points_SD)[-exclude_points_SD:]
        ind_SIT = np.argpartition(d.sit_std, -exclude_points_SIT)[-exclude_points_SIT:]
        
        FYI = [];
        MYI = [];
        
        for i in range(len(d.GR_69_187_V)):
            if d.GR_69_187_V[i]<0.015 and d.GR_187_365_V[i]>-0.03 and d.sigma40[i]<-16.5 and d.sit[i]<2.0 and any(i!=ind_SD) and any(i!=ind_SIT):
                FYI.append(i)
            elif d.GR_69_187_V[i]<0 and d.GR_187_365_V[i]<-0.06 and d.sigma40[i]>-12 and d.sit[i]>2 and any(i!=ind_SD) and any(i!=ind_SIT):
                MYI.append(i)
        
        d_FYI = data(filename, dtype, FYI)
        d_MYI = data(filename, dtype, MYI)
    
        ## Make final class object with FYI and MYI data
        d_final_FYI.concatenate(d_FYI)
        d_final_MYI.concatenate(d_MYI)
        
        d_final_FYI.subset = "FYI"
        d_final_MYI.subset = "MYI"
    if plot==True:
        d_final_FYI.plot('sit')
        d_final_MYI.plot('sit')
    return [d_final_FYI, d_final_MYI]


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

    
    # ## Ocean temperature
    # T_w = (273.15-1.8)  # K
    # ks = 0.25 # Thermal conductivity of snow [w/(m*K)]
    
    # for sit, sd, skt in zip(self.sit, self.sd, self.skt):
    #     ## snow ice interface temperature
    #     T_sii = (self.k * sd * T_w + ks * sit *
    #                 skt) / (self.k * sd + ks * sit)
        
        
    #     ## make temperature profile of snow and ice and concatenate
    #     self.Tp_ice = np.linspace(T_sii, T_w, self.nl_ice)
    #     ## make one extra layer to avoid same temprature in lowest snow layer and ice
    #     self.Tp_snow.append([np.linspace(skt, T_sii,  self.nl_snow+2)[1:-1]])
        
    #     self.T_sii.append(T_sii)

    

        
    ## Ocean temperature
    T_w = (273.15-1.8)  # K
    ks = 0.25 # Thermal conductivity of snow [w/(m*K)]
    
    from sklearn.linear_model import LinearRegression
    for sit, sd, skt in zip(self.sit, self.sd, self.skt):
        ## snow ice interface temperature
        T_sii = (self.k * sd * T_w + ks * sit *
                    skt) / (self.k * sd + ks * sit)
        self.T_sii.append(T_sii)

        x = [0, sd]
        # print(x)
        y = [skt, T_sii]
        model = LinearRegression().fit(np.array(x).reshape(-1,1), np.array(y))
        
        if self.nl_snow == 2:
            sd0 = sd-0.03
            sd1 = 0.03
            xx = np.array([sd0/2, sd0 + sd1/2]).reshape(-1,1)
        else: # snow thinner than 3cm one layer
            xx = np.array([sd/2]).reshape(-1,1)
        Tsnow = model.predict(xx)
        # print(Tsnow)
        self.Tp_snow.append(Tsnow)
        
    if plot==True:
        # ice/snow thickness arrays for showing temperature profiles    
        ice_thick = np.linspace(-self.sit, 0, self.nl_ice + 1)
        snow_thick = np.linspace(0, self.sd, self.nl_snow + 1)
        temp_profile_snow = np.linspace(self.T_sii, self.skt, self.nl_snow + 1)
        
        plt.figure()
        plt.plot(temp_profile_snow, snow_thick,'.-')
        # plt.plot(np.array([np.array(self.Tp_ice).flatten(), self.T_sii]), ice_thick,'.-')
        plt.ylim([-0.1,0.1])
        plt.grid()
        plt.xlabel('Temperature [K]')
        plt.ylabel('Depth [m]')

def make_hist(var, varname, name='OIB'):
    try:
        var = np.concatenate(var)
    except:
        pass
    avgvar=np.nanmean(np.array(var))
    stdvar=np.nanstd(np.array(var))
    fig, ax= plt.subplots()
    if varname=='SIT':
        unit = ' [m]'
        plt.xlim(xmax=7,xmin=0)
        bins = np.linspace(-10, 10, 100)
    elif varname=='SD':
        unit = ' [m]'
        plt.xlim(xmax=1,xmin=0)
        bins = np.linspace(-1, 1, 100)
    elif varname=='rho':
        unit = ' [kg/m^3]'
        plt.xlim(xmax=400,xmin=100)
        bins = np.linspace(100, 400, 50)
    elif varname=='lex':
        unit = ' [mm]'
        plt.xlim(xmax=0.5,xmin=0)
        bins = np.linspace(-1, 1, 100)
    elif varname=='sal':
        unit = ' [PSU]'
        plt.xlim(xmax=20,xmin=0)
        bins = np.linspace(0, 20, 50)
    elif varname=='slope':
        unit = ' '
        plt.xlim(xmax=1,xmin=0)
        bins = np.linspace(-1, 1, 100)    
    ax.hist(var, bins=bins, color='red', alpha=0.5, rwidth=0.85,label=varname + ' avg:'+str(np.round(avgvar,2))+' ,std of: '+str(np.round(stdvar,2)))
    plt.grid(axis='y',alpha=0.75)
    # plt.ylim(ymax=650,ymin=0)
    if varname in ['sal', 'lex', 'rho']:
        plt.title('Simulation values of ' + varname ,fontsize=14,fontweight='bold')
    else:
        plt.title('Distribution of '+ name +' ' + varname + ' measurements',fontsize=14,fontweight='bold')
    plt.xlabel(varname + unit,fontsize=12)
    plt.legend(prop={"size":12}, bbox_to_anchor=(1.0, -0.15))
    plt.ylabel('Count',fontsize=14)
    plt.savefig(os.path.join(pp, 'OE_figures/hist_' + varname + '.png'), bbox_inches='tight')
    plt.show()
    
## Make simulated data
class SMRT_simulations(data):
    
    # init method or constructor
    def __init__(self, icetype, snowpit_data, directory_path):
        
        # directory_path = os.path.dirname(os.path.dirname(os.getcwd()))
        # print(directory_path)
        file = "QSCAT-vs-SMAP-vs-SMOS-vs-ASCAT-vs-AMSR2-vs-ERA-vs-NERSC-OIB-20150403.TEXT"
        ifile = os.path.join(directory_path, "RRDP_v2.2/NERSC_OIB/2013-2019/" + file)
        super().__init__(ifile, float, empty=True)
        self.directory_path     = directory_path
        self.subset             = icetype
        self.nl_snow            = 2
        self.nl_ice             = 1
        self.T_sii              = []
        self.Tp_snow            = []
        self.Tp_ice             = []
        
        self.TbH_AMSR = np.array(self.TbH_AMSR)
        self.TbV_AMSR = np.array(self.TbV_AMSR)

        ## layer specific values
        # choice - two snow layers (one thick and one of 3 cm)
        # bottom layer has high salinity!
        
        ## start with making dataset of 10 values
        self.sd_sim             = []
        ## ice layer
        self.sit_sim             = [] # m
        self.sit_sal_sim         = [7*PSU]
        
        ## Get CB values
        self.SD_CB = np.concatenate([pit.SD for pit in snowpit_data])
        self.lex_CB = np.concatenate([pit.corr for pit in snowpit_data])
        self.rho_CB = np.concatenate([pit.rho for pit in snowpit_data])
        self.sal_CB = np.concatenate([pit.sal for pit in snowpit_data])
        self.SD_CB_tot = np.concatenate([pit.SD_tot for pit in snowpit_data])
        
        ## extract values to be a part of top and bottom indexes
        self.sd_rho_sim = []
        self.sd_lex_sim = []
        self.sd_sal_sim = []
        
    def extract_simulation_values(self, num_simulations=10, plot=True):
        for count, sd in zip(range(num_simulations), self.sd[:num_simulations]):
            ## only extract the number of simulations we want
            # which profiles to include - those within a 10 cm interval of the true depth
            index1 = [self.SD_CB_tot <= sd+0.05]
            index2 = [self.SD_CB_tot >= sd-0.05]
            self.indexx = np.logical_and(index1, index2)[0]
            
            ## calculate cummulated sum of each profile
            self.cumsum = np.concatenate([np.cumsum(sd_cb) for sd_cb in self.SD_CB[self.indexx]])
            
            ## include corr/dens of all values above 3cm of sd in top layer
            self.indexx_t = [cs<=(sd-0.03) for cs in self.cumsum]
            ## include corr/dens of all values below 3cm of snowdepth in bottom layer
            self.indexx_b = [cs>=(sd-0.03) for cs in self.cumsum]
            
            # get correlation lengths
            self.lex_top = np.concatenate(self.lex_CB[self.indexx])[self.indexx_t]
            self.lex_bottom = np.concatenate(self.lex_CB[self.indexx])[self.indexx_b]

            
            # get densities
            self.rho_top = np.concatenate(self.rho_CB[self.indexx])[self.indexx_t]
            self.rho_bottom = np.concatenate(self.rho_CB[self.indexx])[self.indexx_b]       
            
            # get salinities - pick all above lowest value in included profiles
            # for top layer
            self.sal_top = np.concatenate([sal[:-1] for sal in self.sal_CB[self.indexx]])
            self.sal_bottom = np.concatenate([sal[-1:] for sal in self.sal_CB[self.indexx]])
            self.cumsum_sal_t = np.concatenate([np.cumsum(sd_cb)[:-1] for sd_cb in self.SD_CB[self.indexx]])
            self.cumsum_sal_b = np.concatenate([np.cumsum(sd_cb)[-1:] for sd_cb in self.SD_CB[self.indexx]])
            
            # append mean values to simulations
            if self.nl_snow==2:
                self.sd_lex_sim.append([np.mean(self.lex_top), np.mean(self.lex_bottom)])
                self.sd_rho_sim.append([np.mean(self.rho_top), np.mean(self.rho_bottom)])
                self.sd_sal_sim.append([np.mean(self.sal_top), np.mean(self.sal_bottom)])
            elif self.nl_snow > 2:
                self.sd_lex_sim.append(np.concatenate([[np.mean(self.lex_top)]*(self.nl_snow-1), [np.mean(self.lex_bottom)]]))
                self.sd_rho_sim.append(np.concatenate([[np.mean(self.rho_top)]*(self.nl_snow-1), [np.mean(self.rho_bottom)]]))
                self.sd_sal_sim.append(np.concatenate([[np.mean(self.sal_top)]*(self.nl_snow-1), [np.mean(self.sal_bottom)]]))
            
            # plot for each 30th value
            if (count % 30) == 0 and (plot==True):
                plot_simulation_parameters_from_CB(self, sd, count)

        ## distribution based on CB case 0.12-0.45
        # slight tendency to lower lex in top and higher in bottom
        # self.sd_lex_sim         = [0.15e-4] * self.nl_snow
        # ## CB distribution top 0.5-2 PSU, bottom 4-13 PSU
        # self.sd_sal_sim         = [1*PSU, 10*PSU]
        # self.sd_rho_sim         = [300] * self.nl_snow
        
  
        # Results of simulation
        # self.date                = []
        # self.lat                 = []
        # self.lon                 = []
        self.TbV                 = []
        self.TbH                 = []
        self.TbV_AMSR_all        = []
        self.TbH_AMSR_all        = []
        self.res                 = []
        self.sd_all              = []
        self.sd_rho_all          = []
        self.sd_lex_all          = []
        self.sd_sal_all          = []
        self.Tsnow_all           = []
        self.Tice_all            = []
        self.sit_all             = []
        self.slope_all           = []
        
        # print(self.TbV_AMSR[0].shape)
        
    def set_vars(self, plot=False):
        ## GET subset of data
        [d_final_FYI, d_final_MYI] = Get_OIB_data(self.directory_path)

        if self.subset == "FYI":
            # print(d_final_FYI.sit)
            self.concatenate(d_final_FYI)
            self.ice_type       = "firstyear"
            self.k              = 2.10
            self.plot('sit')
            self.TbV_AMSR = self.TbV_AMSR.reshape((-1, 5))
            self.TbH_AMSR = self.TbH_AMSR.reshape((-1, 5))
            # print('mean is:', np.mean(self.TbV_AMSR, axis=0))
            # print('mean is:', np.mean(self.TbH_AMSR, axis=0))
            # print('mean is:', np.mean(self.sd, axis=0))
            # print('mean is:', np.mean(self.sit, axis=0))
            # print(self.sit, self.sd, self.t2m)
            
        elif self.subset == "MYI":
            self.concatenate(d_final_MYI)
            self.ice_type       = "multiyear"
            self.k              = 1.88
            self.plot('sit')
            self.TbV_AMSR = self.TbV_AMSR.reshape((-1, 5))
            self.TbH_AMSR = self.TbH_AMSR.reshape((-1, 5))
            print('mean is:', np.mean(self.TbV_AMSR, axis=0))
            print('mean is:', np.mean(self.TbH_AMSR, axis=0))
        ## fill in temperature information
        # two piece linear temperature profile
        temp_profile(self, plot=plot)
        
    
    def simulated_data(self, num_simulations=10, save=True, plot=True, rough=False):
        
        ## get coor length, density and salinity from CB data
        
        # # create the model/sensor
        m = make_model("iba", "dort")
        
        frequencies=[6.925e9,10.65e9,18.7e9,23.8e9,36.5e9];  # [Hz]
        beamorientation = 55. # [degrees]
        sensor = sensor_list.passive(np.array(frequencies), beamorientation)
        ## make simulations
        # for lav coor lengths, density og salinity som der løbes over
        # print(np.mean(self.TbV_AMSR[0:10], axis=0))
        for num in range(num_simulations):
            # try:
                #print(self.sd_lex_sim[num])
                # self.sd_lex_sim[num] = [0.2*1e-3, 0.7*1e-3]
                #print(self.sd_lex_sim[num])
                nl_snow = self.nl_snow
                if self.sd[num] <= 0.03:
                    print(str(1))
                    nl_snow = 1 # only one snow layer
                    sd_sim = [self.sd[num]]
                    self.sd_sal_sim[num] = [self.sd_sal_sim[num][-1]]
                    self.sd_rho_sim[num] = [self.sd_rho_sim[num][-1]]
                    self.sd_lex_sim[num] = [self.sd_lex_sim[num][-1]]
                    
                    Tsnow = [self.Tp_snow[num][-1]]

                elif self.sd[num] <= 0.06 and nl_snow>2:
                    nl_snow = 2
                    sd_sim = [self.sd[num]-0.03/nl_snow]*(nl_snow-1) +  [0.03]
                    self.sd_sal_sim[num] = self.sd_sal_sim[num][:nl_snow]
                    if np.isnan(self.sd_lex_sim[num]).any():
                        self.sd_lex_sim[num] = [self.sd_lex_sim[num][1]]*nl_snow
                        self.sd_rho_sim[num] = [self.sd_rho_sim[num][1]]*nl_snow
                    else:
                        self.sd_lex_sim[num] = self.sd_lex_sim[num][:nl_snow]
                        self.sd_rho_sim[num] = self.sd_rho_sim[num][:nl_snow]
                    
    
                else:
                    # sd_sim = [self.sd[num]-0.03, 0.03]
                    sd_sim = [self.sd[num]-0.03] +  [0.03]
                    self.sd_sal_sim[num] = self.sd_sal_sim[num][:nl_snow]
                    if np.isnan(self.sd_lex_sim[num]).any():
                        self.sd_lex_sim[num] = [self.sd_lex_sim[num][1]]*nl_snow
                        self.sd_rho_sim[num] = [self.sd_rho_sim[num][1]]*nl_snow
                    
                
                Tsnow = self.Tp_snow[num]
                Tice = self.T_sii[num]
                
                sit = self.sit[num]
    
                # print(sd_sim)
                #print(Tsnow)
                #print(Tice)
                # print([sd_lex*1e3 for sd_lex in self.sd_lex_sim[num]])
                # print((self.sd_rho_sim[num]))
                # print([sd_sal*1e3 for sd_sal in self.sd_sal_sim[num]])
                # self.sd_sal_sim[num] = np.array([0.0, 0.0])
                ## snowpack
                snowpack= make_snowpack(thickness=sd_sim,
                                             microstructure_model="exponential",
                                             ice_permittivity_model=ssp,
                                             density=self.sd_rho_sim[num], 
                                             corr_length=self.sd_lex_sim[num], 
                                             temperature = Tsnow,
                                             salinity=self.sd_sal_sim[num])
        
                ## Make icepack        
                sea_ice_density = 900
                icepack = make_ice_column(ice_type='firstyear',
                        thickness=[sit], temperature=Tice, 
                        microstructure_model='exponential',
                        corr_length=0.15e-3,
                        brine_inclusion_shape='spheres',
                        density=sea_ice_density,
                        salinity=7*PSU,
                        add_water_substrate="ocean")
                
                medium = snowpack + icepack
                
                #run simulations
                res = m.run(sensor, medium)
    
                if rough==True:
                    #% adjust for atmospheric influence
                    # print(self.tcwv[num],res,Tice)
                    [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(self.tcwv[num],res,Tice, m=0.08)
                ## append final values
                self.sd_all.append(sd_sim)
                self.sd_rho_all.append(self.sd_rho_sim[num][:nl_snow])
                self.sd_lex_all.append(self.sd_lex_sim[num][:nl_snow])
                self.sd_sal_all.append(self.sd_sal_sim[num][:nl_snow])
                self.Tsnow_all.append(Tsnow)
                self.Tice_all.append(Tice)
                self.sit_all.append(sit)
                self.slope_all.append(self.slope[num])
                ## append output
                self.TbV.append(TB_AMSR_V)
                self.TbH.append(TB_AMSR_H)
                self.res.append(res)
                self.TbV_AMSR_all.append(self.TbV_AMSR[num])
                self.TbH_AMSR_all.append(self.TbH_AMSR[num])

            # except:
            #     pass
            
        ## save simulated data to .npy file
        if save==True:
            datapath = "C:/Users/Ida Olsen/Documents/Speciale_2022/Optimal_Estimation/OE_ida/Simulated_data"
            with open(os.path.join(datapath,'TB_simulated.npy'), 'wb') as f:
                np.save(f, np.array(self.TbV))
                np.save(f, np.array(self.TbH))
                np.save(f, np.array(self.sd_all))
                np.save(f, np.array(self.sit_all))
                np.save(f, np.array(self.sd_rho_all))
                np.save(f, np.array(self.sd_lex_all))
                np.save(f, np.array(self.sd_sal_all))
                np.save(f, np.array(self.Tsnow_all))
                np.save(f,  np.array(self.t2m))
                np.save(f,  np.array(self.skt))
                np.save(f,  np.array(self.tcwv))
                np.save(f, np.array(self.Tice_all))
                np.save(f, np.array(self.slope_all))
                np.save(f, np.array(self.TbV_AMSR_all))
                np.save(f, np.array(self.TbH_AMSR_all))
                np.save(f, np.array(self.lat[:num]))
                np.save(f, np.array(self.lon[:num]))
                np.save(f, self.date[:num])
        if plot==True:
            freqs = [f*1e-9 for f in frequencies]
            y = np.mean(self.TbV, axis=0)
            err = np.std(self.TbV, axis=0)
            plt.plot(freqs, y, c='orange',label = 'TBV rough')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='r', facecolor='r', alpha=0.3)
            
            y = np.mean(self.TbH, axis=0)
            err = np.std(self.TbH, axis=0)
            plt.plot(freqs,y, c='cyan', label = 'TBH rough')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='b', facecolor='b', alpha=0.3)
            
            y = np.mean(self.TbV_AMSR_all, axis=0)
            err = np.std(self.TbV_AMSR_all, axis=0)
            plt.plot(freqs, y, c='r',label = 'TBV  AMSR')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='r', facecolor='r', alpha=0.3)
            
            y = np.mean(self.TbH_AMSR_all, axis=0)
            err = np.std(self.TbH_AMSR_all, axis=0)
            plt.plot(freqs,y, c='b', label = 'TBH AMSR')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='b', facecolor='b', alpha=0.3)
            
            y = np.mean([r.TbV() for r in self.res], axis=0)
            err = np.std([r.TbV() for r in self.res], axis=0)
            plt.plot(freqs, y, c='k', label = 'TBV')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='k', facecolor='k', alpha=0.3)
            
            y = np.mean([r.TbH() for r in self.res], axis=0)
            err = np.std([r.TbH() for r in self.res], axis=0)
            plt.plot(freqs, y, c='green', label = 'TBH')
            plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='g', facecolor='g', alpha=0.3)
            plt.grid()
            plt.ylim([200,260])
            plt.legend()
            plt.xlabel('Freq GHz')
            plt.ylabel('TB [K]')
            plt.title('Simulated FYI TBs from OIB data 2013-2019')
            plt.savefig(os.path.join(pp, 'OE_figures/Simulated_data.png'), bbox_inches='tight')
            plt.show()

            
            print(self.sd_all)
            make_hist(np.sum(self.sd_all, axis=1), 'SD')
            make_hist(self.slope_all, 'SD')
            make_hist(self.sd_rho_all, 'rho')
            make_hist(np.array(self.sd_sal_all)*1e3, 'sal')
            make_hist(np.array(self.sd_lex_all)*1e3, 'lex')
            make_hist(self.sit_all, 'SIT')
            make_hist(self.slope_all, 'slope')

snowpit_data = fill_data()
directory_path = os.path.dirname(os.path.dirname(os.getcwd()))
SIM_FYI = SMRT_simulations("FYI", snowpit_data=snowpit_data, directory_path=directory_path)
SIM_FYI.set_vars(plot=False)
print('mean is:', SIM_FYI.sd[0])
print('mean is:', np.mean(SIM_FYI.TbV_AMSR, axis=0))
print('mean is:', np.mean(SIM_FYI.TbH_AMSR, axis=0))
ns = len(SIM_FYI.sd)
SIM_FYI.extract_simulation_values(num_simulations=ns, plot=True)
SIM_FYI.simulated_data(num_simulations=ns, save=False, rough=True)
#SIM_FYI.lat = SIM_FYI.lat[:ns]
#SIM_FYI.lon = SIM_FYI.lon[:ns]
#SIM_FYI.plot('sit_all')

def plot_output(self):
    plt.figure(figsize=(7,6))
    frequencies=[6.925e9,10.65e9,18.7e9,23.8e9,36.5e9];
    freqs = [f*1e-9 for f in frequencies]
    # plt.figure(figsize=(6,6))
    # for TbV, TbH, TbV_AMSR, TbH_AMSR in zip(self.TbV, self.TbH, self.TbV_AMSR_all, self.TbH_AMSR_all):
    #     plt.scatter(freqs, TbV, c='r')
    #     plt.scatter(freqs, TbH, c='b')
    #     plt.scatter(freqs, TbV_AMSR, c='k')
    #     plt.scatter(freqs, TbH_AMSR, c='green')
    # print(np.mean(TbV, axis=0).shape)
    # print(np.mean(TbV_AMSR, axis=0).shape)
    # y = np.mean(self.TbV, axis=0)
    for TBV in self.TbV:
        plt.plot(freqs, TBV, c='salmon')
    for TBH in self.TbH:
        plt.plot(freqs, TBH, c='lightblue')
    
    y = np.mean(self.TbH, axis=0)
    err = np.std(self.TbH, axis=0)
    plt.plot(freqs,y, c='b', label = 'TBH')
    plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='b', facecolor='b', alpha=0.3)
    
    y = np.mean(self.TbV, axis=0)
    err = np.std(self.TbV, axis=0)
    plt.plot(freqs, y, c='r',label = 'TBV')
    plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='r', facecolor='r', alpha=0.3)
    
    # y = np.mean(self.TbH, axis=0)
    # err = np.std(self.TbH, axis=0)
    # plt.plot(freqs,y, c='b', label = 'TBH')
    # plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='b', facecolor='b', alpha=0.3)

    
    # y = np.mean(self.TbV_AMSR_all, axis=0)
    # err = np.std(self.TbV_AMSR_all, axis=0)
    # plt.plot(freqs, y, c='k', label = 'TBV AMSR')
    # plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='k', facecolor='k', alpha=0.3)
    
    # y = np.mean(self.TbH_AMSR_all, axis=0)
    # err = np.std(self.TbH_AMSR_all, axis=0)
    
    # plt.plot(freqs, y, c='green', label = 'TBH AMSR')
    # plt.fill_between(freqs, y-err/2, y+err/2, edgecolor='g', facecolor='g', alpha=0.3)
    plt.grid()
    plt.ylim([210,270])
    # plt.legend(['TBV', 'TBH', 'TBV AMSR', 'TBH AMSR'])
    plt.legend()
    plt.xlabel('Frequencies [GHz]')
    plt.ylabel('TB [K]')
    plt.title('Simulated FYI TBs from OIB data 2013-2019')
    plt.savefig(os.path.join(pp, 'OE_figures/Simulated_data.png'), bbox_inches='tight')
    plt.show()

plot_output(SIM_FYI)