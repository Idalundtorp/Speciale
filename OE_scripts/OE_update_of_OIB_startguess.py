# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 23:24:32 2023

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
from sklearn.linear_model import LinearRegression
import os
import pyOptimalEstimation as pyOE
## self made modules
from roughness import roughness
  
class data():
    
    # init method or constructor
    def __init__(self, filename, dtype, ice_type=-1, empty=False):
        if empty == False:
            names = ["latitude","longitude","time","identification","sd_mean","sd_std","SDln","sit_mean","sit_std","SITln","FRB","FRBstd","Tsur","Tair","wsd","wdensity","lat1","lon1","time1","reference-id","69GHzH","69GHzV","73GHzH","73GHzV","107GHzH","107GHzV","187GHzH","187GHzV","238GHzH","238GHzV","365GHzH","365GHzV","890GHzH","890GHzV","Earth Incidence","Earth Azimuth","scanpos","upstreamfile","timediff","latitude2","longitude2","time2","reference-id2","upstreamfile","msl","u10","v10","ws","t2m","skt","istl1","istl2","istl3","istl4","sst","d2m","tcwv","tclw","tciw","ssrd","strd","e","tp","sf","fal","siconc","latitude3","longitude3","time3","reference-id3","upstreamfile","sigma_40","sigma_40_mask","nb_samples","warning","std","latitude4","longitude4","time4","reference-id4","upstreamfile","SMOS_Tbv","SMOS_Tbh","SMOS_RMSE_v","SMOS_RMSE_h","SMOS_nmp","SMOS_data1loss"]
            #print(filename)
            if "-NERSC-" in filename:
                print(filename)
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
            
            try:
                ## AMSR2 V
                self.tb6V     = data1['69V'][index]
                self.tb7V     = data1['73V'][index]
                self.tb10V    = data1['107V'][index]
                self.tb18V    = data1['187V'][index]
                self.tb23V    = data1['238V'][index]
                self.tb36V    = data1['365V'][index]
                self.tb89V    = data1['890V'][index]
                
                ## AMSR2 H
                self.tb6H     = data1['69H'][index]
                self.tb7H     = data1['73H'][index]
                self.tb10H    = data1['107H'][index]
                self.tb18H    = data1['187H'][index]
                self.tb23H    = data1['238H'][index]
                self.tb36H    = data1['365H'][index]
                self.tb89H    = data1['890H'][index]
            except:
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

def Get_OIB_data(directory_path, all_data=True,plot=False):

    # directory_path = os.path.dirname(os.path.dirname(os.getcwd()))
    if all_data==True:
        datafolder= "ALL_DATA" #"2013-2019"
    elif all_data=="old":
        datafolder= "2009-2015"
    elif all_data=="new":
        datafolder= "2016-2019"
    sourcedir = directory_path + "/RRDP_v2.2/NERSC_OIB/"+datafolder+"/*.text"
    file = "QSCAT-vs-SMAP-vs-SMOS-vs-ASCAT-vs-AMSR2-vs-ERA-vs-NERSC-OIB-20150403.TEXT"

    files = os.listdir(directory_path +  "/RRDP_v2.2/NERSC_OIB/"+datafolder+"/")
    ifile = os.path.join(directory_path, "RRDP_v2.2/NERSC_OIB/"+datafolder+"/" + files[0])
    dtype = float
    ## make empty class object
    d_final_FYI = data(ifile, dtype, [])
    d_final_MYI = data(ifile, dtype, [])
    
    ## loop throguh files
    for file in files:
        filename = os.path.join(directory_path, "RRDP_v2.2/NERSC_OIB/"+datafolder+"/" + file)
        data_np = np.genfromtxt(filename, skip_header=1, delimiter=',', encoding='utf8')
        if "-NERSC-" in filename:
            dtype = [float if (e!=2 and e!=3 and e!=21 and e!=22 and e!=23) else object for e in list(range(data_np.shape[1]))]
        else:
            dtype = [float if (e!=2 and e!=3) else object for e in list(range(data_np.shape[1]))]
        print(filename)
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


def SMRT_simulation(X, self, variables):
    
    ## Make icepack           
    sea_ice_density = 900
    icepack = make_ice_column(ice_type='firstyear',
            thickness=[self.sit], 
            temperature=[self.T_sii], 
            microstructure_model='exponential',
            corr_length=[self.sit_lex],
            brine_inclusion_shape='spheres',
            density=sea_ice_density,
            salinity=[self.sit_sal],
            add_water_substrate="ocean")
    if X[3]<=0:
        X[3]=0.1*PSU
    elif X[2]<=0:
        X[2]=0.1*PSU

    try:
        print('updating values')
        self.lex[0]=X[0]
        self.lex[1]=X[1]
        self.lex[2]=self.lex[1]
        print(self.lex)
        
        self.sal[1]=X[2]
        self.sal[2]=X[3]
        print(self.sal)
        
        self.rho[0]=X[4]
        self.rho[1]=X[5]
        self.rho[2]=self.rho[1]
        print(self.rho)
        
        self.calc_temp(plot=False)
        
        snowpack = make_snowpack(self.sd, 
                        microstructure_model='exponential',
                        ice_permittivity_model=ssp, 
                        density=self.rho, 
                        corr_length=self.lex, 
                        temperature = self.Tsnow,
                        salinity=self.sal)
        medium = snowpack + icepack
        
        
        m = make_model("iba", "dort")
        
        frequencies=np.array([6.925e9,10.65e9,18.7e9,23.8e9,36.5e9]);  # [Hz]
        
        beamorientation = 55. # [degrees]
        sensor = sensor_list.passive(frequencies, beamorientation)
        
        #% Run simulations
        res = m.run(sensor, medium)
    
        # #% adjust for atmospheric influence
        self.m = X[6]
        print(self.m)
    
        [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(self.tcwv,res,self.T_sii, m=self.m)
        
        plt.figure(figsize=(8,6))
        frequencies=[6.925e9,10.65e9,18.7e9,23.8e9,36.5e9];
        freqs = [f*1e-9 for f in frequencies]
        
    
        y = TB_AMSR_V
        plt.plot(freqs, y, c='b',label = 'TBV SMRT')
        y = TB_AMSR_H
        plt.plot(freqs, y, c='grey',label = 'TBH SMRT')
        
        y = np.array([np.mean(d_final_FYI.TbV_AMSR.reshape((-1,5)), axis=0)]).flatten()
        plt.plot(freqs, y, c='g',label = 'TBV OIB')
        y = np.array([np.mean(d_final_FYI.TbH_AMSR.reshape((-1,5)), axis=0)]).flatten()
        plt.plot(freqs,y, c='r', label = 'TBH OIB')
        plt.grid()
        plt.ylim([200,270])
        plt.legend(bbox_to_anchor=(1, 0.8))
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('TB [K]')
        plt.title('')
        plt.show()
        return np.concatenate([TB_AMSR_V,TB_AMSR_H])
    except:
        # TB_AMSR_V,TB_AMSR_H = 
        return [170]*10
        
    

class smrt_values():
    # init method or constructor
    def __init__(self, OIB_data, frac=0.25):
        #
        self.pol = 'both'
        self.Tsur = np.mean(OIB_data.skt)
        self.m = 0.19 # roughness variable
        self.tcwv = np.mean(OIB_data.tcwv) #total columnar water vapour
        self.ks = 0.25
        # Snow variables
        sd_tot = np.mean(OIB_data.sd)
        sd0 = sd_tot - sd_tot*frac
        sd1 = sd_tot*frac - 0.01
        sd2 = 0.01
        self.sd = [sd0, sd1, sd2]
        self.lex = [0.24e-3, 0.37e-3, 0.37e-3] ## mm - based on paper by Wisman and Metzler
        self.sal = [0.001*PSU, 0.5*PSU, 13*PSU]  # Based on CB field observations 
        self.rho = [350, 270, 270] # based on paper by Wisman and Metzler
        # Ice varibles
        self.sit = np.mean(OIB_data.sit)
        self.sit_lex = 0.15e-3
        self.sit_sal = 6*PSU
        # Temperatures
    def calc_temp(self, plot=True):
        k = 2.10 # thermal conductivity FYI
        T_w = (273.15-1.8)  # Ocean temperature [K]
        ks = self.ks # Thermal conductivity of snow [w/(m*K)]

        ## snow ice interface temperature
        self.T_sii = (k * np.sum(self.sd) * T_w + ks * self.sit *
                    self.Tsur) / (k * np.sum(self.sd) + ks * self.sit)
        x = [0, np.sum(self.sd)]
        y = [self.Tsur, self.T_sii]
        model = LinearRegression().fit(np.array(x).reshape(-1,1), np.array(y))
        xx = np.array([self.sd[0]/2, self.sd[0] + self.sd[1]/2, self.sd[0]+ self.sd[1] + self.sd[2]/2] ).reshape(-1,1)
        self.Tsnow = model.predict(xx)
        
        #if plot==True:
        # ice/snow thickness arrays for showing temperature profiles    
        ice_thick = np.array([0, -self.sit])
        t_ice = np.array([self.T_sii, T_w])

directory_path = os.path.dirname(os.path.dirname(os.getcwd()))
[d_final_FYI, d_final_MYI] = Get_OIB_data(directory_path, all_data=True,plot=False)
data = smrt_values(d_final_FYI)
data.calc_temp(plot=False)

x_ap = [0.24e-3, 0.37e-3, 0.01*PSU, 13*PSU, 350, 270, 0.19]
x_vars = ["lex0", "lex1", "sal1", "sal2", "rho0", "rho1", "slope"]
x_cov = np.diag([0.001e-3, 0.0001e-3, 1*PSU, 1*PSU, 50, 50, 0.1])

## Define name of y-variables
y_vars = ["6V", "10V", "18V", "23V", "36V", "6H", "10H", "18H", "23H", "36H"]
y_obs  = np.concatenate([np.mean(d_final_FYI.TbV_AMSR.reshape((-1,5)), axis=0)] + [np.mean(d_final_FYI.TbH_AMSR.reshape((-1,5)), axis=0)])
V_sens = np.array([2.0, 2.0, 0.5, 0.5, 0.5])**2
H_sens = np.array([2.0, 2.0, 1.0,1.0,2.0])**2
#y_cov = np.diag([V_sens**2]*5 + [H_sens**2]*5)
y_cov = np.diag(np.concatenate([V_sens, H_sens]))

forwardKwArgs = {"self": data, "variables": x_vars}
# define the optimal estimation
oe = pyOE.optimalEstimation(
        x_vars, x_ap, x_cov, y_vars, y_obs, y_cov, SMRT_simulation
        ,forwardKwArgs=forwardKwArgs
        )
# run the optimization and print outputs
oe.doRetrieval(maxIter=15)

Bias = np.round(y_obs- np.array(oe.y_op),3)
print(Bias)
