# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 09:55:04 2023

@author: Ida Olsen
"""

import netCDF4 as nc
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dt
## Home made modules     
from polar_plots_pred import plot_one as po
from polar_plots_pred import polar_plot as pp
# plot path
pp = "C:/Users/Ida Olsen/Documents/Speciale_2022/figures"
  
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
        datafolder= "2013-2019" #"2013-2019"
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

directory_path = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
[d_final_FYI, d_final_MYI] = Get_OIB_data(directory_path, all_data=True,plot=False)

#%%
def Time_series(self, var):
    plt.figure(figsize=(7,7))
    #dates = [dt.strptime(d.decode('utf8'), '%Y-%m-%dT%H:%M:%SZ') for d in self.date]
    try:
        dates = [d.decode('utf8') for d in self.date] #[dt.strptime(d.decode('utf8')# [dt.date(dt.strptime(d.decode('utf8'), '%Y-%m-%dT%H:%M:%S')) for d in self.date]
    except:
        dates = [d[:-1].decode('utf8') for d in self.date] #[dt.strptime(d.decode('utf8'), '%Y-%m-%dT%H:%M:%SZ') for d in self.date]
        print(dates)
    dates = [d[:10].decode('utf8') for d in self.date]
    plt.figure(figsize=(7, 6))
    plot=plt.scatter(dates, var)
    # plot=plt.scatter(self.date, self.sit[:-1])
    plt.xticks(np.unique(dates)[::2]) # dates[::10])
    plt.grid()
    
    plt.ylim([0.0,0.30])
    plt.xticks(rotation=45, ha='right')
    
    title = 'OIB Time series of SD [m] ' + dates[0][:4] +'-'+ dates[-1][:4]
    plt.title(title,fontsize=16,fontweight='bold')
    plt.ylabel('SD [m]')

    #plt.savefig(savename, bbox_inches='tight')
    plt.show()
    
Time_series(d_final_FYI, d_final_FYI.sd)

import os
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

#%%
def polar_plot(self, title, ylabel='SIT [m]', clim=False, savename=False, saveloc=False, s=8):    

    
    #index = [(TbH<(FYI_OIB_AMSR[5:]-2*FYI_OIB_AMSR_std[5:])).any() for TbH in FYI_data.TbH_AMSR_all][:-1]
    #index1 = [True for TbH in FYI_data.TbH_AMSR_all][:-1]
    plt.figure(figsize=(7, 7))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,65,90],ccrs.PlateCarree())
    # ax.set_extent([np.min(lon)-0.05,np.max(lon)+0.05,np.min(lat)-0.01,np.max(lat)+0.01],ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.OCEAN)        
    ax.add_feature(cfeature.LAND)
    ax.gridlines(draw_labels=False)
    for i, c in zip([index1, index], ['g', 'r']):

        plot=plt.scatter(self.lon[i], self.lat[i], c=c, s=s, transform=ccrs.PlateCarree(),)
    plt.title(title,fontsize=16,fontweight='bold')
    
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(plot, cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel(ylabel);
    cb.ax.tick_params(labelsize=12)
    plt.savefig(savename, bbox_inches='tight')
    plt.show()

def polar_plot_all(self, title, ylabel='SIT [m]', clim=False, savename=False, saveloc=False, s=8):    

    plt.figure(figsize=(7, 7))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,65,90],ccrs.PlateCarree())
    # ax.set_extent([np.min(lon)-0.05,np.max(lon)+0.05,np.min(lat)-0.01,np.max(lat)+0.01],ccrs.PlateCarree())
    ax.coastlines()
    # ax.add_feature(cfeature.OCEAN)        
    ax.add_feature(cfeature.LAND)
    ax.gridlines(draw_labels=False)
    plot=plt.scatter(self.lon, self.lat, c=self.sit, s=s, transform=ccrs.PlateCarree(),)
    plt.title(title,fontsize=16,fontweight='bold')
    
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(plot, cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel(ylabel);
    cb.ax.tick_params(labelsize=12)
    plt.savefig(savename, bbox_inches='tight')
    plt.show()

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
        plt.xlim(xmax=2.2,xmin=0)
        bins = np.linspace(-2, 2, 50)
    elif varname=='SD':
        unit = ' [m]'
        plt.xlim(xmax=0.4,xmin=0)
        bins = np.linspace(-1, 0.4, 100)   
    ax.hist(var, bins=bins, color='red', alpha=0.5, rwidth=0.85,label=varname + ' avg:'+str(np.round(avgvar,2))+' ,std of: '+str(np.round(stdvar,2)))
    plt.grid(axis='y',alpha=0.75)
    # plt.ylim(ymax=650,ymin=0)

    plt.title('Distribution of '+ name +' ' + varname + ' measurements',fontsize=14,fontweight='bold')
    plt.xlabel(varname + unit,fontsize=12)
    plt.legend(prop={"size":12}, bbox_to_anchor=(1.0, -0.15))
    plt.ylabel('Count',fontsize=14)
    plt.savefig(os.path.join(pp, 'OE_figures/hist_' + varname + '.png'), bbox_inches='tight')
    plt.show()  
    
make_hist(d_final_FYI.sd, 'SD', name='OIB')
make_hist(d_final_FYI.sit, 'SIT', name='OIB')

#polar_plot(FYI_data, 'OIB', savename='OIB_outliers')
polar_plot_all(d_final_FYI, 'OIB SIT [m] over FYI', savename='OIB_all')

from Get_Sat_Data import Get_AMSR_ECMWF
[AMSR, ECMWF] = Get_AMSR_ECMWF()


def compare_CB_AMSR_OIB_AMSR(OIB_AMSR2, AMSR2, index=[], index1=[]):
    

    plt.figure(figsize=(7,6))
    color = ['g', 'b', 'r', 'k']
    name = ['V AMSR2 avg. at OIB', 'V AMSR2 avg. at CB', 'H AMSR2 avg. at OIB', 'H AMSR2 avg. at CB']
    ADV = np.array([A.DataSetV for A in AMSR2])
    ADH = np.array([A.DataSetH for A in AMSR2])
    frequencies=[6.925e9,10.65e9,18.7e9,23.8e9,36.5e9];
    freqs = [f*1e-9 for f in frequencies]
    # for i,i1  in zip(index, index1):
    #     print(i1)
    plt.figure(figsize=(7,6))
        #for TBV, TBH in zip(d_final_FYI.TbV_AMSR.reshape(-1,5)[i:i1],d_final_FYI.TbH_AMSR.reshape(-1,5)[i:i1]):
    for TBV, TBH in zip(d_final_FYI.TbV_AMSR.reshape(-1,5),d_final_FYI.TbH_AMSR.reshape(-1,5)):
        plt.plot(freqs, TBV, 'palegreen', 'o',linewidth=0.5) # , label=d_final_FYI.date[i].decode('utf8'))
        plt.plot(freqs, TBH, 'lightsalmon', linewidth=0.5)
    for n, c, A in zip(name, color,[d_final_FYI.TbV_AMSR.reshape(-1,5),ADV,d_final_FYI.TbH_AMSR.reshape(-1,5),ADH]):
         y = np.mean(A, axis=0).flatten()
         err = np.std(A, axis=0).flatten()
         plt.plot(freqs, y, c=c, label = 'TB' + n)
    # labels = []
    labels = [str(i) for i in np.arange(200,261, step=10)]
    plt.yticks(np.arange(200,261, step=10),labels)
        # plt.fill_between(freqs, y-err/2, y+err/2, edgecolor=c, facecolor=c, alpha=0.3)
    # plt.legend(loc='center left',bbox_to_anchor=(1, 0.5),fontsize=12)
    plt.legend(loc='lower left',fontsize=12)
    plt.xlabel('Frequency [GHz]',fontsize=14)
    plt.ylabel('TB [K]',fontsize=14)
    plt.ylim([200,261])
    plt.xlim([6, 37])
    #plt.grid()
    title = 'Comparison of AMSR2 TBs at CB and at OIB FYI'
    plt.title(title,fontsize=16)
    savedic = 'C:/Users/Ida Olsen/Documents/Speciale_2022/Cambridge_Bay/'
    plt.savefig(savedic + title + '.png', bbox_inches='tight')
    plt.show()

dates = [d.decode('utf8') for d in d_final_FYI.date]
years = [d[:4] for d in dates]
index = [i for i in range(len(years)-1) if years[i+1]!=years[i]]
index1 = index[1:] + [len(dates)]

compare_CB_AMSR_OIB_AMSR(d_final_FYI, AMSR, index, index1)