# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 11:06:45 2022

@author: Ida Olsen

Extract relevant AMSR2 subset
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from Get_Sat_Data import Get_AMSR_ECMWF
from prepare_inputs_OE_AMSR2 import OIB_data

[AMSR, ECMWF] = Get_AMSR_ECMWF()

datapath = "C:/Users/Ida Olsen/Documents/Speciale_2022/Optimal_Estimation/OE_ida/Simulated_data/TB_simulated.npy"
# plot path
pp = "C:/Users/Ida Olsen/Documents/Speciale_2022/figures"

with open(datapath, 'rb') as f:
            TbV_sim = np.load(f)
            TbH_sim = np.load(f)
            sd = np.load(f)
            sit = np.load(f)
            sd_rho = np.load(f)
            sd_lex = np.load(f)
            sd_sal = np.load(f)
            Tsnow = np.load(f)
            t2m = np.load(f)
            skt = np.load(f)
            tcwv = np.load(f)
            Tice = np.load(f)
            OIB_surf_rough = np.load(f)
            TbV_AMSR_all = np.load(f)
            TbH_AMSR_all = np.load(f)
            
FYI_data = OIB_data('FYI', new=True)
FYI_OIB_AMSR = np.concatenate([np.mean(FYI_data.TbV_AMSR_all, axis=0), np.mean(FYI_data.TbH_AMSR_all, axis=0)])
FYI_OIB_AMSR_std = np.concatenate([np.std(FYI_data.TbV_AMSR_all, axis=0), np.std(FYI_data.TbH_AMSR_all, axis=0)])


def compare_CB_AMSR_OIB_AMSR(OIB_AMSR2, AMSR2):
    plt.figure(figsize=(7,6))
    color = ['g', 'b', 'r', 'k']
    name = ['V OIB', 'V CB', 'H OIB', 'H CB']
    ADV = np.array([A.DataSetV for A in AMSR2])
    ADH = np.array([A.DataSetH for A in AMSR2])
    frequencies=[6.925e9,10.65e9,18.7e9,23.8e9,36.5e9];
    freqs = [f*1e-9 for f in frequencies]
    for TBV, TBH in zip(FYI_data.TbV_AMSR_all,FYI_data.TbH_AMSR_all):
            plt.plot(freqs, TBV, 'palegreen', 'o',linewidth=0.5)
            plt.plot(freqs, TBH, 'lightsalmon', linewidth=0.5)
    for n, c, A in zip(name, color,[FYI_data.TbV_AMSR_all,ADV,FYI_data.TbH_AMSR_all,ADH]):
        y = np.mean(A, axis=0).flatten()
        err = np.std(A, axis=0).flatten()
        plt.plot(freqs, y, c=c, label = 'TB' + n)
        # plt.fill_between(freqs, y-err/2, y+err/2, edgecolor=c, facecolor=c, alpha=0.3)
    # plt.legend(loc='center left',bbox_to_anchor=(1, 0.5),fontsize=12)
    plt.legend(loc='lower left',fontsize=12)
    plt.xlabel('Frequency [GHz]',fontsize=14)
    plt.ylabel('TB [K]',fontsize=14)
    plt.ylim([200,261])
    plt.xlim([6, 37])
    plt.grid()
    title = 'Comparison of average AMSR2 TBs at CB and at OIB FYI'
    plt.title(title,fontsize=16)
    savedic = 'C:/Users/Ida Olsen/Documents/Speciale_2022/Cambridge_Bay/'
    plt.savefig(savedic + title + '.png', bbox_inches='tight')
    

    
    
compare_CB_AMSR_OIB_AMSR(FYI_data, AMSR)


import os
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

def polar_plot(self, title, ylabel='SIT [m]', clim=False, savename=False, saveloc=False, s=8):    

    
    index = [(TbH<(FYI_OIB_AMSR[5:]-2*FYI_OIB_AMSR_std[5:])).any() for TbH in FYI_data.TbH_AMSR_all][:-1]
    index1 = [True for TbH in FYI_data.TbH_AMSR_all][:-1]
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
    ax.add_feature(cfeature.OCEAN)        
    ax.add_feature(cfeature.LAND)
    ax.gridlines(draw_labels=False)
    plot=plt.scatter(self.lon, self.lat, c=np.mean(self.sd, axis=1)[:-1], s=s, transform=ccrs.PlateCarree(),)
    plt.title(title,fontsize=16,fontweight='bold')
    
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(plot, cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel(ylabel);
    cb.ax.tick_params(labelsize=12)
    plt.savefig(savename, bbox_inches='tight')
    plt.show()
    
#polar_plot(FYI_data, 'OIB', savename='OIB_outliers')
#polar_plot_all(FYI_data, 'OIB', savename='OIB_all')

#%%
def Time_series(self, title):    

    index = 
    plt.figure(figsize=(7, 6))
    plot=plt.scatter(self.date, np.sum(self.sd, axis=1)[:-1])
    # plot=plt.scatter(self.date, self.sit[:-1])
    plt.xticks(self.date[::10])
    plt.grid()
    plt.ylim([0,0.30])
    plt.xticks(rotation=45, ha='right')
    plt.title(title,fontsize=16,fontweight='bold')

    #plt.savefig(savename, bbox_inches='tight')
    plt.show()
    
Time_series(FYI_data, 'time series')