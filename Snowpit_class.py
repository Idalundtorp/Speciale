# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 15:52:58 2022

@author: Ida Olsen

Module for the snowpit class
"""
# python packages
import numpy as np
import matplotlib.pyplot as plt


def plot_vars(var, var_name, unit, site=None):
    avg = np.mean(var)
    std = np.std(var)
    fig, ax= plt.subplots()
    ax.hist(x=var, color='red', alpha=0.5, rwidth=0.85,label= var_name +' avg:'+str(np.round(avg,2))+' and std of: '+str(np.round(std,2)))
    plt.grid(axis='y',alpha=0.75)
    if 'Corr' in var_name:
        plt.xlim(0.1, 0.4)
    elif 'Sal' in var_name:
        plt.xlim(0, 14)
    elif 'SIT' in var_name:
        plt.xlim(0, 4)
    elif 'Dens' in var_name:
        plt.xlim(100, 500)
    else:
        plt.xlim(0, 0.6)
    # if var_name.strip() ==  
    if site is not None:
        plt.title('Distribution of' + var_name + 'at ' + site,fontsize=14,fontweight='bold')
    else:
        plt.title('Distribution of' + var_name + 'at CB',fontsize=14,fontweight='bold')
    
    plt.xlabel(var_name + unit,fontsize=12)
    leg = ax.legend(prop={"size":12})
    plt.ylabel('Count',fontsize=14)
    plt.show()

# Make snowpit class
class Snowpit_info():
    # init method or constructor
    def __init__(self, site):
        self.SD = []
        self.rho = []
        self.corr = []
        self.temp = []
        self.sal = []
        self.Tsur = []
        self.Tsur_t2m = []
        self.Tsur_skt = []
        self.Tcal_snow = []
        self.Tcal_ice = []
        self.SIT = []
        self.site = site
        self.lat = []
        self.lon = []
        self.date = []
    def get_mean(self):
        self.SD_tot = np.array([np.sum(np.array(self.SD[i])) for i in range(len(self.SD))])
        self.rho_m  = np.array([np.mean(np.array(self.rho[i])) for i in range(len(self.rho))])
        self.corr_m = np.array([np.mean(np.array(self.corr[i])) for i in range(len(self.corr))])
        self.temp_m = np.array([np.mean(np.array(self.temp[i])) for i in range(len(self.temp))])
        self.sal_m  = np.array([np.mean(np.array(self.sal[i])) for i in range(len(self.sal))])
        ## remove outliers
        # self.SIT = self.SIT[self.SIT<2.0]
        self.SIT_m = np.mean(self.SIT) # - np.mean(SD_tot)
        
    def plot_var(self, var, var_name, unit):
        plot_vars(var, var_name, unit, site=self.site)
        
    def plot_all_var(self):
        plot_vars(self.SD_tot, ' SD ', 'm', site=self.site)
        plot_vars(self.rho_m, ' Density ', 'kg/m^3', site=self.site)
        plot_vars(self.corr_m*10**3, ' Correlation length ', 'mm', site=self.site)
        plot_vars(self.temp_m, ' Temperature ', 'K', site=self.site)
        plot_vars(self.sal_m*10**3, ' Salinity ', 'PSU', site=self.site)
    
    def plot_indv(self, var, varname, plot_all=False):
        plt.figure()
        plt.gca().invert_yaxis()
        for i in range(len(self.SD)):
            SD_cum = np.array(np.cumsum(self.SD[i]))
            varss = np.array(var[i])
            if plot_all == True:
                plt.plot(np.array(var[i]), SD_cum)
                # np.polyfit( )
        if plot_all == False:
            plt.plot(np.mean(varss), np.mean(SD_cum))
        plt.title('profiles of ' + varname + ' at ' + self.site)
        plt.ylabel('SD [m]')
        plt.xlabel(varname)
        plt.grid()
        plt.show()
        
    def plot_avg_profile(self, var, varname):
        plt.figure()
        plt.gca().invert_yaxis()
        for i in range(len(self.SD)):
            SD_cum = np.array(np.cumsum(self.SD[i]))
            varss = np.array(var[i])
        
        plt.title('profiles of ' + varname + ' at ' + self.site)
        plt.ylabel('SD [m]')
        plt.xlabel(varname)
        plt.grid()
        plt.show()
        
    def plot_all_indv(self, plot_indv):
  
        # plot_indv(self.SD_tot, ' SD [m]')
        plot_indv(self.rho, ' Density [kg/m^3]')
        plot_indv(self.corr*10**3, ' Correlation length [mm]')
        plot_indv(self.temp, ' Temperature [K]')
        plot_indv(self.sal*10**3, ' Salinity [PSU]')   