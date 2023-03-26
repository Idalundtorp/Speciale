#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 12:03:11 2021

@author: s174020

"""
import os
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
def plot_one(lat, lon, z, title, ylabel='SIT [m]', clim=False, savename=False, saveloc=False, s=8):    
    plt.figure(figsize=(7, 7))

    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,65,90],ccrs.PlateCarree())
    # ax.set_extent([np.min(lon)-0.05,np.max(lon)+0.05,np.min(lat)-0.01,np.max(lat)+0.01],ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.OCEAN)        
    ax.add_feature(cfeature.LAND)
    ax.gridlines(draw_labels=False)
    
    plot=plt.scatter(lon, lat, c=z, s=s, transform=ccrs.PlateCarree(),)
    plt.title(title,fontsize=16,fontweight='bold')
    
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(plot, cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel(ylabel);
    cb.ax.tick_params(labelsize=12) 
       
    if clim==False:
        plt.clim(np.nanmin(z),np.nanmax(z))
    else:
        plt.clim(clim)
    if savename:
        if saveloc:
            plt.savefig(os.path.join(saveloc,savename), bbox_inches='tight')
            plt.show()
        else:
            plt.savefig(savename, bbox_inches='tight')


def plot(self, lat,lon,z, title, ylabel='SIT [m]', clim=False, savename=False, saveloc=False, s=8):
    
    plt.figure(figsize=(7, 7))

    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    #ax.set_extent([-180,180,65,90],ccrs.PlateCarree())
    ax.set_extent([np.min(lon)-0.05,np.max(lon)+0.05,np.min(lat)-0.01,np.max(lat)+0.01],ccrs.PlateCarree())
    ax.coastlines()
    # ax.add_feature(cfeature.OCEAN)        
    ax.add_feature(cfeature.LAND)
    ax.gridlines(draw_labels=False)
    
    # c = list(range(len(z)))
    # clim = [0, c[-1]]
    # ylabel = 'count'
    plot=plt.scatter(lon, lat, c=z, alpha=0.5, transform=ccrs.PlateCarree(),)
    
    # plt.scatter(lon[0], lat[0], s=100, c='g',transform=ccrs.PlateCarree(),)
    plt.scatter(self.lon, self.lat, s=50, c='red', transform=ccrs.PlateCarree(),)
    # plot=plt.plot(lat, lon, c=z, transform=ccrs.PlateCarree(),)
    plt.title(title,fontsize=16,fontweight='bold')
    
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(plot, cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel(ylabel);
    cb.ax.tick_params(labelsize=12) 
       
    if clim==False:
        plt.clim(np.nanmin(z),np.nanmax(z))
    else:
        plt.clim(clim)
    if savename:
        if saveloc:
            plt.savefig(os.path.join(saveloc,savename), bbox_inches='tight')
            plt.show()
        else:
            plt.savefig(savename, bbox_inches='tight')

def polar_plot(self, z, title, ylabel='SIT [m]', clim=False, savename=False, saveloc=False, s=8):
    
    plt.figure(figsize=(7, 7))

    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    #ax.set_extent([-180,180,65,90],ccrs.PlateCarree())
    ax.set_extent([np.min(self.lon)-0.05,np.max(self.lon)+0.05,np.min(self.lat)-0.01,np.max(self.lat)+0.01],ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.OCEAN)        
    ax.add_feature(cfeature.LAND)
    ax.gridlines(draw_labels=False)
    plot =plt.scatter(self.lon, self.lat, s=s, c=z, transform=ccrs.PlateCarree(),)
    plt.title(title,fontsize=16,fontweight='bold')
    try:
    	ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    	cb = plt.colorbar(plot, cax=ax_cb, orientation='vertical')
    except:
        cb = plt.colorbar()
    cb.ax.set_ylabel(ylabel);
    cb.ax.tick_params(labelsize=12) 
       
    if clim==False:
        plt.clim(np.nanmin(z),np.nanmax(z))
    else:
        plt.clim(clim)
    if savename:
        if saveloc:
            plt.savefig(os.path.join(saveloc,savename), bbox_inches='tight')
            plt.show()
        else:
            plt.savefig(savename, bbox_inches='tight')