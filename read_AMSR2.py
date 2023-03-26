# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 17:00:00 2022

@author: -
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np

from polar_plots_pred import plot

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import cartopy.crs as ccrs
import numpy as np

def Box(ax, length=None, location=(0.5, 0.05), linewidth=2.0, c='k'):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """

    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
    print('coords', [llx0, llx1, lly0, lly1])
    
    # if loc_type='lat/lon':
    # location with respect to figure
    lx = (location[0]-llx0)/(llx1-llx0)
    ly = (location[1]-lly0)/(lly1-lly0)
    print(location)
    # print(ly)
    
    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    # sbllx = (llx1 + llx0) / 2 
    sbllx = llx0 + (llx1 - llx0) * lx
    # sblly = lly0 + (lly1 - lly0) * location[1]
    sblly = lly0 + (lly1 - lly0) * ly
    tmc = ccrs.TransverseMercator(sbllx, sblly)
    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    #sbx = x0 + (x1 - x0) * location[0]
    #sby = y0 + (y1 - y0) * location[1]
    sbx = x0 + (x1 - x0) * lx
    sby = y0 + (y1 - y0) * ly

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length: 
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length) 

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    bar_ys = [sby - length * 500, sby + length * 500]

    #Plot the scalebar
    #ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
    #ax.plot([sbx, sbx], bar_ys, transform=tmc, color='k', linewidth=linewidth)
    ax.plot(bar_xs, [sby- length * 500, sby- length * 500], transform=tmc, color=c, linewidth=linewidth, alpha=0.8)
    ax.plot([sbx- length * 500, sbx- length * 500], bar_ys, transform=tmc, color=c, linewidth=linewidth, alpha=0.8)
    ax.plot([sbx+ length * 500, sbx+ length * 500], bar_ys, transform=tmc, color=c, linewidth=linewidth, alpha=0.8)
    ax.plot(bar_xs, [sby+ length * 500, sby+ length * 500], transform=tmc, color=c, linewidth=linewidth, alpha=0.8)
    #Plot the scalebar label
    # ax.text(sbx, sby, str(length) + ' km', transform=tmc,
    #         horizontalalignment='center', verticalalignment='bottom')

def scale_bar(ax, length=None, location=(0.5, 0.05), linewidth=2.0, c='k'):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """

    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
    print('coords', [llx0, llx1, lly0, lly1])

    
    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2 
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly)
    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]
    
    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length: 
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length) 

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    bar_ys = [sby - length * 500, sby + length * 500]

    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc, 
            horizontalalignment='center', verticalalignment='bottom')
class AMSR2_data():
    
    # init method or constructor
    def __init__(self, filename):
        h5 = h5py.File(filename,'r')
        self.lat   = np.array(h5['HDFEOS/GRIDS/NpPolarGrid25km'].get('lat'))
        self.lon   = np.array(h5['HDFEOS/GRIDS/NpPolarGrid25km'].get('lon'))
        self.freqs = ['06', '10', '18', '23', '36']
        self.projx = []
        self.projy = []
        self.KD_index = []
        self.namesH = ['6GHzH', '10GHzH', '18GHzH', '23GHzH', '36GHzH']
        self.namesV = ['6GHzV', '10GHzV', '18GHzV', '23GHzV', '36GHzV']
        
        self.DataSetH = np.zeros((self.lat.shape[0], self.lat.shape[1], len(self.freqs) ))
        self.DataSetV = np.zeros((self.lat.shape[0], self.lat.shape[1], len(self.freqs) ))

        self.DataSetH_asc = np.zeros((self.lat.shape[0], self.lat.shape[1], len(self.freqs) ))
        self.DataSetV_asc = np.zeros((self.lat.shape[0], self.lat.shape[1], len(self.freqs) ))        

        self.DataSetH_dec = np.zeros((self.lat.shape[0], self.lat.shape[1], len(self.freqs) ))
        self.DataSetV_dec = np.zeros((self.lat.shape[0], self.lat.shape[1], len(self.freqs) ))
        
        for i in range(len(self.freqs)):
            ## Get mean daily values
            self.DataSetH[:, :, i] = h5['HDFEOS/GRIDS/NpPolarGrid25km/Data Fields'].get('SI_25km_NH_' + self.freqs[i] + 'H_DAY')[...]/10
            self.DataSetV[:, :, i] = h5['HDFEOS/GRIDS/NpPolarGrid25km/Data Fields'].get('SI_25km_NH_' + self.freqs[i] + 'V_DAY')[...]/10
            ## Get ascending
            self.DataSetH_asc[:, :, i] = h5['HDFEOS/GRIDS/NpPolarGrid25km/Data Fields'].get('SI_25km_NH_' + self.freqs[i] + 'H_ASC')[...]/10
            self.DataSetV_asc[:, :, i] = h5['HDFEOS/GRIDS/NpPolarGrid25km/Data Fields'].get('SI_25km_NH_' + self.freqs[i] + 'V_ASC')[...]/10
            ## Get decending
            self.DataSetH_dec[:, :, i] = h5['HDFEOS/GRIDS/NpPolarGrid25km/Data Fields'].get('SI_25km_NH_' + self.freqs[i] + 'H_DSC')[...]/10
            self.DataSetV_dec[:, :, i] = h5['HDFEOS/GRIDS/NpPolarGrid25km/Data Fields'].get('SI_25km_NH_' + self.freqs[i] + 'V_DSC')[...]/10
      
    def plot(self, index, count, polarization, plot_point = [], obs=False, c='red', CB=True):
                
        x = self.lat[index]
        y = self.lon[index]
        if polarization == 'V':
            z = self.DataSetV[:, :, count].flatten()[index]
            name = self.namesV[count]
        else:
            z = self.DataSetH[:, :, count].flatten()[index]
            name = self.namesH[count]

        extent = [-107.5, -100.5, 67, 70]
        if obs==True:
            extent = [-106.5, -104, 68.6, 69.2]
            
        plt.figure(figsize=(7,7))
        plt.axis('off')
        # plt.title('AMSR2 data at:' + name)
        # ax = plt.axes(projection=ccrs.NorthPolarStereo())
        ax = plt.axes(projection=ccrs.Mercator(central_longitude=-50))
        # ax.set_facecolor('white') #white
        ax.coastlines()
        # ax.add_feature(cfeature.LAND, zorder=100, facecolor='black')
        ax.add_feature(cfeature.LAND, alpha=0.5)
        ax.add_feature(cfeature.OCEAN, alpha=0.5) 
        gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlabel_style = {'size': 14, 'color': 'black'}
        gl.ylabel_style = {'size': 14, 'color': 'black'}
        ax.set_extent(extent, ccrs.PlateCarree())

        plott = ax.scatter(y, x, s=40, c='red', alpha=1,transform=ccrs.PlateCarree(),)
        
        if plot_point: # if list is not empty
            try:
                if len(plot_point)>1 and obs==False:
                    xx = [plot_point[i].point[0] for i in range(len(plot_point))]
                    yy = [plot_point[i].point[1] for i in range(len(plot_point))]
                elif len(plot_point)>1 and obs==True:
                    xx = [plot_point[i][0] for i in range(len(plot_point))]
                    yy = [plot_point[i][1] for i in range(len(plot_point))]
            except: # point has no len
                xx = plot_point.point[0]
                yy = plot_point.point[1]
            plt.scatter(yy, xx, s=40, edgecolor='k', linewidth=0.5, facecolor=c, alpha=1,transform=ccrs.PlateCarree(),)
            scale_bar(ax, 25)
            for lo, la,cc in zip(yy,xx,c):
                Box(ax, 25, location=(lo, la), c=cc)
                # Box(ax, 25, location=(lo, la), c=cc)
           
        if CB==True:
            latCB = [69.08897, 68.99985, 68.91067, 68.82158]		
            lonCB = [-105.50422, -105.47468, -105.41325, -105.50422]
            plt.scatter(lonCB, latCB, s=30, c='k', alpha=1,transform=ccrs.PlateCarree(),)
        # plt.colorbar(plott)
        from matplotlib.lines import Line2D
        plt.legend(handles=[Line2D([0], [0], marker='o', color='w',markerfacecolor='k', markersize=10, label='CB pit locations')],prop={"size":12}, framealpha=1.0)


    def plot_asc(self, index, count, polarization):       
        x = self.lat.flatten()[index]
        y = self.lon.flatten()[index]
        if polarization == 'V':
            z = self.DataSetV_asc[:, :, count].flatten()[index]
            name = self.namesV[count]
        else:
            z = self.DataSetH_asc[:, :, count].flatten()[index]
            name = self.namesH[count]

        extent = [-110, -100, 67, 70]
        
        plt.figure()
        plt.axis('off')
        plt.title('AMSR2 data at:' + name)
        ax = plt.axes(projection=ccrs.Mercator(central_longitude=-50))
        ax.set_facecolor('white') #white
        ax.coastlines()
        ax.add_feature(cfeature.LAND, zorder=100, facecolor='black')
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlabel_style = {'size': 14, 'color': 'black'}
        gl.ylabel_style = {'size': 14, 'color': 'black'}
        ax.set_extent(extent, ccrs.PlateCarree())
        
        plot = plt.scatter(y, x, s=40, c=z, alpha=1,transform=ccrs.Geodetic(),)        
        plt.colorbar(plot)


    def plot_dec(self, index, count, polarization):
        
        x = self.lat.flatten()[index]
        y = self.lon.flatten()[index]
        if polarization == 'V':
            z = self.DataSetV_dec[:, :, count].flatten()[index]
            name = self.namesV[count]
        else:
            z = self.DataSetH_dec[:, :, count].flatten()[index]
            name = self.namesH[count]

        extent = [-110, -100, 67, 70]
        
        plt.figure()
        plt.axis('off')
        plt.title('AMSR2 data at:' + name)
        ax = plt.axes(projection=ccrs.Mercator(central_longitude=-50))
        ax.set_facecolor('white') #white
        ax.coastlines()
        ax.add_feature(cfeature.LAND, zorder=100, facecolor='black')
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlabel_style = {'size': 14, 'color': 'black'}
        gl.ylabel_style = {'size': 14, 'color': 'black'}
        ax.set_extent(extent, ccrs.PlateCarree())
        
        plot = plt.scatter(y, x, s=40, c=z, alpha=1,transform=ccrs.Geodetic(),)        
        plt.colorbar(plot)
