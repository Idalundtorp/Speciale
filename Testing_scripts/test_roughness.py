# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 14:52:33 2022

@author: Ida Olsen
"""
import numpy as np
import sys
import os
sys.path.append('C:/Users/Ida Olsen/Documents/Speciale_2022/scripts')
from roughness import roughness
from Make_plots import plot_freq
import matplotlib.pyplot as plt

def test_rough(s, A, E, res):
    '''

    Parameters
    ----------
    s : class holding information of each site
        Contains several attributes.
        This function uses:
            s.site (site name)
            s.temp (average temperature per profile in site)
    res : SMRT simulation output contains simulated TBs
          Data is stored in the res class
    Returns
    -------
    Plots output for roughness values in the range 0.01 to 0.5

    '''
    # tcwv = np.load(s.site + '_ECMWF_tcwv.npy', allow_pickle=True)
    rough = np.linspace(0.008, 0.3, 5)
    data = np.zeros((len(rough)*2,5))
    start = 0
    slut = 2
    
    plt.figure(figsize=(7,6))
    frequencies=np.array([6.925,10.65,18.7,23.8,36.5])
    # colors = [['b', 'orange'] for m in rough]
    colors = ['C' + str(int(num)) for num in range(len(rough))]
    for m, c in zip(rough, colors):
        start += 2 
        slut += 2
        [TB_SMRT_V,TB_SMRT_H, Gamma_rs_H, Gamma_rs_V, QH, QV] = roughness(E.tcwv,res,np.mean(s.Tcal_ice),m=m,several=True,site=s.site)
        
        data[start:slut-1,:] = TB_SMRT_V
        data[start+1:slut,:] = TB_SMRT_H
        
        # data[-2:-1,:] = satdataV_landfree.flatten()
        # data[-1:,:] = satdataH_landfree.flatten()
        # colors.append(['g', 'k'])
        
    # labels = ['TB_SMRT_V','TB_SMRT_H']
    
        plot1,=plt.plot(frequencies, TB_SMRT_V,'.-', linewidth = 1 + m*3, color=c, label = 'ROUGHNESS: ' + str(np.round(m,2)))
        plot2,=plt.plot(frequencies, TB_SMRT_H,'.-', linewidth = 1 + m*3, color=c)
    plot3,=plt.plot(frequencies, A.DataSetV.flatten(),'.-', linewidth = 2, color='k', label='AMSR2 TB')
    plot4,=plt.plot(frequencies, A.DataSetH.flatten(),'.-', linewidth = 2, color='k')
    plt.ylim([200,260])
    plt.grid()
    plt.legend(loc='lower right',fontsize=12) # ,bbox_to_anchor=(1, 0.5))
    plt.title('Roughness sensitivity study at site ' + s.site,fontsize=16)
    plt.xlabel('Frequency [GHz]',fontsize=14)
    plt.ylabel('TB [K]',fontsize=14)
    plt.savefig(os.getcwd() + '/Testing_scripts/figures/' + s.site + '_roughness_test.PNG', bbox_inches='tight')
    plt.show()
    # plot_freq(s, data, labels=labels, cs=np.array(colors).flatten(), title='Roughness test at ' + s.site)
