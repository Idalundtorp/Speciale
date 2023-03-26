# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 17:20:04 2022

@author: Ida Olsen

Load simulation data and do Bias correction
"""
import os
import numpy as np
import matplotlib.pyplot as plt

datapath = "C:/Users/Ida Olsen/Documents/Speciale_2022/Optimal_Estimation/OE_ida/Simulated_data/TB_simulated.npy"
# plot path
pp = "C:/Users/Ida Olsen/Documents/Speciale_2022/figures"

with open(datapath, 'rb') as f:
    TbV = np.load(f)
    TbH = np.load(f)
    b = np.load(f)
    b = np.load(f)
    b = np.load(f)
    b = np.load(f)
    b = np.load(f)
    b = np.load(f)
    b = np.load(f)
    b = np.load(f)
    TbV_AMSR_all = np.load(f)
    TbH_AMSR_all = np.load(f)



frequencies=[6.925e9,10.65e9,18.7e9,23.8e9,36.5e9];  # [Hz]
freqs = [f*1e-9 for f in frequencies]
plt.figure(figsize=(8,6))
for TBV, TBH in zip(TbV, TbH):
    plt.scatter(freqs, TBV, c='r')
    plt.scatter(freqs, TBH, c='b')
plt.legend(['TBV', 'TBH'])
plt.xlabel('Freq GHz')
plt.grid()
plt.ylabel('TB [K]')
plt.title('Simulated FYI TBs from OIB data 2013-2019')
plt.savefig(os.path.join(pp, 'OE_figures/Simulated_data.png'), bbox_inches='tight')
plt.show()



## SIT
avgTbV=np.mean(TbV)
TbVstd=np.std(TbV)

avgTbV_AMSR=np.mean(TbV_AMSR_all)
TbV_AMSRstd=np.std(TbV_AMSR_all)

bins = np.linspace(200, 300, 200)
fig, ax= plt.subplots()
for i, f in zip(range(5), freqs):
    bins = np.linspace(200, 300, 200)
    fig, ax= plt.subplots()
    ax.hist(x=np.concatenate((TbV_AMSR_all))[i::5], bins=bins, color='red', alpha=0.5, rwidth=0.85, 
            label='satelitte AMSR2 TbV  with avg:'+str(np.round(np.mean(np.concatenate((TbV_AMSR_all))[i::5]),2)) + ' and std of: '+str(np.round(np.std(np.concatenate((TbV_AMSR_all))[i::5]),2)))
    ax.hist(np.concatenate((TbV))[i::5], bins, alpha=0.5, rwidth=0.85,
            label='Sim. TbV  with avg:'+str(np.round(np.mean(np.concatenate((TbV))[i::5]),2))+' and std of: '+str(np.round(np.std(np.concatenate((TbV))[i::5]),2)))
    plt.grid(axis='y',alpha=0.75)
    # plt.xlim(xmax=5,xmin=-0.1)
    # plt.ylim(ymax=42,ymin=0)
    plt.title('TBV freq:' + str(np.round(f,2)),fontsize=14,fontweight='bold')
    plt.xlabel('TbV  [K]',fontsize=12)
    plt.legend(prop={"size":12}, bbox_to_anchor=(1.0, -0.15))
    plt.ylabel('Count',fontsize=14)
    # plt.savefig(name + '_hist_TbV .png', bbox_inches='tight')
    plt.show()