# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 14:19:47 2022

@author: Ida Olsen
"""
import matplotlib.pyplot as plt
import numpy as np
import os

basepath =  os.path.dirname(os.path.abspath(os.getcwd()))
# plot path
pp = "C:/Users/Ida Olsen/Documents/Speciale_2022/figures"
def plot_freq(self, y, labels='', cs='r', err=False, title='Tbs simulated at CB', savefig=True):
    plt.figure(figsize=(7,6))
    frequencies = np.array([6.925,10.65,18.7,23.8,36.5]);
    if len(y.shape) > 1:
            for i in range(y.shape[0]):
                try:
                    c = cs[i]
                except:
                    c = cs
                try:
                    label = labels[i]
                except:
                    if len(labels)>1:
                        label = ''
                    else:
                        label = labels
                if err is not False: ## Make plot with errorbars
                    plt.fill_between(frequencies, y[i]-err[i]/2, y[i]+err[i]/2, edgecolor=c, facecolor=c, alpha=0.3)
                plot4,=plt.plot(frequencies, y[i],'.-', color=c, label = label)
    else:
        plot4,=plt.plot(frequencies, y,'.-', color=cs, label = labels)
        if err is not False: ## Make plot with errorbars
            plt.fill_between(frequencies, y-err/2, y+err/2, edgecolor=c, facecolor=c, alpha=0.3)
    plt.ylim([200,260])
    plt.xlim([6, 37])
    plt.grid()
    plt.legend(loc='lower left',fontsize=12)
    plt.title(title,fontsize=16)
    plt.xlabel('Frequency [GHz]',fontsize=14)
    plt.ylabel('TB [K]',fontsize=14)
    if savefig==True:
        plt.savefig(os.path.join(basepath,'figures/CB_figures/Results_obtained_at_' + self.site + '.png'), bbox_inches='tight')
    plt.show()


def plot_simulation_parameters_from_CB(self, sd, count):
        cumsum_ref = np.concatenate([np.cumsum(sd_cb) for sd_cb in self.SD_CB])
        plt.figure(figsize=(7,6))
        plt.gca().invert_yaxis()
        plt.scatter(np.concatenate(self.lex_CB)*1e3, cumsum_ref, s=10, c='blue', label='lex all')
        x = np.linspace(np.min(np.concatenate(self.lex_CB)*1e3),np.max(np.concatenate(self.lex_CB)*1e3), len(np.concatenate(self.lex_CB)))
        y = [sd for lex in np.concatenate(self.lex_CB)]
        yb = [sd-0.05 for lex in np.concatenate(self.lex_CB)]
        yt = [sd+0.05 for lex in np.concatenate(self.lex_CB)]
        plt.title('Correlation lengths for profile with SD, ' + str(np.round(sd*100)) + ' [cm]', fontsize=14)
        plt.plot(x, y, label='SD OIB')
        plt.plot(x, yb, label='lower lim CB SD')
        plt.plot(x, yt, label='upper lim CB SD')
        plt.scatter(self.lex_top*1e3, self.cumsum[self.indexx_t], s=10, c='red', label='lex top layer')
        plt.scatter(self.lex_bottom*1e3, self.cumsum[self.indexx_b], s=10, c='green', label='lex bottom layer')
        plt.xlim([np.min(np.concatenate(self.lex_CB)*1e3),np.max(np.concatenate(self.lex_CB)*1e3)])
        plt.legend(loc='lower right', fontsize=12) #,bbox_to_anchor=(1, 0.5))
        plt.xlabel('Snow correlation length [mm]', fontsize=12)
        plt.ylabel('SD [m]', fontsize=12)
        plt.savefig(os.path.join(pp, 'OE_figures/corr_overview_' + str(count) +'.png'), bbox_inches='tight')
        plt.show()
        
        plt.figure(figsize=(7,6))
        plt.gca().invert_yaxis()
        plt.scatter(np.concatenate(self.rho_CB), cumsum_ref, s=10, c='blue', label='rho all')
        x = np.linspace(np.min(np.concatenate(self.rho_CB)),np.max(np.concatenate(self.rho_CB)), len(np.concatenate(self.rho_CB)))
        y = [sd for rho in np.concatenate(self.rho_CB)]
        yb = [sd-0.05 for rho in np.concatenate(self.rho_CB)]
        yt = [sd+0.05 for rho in np.concatenate(self.rho_CB)]
        plt.title('Densities for profile with SD, ' + str(np.round(sd*100)) + ' [cm]', fontsize=14)
        plt.plot(x, y, label='SD OIB')
        plt.plot(x, yb, label='lower lim CB SD')
        plt.plot(x, yt, label='upper lim CB SD')
        plt.scatter(self.rho_top, self.cumsum[self.indexx_t], s=10, c='red', label='rho top layer')
        plt.scatter(self.rho_bottom, self.cumsum[self.indexx_b], s=10, c='green', label='rho bottom layer')
        plt.xlim([np.min(np.concatenate(self.rho_CB)),np.max(np.concatenate(self.rho_CB))])
        plt.legend(loc='lower right', fontsize=12) #,bbox_to_anchor=(1, 0.5))
        plt.xlabel('Snow density [kg/m3]', fontsize=12)
        plt.ylabel('SD [m]', fontsize=12)
        plt.savefig(os.path.join(pp, 'OE_figures/density_overview_' + str(count) +'.png'), bbox_inches='tight')
        plt.show()
        
        plt.figure(figsize=(7,6))
        plt.gca().invert_yaxis()
        plt.scatter(np.concatenate(self.sal_CB)*1000, cumsum_ref, s=10, c='blue', label='sal all')
        x = np.linspace(np.min(np.concatenate(self.sal_CB))*1000,np.max(np.concatenate(self.sal_CB))*1000, len(np.concatenate(self.sal_CB)))
        y = [sd for sal in np.concatenate(self.sal_CB)]
        yb = [sd-0.05 for sal in np.concatenate(self.sal_CB)]
        yt = [sd+0.05 for sal in np.concatenate(self.sal_CB)]
        plt.title('salinities for profile with SD, ' + str(np.round(sd*100)) + ' [cm]', fontsize=14)
        plt.plot(x, y, label='SD OIB')
        plt.plot(x, yb, label='lower lim CB SD')
        plt.plot(x, yt, label='upper lim CB SD')
        plt.scatter(self.sal_top*1000, self.cumsum_sal_t, s=10, c='red', label='sal top layer')
        plt.scatter(self.sal_bottom*1000, self.cumsum_sal_b, s=10, c='green', label='sal bottom layer')
        plt.xlim([np.min(np.concatenate(self.sal_CB))*1000,np.max(np.concatenate(self.sal_CB))*1000])
        plt.legend(loc='lower right', fontsize=12) #,bbox_to_anchor=(1, 0.5))
        plt.xlabel('salinity [PSU]', fontsize=12)
        plt.ylabel('SD [m]', fontsize=12)
        plt.savefig(os.path.join(pp, 'OE_figures/salinity_overview_' + str(count) +'.png'), bbox_inches='tight')
        plt.show()
