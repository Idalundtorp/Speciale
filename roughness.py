# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 09:08:01 2021

@author: Olsen

Roughness computation is based on paper:
A Parameterized Multifrequency-Polarization Surface Emission Model
by shi et al
"""
# Python packages
import numpy as np
import matplotlib.pyplot as plt
# other packages

# python packages
import numpy as np
 
def Get_slope(site, variance=False, plot=False):
    "Calculation of slope from roughness value and correlation length"
    
    file = "C:/Users/Ida Olsen/Documents/Speciale_2022/Ida_passive_SMRT/Ida_passive_SMRT/Code/smrt_in-out/singleParam_CB2022.csv"
    data = np.genfromtxt(file, delimiter=',', names=True, usecols=(1,2,3,4,5,6,7,8,9,10,11,12), encoding='utf8')
    
    sites = np.array(['AK1', 'AK2', 'AK3', 'AK4'])
    
    # s/l
    slope = data['ice_rms']/data['ice_lc']
    print(slope)
    
    index = np.where(site==sites)[0]
    if variance==True:
        var = np.std(slope)**2
        return [slope[index],var]
    else:
        return slope[index]

def roughness(tcwv,res,avg_Tsli,m=False,site='AK1',several=False):
    '''
    Parameters
    ----------
    tcwv : np.array of floats
        total collumnar water vapour
    m : float
        roughness property = s/l, similar to surface slope (sigma = sqrt(2)*s/l)
        s: rms height
        l: correlation length
    res : SMRT res object
        Result from SMRT simulation contains simulated brightness temperatures
    avg_Tsli : float
        Average snow/ice interface tempature
    several : True/False, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    list
        Modified Brightness temperatures, 
        Polarization mixing constants, 
        roughness parameters (see fig.5 in Shi et al)

    '''
    ## water vapor from ECMWF
    x=tcwv
    
    ## linear regression for Tup=Tdn
    Tdn69=np.mean(0.0297*x+4.0639)
    Tdn107=np.mean(0.0865*x+4.3262)
    Tdn187=np.mean(0.7729*x+5.6202)
    Tdn238=np.mean(2.2411*x+7.2606)
    Tdn365=np.mean(0.8701*x+18.208)
    Tdn=np.array([Tdn69,Tdn107,Tdn187,Tdn238,Tdn365])

    ## linear regressions for transmissivity
    tau69=np.mean(-1e-5*x+0.983)
    tau107=np.mean(-0.0002*x+0.9814)
    tau187=np.mean(-0.0028*x+0.9755)
    tau238=np.mean(-0.0082*x+0.9683)
    tau365=np.mean(-0.0028*x+0.9211)
    tau=np.array([tau69,tau107,tau187,tau238,tau365])
    Tsp=2.7 #K - background radiation
    
    
    # Compute roughness
    if m is False:    
        m=Get_slope(site)
        # if site=='AK1':
        #     m+=0.15


    ## See table IV Shi et al.
    Q_H_log=5.603623+3.095002*np.log(m)-9.377554*m
    Q_Htest=np.exp(Q_H_log)

    Q_V_log=3.216481+2.452773*np.log(m)-6.674134*m
    Q_Vtest=np.exp(Q_V_log)
    
    ## From Ulaby
    Beta_H69=0.958
    Beta_H187=1.013
    Beta_H238=1.015
    Beta_H365=1.026
    
    Q_H69=Beta_H69*Q_Htest
    Q_H187=Beta_H187*Q_Htest
    Q_H238=Beta_H238*Q_Htest
    Q_H365=Beta_H365*Q_Htest
    QH=np.array([Q_H69,Q_Htest,Q_H187,Q_H238,Q_H365]).flatten()
    QV=np.array([Q_Vtest,Q_Vtest,Q_Vtest,Q_Vtest,Q_Vtest]).flatten()
    
    #%% fig, ax = plt.subplots()
    # frequency=[6.9,10.7,18.7,23.8,36.5]
    # plot3,=plt.plot(frequency, QV,'.-m',label='plot3')
    # plot4,=plt.plot(frequency, QH,'.-g',label='plot4')
    # plt.title('Qplot')
    # plt.xlabel('Frequencies')
    # plt.ylabel('Q value')
    # plt.legend(handles=[plot3,plot4])
    # plt.legend([plot3,plot4], ['QV','QH'])
    # plt.show()
    
    #%% calculating Gamma
    if several==True:
        # print('recieving correct')
        TBV = np.array(np.mean(res.TbV(),axis=1))
        TBH = np.array(np.mean(res.TbH(),axis=1))
    elif several==False:
        TBV = res.TbV()
        TBH = res.TbH()
    Gamma_V=(1-np.array(TBV/avg_Tsli))
    Gamma_H=(1-np.array(TBH/avg_Tsli))
    
    
        
    # mixing constants
    Gamma_rs_H=(1-QH)*Gamma_H+QH*Gamma_V
    Gamma_rs_V=(1-QV)*Gamma_V+QV*Gamma_H
    
    TB_AMSR_V=Tdn+tau*(1-Gamma_rs_V)*avg_Tsli+tau*Tdn*Gamma_rs_V+np.square(tau)*Tsp*Gamma_rs_V
    
    TB_AMSR_H=Tdn+tau*(1-Gamma_rs_H)*avg_Tsli+tau*Tdn*Gamma_rs_H+np.square(tau)*Tsp*Gamma_rs_H

    
    return [TB_AMSR_V,TB_AMSR_H, Gamma_rs_H, Gamma_rs_V, QH, QV]
