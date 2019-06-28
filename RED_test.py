#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 13:56:03 2018

@author: aa760
"""

from RED_MODULE import RED
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import os
from datetime import datetime

script_dir = os.path.dirname(__file__)
fig_dir = os.path.join(script_dir,'Output/figures/red_test/')

### Time values
maxtime= 500
timestep = 0.2
T = np.arange(0.0,maxtime,timestep)
ktime = len(T)
dt = np.zeros(ktime)


num_pft = 9



#%%


# Change to fit number of PFTs.
M = np.zeros((ktime,num_pft))
# Change to fit number of mass classes. (currently 10)
N = np.zeros((ktime,num_pft,10))
h_mean = np.zeros((ktime,num_pft))
P_init = [0.0, 0.0, 0.0, 0.0, 0.2565128707298727,0.2176408878894467, 0.1830805822303801,0.4404134877929182, 0.23952657266548183]
P = P_init * np.ones((ktime,num_pft))
nu_tot = np.zeros((ktime,num_pft))

nu_init =   [0.0, 0.0, 0.02296552248299122, 0.0514194592833519,0.599310040473938, 0.15302607417106628, 0.0,0.09776247292757034, 0.0505065992474556, 0.0,0.009994689375162125, 0.014850177802145481, 0.0]
nu_obs = [0.0, 0.0, 0.02296552248299122, 0.0514194592833519,0.599310040473938, 0.15302607417106628, 0.0,0.09776247292757034, 0.0505065992474556, 0.0,0.009994689375162125, 0.014850177802145481, 0.0]
gamma_init = [np.inf, np.inf, np.inf, np.inf, 0.010997714689771407,0.05136042738161767, 0.033836159647465947, 0.12265408740797126, 0.022718759602234066]
N2 = np.zeros((ktime,num_pft,10))
M2 = np.zeros((ktime,num_pft))
h_mean2 = np.zeros((ktime,num_pft))
nu_tot2 = np.zeros((ktime,num_pft))

        
for k in range(0,ktime - 1):
    

    if k == 0:

         m, M[k,:], N[k,:,:], h_mean[k,:], nu_tot[k,:], gamma_init2 = RED(0,nu_tot=nu_init[0:9],P = P_init, init=True,init_typ = 'nu_P')
 #       m, M[k,:], N[k,:,:], h_mean[k,:], nu_tot[k,:], gamma_init2 = RED(0,gamma_init=gamma_init,P = P_init, init=True,init_typ = 'P_gamma_init')

            
        

    year = np.floor(T[k])
    year_prog = T[k] - year


    m, M[k+1,:], N[k+1,:,:], h_mean[k+1,:], nu_tot[k+1,0:9], _ = RED(timestep,N=N[k,:,:],nu_tot=nu_tot[k,0:9],gamma_init=gamma_init2,P = P[k,:],init=False)
    
    _, M2[k+1,:], N2[k+1,:,:], h_mean2[k+1,:], nu_tot2[k+1,0:9], _ = RED(timestep,N=N2[k,:,:],nu_tot=nu_tot2[k,0:9],gamma_init=gamma_init2,P=P[k,:],init=False)

#%%

pft_col = list(np.zeros(10))
pft_name = list(np.zeros(10))
pft_handel = list(np.zeros(10))

pft_name[0] = 'BET-Tr'
pft_name[1] = 'BET-Te'
pft_name[2] = 'BDT'
pft_name[3] = 'NET'
pft_name[4] = 'NDT'
pft_name[5] = 'C3'
pft_name[6] = 'C4'
pft_name[7] = 'Esh'
pft_name[8] = 'Dsh'
pft_name[9] = 'Total'

pft_col[0] = [78.0/255.0,48.0/255.0,132.0/255.0]    # BET-Tr
pft_col[1] = [160.5/255.0,100.0/255.0,150.5/255.0]# BET-Te
pft_col[2] = [223.0/255.0,114.0/255.0,1.0,1.0]  # BDT
pft_col[3] = [150.0/255.0,200.0/255.0,120.0/255.0]  # NET
pft_col[4] = [0.0/255.0,158.0/255.0,96.0/255.0]  # NDT
pft_col[5] = [255.0/255.0,211.0/255.0,0.0/255.0]  # C3
pft_col[6] = [218.0/255.0,165.0/255.0,32.0/255.0]  # C4
pft_col[7] = [255.0/255.0,0.0/255.0,0.0/255.0] # Esh
pft_col[8] = [102.0/255.0,0.0/255.0,0.0/255.0]  # Dsh
pft_col[9] = [0.0/255.0, 0.0/255.0, 0.0/255.0]  # Total

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

#%% Plots
pfts = [4,5,6,7]
with PdfPages(fig_dir+'PFT_outputs_nuP.pdf') as pdf:
    # plot 1
    
    rows = 2
    cols = 1
    
    fig = plt.figure(figsize=(11.69,8.27))
    gs = gridspec.GridSpec(rows,cols)
    
    ax = fig.add_subplot(gs[0,0]) 

    for pft in pfts:
        ax.plot(T,nu_tot[:,pft],color=pft_col[pft])
        ax.plot(T,nu_tot2[:,pft],color=pft_col[pft],label=pft_name[pft],linestyle=':')
        ax.plot([T[0],T[-1]],[nu_obs[pft],nu_obs[pft]],color=pft_col[pft],linestyle='--',marker='|')

        pft_handel[pft] = Patch(color=pft_col[pft], label=pft_name[pft])
    

    ax.legend(pft_handel,pft_name,bbox_to_anchor=(1.005, 1))
    ax.set_xlabel('Stand Age')
    ax.set_ylabel(r'$\nu$')


    ax = fig.add_subplot(gs[1,0]) 
    
    for pft in range(0,num_pft):
        ax.plot(T,M[:,pft],color=pft_col[pft])
        ax.plot(T,M2[:,pft],color=pft_col[pft],linestyle=':')
        
    ax.set_xlabel('Stand Age')
    ax.set_ylabel(r'$M \ (\mathrm{kg \ m^{-2}})$')
    
    fig.suptitle('Gamma Var')
    fig.tight_layout()    
    pdf.attach_note(r'Eq_PFT_outputs')    
    pdf.savefig()
    plt.show()
    plt.close(fig)
    
"""
with PdfPages(fig_dir+'PFT_outputs_Pgam.pdf') as pdf:
    # plot 1
    
    rows = 2
    cols = 1
    
    fig = plt.figure(figsize=(11.69,8.27))
    gs = gridspec.GridSpec(rows,cols)
    
    ax = fig.add_subplot(gs[0,0]) 

    for pft in pfts:
        ax.plot(T,nu_tot[:,pft],color=pft_col[pft])
        ax.plot(T,nu_tot2[:,pft],color=pft_col[pft],label=pft_name[pft],linestyle=':')

        pft_handel[pft] = Patch(color=pft_col[pft], label=pft_name[pft])
    

    ax.legend(pft_handel,pft_name,bbox_to_anchor=(1.005, 1))
    ax.set_xlabel('Stand Age')
    ax.set_ylabel(r'$\nu$')


    ax = fig.add_subplot(gs[1,0]) 
    
    for pft in pfts:
        ax.plot(T,M[:,pft],color=pft_col[pft])
        ax.plot(T,M2[:,pft],color=pft_col[pft],linestyle=':')
    
    ax.set_yscale('log')
    ax.set_xlabel('Stand Age')
    ax.set_ylabel(r'$M \ (\mathrm{kg \ m^{-2}})$')
    
    fig.suptitle('Gamma Fixed')
    fig.tight_layout()    
    pdf.attach_note(r'Eq_PFT_outputs')    
    pdf.savefig()
    plt.show()
    plt.close(fig)
    
"""