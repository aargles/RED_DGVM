#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 13:56:03 2018

In this script we demonstrate how red functions and is initialised with a 
simple run setup. Firstly we declare a observation (i.e. coverage, biomass…) 
and a value, with a paired grid-box assimilate can be used to find a model 
equilibrium. Next we spin up from a coverage to show convergence towards
our diagnosed equilibrium.

@author: aa760
"""

import numpy as np
import matplotlib.pyplot as plt

# RED MODULES
from RED_MODULE import RED, RED_init
from RED_MODULE.RED_SUBROUTINE.MISC_METHODS import LOAD_PFT_VALUES


#%% Initiate RED
help(RED_init) 
# PFTs in element order:
# BET-Tr, Broadleaf Evergreen Tropical Tree
# BET-Te, Broadleaf Evergreen Temperate Tree
# BDT, Broadleaf Deciduous Tree
# NET, Needleleaf Evergreen Tree
# NDT, Needleleaf Deciduous Tree
# C3, Cool Season Grass
# C4, Tropical Grass
# ESh, Evergreen Shrub
# DSh, Deciduous Shrub


observation = 'coverage'
observation_type = 'total'
observation_value = [0.35, 0.3,0.,0.,0.,0.2,0.,0.,0.2]
P_or_gamma = 'P' # Switch 
P_obs = [1.4, 2.0,0.1,0.1,0.,0.3,0.,0.,0.8]  # Observed gridbox assimilate


PFT_competition = True # Adjusts the calculated observed vegetation fraction
                       # to take into account the PFT hierarchy.
overwrite_check = True # Overwrite check, rewrites the PFT values file.


output_key = ('N_m','nu_tot','M_tot','gamma_base','mu_0','g_m')

N_eq_m, nu_eq, M_tot, gamma_base,mu_0_eq, g = RED_init(observation,
                                            observation_value[:],\
                                            observation_type,\
                                            P_or_gamma=P_or_gamma,\
                                            P_tot = P_obs[:],\
                                            PFT_competition=PFT_competition,\
                                            overwrite_check = overwrite_check,\
                                            output_key=output_key)

### Get PFT values
key = ('PFT_name','I','J','J_max','m','alpha')
PFT_name, I, J, J_max, m, alpha = LOAD_PFT_VALUES(key)

# Determine the equilibrium growth of an individual PFT member 
p = np.zeros((I,J_max))

for i in range(0,I):
    
    p[i,:] = g[i,:]/(1.0-alpha[i])

#%% Drive RED over time
dt = 0.5 # Per year
start_year = 0.
end_year = 1000.
year = np.arange(start_year,end_year,dt)
K = len(year)

N = np.zeros((I,J_max,K))
nu_tot = np.zeros((I,K))
M_tot = np.zeros((I,K))
lambda_dem = np.zeros((I,K))
mu_0 = np.zeros((I,K))
F_N = np.zeros((I,J_max,K))
P = np.zeros((I,K))
m_lower = np.ones(I) * 1


output_key = ('N_m','nu_tot','mu_0')
mu_0 = np.zeros((I,K))

#N[:,:,0] = N_eq_m[:,:]  # Uncomment for equilibrium initialisation
#nu_tot[:,0]  = nu_eq[:]

for k in range(0,K-1):
    
    for i in range(0,I):
        
        P[i,k] = sum(N[i,0:J[i],k]*p[i,0:J[i]]) # Determine grid-box assimilate
        
    N[:,:,k+1], nu_tot[:,k+1], mu_0[:,k+1] = \
    \
    RED(dt,P[:,k],N[:,:,k],nu_tot[:,k],gamma_base=gamma_base[:],\
        output_key=output_key)
  
    
#%%

pft_col = ['#1f78b4','#a6cee3','#cab2d6','#33a02c','#b2df8a','#fdbf6f',\
                 '#ff7f00','#e31a1c','#fb9a99','#6e6d6c']                     # BET-Tr

pft_name = ['BET-Tr','BET-Te','BDT','NET','NDT','C3','C4','ESh','DSh']

plt.figure()
lw = 3.0
for i in [0,1,2,3,4,5,6,7,8]:
    
    plt.plot([year[0],year[-1]],[nu_eq[i],nu_eq[i]],color=pft_col[i],marker='x',ls=':',linewidth=2.0)
    plt.plot(year,nu_tot[i,:],label=pft_name[i],color=pft_col[i],linewidth=2.0)
    
plt.xlim(start_year,end_year)
plt.ylim(0,1)
plt.ylabel('coverage (-)')
plt.xlabel('time (year)')
plt.legend(bbox_to_anchor=(1.,0.8))


