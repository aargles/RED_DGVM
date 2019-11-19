"""

Python 2.7, Encoding: UTF-8

Created on Fri Apr 12 14:24:24 2019

@author: A.Argles (aa760@exeter.ac.uk)

Description: Contains the continous form of the equilbrium soultions.
time.
"""


import numpy as np
from scipy.special import gammainc
from scipy.special import gamma
from scipy.optimize import root


def find_mu_nu_continuous(mu_init,nu_obs,alpha,phi_g,tol = None):
    """
    Uses scipy.optimise.root to get the mu value based upon an observed 
    coverage.
    
    Inputs:
        
        - mu_init, the "guessed" turnover ratio for the irretative process,
          defined as: gamma_init * m_init/g_init (-)
        
        - nu_obs, the observed coverage of a PFT. (-)
        
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
       
        - tol, tolerance for termination of the numerical root finding process.
          (-)
        
         
    Outputs:
        
        - mu, the continuous estimate boundary turnover ratio, defined as: 
          gamma_init * m_init/g_init. (-)
                    
    """
    
    args = (alpha,phi_g,nu_obs)
    
    mu = root(Y,mu_init,args = args, jac = dY_dmu, tol = tol).x[0]
    
    return mu



    
#%% Secondary Functions

def Y(mu,alpha,phi_g,nu_obs):    
    """
    Solving for the steady state continuous vegetation fraction, based upon 
    the assumed competition of random overlap for seedlings (seedlings only 
    grow in non-vegetative areas.), while also subtracting nu_obs.
    
    Inputs:
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - J, the number of mass classes.
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - nu_eq, the equilibrium vegetation fraction. (-)
        
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
    
    Outputs:
        
        - Y, the difference between observed and equilibrium fractions. (-)
    
    Equations are derived in technical notes (needs reference)
    
    Called:  "find_mu_nu"
    
    """
    
    
    Y =  nu_eq_analy(mu,alpha,phi_g) - nu_obs
    

    return Y
    

def dY_dmu(mu,alpha,phi_g,nu_obs):
    
    """
    Finds the continuous value for the gradient of vegetation fraction with mu.
    
    Inputs:
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - J, the number of mass classes.
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - nu_eq, the equilibrium vegetation fraction. (-)
        
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
    Outputs:
        
        - dnu_dmu, discrete solution to the vegetation fraction 
          gradient with respect to mu. (-)
          
    Called: - "find_mu_nu"
    
    """    
       
    dY_dmu = nu_eq_analy_dmu(mu,alpha,phi_g)
        
    return dY_dmu


def nu_eq_analy(mu,alpha,phi_g):
    
    """
    Finds the analytical form of the equilibrium fraction from mu and alpha.
    
    Inputs:
          
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - alpha, the reseed fraction. (-)
        
    
    Outputs:
        
        - nu_eq, the equilibrium vegetation fraction of a PFT. (-)
          
    
    """          

    r_alpha = (1. - alpha)/alpha

    x = 1./(1. - phi_g)
    
    gam = Gam(x,x*mu)
    
    E = np.exp(x*mu)

    nu_eq = 1. - r_alpha * mu * (x*mu)**(x-1.0) /(E*gam)
    
    return nu_eq
    
    
    
def nu_eq_analy_dmu(mu,alpha,phi_g):
    """
    Finding the analytical value for the gradient of vegetation fraction based
    upon the analytical solutions.
    
    Inputs:
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
    
    Outputs:
        
        - dnu_anly_dmu, analytical solution to the vegetation fraction 
          gradient. (-)
         
    
    """    
    
    r_alpha = (1. - alpha)/alpha

    x = 1./(1. - phi_g)
    
    gam = Gam(x,x*mu)
    
    E = np.exp(x*mu)

    numerator = E * gam * ((x*mu)**x/mu) - \
                    mu *(x*mu)**(x-1.0) * (x * gam * E - (mu * x)**x/mu) 
    
    denominator = (E * gam ) **2
    
    dnu_anly_dmu = - r_alpha *  numerator / denominator
    
    
    
    return dnu_anly_dmu


def M_tot_eq_analy(mu,nu_eq,m_init,a_init,alpha,phi_g):
    """
    Finding the analytical value for the total biomass.
    
    Inputs:
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - nu_eq, the equilibirum fraction. (-)
        
        - m_init, the mass of a sapling. (kg)
        
        - a_init, the crown area of a sapling. (m^2)
        
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
    
    Outputs:
        
        - M_tot, the toal vegetation carbon from the equilbirum fraction. 
                 (kg m^2)
          
    
    """    
    
    x = 1.0/(1.0 - phi_g)
    
    
    gam = Gam(x + 1.0 , mu * x)
        
    N_eq = N_tot_eq_analy(mu,nu_eq,a_init,alpha,phi_g)
    
    M_eq = m_init * N_eq *(mu*x)**(-x) * np.exp(mu*x) * gam


    return M_eq
    

def N_tot_eq_analy(mu,nu_eq,a_init,alpha,phi_g):
    """
    Finding the analytical value for the total stand density.
    
    Inputs:
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - nu_eq, the equilibirum fraction. (-)
                
        - a_init, the crown area of a sapling. (m^2)
        
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
    
    Outputs:
        
        - N_tot, the toal stand density from the equilbirum fraction. 
                 (kg m^2)
          
    
    """    
    
    x = 1.0/(1.0 - phi_g)
    gam = Gam(phi_g*x+1,mu*x)
        
    N_eq = nu_eq* (x*mu)**(x) /(a_init * np.exp(x*mu)*gam) 
    
    return N_eq
    
    
def Gam(s,x):
    """
    Computes the incomplete gamma function from the scipy functions. 
    
    Inputs:
        
        Y = Integral[t^(s - 1) exp(-t) dt]{0 -> x}
        
        - s, power of the gamma frunction integral Y. (-)
        
        - x, the upperbound to the integral Y. (-)
        
    Outputs:
        
        - Y, the incomplete gamma funtion. (-)
        
    """
    Y = gamma(s) * (1 - gammainc(s,x))
    
    return Y