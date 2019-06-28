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

    nu_eq = 1. - r_alpha * mu * (x * mu) ** (x - 1.0) / (np.exp(x * mu) \
                                 * Gam(x,x*mu))
    
    return nu_eq
    


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