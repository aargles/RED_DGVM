"""

Python 2.7, Encoding: UTF-8

Created on Mon Jan 24 2019

@author: A.Argles (aa760@exeter.ac.uk)
    
Description: Methods used in the Equilibrium Subroutine Libary.

"""

import numpy as np
from scipy.optimize import root


def find_mu_nu(J,mult,nu_eq,mu_init,alpha,phi_g,tol = None):
    """
    Uses scipy.optimise.root to get the mu value based upon the an initial 
    
    Inputs:
        
        - J, the number of mass classes.
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - nu_eq, the equilibrium vegetation fraction. (-)
        
        
        - mu_init, the "guessed" turnover ratio for the irretative process,
          defined as: gamma_init * m_init/g_init (-)
        
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
       
        - tol, tolerance for termination of the numerical root finding process.
          (-)
        
         
    Outputs:
        
        - mu, the numerically estimated turnover ratio, defined as: 
          gamma_init * m_init/g_init (-)
          
    Called: CONVERSION
          
    """
    
    args = (J,mult,nu_eq,alpha,phi_g)
    
    mu = root(Ya,mu_init,args = args, jac = dYa_dmu, tol = tol).x[0]
    
    return mu
    

def find_mu_g(P,J,mult,m_init, a_init, gamma_init, \
                    mu_init,alpha,phi_g,phi_a,tol = None):
    """
    Uses scipy.optimise.root to get the mu value based upon the an initial 
    
    Inputs:

        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
                  
        - J, the number of mass classes.
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - m_init, the initial sapling mass. (kg)
                
        - a_init, the initial sapling area. (m2)

        - gamma_init, the baseline mortalilty rate. (/year) 
        
        - mu_init, the "guessed" turnover ratio for the irretative process,
          defined as: gamma_init * m_init/g_init (-)
                 
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
       
        - tol, tolerance for termination of the numerical root finding process.
          (-)
        
         
    Outputs:
        
        - mu, the numerically estimated turnover ratio, defined as: 
          gamma_init * m_init/g_init (-)
          
          
    Called: CONVERSION
          
    """
    
    args = (P,J,mult,m_init,a_init,gamma_init,alpha,phi_g,phi_a)
    
    mu = root(Yb,mu_init,args = args, jac = dYb_dmu, tol = tol).x[0]
    

    return mu


    
def G_str_discrete(g_init,J,mult,a_init,nu_tot,mu,phi_g,phi_a):
    
    """
    Finding the steady state total growth from a vegetation fraction and intial
    growth rate.
    
    Inputs:
        
        - g_init, the initial sapling growth. (kg/year)
        
        - a_init, the initial crown area. (m2)
        
        - nu_tot, the total fractional coverage. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - phi_a, the crown area scaling parameter. (-)
    
    Outputs:
        
        - G_str, the total structral growth. (kg/m2/year)
          
    Called:  CONVERSION
    
    """      

    _, _, X_G, X_nu = X(J,mult,mu,phi_g=phi_g,phi_a=phi_a)
    

    G_str = g_init *(nu_tot*X_G)/(a_init*X_nu)
    
    
    return G_str
    
def g_init_discrete(G_str,J,mult,a_init,nu_tot,mu,phi_g,phi_a):
    
    """
    Finding the corisponding initial growth from the discrete equation for 
    the steady state.
    
    Inputs:
        
        - G_str, the total structral growth. (kg/m2/year)
        
        - a_init, the initial crown area. (m2)
        
        - nu_tot, the total fractional coverage. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - phi_a, the crown area scaling parameter. (-)
    
    Outputs:
        
        - g_init, the initial sapling growth. (kg/year)
          
    Called: CONVERSION
    
    """      

    _, _, X_G, X_nu = X(J,mult,mu,phi_g=phi_g,phi_a=phi_a)
    
    
    if nu_tot * X_G == 0:
        
        g_init = 0.0
        
    else:
        
        g_init = G_str*(a_init*X_nu)/(nu_tot*X_G)
    
    
    return g_init

    
def nu_eq_discrete(J,mult,mu,nu_min,alpha,phi_g):
    
    """
    Finds the discrete form of the equilibrium fraction from mu and alpha
    
    Inputs:
  
        - J, number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - nu_min, the minimum vegetation fraction allowed. (-)
             
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - alpha, the reseed fraction. (-)
        
    
    Outputs:
        
        - nu_eq, the equilibrium vegetation fraction of a PFT. (-)
          
    Called: CONVERSION
    
    """      

    _, X_N, X_G, _ = X(J,mult,mu,phi_g=phi_g)
    

    nu_eq = max(nu_min,1.0 - mu * (1.0-alpha)/(alpha) * X_N / X_G)
    
    
    return nu_eq
    
    
def number_dist_resolve(J,mult,m,m_init,N,G_seed,g_init,a_init,gamma_init,\
                        nu_tot,nu_eq,nu_min,phi_g):
        
    """
       Solves for the number density distribution of a PFT at equilibrium. RED
       at equilibrium is determined by a tri-diagonal matrix which can be 
       resolved based upon the intial condtions.
       
    Inputs:        
     
    - J, number of mass classes. (-)
    
    - mult, multiplicative bin scaling parameter. (-)
    
    - m, stored carbon biomass of a single tree. (kg)
    
    - m_init, the initial mass. (-)
    
    - N, the number density. (population/m2)
    
    - G_seed, the competitive seedling growth rates. (kg/m2/year)
    
    - a_init, the initial crown area. (m2)
    
    - g_init, the estimated initial growth rates (kg/year)
    
    - gamma_init, the baseline mortality. (/year)
        
    - nu_tot, total vegetation fraction.  (-)

    - nu_eq, equilibrium vegetation fraction.  (-)
    
    - nu_min, the minimum vegetation fraction allowed. (-)
    
    - phi_g, growth-mass scaling power. (-)
    
      Outputs:
       
    - m, the mass range (kg)
    
    - N, the  number density (population/m2)
    
    - Ng_cont, total growth based upon the scaling of the allometry.
      (kg/m2/year)
    
    Called: CONVERSION
    
    """

    Fin = np.zeros(J)
            
    Fout = np.zeros(J)
    
    Ng_cont = 0.0

    for j in range(0,J):
        
        if j == 0:
            
            m[j] = m_init 

            Fin[j] = G_seed/m_init *  (1.0 - nu_eq)

            Fout[j] = g_init/((mult - 1.0)*m[j])            
            
            if (gamma_init != np.inf and Fin[j] > 0.):

                
                N[j] = Fin[j]/(Fout[j] + gamma_init)
                     
            
            else:

                N[j] = nu_min / a_init                        

        else:
            
            m[j] = m[j-1] * mult   
            
            g = g_init*(m[j]/m_init)**phi_g

            Fin[j] = N[j-1] * Fout[j-1]


            if j == J - 1:
                
                Fout[j] = 0
            
            else:
                
                Fout[j] = g/((mult - 1)*m[j])
            
           
            if (gamma_init != np.inf):
            
                N[j] = Fin[j]/(Fout[j] + gamma_init)
            
                
        Ng_cont = Ng_cont + (N[j]*(m[j]/m[0])**phi_g) 
            
        
                
    return m, N, Ng_cont
    

def Ya(mu,J,mult,nu_eq,alpha,phi_g):
    
    """
    Solving for the steady state discretised vegetation fraction, based upon 
    the assumed competition of random overlap for seedlings (seedlings only 
    grow in non-vegetative areas.), while also subtracting nu_eq.
    
    Inputs:
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - J, the number of mass classes.
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - nu_eq, the equilibrium vegetation fraction. (-)
        
        - alpha, the reseed fraction. (-)
        
        - phi_g, the growth mass scaling parameter. (-)
    
    Outputs:
        
        - nu_discrete, the discrete solution to a PFTs vegetation fraction. (-)
    
    Equations are derived in technical notes (needs reference)
    
    Called:  "find_mu_nu"
    
    """
    
    _, X_N, X_G, _ = X(J,mult,mu[0],phi_g=phi_g)
    
    nu_discrete = 1.0 - mu[0] * ((1.0 - alpha)/alpha) * (X_N/X_G) - nu_eq
    
    return nu_discrete
    

def dYa_dmu(mu,J,mult,nu_eq,alpha,phi_g):
    
    """
    Finds the discrete value for the gradient of vegetation fraction with mu.
    
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
   
    _, X_N, X_G,_ = X(J,mult,mu,phi_g=phi_g)
    
    _, dX_N_dmu, dX_G_dmu, _  = dX_dmu(J,mult,mu,phi_g=phi_g)
    
    
    dnu_dmu = - ((1-alpha)/alpha) * ((mu * dX_N_dmu + X_N)*X_G \
                                     - mu*X_N*dX_G_dmu)/X_G**2
    
    
    return dnu_dmu

    
def Yb(mu,P,J,mult,m_init,a_init,gamma_init,alpha,phi_g,phi_a):
    
    """
    Determines mu in relation to a productivity and disturbance rate.
    
    Inputs:
    
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)

        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
          
        - J, the number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
            
        - m_init, the initial mass. (-)
        
        - G_str, the total structral growth. (kg/m2/year)

        - a_init, corresponding initial crown area. (m2)    
        
        - gamma_init, the baseline mortality. (/year)
        
        - alpha, the reseed fraction (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - phi_a, the area mass scaling parameter. (-)
    
    
    Outputs:
        
        - Y, Function which balances nu_eq, P and G_str in the steady state 
             form. (-)
        
    """
    
    
    _, _, X_G, X_nu = X(J,mult,mu,phi_g=phi_g,phi_a=phi_a)
    
    
    Y = X_G / (mu*X_nu) - a_init * (1.0 - alpha) * P /(gamma_init * m_init)
    
    
    return Y
    
def dYb_dmu(mu,P,J,mult,m_init,a_init,gamma_init,alpha,phi_g,phi_a):
    
    """
    Determines mu in relation to a productivity and disturbance rate.
    
    Inputs:
    
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
          
        - J, the number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
            
        - m_init, the initial mass. (-)
        
        - a_init, corresponding initial crown area. (m2)    
        
        - gamma_init, the baseline mortality. (/year)
        
        - alpha, the reseed fraction (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - phi_a, the area mass scaling parameter. (-)
    
    
    Outputs:
        
        - dY_dmu, the gradient of the function which balances nu_eq, P and 
          G_str in the steady state form. (-)
        
    """
    
    
    
    _, _, X_G, X_nu = X(J,mult,mu,phi_g=phi_g,phi_a=phi_a)
    
    _, _, dX_G_dmu, dX_nu_dmu = dX_dmu(J,mult,mu,phi_g=phi_g,\
                                              phi_a=phi_a)
    
    dY_dmu = (mu * X_nu * dX_G_dmu - X_G * (X_nu + mu * dX_nu_dmu))/\
             (X_nu * mu)**2
    
    return dY_dmu
    
    
def X(J,mult,mu,phi_g = 0.0, phi_a = 0.0):
    
    """
    Function X is used to solve for the discretised summation of scaling
    of the variable distribution over mass classes, assuming scaled bin widths.
    
    Inputs:
        
        - J, the number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - phi_a, the area mass scaling parameter. (-)
        
        
    Outputs:

        - X_M, summation of the mass scaling over each mass class. (-)        

        - X_N, summation of the scaling of the number density over each class.
          (-)
         
        - X_G, summation of the scaling of the growth over each class. (-)
        
        
        - X_nu, summation of the vegetation fraction over each mass class. (-)
        
          
    Equations are derived in the technical notes (needs reference)
    
    Called: "G_str_discrete","g_init_discrete","Ya", "dYa_dmu"
        
    """
#todo

    X_M = 1.0
  
    X_N = 1.0

    X_G = 1.0    
    
    X_nu = 1.0
    
    for i in range(1,J):

        s_m = mult ** (i)
        
        s_g = mult ** (i*phi_g)                
        
        s_nu = mult ** (i*phi_a)
        
        prod = 1.0
        
        for j in range(1,i + 1):
            
            if j < J - 1:
                
                kj = k_j(j,mult,mu,phi_g)
            
            elif j == J - 1 :
                
                kj =  (mult**((j-1)*phi_g))/((mult**j - mult**(j-1))*mu)
                
            prod = prod * kj
        
        X_M = X_M + prod*s_m

        X_N = X_N + prod
        
        X_G = X_G + prod * s_g
        
        X_nu = X_nu + prod * s_nu

    return X_M, X_N, X_G, X_nu

    
def dX_dmu(J,mult,mu,phi_g=0.0,phi_a=0.0):
    
    """
    Returns the gradient of the summation of the scaling parameter changes 
    with mu.
    
    Inputs:
        
        - J, the number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - phi, the given variable mass scaling parameter, For Number density
          phi = 0.0. (-)
        
        
    Outputs:
        
        - dX_M_dmu, the gradient value of the summation of the scaling of the 
          biomass over each class with respect to mu. (-)
        
        - dX_N_dmu, the gradient value of the summation of the scaling of the 
          number density over each class with respect to mu. (-)
        
        - dX_G_dmu, the gradient value of the summation of the scaling of the 
          growth over each class with respect to mu. (-)
          
        - dX_nu_dmu, the gradient value of the summation of the scaling of the 
          fractional coverage over each class with respect to mu. (-)
          
          
    Equations are derived in the technical notes (needs reference)
    
    Called: "dnu_discrete_dmu"
        
    """
    
    dX_M_dmu = 0.0
    
    dX_N_dmu = 0.0
    
    dX_G_dmu = 0.0
    
    dX_nu_dmu = 0.0
    
    for i in range(1,J):
        
        prod = 1.0
        series = 0.0
        s_m = mult ** (i)
        s_g = mult ** (i*phi_g)
        s_nu = mult ** (i*phi_a)
        
        for j in range(1,i + 1):
                
            if j < J - 1:
                
                kj = k_j(j,mult,mu,phi_g)
                
                dkdmu =  dk_j_dmu(j,mult,mu,phi_g)

            elif j == J - 1 :
                
                kj =  mult**((j-1)*(phi_g-1))/(mu*(mult-1))

                dkdmu = - mult**((j-1)*(phi_g-1))/((mu**2)*(mult-1))
            
            prod = prod*kj
            
            series = series + dkdmu/kj

        dX_M_dmu = dX_M_dmu + series * prod * s_m
            
        dX_N_dmu = dX_N_dmu + series * prod
        
        dX_G_dmu = dX_G_dmu + series * prod * s_g
        
        dX_nu_dmu = dX_nu_dmu + series * prod * s_nu
    
    return dX_M_dmu,dX_N_dmu, dX_G_dmu,dX_nu_dmu

def k_j(j,mult,mu,phi_g):
    
    """
    Outputs the scaling charateristics of each mass class.
    
    Inputs: 
        
        - j, the given mass class. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
    Outputs:
        
        - k_j, the relative fraction of population in the current class. (-)
        
    Called: "X"

    """

    
    k_j = mult**((j-1)*(phi_g-1))/(mult**(j*(phi_g-1.))+mu*(mult-1.0))

    return k_j

def dk_j_dmu(j,mult,mu,phi_g):
    
    """
    Outputs the scaling charateristics gradient of each mass class against mu.
    
    Inputs: 
        
        - j, the given mass class. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
    Outputs:
        
        - dk_j_dmu, the relative change of fraction population in the current
          class against mu. (-)
        
    Called: "dX_dmu"
    """
    
    dk_j_dmu = -((mult - 1) * mult**(1 + j - (1 - j)*phi_g))\
               /(mu*(mult - 1)*mult**j + mult**(j*phi_g))**2

    return dk_j_dmu