"""

Python 2.7, Encoding: UTF-8

Created on Mon Jan 24 2019

@author: A.Argles (aa760@exeter.ac.uk)
    
Description: Methods used in the Equilibrium Subroutine Libary.

"""

import numpy as np
from .MISC_METHODS import newton
from scipy.optimize import root

def find_mu(observation,observation_value,J,S,m_0,m_scaling,h_0,\
                h_scaling,a_0,a_scaling,g_scaling,alpha,mu_init=0.1,\
                continuous_check=False):
    """
    Uses the newton algorithim to get the mu value required to fit the
    observations based upon an initial guess.
    
    Inputs:
        
                 
        - observation, can be "number density/N_obs", "coverage"/"nu_obs" (-), 
          "carbon mass"/"M_obs" (kgC), "height"/"H_obs" (m), 
          "basal area"/"A_obs" (m2), "assimilate mortality"/"P_gamma". (-) 
         
        - observation_value, the PFT observed value of a total. (-)
        
        - J, the number of mass classes.
        
        - S, multiplicative bin scaling parameter. (-)
                                   
        - m_0, carbon mass at boundary. (kgC)
            
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - h_0, boundary height. (m)
            
        - h_scaling, the height scaling with respect to h_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
                        
        - alpha, the reseed fraction. (-)
        
        - mu_init, the "guessed" turnover ratio for the irretative process,
          defined as: gamma_init * m_init/g_init (-)
                  
    Outputs:
        
        - mu_0, the numerically estimated turnover ratio, defined as: 
          gamma_init * m_init/g_init (-)
          
          
    """

    if observation_value == 0:
        
        mu_0 = float('inf')

    else:
        
        args = (observation,observation_value,J,S,m_0,m_scaling,h_0,\
                h_scaling,a_0,a_scaling,g_scaling,alpha)
        
        mu_0 = newton(diffobs_tot,ddiffobs_tot_dmu,args,mu_init)

    return mu_0
    

def find_mu_P(P,gamma_base,J,S,m_0,m_scaling,a_0,a_scaling,g_scaling,alpha,\
              mu_init=0.1):
    """
    Uses the newton algorithim to get the mu value required to fit growth and
    mortality observations based upon an initial guess.
    
    Inputs:
        
        - P, Net Primary Productivity minus the Local Litterfall. The total 
             carbon assimilate for the gridbox. (kg/m2/year)
 
        - gamma_base, baseline mortalilty. (population/year)
        
        - J, the number of mass classes.
        
        - S, multiplicative bin scaling parameter. (-)
                                   
        - m_0, carbon mass at boundary. (kgC)
            
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
                        
        - alpha, the reseed fraction. (-)
        
        - mu_init, the "guessed" turnover ratio for the irretative process,
          defined as: gamma_init * m_init/g_init (-)
                  
    Outputs:
        
        - mu_0, the numerically estimated turnover ratio, defined as: 
          gamma_init * m_init/g_init (-)
          
          
    """

    if P <= 0:
        
        mu_0 = float('inf')

    else:
        
        args = (P,gamma_base,J,S,m_0,m_scaling,a_0,a_scaling,g_scaling,alpha)
        
        mu_0 = newton(diffP_tot,ddiffP_tot_dmu,args,mu_init)
        
        if mu_0 != mu_0:
            
            mu_0 = float('inf')

    return mu_0    

    
    
def find_mu_sample(observation,observation_value,J,j_lower,j_upper,S,m_scaling,\
                   g_scaling,mu_init=0.3):
    """
    Uses the newton algorithim to get the mu value required to fit the
    observation sample based upon the an initial guess. 
    
    Inputs:
         
        - observation, can be "number density/N_obs", "coverage"/"nu_obs" (-), 
          "carbon mass"/"M_obs" (kgC), "height"/"H_obs" (m), 
          "basal area"/"A_obs" (m2), "assimilate mortality"/"P_gamma". (-) 
         
        - observation_value, the sample distribution. (-)
        
        - J, the number of mass classes.       
        
        - j_lower, the lower boundary class of the sample. (-) 
        
        - j_upper, the higher boundary class of the sample. (-)
            
        - S, multiplicative bin scaling parameter. (-)
        
        - m_scaling, the mass scaling with respect to m_0. (-)

        - g_scaling, the growth scaling with respect to g_0. (-)
        
        - mu_init, the "guessed" turnover ratio for the irretative process,
          defined as: gamma_init * m_init/g_init (-)
 
    Outputs:
        
        - mu_0, the numerically estimated turnover ratio, defined as: 
          gamma_init * m_init/g_init (-)
          
    """
    
    if j_lower < j_upper:
        
        args = (observation,observation_value,J,j_lower,j_upper,S,m_scaling,\
                   g_scaling)
    
        mu_0 = newton(diffsample,ddiffsample_dmu,args,mu_init)

    else:
        
        mu_0 = float('nan')

    return mu_0
    

    
def nu_eq_disc(mu_0,J,S,m_scaling,h_scaling,a_scaling,g_scaling,alpha,nu_min):
    
    """
    Finds the discrete form of the equilibrium fraction from mu and alpha
    
    Inputs:
  
        - mu_0, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - J, the number of mass classes.       
        
        - S, multiplicative bin scaling parameter. (-)
                                               
        - m_scaling, the mass scaling with respect to m_0. (-)
                   
        - h_scaling, the height scaling with respect to h_0. (-)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
                        
        - alpha, the reseed fraction. (-)
        
        - nu_min, minimum vegetation fraction. (-)
               
    Outputs:
        
        - nu_eq, the equilibrium vegetation fraction of a PFT. (-)
          
    
    """      

    _, X_N, X_G, _, _ = X(mu_0,J,S,m_scaling,m_scaling,m_scaling,g_scaling)

    nu_eq = max(nu_min,1.0 - mu_0 * (1.0-alpha)/(alpha) * X_N / X_G)
    
    
    return nu_eq

def nu_eq_star_disc(I,nu_eq,nu_min,C,c_pftg,c_unique,P_or_gamma=None,\
                    P_tot=None,gamma_base=None):
    
    """
    Adjusts the equilibrium to take into account the PFT competition coefficent
    , returning it as "nu_eq_star".
    
    Inputs:
        
        - I, Number of PFTs. (-)
        
        - nu_eq, the equilibrium coverage of a PFT existing wihtin a 
          monoculture. (-)
          
        - nu_min, minimum vegetation fraction. (-)
          
        - C, the number of unique PFT groups. (-)
                  
        - c_pftg, the PFT group dominance hierarchy. (-)
        
        - c_unique, the unique rows of each of the PFT group hierarchy. (-)
        
        - P_or_gamma, method used to employ in splitting mu into it's 
          components, either by knowing a value for P ('P') or a value for
          gamma ('gamma'). (-)
          
        - P_tot, Net Primary Productivity minus the Local Litterfall. The total 
          PFT carbon assimilate for the gridbox. (kg/m2/year)
             
        - gamma_base, the baseline mortalility from the intialisation. 
          (population/year)
          
    Outputs:
        
        - nu_eq_star, the equilibrium coverage that has been adjusted with
          regard to the competitive regime. (-)
          
        - nu_eq_shade, the steady state of the coverage, once PFT shading has
          been calculated. (-)
              
    """
    # Find the number of unique groups
    nu_eq_shade = np.zeros(I)
    nu_eq_star = np.zeros(I)
    nu_shade_group = np.zeros(C)
    nu_above = 0
    groupmax = np.zeros((C,I),dtype=bool)
    zerogrowth = np.zeros(I,dtype=bool)
    groupmax_count = np.zeros(C)
    
    for c in range(0,C):
        
        nu_eq_max = nu_min
        count = len(c_pftg[c,c_pftg[c,:]==True])
        
        for i in range(0,I):

            # If the input growth is zero or input mortality is infinte then
            # the equilibrium fraction must be zero.
            if P_or_gamma == 'P':
                
                if P_tot[i] == 0:
                    
                    zerogrowth[i] = True

            elif P_or_gamma == 'gamma':
                
                if gamma_base[i] == np.inf:
                    
                    zerogrowth[i] = True
            
            if zerogrowth[i] == True:
                
                nu_eq[i] = nu_min
                
            nu_shade_group[c] = nu_shade_group[c] + c_unique[c,i] * nu_eq[i] 

            if (nu_eq_max < nu_eq[i] and c_pftg[c,i] == True):
                
                groupmax[c,groupmax[c,:]] = False
                groupmax[c,i] = True
                nu_eq_max = nu_eq[i]
                groupmax_count[c] = 1.
            
            elif (nu_eq_max == nu_eq[i] and c_pftg[c,i] == True):

                groupmax[c,i] = True
                groupmax_count[c] = groupmax_count[c] + 1.
                
            if nu_shade_group[c] > 1.0 - (I-1)*nu_min:
                
                nu_shade_group[c] = 1.0 - (I-1)*nu_min

        nu_shade_group[c] = nu_shade_group[c] + nu_min*(count-1)
        nu_eq_star[c_pftg[c,:]] = 0.5*(nu_shade_group[c] + nu_eq[c_pftg[c,:]]) 
        nu_eq_star[groupmax[c,:]] = max(nu_min,nu_shade_group[c])                   
        nu_eq_shade[c_pftg[c,:]] = nu_min
        nu_eq_shade[groupmax[c,:]] = max(nu_min,(nu_shade_group[c] -\
                      nu_above - nu_min*(count-groupmax_count[c]))\
                     /groupmax_count[c])
                
        nu_above = sum(nu_eq_shade[:])
        
    return nu_eq_star, nu_eq_shade
   

    
    
def diffsample(mu_0,observation,observation_value,J,j_lower,j_upper,S,\
               m_scaling,g_scaling):
    
    """
    Finds the fractional difference between the estimated equilibrium 
    distribution and that of an observed distribution.
            
        - mu_0, the turnover ratio defined as: gamma_base * m_0/g_0. (-)
            
        - observation, can be "number density/N_obs", "coverage"/"nu_obs" (-), 
          "carbon mass"/"M_obs" (kgC), "height"/"H_obs" (m), 
          "basal area"/"A_obs" (m2), "assimilate mortality"/"P_gamma". (-) 
         
        - observation_value, the sample distribution. (-)
        
        - J, the number of mass classes.       
        
        - j_lower, the lower boundary class of the sample. (-) 
        
        - j_upper, the higher boundary class of the sample. (-)
        
        - S, multiplicative bin scaling parameter. (-)
        
        - m_scaling, the mass scaling with respect to m_0. (-)

        - g_scaling, the growth scaling with respect to g_0. (-)
        
        - mu_init, the "guessed" turnover ratio for the irretative process,
          defined as: gamma_init * m_init/g_init (-)
    
    Outputs:
        
        - diffobs, the difference between eq(mu) and the observations. (-)
                 
    """
    
    if (observation is 'stand density' or observation is 'N_obs'):
            
        prod = 1
        N_a = observation_value[j_lower]
        N_b = observation_value[j_upper]

        if N_a == 0:
            
            return float('inf')
            
        else:
            
            ratio = N_b / N_a
                    
        for j in range(j_lower+1,j_upper+1):
        
            if j < J-1:
                frac_in = g_scaling[j-1]/m_scaling[j-1]
                frac_out = g_scaling[j]/m_scaling[j]
                kj = frac_in/(frac_out + mu_0*(S-1))
                
            elif j == J-1:
                
                frac_in = g_scaling[j-1]/m_scaling[j-1]
                kj = frac_in/(mu_0*(S-1))
                               
            prod = prod * kj
        return ratio - prod        
        
    return
    
    
def ddiffsample_dmu(mu_0,observation,observation_value,J,j_lower,j_upper,S,\
               m_scaling,g_scaling):
        
    """
    Finds the gradient between the fractional ratio between the estimated 
    equilibrium distribution and that of an observed distribution.
            
        - mu_0, the turnover ratio defined as: gamma_base * m_0/g_0. (-)
            
        - observation, can be "number density/N_obs", "coverage"/"nu_obs" (-), 
          "carbon mass"/"M_obs" (kgC), "height"/"H_obs" (m), 
          "basal area"/"A_obs" (m2), "assimilate mortality"/"P_gamma". (-) 
         
        - observation_value, the sample distribution. (-)
        
        - J, the number of mass classes.       
        
        - j_lower, the lower boundary class of the sample. (-) 
        
        - j_upper, the higher boundary class of the sample. (-)
        
        - S, multiplicative bin scaling parameter. (-)
        
        - m_scaling, the mass scaling with respect to m_0. (-)

        - g_scaling, the growth scaling with respect to g_0. (-)
        
        - mu_init, the "guessed" turnover ratio for the irretative process,
          defined as: gamma_init * m_init/g_init (-)
    
    Outputs:
        
        - diffobs, the difference between eq(mu) and the observations. (-)
                 
    """
     
    if (observation is 'stand density' or observation is 'N_obs'):
        
        prod = 1
        series = 0
        
        for j in range(j_lower+1,j_upper+1):
        
            if j < J - 1:
                
                frac_in = g_scaling[j-1]/m_scaling[j-1]
                frac_out = g_scaling[j]/m_scaling[j]
                kj = frac_in/(frac_out + mu_0*(S-1))
                dkjdmu = - (S-1) * frac_in/(frac_out+mu_0*(S-1))**2
                    
            elif j == J - 1 :
                
                frac_in = g_scaling[j-1]/m_scaling[j-1]
                kj = frac_in/(mu_0*(S-1))
                dkjdmu = - kj/mu_0
            
            prod = prod*kj
            series = series + dkjdmu/kj
        
        return - prod * series
        
    
    
    
def diffobs_tot(mu_0,observation,observation_value,J,S,m_0,m_scaling,h_0,\
                h_scaling,a_0,a_scaling,g_scaling,alpha):
    
    """
    Function returns the difference between a mu_0 value and an observation of
    a total gridbox quantity.
    
    Inputs:
        
        - mu_0, the turnover ratio defined as: gamma_base * m_0/g_0. (-)
           
        - observation, can be "coverage"/"nu_obs" (-), "stand density"/"N_obs"
         (/m2) "carbon mass"/"M_obs" (kgC/m2), "height"/"H_obs" (m)
         
        - observation_value, the PFT observed value of a total. (-)
        
        - J, the number of mass classes.
        
        - S, multiplicative bin scaling parameter. (-)
                                   
        - m_0, carbon mass at boundary. (kgC)
            
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - h_0, boundary height. (m)
            
        - h_scaling, the height scaling with respect to h_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
                        
        - alpha, the reseed fraction. (-)
    
    Outputs:
        
        - diffobs, the difference between eq(mu) and the observations. (-)
       
    """
    
    X_M, X_N, X_G, X_nu, X_H = X(mu_0,J,S,m_scaling,h_scaling,a_scaling,\
                            g_scaling)

    nu_eq = 1.0 - mu_0 * (1.0-alpha)/alpha * (X_N/X_G)
    
       
    if (observation == 'coverage' or observation == 'nu_obs'):
        
        return nu_eq - observation_value
        
    else:
        
        N_0 = nu_eq/(a_0 * X_nu)
        
    if (observation == 'stand density' or observation == 'N_obs'):
        
        return N_0 * X_N - observation_value
        
    elif (observation == 'carbon mass' or observation == 'M_obs'):

        return N_0 * m_0 * X_M - observation_value
        
    elif (observation == 'height' or observation == 'H_obs'):
        
        return N_0 * h_0 * X_H - observation_value

    else:

        raise UserWarning('\n Error in calculating the difference'+\
        ' between the  numerical equilibrium and observations: \n'+\
        ' observation not defined')        
        
    return
    
    

def ddiffobs_tot_dmu(mu_0,observation,observation_value,J,S,m_0,m_scaling,h_0,\
                h_scaling,a_0,a_scaling,g_scaling,alpha):
    
    """
    Finds the gradient between the dicrete equilibrium and a observation
    , with respect to mu. 
   
    Inputs:
        
        - mu_0, the turnover ratio defined as: gamma_init * m_init/g_init (-)
           
        - observation, can be "coverage"/"nu_obs" (-), "stand density"/"N_obs"
         (/m2) "carbon mass"/"M_obs" (kgC/m2), "height"/"H_obs" (m)
         
        - observation_value, the PFT observed value of a total. (-)
        
        - J, the number of mass classes.
        
        - S, multiplicative bin scaling parameter. (-)
                                   
        - m_0, carbon mass at boundary. (kgC)
            
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - h_0, boundary height. (m)
            
        - h_scaling, the height scaling with respect to h_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
                        
        - alpha, the reseed fraction. (-)
    
    Outputs:
        
        - diffobs, the difference between eq(mu) and the observations. (-)
    """    
   
    X_M, X_N, X_G, X_nu, X_H = X(mu_0,J,S,m_scaling,h_scaling,a_scaling,\
                                 g_scaling)
    
    dX_M_dmu,dX_N_dmu,dX_G_dmu,dX_nu_dmu,dX_H_dmu  = dX_dmu(mu_0,J,S,
                                                            m_scaling,\
                                                            h_scaling,\
                                                            a_scaling,\
                                                            g_scaling)
    
    dnu_eq_dmu = -(1.0-alpha)/alpha * ((X_N/X_G) \
                   + mu_0* (X_G*dX_N_dmu - X_N*dX_G_dmu)/X_G**2)
    
    if (observation == 'coverage' or observation == 'nu_obs'):
        
        return dnu_eq_dmu
        
    else:
        
        nu_eq = 1.0 - mu_0 * ((1.0-alpha)/alpha) * (X_N/X_G)
        N_0 = nu_eq/(a_0 * X_nu)
        dN_0_dmu = (a_0*X_nu * dnu_eq_dmu - nu_eq*a_0*dX_nu_dmu)/(a_0*X_nu)**2
        
        
    if (observation == 'stand density' or observation == 'N_obs'):
        
        return N_0 * dX_N_dmu + X_N *dN_0_dmu
        
    elif (observation == 'carbon mass' or observation == 'M_obs'):
        
        return  N_0 *m_0*dX_M_dmu + m_0*X_M *dN_0_dmu
        
    elif (observation == 'height' or observation == 'H_obs'):
        
        return N_0 *h_0*dX_H_dmu + h_0*X_H *dN_0_dmu

    else:

        raise UserWarning('\n Error in calculating the gradient in'+\
        ' difference between the  numerical equilibrium and observations: \n'+\
        ' observation not defined')     
        
        
    return 


def diffP_tot(mu_0, P,gamma_base,J,S,m_0,m_scaling,a_0,a_scaling,g_scaling,\
              alpha):
    
    """
    Function returns the difference between a mu_0 P and an observation of
    P, given a known mortality rate. 
    
    Inputs:
        
        - mu_0, the turnover ratio defined as: gamma_init * m_init/g_init (-)
           
        - P, The observed assimilate rate.  Net Primary Productivity minus the 
             Local Litterfall. The total carbon assimilate for the gridbox. 
             (kg/m2/year)
         
        - gamma_base, the baseline mortalilty. (population/year)
        
        - J, the number of mass classes.
        
        - S, multiplicative bin scaling parameter. (-)
                                   
        - m_0, carbon mass at boundary. (kgC)
            
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
                        
        - alpha, the reseed fraction. (-)
    
    Outputs:
        
        - diffP, the difference between the mu assimilate and the observerd 
          assimilate. (kgC/m2/year)
    """
    # mu_0 = gamma_base * m_0 / g
    g_0 = gamma_base * m_0 / mu_0
    _, X_N, X_G, X_nu, _ = X(mu_0,J,S,m_scaling,m_scaling,a_scaling,g_scaling)
    nu_eq = 1.0 - mu_0 * (1.0-alpha)/alpha * (X_N/X_G)
    N_0 = nu_eq/(a_0 * X_nu)
    
    return N_0 * g_0 * X_G - P


def ddiffP_tot_dmu(mu_0, P,gamma_base,J,S,m_0,m_scaling,a_0,a_scaling,\
                  g_scaling,alpha):
    
    """
    Function returns the difference between a mu_0 P and an observation of
    P, given a known mortality rate. 
    
    Inputs:
        
        - mu_0, the turnover ratio defined as: gamma_init * m_init/g_init (-)
           
        - P, The observed assimilate rate.  Net Primary Productivity minus the 
             Local Litterfall. The total carbon assimilate for the gridbox. 
             (kg/m2/year)
         
        - gamma_base, the baseline mortalilty. (population/year)
        
        - J, the number of mass classes.
        
        - S, multiplicative bin scaling parameter. (-)
                                   
        - m_0, carbon mass at boundary. (kgC)
            
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
                        
        - alpha, the reseed fraction. (-)
    
    Outputs:
        
        - diffP, the difference between the mu assimilate and the observerd 
          assimilate. (kgC/m2/year)
    """
     # mu_0 = gamma_base * m_0 / g
    g_0 = gamma_base * m_0 / mu_0  
    _, X_N, X_G, X_nu, _ = X(mu_0,J,S,m_scaling,m_scaling,a_scaling,g_scaling)
    _, dX_N_dmu,dX_G_dmu,dX_nu_dmu,_  = dX_dmu(mu_0,J,S,m_scaling,m_scaling,\
                                                        a_scaling,g_scaling)
    nu_eq = 1.0 - mu_0 * ((1.0-alpha)/alpha) * (X_N/X_G)
    dnu_eq_dmu = -(1.0-alpha)/alpha * ((X_N/X_G) \
               + mu_0* (X_G*dX_N_dmu - X_N*dX_G_dmu)/X_G**2)
    N_0 = nu_eq/(a_0 * X_nu)
    dN_0_dmu = (a_0*X_nu * dnu_eq_dmu - nu_eq*a_0*dX_nu_dmu)/(a_0*X_nu)**2
    
    return N_0 *g_0*dX_G_dmu + g_0*X_G *dN_0_dmu


def N_dist_eq(mu_0,J,S,m_scaling,a_scaling,a_0,g_scaling,nu_shade_eq):
    
    """
    Function to find the equilibrium number density distribution, from a mu,
    and equilibrium shading fractions.
    
    Inputs:
                
        - mu_0, the turnover ratio defined as: gamma_init * m_init/g_init (-)
          
        - J, the number of mass classes.
        
        - S, multiplicative bin scaling parameter. (-)
        
        - m_scaling, the mass scaling with respect to m_0. (-)
                
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - g_scaling, the growth scaling with respect to g_0. (-)
        
        - nu_eq_shade, the steady state of the coverage, once PFT shading has
          been calculated. (-)
          
    Outputs:
        
        - N, the RED equilibrium number density (population/m2)

    """
    
    N = np.zeros(J)
    
    _, _, _, X_nu, _ = X(mu_0,J,S,m_scaling,m_scaling,a_scaling,g_scaling)
    
    N[0] = nu_shade_eq/(a_0*X_nu)
    
    for j in range(1,J):
        
        if j < J-1:
            frac_in = g_scaling[j-1]/m_scaling[j-1]
            frac_out = g_scaling[j]/m_scaling[j]
            kj = frac_in/(frac_out + mu_0*(S-1))
            
        elif j == J-1:
            
            frac_in = g_scaling[j-1]/m_scaling[j-1]
            kj = frac_in/(mu_0*(S-1))
            
        N[j] = N[j-1] * kj
        
    
    return N
    
    
def X(mu_0,J,S,m_scaling,h_scaling,a_scaling,g_scaling):
    
    """
    Function X is used to solve for the discretised summation of scaling
    of the variable distribution over mass classes, assuming scaled bin widths.
    
    Inputs:
        
        - mu_0, the turnover ratio defined as: gamma_init * m_init/g_init (-)
          
        - J, the number of mass classes.
        
        - S, multiplicative bin scaling parameter. (-)
        
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - h_scaling, the height scaling with respect to h_0. (-)
        
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
        
        
    Outputs:

        - X_M, summation of the mass scaling over each mass class. (-)        

        - X_N, summation of the scaling of the number density over each class.
          (-)
         
        - X_G, summation of the scaling of the growth over each class. (-)
        
        - X_nu, summation of the scaling of the vegetation fraction 
          over each mass class. (-)
          
        - X_H, summation of the scaling of the height over each mass clsass. 
          (-)
        
    """
    X_M = 1.0
    X_N = 1.0
    X_G = 1.0    
    X_nu = 1.0
    X_H = 1.0
    
    for i in range(1,J):

        prod = 1.0
        
        for j in range(1,i + 1):
            
            if j < J - 1:
                
                frac_in = g_scaling[j-1]/m_scaling[j-1]
                frac_out = g_scaling[j]/m_scaling[j]
                kj = frac_in/(frac_out + mu_0*(S-1))
            
            elif j == J - 1 :
                
                frac_in = g_scaling[j-1]/m_scaling[j-1]
                kj = frac_in/(mu_0*(S-1))
                
            prod = prod * kj
        
        X_M = X_M + prod*m_scaling[i]
        X_N = X_N + prod
        X_G = X_G + prod*g_scaling[i]
        X_nu = X_nu + prod *a_scaling[i]
        X_H = X_H + prod * h_scaling[i]

    return X_M, X_N, X_G, X_nu, X_H

    
def dX_dmu(mu_0,J,S,m_scaling,h_scaling,a_scaling,g_scaling):
    
    """
    Returns the gradient of the summation of the scaling parameter changes 
    with mu.
    
    Inputs:
        
        - J, the number of mass classes. (-)
        
        - S, multiplicative bin scaling parameter. (-)
        
        - mu_0, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
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
        
        - dX_H_dmu, the gradient value of the summation of the scaling of the 
          height over each class with respect to mu. (-) 
                  
    """
    
    dX_M_dmu = 0.0
    dX_N_dmu = 0.0
    dX_G_dmu = 0.0
    dX_nu_dmu = 0.0
    dX_H_dmu = 0.0
    
    for i in range(1,J):
        
        prod = 1.0
        series = 0.0
        
        for j in range(1,i + 1):
                
            if j < J - 1:
                
                frac_in = g_scaling[j-1]/m_scaling[j-1]
                frac_out = g_scaling[j]/m_scaling[j]
                kj = frac_in/(frac_out + mu_0*(S-1))
                dkjdmu = - (S-1) * frac_in/(frac_out+mu_0*(S-1))**2
                    
            elif j == J - 1 :
                
                frac_in = g_scaling[j-1]/m_scaling[j-1]
                kj = frac_in/(mu_0*(S-1))
                dkjdmu = - kj/mu_0
            
            prod = prod*kj
            
            series = series + dkjdmu/kj

        dX_M_dmu = dX_M_dmu + series * prod * m_scaling[i]
        dX_N_dmu = dX_N_dmu + series * prod
        dX_G_dmu = dX_G_dmu + series * prod * g_scaling[i]
        dX_nu_dmu = dX_nu_dmu + series * prod * a_scaling[i]
        dX_H_dmu = dX_nu_dmu + series * prod * h_scaling[i]
    
    return dX_M_dmu,dX_N_dmu,dX_G_dmu,dX_nu_dmu,dX_H_dmu