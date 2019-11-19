"""

Python 2.7, Encoding: UTF-8

Created on Mon Jan 24 2019

@author: A.Argles (aa760@exeter.ac.uk)
    
Description: Contains the intilisation subroutines called in the main script.

TODO: 
Also contains a module to provide snapshots of the equilibrium states.

"""


from .EQUILIBRIUM_METHODS import * 
from .MISC_METHODS import LOAD_PFT_VALUES, output_key_interpet

import numpy as np


def EQ_get_outputs(output_key,I,J,mu_0,N_eq_m,g_scaling,alpha,gamma_base=None,\
                   P_eq = None):
    
    """
    Sorts the Outputs based upon a requested key.

    Inputs:
        
        - output variables listed in the output key tuple:  
          key, the string containing the output information. 'X_tot', the total
          value of X across all classes. 'X_m', the distribution of X across 
          the mass classes. 'X_mean', the mean value of X across the 
          distribution of trees. 'X_Y_r', the value of X at the r percentile of
          Y. 'X_norm' returns the normalised distribution. Other values can be:
          'mu_0', 'alpha' and anything stored in "PFT_STORE".
    
     - I, Number of PFTs (-)
        
     - J, number of mass classes. (-)      
          
     - mu_0, the turnover ratio a boundary tree defined as: \
          gamma_base * m_0/g_0. (-)
     
     - N_eq_m, the number of trees per unit area over mass at different 
            PFTs, at equilibrium. (population/m2)
            
     - nu_eq_m, the equilibrium coverage of each PFT across mass classes. (-) 

     - g_scaling, the growth-mass scaling. (-) ]
        
     - alpha, the fraction of tree productivity going into seedling 
          production. (-) 
            
     - gamma_base, the baseline mortalility. (population/year).
     
     - P_eq_shade, the total gridbox carbon assimilate at the equilibrium 
       state. (kgC/m2/year) - If P_or_gamma == 'P'
                  
    Returns:
        
        output_dir, dictonary object with the outputs and corresponding 
        strings.
    
    """
    m, h, a = LOAD_PFT_VALUES(('m','h','a'))      
    # compute the N distribution
    output_dict = {}
          

    for key in output_key:

        output_dict.update({key:output_key_interpet(key,I,J,mu_0,N_eq_m,m,a,h,\
                                                    g_scaling,alpha,\
                                                    gamma_base=gamma_base,\
                                                    P_tot=P_eq)})
        

    output_dict_subset = [output_dict[key] for key in output_key]

    return output_dict_subset
            
    
    
    
def EQ_mu(observation,observation_value,observation_type,P_or_gamma,I,J,J_max,\
          S,m_0,m_scaling,h_0,h_scaling,a_0,a_scaling,g_scaling,C,c_pftg,\
          c_unique,alpha,nu_min,PFT_competition,P=None,gamma_base=None,\
          continuous_check=True):
    
    """
    Uses two methods of getting the steady state involving a observation. 
    Returns the equilibrium coverage with regards to the PFT hierarchy.
    
    Inputs:
        
        - observation, can be "coverage"/"nu_obs" (-), "stand density"/"N_obs"
         (/m2) "carbon mass"/"M_obs" (kgC/m2), "height"/"H_obs" (m)
         
        - observation_value, the PFT observed value of a total. (-)
        
        - observation_type, the method employed by RED to find the mu_0 through
          the discrete solutions for the total, average or statistical 
          cummalitve (i.e. meadian or top 20% ). Currently this is just the
          total quantity. Values allowed: ("total","mean" or a tuple with 
          cummlative("cumulative")). (-)
              
              
        - P_or_gamma, method used to employ in splitting mu_0 into it's 
          components, either by knowing a value for P ('P') or a value for
          gamma ('gamma'). (-)
     
        - I, Number of PFTs. (-)
              
        - J, the number of mass classes. (-)
        
        - S, multiplicative bin scaling parameter. (-)
                                   
        - m_0, carbon mass at boundary. (kgC)
            
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - h_0, boundary height. (m)
            
        - h_scaling, the height scaling with respect to h_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)

        - C, the number of unique PFT groups. (-)
                  
        - c_pftg, the PFT group dominance hierarchy. (-)
        
        - c_unique, the unique rows of each of the PFT group hierarchy. (-)
                                          
        - alpha, the reseed fraction. (-)

        - nu_min, minimum vegetation fraction. (-)
        
        - PFT_competition, if True solve under the assumptions of the
          competitive exclusions of RED. Otherwise purely solve to fit the
          observations, therby ignoring the dynamical equilibrium. (-)
          
        - P, Net Primary Productivity minus the Local Litterfall. The total 
             carbon assimilate for the gridbox. (kg/m2/year)

        - gamma_base, the baseline mortalility from the intialisation. 
          (population/year)
             
        - continuous_check, if True solves for the equilibrium using the 
          continuous forms of the RED equations. (-)
        
    Outputs:
    
        
        - mu_0, the turnover ratio a boundary tree defined as: \
          gamma_base * m_0/g_0. (-)
        
        - N, the RED equilibrium number density across the classes
          (population/m2)
          
        - nu_eq_shade, the steady state of the coverage, once PFT shading has
          been calculated. (-) - PFT_competition == True
          
        - nu_eq, the equilibrium vegetation fraction of a PFT. (-) - 
          PFT_competition == False
          
        - P_eq_shade, the total gridbox carbon assimilate at the equilibrium 
          state. (kgC/m2/year) - If P_or_gamma == 'P'
              
    """
    nu_eq = np.zeros(I)
    mu_0 = np.zeros(I)
    N = np.zeros((I,J_max))
    #The total gridbox productivity will be different when compared to a fully
    # saturated grid-box. After adjusting to find the PFT competition gridbox
    # fraction, we must also adjust the total assimilate (P), accordingly. 
    # In this case we assume that the ajusted grid-box productivity is linearly
    # related to the gridbox coverage.
    
    if P_or_gamma == 'P':
    
        P_eq_shade = np.zeros(I)

    if observation_type == 'total':

        
        if (observation == 'coverage' or observation == 'nu_obs'):
            
            # If the observed coverage is equivalent to nu_eq, we derive the 
            # nu_eq_star and nu_eq_adjusted without having to find a mu_0 then 
            # coverage (in comparison to other observations).
            
            if PFT_competition == True:
                
                nu_eq[:] = observation_value[:]
                nu_eq[nu_eq<nu_min] = nu_min
                nu_eq_star,\
                nu_eq_shade = nu_eq_star_disc(I,nu_eq,nu_min,C,c_pftg,\
                                              c_unique,P_or_gamma=P_or_gamma,\
                                              P_tot=P,\
                                              gamma_base=gamma_base)
                
                for i in range(0,I):

                    j = np.arange(0,J[i],dtype=int)                     
                    mu_0[i] = find_mu('coverage',nu_eq_star[i],J[i],S[i],\
                                      m_0[i], m_scaling[i,j],h_0[i],\
                                      h_scaling[i,j],a_0[i],a_scaling[i,j],\
                                      g_scaling[i,j],alpha[i],mu_init=0.1,\
                                      continuous_check=continuous_check)

                    if P_or_gamma == 'P':
                        
                        if observation_value[i] > 0:
                        
                            P_eq_shade[i] = P[i]*nu_eq_shade[i]/\
                                            observation_value[i]
     
                        else:
                            
                            P_eq_shade[i] =  0.0


                                        
            else:
                
                nu_eq[:]  = observation_value[:]
                
                for i in range(0,I):
                    
                    j = np.arange(0,J[i],dtype=int)  
                    mu_0[i] = find_mu('coverage',nu_eq[i],J[i],S[i],\
                                      m_0[i], m_scaling[i,j],h_0[i],\
                                      h_scaling[i,j],a_0[i],a_scaling[i,j],\
                                      g_scaling[i,j],alpha[i],mu_init=0.1,\
                                      continuous_check=continuous_check)
                    
                    if P_or_gamma == 'P':
                        
                        if observation_value[i] != 0:
                        
                            P_eq_shade[i] = P[i]

                        else:
                            
                            P_eq_shade[i] = 0.


                       
        else:
            
            # Assuming the current observations have a close to identical
            # distribution to that of the equilibrium, we first need to find 
            # the equivalent equilibrium coverage that would lend to the 
            # observation.

            for i in range(0,I):
                
                j = np.arange(0,J[i],dtype=int)   
                mu_0[i] = find_mu(observation,observation_value[i],J[i],S[i],\
                                m_0[i],m_scaling[i,j],h_0[i],h_scaling[i,j],\
                                a_0[i],a_scaling[i,j],g_scaling[i,j],alpha[i],\
                                mu_init=0.1,continuous_check=continuous_check)
                nu_eq[i] = nu_eq_disc(mu_0[i],J[i],S[i],m_scaling[i,j],\
                                 h_scaling[i,j],a_scaling[i,j],g_scaling[i,j],\
                                 alpha[i],nu_min)

            # Next we take into account the cross PFT competitive dynamics, the 
            # equilibrium coverage is adjusted according to the model paper. 
            # Space currently occupied by sub-dominant PFTs is given over to 
            # the PFT with the highest current equilibrium coverage. While the 
            # adjusted equilibrium coverages are assumed to be somewhere 
            # between the summation within functional groups and the current
            # observation derived equilibrium coverage.
            
            if PFT_competition == True:
                
                nu_eq_star,\
                nu_eq_shade = nu_eq_star_disc(I,nu_eq,nu_min,C,c_pftg,c_unique,
                                              P_or_gamma=P_or_gamma,\
                                              P_tot=P,\
                                              gamma_base=gamma_base)

        # Refind mu_0, from the derived equilibrium cover. Then find the
        # equilbrium number density distribution.

        if P_or_gamma == 'gamma':
            
            P_eq_shade = np.zeros(I)
        
        if PFT_competition == True:

            for i in range(0,I):
            
                j = np.arange(0,J[i],dtype=int)  
                mu_0[i] = find_mu('coverage',nu_eq_star[i],J[i],S[i],m_0[i],\
                                  m_scaling[i,j],h_0[i],h_scaling[i,j],a_0[i],\
                                  a_scaling[i,j],g_scaling[i,j],alpha[i],\
                                  mu_init=0.1,\
                                  continuous_check=continuous_check)                
                N[i,j] = N_dist_eq(mu_0[i],J[i],S[i],m_scaling[i,j],\
                                   a_scaling[i,j],a_0[i],g_scaling[i,j],\
                                   nu_eq_shade[i])

                if P_or_gamma == 'gamma':
                    
                    P_eq_shade[i] = sum(gamma_base[i] * m_0[i]/mu_0[i]*\
                                        N[i,j] * g_scaling[i,j])/(1.0-alpha[i])

                
                
        elif PFT_competition == False:
            
            for i in range(0,I):
                
                j = np.arange(0,J[i],dtype=int)  
                N[i,j] = N_dist_eq(mu_0[i],J[i],S[i],m_scaling[i,j],\
                                   a_scaling[i,j],a_0[i],g_scaling[i,j],\
                                   nu_eq[i])
            
    if (P_or_gamma == 'P' or P_or_gamma == 'gamma'\
        and PFT_competition == True):
        
        return mu_0, N, nu_eq_shade, P_eq_shade
    
    elif (P_or_gamma == 'P' and PFT_competition == False):

        return mu_0, N, nu_eq, P_eq_shade
    
    elif (PFT_competition == True):
        
        return mu_0, N, nu_eq_shade
        
    elif (PFT_competition == False):
        
        return mu_0, N, nu_eq
        
    return
    
def EQ_sample(observation,observation_value,P_or_gamma,I,J,J_max,S,m_0,\
              m_scaling,m_sample,a_0,a_scaling,g_scaling,C,c_pftg,c_unique,\
              alpha,nu_min,PFT_competition,P_tot,continuous_check):
    
    """
    Function to determine the mu value from the sub sample distribution.
    
    Inputs:
    
      - observation, can be "coverage"/"nu_obs" (-), "stand density"/"N_obs"
        (/m2) "carbon mass"/"M_obs" (kgC/m2), "height"/"H_obs" (m).
    
      - observation_value, the corresponding observation value, must be an 
        array be of the same size of the maximum number of mass classes. (-)
    
       - P_or_gamma, method used to employ in splitting mu_0 into it's 
         components, either by knowing a value for P ('P') or a value for
         gamma ('gamma'). (-)
         
       - I, Number of PFTs. (-)
              
       - J, the number of mass classes. (-)
        
       - S, multiplicative bin scaling parameter. (-)
                                   
       - m_0, carbon mass at boundary. (kgC)
            
       - m_scaling, the mass scaling with respect to m_0. (-)
       
       - m_sample, the boundary masses, for each PFT contains the lower and 
         upper boundaries for the sample mass. (-)
        
       - a_0, boundary crown area. (m2)
            
       - a_scaling, the crown area scaling with respect to a_0. (-)
                   
       - g_scaling, the growth scaling with respect to g_0. (-)
    
       - C, the number of unique PFT groups. (-)
                  
       - c_pftg, the PFT group dominance hierarchy. (-)
        
       - c_unique, the unique rows of each of the PFT group hierarchy. (-)
                                          
       - alpha, the reseed fraction. (-)
    
       - nu_min, minimum vegetation fraction. (-)
        
       - PFT_competition, if True solve under the assumptions of the
         competitive exclusions of RED. Otherwise purely solve to fit the
         observations, therby ignoring the dynamical equilibrium. (-)
          
       - P, Net Primary Productivity minus the Local Litterfall. The total 
             carbon assimilate for the gridbox. (kg/m2/year)
        
       - continuous_check, if True solves for the equilibrium using the 
         continuous forms of the RED equations. (-)
         
    """
    # Firstly, check to see if the dimensions are correct
    
    if (len(observation_value[0,:]) is not J_max or\
        len(observation_value[:,0]) is not I):
        
        raise UserWarning('\nSample Error:\n Observation do not have the'+\
                          ' same size when compared to the'+\
                          ' the maximum number of classes.'+\
                          ' Please input an array of observations of size:'+\
                          str(I)+'x'+str(J_max))
    
    # Secondly, find the closet mass classes nearests to the chosen bounds.   
    # We choose this to be the nearest mass left nearest class or zero for the
    # left boundary. or the right nearest or top class for the right boundary.
    j_lower = np.zeros(I,dtype=int)
    j_upper = np.zeros(I,dtype=int)
    mu_0 = np.zeros(I)     
    nu_eq = np.zeros(I)
    N = np.zeros((I,J_max))
    
    for i in range(0,I):
        
        m_lower_diff = float('inf')
        m_upper_diff = float('inf')

        for j in range(0,J[i]):

            m = m_0[i] * m_scaling[i,j]
            
            down_diff = abs(m - m_sample[0][i]) 
            up_diff = abs(m - m_sample[1][i])
            
            if m_lower_diff > down_diff and m <= m_sample[0][i]:
                
                m_lower_diff = down_diff
                j_lower[i] = j

            elif m_sample[0][i] < m_0[i]:
                
                j_lower[i] = 0

            if m_upper_diff > up_diff and m >= m_sample[1][i]:
                
                m_upper_diff = up_diff
                j_upper[i] = j
    
            elif m_sample[1][i] > m_0[i] * m_scaling[i,J[i]-1]:
                
                j_upper[i] = J[i] - 1
                        
        if j_lower[i] > j_upper[i] and J[i] > 1:
            
            raise UserWarning('\n\nSample Error:\n\nThe lower boundary'+\
                              ' class is equal to the upper boundary'+\
                              ' class and the total number of mass '+\
                              'classes is greater than one. Please '+\
                              'pick a wider sample.')
        
        # Assuming the current observations have a close to identical
        # distribution to that of the equilibrium, we first need to find 
        # the equivalent equilibrium coverage that would lend to the 
        # observation.
        j = np.arange(0,J[i],dtype=int)
        mu_0[i] = find_mu_sample(observation,observation_value[i,j],J[i],\
                                 j_lower[i],j_upper[i],S[i],m_scaling[i,j],\
                                 g_scaling[i,j],mu_init=0.3)      
        nu_eq[i] = nu_eq_disc(mu_0[i],J[i],S[i],m_scaling[i,j],\
                 m_scaling[i,j],a_scaling[i,j],g_scaling[i,j],\
                 alpha[i],nu_min)
        
    # Next we take into account the cross PFT competitive dynamics, the 
    # equilibrium coverage is adjusted according to the model paper. 
    # Space currently occupied by sub-dominant PFTs is given over to 
    # the PFT with the highest current equilibrium coverage. While the 
    # adjusted equilibrium coverages are assumed to be somewhere 
    # between the summation within functional groups and the current
    # observation derived equilibrium coverage.
    
    if PFT_competition == True:
        
        nu_eq_star, nu_eq_shade = nu_eq_star_disc(I,nu_eq,nu_min,\
                                                 C,c_pftg,c_unique)
        
        for i in range(0,I):
            
            j = np.arange(0,J[i],dtype=int) 
            mu_0[i] = find_mu('coverage',nu_eq_star[i],J[i],S[i],m_0[i],\
                              m_scaling[i,j],m_0[i],m_scaling[i,j],a_0[i],\
                              a_scaling[i,j],g_scaling[i,j],alpha[i],\
                              mu_init=0.1)
            N[i,j] = N_dist_eq(mu_0[i],J[i],S[i],m_scaling[i,j],\
                               a_scaling[i,j],a_0[i],g_scaling[i,j],\
                               nu_eq_shade[i])        
            
                            
        if P_or_gamma == 'P':
            
            if nu_eq[i] != 0:
            
                P_eq_shade[i] = P[i]*nu_eq_shade[i]/\
                                nu_eq[i]
                
            else:
                
                P_eq_shade[i] = 0.
                
        elif PFT_competition == False:
            
            for i in range(0,I):
                
                j = np.arange(0,J[i],dtype=int)  
                N[i,j] = N_dist_eq(mu_0[i],J[i],S[i],m_scaling[i,j],\
                                   a_scaling[i,j],a_0[i],g_scaling[i,j],\
                                   nu_eq[i])
            
    if (P_or_gamma == 'P' and PFT_competition == True):
        
        return mu_0, N, nu_eq_shade, P_eq_shade
    
    elif (P_or_gamma == 'P' and PFT_competition == False):

        return mu_0, N, nu_eq, P_eq_shade
    
    elif (PFT_competition == True):
        
        return mu_0, N, nu_eq_shade
        
    elif (PFT_competition == False):
        
        return mu_0, N, nu_eq
        
    return
        
    
def EQ_P_gamma(observation_type,P,gamma_base,I,J,J_max,S,m_0,m_scaling,a_0,\
               a_scaling,g_scaling,C,c_pftg,c_unique,alpha,nu_min,\
               PFT_competition):
    
    """
    From a total gridbox assmilate rate and a baseline mortality derive the
    RED equilibrium distribution and mu_0 values.
    
    Inputs:
        
        - P, Net Primary Productivity minus the Local Litterfall. The total 
             carbon assimilate for the gridbox. (kg/m2/year)
 
        - gamma_base, baseline mortalilty. (population/year)
            
        - I, Number of PFTs. (-)
              
        - J, the number of mass classes. (-)
        
        - S, multiplicative bin scaling parameter. (-)
                                   
        - m_0, carbon mass at boundary. (kgC)
            
        - m_scaling, the mass scaling with respect to m_0. (-)
        
        - a_0, boundary crown area. (m2)
            
        - a_scaling, the crown area scaling with respect to a_0. (-)
        
        - g_scaling, the growth scaling with respect to g_0. (-)

        - C, the number of unique PFT groups. (-)
                  
        - c_pftg, the PFT group dominance hierarchy. (-)
        
        - c_unique, the unique rows of each of the PFT group hierarchy. (-)
                                          
        - alpha, the reseed fraction. (-)

        - nu_min, minimum vegetation fraction. (-)  
        
    Outputs:
    
        - N, the RED equilibrium number density (population/m2)
          
        - mu_0, the turnover ratio a boundary tree defined as: \
          gamma_base * m_0/g_0. (-)        
    
    """
    mu_0 = np.zeros(I)
    nu_eq = np.zeros(I)
    N = np.zeros((I,J_max))
    
    if observation_type == 'total':
         
        # Assuming the current observations have a close to identical
        # distribution to that of the equilibrium, we first need to find 
        # the equivalent equilibrium coverage that would lend to the 
        # observation.
        
        for i in range(0,I):
            
            j = np.arange(0,J[i],dtype=int)  
            mu_0[i] = find_mu_P(P[i],gamma_base[i],J[i],S[i],m_0[i],\
                              m_scaling[i,j],a_0[i],a_scaling[i,j],\
                              g_scaling[i,j],alpha[i],mu_init=0.2)
            nu_eq[i] = nu_eq_disc(mu_0[i],J[i],S[i],m_scaling[i,j],\
                 m_scaling[i,j],a_scaling[i,j],g_scaling[i,j],\
                 alpha[i],nu_min)
        
        # Next we take into account the cross PFT competitive dynamics, the 
        # equilibrium coverage is adjusted according to the model paper. 
        # Space currently occupied by sub-dominant PFTs is given over to 
        # the PFT with the highest current equilibrium coverage. While the 
        # adjusted equilibrium coverages are assumed to be somewhere 
        # between the summation within functional groups and the current
        # observation derived equilibrium coverage.
        
        if PFT_competition == True:
            
            nu_eq_star, nu_eq_shade = nu_eq_star_disc(I,nu_eq,nu_min,\
                                                     C,c_pftg,c_unique)
            
            
        # Refind mu_0, from the derived equilibrium cover. Then find the
        # equilbrium number density distribution.
        
        if PFT_competition == True:
            
            for i in range(0,I):
                
                j = np.arange(0,J[i],dtype=int) 
                mu_0[i] = find_mu('coverage',nu_eq_star[i],J[i],S[i],m_0[i],\
                                  m_scaling[i,j],m_0[i],m_scaling[i,j],a_0[i],\
                                  a_scaling[i,j],g_scaling[i,j],alpha[i],\
                                  mu_init=0.1)

                N[i,j] = N_dist_eq(mu_0[i],J[i],S[i],m_scaling[i,j],\
                                   a_scaling[i,j],a_0[i],g_scaling[i,j],\
                                   nu_eq_shade[i])
                
        elif PFT_competition == False:
            
            for i in range(0,I):
                
                j = np.arange(0,J[i],dtype=int)  
                N[i,j] = N_dist_eq(mu_0[i],J[i],S[i],m_scaling[i,j],\
                                   a_scaling[i,j],a_0[i],g_scaling[i,j],\
                                   nu_eq[i])
            
            
    return mu_0, N
    
    
def get_gamma(J,mu_0,N_eq,P_eq,m_0,g_scaling,alpha):
    
    """
    For a given mu_0, equilibrium coverage and total carbon assimilate, 
    find the required mortality rate to fit these values.    
    
    Inputs:
        
        - J, the number of mass classes. (-)
        
        - mu_0, the turnover ratio a boundary tree defined as: \
          gamma_base * m_0/g_0. (-)
          
        - N_eq, the RED equilibrium number density across the classes
          (population/m2)
          
        - P, Net Primary Productivity minus the Local Litterfall. The total 
             carbon assimilate for the gridbox. (kg/m2/year)
          
        - m_0, carbon mass at boundary. (kgC)
        
        - g_scaling, the growth scaling with respect to g_0. (-)
        
        - alpha, the reseed fraction. (-)
        
    Outpus:
        
        - gamma_base,  the baseline mortalility. (population/year)
    
    """
    
    # The scaling sum of N with respect to G
    N_g_scaling_sum = 0.
    
    
    for j in range(0,J):
        
        N_g_scaling_sum = N_g_scaling_sum + N_eq[j] * g_scaling[j] 

    if N_g_scaling_sum != 0:
        
        g_0 = (1-alpha)*P_eq/N_g_scaling_sum

        gamma_base = mu_0 * g_0 / m_0
        
        # If the baseline mortality rate is bellow or equal to zero set to 
        # infinity
        
        if gamma_base <= 0:
            
            gamma_base = float('inf')

    else:
        
        gamma_base = float('inf')
        
        
    return gamma_base