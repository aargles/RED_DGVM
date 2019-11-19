"""

Python 2.7, Encoding: UTF-8

Created on Mon Jul 3 2018

@author: A.Argles (aa760@exeter.ac.uk)
    
Description: Describes the dynamical methods used drive RED forwards through
time.

"""

from .DYNAMICAL_METHODS import * 
from .MISC_METHODS import output_key_interpet

def DYN_get_outputs(output_key,I,J,mu_0,N_m,m,a,h,g_scaling,alpha,\
                        P_tot = None,gamma_base=None,gamma=None,\
                        lambda_dem = None,f_par=None,f_par_seed=None,\
                        F_s=None,F_in=None,F_out=None):
    
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
         
        - N_m, the number of trees per unit area over mass at different 
                PFTs. (population/m2)
                
        - g_scaling, the growth-mass scaling. (-) 
         
        - alpha, the fraction of tree productivity going into seedling 
              production. (-) 
         
        - P_tot, the total gridbox carbon assimilate at the state.
          (kgC/m2/year)
              
        - gamma_base, the baseline mortalility. (population/year)
        
        - gamma, the total mortalility across mass classes. (population/year)    
    
        - lambda_dem, the total demographic litterfall, from compeitition and 
          mortality. (-)
                   
        - f_par, fraction of photosyntheticly active radiation across mass 
          classes. (-)
          
        - f_par_seed, fraction of photosyntheticly active radiation across on
          seeds. (-)
          
        - F_s, seed flux into the inital mass class. (population/m2/year)
        
        - F_in, the population flux entering into a mass class. 
          (population/m2/year)
        
        - F_out, the population flux leaving a mass class. (population/m2/year)
        
       state. (kgC/m2/year) - If P_or_gamma == 'P'
                  
    Returns:
        
        output_dir, dictonary object with the outputs and corresponding 
        strings.
    
    """
     
    # compute the N distribution
    output_dict = {}
          

    for key in output_key:
        
        if key is 'lambda_dem':
            
            output_dict.update({key:lambda_dem})
            
        else:
            output_dict.update({key:output_key_interpet(key,I,J,mu_0,N_m,m,a,\
                                                        h,g_scaling,alpha,\
                                                        P_tot=P_tot,\
                                                        gamma_base=gamma_base,\
                                                        gamma=gamma,\
                                                        f_par=f_par,\
                                                        f_par_seed=f_par_seed,\
                                                        F_s=F_s,F_in=F_in,\
                                                        F_out=F_out)})
        

    output_dict_subset = [output_dict[key] for key in output_key]
            
    return output_dict_subset