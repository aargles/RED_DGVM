"""
Top level routine definitions of the RED MODULE
Python 2.7, Encoding: UTF-8

Created on Mon Jul  2

@author: A.Argles (aa760@exeter.ac.uk)

"""

from .RED_SUBROUTINE import *
import numpy as np
import os

def RED_init(observation=None, observation_value=None,observation_type=None,\
             P_or_gamma=None,P_tot=None,gamma_base=None,overwrite_check=None,
             PFT_competition=True,continuous_check=False,output_key='default'):
    
    """
    RED_init, employs the equilibrium solutions to the discrete form of RED to 
    initialise the model, as outlined in (Argles, et al 2019.) . There are a 
    few processes, firstly if there is an observation we can diagnose the 
    required mu, and the equilibrium number density distribution. Secondly if 
    either the total carbon assimilate or baseline mortality is known we can 
    diagnose the required mortality or growth (respectfully) needed to meet the
    observation. RED_init then provides the necessary outputs required to drive
    the dynamical function (RED).
    
    Inputs:
                
    - observation, can be "stand density/N_obs", "coverage"/"nu_obs" (-), 
      "carbon mass"/"M_obs" (kgC), "height"/"H_obs" (m), 
      "basal area"/"A_obs" (m2), "assimilate mortality"/"P_gamma". (-) 
     
    - observation_value, the PFT observed value of a total or average
      quantity (todo) - apart for the observation "assimilate mortality"
      /"P_gamma", which has obervered values in the inputs of P and
      gamma_base. (-)
    
    - observation_type, the method employed by RED to find the mu through
      the discrete solutions for the total, average or statistical 
      cummalitve (i.e. meadian or top 20% ). Currently this is just the
      total quantity. Values allowed: ("total","mean" or a tuple with 
      cummlative("cumulative",a fraction) or a sample of the distribution 
      ("sample",truncation_points) where truncation_points is a list of single
      values related to the observation).
      
    - P_or_gamma, method used to employ in splitting mu into it's 
      components, either by knowing a value for P ('P') or a value for
      gamma ('gamma'). (-)
      
    - P_tot, Net Primary Productivity minus the Local Litterfall. The total PFT 
          carbon assimilate for the gridbox. (kg/m2/year)
         
    - gamma_base, the baseline mortalility from the intialisation. 
      (population/year)
      
    - overwrite_check, if True overwrites PFT constants and values stores
      in RED_MODULE/PDT_values.pk1. Otherwise Loads in the file.
      
    - PFT_competition, if True solve under the assumptions of the
      competitive exclusions of RED. Otherwise purely solve to fit the
      observations, therby ignoring the dynamical equilibrium. (-)
      
    - continuous_check, if True solves for the equilibrium using the continuous
      forms of the RED equations. (-)
      
    - output_key, either a tuple containing specific request for outputs:
                In the format ("var1","var2",...):
                
      
    Outputs:
        
    - N_eq_m, the number of trees per unit area over mass at different 
            PFTs, at equilibrium. (population/m2)
            
    - nu_eq, the equilibrium total coverage of each PFT. (-) 

    - gamma_base, the baseline mortalility. (population/year) - When 
                 P_or_gamma == 'P'
                 
    - output variables listed in the output key tuple:  
      key, the string containing the output information. 'X_tot', the total
      value of X across all classes. 'X_m', the distribution of X across 
      the mass classes. 'X_mean', the mean value of X across the 
      distribution of trees. 'X_Y_r', the value of X at the r percentile of
      Y. 'X_norm' returns the normalised distribution. Other values can be:
      'mu_0', 'alpha' and anything stored in "PFT_STORE".
    
    """    
    ### Load PFTs physiological and allometric parameters.
    #   Location: ~/RED_SUBROUTINE/STORE.py
    #
    #   If overwrite_check = True or there is no file "PFT_values.pk1" in
    #   in ~/RED_SUBROUTINE/MISC_METHODS/ call SAVE_PFT_VALUES.
    
    if overwrite_check is None:
        
        file_path = (os.path.join(os.path.dirname( __file__ ),'PFT_STORE/'))
        
        value_path = file_path+'PFT_values.pk1'
        
        if os.path.isfile(value_path) != True:
        
            overwrite_check = True
        
        else:
            
            overwrite_check = False
            
    if overwrite_check == True:
        
        SAVE_PFT_VALUES()           
        
    #   Declare into local namespace, values which determine PFTs 
    #   characteristics:  
    
    if ((P_or_gamma == 'gamma' and gamma_base is None) or \
        ((observation == 'assimilate mortality' or observation == 'P_gamma') \
         and gamma_base is None)):
        
        equilibrium_PFT_key = ('I','J','J_max','S','m_0','m_scaling','h_0',\
                               'h_scaling','a_0','a_scaling','g_scaling',\
                               'gamma_base','C','c_pftg','c_unique',\
                               'alpha','phi_g','phi_h','phi_a','nu_min')
        I,J,J_max,S,m_0,m_scaling, h_0,h_scaling,a_0,a_scaling,g_scaling,\
        gamma_base,C,c_pftg,c_unique,alpha,phi_g,phi_h,phi_a,nuzz_min = \
        LOAD_PFT_VALUES(load=equilibrium_PFT_key)

        
    else:
        
        equilibrium_PFT_key = ('I','J','J_max','S','m_0','m_scaling','h_0',\
                               'h_scaling','a_0','a_scaling','g_scaling',\
                               'C','c_pftg','c_unique','alpha',\
                               'phi_g','phi_h','phi_a','nu_min')
        I,J,J_max,S,m_0,m_scaling, h_0,h_scaling,a_0,a_scaling,g_scaling,\
        C,c_pftg,c_unique,alpha,phi_g,phi_h,phi_a,nu_min = \
        LOAD_PFT_VALUES(load=equilibrium_PFT_key)
        
    ### Check Observations and methods are valid.
    #   Location: ~RED_SUBROUTINE/EXCEPTION.py        
    INITILISATION_CHECK(I,J_max,observation,observation_value,\
                        observation_type,P_or_gamma,P_tot,gamma_base)
    ### Solving for the equilibrium number density at a given mu
    #  Location: /RED_SUBROUTINE/Equilibrium.py
    #
    #  Find the equilibrium distribution and corresponding mu value. To ignore
    #  PFT competition set PFT_competition == False
    #
    #  Note, from now on "X_eq" stands for the equilibrium value irrespective 
    #  of the PFT_competition. However this is not the case within the
    #  equilibrium subroutines.

    if (observation == 'assimilate mortality' or observation == 'P_gamma'): 

         mu_0, N_eq_m = EQ_P_gamma(observation_type,P_tot,gamma_base,I,J,\
                                   J_max,S,m_0,m_scaling,a_0,a_scaling,\
                                   g_scaling,C,c_pftg,c_unique,alpha,nu_min,\
                                     PFT_competition)
    
    elif type(observation_type) == tuple:
        
        if observation_type[0] is 'sample':

            m_sample = observation_type[1] # The bounds of the sample, must be
                                           # an array with lower and upper 
                                           # bounds for the sub sample for each
                                           # PFT
            if len(m_sample[0]) is not I:
                
                raise UserWarning('\n Error: mu_method: sample \n'+\
                                 '- The number of sample bounds must be'+\
                                 'equal to the number of PFTs, array must be'+\
                                 ' of size Ix2.')
            
            if len(m_sample) is not 2:
            
                raise UserWarning('\n Error: mu_method: sample \n'+\
                                  '- Two bounds must be provided, input arry'+\
                                  'must be of size Ix2.')
            
                
            
            if P_or_gamma == 'P':
                 
                mu_0, N_eq_m, nu_eq, P_eq  = EQ_sample(observation,\
                                                       observation_value,\
                                                       P_or_gamma,I,J,J_max,S,\
                                                       m_0,m_scaling,m_sample,\
                                                       a_0,a_scaling,\
                                                       g_scaling,C,c_pftg,\
                                                       c_unique,alpha,nu_min,\
                                                       PFT_competition,P_tot,\
                                                       continuous_check)      
            else:

                 mu_0, N_eq_m, nu_eq  = EQ_sample(observation,\
                                                  observation_value,\
                                                  P_or_gamma,I,J,J_max,S,m_0,\
                                                  m_scaling,m_sample,a_0,\
                                                  a_scaling,g_scaling,C,\
                                                  c_pftg,c_unique,alpha,\
                                                  nu_min,PFT_competition,\
                                                  P_tot,continuous_check)    
                 
                 
    else:
        
        if (P_or_gamma == 'P' or P_or_gamma == 'gamma'):

            mu_0, N_eq_m, nu_eq, P_eq = EQ_mu(observation,observation_value,\
                                              observation_type,P_or_gamma,I,\
                                              J,J_max,S,m_0,m_scaling,h_0,\
                                              h_scaling,a_0,a_scaling,\
                                              g_scaling,C,c_pftg,c_unique,\
                                              alpha,nu_min,PFT_competition,\
                                              P=P_tot,\
                                              gamma_base=gamma_base,\
                                              continuous_check=\
                                              continuous_check)   
                    
        else:
            

            mu_0, N_eq_m, nu_eq = EQ_mu(observation,observation_value,\
                                         observation_type,P_or_gamma,I,J,\
                                         J_max,S,m_0,m_scaling,h_0,h_scaling,\
                                         a_0,a_scaling,g_scaling,C,c_pftg,\
                                         c_unique,alpha,nu_min,\
                                         PFT_competition,\
                                         continuous_check=continuous_check)   
            
    if P_or_gamma == 'P':

        gamma_base = np.zeros(I)

        for i in range(0,I):
          
            j = np.linspace(0,J[i] - 1,J[i],dtype=int)  
            gamma_base[i] = get_gamma(J[i],mu_0[i],N_eq_m[i,j],P_eq[i],m_0[i],
                                        g_scaling[i,j],alpha[i])

    ### If outputs are requested through output_key, compute and return 
    #   requested values.

    if type(output_key) == tuple:

        # Check what variables are contained (total,average,over mass, 
        # physiological or percentile)
        #
        # Under /RED_SUBROUTINE/EQUILIBRIUM.py


        if (P_or_gamma == 'P' or P_or_gamma == 'gamma' or\
            observation == 'assimilate mortality'\
            or observation == 'P_gamma'):

            return EQ_get_outputs(output_key,I,J,mu_0,N_eq_m,g_scaling,alpha,\
                                  gamma_base=gamma_base,P_eq=P_eq)

        else:

            return EQ_get_outputs(output_key,I,J,mu_0,N_eq_m,g_scaling,alpha)

    elif output_key is 'default':
 
        if P_or_gamma == 'P':
            
            return N_eq_m, nu_eq, P_eq, gamma_base
        
        else:
        
            return N_eq_m, nu_eq
                
    return
    
    
def RED(dt=None,P_tot = None,N=None,nu_tot=None,gamma_base=None,gamma_add=None,\
        output_key='default'):
    
    """
    Robust Ecosystem Demography (RED) dynamic vegetation model. Simulates
    changes to the vegetation demography and total biomass at a given timestep.

    Inputs:

    - dt relative to a year. Time-dependent parameters, such as growth,
      are measured at the yearly rate. (year)
    
    - P_tot, Net Primary Productivity minus the Local Litterfall. The total 
      carbon assimilate for the gridbox. (kg/m2/year)

    - N, the number of trees per unit area of sizes and of different PFTs.
      (population/m2)
    
    - nu_tot, The total fraction of a given plant functional type. (-)
    
    - gamma_base, the baseline mortalility. (/year)
      
    - gamma_add, any additional mortality from external impacts. (/year)
    
          
    - output_key, either a tuple containing specific request for outputs:
                In the format ("var1","var2",...):
           
      
      
    Outputs:
    
    output_key='default':
        
    - N, the number of trees per unit area of sizes and of different PFTs.
      (population/m2)
    
    - nu_tot, The total fraction of a given plant functional type. (-)
    
    - lambda_dem, the total demographic litterfall, from compeitition and 
      mortality. (-)
        
    output_key=(var1,var2,..)
    
    Additional Information:
    
    - If either PROD = 0, N = 0, nu_tot = 0 or gamma_add = 0, then the
      variables will be made into a zeros array such that the script can still
      be compiled.
   

      @author: A.Argles (aa760@exeter.ac.uk)

      """
    
    
    ### Load PFTs physiological and allometric parameters.
    #   Location: ~/RED_SUBROUTINE/STORE.py
    #
    #   If overwrite_check = True or there is no file "PFT_values.pk1" in
    #   in ~/RED_SUBROUTINE/MISC_METHODS/ call SAVE_PFT_VALUES.   
    file_path = (os.path.join(os.path.dirname( __file__ ),'PFT_STORE/'))
    value_path = file_path+'PFT_values.pk1'
    
    if os.path.isfile(value_path) != True:
    
        overwrite_check = True
    
    else:
        
        overwrite_check = False
        
    if overwrite_check == True:
        
        SAVE_PFT_VALUES()           
   
    #   Declare into local namespace, values which determine PFTs 
    #   characteristic required for the operation
    
    if gamma_base is None:
        
        dynamical_PFT_key = ('I','J','J_max','m_0','m','h_0','h','h_ordered',\
                             'a_0','a','g_scaling','gamma_base','c_pft',
                             'alpha','nu_min','i_ordered','j_ordered')
        I, J, J_max, m_0, m, h_0, h, h_ordered, a_0, \
        a, g_scaling, gamma_base, \
        c_pft, alpha, nu_min, \
        i_ordered, j_ordered = LOAD_PFT_VALUES(dynamical_PFT_key)
        
    else:

        dynamical_PFT_key = ('I','J','J_max','m_0','m','h_0','h','h_ordered',\
                             'a_0','a','g_scaling','c_pft',
                             'alpha','nu_min','i_ordered','j_ordered')
        I, J, J_max, m_0, m, h_0, h, h_ordered, a_0, \
        a, g_scaling, c_pft, alpha, nu_min, \
        i_ordered, j_ordered = LOAD_PFT_VALUES(dynamical_PFT_key)
        

    ### Expectation Handling
    #  From the initial condition check to make sure that inputs are correctly
    #  inputted, otherwise inform user of input errors:
    #
    #  If no additional mortalility is added we create an array of zeros
    
    DYNAMICAL_CHECK(dt,I,J_max,P_tot,N,nu_tot,gamma_base,gamma_add)
        
    if gamma_add is None:
        
        gamma_add = np.zeros((I,J_max))
    
    ### Allometry
    #  From the allometric growth scaling relationship for G
    #  (g_scaling) and the relatationship between P_tot and G_tot,
    #  derive g0.
    #  inputs mutable
    G_tot = np.zeros(I) # Total growth
    P_s = np.zeros(I) # Tot assimilate avaliable for seedlings
    g_0 = np.zeros(I) # lower boundary growth value.
    g = np.zeros((I,J_max)) # tree growth.
    gamma = np.zeros((I,J_max)) # total mortality    
    lambda_dem = np.zeros(I)      # Demographic Litter Term
    mu_0 = np.zeros(I)
    
    for i in range(0,I):
        
        P_s[i] = alpha[i] * P_tot[i]       
        lambda_dem[i] = P_s[i] * nu_tot[i] # Add lost seedlings to demographic
                                           # litter term.
        G_tot[i] = (1-alpha[i]) * P_tot[i]
        N_g_scaling_sum = 0
        
        for j in range(0,J[i]):
            
            N_g_scaling_sum = N_g_scaling_sum + N[i,j] * g_scaling[i,j]
                    
        if N_g_scaling_sum > 0:
            
            g_0[i] = G_tot[i]/N_g_scaling_sum
            
        else:
            
            g_0[i] = 0.
        
        if g_0[i] == 0:
            
            mu_0[i] = float('inf')
            
        else:
            
            mu_0[i] = gamma_base[i] * m_0[i] / g_0[i]
        
        for j in range(0,J[i]):
            
            g[i,j] = g_0[i] * g_scaling[i,j]
            gamma[i,j] = gamma_add[i,j] + gamma_base[i]    
    
    ### Competition
    # Find shading on the seedling fraction. If there is canopy saturation, i.e
    # nu_shade >= 1, we apply the fpar_cond
    #
    # If fpar_cond == True, then we go through the size classes across all PFTs
    # that are sorted in terms of height. Once we find which height canopy
    # saturation occurs we then we set f_par = 0 for the class and lower 
    # classes
    nu_shade = np.zeros(I) # Fraction of shading at the seedling level
    f_par_cond = False     
    
    for i in range(0,I):
        
        for j in range(0,I):
            
            nu_shade[i] = min(1,nu_shade[i] + c_pft[i,j] * nu_tot[j])

            if (nu_shade[i] == 1):
                
                f_par_cond = True

                break
    
    if f_par_cond == True:

        f_par = np.ones((I,J_max))
        nu_shade_m = np.zeros((I,I,J_max))
        top = I*J_max - 1
        
        for l in range(top,-1,-1):
            
            # Making sure that the loop does not deal with any non classes.
            
            if h_ordered[l] == 0:
                
                break
            
            i = i_ordered[l]
            j = j_ordered[l]
            
            if l != top:
                
                i_up = i_ordered[l+1]
                j_up = j_ordered[l+1]

                for i2 in range(0,I):
                    
                    nu_shade_m[i2,i,j] = nu_shade_m[i2,i_up,j_up] + \
                                         c_pft[i2,i]*N[i,j]*a[i,j]
                
                    if nu_shade_m[i2,i,j] >= 1:
                        
                        f_par[i2,(h[i2,:]<=h_ordered[l])] = 0 # Apply canopy
                                                              # Saturation
                    
            
            else:
                
                for i2 in range(0,I):
                    
                    nu_shade_m[i2,i,j] = c_pft[i2,i]*N[i,j]*a[i,j]

        # If f_par is true, we need to estimate the growth lost, and add it on
        # to the demographic litter term (lamdba_dem), in addition this impacts
        # the total growth avaliable for spreading seeds, therefore P_s needs
        # to be adjusted.
        
        for i in range(0,I):

            G_tot[i] = 0.            
            
            for j in range(0,J[i]):
                
                lambda_dem[i] = (1-f_par[i,j]) * g[i,j] # Add lost growth to 
                                                        # demographic litter.
                g[i,j] = f_par[i,j] * g[i,j]
                G_tot[i] = G_tot[i] + N[i,j] * g[i,j]

            P_s[i] = alpha[i] * G_tot[i]/(1-alpha[i])
            lambda_dem[i] = lambda_dem[i] + P_s[i] * nu_tot[i]
    
    ### Demography
    # Applying an upwind method, firstly calculate the fraction of seedlings 
    # entering into the lowest mass class. The work through the mass classes
    # estimating the overall flux.
    F_s = np.zeros(I)
    F_in = np.zeros((I,J_max))
    F_out = np.zeros((I,J_max))
    dN_dt = np.zeros((I,J_max))
    N_new = np.zeros((I,J_max))
    nu_tot_new = np.zeros(I)

    for i in range(0,I):

        F_s[i] = P_s[i]/m_0[i] * (1-nu_shade[i])
    
        for j in range(0,J[i]):
            
            if j < J[i] - 1:

                F_out[i,j] = N[i,j]*g[i,j]/(m[i,j+1]-m[i,j])
                
            if j == 0:
                
                F_in[i,j] = F_s[i]
        
            else:
                
                F_in[i,j] = F_out[i,j-1]

            # Upwind scheme
            if gamma[i,j] != float('inf'):
                dN_dt[i,j] = F_in[i,j] - F_out[i,j] - gamma[i,j] * N[i,j]

            if dN_dt[i,j]*dt < -N[i,j] or gamma[i,j] == float('inf'):
                                  # Keep number density positive if rate of
                                  # change is greater than the current 
                                  # population in the current class.
            
                dN_dt[i,j] = -N[i,j]/dt
                

            lambda_dem[i] = (F_in[i,j] - dN_dt[i,j] - F_out[i,j])*m[i,j]

            N_new[i,j] = N[i,j] + dN_dt[i,j]*dt
            nu_tot_new[i] = nu_tot_new[i] + N_new[i,j] * a[i,j]
        
        if nu_tot_new[i] < nu_min:

            nu_remainder = nu_min - nu_tot_new[i]
            nu_tot_new[i] = nu_tot_new[i] + nu_remainder
            N_remainder = nu_remainder/a_0[i]
            N_new[i,0] = N_new[i,0] + N_remainder
            lambda_dem[i] = lambda_dem[i] - N_remainder * m_0[i] #Subtract 
                                                                 # nu_min cost
                                                                 # from litter
                                                                 # to conserve
                                                                 # carbon. 
                                                                

    if type(output_key) == tuple:
        
        # Check what variables are contained (total,average,over mass, 
        # physiological or percentile)
        #
        # Under /RED_SUBROUTINE/DYNAMICAL.py

        return DYN_get_outputs(output_key,I,J,mu_0,N_new,m,a,h,g_scaling,\
                               alpha,P_tot=P_tot,gamma_base=gamma_base,\
                               gamma=gamma,lambda_dem=lambda_dem,F_s=F_s,\
                               F_in=F_in,F_out=F_out)
        

    elif output_key is 'default':
       
        return N_new, nu_tot_new, lambda_dem
        
        
    return