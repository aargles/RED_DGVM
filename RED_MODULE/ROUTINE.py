"""
Top level routine definitions of the RED MODULE
Python 2.7, Encoding: UTF-8

Created on Mon Jul  2

@author: A.Argles (aa760@exeter.ac.uk)

Routine Tree:
    
    Routine Tree:
    
        Functions ((#) denotes call order):
        
            [RED]
              |
              |
              |
            [RED_SUBROUTINE]
                     |           
              |----------------------------|
              |                            |
            [STORE: PFT_VALUES (1)]  [SUBROUTINES: CONVERSION(2), ALLOMETRY(3),
                                      COMPETITION(4),DEMOGRAPHY(5)]
             

"""

from .RED_SUBROUTINE import *
import numpy as np



def RED(dt=None,P = None,N = None,nu_tot = None,gamma_init=None,\
        gamma_add = None,init = False, init_typ = 'nu_P'):
    
    """
    Robust Ecosystem Demography (RED) dynamic vegetation model. Simulates
    changes to the vegetation demography and total biomass at a given timestep.

    Inputs:

    - dt relative to a year. Time-dependent parameters, such as growth,
      are measured at the yearly rate. (year)
    
    - P, is the amount of primary productivity retained, and thus not lost 
      by litter, in the vegetative biomass. It is given as the net primary 
      productivity minus the litterfall. (kg/m2/year)

    - N, the number of trees per unit area of sizes and of different PFTs.
      (population/m2)
    
    - nu_tot, The total fraction of a given plant functional type. (-)
    
    - gamma_init, the baseline mortalility from the intialisation. 
      (population/year)
      
    - gamma_add, any additional mortality from external impacts 
      (population/year)
    
    - init, If Initiate is equal to true, then the initial number density is 
      assumed to be given by the steady state which the analytical solution
      produces. (-)
      
    - init_typ, Only used if init == True. Picks the method used to initialise 
      the model. init_typ == 'nu_P', fits gamma_init to nu_tot and 
      productivity. init_typ == 'nu_gamma_init', fits a growth rate from a 
      vegetation fraction and gamma_init. init_typ == 'P_gamma_init', fits a 
      vegetation fraction from a gamma_init and growth rate. (-)
      
    Outputs:
    
    - N, The updated number density. (population/m2)
    
    - h_mean, The mean Height of the grid box PFT from the top 10% of
      population. (m)
    
    - M, The Biomass distribution. (kg/m2)
    
    - gamma_init, the baseline mortalility from the intialisation. 
      (population/year)
    
    Additional Information:
    
    - If either PROD = 0, N = 0, nu_tot = 0 or gamma_add = 0, then the
      variables will be made into a zeros array such that the script can still
      be compiled.
   

      @author: A.Argles (aa760@exeter.ac.uk)

      """
    
    
    ### Saved PFT values
    #   Location: ~/RED_SUBROUTINE/STORE.py
    #  
    #   Declare into local namespace, values which determine PFTs 
    #   characteristics:
    #
    #   I, number of PFTs. PFT_name - Corresponding name, PFT_group -
    #   plant functional group (trees,shrub,grass). 
    #   
    #   J, Number of mass sizes for a given PFT. mult - Multipliticative
    #   bin scaling parameter. J_max - Maximum number of mass classes.
    #   
    #   m_init, h_init, a_init, Initial values variable values, respectivley -
    #   mass, height, crown area and for PFTs. 
    #   alpha - reseed fraction. 
    #
    #   gamma_init, background moratlility if init == True, 
    #   init_type == 'nu_gamma_init' or 'P_gamma_init'
    #
    #   From N&S allometry 2004, phi_g growth - mass scaling. phi_h height
    #   scaling & phi_a the crown-mass scaling. 
    #
    #   nu_min, minimum vegetation fraction.
    
    
    if (init == True and gamma_init is None and (init_typ == 'nu_gamma_init' \
        or init_typ == 'P_gamma_init')):

        I, PFT_name, PFT_group, J, mult, J_max, m_init, h_init, a_init, alpha,\
        gamma_init, phi_g, phi_h, phi_a,nu_min = PFT_VALUES()
        
    else:
        
        I, PFT_name, PFT_group, J, mult, J_max, m_init, h_init, a_init, alpha,\
        _, phi_g, phi_h, phi_a,nu_min = PFT_VALUES()
 
    
    ### Expectation Handling
    #  From the initial condition check to make sure that inputs are correctly
    #  inputted, otherwise inform user of input errors:
    #
    #  If no additional mortalility is added we create an array of zeros
    
    if init == True:
        
        INITILISATION_CHECK(I,P,nu_tot,gamma_init,init_typ)
        
        if N == None:
            
            N = np.zeros((I,J_max))
            
        if init_typ == 'nu_P' and gamma_init is None:
            
            gamma_init = np.zeros(I)
            
        elif init_typ == 'nu_gamma_init' and P is None:
            
            P = np.zeros(I)
            
        elif init_typ == 'P_gamma_init' and nu_tot is None:
            
            nu_tot = np.zeros(I)
        
    else:
        
        DYNAMICAL_CHECK(dt,I,J_max,P,N,nu_tot,gamma_init,gamma_add)
        
     
    if gamma_add is None:
        
        gamma_add = np.zeros((I,J_max))
        

    ### Initilisation
    #  From inputs & saved variables, estimates the equilibrium number density
    #  using the discretised steady state solutions to avoid "spinning up".
    #
    #  Location: ~ /RED_SUBROUTINE/EQUILIBRIUM.py
    #
    #  PFT_HIERARCHY - Adds all the vegetation fractions of different PFT
    #  catergories. Consining dominance and subdominance.
    #
    # If init_type == 'nu_gamma_init' or 'nu_P' Loops runs through
    # each of the PFTs:
    #
    #   
    #   Finds the maximum vegtation fractions of different PFT catergories
    #
    #   Find the indices of the PFT mass classes.
    #
    #   Find the initial growth rate, and initiate == true estimate the number
    #   distribution.
    #
    #   Estimate the physiological variation of height, crown area and growth 
    #   for an individual.
    #
    #   Estimate how the vegetation fraction varies with mass.
    #
    #   Derive how the mortality varies across the mass classes.
 
    
    # Declared Variables:
    #    Allometric Vars:
    #      m, PFT masses.
    #      g_init, PFT initial growth rate.      
    #      g_iso, Non-competitive growth rates.
    #      h, PFT height
    #      a, PFT crown area
    #           
    #   Demographic Vars:
    #      nu, the vegetation fraction across masses.       
    #      gamma, the total mortality across masses.
    #      M, the total biomass across for a PFT, across the mass range.
    #      h_mean, the mean height of a PFT within the population.
    #      G, the total structural growth of the population.
    

    m = np.zeros((I,J_max))
    
    g_init = np.zeros(I)
            
    g_iso = np.zeros((I,J_max))

    h = np.zeros((I,J_max))
    
    a = np.zeros((I,J_max))
    
    nu = np.zeros((I,J_max))
    
    gamma = np.zeros((I,J_max))
    
    M = np.zeros(I)
    
    h_mean = np.zeros(I)
    
    G = np.zeros(I)
    
    
    
    if init == True:
        
        
        nu_eq = np.zeros(I)
        
        nu_max_group = np.zeros(3)
    
        max_group_count = np.zeros(3)
    
        nu_sum_group = np.zeros(3)
    
        num_pft_group = np.zeros(3)
      
        
        if (init_typ == 'nu_gamma_init' or init_typ == 'nu_P'):

            nu_tot, nu_eq = \
            \
            PFT_HIERARCHY(P,PFT_group,I,nu_tot,nu_eq,nu_min,init_typ,\
                          num_pft_group,nu_sum_group,max_group_count,\
                          nu_max_group)
        
        
            for i in range(0,I):
                
                j = np.linspace(0,J[i] - 1,J[i],dtype=int)  
                m[i,j], N[i,j], g_init[i], gamma_init[i] = \
                \
                INTILISATION_nu(P[i],J[i],mult[i],m[i,j],m_init[i],N[i,j],\
                                nu_tot[i],nu_eq[i],nu_min,a_init[i],\
                                gamma_init[i],alpha[i],phi_g[i],phi_a[i],\
                                init_typ)
                
                
                M[i], N[i,j], h_mean[i], G[i] = \
                \
                INTILISATION_OUTPUT(PFT_group[i],J[i],m[i,j],M[i],N[i,j],\
                                    h[i,j],h_mean[i],h_init[i],a[i,j],\
                                    a_init[i],g_init[i],G[i],gamma_init[i],\
                                    nu[i,j],nu_tot[i],phi_g[i],phi_h[i],\
                                    phi_a[i])
                

        elif init_typ == 'P_gamma_init':
           
            for i in range(0,I):
                
                # Truncate Growth to zero if negative
                
                g_init[i], nu_eq[i] = \
                \
                EQUIL_FRAC(P[i],J[i],mult[i],m_init[i],nu_min,a_init[i],\
                           gamma_init[i],alpha[i],phi_g[i],phi_a[i])
        
                
            nu_tot, nu_eq = \
            \
            PFT_HIERARCHY(P,PFT_group,I,nu_tot,nu_eq,nu_min,init_typ,\
                          num_pft_group,nu_sum_group,max_group_count,\
                          nu_max_group)


            for i in range(0,I):
                
                
                j = np.linspace(0,J[i] - 1,J[i],dtype=int)  
                
                m[i,j], N[i,j], g_init[i] = \
                \
                INTILISATION_P_gamma_init(P[i],J[i],mult[i],m[i,j],m_init[i],\
                                          N[i,j],g_init[i],nu_tot[i],nu_eq[i],\
                                          nu_min,a_init[i],gamma_init[i],\
                                          alpha[i],phi_g[i])
                
                M[i], N[i,j], h_mean[i], G[i] = \
                \
                INTILISATION_OUTPUT(PFT_group[i],J[i],m[i,j],M[i],N[i,j],\
                                    h[i,j],h_mean[i],h_init[i],a[i,j],\
                                    a_init[i],g_init[i],G[i],gamma_init[i],\
                                    nu[i,j],nu_tot[i],phi_g[i],phi_h[i],\
                                    phi_a[i])

        
    
    ### Allometric
    #  From inputs & saved variables, estimates physiological variation across 
    #  masses:
    #
    #  Location: ~/RED_SUBROUTINE/DYNAMICAL.py
    #
    #
    #  CONVERSION - From the initial (or estimate) distribution and the total
    #  growth rate estimate the initial growth rate
    #
    #  ALLOMETRY - Using N&S for woody PFTs get values for how the growth,
    #  height and crown area vary with mass. Grass based PFTs are assumed to 
    #  follow a different allometry.
    
    if init == False:
        
        
        for i in range(0,I):
            
            j = np.linspace(0,J[i] - 1,J[i],dtype=int)  
    
            m[i,j], N[i,j], g_init[i] = \
            \
            CONVERSION(P[i],J[i],mult[i],m[i,j],m_init[i],N[i,j],nu_tot[i],\
                       nu_min,a_init[i],gamma_init[i],alpha[i],phi_g[i])
                
            h[i,j],a[i,j],g_iso[i,j],nu[i,j],gamma[i,j] = \
            \
            ALLOMETRY(PFT_group[i],J[i],m[i,j],N[i,j],h[i,j],h_init[i],a[i,j],\
                      a_init[i],g_iso[i,j],g_init[i],nu[i,j],gamma[i,j],\
                      gamma_init[i],gamma_add[i,j],phi_g[i],phi_h[i],phi_a[i])   

    ### Competition
    #  Location: ~/RED_SUBROUTINE/DYNAMICAL.py
    #    Competition Vars:
    #      g, Competitive growth rates
    #      G_seed, the compeititve growth going into seedling production.
    #      Seeds_in, the number of seedlings for a PFT, used as the initial
    #      boundary condition for the DEMOGRAPHY function.
    #
    #  COMPETITION - sorting all PFTs into order of heights find the
    #  respective fraction of light. Returns the adjusted growth rates. 
        

        g = np.zeros((I,J_max))
        
        G_seed = np.zeros(I)
        
        Seeds_in = np.zeros(I)
        
        
        g, G_seed, Seeds_in =\
        \
        COMPETITION(PFT_group,I,J_max,m_init,N,Seeds_in,h,g,g_iso,G_seed,nu,\
                    alpha)
        

    ### Demography
    #  Location: ~/RED_SUBROUTINE/DYNAMICAL.py
    #
    #  DEMOGRAPHY - updates the demographic profile of the grid box PFTs, while
    #  returning estimates for the total vegetative carbon, mean height and 
    #  vegetative fraction.
    

        for i in range(0,I):
            
                
            j = np.linspace(0,J[i] - 1,J[i],dtype=int)  
            
            M[i], N[i,j], h_mean[i], G[i], nu_tot[i] = \
            \
            DEMOGRAPHY(dt,J[i],mult[i],m[i,j],M[i],N[i,j],Seeds_in[i],h[i,j],\
                       h_mean[i],a[i,j],g[i,j],G[i],nu[i,j],nu_tot[i],nu_min,\
                       gamma[i,j])
            
    

    return m, M, N, h_mean, nu_tot, gamma_init