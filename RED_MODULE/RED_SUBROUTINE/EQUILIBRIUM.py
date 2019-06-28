"""

Python 2.7, Encoding: UTF-8

Created on Mon Jan 24 2019

@author: A.Argles (aa760@exeter.ac.uk)
    
Description: Contains the intilisation subroutines called in the main script.

TODO: 
Also contains a module to provide snapshots of the equilibrium states.

"""


from .EQUILIBRIUM_METHODS import * 


def PFT_HIERARCHY(P,PFT_group,I,nu_tot,nu_eq,nu_min,init_typ,num_pft_group,\
                  nu_sum_group,max_group_count,nu_max_group):

    """
    Adding the vegetation fraction of each PFT group, estimates the dominant 
    and sub-dominant PFTs through shading or the equilibrium fraction. 
    
    Inputs:
        
        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
        
        - PFT_group, PFT physiology woody vs grass. (-)
        
        - I, Number of PFTs (-)
    
        - J, number of mass classes. (-)
        
        - nu_tot, total vegetation fraction.  (-)
        
        - nu_eq, the equilibrium fraction. (-)
                
        - nu_min, the minimum vegetation fraction allowed. (-)
        
        - init_typ, Only used if init == True. Picks the method used to 
          initialise the model. init_typ == 'nu_P', fits gamma_init to
          nu_tot and productivity. init_typ == 'nu_gamma_init', fits a growth
          rate from a vegetation fraction and gamma_init. init_typ == 
          'P_gamma_init', fits a vegetation fraction from a gamma_init and
          growth rate. (-)
          
        - num_pft_group, contains the number of PFTs in each PFT group. (-)
          
        - nu_sum_group, the total vegetation fraction of PFTs in a given group.
          (-)
        
        - max_group_count, the total number of grouped PFTs with the maximum 
          vegetationg fraction of the group. (-)
        
        - nu_max_group, the maximum vegetation fraction of a group. (-)
        
        
    """

    # If init_type == 'nu_gamma_init' or 'nu_P' Loops runs through
    # each of the PFTs:
    #
    #   Adds all the vegetation fractions of different PFT catergories to
    #   estimate the corresponding equilibrium fractions.

    nu_max_group[0] = nu_min
    
    nu_max_group[1] = nu_min
    
    nu_max_group[2] = nu_min
           

    if (init_typ == 'nu_gamma_init' or init_typ == 'nu_P'):
        
        for i in range(0,I):
            
            
            
            if nu_tot[i] > 1.0:
                raise UserWarning('Vegetation Fraction Modified\n'+\
                      'Vegetation fraction is above the maxima, therefore '+\
                      'the vegetation fraction is modified to be '+\
                      'less than one. Such that the analytical '+\
                      'solutions are solvable')
            
                
            elif nu_tot[i] < nu_min:
                
                nu_tot[i] = nu_min

            if init_typ == 'nu_P' and P[i] <= 0.0:
                
                nu_tot[i] = nu_min


            if PFT_group[i] == 'T':
    
                gp = 0
                
            elif PFT_group[i] == 'S':
                
                gp = 1
                
            elif PFT_group[i] == 'G':
                
                gp = 2
                
            if nu_tot[i] > nu_max_group[gp]:
                
                nu_max_group[gp] = nu_tot[i]
    
                max_group_count[gp] = 1.0
                
            elif nu_tot[i] == nu_max_group[gp]:
                
                max_group_count[gp] =  max_group_count[gp] + 1.0
    
            nu_sum_group[gp] = nu_sum_group[gp] + nu_tot[i]

            if nu_sum_group[gp] >= 1.0:
                
                nu_sum_group[gp] = 1.0 - (I-1)*nu_min
    
            num_pft_group[gp] = num_pft_group[gp] + 1.0


    
        for i in range(0,I):
            

            
            if PFT_group[i] == 'T':
    
                gp = 0
                
            elif PFT_group[i] == 'S':
                
                gp = 1
                
            elif PFT_group[i] == 'G':
                
                gp = 2
            
            
            nu_eq[i] = sum(nu_sum_group[0:gp+1]) - (nu_max_group[gp] - \
                       nu_tot[i])
        

        
            if nu_tot[i] == nu_max_group[gp]:
                

                    
                nu_tot[i] = nu_eq[i] - sum(nu_sum_group[0:gp]) - (num_pft_group[gp] - 1.0 ) * nu_min
                

                
                
        
                nu_tot[i] = max(nu_min,nu_tot[i] / max_group_count[gp])
      
            
            else:
                    
                nu_tot[i] = nu_min

    
            if nu_tot[i] >= 1.0:
                
                nu_tot[i] = 1.0 - (I-1)*nu_min

            if nu_eq[i] >= 1.0:
                
                nu_eq[i] = 1.0 - (I-1)*nu_min


    

    elif init_typ == 'P_gamma_init': 
        
 
        for i in range(0,I):
            
            if PFT_group[i] == 'T':
    
                gp = 0
                
            elif PFT_group[i] == 'S':
                
                gp = 1
                
            elif PFT_group[i] == 'G':
                
                gp = 2
                
            if nu_eq[i] > nu_max_group[gp]:
                
                nu_max_group[gp] = nu_eq[i]

                max_group_count[gp] = 1.0

            elif nu_eq[i] == nu_max_group[gp]:
                
                max_group_count[gp] = max_group_count[gp] + 1.0

            
            num_pft_group[gp] = num_pft_group[gp] + 1.0

        nu_shade = [0.,0.,0.]

        nu_shade[0] = (num_pft_group[0] - 1.0)* nu_min

        nu_shade[1] =  nu_shade[0] + (num_pft_group[1] - 1.0)* nu_min + max(nu_min,nu_max_group[0]-nu_shade[0])
        
        nu_shade[2] =  nu_shade[1] + (num_pft_group[2] - 1.0)* nu_min + max(nu_min,nu_max_group[1]-nu_shade[1])

        for i in range(0,I):
            
            if PFT_group[i] == 'T':
    
                gp = 0
                
            elif PFT_group[i] == 'S':
                
                gp = 1
                
            elif PFT_group[i] == 'G':
                
                gp = 2
                
                

            if nu_eq[i] == nu_max_group[gp]:
                
                nu_tot[i] = max(nu_min,(nu_max_group[gp] - nu_shade[gp])\
                            /max_group_count[gp])
                
            else:
                
                nu_eq[i] = nu_max_group[gp]
                
                nu_tot[i] = nu_min


  
            if nu_tot[i] >= 1.0:
                
                nu_tot[i] = 1.0 - (I-1)*nu_min

            if nu_eq[i] >= 1.0:
                
                nu_eq[i] = 1.0 - (I-1)*nu_min

        
    return nu_tot, nu_eq
            
            

def INTILISATION_nu(P,J,mult,m,m_init,N,nu_tot,nu_eq,nu_min,a_init,gamma_init,\
                    alpha,phi_g,phi_a,init_typ):
    
    """
    Uses two methods of getting the steady state involving a observational 
    coverage. Returns the demographic profile of the equilibrium.
    
    Inputs:
                
        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
        
        - J, number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - m, stored carbon biomass of a single tree. (kg)
        
        - m_init, the initial mass. (-)
        
        - N, the number density. (population/m2)
        
        - nu_tot, total vegetation fraction.  (-)
        
        - nu_eq, the equilibrium vegetation fraction allowed. (-)

        - nu_min, the minimum vegetation fraction allowed. (-)
        
        - a_init, the initial area. (m2)
        
        - gamma_init, the baseline mortality. (/year)
        
        - alpha, the fraction of P going into seedling production. (-)
        
        - phi_g, growth-mass scaling power. (-)
        
        - phi_a, crown area - mass scaling power. (-)
          
        - init_typ, Only used if init == True. Picks the method used to 
          initialise the model. init_typ == 'nu_P', fits gamma_init to
          nu_tot and productivity. init_typ == 'nu_gamma_init', fits a growth
          rate from a vegetation fraction and gamma_init. init_typ == 
          'P_gamma_init', fits a vegetation fraction from a gamma_init and
          growth rate. (-)
      
        
    Outputs:
       
        - m, the mass range (kg)

        - N, the  number density (population/m2)
        
        - g_init, the estimated initial growth rates (kg/m2/year)
    
        - gamma_init, the baseline mortality. (/year)
    """
        
        
    ### Initiate == True and init_typ == 'nu_P' or 'nu_gamma_init'.
    #
    #    If nu_tot equal or greater than 1.0, change to 1.0 - I*nu_min, such 
    #    that it compiles. This will also raise a user warning.
    #    
    #    From the discrete solutions for the vegetation fraction. Numerically
    #    solve for "mu" given the vegetation fraction.
    #
    #    From the initial growth rate we can estimate the baseline mortalilty,
    #    using the steady state solutions. This is done by estimating an growth 
    #    rate from the total structral growth. Then using this growth rate and
    #    the solved mu to estimate the baseline mortalility.
    #
    #    If input P is not provided, but a gamma initial is, we estimate the 
    #    steady state productivity.
    #
    #   * The number of the population dying is equal to the number of seedlings
    #     and seedlings can only grow in free space (see write up for 
    #     derivation of discrete equilbirium solutions). 
        
        
    if J != 1:
        
        # Estimate mu for the newton root finding method
        
        mu_init = 0.5     # Guessed mu

        mu = find_mu_nu(J,mult,nu_eq,mu_init,alpha,phi_g)


        # Finding the baseline mortalilty or total growth
        
        if init_typ == 'nu_gamma_init':

            g_init_eq = (gamma_init * m_init) / mu
            
            G_str_eq = G_str_discrete(g_init_eq,J,mult,a_init,nu_eq,mu,\
                                    phi_g,phi_a)
            
            if nu_eq != 0.0:
                G_seed = (alpha)/(1.0-alpha)* G_str_eq * nu_tot/nu_eq
            else:
                G_seed = 0.0
                
            G_str = G_seed/alpha * (1.0-alpha)

        elif init_typ == 'nu_P':
            
            if P > 0.0:
                
                G_str = (1.0 - alpha) * P * nu_tot

                G_seed = alpha * P * nu_tot
            
                G_str_eq = (1.0 - alpha) * P * nu_eq
            
                g_init_eq = g_init_discrete(G_str_eq,J,mult,a_init,nu_eq,\
                                            mu,phi_g,phi_a)
            
                gamma_init = mu * g_init_eq / m_init
                
                
            else:
                                
                gamma_init = float('inf')
                
                G_seed = alpha * P * nu_tot
                
                g_init_eq = 0.0
                
                G_str = 0.0

        else:
            
            g_init_eq = 0.0

            gamma_init = float('inf')
            
            G_str = 0.0

            
        
        m, N, Ng_cont = number_dist_resolve(J,mult,m,m_init,N,G_seed,\
                                            g_init_eq,a_init,gamma_init,\
                                            nu_tot,nu_eq,nu_min,phi_g)

    # One mass class initialisation
    
    elif J == 1:
        
        m[0] = m_init
        
        N[0] = nu_tot / a_init

        Ng_cont = N[0] 

        if init_typ == 'nu_gamma_init':
            
            mu = (1.0 - nu_eq) * alpha/(1.0-alpha)
            
            g_init_eq = gamma_init * m_init / mu

            G_str = g_init_eq * N[0]
            
        elif init_typ == 'nu_P':

            if P > 0.0:
                
                G_seed_eq = P * nu_eq * alpha
                
                N_eq = nu_eq / a_init
             
                if N_eq == 0:
                    
                    gamma_init = float('inf')
                    
                else:
                    
                    
                    gamma_init = G_seed_eq *(1-nu_eq) / (N_eq * m_init)  
                
                
                G_str = P * nu_tot * (1.0-alpha)

            
            else:
                
                gamma_init = float('inf')
                
                mu = float('inf')
                
                G_str = 0.0
    
                
        else:
            
            gamma_init = float('inf')
            
            G_str = 0.0
            
        
    if Ng_cont == 0.0:

        g_init = 0.0

    else:        
        
        g_init = G_str/Ng_cont
    
    if J ==1:
        
        mu = gamma_init * m_init / g_init
            

        
    return  m, N, g_init, gamma_init
    
    
    
def EQUIL_FRAC(P,J,mult,m_init,nu_min,a_init,gamma_init,alpha,phi_g,phi_a):
    
    """
    Using the discrete form of the steady state solutions, uses one method of 
    Returns the equilibrium fraction from a productivity and baseline
    mortalilty.
    
    Inputs:
                
        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
        
        - J, number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - m_init, the initial mass. (-)
        
        - nu_min, the minimum vegetation fraction allowed. (-)
        
        - a_init, the initial area. (m2)
        
        - gamma_init, the baseline mortality. (/year)
        
        - alpha, the fraction of P going into seedling production. (-)
        
        - phi_g, growth-mass scaling power. (-)
        
        - phi_a, crown area - mass scaling power. (-)
               
        
    Outputs:
       
        - m, the mass range (kg)

        - N, the  number density (population/m2)
        
        - g_init, the estimated initial growth rates (kg/m2/year)
        
        - nu_eq, the equilibrium vegetation fraction. (-)
    
        - gamma_init, the baseline mortality. (/year)
    """
        


    if J != 1:
    
    # Estimate mu for the newton root finding method
    
        mu_init = 0.2     # Guessed mu
        
        
        if P > 0.0:
    
            mu = find_mu_g(P,J,mult,m_init, a_init,gamma_init,mu_init,\
                           alpha,phi_g,phi_a,tol = None)
        
            g_init = gamma_init * m_init / mu
            
            nu_eq = nu_eq_discrete(J,mult,mu,nu_min,alpha,phi_g)

        else:
            
            g_init = 0.0
            
            nu_eq = nu_min

    elif J==1:


        if P > 0.0:
            
            g_init = (1.0 - alpha) * P * a_init
    
            mu = gamma_init * m_init / g_init
            
            nu_eq = max(nu_min,1.0 - mu * (1.0 - alpha)/alpha)
                        
         
            
        else:
            
            g_init = 0     
            
            nu_eq = nu_min
    
    
    return g_init, nu_eq
    

def INTILISATION_P_gamma_init(P,J,mult,m,m_init,N,g_init,nu_tot,nu_eq,nu_min,\
                              a_init,gamma_init,alpha,phi_g):
    """
    Returns the number desnity of the steady state based upon a current 
    vegetation fraction and the corresponding equilibrium fraction.
    
        Inputs:
                
        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
        
        - J, number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - m, stored carbon biomass of a single tree. (kg)
        
        - m_init, the initial mass. (-)
        
        - N, the number density. (population/m2)
        
        - g_init, the estimated initial growth rate. (kg/m2/year)
        
        - nu_tot, total vegetation fraction.  (-)
        
        - nu_eq, the equilibrium vegetation fraction allowed. (-)
        
        - nu_min, the minimum vegetation fraction allowed. (-)
        
        - a_init, the initial area. (m2)
        
        - gamma_init, the baseline mortality. (/year)
        
        - alpha, the fraction of P going into seedling production. (-)
        
        - phi_g, growth-mass scaling power. (-)
        
        - phi_a, crown area - mass scaling power. (-)
          
        - init_typ, Only used if init == True. Picks the method used to 
          initialise the model. init_typ == 'nu_P', fits gamma_init to
          nu_tot and productivity. init_typ == 'nu_gamma_init', fits a growth
          rate from a vegetation fraction and gamma_init. init_typ == 
          'P_gamma_init', fits a vegetation fraction from a gamma_init and
          growth rate. (-)
      
        
    Outputs:
       
        - m, the mass range (kg)

        - N, the  number density (population/m2)
    
    """
    
    if J != 1:
        
        G_seed = alpha * nu_tot * P
    
        m, N, NG_cont = number_dist_resolve(J,mult,m,m_init,N,G_seed,g_init,\
                                            a_init,gamma_init,nu_tot,nu_eq,\
                                            nu_min,phi_g)
    
    elif J == 1:
        
        m[0] = m_init
        
        N[0] = nu_tot / a_init

        g_init = alpha * nu_tot * P / N[0]
    
    
    return m, N, g_init
    
    
    
def INTILISATION_OUTPUT(PFT_group,J,m,M,N,h,h_mean,h_init,a,a_init,g_init,G,\
                        gamma_init,nu,nu_tot,phi_g,phi_h,phi_a):
    """
    Function returns the numerical integral of various outputs across the 
    number density. Assumes Niklas & Spatz* allometry for multiple mass classed
    pfts:
    
              
    Inputs:
        
            - PFT_group, Defines the physiology of the PFT, Tree, Shrub or
              Grass. (-)
            
            - J, The number of mass classes within a PFT
            
            - m, stored carbon biomass. (kg)
            
            - h, the height of the tree across mass. (m)
                        
            - h_init, corresponding initial height. (m)

            - h_mean, the mean vegetation height in the top 10% of the pop. (m)

            - a, the crown area of an individual across mass. (m2)  
                
            - a_init, corresponding initial crown area. (m2)    

            - g_init, the initial growth rate. (kg/year)

            - G, the total structural growth of the cohort. (kg / year / m2)
            
            - gamma_init, the baseline mortality. (population/year)
        
            - nu, the vegetation fraction across masses. (-)        
    
            - nu_tot, the total vegetation fraction. (-)
        
            - phi_g, growth-mass scaling power. (-)
            
            - phi_h, height-mass scaling power. (-)
            
            - phi_a, crown area-mass scaling power. (-)
            
    Outputs:
    
            - M, The Biomass distribution. (kg/m2)
            
            - N, The updated number density. (population/m2)
            
            - h_mean, The mean Height of the grid box PFT from the top 10% of
              population. (m)
            
            - G, the total structural growth of the cohort. (kg / year / m2)
                                    
    
    *Niklas, Karl J., and Hanns-Christof Spatz. "Growth and hydraulic (not
    mechanical) constraints govern the scaling of tree height and mass." 
    Proceedings of the National Academy of Sciences 101.44 (2004): 15661-15663.
    """
       
    ### Demographic loop
    #
    # Loop uses allometric scaling for trees, grasses and shrubs and outputs
    # the sum product of total biomass, growth and fraction.    

    for j in range(0,J):
        
        if PFT_group == 'T':
            
            h[j] = h_init * (m[j]/m[0]) ** phi_h

            a[j] = a_init * (m[j]/m[0]) ** phi_a

            g = g_init * (m[j]/m[0]) ** phi_g

        elif PFT_group == 'S':
            
            h[j] = h_init * (m[j]/m[0]) ** phi_h

            a[j] = a_init * (m[j]/m[0]) ** phi_a

            g = g_init * (m[j]/m[0]) ** phi_g

        elif PFT_group == 'G':

            h[j] = h_init
            
            a[j] = a_init
            
            g = g_init

        
        M = M + N[j] * m[j]
        
        nu[j] = N[j] * a[j]

        G = G + N[j] * g

            
    ### Height Loop
    #
    #  Loop runs from top to bottom adding the total height density and number
    #  density for the process of finding the average height within the top 10%
    #  of the population.      
    #
    #  nu must be normalilsed to isolate the fractional distribution from the 
    #  rest of the grid box.

    
    if J > 1:
        
        N_tot = 0.0 
        
        h_tot = 0.0
        
        nu_h_tot = 0.0
        
        nu_h_check = 0.1
        
        
        
        for j in range(J-1,-1,-1):
        
        
            nu_h_tot = nu_h_tot + nu[j]/nu_tot
            
            if nu_h_tot < nu_h_check:
            
                N_tot = N_tot + N[j]
                
                h_tot = h_tot + N[j] * h[j]
                
                
            elif nu_h_tot >= nu_h_check:
                            
                N_remainder = (nu_h_tot - nu_h_check)/ a[j]
            
                N_tot = N_tot + N_remainder
                
                h_tot = h_tot + N_remainder * h[j]
        
                break
            
            if N_tot > 0.0:
                
                h_mean = h_tot / N_tot
        
            else:
                
                h_mean = h[0]
                
        else:
            
            h_mean = h[0]

    return M, N, h_mean, G