"""

Python 2.7, Encoding: UTF-8

Created on Mon Jul 3 2018

@author: A.Argles (aa760@exeter.ac.uk)
    
Description: Describes the dynamical methods used drive RED forwards through
time.

"""


import numpy as np
from scipy.optimize import root
from .DYNAMICAL_METHODS import * 



def CONVERSION(P,J,mult,m,m_init,N,nu_tot,nu_min,a_init,gamma_init,alpha,\
               phi_g):
    
    """
    Estimates of the initial growth rate from an inputted number density and 
    vegetative growth.
    
    Inputs:
                       
        - P, productivity minus litterfall. The amount of biomass grown into the
          structure. (kg/m2/year)
        
        - J, number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - m, stored carbon biomass of a single tree. (kg)
        
        - m_init, the initial mass. (-)
        
        - N, the number density. (population/m2)
    
        - nu_tot, the total vegetation fraction. (-)
        
        - nu_min, the minimum vegetation fraction allowed. (-)
        
        - a_init, the initial area. (m2)
        
        - alpha, the fraction of P going into seedling production. (-)

        - phi_g, growth-mass scaling power. (-)
                
    Outputs:
       
        - m, the mass range (kg)

        - N, the  number density (population/m2)
        
        - g_init, the estimated initial growth rates (kg/m2/year)
        
        - nu_eq, the equilibrium vegetation fraction. (-)
    
        - gamma_init, the baseline mortality. (/year)
    """
    
    ### Declare vars:
    #
    #  G_str, Estimatation of the total structral growth of a PFT 
    #   
    #  NG_cont, The contribtution of the given number density towards the 
    #  structral growth.

    G_str = (1 - alpha) * P * nu_tot

    #  Loop runs through the mass ranges. While getting the contribution of
    #  each mass class tot he total growth.
    #  
    #  P, the productivity ideally should be a mean measurement over a period
    #  as the variable is volitile.
    #
    #  Estimate the number density in mass space.
    #
    #  If by the end of the loop, Ng_tot is equal to zero. Adjust the number
    #  density within the first class to equal the minimum vegetation fraction.
    #  This is because N is all zeros.
    #
    #  Find the initial growth rate.

    
    Ng_cont = 0.0
    
    for j in range(0,J):
        
        if j == 0:
            
            m[0] = m_init 
        
        else:
            
            m[j] = m[j-1] * mult

        Ng_cont = Ng_cont + (N[j]*(m[j]/m[0])**phi_g) 
        
        
    if Ng_cont == 0:
        
        N[0] = nu_min/a_init
        
        Ng_cont = N[0]

    if G_str > 0.:
        g_init = G_str/Ng_cont
    else:
        g_init = 0.
    
    return  m, N, g_init


def ALLOMETRY(PFT_group,J,m,N,h,h_init,a,a_init,g_iso,g_init,nu,gamma,\
              gamma_init,gamma_add,phi_g,phi_h,phi_a):
    """
    Function to return estimates for the variations of a PFTs height, crown 
    area and growth. There is an assumed differnce between grasses and tree
    based PFTs:
        
            - Woody PFTs follow the theory set out by Niklas & Spatz 2004*, 
              that the relationships are determined by a power law. 
              (phi_g = 3/4,phi_h = 1/4, phi_a = 1/2)
              
            - Grasses are assumed to have a similar power scaling although
              there is a lack of information to do with size scaling. Currently
              we assume one mass class.
              
    Inputs:
        
            - PFT_group, Defines the physiology of the PFT, Tree, Shrub or
              Grass. (-)
            
            - J, The number of mass classes within a PFT
            
            - m, stored carbon biomass. (kg)
            
            - h, the height of the tree across mass. (m)
                        
            - h_init, corresponding initial height. (m)
            
            - a, the crown area of a tree across masses. (m2)
        
            - a_init, corresponding initial crown area. (m2)    
            
            - g_iso, the non-competitive growth rate of an individual. 
              (kg/year)
              
            - g_init, the initial growth rate. (kg/year)
            
            - nu, the vegetation fraction across masses (-)
    
            - gamma, the total mortality across masses (population/year)
            
            - gamma_init, the baseline mortality. (population/year)
            
            - gamma_add, any additional mortality. (population/year)
                    
            - phi_g, growth-mass scaling power. (-)
            
            - phi_h, height-mass scaling power. (-)
            
            - phi_a, crown area-mass scaling power. (-)
            
    Out:
        
            - h, the height of the tree across mass. (m)
    
            - a, the crown area of a tree across mass. (m2)   
    
            - g_iso, the non-competitive growth rate of an individual 
              across mass. (kg/year)
              
            - nu, the vegetation fraction across masses. (-)
            
            - gamma, the total mortality across masses. (population/year)
            
    
    *Niklas, Karl J., and Hanns-Christof Spatz. "Growth and hydraulic (not
    mechanical) constraints govern the scaling of tree height and mass." 
    Proceedings of the National Academy of Sciences 101.44 (2004): 15661-15663.
    """
    
    ### Tree PFTs:
    #
    #  Run through the mass classes and apply N&S allometry.
    #
    #  Add up the baseline and additional mortality to find the total.
    #
    #  Find the corresponding vegetation fractions across mass classes.
    
    if PFT_group == 'T':
        
        for j in range(0,J):
            
            h[j] = h_init * (m[j]/m[0])**phi_h
            
            a[j] = a_init * (m[j]/m[0])**phi_a

            g_iso[j] = g_init * (m[j]/m[0])**phi_g

            nu[j] = N[j] * a[j]

            gamma[j] = gamma_init + gamma_add[j]


    ### Shrub PFTs:
    
    elif PFT_group == 'S':
        
        for j in range(0,J):
            
            h[j] = h_init * (m[j]/m[0])**phi_h
            
            a[j] = a_init * (m[j]/m[0])**phi_a

            g_iso[j] = g_init * (m[j]/m[0])**phi_g

            nu[j] = N[j] * a[j]

            gamma[j] = gamma_init + gamma_add[j]

    ### Grass PFTs:
    #
    #  Currently grasses only have one mass class, so this is only for the 
    #  potenial for multiple grasses
    
    elif PFT_group == 'G':
        
        for j in range(0,J):
            
            h[j] = h_init * (m[j]/m[0])**phi_h
            
            a[j] = a_init * (m[j]/m[0])**phi_a

            g_iso[j] = g_init * (m[j]/m[0])**phi_g

            nu[j] = N[j] * a[j]
            
            gamma[j] = gamma_init + gamma_add[j]





    
    return h, a, g_iso, nu, gamma

    
    
    
def COMPETITION(PFT_group,I,J_max,m_init,N,Seeds_in,h,g,g_iso,G_seed,nu,alpha):
    """
    Simulates light competition between PFTs, returns the competitive growth
    rates assumed to be reduced according to by minimum overlap (Perfect
    plasticity assumption) for cohort members. While the growth rates for
    seedlings are assumed to follow random overlap and have no plasticity.
    
    Inputs:
        
        - I, the number of PFTs. (-)
        
        - m_init, the initial mass, which is assumed to be the seedling mass 
          (kg)
        
        - J_max, the maximum number of mass classes. (-)
        
        - N, the number density of the PFT. (population/m2)
        
        - Seeds_in, The seed flux, how many seedlings are grown within a year.
          (Population/m2/year)
        
        - h, the height of the tree across mass. (m)
        
        - g, the competitive growth rate. (kg/year)
        
        - G_seed, the competitive seedling growth rates. (kg/year/m2)
                
        - nu, the vegetation fraction of a PFT across each size class. (-)
        

        
        - alpha, the fraction of competitive productivity going into seedling
          production.
        
        
    Outputs:
        
        - g, the competitive growth rate. (kg/year)
        
        - G_seed, the competitive seedling growth rates. (kg/year/m2)
        
        - Seeds_in, The seed flux, how many seedlings are grown within a year.
          (Population/m2/year)
        
    """
    ### Flattening
    #
    #  Find the total length of the flattened height matrix, L.
    #
    #  Define a flattened height, pft and mass class matrix, h_flat, i_flat,
    #  j_flat.
    #
    #  Declare integer array to record the PFT and mass class of the element,
    #  pft_flat,m_flat. As python is row-major, i -> i + 1 until i = I - 1 we 
    #  reset with i = 0 and j -> j + 1.
    #
    #
    #  Calls a quicksort algorithm which sorts from smallest to tallest trees
    #  while recording the current PFT, i, and mass class, j, to retain matrix
    #  position cordinates.
    
    L = I * J_max
    
    h_flat = np.zeros(L)
    
    i_flat = np.zeros(L,dtype=int)
    
    j_flat = np.zeros(L,dtype=int)
    
    
    for l in range(0,L):
        
        if l == 0:
            
            i = 0
            j = 0
        
        h_flat[l] = h[i,j]

        i_flat[l] = i

        j_flat[l] = j

        i = i + 1
        
        if i > I - 1:
            
            i = 0
            
            j = j +1
            
    ### Quicksort
    #
    #  Calls a quicksort algorithm which sorts from smallest to tallest trees
    #  while recording the current PFT, i, and mass class, j, to retain matrix
    #  position cordinates.
    #
    #  We define the sorted cords:
        
    h_ordered, i_ordered, j_ordered = quicksort_record_row_col(h_flat,i_flat,\
                                                               j_flat,0,L-1)
    
    
    ### frac_par & frac_par_seed
    #
    #  For loop runs through each of the mass classes to find the fraction of 
    #  photosynthetic radiation at a paticular height(frac_par) and for 
    #  seedlings, frac_par_seed.
    #
    #
    #  This is done by adding the cummaltlative vegetation fractions from the
    #  tallest class to the smallest.
    #
    #  We assume frac_par_top = 1.0, if the nu_above > 1, then frac_par -> 0.0,
    #  otherwise frac_par -> 1.0
    #
    #  Seedlings follow random overlap, the frac_par is proportional to one, 
    #  minus the total growth.
    #
    #  The competitive growth is defined as the product of the fraction of
    #  photosyntheically active radiation and isolated growth.
    #
    #  Competitive growth rates are calcualted within the next for loop.
    #
    #  Seedling growth rates area assumed to be attenuated by the ammount
    #  across PFTs.
    #
    #  Also estimates the initial boundary condition (Seeds_in) for the
    #  demography function.

    
    frac_par = np.zeros((I,J_max))
    
    nu_above_grass = 0.0
    
    nu_above_shrub = 0.0
    
    nu_above_tree = 0.0
        
    frac_par_top = 1.0
    
    
    top = L - 1
    
    for l in range(top,-1,-1):
        
        # Making sure that the loop does not deal with any non-used classes.
        if h_ordered[l] == 0:
            
            break
        i = i_ordered[l]
        j = j_ordered[l]

        if l != top:
            
            i_up = i_ordered[l+1]
            j_up = j_ordered[l+1]
            

            if PFT_group[i_up] ==  'T':
                
                nu_above_tree = nu_above_tree + nu[i_up,j_up]

                nu_above_shrub = nu_above_shrub + nu[i_up,j_up]

                nu_above_grass = nu_above_grass + nu[i_up,j_up]

            elif PFT_group[i_up] == 'S':
                
                nu_above_shrub = nu_above_shrub + nu[i_up,j_up]
    
                nu_above_grass = nu_above_grass + nu[i_up,j_up]

            elif PFT_group[i_up] == 'G':
            
                nu_above_grass = nu_above_grass + nu[i_up,j_up]


        frac_par[i,j] = frac_par_top * np.ceil(max(0.0,1.0 - nu_above_grass))
        
        g[i,j] = g_iso[i,j] * frac_par[i,j]

        G_seed[i] = G_seed[i] + alpha[i]/(1-alpha[i])*N[i,j]*g[i,j]


    # Adding up the final height class before the seedling level.
    i_up = i_ordered[l+1]
    j_up = j_ordered[l+1]

    
    if PFT_group[i_up] ==  'T':
        
        nu_above_tree = nu_above_tree + nu[i_up,j_up]

        nu_above_shrub = nu_above_shrub + nu[i_up,j_up]

        nu_above_grass = nu_above_grass + nu[i_up,j_up]

    elif PFT_group[i_up] == 'S':
        
        nu_above_shrub = nu_above_shrub + nu[i_up,j_up]

        nu_above_grass = nu_above_grass + nu[i_up,j_up]

    elif PFT_group[i_up] == 'G':
    
        nu_above_grass = nu_above_grass + nu[i_up,j_up]



    frac_par_tree_seed = frac_par_top * max(0.0,(1.0 - nu_above_tree))
    
    frac_par_shrub_seed = frac_par_top * max(0.0,(1.0 - nu_above_shrub))
    
    frac_par_grass_seed = frac_par_top * max(0.0,(1.0 - nu_above_grass))

    
    for i in range(0,I):
        
        if PFT_group[i] ==  'T':
            
            
            G_seed[i] = frac_par_tree_seed * (G_seed[i])
                    
            Seeds_in[i] = max(0.0,G_seed[i]/m_init[i])        
            
        


        elif PFT_group[i] == 'S':

            G_seed[i] = frac_par_shrub_seed * G_seed[i]
        
            Seeds_in[i] = max(0.0,G_seed[i]/m_init[i])        

        elif PFT_group[i] == 'G':        

            G_seed[i] = frac_par_grass_seed * G_seed[i]

            Seeds_in[i] = max(0.0,G_seed[i]/m_init[i])   

          
            
    
    return g, G_seed, Seeds_in
    
    
    
def DEMOGRAPHY(dt,J,mult,m,M,N,Seeds_in,h,h_mean,a,g,G,nu,nu_tot,nu_min,gamma):
    """
    Updates the demographic profile of the gridbox, runs through each of the
    mass classes estimating the population flux through a given PFT, while
    reducing the population by the mortaility term.
    
    Returns the totals of the allometric parameters (biomass, growth, 
    vegetation fraction.), the mean height and the updated number density.
    
    Inputs:
        
        - dt,timestep used for integration. (year)
        
        - J, number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)

        - m, stored carbon biomass of a single tree. (kg)
        
        - M, the total biomass of the cohort. (kg)
        
        - N, the number density. (population/m2)   
        
        - Seeds_in, The seed flux, how many seedlings are grown within a year.
          (Population/m2/year)
        
        - h, the height of the tree across mass. (m)
        
        - h_mean, the mean vegetation height in the top 10% of the pop. (m)
        
        - a, the crown area of an individual across mass. (m2)  
        
        - g, the competitive growth rate. (kg/year)
        
        - G, the total structural growth of the cohort. (kg / year / m2)
        
        - nu, the vegetation fraction across masses. (-)
    
        - nu_tot, the total vegetation fraction. (-)
        
        - nu_min, the minimum vegetation fraction. (-)
        
        - gamma, the total mortality across masses. (population/year)
        
    Outputs:

        - M, The Biomass distribution. (kg/m2)
        
        - N, The updated number density. (population/m2)
        
        - h_mean, The mean Height of the grid box PFT from the top 10% of
          population. (m)
        
        - G, the total structural growth of the cohort. (kg / year / m2)
          
        - gamma_init, the baseline mortalility from the intialisation. 
          (population/year)
        
        
    """
    
    ### Demographic loop
    #
    #  Loop runs through each of the mass classes, determining the flux out of
    #  the population (pop_out). While the population entering the class 
    #  (pop_in) equivalent to the previous population leaving. The boundary
    #  condition of pop_in is determined by the number of seeds coming into the
    #  lowest class, while the population is not allowed to leave the highest
    #  mass class. We assume population does not move backward because of
    #  negative growth rates, instead this provides additional mortality.
    #
    #  Loop uses the RED demographic differential equation to change the number
    #  density, if the resultant number density is negative we truncate to
    #  zero. Once the updated number density is found, the totals of the 
    #  biomass, growth and vegetation fraction are found.
    # 
    #  If the resultant vegetation fraction is below the minimum adjustments
    #  are made to the resultant outputs.

    nu_tot = 0.0
    
    dNdt = np.zeros(J)
    
    pop_in = np.zeros(J)
    
    pop_in[0] = Seeds_in
    
    pop_out = np.zeros(J)

    nu = np.zeros(J)
    for j in range(0,J):        

        if j == J - 1:

            m_inc = m[j] * (mult - 1)
            
            pop_out[j] = 0
            
        else:
            
            m_inc = m[j] * (mult - 1)
   
            pop_out[j] = N[j]*g[j]/m_inc
        
            
            if pop_out[j] > 0.0:
                
                pop_in[j+1] = pop_out[j]
            
            else:
        
                pop_in[j+1] = 0.0     


        if gamma[j] != np.inf:
            
            # The primary differential equation.
            dNdt[j] = pop_in[j] - pop_out[j] - N[j] * gamma[j]   
            N[j] = max(0.0,N[j] + dNdt[j] * dt)
        
        else:
            
            N[j] = 0.0

        M = M + N[j] * m[j]

        G = G + N[j] * g[j]
        
        nu[j] = N[j]*a[j]

        nu_tot = nu_tot + N[j] * a[j]




    if nu_tot < nu_min:
        
        dN_min = (nu_min - nu_tot)/a[0]
        
        N[0] = N[0] + dN_min

        
        M = M + dN_min * m[0]

        G = G + dN_min * g[0]

        nu_tot = nu_tot + dN_min * a[0]
        
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

    return M, N, h_mean, G, nu_tot