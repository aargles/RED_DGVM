"""
Bottom Level python script used to store allometric parameters for RED

Python 2.7, Encoding: UTF-8

Created on Mon Jul  2

@author: A.Argles (aa760@exeter.ac.uk)
"""
import numpy as np
import pickle as pk
import os
from .MISC_METHODS import quicksort_record_row_col

def SAVE_PFT_VALUES():
    """
    Repository for the PFTs used within RED, saves the values in PFT_STORE,
    currently overwrites any saved file.
    
    Inputs: 
        
        - N/A
    
    Outputs: 
        
       save output file containing:
        
            - PFT_name, Name of the PFT (-)
            
            - I, Number of PFTs (-)
            
            - J, number of mass classes. (-)
            
            - J_max, the maximum number of mass classes across all PFTs. (-)
            
            - S, bin width scaling parameter. (-)
       
            - m, stored carbon biomass. (kgC)
                           
            - m_0, carbon mass at boundary. (kgC)
            
            - m_scaling, the mass scaling with respect to m_0. (-)
            
            - h, the height of the tree across mass. (m)
            
            - h_0, boundary height. (m)
            
            - h_scaling, the height-mass scaling. (-)
            
            - h_ordered, the flattend height array. (m)
    
            - a, the crown area of a tree across mass. (m2)
    
            - a_0, boundary crown area. (m2)
            
            - a_scaling, the crown area-mass scaling. (-)
            
            - g_scaling, the growth-mass scaling. (-) 
            
            - gamma_base, baseline mortalilty. (population/year)
            
            - C, the number of unique PFT groups. (-)
         
            - c_pft, competition coefficient. How much shading among plant 
              functional groups. (-)
              
            - c_pftg, the PFT group dominance hierarchy. (-)
            
            - c_unique, the unique rows of each of the PFT group hierarchy. (-)
            
            - alpha, the fraction of tree productivity going into seedling 
              production. (-) 
            
            - phi_g, growth-mass scaling power. (-)
            
            - phi_h, height-mass scaling power. (-)
            
            - phi_a, crown area-mass scaling power. (-)
            
            - nu_min, minimum vegetation fraction. (-)
                        
            - i_ordered, ordered relative from smallest to largest height. (-)
              
            - j_ordered, ordered relative from smallest to largest height. (-)
              
    
    """
    
    #  Follows the JULES9 scheme (A.B.Harper 2016)
    # PFT description: 
    #  BET-Tr, Tropical Broadleaf Evergreen Trees
    #  BET-Te, Temperate Broadleaf Evergreen Trees
    #  BDT, Broadleaf Deciduous Tree
    #  NET, Needle-leaf Evergreen Tree
    #  NDT, Needle-leaf Deciduous Tree
    #  C3, Cool Season Grasses
    #  C4, Adapted, Tropical/Dry Grasses
    PFT_name = ['BET-Tr', 'BET-Te', 'BDT', 'NET', 'NDT','C3','C4','ESh','DSh']
    PFT_type = ['T','T','T','T','T','G','G','S','S'] # Array used in
                                                     # detailing the allometry

    ### Number of PFTs used in Simulation
    I = len(PFT_name)
   
    ### Boundary Constants
    #  Boundary Carbon Mass
    #  Boundary Sapling Height
    #  Boundaru Crown Area 
    m_0 = [1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.15, 0.15, 0.5] 
    a_0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25]      
    h_0 = [3.0, 3.0, 3.0, 3.0, 3.0, 0.05, 0.05, 3.0, 3.0]

    ### Demorgraphic
    #  Baseline mortality
    #  Reseed Fraction
    #  Competition Coefficent.
    #  Number of unique functional groups
    #  The PFT group dominance hierarchy 
    gamma_base = [0.06,0.05,0.05,0.075,0.065,0.2,0.2,0.15,0.15]
    alpha = [0.1,0.1,0.1,0.1,0.1,0.6,0.6,0.35,0.35] 
    c_pft = np.zeros((I,I))
    c_pft[0,:] = [1,1,1,1,1,0,0,0,0] # BET-Tr Competition coef --|
    c_pft[1,:] = [1,1,1,1,1,0,0,0,0] # BET-Te Competition coef   |
    c_pft[2,:] = [1,1,1,1,1,0,0,0,0] # BET-Te Competition coef   |--Trees
    c_pft[3,:] = [1,1,1,1,1,0,0,0,0] # NET Competition coef      |
    c_pft[4,:] = [1,1,1,1,1,0,0,0,0] # NDT Competition coef -----|
    
    c_pft[5,:] = [1,1,1,1,1,1,1,1,1] # C3 Competition coef  -----|-- Grasses
    c_pft[6,:] = [1,1,1,1,1,1,1,1,1] # C4 Competition coef  -----|
    
    c_pft[7,:] = [1,1,1,1,1,0,0,1,1] # Esh Competition coef -----|-- Shrubs
    c_pft[8,:] = [1,1,1,1,1,0,0,1,1] # Esh Competition coef -----|
    c_unique = np.vstack(np.vstack({tuple(row) for row in c_pft}))
    C = len(c_unique) # Number of unique rows
                      # Equals number of
                      # unique groups
    # Reorder c_unique to match the order of dominance (dominant -> 
    # subdominant). The dominant PFT functional group is defined as the PFT 
    # group with the most amount of zero's.
    pftg_order = np.zeros(C,dtype=int)
    
    for c in range(0,C):
        pftg_order[c] = I - len(c_unique[c,c_unique[c,:]==0])
    
    c_unique[:,:] = c_unique[np.argsort(pftg_order),:]
    c_pftg = np.zeros((C,I),dtype=bool)   
    for c in range(0,C):
                
        for i in range(0,I):
            
            if (c_unique[c,:] == c_pft[i,:]).all():
                
                c_pftg[c,i] = True
    


    ### Allometry Scaling Parameters
    #  Growth Power 
    #  Height Power
    #  Crown Area Power
    phi_g = [0.75,0.75,0.75,0.75,0.75,0.5,0.5,0.75,0.75]
    phi_h = [0.25,0.25,0.25,0.25,0.25,0.2,0.2,0.25,0.25] 
    phi_a = [0.5,0.5,0.5,0.5,0.5,0.4,0.4,0.5,0.5]        
    
    ### Mass Ranges and Binsizes
    #  Number of mass classes
    #  Bin scaling parameter
    #  Largest number of mass classes
    #  Minimum veg fraction    
    J = [10,10,10,10,10,1,1,8,8]                      
    S = [2.32,2.32,2.35,2.35,2.32,1.5,1.5,2.80,2.80]            
    J_max = max(J[0:I])
    nu_min = 1e-3
    
    ### Compute constant allometry, and allometric scaling.
    # m, carbon mass
    # h, height across m
    # a, the crown area across m
    #
    # Scaling terms save on memory, avoids multiple power operations when
    # finding the equilibrium.
    m = np.zeros((I,J_max))
    h = np.zeros((I,J_max))
    a = np.zeros((I,J_max))
    m_scaling = np.zeros((I,J_max))
    h_scaling = np.zeros((I,J_max))
    a_scaling = np.zeros((I,J_max))
    g_scaling = np.zeros((I,J_max))
    
    for i in range(0,I):
    
        m[i,0] = m_0[i]
        h[i,0] = h_0[i]
        a[i,0] = a_0[i]
        m_scaling[i,0] = 1.
        h_scaling[i,0] = 1.
        a_scaling[i,0] = 1.
        g_scaling[i,0] = 1.

        for j in range(1,J[i]):
            
            ### Functional Group Allometric Relationships 
            #  can be expanded upon             
            if PFT_type[i] == 'T': # Trees follow Niklas & Spatz 2004
                
                m[i,j] = m[i,j-1] * S[i]
                h[i,j] = h[i,0] * (m[i,j]/m[i,0])**phi_h[i]
                a[i,j] = a[i,0] * (m[i,j]/m[i,0])**phi_a[i]
                m_scaling[i,j] = m[i,j]/m[i,0]
                h_scaling[i,j] = h[i,j]/h[i,0]
                a_scaling[i,j] = a[i,j]/a[i,0]
                g_scaling[i,j] = (m[i,j]/m[i,0])**phi_g[i]

            elif PFT_type[i] == 'S': # Shrubs follow Niklas & Spatz 2004
            
                m[i,j] = m[i,j-1] * S[i]
                h[i,j] = h[i,0] * (m[i,j]/m[i,0])**phi_h[i]
                a[i,j] = a[i,0] * (m[i,j]/m[i,0])**phi_a[i]
                m_scaling[i,j] = m[i,j]/m[i,0]
                h_scaling[i,j] = h[i,j]/h[i,0]
                a_scaling[i,j] = a[i,j]/a[i,0]
                g_scaling[i,j] = (m[i,j]/m[i,0])**phi_g[i]

            elif PFT_type[i] == 'G': # Grasses are currently assumed to posses
                                     # only one mass class
                m[i,j] = m[i,j-1] * S[i]
                h[i,j] = h[i,0] * (m[i,j]/m[i,0])**phi_h[i]
                a[i,j] = a[i,0] * (m[i,j]/m[i,0])**phi_a[i]
                m_scaling[i,j] = m[i,j]/m[i,0]
                h_scaling[i,j] = h[i,j]/h[i,0]
                a_scaling[i,j] = a[i,j]/a[i,0]
                g_scaling[i,j] = (m[i,j]/m[i,0])**phi_g[i]             
    ### Find height order of PFT Classes
    #  h_flat the flattened array of height
    #  i_flat the flattened array of PFTs with respect to h_flat
    #  j_flat the flattened array of classes with respect to h_flat
    
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
            
    # Apply quicksort alogrithim to h_flat to return the ordered array
    h_ordered, i_ordered, j_ordered = quicksort_record_row_col(h_flat,i_flat,\
                                                               j_flat,0,L-1)
       
    ### Save values as PFT
    script_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),\
                                              '..', 'PFT_STORE'))
    script_name = '/PFT_values.pk1'
    PFT_dict = {'PFT_name':PFT_name,'I':I,'J':J,'J_max':J_max,'S':S,'m':m,\
                'm_0':m_0,'m_scaling':m_scaling,'h':h,'h_0':h_0,\
                'h_scaling':h_scaling,'h_ordered':h_ordered,'a':a,'a_0':a_0,\
                'a_scaling':a_scaling,'g_scaling':g_scaling,\
                'gamma_base':gamma_base,'C':C,'c_pft':c_pft,'c_pftg':c_pftg,
                'c_unique':c_unique,'alpha':alpha,'phi_g':phi_g,'phi_h':phi_h,\
                'phi_a':phi_a,'nu_min':nu_min,'i_ordered':i_ordered,\
                'j_ordered':j_ordered}
                
    with open(script_dir+script_name,'wb') as f:
        
        pk.dump(PFT_dict,f)
    
    return
    
    

        
        
def PFT_VALUES():
    """
    Repository for the PFTs used within RED, returns values for use concerning
    physical scaling and baseline death rate.
    
    Inputs: 
        
        - N/A 
    
    Outputs: 
        
        - I, Number of PFTs (-)
        
        - PFT_name, Name of the PFT (-)
        
        - PFT_group, PFT physiology woody vs grass. (-)
        
        - J, number of mass classes. (-)
        
        - mult, bin width scaling parameter. (-)
        
        - J_max, the maximum number of mass classes across all PFTs. (-)
        
        - m_init, initial mass boundary condition. (kg)
        
        - h_init, corresponding initial height. (m)
        
        - a_init, corresponding initial crown area. (m2)
        
        - gamma_init, baseline mortalilty. (population/year)
        
        - alpha, the fraction of tree productivity going into seedling 
          production. (-) 
        
        - phi_g, growth-mass scaling power. (-)
        
        - phi_h, height-mass scaling power. (-)
        
        - phi_a, crown area-mass scaling power. (-)
        
        - nu_min, minimum vegetation fraction. (-)
    
    """
    
    ### Number of PFTs used in Simulation
    #  Follows the JULES9 scheme (A.B.Harper 2016)
    
    I = 9
    
    ### PFT description 
    #  BET-Tr, Tropical Broadleaf Evergreen Trees
    #  BET-Te, Temperate Broadleaf Evergreen Trees
    #  BDT, Broadleaf Deciduous Tree
    #  NET, Needle-leaf Evergreen Tree
    #  NDT, Needle-leaf Deciduous Tree
    #  C3, Cool Season Grasses
    #  C4, Adapted, Tropical/Dry Grasses
    #  'W' = woody, 'G' = Grass - Allometry
    
    PFT_name = ['BET-Tr', 'BET-Te', 'BDT', 'NET', 'NDT','C3','C4','ESh','DSh']
    PFT_group = ['T','T','T','T','T','G','G','S','S']
    
    ### Drivers
    #  Sapling Mass
    #  Sapling Height
    #  Sapling Crown Area 

    m_init = [1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.15, 0.15, 0.5] 
    a_init = [0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25]      
    h_init = [3.0, 3.0, 3.0, 3.0, 3.0, 0.05, 0.05, 3.0, 3.0]
    
    ### Demorgraphic
    #  Gamma_init
    #  Reseed Fraction
    
    gamma_init = [0.06,0.05,0.05,0.075,0.065,0.2,0.2,0.15,0.15]
    alpha = [0.05,0.1,0.1,0.1,0.1,0.6,0.6,0.35,0.35] 
    
    ### Allometry Scaling Parameters
    #  Growth Power 
    #  Height Power
    #  Crown Area Power

    phi_g = [0.75,0.75,0.75,0.75,0.75,0.5,0.5,0.75,0.75]
    phi_h = [0.25,0.25,0.25,0.25,0.25,0.2,0.2,0.25,0.25] 
    phi_a = [0.5,0.5,0.5,0.5,0.5,0.4,0.4,0.5,0.5]        
    
    ### Mass Ranges and Binsizes
    #  Number of mass classes
    #  Bin scaling parameter
    #  Largest number of mass classes
    #  Minimum veg fraction
    
    J = [10,10,10,10,10,1,1,8,8]                      
    mult = [2.32,2.32,2.35,2.35,2.32,1.5,1.5,2.80,2.80]            
    J_max = max(J[0:I])                     
    nu_min = 0.001
    
    return     I, PFT_name, PFT_group, J, mult, J_max, m_init,h_init,\
               a_init, alpha,gamma_init, phi_g, phi_h, phi_a, nu_min