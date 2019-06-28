"""
Bottom Level python script used to store allometric parameters for RED

Python 2.7, Encoding: UTF-8

Created on Mon Jul  2

@author: A.Argles (aa760@exeter.ac.uk)
"""

#%% Primary Store

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
    a_init = [0.5, 0.5, 0.75, 0.75, 1.0, 0.5, 0.25, 0.25, 0.25]      
    h_init = [3.0, 3.0, 3.0, 3.0, 3.0, 0.05, 0.05, 3.0, 3.0]
    
    ### Demorgraphic
    #  Gamma_init
    #  Reseed Fraction
    
    gamma_init = [0.06,0.05,0.05,0.075,0.065,0.2,0.2,0.15,0.15]
    alpha = [0.1,0.1,0.1,0.3,0.3,0.5,0.45,0.35,0.35] 
    
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
    mult = [4.0,4.0,4.0,4.0,4.0,1.5,1.5,4.0,4.0]            
    J_max = max(J[0:I])                     
    nu_min = 1.0e-3
    
    return     I, PFT_name, PFT_group, J, mult, J_max, m_init,h_init,\
               a_init, alpha,gamma_init, phi_g, phi_h, phi_a, nu_min
               
#%% Secondary Store
    