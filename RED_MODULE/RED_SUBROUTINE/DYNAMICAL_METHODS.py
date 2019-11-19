"""

Python 2.7, Encoding: UTF-8

Created on Mon Jan 24 2019

@author: A.Argles (aa760@exeter.ac.uk)
    
Description: Methods used in the Dynamical Subroutine Libary.

"""  
    
def CFL_cond(dt,m_inc,g,C_max=1.0):
    """
    Computes the Courant-Friedrichs-Lewy Condition for the Fokker-Planck 
    equation within the dynamical runs that update the demographic equilibrium.
    If the condition is broken, then the numerical simulation will not converge
    towards the non-numerical solutions.
    
    Inputs:
            
        - dt,timestep used for integration. (year)
            
        - m_inc, the mass class bin width. (kg)
        
        - g, the competitive growth rate. (kg/year)
              
        - C_max, the maximum allowed Courant number. (-)
              
    Outputs:
            
        - C, the Courant number. (-)
            
    """
    
    if m_inc != 0.0:
        
        C = abs(g) * dt/m_inc
        
    else:
        
        C = float('inf')
    
    return C
    
    
    