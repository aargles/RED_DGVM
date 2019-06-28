#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 10:53:55 2019

@author: aa760
"""

def INITILISATION_CHECK(I,P,nu_tot,gamma_init,init_typ):
    
    """
    Checks to see if the initilisation method has valid inputs with the same 
    size as the number of PFTs.
    
    Inputs:
        

        - I, Number of PFTs (-)
        
        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
         
        - nu_tot, total vegetation fraction.  (-)
        
        - gamma_init, baseline mortalilty. (population/year)
  
        - init_typ, Only used if init == True. Picks the method used to 
          initialise the model. init_typ == 'nu_P', fits gamma_init to
          nu_tot and productivity. init_typ == 'nu_gamma_init', fits a growth
          rate from a vegetation fraction and gamma_init. init_typ == 
          'P_gamma_init', fits a vegetation fraction from a gamma_init and
          growth rate. (-)
          
    
    Outputs:
        
        - N/A
        
    """
    
    
    ere_check = False

    
    if init_typ == 'nu_P':

        ere_msg = 'Intilisation Method: ' + \
                    'Observational Coverage with Productivity.\n'        
                  
            
        if nu_tot is None:
            
            ere_check = True
            
            ere_msg = ere_msg+ '\n - Vegetation coverage (nu_tot) has not'+\
                      ' been entered.'

            if I != len(nu_tot):
                
                ere_check = True
                
                ere_msg = ere_msg + '\n - Vegetation coverage (nu_tot) has' +\
                          ' size mismatch with the number of PFTs (I).'

        if P is None:
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Productivity (P) has not been entered.'
            
        else:
            
            if I != len(P):
                
                ere_check = True
                
                ere_msg = ere_msg + '\n - Productivity (P) has size mismatch'+\
                          ' with the number of PFTs (I).'
    
    
    elif init_typ == 'nu_gamma_init':

        ere_msg = 'Intilisation Method: ' + \
                    'Observational Coverage with fixed mortalility.\n'        
                  
            
        if nu_tot is None:
            
            ere_check = True
            
            ere_msg = ere_msg+ '\n - Vegetation coverage (nu_tot) has not'+\
                      ' been entered.'

            if I != len(nu_tot):
                
                ere_check = True
                
                ere_msg = ere_msg + '\n - Vegetation coverage (nu_tot) has' +\
                          ' size mismatch with the number of PFTs (I).'            


        if gamma_init is None:
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Baseline mortalility (gamma_init)'+\
                     ' has not been entered.'
            
        else:
            
            if I != len(gamma_init):
                
                ere_check = True
                
                ere_msg = ere_msg + '\n - Baseline mortalility (gamma_init)'+\
                          ' has size mismatch with the number of PFTs (I).'
                          

    elif init_typ == 'P_gamma_init':
        
        ere_msg = 'Intilisation Method: ' + \
            'Productivity with fixed mortalility.\n'        
          
            
        if P is None:
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Productivity (P) has not been entered.'
            
        else:
            
            if I != len(P):
                
                ere_check = True
                
                ere_msg = ere_msg + '\n - Productivity (P) has size mismatch'+\
                          ' with the number of PFTs (I).'

        if gamma_init is None:
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Baseline mortalility (gamma_init)'+\
                     ' has not been entered.'
            
        else:
            
            if I != len(gamma_init):
                
                ere_check = True
                
                ere_msg = ere_msg + '\n - Baseline mortalility (gamma_init)'+\
                          ' has size mismatch with the number of PFTs (I).'
                          

    else:
        
        ere_check = True
        ere_msg = 'Invalid Intilisation method inputted.\n'+\
                  'Please choose: \n - init_typ = nu_P, \n - init_typ'+\
                  ' = nu_gamma_init \n - init_typ = P_gamma_init '
                          


            
    if ere_check == True:
        
        raise UserWarning('Error(s) in the intialisation:\n'+ere_msg)
                    
            
    return
    
    
    
    
def DYNAMICAL_CHECK(dt,I,J_max,P,N,nu_tot,gamma_init,gamma_add):
    
    """
    Checks to see if the dynamical inputs area valid and posses the same size
    as the number of mass classes and the number of PFTs.
    
    Inputs:
        
        - dt,timestep used for integration. (year)
     
        - I, Number of PFTs (-)
        
        - J_max, the maximum number of mass classes across all PFTs. (-)
                        
        - P, productivity minus litterfall. The amount of biomass grown into 
          the structure. (kg/m2/year)
         
        - N, the number density. (population/m2)
    
        - nu_tot, total vegetation fraction.  (-)
        
        - gamma_init, the baseline mortality. (population/year)
        
        - gamma_add, any additional mortality. (population/year)

    Outputs:
        
        - N/A
        
    """
    
    
    ere_check = False

    ere_msg = ''
    
    if dt is None:
        
        ere_check = True
        
        ere_msg = ere_msg + '\n - Timestep (dt) has not been defined.'
    
    if P is None:
        
        ere_check = True
        
        ere_msg = ere_msg + '\n - Productivity (P) has not been entered.'
        
    else:
        
        if I != len(P):
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Productivity (P) has size mismatch'+\
                      ' with the number of PFTs (I).'

    if N is None:
        
        ere_check = True
        
        ere_msg = ere_msg + '\n - Number denisty (N)'+\
                 ' has not been entered.'
                 
        if I != len(N[:,0]):
            
            ere_msg = ere_msg + '\n - Number density (N) has size mismatch'+\
                      ' with the number of PFTs (I).'
                      
        if J_max != len(N[0,:]):
            
            ere_msg = ere_msg + '\n - Number density (N) has size mismatch'+\
                      ' with the maximum mass class (J_max).'
        
    else:
        
        if I != len(gamma_init):
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Baseline mortalility (gamma_init)'+\
                      ' has size mismatch with the number of PFTs (I).'    
            
    if nu_tot is None:
            
        ere_check = True
        
        ere_msg = ere_msg+ '\n - Vegetation coverage (nu_tot) has not'+\
                  ' been entered.'

        if I != len(nu_tot):
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Vegetation coverage (nu_tot) has' +\
                      ' size mismatch with the number of PFTs (I).'
    

    if gamma_init is None:
        
        ere_check = True
        
        ere_msg = ere_msg + '\n - Baseline mortalility (gamma_init)'+\
                 ' has not been entered.'
        
    else:
        
        if I != len(gamma_init):
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Baseline mortalility (gamma_init)'+\
                      ' has size mismatch with the number of PFTs (I).'
                          
    if gamma_add != None:
        
        if I != len(N[:,0]):
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Additional Mortalility (gamma_add) has'+\
                      ' size mismatch with the number of PFTs (I).'
                      
        if J_max != len(N[0,:]):
            
            ere_check = True
            
            ere_msg = ere_msg + '\n - Additional Mortalility (gamma_add) has'+\
                     ' size mismatch with the maximum mass class (J_max).'
        
            
    if ere_check == True:
        
        raise UserWarning('Error(s) in the dynamical run input:\n'+ere_msg)
                    
            
    return