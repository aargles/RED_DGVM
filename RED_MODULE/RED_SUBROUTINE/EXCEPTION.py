#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 10:53:55 2019

@author: aa760
"""

def INITILISATION_CHECK(I,J,observation,observation_value,observation_type,\
                        P_or_gamma,P,gamma_base):
    """
    Checks to see if the initilisation method has valid inputs with the same 
    size as the number of PFTs and valid method employed 
    
    Inputs:
       
        - I, Number of PFTs (-)
            
        - J, number of mass classes. (-)
            
        - observation, can be "stand density"/"nu_obs" (-), "carbon mass"
          /"M_obs" (kgC), "height"/"H_obs" (m), "basal area"/"A_obs" (m2) or 
         "assimilate mortality"/"P_gamma". (-) 
         
        - observation_value, the PFT observed value of a total or average
          quantity (todo) - apart for the observation assimilate mortality"
          /"P_gamma", which has obervered values in the inputs of P and
          gamma_base. (-)
        
        - observation_type, the method employed by RED to find the mu through
          the discrete solutions for the total, average or statistical 
          cummalitve (i.e. meadian or top 20% ). Currently this is just the
          total quantity. Values allowed:
              ("total","mean" or a tuple with cummlative("cumulative")).
          
        - P_or_gamma, method used to employ in splitting mu into it's 
          components, either by kntowing a value for P ('P') or a value for
          gamma ('gamma'). (-)
          
        - P, Net Primary Productivity minus the Local Litterfall. The total 
             carbon assimilate for the gridbox. (kg/m2/year)
             
        - gamma_base, the baseline mortalility from the intialisation. 
          (population/year)    
          
    Outputs:
        
        - N/A
        
    """
    ### Run model with values for 
    # If error check = True, raise error
    ere_check = False
    ere_check_primary = False  
    ere_check_secondary = False
    ere_msg = ''
    ere_msg_primary = ''
    ere_msg_secondary = ''

    if (observation is 'N_obs' or observation is 'stand density'):
        
        if observation_value is None:
            
            ere_check_primary = True
            ere_msg_primary = ere_msg_primary+ '\n - Stand Density has not'+\
                           ' been entered.'

        elif (observation_type is 'total' or observation_type is 'mean'):
            
            if len(observation_value) is not I:
                    
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary + '\n - Stand density'+\
                                  ' has size mismatch with the number of'+\
                                  ' PFTs ('+str(I)+').'
                          
            else:
                
                for obs in observation_value:
                    
                    if obs < 0:
                        
                        ere_check_primary = True
                        ere_msg_primary = ere_msg_primary + '\n - Stand'+\
                        ' Density has a negative value.'
                        
                        break                   
                    
        if ere_check_primary == True:

            ere_msg_primary = 'mu method: Observational Stand Density'+\
            ' errors:\n' + ere_msg_primary
            
            
    elif (observation is 'nu_obs' or observation is 'coverage'):
        
        if observation_value is None:
            
            ere_check_primary = True
            ere_msg_primary = ere_msg_primary+ '\n - Coverage has not been '+\
                              'entered.'

        elif (observation_type is 'total' or observation_type is 'mean'):
            
            if len(observation_value) is not I:
                    
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary + '\n - Coverage has size '+\
                                  'mismatch with the number of PFTs ('+str(I)\
                                                                     +').'
                          
            else:
                
                for obs in observation_value:
                    
                    if obs < 0:
                        
                        ere_check_primary = True
                        ere_msg_primary = ere_msg_primary + '\n - Coverage '+\
                        'has a negative value.'
                        
                        break                   
                    
        if ere_check_primary == True:

            ere_msg_primary = 'mu method: Observational Coverage errors: \n' +\
            ere_msg_primary
                    
            
    elif (observation == 'M_obs' or observation == 'carbon mass'):
            
        if observation_value is None:
            
            ere_check_primary = True
            ere_msg_primary = ere_msg_primary+ '\n - Mass has not'+\
                      ' been entered.'
                      
        elif (observation_type == 'total' or observation_type == 'mean'):
            
            if len(observation_value) is not I:
                    
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary + '\n - Mass has' +\
                          ' size mismatch with the number of PFTs (I).'
                          
            else:
                
                for obs in observation_value:
                    
                    if obs < 0:
                        
                        ere_check_primary = True
                        ere_msg_primary = ere_msg_primary + '\n - Mass has a'+\
                        ' negative value.'
                        
                        break
        if ere_check_primary == True:

            ere_msg_primary ='mu method: Observational Carbon Mass errors:'+\
                             ' \n'   +  ere_msg_primary
            
    elif (observation == 'H_obs' or observation == 'height'):

        if observation_value is None:
            
            ere_check_primary = True
            ere_msg_primary = ere_msg_primary+ '\n - Height has not'+\
                      ' been entered.'
                      
        elif (observation_type == 'total' or observation_type == 'mean'):

            if len(observation_value) is not I:
                    
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary + '\n - Height has' +\
                          ' size mismatch with the number of PFTs (I).'
                                             
            else:
                
                for obs in observation_value:
                    
                    if obs < 0:
                        
                        ere_check_primary = True
                        ere_msg_primary = ere_msg_primary + '\n - Height has'+\
                        ' a negative value.'
                        
                        break
                    
        if ere_check_primary == True:

            ere_msg_primary ='mu method: Observational Height errors:'+\
                             ' \n'   +  ere_msg_primary
                
    elif (observation == 'A_obs' or observation == 'basal area'):
        
        if observation_value is None:
            
            ere_check_primary = True
            ere_msg_primary = ere_msg_primary+ '\n - Basal Area has not'+\
                      ' been entered.'
                      
        elif (observation_type == 'total' or observation_type == 'mean'):
            
            if len(observation_value) is not I:
                    
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary + '\n - Basal Area  has' +\
                          ' size mismatch with the number of PFTs (I).'
                                                
            else:
                
                for obs in observation_value:
                    
                    if obs < 0:
                        
                        ere_check_primary = True
                        ere_msg_primary = ere_msg_primary + '\n - Basal Area'+\
                        ' has a negative value.'
                        
                        break
                    
        if ere_check_primary == True:

            ere_msg_primary ='mu method: Observational Basal Area errors:'+\
                             ' \n'   +  ere_msg_primary
                             
    elif (observation == 'P_gamma' or observation == 'assimilate mortality'):
                
        if P is None:
            
            ere_check_primary = True
            ere_msg_primary = ere_msg_primary + \
                      '\n - Total carbon assimilate/P has not been' +\
                      ' been entered.'
                      
        elif (observation_type == 'total' or observation_type == 'mean'):
                      
            if len(P) is not I:
            
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary + '\n - Total carbon '+\
                'assimilate/P has size mismatch with the number of PFTs (I).'
                                                         
            for obs in P:
                
                if obs < 0:
                    
                    ere_check_primary = True
                    ere_msg_primary = ere_msg_primary + '\n - Total carbon'+\
                    ' assimilate/P has a negative value.'
                    break
                        
        if gamma_base is None:
            
            ere_check_primary = True
            ere_msg_primary = ere_msg_primary + \
                      '\n - Baseline mortality/gamma_base has not been' +\
                      ' been entered.'
                      
        elif (observation_type == 'total' or observation_type == 'mean'):
            
            if len(gamma_base) is not I:
                
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary + '\n - Baseline '+\
                'mortality/gamma_base has size mismatch with the number'+\
                ' of PFTs (I).'
                                                             
            for obs in gamma_base:
                
                if obs < 0:
                    
                    ere_check_primary = True
                    ere_msg_primary = ere_msg_primary + '\n - Baseline'+\
                    ' mortality/gamma_base has a negative value.'
                    break
                
                
        if ere_check_primary == True:
            
              ere_msg_primary =  'mu_method: ' + \
             'Assimilate carbon and Mortality errors: \n' + ere_msg_primary
                                          
    
    elif observation is None:
        
        ere_check_primary = True        
        ere_msg_primary = ere_msg_primary + 'mu_method: Invalid \n'+\
                            '\n - Observation variable has not been chosen.'+\
                            ' To use RED_init please pick out of observation'+\
                            '= "coverage"/"nu_obs" (-), "carbon mass"/"M_obs'+\
                            '" (kgC), "height"/"H_obs" (m)'+\
                            '  or "basal area"/"A_obs" (m2).'                         
    else:
        
        ere_check_primary = True
        
        ere_msg_primary = ere_msg_primary + 'mu_method: "'+str(observation)+\
                            '" errors: \n\n - Observation picked is not'+\
                            ' modelled. To use RED_init please pick out of '+\
                            'observation = "coverage"/"nu_obs" (-), "carbon'+\
                            ' mass"/"M_obs" (kgC), "height"/"H_obs" (m)'+\
                            '  or "basal area"/"A_obs" (m2).'
                            
    if observation_type is None:
    
        ere_check_primary = True
               
        ere_msg_primary = ere_msg_primary + '\n\n' +'observation_type: '+\
                          'Invalid: \n\n - Observation type has not been '+\
                          'entered.- Input must be'+\
                         ' a valid method. observation_type = "total",'+\
                         ' "mean", or ("cummlative",fraction between 0'+\
                         ' and 1)'
        
    elif observation_type is 'total':
        
        pass
    
    elif observation_type is 'mean':
        
        pass
        
    elif type(observation_type) == tuple:
   
        if observation_type[0] is 'cummlative':
             
            if len(observation_type) is not 2:
                
                ere_check_primary = True    
                ere_msg_primary = ere_msg_primary + '\n\n'+\
                                  'observation_type: cummlative: \n\n'+\
                                  ' - Input must be a tuple of length 2.'+\
                                  'i.e ("cummlative",fraction between 0'+\
                                  ' and 1)'
                                 
        elif observation_type[0] is 'sample':
    
            if len(observation_type) is not 2:
                
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary + '\n\n'+\
                                  'observation_type: sample: \n -'+\
                                  '- Input must be a tuple of length 2.'+\
                                  'i.e. ("sample",[truncation points])'
   
            if len(observation_type[1][0]) is not I:
                
                ere_check_primary = True
                ere_msg_primary = ere_msg_primary+ '\n\n'+\
                                  'observation_type: sample: \n -'+\
                                  '- Input must contain truncation points '+\
                                  'for all PFTs.'
                
    elif observation_type is 'cummlative':
        
        ere_check_primary = True
        ere_msg_primary = ere_msg_primary + '\n\n'+\
                              'observation_type: cummlative: \n\n'+\
                              ' - Input must be a tuple of length 2.'+\
                              'i.e ("cummlative",fraction between 0'+\
                              ' and 1)'
                              
                                    
    elif observation_type is 'sample':
        
        ere_check_primary = True
        ere_msg_primary = ere_msg_primary + '\n\n'+\
                              'observation_type: sample: \n\n'+\
                              ' - Input must be a tuple of length 2.'+\
                              'i.e ("sample",[truncation points])'
            
                         
    else:
        
        ere_check_primary = True
        ere_msg_primary = ere_msg_primary+'\n\n'+'observation_type: '+\
                         str(observation_type) + ': \n - Input must be'+\
                         ' a valid method. observation_type = "total",'+\
                         ' "mean", ("cummlative",fraction between 0'+\
                         ' and 1) or ("sample",[truncation points])'
                         
                                    
    if P_or_gamma is not None: 
        
        if P_or_gamma == 'P':
            
            ere_msg_secondary = ere_msg_secondary + '\n\n'+'P_or_gamma: '+\
            'P errors: \n'
            
            if P is None:
                
                ere_check_secondary = True
                ere_msg_secondary = ere_msg_secondary + \
                          '\n - Total carbon assimilate/P has not been' +\
                          ' been entered.'
                          
            else:
                
                if len(P) is not I:
                    
                    ere_check_secondary = True
                    ere_msg_secondary = ere_msg_secondary + '\n - Total '+\
                    'carbon  assimilate/P has size mismatch with the number'+\
                    ' of PFTs (I).'
                                                                 
                for obs in P:
                    
                    if obs < 0:
                        
                        ere_check_primary = True
                        ere_msg_secondary = ere_msg_secondary + '\n - Total '+\
                        'carbon assimilate/P has a negative value.'
                        break
                        
                
            if gamma_base is not None:
                
                ere_check_secondary = True
                ere_msg_secondary = ere_msg_secondary + '\n - The baseline'+\
                ' mortality/gamma_base has been entered, despite not being'+\
                'applicable for this method.'
                          
        elif P_or_gamma == 'gamma':
            

            if gamma_base is None:
                
                ere_msg_secondary = ere_msg_secondary + '\n\n'+\
                'P_or_gamma method: gamma errors: \n'
                ere_check_secondary = True
                ere_msg_secondary = '\n - Baseline mortality/gamma_base has'+\
                ' not been been entered.'
            else:
                
                if len(gamma_base) is not I:
                    
                    ere_check_secondary = True
                    ere_msg_secondary = ere_msg_secondary + '\n - Baseline'+\
                    ' mortality/gamma_base has size mismatch with the number'+\
                    ' of PFTs (I).'
                                                                 
                for obs in gamma_base:
                    
                    if obs < 0:
                        
                        ere_check_secondary = True
                        ere_msg_secondary = ere_msg_secondary +\
                        '\n - Baseline mortality/gamma_base has a negative'+\
                        ' value.'
                        break                
                
            if P is not None:
                
                ere_check_secondary = True
                ere_msg_secondary = ere_msg_secondary + '\n - The total'+\
                ' carbon assimilate has been entered has been entered,'+\
                ' despite not being applicable for this method.'  

        else:
            
            ere_check_secondary = True
            ere_msg_secondary = ere_msg_secondary + '\n\n'+'P_or_gamma"'+\
            str(observation)+'"'+' errors: \n - is not the total carbon'+\
            ' assimilate/"P" or the baseline mortality/"gamma".'

    if ere_check_primary == True:
        
        ere_check = True
        ere_msg = ere_msg + ere_msg_primary
                      
    if ere_check_secondary == True:
        
        ere_check = True
        ere_msg = ere_msg + ere_msg_secondary
        
    if ere_check == True:
        
        raise UserWarning('Error(s) in the intialisation:\n'+ere_msg)
                    
    return
    
    
    
    
def DYNAMICAL_CHECK(dt,I,J_max,P,N,nu_tot,gamma_base,gamma_add):
    
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
        
        - gamma_base, the baseline mortality. (/year)
        
        - gamma_add, any additional mortality. (/year)

    Outputs:
        
        - N/A
        
    """
    
    
    ere_check_primary = False

    ere_msg = ''
    
    if dt is None:
        
        ere_check_primary = True
        
        ere_msg = ere_msg + '\n - Timestep (dt) has not been defined.'
    
    if P is None:
        
        ere_check_primary = True
        
        ere_msg = ere_msg + '\n - Productivity (P) has not been entered.'
        
    else:
        
        if I != len(P):
            
            ere_check_primary = True
            
            ere_msg = ere_msg + '\n - Productivity (P) has size mismatch'+\
                      ' with the number of PFTs (I).'

    if N is None:
        
        ere_check_primary = True
        
        ere_msg = ere_msg + '\n - Number denisty (N)'+\
                 ' has not been entered.'
                 
        if I != len(N[:,0]):
            
            ere_msg = ere_msg + '\n - Number density (N) has size mismatch'+\
                      ' with the number of PFTs (I).'
                      
        if J_max != len(N[0,:]):
            
            ere_msg = ere_msg + '\n - Number density (N) has size mismatch'+\
                      ' with the maximum mass class (J_max).'
        
    else:
        
        if I != len(gamma_base):
            
            ere_check_primary = True
            
            ere_msg = ere_msg + '\n - Baseline mortalility (gamma_base)'+\
                      ' has size mismatch with the number of PFTs (I).'    
            
    if nu_tot is None:
            
        ere_check_primary = True
        
        ere_msg = ere_msg+ '\n - Vegetation coverage (nu_tot) has not'+\
                  ' been entered.'

        if I != len(nu_tot):
            
            ere_check_primary = True
            
            ere_msg = ere_msg + '\n - Vegetation coverage (nu_tot) has' +\
                      ' size mismatch with the number of PFTs (I).'
    

    if gamma_base is None:
        
        ere_check_primary = True
        
        ere_msg = ere_msg + '\n - Baseline mortalility (gamma_base)'+\
                 ' has not been entered.'
        
    else:
        
        if I != len(gamma_base):
            
            ere_check_primary = True
            
            ere_msg = ere_msg + '\n - Baseline mortalility (gamma_base)'+\
                      ' has size mismatch with the number of PFTs (I).'
                          
    if gamma_add is not None:
        
        if I != len(N[:,0]):
            
            ere_check_primary = True
            
            ere_msg = ere_msg + '\n - Additional Mortalility (gamma_add) has'+\
                      ' size mismatch with the number of PFTs (I).'
                      
        if J_max != len(N[0,:]):
            
            ere_check_primary = True
            
            ere_msg = ere_msg + '\n - Additional Mortalility (gamma_add) has'+\
                     ' size mismatch with the maximum mass class (J_max).'
        
            
    if ere_check_primary == True:
        
        raise UserWarning('Error(s) in the dynamical run input:\n'+ere_msg)
                    
            
    return