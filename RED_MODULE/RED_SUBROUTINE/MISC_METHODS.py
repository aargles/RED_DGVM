"""
Python 2.7, Encoding: UTF-8

Created on Mon Jul  2

@author: A.Argles (aa760@exeter.ac.uk)
"""


#%% Script hold miscellaneous methods. 

#from EQUILIBRIUM_METHODS import nu_eq
from .EQUILIBRIUM_ANALYTICAL import nu_eq_analy
#from EQUILIBRIUM_METHODS import X, k_j
from collections import Sequence
import numpy as np
import pickle as pk
import os


def output_key_interpet(key,I,J,mu_0,N_m,m,a,h,g_scaling,alpha,\
                        P_tot = None,gamma_base=None,gamma=None,f_par=None,\
                        f_par_seed=None,F_s=None,F_in=None,F_out=None):
    """
    From a output string key, derive the output required. 
    
    Inputs: 
        
        - key, the string containing the output information. 'X_tot', the total
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
                   
        - f_par, fraction of photosyntheticly active radiation across mass 
          classes. (-)
          
        - f_par_seed, fraction of photosyntheticly active radiation across on
          seeds. (-)
          
        - F_s, seed flux into the inital mass class. (population/m2/year)
        
        - F_in, the population flux entering into a mass class. 
          (population/m2/year)
        
        - F_out, the population flux leaving a mass class. (population/m2/year)
        
        
         
    Outputs:
        
        - Output variable value requested. For a variable X and corresponding
          method.
    
    """
    nonpercentile_output_list = ['tot','m','norm','mean']
    output_check = False
    output_type = ''
    spl = key.split('_')
    
    for out in nonpercentile_output_list: # Find the type of output
        
        if spl[-1] == out:
            
            output_check = True
            output_type = out

            break

    if (key[-1:].isdigit() == True and key != 'mu_0'): # percentile could be requested 
        
        percentile = float(spl[-1])
        
        if ( percentile <= 0. and percentile <= 100.):
            
            raise UserWarning('\n Output Key Error: \n  percentile value'+\
            'wihtin chosen method must be between or equal to 0 or 100,'+\
            ' such that: "X_Y_0"<= "X_Y_r" <= "X_Y_100".')

        output_check = True
        output_type = 'r'
        
    if key == 'mu_0':
        
        output_check = True
        output_type = 'other'        
        
    elif key == 'gamma_base':
        
        if gamma_base is None:
            
            raise UserWarning('\n Output Key Error: Invalid output requested'+\
                              ' \n - The key: "gamma_base", is not'+\
                              ' determined.')
        
        else:
            
            output_check = True
            output_type = 'other'
            
    elif key == 'f_par_seed':
        
        if gamma_base is None:
            
            raise UserWarning('\n Output Key Error: Invalid output requested'+\
                              ' \n - The key: "f_par_seed", is not'+\
                              ' determined.')
            
        else:
            
            output_check = True
            output_type = 'other'
            
    elif key == 'F_s':
        
        if F_s is None:
            
            raise UserWarning('\n Output Key Error: Invalid output requested'+\
                              ' \n - The key: "F_s", is not'+\
                              ' determined.')
            
        else:
            
            output_check = True
            output_type = 'other'   
        
        
    if output_check == False: # or a PFT specific value or mu or gamma_base
        
        pft_store_list = ['PFT_name','I','J','J_max','S','m_0','m_scaling',\
        'h_0','h_scaling','h_ordered','a_0','a_scaling','g_scaling','C',\
        'c_pft','c_pftg','c_unique','alpha','phi_g','phi_h','phi_a',\
        'nu_min','i_ordered','j_ordered']
        
        for out in pft_store_list:
            
            if key == out:
            
                output_check = True
                output_type = 'store'
                
                break
            
    if output_check == False: 
                
        raise UserWarning('\n Output Key Error: Invalid output requested: \n'+\
                          ' - The key: '+key+' is invalid. Please put in:'+\
                          ' "X_tot", "X_m", "X_mean", "X_Y_r", "X_norm"'+\
                          ' or a Variable defined in PFT_STORE')
        

    if output_type == 'store':
        
        return LOAD_PFT_VALUES(tuple([key]))[0]
        
    elif output_type == 'other':
        
        if key == 'gamma_base':
                    
            return gamma_base
            
        elif key == 'mu_0':
            
            return mu_0

        elif key == 'f_par_seed':
            
            return f_par_seed
            
        elif key == 'F_s':
            
            return F_s
        
    else:
        
        variable_list = ['N','M','P','G','P_s','H','nu','h','p','g','p_s',\
                         'a','gamma','lambda_dem','f_par','F_in','F_out','m']

        if output_type == 'r':
            
            sep = key.split('_')
            var1 = sep[0]
            var2 = '-'.join(sep[1:])
            first_variable_check = False
            second_variable_check = False
            
        else:
            
            variable_check = False
            
        for var in variable_list:
            
            if output_type == 'r':

                if var1.find(var) != -1:
                    
                    if ((var=='P' or var == 'G' or var == 'P_s' or \
                           var == 'p' or var == 'g' or var == 'p_s' or \
                           var == 'lambda_dem')
                        and P_tot is None):
                        
                        return UserWarning('\n Output Key Error: No '+\
                                  'Assimilate \n the method: '+key+\
                                  ' requires P to be estimated, if RED_init'+\
                                  ' is being used employ a P_or_gamma method'+\
                                  ' or an applicable "assimilate '+\
                                  'mortality"/"P_gamma" observation.')
                        
                    elif ((var=='P' or var == 'G' or var == 'P_s' or \
                           var == 'p' or var == 'g' or var == 'p_s')\
                          and P_tot is not None):
                        
                        N_p_sum = np.sum(N_m[:,:]*g_scaling[:,:],axis=1)
                        p_0 = np.divide(P_tot[:],N_p_sum[:]) 
                        p = p_0[:,None]* g_scaling[:]
                        g = (1-np.array(alpha[:])[:,None])*p[:,:]
                        p_s = np.array(alpha[:])[:,None] * p[:,:]

                        
                    first_variable = var
                    first_variable_check = True

                    if  second_variable_check == True:
                        
                        break
                    
                    
                if var2.find(var) != -1:
                                        
                    if ((var=='P' or var == 'G' or var == 'P_s' or \
                           var == 'p' or var == 'g' or var == 'p_s' or \
                            var == 'lambda_dem')
                        and P_tot is None):
                        
                        return UserWarning('\n Output Key Error: No '+\
                                  'Assimilate \n the method: '+key+\
                                  ' requires P to be estimated, if RED_init'+\
                                  ' is being used employ a P_or_gamma method'+\
                                  ' or an applicable "assimilate '+\
                                  'mortality"/"P_gamma" observation.')
                        
                    elif ((var=='P' or var == 'G' or var == 'P_s' or \
                           var == 'p' or var == 'g' or var == 'p_s')\
                          and P_tot is not None):
                        
                        N_p_sum = np.sum(N_m[:,:]*g_scaling[:,:],axis=1)
                        p_0 = np.divide(P_tot[:],N_p_sum[:]) 
                        p = p_0[:,None]* g_scaling[:]
                        g = (1-np.array(alpha[:])[:,None])*p[:,:]
                        p_s = np.array(alpha[:])[:,None] * p[:,:]

                    second_variable = var
                    second_variable_check = True
                    
                    if  first_variable_check == True:
                        
                        break
                    
            else:
                
                if key.find(var) != -1:
                    
                    if ((var=='P' or var == 'G' or var == 'P_s' or \
                           var == 'p' or var == 'g' or var == 'p_s' or\
                           var == 'lambda_dem') and P_tot is None):
                        
                        raise UserWarning('\n Output Key Error: No '+\
                                  'Assimilate \n the method: '+key+\
                                  ' requires P to be estimated, if RED_init'+\
                                  ' is being used employ a P_or_gamma method'+\
                                  ' or an applicable "assimilate '+\
                                  'mortality"/"P_gamma" observation.')
                        
                    elif ((var=='P' or var == 'G' or var == 'P_s' or \
                           var == 'p' or var == 'g' or var == 'p_s' or\
                           var == 'lambda_dem') and P_tot is not None):
                        
                        N_p_sum = np.sum(N_m[:,:]*g_scaling[:,:],axis=1)
                        p_0 = np.divide(P_tot[:],N_p_sum[:]) 
                        p = p_0[:,None]* g_scaling[:]
                        g = (1-np.array(alpha[:])[:,None])*p[:,:]
                        p_s = np.array(alpha[:])[:,None] * p[:,:]

                    variable = var
                    variable_check = True
                    
                    break

        if output_type == 'r':
        
            if (first_variable_check == False or \
                second_variable_check == False):
                
                raise UserWarning('\n Output Key Error: Invalid output '+\
                                  'requested: \n - if percentile output is'+\
                                   ' choosen there must be two valid'+\
                                   ' variables choosen; "N", "M", "P", "H",'+\
                                   ' "nu", "m", "h", "p", "g", "p_s", "a"'+\
                                   'and arrange accordingly: "X_Y_r", where'+\
                                   ' "r" is the digits of the percentile')
                
            if first_variable == 'N':
                
                X = N_m[:,:]
                
            elif first_variable == 'M':
                
                X  = N_m[:,:] * m[:,:]

            elif first_variable == 'P':
            
                X = N_m[:,:] * p[:,:]

            elif first_variable == 'G':
                
                X = N_m[:,:] * g[:,:]

            elif first_variable == 'P_s':
                
                X = N_m[:,:] * p_s[:,:]

            elif first_variable == 'gamma':
                
                if (gamma is None and gamma_base is None):
                    
                    raise UserWarning('\n Output Key Error: No '+\
                                  'Mortality \n the method: '+key+\
                                  ' requires gamma or gamma_base '+\
                                  'to be provided.')

                elif (gamma is None and gamma_base is not None):

                    gamma = np.zeros((I,max(J)))
                    
                    for i in range(0,I):
                        
                        gamma[i,:] = gamma_base[i]

                X = gamma[:,:]

            elif first_variable == 'H':
                
                X = N_m[:,:] * h[:,:]

            elif first_variable == 'nu':

                X = N_m[:,:] * a[:,:]

            elif first_variable == 'f_par':
                
                if f_par is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'f_par \n the method: '+key+\
                          ' requires f_par'+\
                          'to be provided within the input.')
                    
                X = f_par[:,:]

            elif first_variable == 'F_in':
                
                if F_in is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'F_in \n the method: '+key+\
                          ' requires F_in'+\
                          'to be provided within the input.')
            
                X = F_in[:,:]

            elif first_variable == 'F_out':
                
                if F_out is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'F_out \n the method: '+key+\
                          ' requires F_out'+\
                          'to be provided within the input.')
            
                X = F_out[:,:]
            
            elif first_variable == 'm':

                X = m[:,:]

            elif first_variable == 'h':
                
                X = h[:,:]

            elif first_variable == 'p':

                X = p[:,:]

            elif first_variable == 'g':
                
                X = g[:,:]
                    
            elif first_variable == 'p_s':
                
                X = p_s[:,:]

            elif first_variable == 'a':
                
                X = a[:,:]

            if second_variable == 'N':
                
                Y = N_m[:,:]
                
            elif second_variable == 'M':
                
                Y  = N_m[:,:] * m[:,:]

            elif second_variable == 'P':
            
                Y = N_m[:,:] * p[:,:]

            elif second_variable == 'G':
                
                Y = N_m[:,:] * g[:,:]

            elif second_variable == 'P_s':
                
                Y = N_m[:,:] * p_s[:,:]

            elif second_variable == 'gamma':
                
                if (gamma is None and gamma_base is None):
                    
                    raise UserWarning('\n Output Key Error: No '+\
                                  'Mortality \n the method: '+key+\
                                  ' requires gamma or gamma_base '+\
                                  'to be provided.')

                elif (gamma is None and gamma_base is not None):

                    gamma = np.zeros((I,max(J)))
                    
                    for i in range(0,I):
                        
                        gamma[i,:] = gamma_base[i]

                Y = gamma[:,:]

            elif second_variable == 'H':
                
                Y = N_m[:,:] * h[:,:]

            elif second_variable == 'nu':

                Y = N_m[:,:] * a[:,:]
 
            elif second_variable == 'f_par':
                
                if f_par is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'f_par \n the method: '+key+\
                          ' requires f_par'+\
                          'to be provided within the input.')
                    
                Y = f_par[:,:]

            elif second_variable == 'F_in':
                
                if F_in is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'F_in \n the method: '+key+\
                          ' requires F_in'+\
                          'to be provided within the input.')
            
                Y = F_in[:,:]

            elif second_variable == 'F_out':
                
                if F_out is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'F_out \n the method: '+key+\
                          ' requires F_out'+\
                          'to be provided within the input.')
            
                Y = F_out[:,:]
            

            elif second_variable == 'm':
                
                Y = m[:,:]

            elif second_variable == 'h':
                
                Y = h[:,:]

            elif second_variable == 'p':
                
                Y = p[:,:]

            elif second_variable == 'g':
                
                Y = g[:,:]
                    
            elif second_variable == 'p_s':
                
                Y = p_s[:,:]

            elif second_variable == 'a':
                
                Y = a[:,:]         
    
        else:
            
            if variable_check == False:
                
                raise UserWarning('\n Output Key Error: Invalid output '+\
                                    'requested: \n - if percentile output is'+\
                                    ' choosen there must be one valid'+\
                                    ' variable; "N", "M", "P", "H",'+\
                                    ' "nu", "m", "h", "p", "g", "p_s", "a"'+\
                                    'and arrange accordingly: "X_Y_r", where'+\
                                    ' "r" is the digits of the percentile.')
    
            if variable == 'N':
                
                Y = N_m[:,:]
                
            elif variable == 'M':
                
                Y  = N_m[:,:] * m[:,:]
    
            elif variable == 'P':
                
                Y = N_m[:,:] * p[:,:]
    
            elif variable == 'G':
                
                Y = N_m[:,:] * g[:,:]
    
            elif variable == 'P_s':
                
                Y = N_m[:,:] * p_s[:,:]

            elif variable == 'gamma':
                
                if (gamma is None and gamma_base is None):
                    
                    raise UserWarning('\n Output Key Error: No '+\
                                  'Mortality \n the method: '+key+\
                                  ' requires gamma or gamma_base '+\
                                  'to be provided.')

                elif (gamma is None and gamma_base is not None):

                    gamma = np.zeros((I,max(J)))
                    
                    for i in range(0,I):
                        
                        gamma[i,:] = gamma_base[i]

                Y = gamma[:,:]
    
            elif variable == 'H':
                
                Y = N_m[:,:] * h[:,:]
    
            elif variable == 'nu':
    
                Y = N_m[:,:] * a[:,:]

 
            elif variable == 'f_par':
                
                if f_par is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'f_par \n the method: '+key+\
                          ' requires f_par'+\
                          'to be provided within the input.')
                    
                Y = f_par[:,:]

            elif variable == 'F_in':
                
                if F_in is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'F_in \n the method: '+key+\
                          ' requires F_in'+\
                          'to be provided within the input.')
            
                Y = F_in[:,:]

            elif variable == 'F_out':
                
                if F_out is None:
                    
                    raise UserWarning('\n Output Key Error: No '+\
                          'F_out \n the method: '+key+\
                          ' requires F_out'+\
                          'to be provided within the input.')
            
                Y = F_out[:,:]
    
            elif variable == 'm':
                
                Y = m[:,:]
    
            elif variable == 'h':
                
                Y = h[:,:]
    
            elif variable == 'p':
                
                Y = p[:,:]
    
            elif variable == 'g':
                
                Y = g[:,:]
                    
            elif variable == 'p_s':
                
                Y = p_s[:,:]
    
            elif variable == 'a':
                
                Y = a[:,:]

    if output_type == 'tot':
        
        return get_sum(I,J,Y)
        
    elif output_type == 'm':
        
        return Y
        
    elif output_type == 'mean':
   
        return get_mean(I,J,N_m,Y)
        
    elif output_type == 'norm':

        return get_normalised(I,J,Y)
        
        
    elif output_type == 'r':
        
        return get_percentile(I,J,X,Y,percentile)

    
    


def get_percentile(I,J,X,Y,percentile):
    """
    Returns the percentile value of a variable across each PFT with respect
    to another variables distribution.
    
    Inputs:
        
         - I, Number of PFTs (-)
            
         - J, number of mass classes. (-)
          
         - Y, the variable which the distribution is taken 
           (Normaly number density). (-)
         
         - X, variable which the percentile value is taken, which has length 
           equal to the maximum number mass classes. (-)
           
         - percentile, value between 0.0 and 100.0. (%)

    Outputs:
        
        - Y_percentile, the percentile of Y with respect to the distriubiton
          of the percentile of X. (-)
    
    """
    X_percentile = np.zeros(I)

    for i in range(0,I):
        
        j = np.arange(0,J[i],dtype=int)
        Y_percentile = np.percentile(Y[i,j],percentile) # Assumes linear
        X_percentile[i] = np.interp(Y_percentile,Y[i,j],X[i,j]) # Interpolate
        
    return X_percentile
   
def get_sum(I,J,Y):
    
    """
    From the variable number of classes across each PFT, return the summation 
    of a particular value across the mass classes.
    
    Inputs:
        
         - I, Number of PFTs (-)
            
         - J, number of mass classes. (-)
          
         - Y, the variable being summed, which has length equal to the maximum
              number mass classes. (-)
         
         
    Outputs:
        
         - Y_sum, the summation of the Y variable across each of the mass
           classes. (-)
           
    """    
    Y_sum = np.zeros(I)
    
    for i in range(0,I):
        
        j = np.arange(0,J[i],dtype=int)
        Y_sum[i] = sum(Y[i,j])

    return Y_sum
    

def get_mean(I,J,X,Y):
    
    """
    Get the mean of one variable with respect to the distribution of another,
    across all PFTs.
    
    Inputs:
        
         - I, Number of PFTs (-)
            
         - J, number of mass classes. (-)
          
         - X, the variable which the distribution is taken 
           (Normaly number density). (-)
         
         - Y, variable which the mean is taken, which has length equal to the 
           maximum number mass classes. (-)

    Outputs:
        
        - Y_mean, the mean of Y with respect to the distriubiton of X. (-)
    
    """
    Y_mean = np.zeros(I)
    
    for i in range(0,I):
        
        j = np.arange(0,J[i],dtype=int)
        X_sum = sum(X[i,j])

        if X_sum == 0:
            
            f_X = float('nan')
            
        else:
            
            f_X = X[i,j]/X_sum
        
        Y_mean[i] = sum(f_X * Y[i,j]) 
    
    return Y_mean
    
def get_normalised(I,J,Y):

    """
    Return the normalised array over the variable Y.
    
    Inputs:
        
         - I, Number of PFTs (-)
            
         - J, number of mass classes. (-)
          
         - Y, the variable to normalise. (-)
         
    Outputs:
        
        -  Y_norm, the mean of Y with respect to the distriubiton of X. (-)
    
    """
    Y_norm = np.zeros((I,max(J)))
    
    for i in range(0,I):
        
        j = np.arange(0,J[i],dtype=int)
        Y_sum = sum(Y[i,j])
        
        if Y_sum == 0:
            
            Y_norm[i] = float('nan')
            
        else:

            Y_norm[i,j] = Y[i,j]/Y_sum
    
    return Y_norm

def quicksort_record_row_col(A_flat,i_cord,j_cord,start,end):
    """
    Quicksort algrorithim designed to record the ith row and jth colunm and of a
    flattened matrix. 
    
    Inputs:
            
        - A_flat, a flattened matrix of IxJ length.
            
        - i_cord, records the integer row position of the corresponding 
          element in A_flat. *
              
        - j_cord, records the integer column position of the corresponding
          element in A_flat.*
              
    Outputs:
            
            - A_flat, ordered from smallest to largest.
            
            - i_cord, ordered relative from smallest to largest A_flat.
              
            - j_cord, ordered relative from smallest to largest A_flat. 
              
              
    Called:  "COMPETITION", Overloaded function - Calls itself.
          
    
    * Such that when element l, is called in A_flat we get the same value
      from the non-flatten matrix, A -
          
          A_flat[l] = A[i_cord[l],j_cord[l]]
          
          
    
    """
    
    pivot = A_flat[int((start+end)/2)]
              
    l_up = start
    
    l_down = end
    
    while True:
        
        while A_flat[l_up] < pivot:
            
            l_up = l_up + 1
            
        
        while pivot < A_flat[l_down]:
            l_down = l_down - 1
        
            
        if l_up >= l_down:
            
            break

        A_flat[l_up], A_flat[l_down] =  A_flat[l_down], A_flat[l_up]

        i_cord[l_up], i_cord[l_down] =  i_cord[l_down],i_cord[l_up]

        j_cord[l_up], j_cord[l_down] =  j_cord[l_down],j_cord[l_up]

        l_up = l_up + 1
        l_down = l_down - 1
    
    if start < l_up - 1:
        
        A_flat, i_cord, j_cord = quicksort_record_row_col(A_flat,i_cord,\
                                                          j_cord,start,l_up-1)
    
    if l_down + 1 < end:
        
        A_flat, i_cord, j_cord = quicksort_record_row_col(A_flat,i_cord,\
                                                          j_cord,l_down+1,end)

    return A_flat,i_cord,j_cord

#
#def find_optimum_mult(mu=0.25,method='nu',tol=1.48e-08,pft_params=None,mult_init=[20.,1.05]):
#
#    """
#    For given input parameters and find the optimal bin logarithmic scaling to 
#    minimise the difference between analytical and discrete forms.
#    
#    Inputs:
#        
#        - mu, the estimated turnover ratio, defined as 
#          gamma_init * m_init/g_init (-)
#          
#        - method, the quantity we are trying to find the optimal bin scaling 
#          for, currently coverage (nu) is completed. (-)
#         
#        - pft_params, a tuple containing the arguments for optimising a 
#                      specific PFT at a specific method. (-)
#          
#        - tol, tolerance for termination of the numerical root finding process.
#          (-)
#    Output:
#        
#        - S, the optimal class scaling (-).
#    
#    """
#    if method == 'nu':
#        
#        if pft_params == None:
#            I, _, _, J, mult, _, _,_,\
#                       _, alpha,gamma_init, phi_g, phi_h, phi_a, nu_min = \
#                       PFT_VALUES()
#        
#                
#                       
#            
#            for i in range(0,I):
#                
#                        
#                args = (mu,J[i],phi_g[i])
#                
#                mult[i] = newton(Za,mult_init,args = args, fprime = dZa_dS, tol = tol)
#                
#                if mult[i] != mult[i]:
#                    mult_bounds = [mult_init[-1],mult_init[0]]
#                    mult[i] = gss(Za,args,mult_bounds)
#                    
#        
#        else:    
#            
#            J = pft_params[0]
#            phi_g = pft_params[1]
#            args = (mu,J,phi_g)
#            mult = newton(Za,dZa_dS,args,mult_init,tol = tol)
#            if mult > 1e10:
#                mult_bounds = [mult_init[-1],mult_init[0]]
#                mult = gss(Za,args,mult_bounds)
#
#        
#    return mult
              
def newton(function,derivative,args,x_init,tol=1.48e-08):
    """
    Modified Newton root finding algorithim, customised to ouput nan value if 
    no root is found. Uses multiple x_init.
    
    Inputs:
        
        - function, a definition of a process with args and xinit as intputs, 
          i.e. function(x,args[0],...,args[-1]). (-)
          
        - derivative, the derivative of the function with respect to x. Takes 
          identicial inouts as the function, i.e. derivative(x,args[0],...
          ,args[-1]). (-)
          
        - args, the other arguments of the function besides x. A tuple of n 
          arguements. (-)
          
        - x_init, the initial guess of the root, can be an array of guesses.
         (-)
        
        - tol, the required tolerance of error. (-)
        
    Outputs:
        
        - x, the numerical approximation of the root.
        
    
    """
    

    if isinstance(x_init,Sequence) == False:
        
        array_of_guesses = False
        
    else:
        num_guess = len(x_init)
    
        array_of_guesses = True

        

    if array_of_guesses == True:
        
        for i in range(0,num_guess):

            x = x_init[i] 

            error_abs = float("inf")
    
            error_check = False
       
            convergence_count = 0
            
            while error_check == False:

                error_prev = error_abs
                
                myargs = (x,) + args
            
                f_x = function(*myargs)
            
                df_dx = derivative(*myargs)
                
                if df_dx != 0:
                    
                    x = x - f_x/df_dx
                
                elif f_x < 0:
                    
                    x = float('inf')
                    
                    break
                    
                elif f_x > 0:
                    
                    x = -float('inf')
                    
                    break
                    
                else:
                    
                    x = float('nan')
                    
                    break
                
                
                
                myargs = (x,) + args
             
                error_abs = abs(function(*myargs) - 0.)

                if error_abs <= tol:
                    
                    error_check = True
                    
                    break
                
                if error_prev <= error_abs:
                    
                    
                    convergence_count = convergence_count + 1
                    
                    if convergence_count == 50:
                        
                        x = float("nan")
                        
                        break                    
                    
                else:
                    
                    convergence_count = 0
                    
                    
                    
            if error_check == True:
                
                break
            
    else:
        
        x = x_init
        error_abs = float('inf')
        error_check = False
        convergence_count = 0
        count = 0
        
        while error_check == False:

            if count > 0:
                
                error_prev2 = error_prev

            else:
                
                error_prev2 = float('inf')
                
            error_prev = error_abs
            myargs = (x,) + args
            f_x = function(*myargs)
            df_dx = derivative(*myargs)

            if df_dx != 0:
                
                x = x - f_x/df_dx
            
            elif f_x < 0:
                
                x = float('inf')
                
                break
                
            elif f_x > 0:
                
                x = -float('inf')
                
                break
                
            else:
                x = float('nan')

                break
            
            myargs = (x,) + args
            error_abs = abs(function(*myargs) - 0.)
            
            if error_abs <= tol:
                
                error_check = True
                
                break
            
            if (error_prev <= error_abs or error_prev2 <= error_abs or x != x):
                
                convergence_count = convergence_count + 1
                
                if convergence_count == 50:
                    
                    x = float("nan")
                    
                    break                    
                
            else:
                
                convergence_count = 0
            
                
            if count == 1000:
                
                x = float('nan')
                
                break
                
            count = count + 1 
            
            
                        
    return x
        
        
        


def gss(function,args,x_bounds,tol=1.48e-08):
    """
    Golden-section search  (GSS) Algorithim for solving for the local extremum
    of a function. Used when the newton method returns a impractical result 
    (mult>>1) for optimising the bin scalin0g.
            
    Inputs:
        
        - function, a definition of a process with args and xinit as intputs, 
          i.e. function(x,args[0],...,args[-1]). (-)
          
        - args, the other arguments of the function besides x. A tuple of n 
          arguements. (-)

        - x_bounds, where to check for the root between x. An array with [x_1
          ,x_2], where x_1 < x_2. (-)
                
        - tol, the required tolerance of error. (-)
        
    Outputs:
        
        - x, the numerical approximation of the root.

        
    """
    
    error_prev = float("inf")

    error_check = False
   
    convergence_count = 0
    
    gr = (5.0**(0.5) + 1.0)/2.0
    
    x_1 = x_bounds[0]

    x_2 = x_bounds[1]
    x_a = x_2 - (x_2 - x_1)/gr
    x_b = x_2 + (x_2 - x_2)/gr 
    
    while error_check == False:
        
        error_abs = abs(x_a - x_b)

        if error_prev <= error_abs:
            convergence_count = convergence_count + 1           
            if convergence_count == 50:               
                x_a = float("nan")
                x_b = float("nan")
                break                                
        else:
            
            convergence_count = 0
                       
        error_prev = error_abs
                
        if error_abs <= tol:
            
            error_check = True
            
            break
        
        
        error_prev = error_abs
        
        myargs_x1 = (x_a,) + args
        myargs_x2 = (x_b,) + args
        
        f_xa = function(*myargs_x1)
        f_xb = function(*myargs_x2)
        if f_xa < f_xb:
            x_2 = x_b            
        else:
            x_1 = x_a              

    return (x_a + x_b)/2.0

def LOAD_PFT_VALUES(load='all'):
    """
    Based on input conditions loads and returns warrented values from 
    'PFT_values.pk1' held in 'RED_MODULE/PFT_STORE/' 
    
    Inputs: 
        
        - load, what is needed when called? load = 'all', or a tuple of 
          keynames from the dictonary of PFT_values.pk1 (-)
    
    Outputs: 
        
        if load == 'all', returns all variables saved in 'PFT_values.pk1':
            
            Keynames:
            
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
            
            - c_pft, competition coefficient. How much shading among plant 
              functional groups. (-)
            
            - alpha, the fraction of tree productivity going into seedling 
              production. (-) 
            
            - phi_g, growth-mass scaling power. (-)
            
            - phi_h, height-mass scaling power. (-)
            
            - phi_a, crown area-mass scaling power. (-)
            
            - nu_min, minimum vegetation fraction. (-)
                        
            - i_ordered, ordered relative from smallest to largest height. (-)
              
            - j_ordered, ordered relative from smallest to largest height. (-)
            

        
    """
    
    script_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),\
                                              '..', 'PFT_STORE'))
    script_name = '/PFT_values.pk1'
    
    with open(script_dir+script_name,'rb') as f:
        pft_dict = pk.load(f)
    
    if type(load) is not tuple:
        
        if load is 'all':
            
            return pft_dict.values()
            
        else:
            
            raise UserWarning('\n Error in Loading PFT Values: \n\n - '+\
                              'load must equal either a tuple of keywords or'+\
                              ' must be a string of "all".')
    else:             
        
        pft_dict_subset = [pft_dict[key] for key in load]
        
        return pft_dict_subset

        
        
        
#TODO - Re add in optimum method for finding the bin scaling. 
#def Za(mult,mu,J,phi_g):
#    
#    """
#    Solving for the optimum bin scaling, S, for a given mu to minimise the 
#    difference between analytical and numerical forms of the equilibrium
#    solution. The function describes the condition needed to be true for a 
#    maximum. d_[nu_analytical - nu_discrete]/dS = 0.
#    
#    Inputs:
#        
#        - mult, multiplicative bin scaling parameter. (-)
#        
#        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
#        
#        - J, the number of mass classes.
#        
#        - phi_g, the growth mass scaling parameter. (-)
#    
#    Outputs:
#        
#        - 
#       
#    
#    """
#    
#    _, X_N, X_G, _ = X(J,mult,mu,phi_g=phi_g)
#    _, dX_N_dS, dX_G_dS, _ =  dX_dS(J,mult,mu,phi_g=phi_g)
#    
#    Z = X_G * dX_N_dS - X_N * dX_G_dS
#    
#    return Z
#    
#
#def dZa_dS(mult,mu,J,phi_g):
#    
#    """
#    Solving for the optimum bin scaling, S, for a given mu to minimise the 
#    difference between analytical and numerical forms of the equilibrium
#    solution. The function describes the condition needed to be true for a 
#    maximum. d_[nu_analytical - nu_discrete]/dS = 0.
#    
#    Inputs:
#        
#        - mult, multiplicative bin scaling parameter. (-)
#        
#        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
#        
#        - J, the number of mass classes.
#        
#        
#        - phi_g, the growth mass scaling parameter. (-)
#    
#    Outputs:
#        
#        - 
#       
#    
#    """
#   
#    _, X_N, X_G,_ = X(J,mult,mu,phi_g=phi_g)
#    
#    _, dX_N2_dS2, dX_G2_dS2, _  = dX2_dS2(J,mult,mu,phi_g=phi_g)
#    
#    
#    dZ_dS = X_G*dX_N2_dS2 - dX_G2_dS2*X_N
#    
#    return dZ_dS
#
#
#
#
#
#def dX_dS(J,mult,mu,phi_g = 0.0, phi_a = 0.0):
#    
#    """
#    Function dX_dS is used to determine the optimal scaling to fit the 
#    continous solutions.
#    
#    Inputs:
#        
#        - J, the number of mass classes. (-)
#        
#        - mult, multiplicative bin scaling parameter. (-)
#        
#        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
#        
#        - phi_g, the growth mass scaling parameter. (-)
#        
#        - phi_a, the area mass scaling parameter. (-)
#        
#        
#    Outputs:
#
#        - dX_M_dS, summation of the mass scaling bin scaling differential 
#                   over each mass class. (-)        
#
#        - dX_N_dS, summation of the scaling of the number density bin scaling
#                   differential over each mass class. (-)   
#         
#        - dX_G_dS, summation of the scaling of the growth bin scaling 
#                   differential over each mass class. (-)   
#        
#    
#        - dX_nu_dS, summation of the vegetation fraction bin scaling 
#                    differential over each mass class. (-)  
#        
#    """
#    
#    dX_M_dS = 0.
#    
#    dX_N_dS = 0.
#    
#    dX_G_dS = 0.
#    
#    dX_nu_dS = 0.
#    
#    
#    for i in range(0,J):
#        
#        prod = 1
#        series = 0       
#        s_m = mult ** (i)
#        ds_m_dS = i * mult ** (i-1)
#        s_g = mult ** (i*phi_g)
#        ds_g_dS = i*phi_g *mult**(- 1)
#        s_nu = mult ** (i*phi_a)
#        ds_nu_dS = i*phi_a *mult**(- 1)
#
#        for j in range(1,i + 1):
#            
#            if j < J - 1:
#                
#                kj = k_j(j,mult,mu,phi_g)
#                dkjdS = dk_j_dS(j,mult,mu,phi_g)
#            
#            elif j == J - 1 :
#                
#                kj =  mult**((j-1)*(phi_g-1))/(mu*(mult-1))
#                dkjdS =  ((j-1)*(phi_g-1) * (1-mult**(-1)) - 1)/(mult-1.) * kj
#           
#                
#            prod = prod * kj            
#            series = series + dkjdS/kj
#            
#        dX_M_dS = dX_M_dS + prod * (series+ds_m_dS) * s_m 
#        dX_N_dS = dX_N_dS + prod * series
#        dX_G_dS = dX_G_dS + prod * (series+ds_g_dS)* s_g
#        dX_nu_dS = dX_nu_dS + prod * (series+ds_nu_dS)* s_nu
#        
#        
#        
#    return dX_M_dS, dX_N_dS, dX_G_dS, dX_nu_dS
#        
#    
#    
#def dX2_dS2(J,mult,mu,phi_g = 0.0, phi_a = 0.0):
#    
#    """
#    Function dX2_dS2 is used to determine the optimal scaling to fit the 
#    continous solutions. It is the double differntial of X with respect to the 
#    bin scaling.
#    
#    Inputs:
#        
#        - J, the number of mass classes. (-)
#        
#        - mult, multiplicative bin scaling parameter. (-)
#        
#        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
#        
#        - phi_g, the growth mass scaling parameter. (-)
#        
#        - phi_a, the area mass scaling parameter. (-)
#        
#        
#    Outputs:
#
#        - dX_M2_dS2, summation of the mass scaling bin scaling differential 
#                   over each mass class. (-)        
#
#        - dX_N2_dS2, summation of the scaling of the number density bin scaling
#                   differential over each mass class. (-)   
#         
#        - dX_G2_dS2, summation of the scaling of the growth bin scaling 
#                   differential over each mass class. (-)   
#        
#    
#        - dX_nu2_dS2, summation of the vegetation fraction bin scaling 
#                    differential over each mass class. (-)  
#        
#    """
#    
#    dX_M2_dS2 = 0.
#    
#    dX_N2_dS2 = 0.
#    
#    dX_G2_dS2 = 0.
#    
#    dX_nu2_dS2 = 0.
#    
#    
#    for i in range(0,J):
#        
#        prod = 1
#        series_M = 0
#        series = 0       
#        series_g = 0
#        series_nu = 0
#        
#        s_m = mult ** (i)
#        ds_m_dS = i * mult ** (i-1)
#        s_g = mult ** (i*phi_g)
#        ds_g2_dS2 = i*phi_g*mult**(-2)*(i*phi_g-1)
#        s_nu = mult ** (i*phi_a)
#        ds_nu_dS = i*phi_a *mult**(i*phi_a - 1)
#
#        for j in range(1,i + 1):
#            
#            if j < J - 1:
#                
#                kj = k_j(j,mult,mu,phi_g)
#                dkjdS = dk_j_dS(j,mult,mu,phi_g)
#                dkj2dS2 = dk_j2_dS2(j,mult,mu,phi_g)
#            
#            elif j == J - 1 :
#                
#                kj =  mult**((j-1)*(phi_g-1))/(mu*(mult-1))
#                dkjdS =  ((j-1)*(phi_g-1) * (1-mult**(-1)) - 1)/(mult-1.) * kj
#                dkj2dS2 = dkjdS**2/kj + (kj/(mult-1.))*((j-1)*(phi_g-1.)\
#                                         /(mult**2)-dkjdS/kj)
#
#            prod = prod * kj            
#            series = series + dkj2dS2/kj
#            series_g = series_g + (2*j*phi_g*mult**(-1)*dkjdS - dkj2dS2)/kj 
#            
#        dX_M2_dS2 = dX_M2_dS2 + prod * series * s_m + ds_m_dS *prod
#        dX_N2_dS2 = dX_N2_dS2 + prod * series
#        dX_G2_dS2 = dX_G2_dS2 + s_g * prod *(ds_g2_dS2 + series_g) 
#        dX_nu2_dS2 = dX_nu2_dS2 + prod * series * s_nu + ds_nu_dS *prod
#        
#        
#        
#    return dX_M2_dS2, dX_N2_dS2, dX_G2_dS2, dX_nu2_dS2
#        
#    
#            
#        
#def dk_j_dS(j,mult,mu,phi_g):
#    """
#    Outputs the scaling charateristics gradient of each mass class against S.
#    
#    Inputs: 
#    
#    - j, the given mass class. (-)
#    
#    - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
#    
#    - phi_g, the growth mass scaling parameter. (-)
#    
#    Outputs:
#    
#    - dk_j_dS, the relative change of fraction population in the current
#      class against S. (-)
#    
#    """
#    dk_j_dS = k_j(j,mult,mu,phi_g) * ((j-1.) * (phi_g-1.)/mult - \
#                  k_j(j,mult,mu,phi_g)*(j*(phi_g-1)*mult**(phi_g-2)+\
#                      mu*mult**((j-1)*(1-phi_g))))
#
#    return dk_j_dS
#    
#    
#def dk_j2_dS2(j,mult,mu,phi_g):
#    """
#    The double derivative of the scaling charateristic function with respect to
#    the bin scaling.
#        Inputs: 
#    
#    - j, the given mass class. (-)
#    
#    - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
#    
#    - phi_g, the growth mass scaling parameter. (-)
#    
#    Outputs:
#    
#    - dk_j2_dS2, the relative change of fraction population in the current
#      class against the double derivative of S. (-)
#    
#    """
#    kj = k_j(j,mult,mu,phi_g)
#    dkjdS = dk_j_dS(j,mult,mu,phi_g)
#    r = dkjdS/kj
#    dk_j2_dS2 = kj * (-r*(j-1.)*(phi_g-1.)*mult**(-1)-(j-1.)*(phi_g\
#                      -1.)*mult**(-2.)-kj*(phi_g-1.)*mult**(-1.)*(j*(phi_g-1.)\
#                      *mult**(phi_g-2.)-mu*(j-1.)*mult**((j-1.)*(1.-phi_g))))
#    
#    return dk_j2_dS2