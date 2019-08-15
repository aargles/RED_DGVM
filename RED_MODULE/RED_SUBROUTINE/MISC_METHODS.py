"""
Python 2.7, Encoding: UTF-8

Created on Mon Jul  2

@author: A.Argles (aa760@exeter.ac.uk)
"""


#%% Script hold miscellaneous methods. 
from STORE import PFT_VALUES
from EQUILIBRIUM import nu_eq_discrete
from EQUILIBRIUM_ANALYTICAL import nu_eq_analy
from EQUILIBRIUM_METHODS import X, k_j
from collections import Sequence

def find_optimum_mult(mu=0.25,method='nu',tol=1.48e-08,pft_params=None,mult_init=[20.,1.05]):

    """
    For given input parameters and find the optimal bin logarithmic scaling to 
    minimise the difference between analytical and discrete forms.
    
    Inputs:
        
        - mu, the estimated turnover ratio, defined as 
          gamma_init * m_init/g_init (-)
          
        - method, the quantity we are trying to find the optimal bin scaling 
          for, currently coverage (nu) is completed. (-)
         
        - pft_params, a tuple containing the arguments for optimising a 
                      specific PFT at a specific method. (-)
          
        - tol, tolerance for termination of the numerical root finding process.
          (-)
    Output:
        
        - S, the optimal class scaling (-).
    
    """
    if method == 'nu':
        
        if pft_params == None:
            I, _, _, J, mult, _, _,_,\
                       _, alpha,gamma_init, phi_g, phi_h, phi_a, nu_min = \
                       PFT_VALUES()
        
                
                       
            
            for i in range(0,I):
                
                        
                args = (mu,J[i],phi_g[i])
                
                mult[i] = newton(Za,mult_init,args = args, fprime = dZa_dS, tol = tol)
                
                if mult[i] != mult[i]:
                    mult_bounds = [mult_init[-1],mult_init[0]]
                    mult[i] = gss(Za,args,mult_bounds)
                    
        
        else:    
            
            J = pft_params[0]
            phi_g = pft_params[1]
            args = (mu,J,phi_g)
            mult = newton(Za,dZa_dS,args,mult_init,tol = tol)
            if mult > 1e10:
                mult_bounds = [mult_init[-1],mult_init[0]]
                mult = gss(Za,args,mult_bounds)

        
    return mult
              
               
#%% Secondary Functions
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




def Za(mult,mu,J,phi_g):
    
    """
    Solving for the optimum bin scaling, S, for a given mu to minimise the 
    difference between analytical and numerical forms of the equilibrium
    solution. The function describes the condition needed to be true for a 
    maximum. d_[nu_analytical - nu_discrete]/dS = 0.
    
    Inputs:
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - J, the number of mass classes.
        
        - phi_g, the growth mass scaling parameter. (-)
    
    Outputs:
        
        - 
       
    
    """
    
    _, X_N, X_G, _ = X(J,mult,mu,phi_g=phi_g)
    _, dX_N_dS, dX_G_dS, _ =  dX_dS(J,mult,mu,phi_g=phi_g)
    
    Z = X_G * dX_N_dS - X_N * dX_G_dS
    
    return Z
    

def dZa_dS(mult,mu,J,phi_g):
    
    """
    Solving for the optimum bin scaling, S, for a given mu to minimise the 
    difference between analytical and numerical forms of the equilibrium
    solution. The function describes the condition needed to be true for a 
    maximum. d_[nu_analytical - nu_discrete]/dS = 0.
    
    Inputs:
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - J, the number of mass classes.
        
        
        - phi_g, the growth mass scaling parameter. (-)
    
    Outputs:
        
        - 
       
    
    """
   
    _, X_N, X_G,_ = X(J,mult,mu,phi_g=phi_g)
    
    _, dX_N2_dS2, dX_G2_dS2, _  = dX2_dS2(J,mult,mu,phi_g=phi_g)
    
    
    dZ_dS = X_G*dX_N2_dS2 - dX_G2_dS2*X_N
    
    return dZ_dS





def dX_dS(J,mult,mu,phi_g = 0.0, phi_a = 0.0):
    
    """
    Function dX_dS is used to determine the optimal scaling to fit the 
    continous solutions.
    
    Inputs:
        
        - J, the number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - phi_a, the area mass scaling parameter. (-)
        
        
    Outputs:

        - dX_M_dS, summation of the mass scaling bin scaling differential 
                   over each mass class. (-)        

        - dX_N_dS, summation of the scaling of the number density bin scaling
                   differential over each mass class. (-)   
         
        - dX_G_dS, summation of the scaling of the growth bin scaling 
                   differential over each mass class. (-)   
        
    
        - dX_nu_dS, summation of the vegetation fraction bin scaling 
                    differential over each mass class. (-)  
        
    """
    
    dX_M_dS = 0.
    
    dX_N_dS = 0.
    
    dX_G_dS = 0.
    
    dX_nu_dS = 0.
    
    
    for i in range(0,J):
        
        prod = 1
        series = 0       
        s_m = mult ** (i)
        ds_m_dS = i * mult ** (i-1)
        s_g = mult ** (i*phi_g)
        ds_g_dS = i*phi_g *mult**(- 1)
        s_nu = mult ** (i*phi_a)
        ds_nu_dS = i*phi_a *mult**(- 1)

        for j in range(1,i + 1):
            
            if j < J - 1:
                
                kj = k_j(j,mult,mu,phi_g)
                dkjdS = dk_j_dS(j,mult,mu,phi_g)
            
            elif j == J - 1 :
                
                kj =  mult**((j-1)*(phi_g-1))/(mu*(mult-1))
                dkjdS =  ((j-1)*(phi_g-1) * (1-mult**(-1)) - 1)/(mult-1.) * kj
           
                
            prod = prod * kj            
            series = series + dkjdS/kj
            
        dX_M_dS = dX_M_dS + prod * (series+ds_m_dS) * s_m 
        dX_N_dS = dX_N_dS + prod * series
        dX_G_dS = dX_G_dS + prod * (series+ds_g_dS)* s_g
        dX_nu_dS = dX_nu_dS + prod * (series+ds_nu_dS)* s_nu
        
        
        
    return dX_M_dS, dX_N_dS, dX_G_dS, dX_nu_dS
        
    
    
def dX2_dS2(J,mult,mu,phi_g = 0.0, phi_a = 0.0):
    
    """
    Function dX2_dS2 is used to determine the optimal scaling to fit the 
    continous solutions. It is the double differntial of X with respect to the 
    bin scaling.
    
    Inputs:
        
        - J, the number of mass classes. (-)
        
        - mult, multiplicative bin scaling parameter. (-)
        
        - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
        
        - phi_g, the growth mass scaling parameter. (-)
        
        - phi_a, the area mass scaling parameter. (-)
        
        
    Outputs:

        - dX_M2_dS2, summation of the mass scaling bin scaling differential 
                   over each mass class. (-)        

        - dX_N2_dS2, summation of the scaling of the number density bin scaling
                   differential over each mass class. (-)   
         
        - dX_G2_dS2, summation of the scaling of the growth bin scaling 
                   differential over each mass class. (-)   
        
    
        - dX_nu2_dS2, summation of the vegetation fraction bin scaling 
                    differential over each mass class. (-)  
        
    """
    
    dX_M2_dS2 = 0.
    
    dX_N2_dS2 = 0.
    
    dX_G2_dS2 = 0.
    
    dX_nu2_dS2 = 0.
    
    
    for i in range(0,J):
        
        prod = 1
        series_M = 0
        series = 0       
        series_g = 0
        series_nu = 0
        
        s_m = mult ** (i)
        ds_m_dS = i * mult ** (i-1)
        s_g = mult ** (i*phi_g)
        ds_g2_dS2 = i*phi_g*mult**(-2)*(i*phi_g-1)
        s_nu = mult ** (i*phi_a)
        ds_nu_dS = i*phi_a *mult**(i*phi_a - 1)

        for j in range(1,i + 1):
            
            if j < J - 1:
                
                kj = k_j(j,mult,mu,phi_g)
                dkjdS = dk_j_dS(j,mult,mu,phi_g)
                dkj2dS2 = dk_j2_dS2(j,mult,mu,phi_g)
            
            elif j == J - 1 :
                
                kj =  mult**((j-1)*(phi_g-1))/(mu*(mult-1))
                dkjdS =  ((j-1)*(phi_g-1) * (1-mult**(-1)) - 1)/(mult-1.) * kj
                dkj2dS2 = dkjdS**2/kj + (kj/(mult-1.))*((j-1)*(phi_g-1.)\
                                         /(mult**2)-dkjdS/kj)

            prod = prod * kj            
            series = series + dkj2dS2/kj
            series_g = series_g + (2*j*phi_g*mult**(-1)*dkjdS - dkj2dS2)/kj 
            
        dX_M2_dS2 = dX_M2_dS2 + prod * series * s_m + ds_m_dS *prod
        dX_N2_dS2 = dX_N2_dS2 + prod * series
        dX_G2_dS2 = dX_G2_dS2 + s_g * prod *(ds_g2_dS2 + series_g) 
        dX_nu2_dS2 = dX_nu2_dS2 + prod * series * s_nu + ds_nu_dS *prod
        
        
        
    return dX_M2_dS2, dX_N2_dS2, dX_G2_dS2, dX_nu2_dS2
        
    
            
        
def dk_j_dS(j,mult,mu,phi_g):
    """
    Outputs the scaling charateristics gradient of each mass class against S.
    
    Inputs: 
    
    - j, the given mass class. (-)
    
    - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
    
    - phi_g, the growth mass scaling parameter. (-)
    
    Outputs:
    
    - dk_j_dS, the relative change of fraction population in the current
      class against S. (-)
    
    """
    dk_j_dS = k_j(j,mult,mu,phi_g) * ((j-1.) * (phi_g-1.)/mult - \
                  k_j(j,mult,mu,phi_g)*(j*(phi_g-1)*mult**(phi_g-2)+\
                      mu*mult**((j-1)*(1-phi_g))))

    return dk_j_dS
    
    
def dk_j2_dS2(j,mult,mu,phi_g):
    """
    The double derivative of the scaling charateristic function with respect to
    the bin scaling.
        Inputs: 
    
    - j, the given mass class. (-)
    
    - mu, the turnover ratio defined as: gamma_init * m_init/g_init (-)
    
    - phi_g, the growth mass scaling parameter. (-)
    
    Outputs:
    
    - dk_j2_dS2, the relative change of fraction population in the current
      class against the double derivative of S. (-)
    
    """
    kj = k_j(j,mult,mu,phi_g)
    dkjdS = dk_j_dS(j,mult,mu,phi_g)
    r = dkjdS/kj
    dk_j2_dS2 = kj * (-r*(j-1.)*(phi_g-1.)*mult**(-1)-(j-1.)*(phi_g\
                      -1.)*mult**(-2.)-kj*(phi_g-1.)*mult**(-1.)*(j*(phi_g-1.)\
                      *mult**(phi_g-2.)-mu*(j-1.)*mult**((j-1.)*(1.-phi_g))))
    
    return dk_j2_dS2