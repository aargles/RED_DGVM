"""

Python 2.7, Encoding: UTF-8

Created on Mon Jan 24 2019

@author: A.Argles (aa760@exeter.ac.uk)
    
Description: Methods used in the Dynamical Subroutine Libary.

"""

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
    
    pivot = A_flat[(start+end)/2]
              
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
    
    
    