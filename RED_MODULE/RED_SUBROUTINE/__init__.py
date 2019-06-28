"""
Initialise top level RED 
Python 2.7, Encoding: UTF-8
Imports RED subroutines into namespace

Created on Mon Jul  2

@author: A.Argles (aa760@exeter.ac.uk)
"""

from .EXCEPTION import INITILISATION_CHECK, DYNAMICAL_CHECK
from .STORE import PFT_VALUES
from .DYNAMICAL import CONVERSION, ALLOMETRY, COMPETITION, DEMOGRAPHY
from .EQUILIBRIUM import PFT_HIERARCHY, INTILISATION_nu,EQUIL_FRAC,\
                         INTILISATION_P_gamma_init,INTILISATION_OUTPUT
                         
                         