"""
Initialise top level RED 
Python 2.7, Encoding: UTF-8
Imports RED subroutines into namespace

Created on Mon Jul  2

@author: A.Argles (aa760@exeter.ac.uk)
"""

from .EXCEPTION import INITILISATION_CHECK, DYNAMICAL_CHECK
from .STORE import SAVE_PFT_VALUES
from .MISC_METHODS import LOAD_PFT_VALUES
from .DYNAMICAL import DYN_get_outputs
from .EQUILIBRIUM import EQ_get_outputs,EQ_mu,EQ_P_gamma,EQ_sample,get_gamma
                         