# -*- coding: utf-8 -*-
"""
Default and constant values in the system
"""
import numpy as np

#%% PCS
ANGLE = np.pi
MVOL_CC = 3.69e-5 #calcite molar volume (m3/mol)
GAS_CONSTANT = 8.314 # gas constant
TEMPERATURE = 298.3 # temperature in kelvin
N_IONS = 1 # number of ions per unit salt

INTERFACIAL_ENERGY = 0.485
PORE_SIZE = 2E-7
PORE_DENSITY = 2E+5
CRYSTAL_SIZE = 5E-7
                   
#%%                
MOLAR_VOLUME = {'CH': 33.10e-3, # l/mol
                'CC': 36.90e-3,
                'CSH_TobH': 55.30e-3,
                'CSH_TobD': 47.95e-3,
                'CSH_JenH': 75.63e-3,
                'CSH_JenD': 80.58e-3,
                'CaO':21.0e-3,
                'SiO2':27.5e-3}