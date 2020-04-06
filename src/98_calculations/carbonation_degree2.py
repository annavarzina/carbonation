# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
#%% DATA

alpha = np.array([0.,0.3,0.55,1.05,1.25,1.55])
time = np.array([0.,2.,4.5,10.,12.,20.]) #days
time_new = np.arange(0.,100.,1.)

mass0 = 22.977 *1e-3 #g
dmass = alpha*mass0/100 
dmolmass = 25.993 #26.812 g/mol
xmol = dmass/dmolmass #constant = rho_CC/C_CC - rho_CH/C_CH

sa  = 29.19e+6 #m2
d = xmol*0.0331*1e+15/sa#m

msa = dmass/sa #delta mass per 1 um2
molsa = xmol/sa