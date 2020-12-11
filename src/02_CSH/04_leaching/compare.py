# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 01:24:16 2020

@author: avarzina
"""

import numpy as np
import matplotlib.pylab as plt

fraction = np.array([1., 0.1, 0.01, 0.001])
casi1 = np.array([-9.56, -9.67, -10.2, -11.09])
casi2 = np.array([-9.74, -9.86, -10.38, -11.27])
casi3 = np.array([-9.87, -9.99, -10.52, -11.42,])
#%%
plt.figure()
plt.plot(fraction, casi1, '-*', label = "Ca/Si 1.67")
plt.plot(fraction, casi2, '-o', label = "Ca/Si 1.25")
plt.plot(fraction, casi3, '-x', label = "Ca/Si 0.83")
plt.plot(np.array([0.5,0.5,0.5]), np.array([-8.40, -9.98, -10.99]), '+', label = "T.-B.")
plt.legend()
plt.xlabel(r'Mixing fraction $\sigma$')
plt.ylabel(r'$log K$ ($mol/s/m^2$)')
plt.show()