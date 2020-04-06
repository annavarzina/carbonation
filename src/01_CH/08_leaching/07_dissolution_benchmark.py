# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 14:31:38 2019

@author: mat
"""

import numpy as np

from scipy.optimize import root
#from scipy.special import erfc
import matplotlib.pylab as plt
from math import exp
from math import erfc

def equation(x, ceq, cm, c0):
    a = x*exp(x**2)*erfc(-x)
    b = (ceq-c0)/(cm-ceq)
    return (np.sqrt(np.pi)*a-b)
#%% Benchmark from thesis of Ravi
cm = 1 #mol/l
ceq = 0.4 #mol/l
c0 = 0.1 #mol/l
k_range = np.arange(-2,1,0.01)
f = []
'''
for k in k_range:
    f.append(equation(k, ceq, cm, c0))
plt.plot(k_range, np.array(f))
plt.show()
'''
res = root(equation, 0, args = (ceq, cm, c0))
k = res.x[0]
print(k)

time = np.arange(0, 8 ,0.1) #hours
time_s = time*3600 #seconds
D = 1e-3 # mm2/s = 1e-9 m2/s
plt.figure()
plt.plot(time, np.abs(2*k*np.sqrt(D*time_s)))
plt.xlabel('Time (h)')
plt.ylabel('dr (mm)')
plt.show()

#%%

cm = 0.199999 #mol/l
ceq = 0.01949 #mol/l
c0 = 0 #mol/l
res = root(equation, 0, args = (ceq, cm, c0))
k = res.x[0]
print(k)
time = np.arange(0, 1.0 ,0.001) #minutes
time_s = time#*60#seconds
D = 1e-9*1e+12 #um2/s
#D = 1e-9*1e+6 # mm2/s 
plt.figure()
plt.plot(time, np.abs(2*k*np.sqrt(D*time_s)))
plt.xlabel('Time (min)')
plt.ylabel('dr (um)')
plt.show()


plt.figure()
plt.plot(time_s, np.abs(2*k*np.sqrt(D*time_s)))
plt.xlabel('Time (sec)')
plt.ylabel('dr (um)')
plt.show()
#%%
print(np.abs(2*k*np.sqrt(D*0.09074999999999936)))
print(np.abs(2*k*np.sqrt(D*0.3688500000000422)))
print(np.abs(2*k*np.sqrt(D*0.8334000000002864)))
dissolution_time = np.array([0.,0.09074999999999936, 0.3688500000000422, 0.8334000000002864])
results = np.array([0, 1.,2.,3.])
plt.figure()
plt.plot(time, np.abs(2*k*np.sqrt(D*time_s)))
plt.plot(dissolution_time, results)
plt.xlabel('Time (min)')
plt.ylabel('dr (um)')
plt.show()