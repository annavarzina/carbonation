# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 22:06:22 2019

@author: avarzina
"""

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
#fraction = np.array([0.001, 0.002, 0.003, 0.005, 0.01, 0.05,  0.1, 0.2, 0.3, 0.5, 1 ])
#diffusivity = np.array([9.89E-10, 9.68E-10, 9.45E-10, 9.16E-10, 8.42E-10, 7.30E-10, 3.06E-10, 1.47E-10, 7.54E-11, 4.60E-11, 1.78E-11])
#%% LOAD DATA
x = pd.read_csv('subgrid.csv')
fraction = x['Fraction'].values
diffusivity = x['D'].values

fr = np.arange(0.0005,1,0.0001)
#%% PREDICTION 1
#'''
def diff_predict(f, c0, c1,c2,c3):
    d = c0 - c1* np.exp(-f*c2) + c3*f#**(1./2.)
    return d

c, cov = curve_fit(diff_predict, fraction, diffusivity)
print(c)

diff_opt = diff_predict(fraction, c[0],c[1],c[2],c[3])
diff_opt2=diff_predict(fr, c[0],c[1],c[2],c[3])
print('R2: ', r2_score(diff_opt, diffusivity))
#'''
'''
def diff_predict(f, c0, c1,c2,c3):
    #d = dw*(1 - np.exp(-f*c))
    d = c0 + c1*np.log(f*c2)
    return d

c, cov = curve_fit(diff_predict, fraction, diffusivity)

diff_opt = diff_predict(fraction, c[0], c[1], c[2], c[3])
'''
'''
def diff_predict(f, c0, c1):
    d = c0*(1 - np.exp(-f*c1))
    #d = c0*(1-np.exp(-f*c1))
    return d

c, cov = curve_fit(diff_predict, fraction, diffusivity)

diff_opt = diff_predict(fraction, c[0],c[1])
diff_opt2=diff_predict(fr, c[0],c[1])
'''

plt.loglog(fraction, diffusivity, 'r.')
plt.loglog(fr, diff_opt2, 'b-')
plt.xlabel('Fraction')
plt.ylabel('D')
plt.show()

plt.plot(fraction, diffusivity/1e-9, 'r.')
plt.plot(fr, diff_opt2/1e-9, 'b-')
plt.xlabel('Fraction')
plt.ylabel('D')
plt.show()

#%% PREDICTION 2

def diff_predict(f, c0, c1):
    d = c0*(1 - np.exp(-f*c1))
    #d = c0*(1-np.exp(-f*c1))
    return d

c, cov = curve_fit(diff_predict, fraction, diffusivity)
print(c)

diff_opt = diff_predict(fraction, c[0],c[1])
diff_opt2=diff_predict(fr, c[0],c[1])
print('R2: ', r2_score(diff_opt, diffusivity))


plt.loglog(fraction, diffusivity, 'r.', label = "Result")
plt.loglog(fr, diff_opt2, 'b-', label = "Fitting")
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.legend()
plt.show()

plt.plot(fraction, diffusivity, 'r.', label = "Result")
plt.plot(fr, diff_opt2, 'b-', label = "Fitting")
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.yscale("log")
plt.legend()
plt.show()
#%% PREDICTION 2

def diff_predict(f, c0, c1, c2):
    d = c0 + c1*np.exp(-f*c2)
    #d = c0*(1-np.exp(-f*c1))
    return d

c, cov = curve_fit(diff_predict, fraction, diffusivity)
print(c)

diff_opt = diff_predict(fraction, c[0],c[1], c[2])
diff_opt2=diff_predict(fr, c[0],c[1], c[2])
print('R2: ', r2_score(diff_opt, diffusivity))


plt.loglog(fraction, diffusivity, 'r.', label = "Result")
plt.loglog(fr, diff_opt2, 'b-', label = "Fitting")
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.legend()
plt.show()

plt.plot(fraction, diffusivity, 'r.', label = "Result")
plt.plot(fr, diff_opt2, 'b-', label = "Fitting")
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.yscale("log")
plt.legend()
plt.show()
#%% PREDICTION 3
def diff_predict(f, c0,c1,c2,c3):
    d = c0*(1 - c1*np.log((f+c3)*c2))#**(1./2.)
    return d
D = 1e-9
c, cov = curve_fit(diff_predict, fraction, diffusivity/D)
print(c)

diff_opt = diff_predict(fraction, c[0], c[1], c[2], c[3])
diff_opt2=diff_predict(fr, c[0], c[1], c[2], c[3])
print('R2: ', r2_score(diff_opt, diffusivity/D))

plt.loglog(fraction, diffusivity/D, 'r.')
plt.loglog(fr, diff_opt2, 'b-')
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.show()

plt.plot(fraction, diffusivity/D, 'r.')
plt.plot(fr, diff_opt2, 'b-')
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.show()
#%% PREDICTION 3
def diff_predict(f, c0,c1,c2):
    d = c0-c1/f**c2#**(1./2.)
    return d

test  = diff_predict(fr, 1.0e-9, 1e-12, 0.9)
plt.plot(fraction, diffusivity, 'r.')
plt.plot(fr, test, 'b-')
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.show()

plt.loglog(fraction, diffusivity, 'r.')
plt.loglog(fr, test, 'b-')
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.show()

#%%
D = 1e-9
c, cov = curve_fit(diff_predict, fraction, diffusivity)
print(c)

diff_opt = diff_predict(fraction, c[0], c[1],c[2])
diff_opt2=diff_predict(fr, c[0], c[1],c[2])
print('R2: ', r2_score(diff_opt, diffusivity))

plt.loglog(fraction, diffusivity, 'r.')
plt.loglog(fr, diff_opt2, 'b-')
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.show()

plt.plot(fraction, diffusivity, 'r.')
plt.plot(fr, diff_opt2, 'b-')
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.show()
#%%
plt.plot(fraction, diffusivity, 'r-')
plt.xlabel('Mixing fraction')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.yscale("log")
plt.show()