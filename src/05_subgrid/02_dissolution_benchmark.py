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

def opt_fun(x, ceq, cm, c0):
    a = x*exp(x**2)*erfc(-x)
    b = (ceq-c0)/(cm-ceq)
    return (np.sqrt(np.pi)*a-b)
#%%
cm = 1#0.95
ceq = 0.4#0.01949
c0 = 0.1# 0
crange = np.arange(-2,1,0.01)
f = []
for v in crange:
    f.append(opt_fun(v, ceq, cm, c0))
plt.plot(crange, np.array(f))
plt.show()

res = root(opt_fun, 0, args = (ceq, cm, c0))
#res = minimize(lambda x: opt_fun(x, ceq, cm, c0), x0=-0.25)
#print(res.success)
#print(res.x[0])

k = res.x[0]
print(k)
time = np.arange(0, 8 ,0.1)
plt.plot(time, np.abs(2*k*np.sqrt(1e-3*time*3600)))
plt.show()

#%%

cm = 0.95
ceq = 0.01949
c0 = 0
crange = np.arange(-2,1,0.01)
f = []
for v in crange:
    f.append(opt_fun(v, ceq, cm, c0))
plt.plot(crange, np.array(f))
plt.show()

res = root(opt_fun, 0, args = (ceq, cm, c0))
#res = minimize(lambda x: opt_fun(x, ceq, cm, c0), x0=-0.25)
#print(res.success)
#print(res.x[0])

k = res.x[0]
print(k)
time = np.arange(0, 10 ,0.01)
plt.plot(time, np.abs(2*k*np.sqrt(1e+3*time)))
plt.show()