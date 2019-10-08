# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:11:51 2019

@author: avarzina
"""

import numpy as np 
import matplotlib.pylab as plt

Dw = 1e-9
Dch = np.array([1e-11,1e-12,1e-13])

Mw = np.arange(0,1.01,0.01)
Mch = 1 - Mw

def deff(dw, dch, mw, mch):
    return 1/(mw/dw+mch/dch)
De = deff(Dw, Dch[2], Mw, Mch)
#%%
for i in range(0,len(Dch)):
    plt.plot(Mch, deff(Dw, Dch[i], Mw, Mch), label = str(Dch[i]))
plt.legend()
plt.yscale('log')
plt.ylabel('Deff')
plt.xlabel('Mass fraction of CH')
plt.show()
#%%

def deff2(dw, dch, mw, mch):
    return mw*dw+mch*dch

for i in range(0,len(Dch)):
    plt.plot(Mch, deff2(Dw, Dch[i], Mw, Mch), label = str(Dch[i]))
plt.legend()
plt.yscale('log')
plt.ylabel('Deff')
plt.xlabel('Mass fraction of CH')
plt.show()