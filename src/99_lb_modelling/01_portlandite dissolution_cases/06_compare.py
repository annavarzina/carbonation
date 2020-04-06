# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:57:33 2019

@author: avarzina
"""

from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
import numpy as np
import pickle
np.set_printoptions(precision=5, threshold=np.inf)

namesCa = ['Ca_D1e-09',
         'Ca_D1e-10',
         'Ca_D1e-11',
         'Ca_D1e-12']
resultsCa = {}
for i in range(0,len(namesCa)):
    with open(namesCa[i] + '.pkl', 'rb') as f:
        resultsCa.update({str(i): pickle.load(f)})
        
        
namesCH = ['CH_D1e-09',
         'CH_D1e-10',
         'CH_D1e-11',
         'CH_D1e-12']
resultsCH = {}
for i in range(0,len(namesCH)):
    with open(namesCH[i] + '.pkl', 'rb') as f:
        resultsCH.update({str(i): pickle.load(f)})
        
#%%        
plt.figure()
for i in range(0,len(namesCa)):
    plt.plot(resultsCa[str(i)], label=namesCa[i])
plt.xlabel('Iterations')
plt.ylabel('Average Ca  [mM]')
plt.title('Average Ca concentration in aqeuous phase')
plt.legend()
plt.show()
    
plt.figure()
for i in range(0,len(namesCH)):
    plt.plot(resultsCH[str(i)], label=namesCH[i])
plt.xlabel('Iterations')
plt.ylabel('Total CH [mM]')
plt.title('Total portlandite mass')
plt.legend()
plt.show()