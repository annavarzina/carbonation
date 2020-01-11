# -*- coding: utf-8 -*-
"""
...
"""

#%% PYTHON MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=3, threshold=np.inf)
scale = 100
#import phrqc
fractions = np.array([1., 0.7, 0.5, 0.3, 0.2, 0.1,
                      0.07, 0.05, 0.03,0.02, 0.01,
                      0.007, 0.005, 0.003, 0.002, 0.001 ])
time = {}
ch = {}
ca ={}
dch = {}
dissolution = {}
nn="03_subgrid_leaching_depth"#os.path.basename(__file__)[:-3] 
for f in fractions:
    path = root_dir+'\\results\\output\\10_subgrid_leaching\\' + nn + str(f)+'\\'
    dissolution[nn+ str(f)] = np.load(path +'dis_time'+'.npy')
    ch[nn+ str(f)] = np.load(path +'CH'+'.npy')
    ca[nn+ str(f)] = np.load(path +'Ca_prof'+'.npy')
    dch[nn+ str(f)] = np.load(path +'dCH'+'.npy')
    time[nn+ str(f)] = np.load(path +'time'+'.npy')
#%%
plt.figure()
for f in fractions:    
    t = np.array(dissolution[nn + str(f)])*scale
    x = np.arange(1,len(dissolution[nn + str(f)])+1)
    plt.plot(t, x, label = f)
plt.ylabel("Dissolved length (um)")
plt.xlabel("Time")
plt.legend()
plt.show()

plt.figure()
for f in fractions:  
    t = np.array(dissolution[nn + str(f)])*scale
    x = np.arange(1,len(dissolution[nn + str(f)])+1)*1e-6
    D = (x**2)/2/t
    plt.plot(t, D, label = f)
plt.ylabel("D (m2/s)")
plt.xlabel("Time (s)")
plt.yscale("log")
plt.legend()
plt.show()

d1 = []
dn = []
for f in fractions:
    t = np.array(dissolution[nn + str(f)])*scale
    x = np.arange(1,len(dissolution[nn + str(f)])+1)*1e-6
    D = (x**2)/2/t
    d1.append(D[0])
    dn.append(D[-1])
    
plt.figure()
plt.plot( d1, fractions, label = "1")
#plt.plot(fractions, dn, label = "20")
plt.xlabel("D (m2/s)")
plt.ylabel("Fraction")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()