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

#%% rate
scale = 100
dim = 10**(-15)*10**3*10**8 # convert to mmol/l/s/cm2
plt.figure(figsize = (6,6))
for f in fractions:
    plt.plot(np.abs(dch[nn+ str(f)][1:])*scale*dim, label = f)
plt.ylabel(r"Rate of dissolution $k$ $(mmol / (l\cdot s\cdot cm^2)$")
plt.xlabel(r"Dissolved portlandite length $(\mu)$")
plt.yscale("log")
#plt.xscale("log")
plt.legend()
plt.show()

#%% CH rate
r1 = []
rn = []
rmean = []
for f in fractions:    
    r1.append(dch[nn + str(f)][1]*scale*dim)
    rn.append(dch[nn + str(f)][-11]*scale*dim)
    rmean.append(np.mean(dch[nn + str(f)][:])*scale*dim)
    
    
plt.figure()
plt.plot(fractions, np.abs(r1), label = "1")
plt.plot(fractions, np.abs(rn), label = "20")
plt.plot(fractions, np.abs(rmean), label = "mean")
plt.ylabel(r"Rate of dissolution $k$ $(mmol / (l\cdot s\cdot cm^2)$")
plt.xlabel(r"Fraction of equilibrated water $\sigma$")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.show()

#%%
r = rn#rmean
idx = np.where(np.logical_and(np.abs(r)>=0.39e-5, np.abs(r)<=6.2e-5))[0]
plt.figure(figsize = (8,4), dpi = 500)
plt.plot(np.abs(r), fractions)
plt.fill_between(np.abs(r)[idx],fractions[idx], color = "#6f8191", alpha=.5,label = "Johannsen K., et.al. (1999)")
plt.xlabel(r"Rate of dissolution $R$ $(mmol / (l\cdot s\cdot cm^2)$",fontsize=14)
plt.ylabel(r"Mixing parameter $\sigma$ ", fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.yscale("log")
plt.xscale("log")
plt.legend(fontsize=12)
plt.show()
#%% ch
scale = 100
chsc = np.array(ch[nn+ str(fractions[0])]*scale)
print(chsc)
