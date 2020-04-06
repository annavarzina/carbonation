# -*- coding: utf-8 -*-
'''
Compare the results for different water layer
'''
#%% MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import misc_func as fn
#%% SETTINGS
Ts =200.
scale=50.
fname = 'subrid'
fpath = root_dir+'\\results\\temp\\'
fn.make_output_dir(fpath)

#names = np.array(['01_ll1_p005_subgrid', '01_ll3_p005_subgrid', '01_ll5_p005_subgrid'])
names = np.array(['03_subgrid_leaching_depth_05', '04_subgrid_leaching_depth_005', 
                  '05_subgrid_leaching_depth_0005'])
label = np.array(['0.5', '0.05','0.005'])
linetype = np.array(['-', '--', '-.'])
#%%  LOAD
    
time = {}
ch = {}
ca ={}
dch = {}
dissolution = {}
for nn in names:
    path = root_dir+'\\results\\temp\\' + nn + '\\'
    dissolution[nn] = np.load(path +'dis_time'+'.npy')
    ch[nn] = np.load(path +'CH'+'.npy')
    ca[nn] = np.load(path +'Ca_prof'+'.npy')
    dch[nn] = np.load(path +'dCH'+'.npy')
    time[nn] = np.load(path +'time'+'.npy')
#%%DISSOLUTION

plt.figure()
for i in range(0,len(names)):    
    t = np.array(dissolution[names[i]])*scale
    x = np.arange(1,len(dissolution[names[i]])+1)
    plt.plot(t, x, label = label[i])
plt.ylabel("Dissolved length (um)")
plt.xlabel("Time")
plt.legend()
plt.show()
#%%

plt.figure()
for i in range(0,len(names)):
    t = np.array(dissolution[names[i]])*scale
    x = np.arange(1,len(dissolution[names[i]])+1)*1e-6
    D = x**2/2/t
    plt.plot(t, D, label = label[i])
plt.ylabel("D (m2/s)")
plt.xlabel("Time (s)")
plt.legend()
plt.show()

#%% CH    
plt.figure()
for i in range(0,len(names)):
    plt.plot(np.array(time[names[i]]), ch[names[i]], label = label[i])
plt.ylabel("CH (mol/l)")
plt.xlabel("Time")
plt.legend()
plt.show()

#%% CH rate    
plt.figure()
for i in range(0,len(names)):
    plt.plot(np.array(time[names[i]]), np.abs(dch[names[i]]), label = label[i])
plt.ylabel("CH rate (*10^-15 mol/s)")
plt.xlabel("Time")
plt.ylim(0,10)
plt.legend()
plt.show()
#%% Ca    
plt.figure()
for i in range(0,len(names)):
    plt.plot( ca[names[i]], label = label[i])
plt.ylabel("CH (mol/l)")
plt.xlabel("X (um)")
plt.legend()
plt.show()