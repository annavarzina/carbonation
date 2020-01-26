# -*- coding: utf-8 -*-
'''
The default case
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
import func as cf
#%% SETTINGS

fname = 'dc'
fpath = root_dir+'\\results\\output\\00_default\\'
fn.make_output_dir(fpath)
#names = np.array([ '01_ie01_p05', '02_ie05_p05', '03_ie1_p05'])
name = '02_c05_p005_Db09'

linetype = '-'
scale = 100

results = {}
path = root_dir+'\\results\\output\\09_crystal_size\\' + name + '\\'
results = fn.load_obj(path + name +'_results')

for n in ['time', 'portlandite', 'calcite']:
    temp = np.array(results[n])
    temp *= scale 
    results[n]= temp.tolist()

concCA = np.load(path +'Ca.npy')
concC = np.load(path +'C.npy')
ph = np.load(path +'pH.npy')
de = np.load(path +'De.npy')
poros = np.load(path +'poros.npy')
si = np.load(path +'SI.npy')

#%% MAIN PROPERTIES
r1 = 1000
r2 = 10000
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
          'Average pH', 'Input C', 'Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'pH', 
         'C (1, 0)', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_average ph', '_input_c', '_poros']
ylabel = [r'Portlandite $\cdot 10^{-15}$ (mol)', r'Calcite  $\cdot 10^{-15}$ mol', 
          r'Dissolved Ca $\cdot 10^{-15}$ (mol)', r'Dissolved C $\cdot 10^{-15}$ (mol)' ,
          'Average pH', 'Boundary C (mol/l)', 'Porosity']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    #for i in range(0, len(names)):
    plt.plot(np.array(results['time'][:r2])/3600, results[comp[k]][:r2],
             ls=linetype)
    plt.ylabel(ylabel[k])
    plt.xlabel('Time (h)')
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 

#%% pH profile
nch = results['portlandite_cells'][-1]
ncc = results['calcite_cells'][-1]
plt.figure(figsize=(8,4))
plt.plot(range(0,len(ph[1,:])-1),ph[1,0:-1], color = "#f2630d", label = 'pH')
plt.fill_between(np.arange(1,7,1),ph[1,1:7], color = "#6f8191", alpha=.5,label = "calcite")
plt.fill_between(np.arange(6,11,1),ph[1,6:11], color = "#cbd2d9", alpha=.5, label = "depleted portlandite")
plt.fill_between(np.arange(10,30,1),ph[1,10:30], color = "#3d4249", alpha=.5, label = "portlandite")
plt.ylabel('pH')
plt.legend()
plt.ylim(8.1,12.6)
plt.xlabel(r'Distance ($\mu m$)')

#%% porosity profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(ph[1,:])-1),poros[1,0:-1])
plt.ylabel('Porosity')
plt.yscale("log")
plt.xlabel(r'Distance ($\mu m$)')

#%% Ca profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(concCA[1,:])-1),concCA[1,0:-1], label = 'Ca')
plt.fill_between(np.arange(1,7,1),concCA[1,1:7], color = "#6f8191", alpha=.5,label = "calcite")
plt.fill_between(np.arange(6,11,1),concCA[1,6:11], color = "#cbd2d9", alpha=.5, label = "depleted portlandite")
plt.fill_between(np.arange(10,30,1),concCA[1,10:30], color = "#3d4249", alpha=.5, label = "portlandite")
plt.ylabel('Dissolved Ca (mol/l)')
plt.xlabel(r'Distance ($\mu m$)')
plt.legend()
#%% C profile 
plt.figure(figsize=(8,4))
plt.plot(range(1,len(concC[1,:])-1),concC[1,1:-1])
plt.ylabel('C')
plt.xlabel(r'Distance ($\mu m$)')

#%% SI profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(si[1,:])-1),si[1,0:-1])
plt.ylabel('SI')
plt.xlabel(r'Distance ($\mu m$)')

#%% CC profile
cc_mass = 100.09
ch_mass = 74.09
cc = []
ch = []
for i in range(1,10):
    cc.append(results['calcite (1, ' + str(i) +')' ][-1]*cc_mass)
    ch.append(results['portlandite (1, ' + str(i) +')' ][-1]*ch_mass)
chmax = results['portlandite (1, ' + str(i) +')' ][1]*ch_mass
for i in range(10,30):
    cc.append(0)
    ch.append(chmax)
    
plt.figure(figsize=(8,4))
plt.plot(range(1,30),cc[:])
plt.plot(range(1,30),ch[:])
plt.fill_between(np.arange(1,7,1),cc[0:6], color = "#6f8191", alpha=.5,label = "calcite")
plt.fill_between(np.arange(8,30,1),ch[7:30], color = "#3d4249", alpha=.5, label = "portlandite")
plt.ylabel(r'Mineral mass  $\cdot 10^{-15}$ (g) in 1 $\mu m^3$ ')
plt.xlabel(r'Distance ($\mu m$)')
plt.legend()

#%% Points
r1 = 1
r3 = 10000
p = 7
f = "Ca"
comp = ['%s (1, %s)'%(f,i) for i in range(1,p)]
title = ['%s in %s'%(f,i) for i in range(1,p)]

plt.figure(figsize=(8,4))    
for k in range(0, len(comp)):
    plt.plot(np.array(results['time'][r1:r2])/3600, results[comp[k]][r1:r2],
             ls=linetype, label = title[k])
plt.xlabel('Time (h)')
plt.ylabel('Dissolved Ca (mol/l)')
plt.legend(loc = "upper right")
plt.show() 