# -*- coding: utf-8 -*-
'''
Compare the results for different fraction
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
Ts =3600.
fname = 'fraction'
fpath = root_dir+'\\results\\output\\03_fraction\\'
fn.make_output_dir(fpath)

#names = np.array(['01_ll1_p005_subgrid', '01_ll3_p005_subgrid', '01_ll5_p005_subgrid'])
#names = np.array(['01_example_1D_constD_fr04', '02_example_1D_constD_mlvl']) #,'02_example_1D_constD_fr1', 
#names = np.array(['01_example_1D_constD_fr04','05_neq_constD_fr04'])
#label = np.array(['eq', 'neq'])
names = np.array(['05_neq_constD_fr04', '06_neq_constD_mlvl'])
scale = 100
names = np.array(['06_neq_Di_fr04', '06_neq_Di_mlvl'])
scale = 50
label = np.array([r'$\sigma=0.004$', r'$\sigma=1.0$'])
p = '\\results\\output\\03_fraction\\'
linetype = np.array(['-', '--', '-.', ':'])


results = {}
for nn in names:
    path = root_dir + p+ nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

for i in range(0, len(names)):
    for n in ['time', 'portlandite', 'calcite']:
        temp = np.array(results[names[i]][n])
        temp *= scale
        results[names[i]][n]= temp.tolist()
#%% PARAMS
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
          'Average pH', 'Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_average ph', 'avg_poros']
ylabel = [r'Portlandite $\cdot 10^{-15}$ (mol)', r'Calcite  $\cdot 10^{-15}$ mol', 
          r'Dissolved Ca $\cdot 10^{-15}$ (mol)', r'Dissolved C $\cdot 10^{-15}$ (mol)' ,
          'Average pH', 'Boundary C (mol/l)', 'Porosity']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(np.array(results[names[i]]['time'])/3600,
                 results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.ylabel(ylabel[k], fontsize = 14)
    plt.xlabel('Time (h)', fontsize = 14)
    plt.legend(fontsize = 12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 

#%% DISSOLUTION/PRECIPITATION RATE

titles = ['Dissolution rate', 'Precipitation rate' ]
comp =  ['portlandite', 'calcite']
suffix = ['_CH_rate', '_CC_rate' ]
s =2
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        rate = np.abs(cf.get_rate(results[names[i]][comp[k]],
                             results[names[i]]['time'][2] - results[names[i]]['time'][1],
                             step = s ))
        #print(len(rate))
        plt.plot(np.array(results[names[i]]['time'][::s])/3600, 
                 rate,
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (h)')
    plt.ylabel('Rate (mol/s)')
    plt.yscale("log")
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show()
#plt.savefig(fpath + fname + '_CH_rate')
#%% POINTS
f = ('C', 'C') # ('Ca', 'Ca') ('Volume CH', 'vol_CH') ('Volume CC', 'vol_CC') ('pH', 'pH') ('De', 'De')('Porosity', 'poros')
title = ['%s in %s'%(f[0],i) for i in range(1,9)]
comp = ['%s (1, %s)'%(f[1],i) for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
    plt.title(label[i])
    plt.legend()
    plt.show() 

#%% C profile
c_conc = np.load(path + 'C.npy')
plt.figure(figsize = (8,4))
plt.plot(c_conc[1,0:-2])
plt.xlabel(r'Length ($\mu$m)')
plt.ylabel('pH')
#%% pH profile
pH = np.load(path + 'pH.npy')
plt.figure(figsize = (8,4))
plt.plot(pH[1,0:-2])
plt.xlabel(r'Length ($\mu$m)')
plt.ylabel('pH')

