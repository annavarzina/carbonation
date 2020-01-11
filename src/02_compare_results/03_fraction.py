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
names = np.array(['05_neq_constD_fr04', '06_neq_constD_mlvl'])
label = np.array([r'$\sigma=0.004$', r'$\sigma=1.0$'])
#names = np.array(['01_example_1D_constD_fr04','05_neq_constD_fr04'])
#label = np.array(['eq', 'neq'])
linetype = np.array(['-', '--', '-.', ':'])
results = {}
for nn in names:
    path = root_dir+'\\results\\output\\03_fraction\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')
    #%%
scale = 100
for i in range(0, len(names)):
    temp = np.array(results[names[i]]['time'])
    temp *= scale
    results[names[i]]['time']= temp.tolist()
#%% CH DISSOLUTION 
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
          'Average pH', 'Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_average ph', 'avg_poros']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 

#%% DISSOLUTION RATE

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
        plt.plot(results[names[i]]['time'][::s], 
                 rate,
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.ylabel('Rate (mol/s)')
    plt.yscale("log")
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show()
#plt.savefig(fpath + fname + '_CH_rate')
#%%
    
titles = [ 'Volume CH in (1,2)', 'Volume CH in (1,3)','Volume CH in (1,4)',
          'Volume CH in (1,5)', 'Volume CH in (1,6)','Volume CH in (1,7)',]
comp =  [ 'vol_CH (1, 2)', 'vol_CH (1, 3)',  'vol_CH (1, 4)','vol_CH (1, 5)', 
         'vol_CH (1, 6)',  'vol_CH (1, 7)']

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = titles[k])
    plt.title(label[i])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.show() 
#%%
titles = [ 'Volume CC in (1,1)','Volume CC in (1,2)', 'Volume CC in (1,3)','Volume CC in (1,4)',
          'Volume CC in (1,5)', 'Volume CC in (1,6)','Volume CC in (1,7)']
comp =  [ 'vol_CC (1, 1)','vol_CC (1, 2)', 'vol_CC (1, 3)',  'vol_CC (1, 4)',
         'vol_CC (1, 5)', 'vol_CC (1, 6)',  'vol_CC (1, 7)']

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = titles[k])
    plt.title(label[i])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.show() 
    
#%%
titles = [ 'pH in 1','pH in 2', 'pH in 3','pH in 4',
          'pH in 5', 'pH in 6','pH in 7']
comp =  [ 'pH (1, 1)','pH (1, 2)', 'pH (1, 3)',  'pH (1, 4)',
         'pH (1, 5)', 'pH (1, 6)',  'pH (1, 7)']

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = titles[k])
    plt.title(label[i])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.show() 
#%%
comp =  [ 'De (1, 1)','De (1, 2)', 'De (1, 3)',  'De (1, 4)',
         'De (1, 5)', 'De (1, 6)',  'De (1, 7)']

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = titles[k])
    plt.xlabel('Time (s)')
    plt.yscale("log")
    plt.legend()
    plt.show()  
    
#%%
comp =  [ 'poros (1, 1)','poros (1, 2)', 'poros (1, 3)',  'poros (1, 4)',
         'poros (1, 5)', 'poros (1, 6)',  'poros (1, 7)']

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = titles[k])
    plt.xlabel('Time (s)')
    plt.yscale("log")
    plt.legend()
    plt.show()  
    
       
#%% Ca
comp = ['Ca (1, %s)'%i for i in range(1,9)]
title = ['Ca in %s'%i for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
    plt.title(label[i])
    plt.legend()
    plt.show() 
    #%% Ca
comp = ['C (1, %s)'%i for i in range(1,9)]
title = ['C in %s'%i for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
    plt.title(label[i])
    plt.legend()
    plt.show() 
#%% C
c_conc = np.load(path + 'C.npy')
plt.figure(figsize = (8,4))
plt.plot(c_conc[1,0:-2])
plt.xlabel(r'Length ($\mu$m)')
plt.ylabel('pH')
#%%
pH = np.load(path + 'pH.npy')
plt.figure(figsize = (8,4))
plt.plot(pH[1,0:-2])
plt.xlabel(r'Length ($\mu$m)')
plt.ylabel('pH')
#%%
'''
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['co2_uptake'],
             ls=linetype[i], label = label[i])
plt.xlabel('Time (s)')
plt.legend()
plt.show()
''' 
#results[names[0]]['C (1, 1)'][0:5]

