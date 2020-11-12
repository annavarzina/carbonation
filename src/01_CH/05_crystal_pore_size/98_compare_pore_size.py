# -*- coding: utf-8 -*-
'''
Compare the results for different pore sizes (PS)
'''
#%% MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import misc_func as fn
import plot_func as cf
#%% SETTINGS
Ts =100.
fname = 'ps'
fpath = root_dir+'\\results\\output\\05_porosity\\'
#fpath = root_dir+'\\results\\output\\05_crystal_pore_size\\compare\\'
fn.make_output_dir(fpath)
#names = np.array([ '01_ie01_p05', '02_ie05_p05', '03_ie1_p05'])
#names = np.array([ '01_p0005_p005_Db09', '02_p001_p005_Db09', '03_p005_p005_Db09'])

names = np.array(['02_c07_p0005',  '01_c07_p001', '03_c07_p005'])
names = np.array(['02_c05_p0005', '01_c05_p001', '03_c05_p005'])
names = np.array(['02_co3_p0005', '01_co3_p001', '03_co3_p005'])
label = np.array(['PS 0.005','PS 0.01', 'PS 0.05'])
linetype = np.array(['-', '--', '-.',':'])



names = np.array(['01_porosity_0.05', '02_porosity_0.1', '03_porosity_0.25', '04_porosity_0.5'])
label = np.array(['pore size 0.005','pore size 0.01', 'pore size 0.025', 'pore size 0.05'])

results = {}
for nn in names:
    #path = root_dir+'\\results\\output\\05_crystal_pore_size\\' + nn + '\\'
    path = root_dir+'\\results\\output\\05_porosity\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

#%% SCALE    
scale = 50
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
    plt.ylabel(ylabel[k])
    plt.xlabel('Time (h)')
    plt.legend()
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
f = ('Porosity', 'poros')# ('Ca', 'Ca') ('C', 'C') ('Volume CH', 'vol_CH') ('Volume CC', 'vol_CC') ('pH', 'pH') ('De', 'De')('Porosity', 'poros')
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
#%% POINTS log
f = ('Porosity', 'poros')#   ('De', 'De')('Porosity', 'poros')
title = ['%s in %s'%(f[0],i) for i in range(1,9)]
comp = ['%s (1, %s)'%(f[1],i) for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
    plt.yscale("log")
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
