# -*- coding: utf-8 -*-
'''
Compare the results for different crystal sizes (CS)
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
Ts =100.
fname = 'cs'
fpath = root_dir+'\\results\\output\\05_crystal_pore_size\\compare\\'
fn.make_output_dir(fpath)
#names = np.array([ '01_ie01_p05', '02_ie05_p05', '03_ie1_p05'])
#names = np.array([ '01_c03_p005_Db09', '02_c05_p005_Db09', '03_c07_p005_Db09'])
#names = np.array([ '03_p0005_c07', '03_p001_c07', '03_p0025_c07', '03_p005_c07'])
#label = np.array(['0.005','0.01', '0.025', '0.05'])

names = np.array(['02_co3_p0005', '02_c05_p0005', '02_c07_p0005'])
names = np.array(['03_co3_p005', '03_c05_p005', '03_c07_p005'])
names = np.array(['01_co3_p001', '01_c05_p001', '01_c07_p001'])
label = np.array(['CS 0.3', 'CS 0.5', 'CS 0.7'])
linetype = np.array(['-', '--', '-.',':'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\05_crystal_pore_size\\' + nn + '\\'
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
'''
f = ('Ca', 'Ca') # ('C', 'C') ('Volume CH', 'vol_CH') ('Volume CC', 'vol_CC') ('pH', 'pH') ('De', 'De')('Porosity', 'poros')
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
'''
#%% COMPARE POINT for cases
'''
f = 'poros (1, 5)'
plt.figure(figsize=(8,4))  
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]][f],
             ls=linetype[i], label = label[i])
plt.xlabel('Time (s)')
plt.ylabel('Porosity in voxel 5')
plt.yscale("log")
plt.legend()
plt.show() 

for i in range(0, len(names)):
    print(results[names[i]][f][-1])
'''
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

#%%
if(False):
    names = np.array([ '03_p0005_c07', '03_p001_c07', '03_p0025_c07', '03_p005_c07'])
    tdict = {}
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        print(i)
        a = [0]
        comp =  ['poros (1, ' + str(m) + ')' for m in range(1,5) ]
        thres = 1e-4
        for p in comp:
            for j in range(1,len(results[names[i]][p])-1):
                if np.abs(results[names[i]][p][j] - results[names[i]][p][j-1]) < thres:
                    if np.abs(results[names[i]][p][j+1] - results[names[i]][p][j]) > thres:
                        a.append(results[names[i]]['time'][j] )
        print(np.sort(a))
        x = np.arange(0, len(a[:-2]))
        plt.plot(np.sort(a)[:-2]/3600, x, label = 'pore size ' + label[i], ls = linetype[i])
    plt.ylabel(r"Thickness ($\mu m$)",fontsize = 14)
    plt.xlabel("Time (h)",fontsize = 14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.legend(fontsize = 12)
    plt.show()