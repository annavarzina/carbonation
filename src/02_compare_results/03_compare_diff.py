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
import func as cf
#%% SETTINGS
Ts =1000.
fname = 'diff_p05_ll'
fpath = root_dir+'\\results\\output\\03_diffusivity\\'
fn.make_output_dir(fpath)
#names = np.array(['04_fixD_p005_arch', '01_fixD_p005', '02_fixD_p005_D11', '03_fixD_p005_D12'])#, '01_fixD_005p_D11'])
#names = np.array(['04_fixD_p05_arch', '01_fixD_p05', '02_fixD_p05_D11', '03_fixD_p05_D12'])#, '01_fixD_005p_D11'])
label = np.array(['Archies', '1e-10', '1e-11', '1e-12'])
linetype = np.array(['-', '--', '-.', ':'])
names = np.array(['03_fixD_p05_D12_ll_all', '05_fixD_p05_D15_ll_all'])
label = np.array(['1e-12', '1e-15'])
linetype = np.array(['-', '--'])
results = {}
for nn in names:
    path = root_dir+'\\results\\output\\03_diffusivity\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

#%% SCALE
    
scale = 50
for i in range(0, len(names)):
    temp = np.array(results[names[i]]['time'])
    temp *= scale
    results[names[i]]['time']= temp.tolist()
    
#%% CH DISSOLUTION 

titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
          'Average pH', 'Input C', 'Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'C (1, 0)', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_average ph', '_input_c', '_poros']
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
'''
titles = ['Dissolution rate', 'Precipitation rate' ]
comp =  ['portlandite', 'calcite']
suffix = ['_CH_rate', '_CC_rate' ]
rstart = 500
rend = len(results[names[1]]['time']) - 1

for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'][rstart:rend], 
                 cf.get_rate(results[names[i]][comp[k]][rstart:rend],
                             results[names[i]]['time'][2] - \
                             results[names[i]]['time'][1]),
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.ylabel('Rate (mol/s)')
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show()
#plt.savefig(fpath + fname + '_CH_rate')
#%%
k = 'calcite'# 'avg_poros' #'portlandite''calcite'
er = np.array([])
for nn in names:
    d_ch = results[nn][k][0]- results[nn][k][-2]
    er = np.append(er, d_ch)
print(er)
for i in range(0,len(er)):
    e = (er[i]-er[-1])/er[-1] * 100    
    print('%s, %s ' %(names[i], str(e)))
'''
#%%
n = np.sort(results[names[1]].keys())
points = ['(1, 0)', '(1, 1)', '(1, 2)']
titles = ['Carbon ', 'Calcium ', 'Effective diffusivity ', 
          'CC volume ', 'CH volume ', 'Porosity ']
comp =  [ 'C ', 'Ca ', 'De ', 'vol_CC ', 'vol_CH ', 'poros ']

tstart = 0
tend = 1*scale

for k in range(0, len(comp)):
    for j in range(0, len(points)):
        plt.figure(figsize=(8,4))
        for i in range(0, len(names)):
            indices = np.where(np.logical_and(np.array(results[names[i]]['time'])<=tend,
                          np.array(results[names[i]]['time'])>=tstart))
            plt.plot(np.array(results[names[i]]['time'])[indices], 
                     np.array(results[names[i]][comp[k]+points[j]])[indices],
                     ls=linetype[i], label = label[i])
        plt.title(titles[k]+points[j])
        plt.xlabel('Time (s)')
        plt.legend()
        plt.show() 