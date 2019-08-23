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
fname = 'diff'
fpath = root_dir+'\\results\\output\\compare_diffusivity\\'
fn.make_output_dir(fpath)
#names = np.array(['05_mvol_40', '01_reference', '05_mvol_10', '05_mvol_2', '05_mvol_1'])
#label = np.array(['0.331*40', '0.331*20','0.331*10', '0.331*2', '0.331'])
#linetype = np.array(['-', '--', '-.', ':', '-'])
names = np.array(['01_fixD_p005_arch', '01_dixD_p005', '01_dixD_p005_D11', '01_dixD_p005_D12'])#, '01_fixD_005p_D11'])
label = np.array(['Archies', '1e-10', '1e-11', '1e-12'])
linetype = np.array(['-', '--', '-.', ':'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\diffusivity\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

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
