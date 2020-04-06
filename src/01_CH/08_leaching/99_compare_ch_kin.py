# -*- coding: utf-8 -*-
'''
Compare the results:
        poros  D     liqlayer
    1)  0.05 Archie     0
    2)  0.05 Fixed(-11) 0
    1)  0.5 Archie      0
    2)  0.5 Archie      4
    1)  0.5 Fixed(-11   0
    2)  0.5 Fixed(-11)  4
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
scale = 50
fname = '10_compare_ch_kin'
fpath = root_dir+'\\results\\output\\10_ch_kinetics\\'
fn.make_output_dir(fpath)
linetype = np.array(['-', '--', '-.', ':', '-', '--', '-.', ':'])

paths = np.repeat(root_dir+'\\results\\output\\10_ch_kinetics\\', 6) 
#names = np.array(['01_p005_ll0_10s', '02_p005_ll4_10s',
#                  '03_p05_ll0_10s', '04_p05_ll4_10s'])
#names = np.array(['01_p005_ll0_100s', '02_p005_ll4_100s',
#                  '03_p05_ll0_100s', '04_p05_ll4_100s'])
names = np.array(['01_p005_ll0_1000s', '01_p005_ll1_1000s',
                  '02_p005_ll4_1000s', '03_p05_ll0_1000s',
                  '03_p05_ll1_1000s', '04_p05_ll4_1000s'])
porosity = np.array([ '0.05','0.05', '0.05', '0.5', '0.5', '0.5'])
liqlayer =  np.array(['0','1','4','0', '1', '4'])
label = np.array(['p0.05; l0', 'p0.05; l1',
                  'p0.05; l4', 'p0.5; l0',
                  'p0.5; l1', 'p0.5; l4'])
results = {}
for i in range(0, len(names)):
    
    path = paths[i] + names[i] + '\\'
    results[names[i]] = fn.load_obj(path +  names[i] +'_results')
    
for i in range(0, len(names)):
    temp = np.array(results[names[i]]['time'])
    temp *= scale
    results[names[i]]['time']= temp.tolist()
    
#%% PROPERTIES
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
           'Input C', 'Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'C (1, 0)', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_input_c', '_poros']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend(loc='upper right')
    #plt.savefig(fpath + fname + suffix[k])
    plt.show() 
#%%
n = np.sort(results[names[1]].keys())

tstart = 0
tend = 20*scale

linetype = np.array(['-', '--', '-.', ':', '-', '--', '-.', ':',
                     '-', '--', '-.', ':', '-', '--', '-.', ':'])
points = [(1, 0), (1, 1), (1, 2), (1,3), (1,4), (1,5), (1,6), (1,7),
          (1,8)]#, (1,9), (1,10), (1,11), (1,12), (1,13), (1,14)]
label = [str(points[i][1]) for i in range(0,len(points))]

titles = ['Carbon ', 'Calcium ','Effective diffusivity ',
          'CC volume ', 'CH volume ']
comp =  [ 'C ', 'Ca ','De ', 'vol_CC ', 'vol_CH ']
ylims = [(0,0.1), (0, 0.08), (5e-13,5e-9),
         (0, 1), (0, 1)]
# Carbon 


#for k in range(0, len(comp)):    
    #for i in range(0, len(names)):
k=2

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))
    for j in range(0, len(points)):
        indices = np.where(np.logical_and(np.array(results[names[i]]['time'])<=tend,
                      np.array(results[names[i]]['time'])>=tstart))
    
        plt.plot(np.array(results[names[i]]['time'])[indices], 
                 np.array(results[names[i]][comp[k]+str(points[j])])[indices],
                 ls=linetype[j], label = label[j])
        if k==2:
        
            plt.yscale('log')
    plt.title(titles[k]+'. P:' + str(porosity[i]) +\
              #'. D:' + str(diffusivity[i]) +\
              '. L:' + str(liqlayer[i])) #names[i]
    plt.ylim(ylims[k])
    plt.xlabel('Time (s)')
    plt.legend(loc='upper right')
    plt.show() 
