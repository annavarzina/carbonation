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
fname = '07_set1_compare_100s'
fpath = root_dir+'\\results\\output\\09_set_1\\'
fn.make_output_dir(fpath)
'''
paths = np.array([root_dir+'\\results\\output\\05_porosity\\',
                  root_dir+'\\results\\output\\05_porosity\\',
                  root_dir+'\\results\\output\\05_porosity\\',
                  root_dir+'\\results\\output\\05_porosity\\',
                  root_dir+'\\results\\output\\04_liquid_layer\\',
                  root_dir+'\\results\\output\\04_liquid_layer\\'])
names = np.array(['01_porosity_0.05archie', '01_porosity_0.05', 
                  '04_porosity_0.5archie', '04_porosity_0.5', 
                  '01_ll4_p05', '01_ll4_p05_Dfix'])

p005_A_ll0 = root_dir+'\\results\\output\\05_porosity\\01_porosity_0.05archie'
p005_F_ll0 = root_dir+'\\results\\output\\05_porosity\\01_porosity_0.05'
p05_A_ll0 = root_dir+'\\results\\output\\05_porosity\\04_porosity_0.5archie'
p05_F_ll0 = root_dir+'\\results\\output\\05_porosity\\04_porosity_0.5'
p05_A_ll4 = root_dir+'\\results\\output\\04_liquid_layer\\01_ll4_p05'
p05_F_ll4 = root_dir+'\\results\\output\\04_liquid_layer\\01_ll4_p05_Dfix'
'''
    
paths = np.repeat(root_dir+'\\results\\output\\09_set_1\\', 8) 
names = np.array(['01_p005_A_ll0_1000s','02_p005_F_ll0_1000s',
                  '03_p05_A_ll0_1000s','04_p05_F_ll0_1000s',
                  '05_p05_A_ll4_1000s','06_p05_F_ll4_1000s',
                  '07_p005_A_ll4_1000s', '08_p005_F_ll4_1000s'])
porosity = np.array(['0.05', '0.05', '0.5', '0.5', '0.5', '0.5', '0.05', '0.05'])
diffusivity =  np.array(['A', 'F', 'A', 'F', 'A', 'F', 'A', 'F'])
liqlayer =  np.array(['0', '0', '0', '0', '4', '4', '4', '4'])
label = np.array(['p0.05; A; l0', 'p0.05; F; l0', 'p0.5;  A; l0', 
                  'p0.5;  F; l0', 'p0.5;  A; l4', 'p0.5;  F; l4',
                  'p0.05;  A; l4', 'p0.05;  F; l4'])
linetype = np.array(['-', '--', '-.', ':', '-', '--', '-.', ':'])

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
    plt.legend(loc='upper right')
    plt.savefig(fpath + fname + suffix[k])
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
          'CC volume ', 'CH volume ', 'Porosity ', 'pH']
comp =  [ 'C ', 'Ca ','De ', 'vol_CC ', 'vol_CH ', 'poros ', 'pH ']
ylims = [(0,0.1), (0, 0.08), (1e-14,1e-8),
         (0, 1), (0, 1), (0,1), (0,14)]
# Carbon 


#for k in range(0, len(comp)):    
    #for i in range(0, len(names)):
k=3
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
              '. D:' + str(diffusivity[i]) +\
              '. L:' + str(liqlayer[i])) #names[i]
    plt.ylim(ylims[k])
    plt.xlabel('Time (s)')
    plt.legend(loc='upper right')
    plt.show() 
