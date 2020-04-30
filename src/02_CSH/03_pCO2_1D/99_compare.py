# -*- coding: utf-8 -*-
'''
Compare the results for different CSH
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
#%% SETTINGS
Ts =1000.
fname = 'compare'
fpath = root_dir+'\\results\\temp\\01_develop\\'
fn.make_output_dir(fpath)

names = np.array(['01_1D_1voxel_3.4', '01_1D_1voxel_2.0',
                  '01_1D_1voxel_1.0', '02_1D_3.4'])
label = np.array([ 'CO2 0.03% 1vox', 'CO2 1% 1vox',
                  'CO2 10% 1vox', 'CO2 0.03% 10vox'])
scale = 10
linetype = np.array(['-', '--', '-.', ':', '-'])
results = {}
for nn in names:
    path = root_dir+'\\results\\temp\\01_develop\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

print(np.sort(results[names[0]].keys()))
#%% SCALE
#dtime = [350,350,350,350,200,0]
keys = ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'time', 'avg_poros']
sres = {}
for i in range(0, len(names)):
    temp = {}
    a =np.array(results[names[i]]['time']) <= Ts/scale
    s = np.size(results[names[i]]['time'])
    for k in keys:
        if(np.size(results[names[i]][k])==s):
            temp[k] = np.array(results[names[i]][k])[a]
    temp['time'] *= scale 
    #temp['time']+=dtime[i]
    temp['portlandite'] *=scale
    temp['calcite'] *=scale
    sres[names[i]] = temp
#%% MAIN PROPERTIES
titles = ['Calcite', 'Calcium', 'Carbon','Average pH']
comp =  ['calcite', 'Ca', 'C', 'pH']
suffix = ['_calcite', '_calcium', '_carbon','_average ph']
ylabel = [r'Calcite  $\cdot 10^{-15}$ mol', 
          r'Dissolved Ca $\cdot 10^{-15}$ (mol)', 
          r'Dissolved C $\cdot 10^{-15}$ (mol)' ,
          'Average pH']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    #for i in range(0, len(names)):
    for i in range(0, len(names)):
        plt.plot(sres[names[i]]['time'], sres[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.ylabel(ylabel[k], fontsize=14)
    #plt.title(titles[k])
    plt.xlabel('Time (s)', fontsize=14)
    plt.legend(fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 
#%%
titles = ['CSHQ_JenD','CSHQ_JenH','CSHQ_TobD','CSHQ_TobH','Ca_solid','Si_solid' ]
comp =  ['CSHQ_JenD','CSHQ_JenH','CSHQ_TobD','CSHQ_TobH','Ca_solid','Si_solid' ]
suffix = ['CSHQ_JenD','CSHQ_JenH','CSHQ_TobD','CSHQ_TobH','Ca_solid','Si_solid' ]
ylabel = ['CSHQ_JenD','CSHQ_JenH','CSHQ_TobD','CSHQ_TobH','Ca_solid','Si_solid']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    #for i in range(0, len(names)):
    for i in [0,1,2]:
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.ylabel(ylabel[k], fontsize=14)
    #plt.title(titles[k])
    plt.xlabel('Time (s)', fontsize=14)
    plt.legend(fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig(fpath + fname + suffix[k])
    plt.show()
#%%
titles =  ['C (1, 3)','Ca (1, 3)','C (1, 4)','Ca (1, 4)','C (1, 5)','Ca (1, 5)','Si (1, 5)' ]
comp =  ['C (1, 3)','Ca (1, 3)','C (1, 4)','Ca (1, 4)','C (1, 5)','Ca (1, 5)','Si (1, 5)' ]
suffix =  ['C (1, 3)','Ca (1, 3)','C (1, 4)','Ca (1, 4)','C (1, 5)','Ca (1, 5)','Si (1, 5)' ]
ylabel =  ['C (1, 3)','Ca (1, 3)','C (1, 4)','Ca (1, 4)','C (1, 5)','Ca (1, 5)','Si (1, 5)' ]
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    #for i in range(0, len(names)):
    for i in [0,1,2]:
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.ylabel(ylabel[k], fontsize=14)
    #plt.title(titles[k])
    plt.xlabel('Time (s)', fontsize=14)
    plt.legend(fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 
#%% plot ca/si against density
plt.figure()
for i in [0,1,2]:
    plt.plot(results[names[i]]['Ca_Si'], results[names[i]]['csh_density'],
                 ls=linetype[i], label = label[i])
plt.legend()
plt.ylabel('CSH density')
plt.xlabel('C/S')
plt.show()

#%% plot CSH
plt.figure()
for i in [0,1,2,3]:
    plt.plot(results[names[i]]['time'], results[names[i]]['CSHQ_JenD (1, 5)'],
                 ls=linetype[0], label = label[i])
plt.legend()
plt.ylabel('CSHQ_JenD')
plt.xlabel('Time')
plt.show()

plt.figure()
for i in [0,1,2,3]:
    plt.plot(results[names[i]]['time'], results[names[i]]['CSHQ_JenH (1, 5)'],
                 ls=linetype[0], label = label[i])
plt.legend()
plt.ylabel('CSHQ_JenH')
plt.xlabel('Time')
plt.show()

plt.figure()
for i in [0,1,2,3]:
    plt.plot(results[names[i]]['time'], results[names[i]]['CSHQ_TobD (1, 5)'],
                 ls=linetype[0], label = label[i])
plt.legend()
plt.ylabel('CSHQ_TobD')
plt.xlabel('Time')
plt.show()


plt.figure()
for i in [0,1,2,3]:
    plt.plot(results[names[i]]['time'], results[names[i]]['CSHQ_TobH (1, 5)'],
                 ls=linetype[0], label = label[i])
plt.legend()
plt.ylabel('CSHQ_TobH')
plt.xlabel('Time')
plt.show()