# -*- coding: utf-8 -*-
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
#%% SETTINGS
Ts =1000.
scale = 50
fname = '01_compare_sg'
fpath = root_dir+'\\results\\temp\\01_subgrid\\'
fn.make_output_dir(fpath)
linetype = np.array(['-', '-', '-',# ':', '-', '--', '-.', ':',
                     '--', '--', '--', #':', '-', '--', '-.', ':'])
                     '-.', '-.', '-.',])
color = np.array(['r']) 
names = np.array(['01_test_sg_p01'])
#names = np.array(['01_test_D_p001_10', '01_test_D_p005_10', '01_test_D_p01_10',
#                  '01_test_Db_p001_10', '01_test_Db_p005_10', '01_test_Db_p01_10',
#                  '01_test_sg_p001', '01_test_sg_p005', '01_test_sg_p01'])
#names = np.array(['01_test_D_p001_Caeq', '01_test_D_p005_Caeq', '01_test_D_p01_Caeq',
#                  '01_test_Db_p001_Caeq', '01_test_Db_p005_Caeq', '01_test_Db_p01_Caeq',
#                  '01_test_sg_p001_Caeq', '01_test_sg_p005_Caeq', '01_test_sg_p01_Caeq'])
porosity = np.array([ '0.1'])
case =  np.array(['Subgrid'])
label = np.array(['Subgrid p0.1',])
results = {}
for i in range(0, len(names)):
    
    path = fpath + names[i] + '\\'
    results[names[i]] = fn.load_obj(path +  names[i] +'_results')
    
for i in range(0, len(names)):
    temp = np.array(results[names[i]]['time'])
    temp *= scale
    results[names[i]]['time']= temp.tolist()
#%%    
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon','Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon','_poros']
position = ['lower left','upper left','lower right' , 'upper right', 'upper right']
for k in range(0, len(comp)):
    plt.figure(figsize=(10,6))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i], color = color[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend(loc=position[k])
    #plt.savefig(fpath + fname + suffix[k])
    plt.show() 
    
plt.figure(figsize=(10,6))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], 
             (results[names[i]]['avg_poros'] - results[names[i]]['avg_poros'][0]),
             ls=linetype[i], label = label[i], color = color[i])
plt.title('Porosity')
plt.xlabel('Time (s)')
plt.ylabel('Porosity change (%)')
plt.legend(loc='lower left')
#plt.savefig(fpath + fname + suffix[k])
plt.show() 

  
#%% PROPERTIES in node (1,1)
titles = [ 'Calcite', 'Calcium', 'Carbon', 'pH']
comp =  [ 'calcite (1, 1)', 'Ca (1, 1)', 'C (1, 1)','pH (1, 1)']
suffix = [ '_calcite', '_calcium', '_carbon','_poros','_pH']
position = ['upper left','lower right' , 'upper right','lower right']
for k in range(0, len(comp)):
    plt.figure(figsize=(10,6))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i], color = color[i])
    plt.title(titles[k] + ' (1, 1)')
    plt.xlabel('Time (s)')
    plt.legend(loc=position[k])
    #plt.savefig(fpath + fname + suffix[k])
    plt.show() 
    
#%% PROPERTIES in node (1,2)
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon', 'pH']
comp =  ['portlandite (1, 2)', 'calcite (1, 2)', 'Ca (1, 2)', 'C (1, 2)','pH']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon','_poros','_pH']
position = ['lower left','upper left','lower right' , 'upper right','lower right']
for k in range(0, len(comp)):
    plt.figure(figsize=(10,6))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i], color = color[i])
    plt.title(titles[k] + ' (1, 2)')
    plt.xlabel('Time (s)')
    plt.legend(loc=position[k])
    #plt.savefig(fpath + fname + suffix[k])
    plt.show() 
    
#%% PROPERTIES in node (1,3)
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon']
comp =  ['portlandite (1, 3)', 'calcite (1, 3)', 'Ca (1, 3)', 'C (1, 3)']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon','_poros']
position = ['lower left','upper left','lower right' , 'upper right']
for k in range(0, len(comp)):
    plt.figure(figsize=(10,6))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i], color = color[i])
    plt.title(titles[k] + ' (1, 3)')
    plt.xlabel('Time (s)')
    plt.legend(loc=position[k])
    #plt.savefig(fpath + fname + suffix[k])
    plt.show() 
    
#%% PROPERTIES in node (1,0)
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon']
comp =  ['portlandite (1, 0)', 'calcite (1, 0)', 'Ca (1, 0)', 'C (1, 0)']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon','_poros']
position = ['lower left','upper left','lower right' , 'upper right']
for k in range(0, len(comp)):
    plt.figure(figsize=(10,6))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i], color = color[i])
    plt.title(titles[k] + ' (1, 3)')
    plt.xlabel('Time (s)')
    plt.legend(loc=position[k])
    #plt.savefig(fpath + fname + suffix[k])
    plt.show() 