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
Ts =100.
fname = 'poros_Dfix'
fpath = root_dir+'\\results\\output\\05_porosity\\'
fn.make_output_dir(fpath)
'''
names = np.array(['01_poros001_subgrid', '02_poros005_subgrid',
                  '03_poros01_subgrid', '04_poros015_subgrid'])
names = np.array(['01_poros001_subgrid_poros', '02_poros005_subgrid_poros',
                  '03_poros01_subgrid_poros', '04_poros015_subgrid_poros'])
names = np.array([ '01_poros001_subgrid_lowfr','02_poros005_subgrid_lowfr',#
                  '03_poros01_subgrid_lowfr', '04_poros015_subgrid_lowfr'])
names = np.array([ '01_poros001_subgrid_lowfr_archie','02_poros005_subgrid_lowfr_archie',
                  '03_poros01_subgrid_lowfr_archie', '04_poros015_subgrid_lowfr_archie'])

names = np.array([ '02_poros005_subgrid_lowfr', '02_poros005_subgrid_lowfr_archie',
                   '02_poros005_subgrid', '02_poros005_subgrid_poros',
                   '02_poros005_Dborder11'])
label = np.array(['0.01','0.05', '0.1', '1.'])

names = np.array([ '03_poros01_subgrid_lowfr', '03_poros01_subgrid_lowfr_archie',
                   '03_poros01_subgrid', '03_poros01_subgrid_lowfr04'])
names = np.array([ '04_poros015_subgrid_lowfr', '04_poros015_subgrid_lowfr_archie',
                   '04_poros015_subgrid', '04_poros015_subgrid_lowfr04'])
names = np.array(['01_poros001_subgrid', '02_poros005_subgrid',
                  '03_poros01_subgrid', '04_poros015_subgrid'])
names = np.array([ '03_poros01_subgrid_lowfr',#, '03_poros01_subgrid_lowfr_archie',
                   '03_poros01_subgrid', '03_poros01_subgrid_lowfr04'])
names = np.array([ '02_poros005_subgrid_lowfr',# '02_poros005_subgrid_lowfr_archie',
                   '02_poros005_subgrid', '02_poros005_subgrid_lowfr04'])
names = np.array([ '01_poros001_subgrid_lowfr','02_poros005_subgrid_lowfr',#
                  '03_poros01_subgrid_lowfr', '04_poros015_subgrid_lowfr'])
label = np.array(['fr 0.01','fr 1.0', 'fr 0.004' ])
#label = np.array(['fr 0.01', 'fr 0.01 + A','fr 1.0', 'dyn fr', 'Dbord'])
label = np.array(['fr 0.01', 'fr 0.01 + A','fr 1.0', 'fr 0.004' ])
label = np.array(['0.01','0.05', '0.1', '0.15'])
names = np.array([ '02_poros005_subgrid', '02_poros005_subgrid_lowfr',
                   '02_poros005_Dborder09', '02_poros005_Dborder'])
label = np.array(['fr 1.0','fr 0.01', 'D 09', 'D 11' ])
names = np.array([ '02_poros005_Dborder09', #'01_poros001_Dborder09', 
                    '03_poros01_Dborder09','04_poros015_Dborder09', ])
'''
names = np.array(['01_p005_fr04_constD12','02_p01_fr04_constD12', '03_p015_fr04_constD12' ])
label = np.array(['0.05', '0.1', '0.15'])
linetype = np.array(['-', '--', '-.', ':', '-'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\05_porosity\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')
#%% SCALE
    
scale = 100
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
#%% POROSITY
field = 'avg_poros'
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'],
             np.array(results[names[i]][field])-results[names[i]][field][0],
             ls=linetype[i], label = label[i])
plt.ylabel('Change in porosity (%)')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_poros')
plt.show() 

#%% DISSOLUTION RATE

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
titles = [ 'Porosity in (1,1)','Porosity in (1,2)', 'Porosityin (1,3)','Porosity in (1,4)']
comp =  [ 'poros (1, 1)','poros (1, 2)', 'poros (1, 3)',  'poros (1, 4)']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.yscale("log")
    plt.legend()
    plt.show() 
#%%
titles = [ 'De in (1,1)','De in (1,2)', 'De in (1,3)','De in (1,4)']
comp =  [ 'De (1, 1)','De (1, 2)', 'De (1, 3)',  'De (1, 4)']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.yscale("log")
    plt.legend()
    plt.show() 
#%%
titles = [ 'Volume CH in (1,2)', 'Volume CH in (1,3)','Volume CH in (1,4)']
comp =  [ 'vol_CH (1, 2)', 'vol_CH (1, 3)',  'vol_CH (1, 4)']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.show() 
#%%
titles = [ 'Volume CC in (1,1)','Volume CC in (1,2)', 'Volume CC in (1,3)','Volume CC in (1,4)']
comp =  [ 'vol_CC (1, 1)','vol_CC (1, 2)', 'vol_CC (1, 3)',  'vol_CC (1, 4)']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.show() 
       
#%%
#'''
titles = [' Ca in (1,0)',' Ca in (1,1)' ,' Ca in (1,2)', 'Ca in (1,3)']
comp =  ['Ca (1, 0)','Ca (1, 1)' ,'Ca (1, 2)', 'Ca (1, 3)']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.show() 
    
#%% 
titles = [' C in (1,0)',' C in (1,1)' ,' C in (1,2)', 'C in (1,3)', 'C in (1,4)']
comp =  ['C (1, 0)','C (1, 1)' ,'C (1, 2)', 'C (1, 3)', 'C (1, 4)']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.show() 
#'''
#%%
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][0:500], results[names[i]]['vol_CC (1, 1)'][0:500],
             ls=linetype[i], label = label[i])
plt.title('CC (1, 1)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][0:500], results[names[i]]['vol_CC (1, 2)'][0:500],
             ls=linetype[i], label = label[i])
plt.title('CC (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][0:500], results[names[i]]['Ca (1, 1)'][0:500],
             ls=linetype[i], label = label[i])
plt.title('Ca (1, 1)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][0:500], results[names[i]]['Ca (1, 2)'][0:500],
             ls=linetype[i], label = label[i])
plt.title('Ca (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][0:500], results[names[i]]['C (1, 1)'][0:500],
             ls=linetype[i], label = label[i])
plt.title('C (1, 1)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


mvol = np.array([0.01, 0.05, 0.1, 0.15])
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][0:500], results[names[i]]['vol_CH (1, 2)'][0:500],
             ls=linetype[i], label = label[i])
plt.title('CH (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

#%%
for i in range(0, len(names)):
    print(results[names[i]]['De (1, 1)'][-1])
for i in range(0, len(names)):
    print(results[names[i]]['poros (1, 1)'][-1])

#%% ADDITIONAL
r1= 20
r2 =100
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['C (1, 0)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('C (1, 0)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

#%%
f = 'C (1, 1)'
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]][f][r1:r2],
             ls=linetype[i], label = label[i])
plt.title(f)
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

#%%
r1=0
r2 = 500
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], 
             np.array(results[names[i]]['vol_CH (1, 2)'][r1:r2])-results[names[i]]['vol_CH (1, 2)'][0],
             ls=linetype[i], label = label[i])
plt.title('CH volume (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['portlandite (1, 2)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('CH (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

#%%
r1=0
r2 = 500
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['calcite (1, 1)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('CC (1, 1)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['vol_CC (1, 1)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('CC volume (1, 1)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['vol_CC (1, 2)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('CC volume (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['Ca (1, 2)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('Ca (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['C (1, 2)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('C (1, 2)')
plt.ylim(0,1e-5)
plt.xlabel('Time (s)')
plt.legend()
plt.show() 
#%%
r1=0
r2 = 20
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['poros (1, 2)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('Porosity (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 


plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], 
             np.array(results[names[i]]['poros (1, 2)'][r1:r2])-results[names[i]]['poros (1, 2)'][0],
             ls=linetype[i], label = label[i])
plt.title('Porosity (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][r1:r2], results[names[i]]['poros (1, 1)'][r1:r2],
             ls=linetype[i], label = label[i])
plt.title('Porosity (1, 1)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 