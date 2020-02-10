# -*- coding: utf-8 -*-
'''
The default case
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

fname = 'dc'
fpath = root_dir+'\\results\\output\\00_default\\'
fn.make_output_dir(fpath)
#names = np.array([ '01_ie01_p05', '02_ie05_p05', '03_ie1_p05'])
name = '02_c05_p005_Db09'
name = '04_ccD13_Di'
#name = '05_ccD13_Di_PSCS'
linetype = '-'
scale = 100

results = {}
#path = root_dir+'\\results\\output\\09_crystal_size\\' + name + '\\' 
path = root_dir+'\\results\\output\\13_validation\\' + name + '\\'
results = fn.load_obj(path + name +'_results')
for n in ['time', 'portlandite', 'calcite']:
    temp = np.array(results[n])
    temp *= scale 
    results[n]= temp.tolist()

concCA = np.load(path +'Ca.npy')
concC = np.load(path +'C.npy')
ph = np.load(path +'pH.npy')
de = np.load(path +'De.npy')
poros = np.load(path +'poros.npy')
si = np.load(path +'SI.npy')

#%% MAIN PROPERTIES
r1 = 1000
r2 = 10000
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
          'Average pH', 'Input C', 'Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'pH', 
         'C (1, 0)', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_average ph', '_input_c', '_poros']
ylabel = [r'Portlandite $\cdot 10^{-15}$ (mol)', r'Calcite  $\cdot 10^{-15}$ mol', 
          r'Dissolved Ca $\cdot 10^{-15}$ (mol)', r'Dissolved C $\cdot 10^{-15}$ (mol)' ,
          'Average pH', 'Boundary C (mol/l)', 'Porosity']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    #for i in range(0, len(names)):
    plt.plot(np.array(results['time'][:r2])/3600, results[comp[k]][:r2],
             ls=linetype)
    plt.ylabel(ylabel[k])
    plt.xlabel('Time (h)')
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 

#%% pH profile
nch = results['portlandite_cells'][-1]
ncc = results['calcite_cells'][-1]
plt.figure(figsize=(8,4))
plt.plot(range(0,len(ph[1,:])-1),ph[1,0:-1], color = "#f2630d", label = 'pH')
plt.fill_between(np.arange(1,7,1),ph[1,1:7], color = "#6f8191", alpha=.5,label = "calcite")
plt.fill_between(np.arange(6,11,1),ph[1,6:11], color = "#cbd2d9", alpha=.5, label = "depleted portlandite")
plt.fill_between(np.arange(10,30,1),ph[1,10:30], color = "#3d4249", alpha=.5, label = "portlandite")
plt.ylabel('pH')
plt.legend()
plt.ylim(6.0,12.6)
plt.xlabel(r'Distance ($\mu m$)')
#%% DISSOLUTION RATE
titles = ['Dissolution rate', 'Precipitation rate' ]
comp =  ['portlandite', 'calcite']
suffix = ['_CH_rate', '_CC_rate' ]
s =2
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))    
    rate = np.abs(cf.get_rate(results[comp[k]],
                         results['time'][2] - results['time'][1],
                         step = s ))
    plt.plot(np.array(results['time'][::s] )/3600,
             rate)
    plt.title(titles[k])
    plt.xlabel('Time (h)')
    plt.ylabel('Rate (mol/s)')
    plt.yscale("log")
    plt.show()
#plt.savefig(fpath + fname + '_CH_rate')

#%% Change in porosity
plt.figure(figsize=(8,4))
cp =  np.array(results['avg_poros']) -results['avg_poros'][0]
print(cp[-1]*100)
plt.plot(results['time'],cp)
plt.xlabel('Time')
plt.ylabel('Change in porosity')
plt.legend()
plt.show
#%%
plt.figure(figsize=(8,4))
plt.plot(results['time'], np.array(results['portlandite']) -
            results['portlandite'][0])
plt.xlabel('Time')
plt.ylabel(r'Change in Portlandite ($\cdot 10^{-15}$ mol/l)')
plt.legend()
plt.show

plt.figure(figsize=(8,4))
plt.plot(results['time'], np.array(results['calcite']))
plt.xlabel('Time')
plt.ylabel(r'Change in Calcite ($\cdot 10^{-15}$ mol/l)')
plt.legend()
plt.show
print(np.array(results['portlandite'][-1]) -
            results['portlandite'][0])
print(np.array(results['calcite'][-1]) -
            results['calcite'][0])


#%% porosity profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(ph[1,:])-1),poros[1,0:-1])
plt.ylabel('Porosity')
plt.yscale("log")
plt.xlabel(r'Distance ($\mu m$)')

#%% Ca profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(concCA[1,:])-1),concCA[1,0:-1], label = 'Ca')
plt.fill_between(np.arange(1,7,1),concCA[1,1:7], color = "#6f8191", alpha=.5,label = "calcite")
plt.fill_between(np.arange(6,11,1),concCA[1,6:11], color = "#cbd2d9", alpha=.5, label = "depleted portlandite")
plt.fill_between(np.arange(10,30,1),concCA[1,10:30], color = "#3d4249", alpha=.5, label = "portlandite")
plt.ylabel('Dissolved Ca (mol/l)')
plt.xlabel(r'Distance ($\mu m$)')
plt.legend()
#%% C profile 
plt.figure(figsize=(8,4))
plt.plot(range(1,len(concC[1,:])-1),concC[1,1:-1])
plt.ylabel('C')
plt.xlabel(r'Distance ($\mu m$)')

#%% SI profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(si[1,:])-1),si[1,0:-1])
plt.ylabel('SI')
plt.xlabel(r'Distance ($\mu m$)')

#%% CC profile
cc_mass = 100.09
ch_mass = 74.09
cc = []
ch = []
for i in range(1,10):
    cc.append(results['calcite (1, ' + str(i) +')' ][-1]*cc_mass)
    ch.append(results['portlandite (1, ' + str(i) +')' ][-1]*ch_mass)
chmax = results['portlandite (1, ' + str(i) +')' ][1]*ch_mass
for i in range(10,30):
    cc.append(0)
    ch.append(chmax)
    
plt.figure(figsize=(8,4))
plt.plot(range(1,30),cc[:])
plt.plot(range(1,30),ch[:])
plt.fill_between(np.arange(1,7,1),cc[0:6], color = "#6f8191", alpha=.5,label = "calcite")
plt.fill_between(np.arange(8,30,1),ch[7:30], color = "#3d4249", alpha=.5, label = "portlandite")
plt.ylabel(r'Mineral mass  $\cdot 10^{-15}$ (g) in 1 $\mu m^3$ ')
plt.xlabel(r'Distance ($\mu m$)')
plt.legend()

#%% Points
r1 = 1
r3 = 10000
p = 7
f = "Ca"
comp = ['%s (1, %s)'%(f,i) for i in range(1,p)]
title = ['%s in %s'%(f,i) for i in range(1,p)]

plt.figure(figsize=(8,4))    
for k in range(0, len(comp)):
    plt.plot(np.array(results['time'][r1:r2])/3600, results[comp[k]][r1:r2],
             ls=linetype, label = title[k])
plt.xlabel('Time (h)')
plt.ylabel('Dissolved'+ f +' (mol/l)')
plt.legend(loc = "upper right")
plt.show() 


#%%
a = []
#p = 'poros (1, 8)'
comp =  [ 'poros (1, ' + str(i) + ')' for i in range(1,15) ]
thres = 1e-4
for p in comp:
    for i in range(1,len(results[p])-1):
        if np.abs(results[p][i] - results[p][i-1]) < thres:
            if np.abs(results[p][i+1] - results[p][i]) > thres:
                a.append(results['time'][i] )
#t = carb_rt.transition_time
x = np.arange(0, len(a))
plt.figure(figsize=(8,4))
plt.plot(sorted(a), x)
plt.ylabel("Dissolved length (um)")
plt.xlabel("Time")
#plt.ylim([0,10])
#plt.xlim([0,140])
plt.legend()
plt.show()

#%%
dtime = np.array(sorted(a))
indices = np.array([np.where(results['time']==dtime[i])[0][0] for i in range(0, len(dtime))])
dtime = np.insert(dtime, 0, 0)
indices = np.insert(indices, 0, 0)

titles = ['Dissolution rate', 'Precipitation rate' ]
comp =  ['portlandite', 'calcite']
suffix = ['_CH_rate', '_CC_rate' ]
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4)) 
    rate = np.zeros(len(indices)-1)
    for i in range(1, len(indices)):
        rate[i-1] = np.abs(results[comp[k]][indices[i]] - results[comp[k]][indices[i-1]])/\
                         (results['time'][indices[i]]  - results['time'][indices[i-1]])  
    plt.plot(dtime[1:]/3600,  rate)
    plt.title(titles[k])
    plt.xlabel('Time (h)')
    plt.ylabel('Rate (mol/l/s/um2)')
    plt.yscale("log")
    plt.show()