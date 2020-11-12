# -*- coding: utf-8 -*-
'''
The default case
'''
#%% MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import misc_func as fn
import plot_func as cf
#%% SETTINGS

fname = 'dc'
fpath = root_dir+'\\results\\output\\00_examples\\'#'\\results\\output\\00_default\\'
fn.make_output_dir(fpath)
name = '01_example_default'
path = root_dir+'\\results\\output\\00_examples\\' + name + '\\'
#name = '01_ie01_p05'
#path = root_dir+'\\results\\output\\07_internal_energy\\' + name + '\\' 
#path = root_dir+'\\results\\output\\13_validation\\' + name + '\\'
#name = '04_ccD13_Di'
#name = '05_ccD13_Di_PSCS'
linetype = '-'
scale = 50

results = {}
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

phaseCC = np.load(path +'CC.npy')
phaseCH = np.load(path +'CH.npy')
#%% MAIN PROPERTIES
r1 = 1000
r2 = 10000
titles = ['Portlandite', 'Calcite']#, 'Calcium', 'Carbon', 'Average pH',  'Porosity']
comp =  ['portlandite', 'calcite']#, 'Ca', 'C', 'pH',  'avg_poros']
suffix = ['_portlandite', '_calcite']#, '_calcium', '_carbon', '_average ph', '_input_c', '_poros']
ylabel = [r'Portlandite $\cdot 10^{-15}$ (mol)', r'Calcite  $\cdot 10^{-15}$ mol', 
          r'Dissolved Ca $\cdot 10^{-15}$ (mol)', r'Dissolved C $\cdot 10^{-15}$ (mol)' ,
          'Average pH',  'Porosity']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    #for i in range(0, len(names)):
    plt.plot(np.array(results['time'][:r2])/3600, results[comp[k]][:r2],
             ls=linetype)
    plt.ylabel(ylabel[k],fontsize = 14)
    plt.xlabel('Time (h)',fontsize = 14)
    plt.legend(fontsize = 12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 
    
    
#%% CC and CH together
fig, ax1 = plt.subplots(figsize=(8,4), dpi = 500)

color = 'tab:orange'
ax1.set_xlabel('Time (h)',fontsize = 14)
ax1.set_ylabel(r'Portlandite $\cdot 10^{-15}$ (mol)', color=color,fontsize = 14)
ax1.plot(np.array(results['time'])/3600, results['portlandite'], color=color)
ax1.tick_params(axis='y', labelcolor=color, which='major', labelsize=12)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel(r'Calcite  $\cdot 10^{-15}$ (mol)', color=color,fontsize = 14)  # we already handled the x-label with ax1
ax2.plot(np.array(results['time'])/3600, results['calcite'], color=color)
ax2.tick_params(axis='y', labelcolor=color, which='major', labelsize=12)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
#%% CC and CH profile

phaseCC = np.load(path +'CC.npy')
phaseCH = np.load(path +'CH.npy')
cc_mass = 100.09
ch_mass = 74.09
phaseCC[1,1:3] = phaseCC[1,3]
phaseCC[1,0:-1]=phaseCC[1,0:-1]*cc_mass*50/1000
phaseCH[1,0:-1]=phaseCH[1,0:-1]*ch_mass*50/1000
plt.figure(figsize=(8,4), dpi = 500)
plt.plot(range(0,len(phaseCC[1,:])-1),phaseCC[1,0:-1], label = 'CC')
plt.plot(range(0,len(phaseCH[1,:])-1),phaseCH[1,0:-1], label = 'CH')
plt.fill_between(np.arange(0,7,1),phaseCC[1,0:7], color = "#6f8191",
                 alpha=.5,label = "calcite")
plt.fill_between(np.arange(9,30,1),phaseCH[1,9:30], color = "#3d4249",
                 alpha=.5, label = "portlandite")
plt.ylabel(r'Mineral mass  $\cdot 10^{-12}$ (g) in 1 $\mu m^3$ ', fontsize=14)
plt.xlabel(r'Distance ($\mu m$)', fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend(fontsize=12)
#%% pH profile
nch = results['portlandite_cells'][-1]
ncc = results['calcite_cells'][-1]
plt.figure(figsize=(8,4), dpi = 500)
plt.plot(range(0,len(ph[1,:])-1),ph[1,0:-1], color = "#f2630d", label = 'pH')
plt.fill_between(np.arange(1,7,1),ph[1,1:7], color = "#6f8191", alpha=.5,label = "calcite")
plt.fill_between(np.arange(6,11,1),ph[1,6:11], color = "#cbd2d9", alpha=.5, label = "depleted portlandite")
plt.fill_between(np.arange(10,30,1),ph[1,10:30], color = "#3d4249", alpha=.5, label = "portlandite")
plt.ylabel('pH',fontsize = 14)
plt.legend(fontsize = 12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylim(6.0,12.6)
plt.xlabel(r'Distance ($\mu m$)',fontsize = 14)

#%% Ca and C profile 
concC[1, 10:16] = concC[1, 16]
plt.figure(figsize=(8,4), dpi = 500)
plt.plot(range(0,len(concCA[1,:])-1),concCA[1,0:-1], label = 'Ca')
plt.plot(range(0,len(concC[1,:])-1),concC[1,0:-1], label = 'C')
plt.fill_between(np.arange(1,7,1),concCA[1,1:7], color = "#6f8191",
                 alpha=.5,label = "calcite")
plt.fill_between(np.arange(6,11,1),concCA[1,6:11], color = "#cbd2d9", 
                 alpha=.5, label = "depleted portlandite")
plt.fill_between(np.arange(10,30,1),concCA[1,10:30], color = "#3d4249",
                 alpha=.5, label = "portlandite")
plt.ylabel('Ca and C molarity (mol/l)',fontsize = 14)
plt.xlabel(r'Distance ($\mu m$)',fontsize = 14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.yscale('log')
plt.legend(fontsize=12)
#%%
np.savetxt('case_study_profile.txt', 
           np.transpose([range(0,len(phaseCC[1,:])-1), phaseCH[1,0:-1], phaseCC[1,0:-1], ph[1,0:-1], concCA[1,0:-1], concC[1,0:-1]]), 
           header = 'distance;\t portlandite;\t calcite;\t ph;\t Ca;\t C', delimiter = '; ', fmt ='%f')


#%% C profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(concC[1,:])-1),concC[1,0:-1])
plt.plot(range(0,len(concCA[1,:])-1),concCA[1,0:-1], label = 'Ca')
plt.ylabel('C')
plt.xlabel(r'Distance ($\mu m$)')
plt.yscale('log')

print(concC[1,0])
print(concCA[1,0])

#%% SI profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(si[1,:])-1),si[1,0:-1])
plt.ylabel('SI')
plt.xlabel(r'Distance ($\mu m$)')

#%% DISSOLUTION RATE
titles = ['Dissolution rate', 'Precipitation rate' ]
comp =  ['portlandite', 'calcite']
suffix = ['_CH_rate', '_CC_rate' ]
s =1
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))    
    rate = np.abs(cf.get_rate(results[comp[k]],
                         results['time'][2] - results['time'][1],
                         step = s ))
    plt.plot(np.array(results['time'][::s] )/3600,
             rate*1e-3)
    plt.title(titles[k])
    plt.xlabel('Time (h)')
    plt.ylabel('Rate (mmol/l/s)')
    plt.yscale("log")
    plt.show()
#plt.savefig(fpath + fname + '_CH_rate')

#%% Change in porosity
plt.figure(figsize=(8,4))
cp =  np.array(results['avg_poros']) -results['avg_poros'][0]
print(cp[-1]*100)
plt.plot(results['time'],cp)
plt.xlabel('Time', fontsize = 14)
plt.ylabel('Change in porosity', fontsize = 14)
plt.legend(fontsize = 12)
plt.show
#%% smooth change in porosity
idx = np.array([0, 457,1905, 4605, 8833, 11748])
plt.figure(figsize=(8,4))
cp=[]
t = [] 
for i in range(0, len(idx)):
    cp.append((np.array(results['avg_poros'][idx[i]]) -results['avg_poros'][0])*100)
    t.append(results['time'][idx[i]]/3600)
print(cp[-1]*100)
plt.plot(t,cp)
plt.xlabel('Time (h)', fontsize = 14)
plt.ylabel('Change in porosity (%)', fontsize = 14)
plt.legend(fontsize = 12)
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

'''
plt.figure(figsize=(8,4))
plt.plot(results['time'], np.array(results['portlandite']) -
            results['portlandite'][0], label = 'CH')
plt.plot(results['time'], np.array(results['calcite']), label = 'CC')
plt.xlabel('Time')
plt.ylabel(r'Change in Portlandite ($\cdot 10^{-15}$ mol/l)')
plt.legend()
plt.show
'''
#%% porosity profile 
plt.figure(figsize=(8,4))
plt.plot(range(0,len(ph[1,:])-1),poros[1,0:-1])
plt.ylabel('Porosity')
plt.yscale("log")
plt.xlabel(r'Distance ($\mu m$)')


#%% Points
r1 = 1
r3 = 10000
p = 7
f = "C"
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
a = [0]
#p = 'poros (1, 8)'
comp =  [ 'poros (1, ' + str(i) + ')' for i in range(1,8) ]
thres = 1e-4
for p in comp:
    for i in range(1,len(results[p])-1):
        if np.abs(results[p][i] - results[p][i-1]) < thres:
            if np.abs(results[p][i+1] - results[p][i]) > thres:
                a.append(results['time'][i] )
b = np.sort(a)

a = np.delete(b,2) #a = np.delete(b,1)
#t = carb_rt.transition_time
x = np.arange(0, len(a))
plt.figure(figsize=(8,4))
plt.plot(np.sort(a)/3600, x)
plt.ylabel(r"Thickness ($\mu$m)",fontsize = 14)
plt.xlabel("Time (h)",fontsize = 14)
plt.tick_params(axis='both', which='major', labelsize=12)
#plt.ylim([0,10])
#plt.xlim([0,140])
plt.legend(fontsize = 12)
plt.show()

idx = np.array([0, 457,1905, 4605, 8833])
plt.figure(figsize=(8,4))
cp=[]
t = [] 
for i in range(0, len(idx)):
    cp.append(np.array(results['avg_poros'][idx[i]]*100))
    t.append(results['time'][idx[i]]/3600)
print(cp[-1]*100)
plt.plot(t,cp)
plt.xlabel('Time (h)', fontsize = 14)
plt.ylabel('Change in porosity (%)', fontsize = 14)
plt.legend(fontsize = 12)
plt.show
#%% porosity and thickness together
fig, ax1 = plt.subplots(figsize=(8,4), dpi = 500)

color = 'tab:orange'
ax1.set_xlabel('Time (h)',fontsize = 14)
ax1.set_ylabel('Average porosity (%)', color=color,fontsize = 14)
ax1.plot(t,cp, color = color)
ax1.tick_params(axis='y', labelcolor=color, which='major', labelsize=12)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel(r"Thickness ($\mu$m)", color=color,fontsize = 14)  # we already handled the x-label with ax1
ax2.plot(np.sort(a)/3600, x, color=color)
ax2.tick_params(axis='y', labelcolor=color, which='major', labelsize=12)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


#%%
dtime = np.array(sorted(a))[1:]
dtime = np.insert(dtime, 4, results['time'][-1])
indices = np.array([np.where(results['time']==dtime[i])[0][0] for i in range(0, len(dtime))])
dtime = np.insert(dtime, 0, 0)
indices = np.insert(indices, 0, 0)

titles = ['Dissolution rate', 'Precipitation rate' ]
comp =  ['portlandite', 'calcite']
suffix = ['_CH_rate', '_CC_rate' ]
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4), dpi = 500) 
    rate = np.zeros(len(indices)-1)
    for i in range(1, len(indices)):
        rate[i-1] = np.abs(results[comp[k]][indices[i]] - results[comp[k]][indices[i-1]])/\
                         (results['time'][indices[i]]  - results['time'][indices[i-1]]) *1e-3 
    plt.plot(dtime[0:-1]/3600,  rate)
    #plt.title(titles[k])
    plt.xlabel('Time (h)',fontsize = 14)
    plt.ylabel(r"Rate $(mmol / (l\cdot s\cdot cm^2)$",fontsize = 14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.yscale("log")
    plt.show()
    #%%
np.savetxt('case_study_rates_and_porosity.txt', 
           np.transpose([t, cp, x, rate]), 
           header = 'time;\t porosity;\t thickness;\t reaction rate', delimiter = '; ', fmt ='%f')
#%% Degree of carbonation
mm_CH = 74
mm_CC = 100
mass0 = results['portlandite'][0]*mm_CH + results['calcite'][0]*mm_CC
mass_f = results['portlandite'][-1]*mm_CH + results['calcite'][-1]*mm_CC

plt.figure(figsize=(8,4))
d = ((np.array(results['portlandite'])*mm_CH + \
              np.array(results['calcite'])*mm_CC)- \
              results['portlandite'][0]*mm_CH)*1e-15
plt.plot(np.array(results['time'])/24./3600, d)
plt.title('Degree of carbonation')
plt.xlabel('Time (d)')
plt.ylabel('Mass increase (g)')
plt.show

#%%
ch = 28.7
gain_per_voxel =(ch*0.9*mm_CC-ch*mm_CH)*1e-15
print(gain_per_voxel)
max_gain = (ch*mm_CC-ch*mm_CH)*1e-15
print(max_gain)

#%%

np.savetxt('case_study_time_series.txt', 
           np.transpose([np.array(results['time'])/3600, np.array(results['portlandite']), np.array(results['calcite']), d]), 
           header = 'time;\t portlandite;\t calcite;\t mass increase', delimiter = '; ')


