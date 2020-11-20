# -*- coding: utf-8 -*-
'''
Compare the results for different PCO2
'''
#%% MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(root_dir)
sys.path.append(src_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import misc_func as fn
import plot_func as cf
#%% CONSOLE DEBUGGING
# In IPython shell use this:
'''
import sys,os
root_dir = os.getcwd()
src_dir = root_dir + '\\src'
csh_dir = src_dir +'\\02_CSH'
sys.path.append(root_dir)
sys.path.append(src_dir)
sys.path.append(csh_dir)
'''
#%% SETTINGS
Ts =3600. * 3.
fname = 'csh'
fpath = root_dir+'\\results\\output_csh\\02_24h\\compare\\'
fn.make_output_dir(fpath)
names = np.array(['02_D11', '01_D12', '02_D13','03_D12_CSH-archie'])
label = np.array(['D_CC-11', 'D_CC-12', 'D_CC-13', 'D_CC-12_D_CSH-A'])
linetype = np.array(['-', '--', '-.', ':', '-', '--','-', '--', '-.', ':', '-', '--'])

scale = 50
results = {}
for nn in names:
    path = root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')
#%% SCALE
    
for i in range(0, len(names)):
    temp = np.array(results[names[i]]['time'])
    temp *= scale/3600
    results[names[i]]['time']= temp.tolist()
    #temp = np.array(results[names[i]]['portlandite'])
    #temp *= scale
    #results[names[i]]['portlandite']= temp.tolist()
    #temp = np.array(results[names[i]]['calcite'])
    #temp *= scale
    #results[names[i]]['calcite']= temp.tolist()
    

#%% PARAMS
s = 5
titles = ['CSHQ TobD', 'CSHQ TobH', 'CSHQ JenD','CSHQ JenH',
          'Calcite', 'Calcium', 'Carbon', 'Silicium',
          'Average pH', 'Porosity']
comp =  ['CSHQ_TobD', 'CSHQ_TobH', 'CSHQ_JenD','CSHQ_JenH', 
         'calcite', 'Ca', 'C', 'Si', 'pH', 'avg_poros']
suffix = ['_csh_tobd', '_csh_tobh', '_csh_jend', '_csh_jenh', 
          '_calcite', '_calcium', '_carbon', '_silicium',
          '_average ph', 'avg_poros']
ylabel = [r'CSHQ_TobD $\cdot 10^{-15}$ (mol)', 
          r'CSHQ_TobH $\cdot 10^{-15}$ (mol)', 
          r'CSHQ_JenD $\cdot 10^{-15}$ (mol)', 
          r'CSHQ_JenH $\cdot 10^{-15}$ (mol)', 
          r'Calcite  $\cdot 10^{-15}$ mol', 
          r'Dissolved Ca $\cdot 10^{-15}$ (mol)', 
          r'Dissolved C $\cdot 10^{-15}$ (mol)' ,
          r'Dissolved Si $\cdot 10^{-15}$ (mol)' ,
          'Average pH',  'Porosity']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(np.array(results[names[i]]['time'])[s:],
                 results[names[i]][comp[k]][s:],
                 ls=linetype[i], label = label[i])
    plt.ylabel(ylabel[k], fontsize = 14)
    plt.xlabel('Time (d)', fontsize = 14)
    plt.legend(fontsize = 12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 

#np.savetxt('pco2.txt', np.c_(np.array(results[names[0]]['time']), results[names[0]]['portlandite']) )

#%% POROSITY DIFFERENCE
pt = []
dt = results[names[1]]['time'][2] - results[names[1]]['time'][1] 
t1 = 1
#print(get_porosity_val(results[names[1]], t1, dt))
text2 = ''
for i in range(0, len(names)): 
    nn = names[i]
    p = cf.get_porosity_val(results[nn], t1, dt)
    pt.append(p)
    text2 += ' \nPorosity at time ' + str(t1) + ' (р) is ' + \
              str(p) + ' for ' + str(label[i]) + ' CO2'
plt.figure(figsize=(8,4))
plt.plot(label, pt)
plt.title('Porosity at time %s hour' %t1)
plt.xlabel('CO2')
plt.ylabel('Porosity')
plt.tick_params(axis='both', which='major', labelsize=12)
plt.savefig(fpath + fname + '_poros_profile')
plt.show

  

#%% COMPARE POINT for cases
f = 'Ca (1, 3)'
plt.figure(figsize=(8,4))  
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]][f],
             ls=linetype[i], label = label[i])
plt.xlabel('Time (d)')
#plt.yscale("log")
plt.ylabel(f)
plt.legend()
plt.show() 
#%% Ca
comp = ['Ca (1, %s)'%i for i in range(3,8)]
title = ['Ca in %s'%i for i in range(3,8)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (d)')
    plt.yscale("log")
    plt.legend()
    plt.show() 
    
#%% Points
    
f = ('Ca', 'Ca') # ('Ca', 'Ca') ('Volume CH', 'vol_CH') ('Volume CC', 'vol_CC') ('pH', 'pH') ('De', 'De')('Porosity', 'poros')
title = ['%s in %s'%(f[0],i) for i in range(1,9)]
comp = ['%s (1, %s)'%(f[1],i) for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (d)')
    plt.title(label[i])
    plt.legend()
    plt.show() 


#%% pH
      
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
plt.figure(figsize = (8,4))
for i in range(0, len(names)): 
    pH = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + names[i] + '\\' + 'pH.npy')
    plt.plot(pH[1,1:-2], label=label[i], ls=linetype[i], linewidth = linewidths[i],
             color = 'tab:blue' )
plt.xlabel(r'Distance ($\mu$m)', fontsize = 14)
plt.ylabel('pH', fontsize = 14)
plt.legend(loc= "right", fontsize = 12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()

#%% C and Ca profiles together
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
C= {}
Ca = {}
for nn in names:
    C[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'C.npy')
    Ca[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'Ca.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
color = 'tab:orange'
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('Carbon (mol/l)', color=color,fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, C[names[i]][1,1:-2], 
             color=color, ls=linetype[i], linewidth = linewidths[i])

ax1.tick_params(axis='y', labelcolor=color, which='major', labelsize=12)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Calcium (mol/l)', color=color,fontsize = 14)  # we already handled the x-label with ax1

for i in range(0, len(names)):
    ax2.plot(x, Ca[names[i]][1,1:-2], 
             color=color, ls=linetype[i], linewidth = linewidths[i],
             label = label[i])
ax2.tick_params(axis='y', labelcolor=color, which='major', labelsize=12)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.legend(loc = "center right")
#plt.yscale("log")
plt.show()

#%% CC profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
CC= {}
for nn in names:
    CC[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'CC.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('Calcite',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, CC[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center right")
plt.show()

#%% CSH TobD profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
CSHQ_TobD= {}
for nn in names:
    CSHQ_TobD[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'CSHQ_TobD.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('CSHQ TobD',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, CSHQ_TobD[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center left")
plt.show()

#%% CSH TobH profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
CSHQ_TobH= {}
for nn in names:
    CSHQ_TobH[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'CSHQ_TobH.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('CSHQ TobH',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, CSHQ_TobH[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center left")
plt.show()


#%% CSH JenD profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
CSHQ_JenD= {}
for nn in names:
    CSHQ_JenD[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'CSHQ_JenD.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('CSHQ JenD',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, CSHQ_JenD[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center left")
plt.show()

#%% CSH JenH profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
CSHQ_JenH= {}
for nn in names:
    CSHQ_JenH[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'CSHQ_JenH.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('CSHQ JenH',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, CSHQ_JenH[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center left")
plt.show()


#%% Ca profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
Ca= {}
for nn in names:
    Ca[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'Ca.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('Ca',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, Ca[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center left")
plt.show()
#%% Si profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
Si= {}
for nn in names:
    Si[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'Si.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('Si',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, Si[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center right")
plt.show()
#%% C profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
C= {}
for nn in names:
    C[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'C.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('C',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, C[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center left")
plt.show()
#%% pH profile
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
linewidths = np.array([1.5, 1, 1, 1, 1])
pH= {}
for nn in names:
    pH[nn] = np.load(root_dir+'\\results\\output_csh\\02_24h\\' + nn + '\\' + 'pH.npy')
    
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('pH',fontsize = 14)
for i in range(0, len(names)):
    ax1.plot(x, pH[names[i]][1,1:-2], 
             label=label[i], ls=linetype[i], linewidth = linewidths[i])
plt.legend(loc = "center left")
plt.show()

#%% CSH in point 5

csh_TobD = 'CSHQ_TobD (1, 5)'
csh_TobH = 'CSHQ_TobH (1, 5)'
csh_JenD = 'CSHQ_JenD (1, 5)'
csh_JenH = 'CSHQ_JenH (1, 5)'
plt.figure(figsize=(8,4))  
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]][csh_TobD],
             ls=linetype[i], label = label[i])
    plt.plot(results[names[i]]['time'], results[names[i]][csh_TobH],
             ls=linetype[i], label = label[i])
    plt.plot(results[names[i]]['time'], results[names[i]][csh_JenD],
             ls=linetype[i], label = label[i])
    plt.plot(results[names[i]]['time'], results[names[i]][csh_JenH],
             ls=linetype[i], label = label[i])
plt.xlabel('Time (d)')
plt.ylabel('CSH')
plt.legend()
plt.show() 
    