# -*- coding: utf-8 -*-
'''
Compare the results for different PCO2
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
fname = 'pco2'
fpath = root_dir+'\\results\\output\\06_pco2\\compare\\'
fn.make_output_dir(fpath)
#names = np.array([ '02_pco2_0', '02_pco2_1', 
#                  '02_pco2_2', '02_pco2_3', '01_reference'])
#names = np.array(['02_pco2_1_p005', '03_pco2_2_p005','04_pco2_3_p005', '05_pco2_34_p005'])
#names = np.array(['02_pco2_1_p005', '03_pco2_2_p005', '04_pco2_3_p005', '05_pco2_34_p005'])
#names = np.array(['01_p005_c34m_ll1', '03_p005_c252m_ll1', '06_p005_c152m_ll1', '05_p005_c1m_ll1'])
#label = np.array(['0.04%', '0.3%', '3%', '10%'])
#names = np.array(['01_p005_c34', '03_p005_c252', '06_p005_c152', '05_p005_c1'])
scale = 100
label = np.array(['0.04%',  '0.3%',  '3%','10%'])
#names = np.array(['01_p005_c34m', '02_p005_c3m', '03_p005_c252m', '04_p005_c2m', '06_p005_c152m', '05_p005_c1m'])
names = np.array(['01_p005_c34m', '03_p005_c252m', '06_p005_c152m', '05_p005_c1m'])

label = np.array(['0.04%',  '0.3%',  '3%','10%'])

names = np.array(['08_p0_D13', '08_p052_D13', '08_p1_D13', '08_p152_D13', '08_p2_D13', '08_p34_D13'])

names = np.array(['08_p052_D13', '08_p1_D13', '08_p152_D13', '08_p252_D13', '08_p34_D13'])
label = np.array(['30%', '10%',  '3%',  '0.3%', '0.04%'])
names = np.array(['08_p1_D13', '08_p152_D13', '08_p252_D13', '08_p34_D13'])
label = np.array(['10%',  '3%',  '0.3%', '0.04%'])
linetype = np.array(['-', '--', '-.', ':', '-', '--'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\06_pco2\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')
#%% SCALE
    
for i in range(0, len(names)):
    temp = np.array(results[names[i]]['time'])
    temp *= scale/3600
    results[names[i]]['time']= temp.tolist()
    temp = np.array(results[names[i]]['portlandite'])
    temp *= scale
    results[names[i]]['portlandite']= temp.tolist()
    temp = np.array(results[names[i]]['calcite'])
    temp *= scale
    results[names[i]]['calcite']= temp.tolist()
    

#%% PARAMS
s = 5
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
          'Average pH', 'Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_average ph', 'avg_poros']
ylabel = [r'Portlandite $\cdot 10^{-15}$ (mol)', r'Calcite  $\cdot 10^{-15}$ mol', 
          r'Dissolved Ca $\cdot 10^{-15}$ (mol)', r'Dissolved C $\cdot 10^{-15}$ (mol)' ,
          'Average pH',  'Porosity']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(np.array(results[names[i]]['time'])[s:],
                 results[names[i]][comp[k]][s:],
                 ls=linetype[i], label = label[i])
    plt.ylabel(ylabel[k], fontsize = 14)
    plt.xlabel('Time (h)', fontsize = 14)
    plt.legend(fontsize = 12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 

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
#%%
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], np.array(results[names[i]]['portlandite']) -
            results[names[i]]['portlandite'][0], label = label[i])
plt.xlabel('Time')
plt.ylabel(r'Change in Portlandite ($\cdot 10^{-15}$ mol/l)')
plt.legend()
plt.show



plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], np.array(results[names[i]]['calcite']), label = label[i])
plt.xlabel('Time')
plt.ylabel(r'Change in Calcite ($\cdot 10^{-15}$ mol/l)')
plt.legend()
plt.show

#%% Change in porosity
s = 500
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    cp =  np.array(results[names[i]]['avg_poros']) -results[names[i]]['avg_poros'][s]
    print(cp[-1]*100)
    plt.plot(results[names[i]]['time'][s:],cp[s:], label = label[i])
plt.xlabel('Time')
plt.ylabel('Change in porosity')
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend()
plt.show
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
#%% TEXT  
text1 = ''
for i in range(0, len(names)): 
    nn = names[i]
    #print('\nPortlandite dissolved in %s for %s' %(get_CH_dissolution_time(results[nn], Ts), nn))
    text1 += ' \nPortlandite dissolved in ' + \
             cf.get_CH_dissolution_time(results[nn], Ts) + '. for ' + label[i] + ' CO2'
             
text3 = ''
for i in range(0, len(names)): 
    text3 += ' \nPorosity in a node (1,1) ' + \
             str(results[names[i]]['poros (1, 1)'][10000]) + \
             '. for ' + label[i] + ' CO2'
text4 = ''
for i in range(0, len(names)): 
    text4 += ' \npH in a node (1,1) ' + \
             str(results[names[i]]['pH (1, 1)'][10000]) + \
             '. for ' + label[i] + ' CO2'
text5 = ''
for i in range(0, len(names)): 
    text5 += ' \nDe in a node (1,1) ' + \
             str(results[names[i]]['De (1, 1)'][10000]) + \
             '. for ' + label[i] + ' CO2'

text = ''
for i in range(1,6):
    text += eval('text'+str(i))
np.save(fpath +fname, text)

#%% Degree of carbonation
mm_CH = 74
mm_CC = 100
mass0 = results[names[0]]['portlandite'][0]*mm_CH + results[names[0]]['calcite'][0]*mm_CC
mass_f = results[names[0]]['portlandite'][-1]*mm_CH + results[names[0]]['calcite'][-1]*mm_CC

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    d = ((np.array(results[names[i]]['portlandite'])*mm_CH + \
              np.array(results[names[i]]['calcite'])*mm_CC)- \
              results[names[i]]['portlandite'][0]*mm_CH)*1e-15
    print(d[-1])
    plt.plot(results[names[i]]['time'], d,
             label = label[i], ls = linetype[i])
#plt.title('Degree of carbonation')
plt.xlabel('Time (h)', fontsize = 14)
plt.ylabel('Mass increase (g)', fontsize = 14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend(fontsize = 12)
plt.show
#%% COMPARE POINT for cases
f = 'Ca (1, 1)'
plt.figure(figsize=(8,4))  
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]][f],
             ls=linetype[i], label = label[i])
plt.xlabel('Time (s)')
#plt.yscale("log")
plt.ylabel(f)
plt.legend()
plt.show() 
#%% De
comp = ['De (1, %s)'%i for i in range(1,7)]
title = ['De in %s'%i for i in range(1,7)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
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
    plt.xlabel('Time (s)')
    plt.title(label[i])
    plt.legend()
    plt.show() 
#%%
      
plt.figure(figsize = (8,4))
for i in range(0, len(names)): 
    pH = np.load(root_dir+'\\results\\output\\06_pco2\\' + names[i] + '\\' + 'pH.npy')
    plt.plot(pH[1,1:-2], label=label[i], ls=linetype[i])
plt.xlabel(r'Distance ($\mu$m)', fontsize = 14)
plt.ylabel('pH', fontsize = 14)
plt.legend(loc= "right", fontsize = 12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()
    #%%
      
plt.figure(figsize = (8,4))
for i in range(0, len(names)): 
    C = np.load(root_dir+'\\results\\output\\06_pco2\\' + names[i] + '\\' + 'C.npy')
    plt.plot(C[1,1:-2], label=label[i], ls=linetype[i])
plt.xlabel(r'Distance ($\mu$m)', fontsize = 14)
plt.ylabel('C (mol/l)', fontsize = 14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend(loc= "right", fontsize = 12)
plt.show()
#%%
      
plt.figure(figsize = (8,4))
for i in range(0, len(names)): 
    Ca = np.load(root_dir+'\\results\\output\\06_pco2\\' + names[i] + '\\' + 'Ca.npy')
    plt.plot(Ca[1,1:-2], label=label[i], ls=linetype[i])
plt.xlabel(r'Distance ($\mu$m)', fontsize = 14)
plt.ylabel('Ca (mol/l)', fontsize = 14)
plt.legend(loc= "right", fontsize = 12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()


#%% CC profile
#c = np.array(['#C22500','#FFEA5E','#9CDFFF','#0172D8'])
#c = np.array(['b','r','g','y'])
c = plt.rcParams['axes.prop_cycle'].by_key()['color']
cc_mass = 100.09
ch_mass = 74.09
cc = {}
ch = {}
for n in range(0, len(names)):
    cc[names[n]] = []
    ch[names[n]] = []
    for i in range(1,10):
        cc[names[n]].append(results[names[n]]['calcite (1, ' + str(i) +')' ][-1]*cc_mass)
        ch[names[n]].append(results[names[n]]['portlandite (1, ' + str(i) +')' ][-1]*ch_mass)
    chmax = results[names[n]]['portlandite (1, ' + str(i) +')' ][1]*ch_mass
    for i in range(10,30):
        cc[names[n]].append(0)
        ch[names[n]].append(chmax)
        
plt.figure(figsize=(8,4))
for n in range(0, len(names)):
    plt.plot(range(1,30),cc[names[n]][:], ls=linetype[n], color = c[n], label = label[n])
    plt.plot(range(1,30),ch[names[n]][:], ls=linetype[n], color = c[n])
plt.ylabel(r'Mineral mass  $\cdot 10^{-15}$ (g) in 1 $\mu m^3$ ')
plt.xlabel(r'Distance ($\mu m$)')
plt.legend()
