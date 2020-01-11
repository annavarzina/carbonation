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
names = np.array(['02_pco2_1_p005', '03_pco2_2_p005', 
                  '04_pco2_3_p005', '05_pco2_34_p005'])
names = np.array(['02_pco2_1_p005', '03_pco2_2_p005', 
                  '04_pco2_3_p005', '05_pco2_34_p005'])
names = np.array(['01_p005_c34', '02_p005_c3', '03_p005_c252', '04_p005_c2', '05_p005_c1'])
label = np.array(['0.04%', '0.1%', '0.3%', '1%', '10%'])
names = np.array(['01_p005_c34m_ll1', '03_p005_c252m_ll1', '06_p005_c152m_ll1', '05_p005_c1m_ll1'])
label = np.array(['0.04%', '0.3%', '3%', '10%'])
names = np.array(['01_p005_c34m','03_p005_c252m', '04_p005_c2m', '05_p005_c1m'])
label = np.array(['0.04%', '0.3%', '1%','10%'])
linetype = np.array(['-', '--', '-.', ':', '-'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\06_pco2\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')
#%% SCALE
    
scale = 100
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
    

#%% CH DISSOLUTION 
#label = np.array(['0.03%', '10%', '1%', '0.1%', '0.01%'])) 
cf.plot_results(results, names, 'time', 'portlandite', label, linetype,
             'Portlandite', 'Time (h)', 'Portlandite*1e-12 (mol)', fpath, fname, '_CH',
             fsize = (8,4))
'''
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['portlandite'], 
             ls=linetype[i], label = label[i])      
plt.title('Portlandite')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_CH')
plt.show() 
'''

#%% CC PRECIPITATION
cf.plot_results(results, names, 'time', 'calcite', label, linetype,
             'Calcite', 'Time (h)', 'Calcite*1e-12 (mol)', fpath, fname, '_CC',
             fsize = (8,4))

#%% FULL PORTLANDITE TIME DISSOLUTION    
text1 = ''
for i in range(0, len(names)): 
    nn = names[i]
    #print('\nPortlandite dissolved in %s for %s' %(get_CH_dissolution_time(results[nn], Ts), nn))
    text1 += ' \nPortlandite dissolved in ' + \
             cf.get_CH_dissolution_time(results[nn], Ts) + '. for ' + label[i] + ' CO2'
             
#%% LIQUID CA 
cf.plot_results(results, names, 'time', 'Ca', label, linetype,
             'Ca in liquid', 'Time (h)', 'Ca*1e-12 (mol)', fpath, fname, '_Ca',
             fsize = (8,4))

#%% LIQUID C
cf.plot_results(results, names, 'time', 'C', label, linetype,
             'C in liquid', 'Time (h)', 'C*1e-12 (mol)', fpath, fname, '_C',
             fsize = (8,4))

#%% AVG PH
cf.plot_results(results, names, 'time', 'pH', label, linetype,
             'Average pH', 'Time (h)', 'pH', fpath, fname, '_pH',
             fsize = (8,4))

#%% INPUT C
cf.plot_results(results, names, 'time', 'C (1, 0)', label, linetype,
             'Input C', 'Time (h)', 'C (mol/l)', fpath, fname, '_inC',
             fsize = (8,4))

#%% POROSITY
cf.plot_results(results, names, 'time', 'avg_poros', label, linetype,
             'Porosity', 'Time (h)', 'Porosity', fpath, fname, '_poros',
             fsize = (8,4))

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
plt.savefig(fpath + fname + '_poros_profile')
plt.show
#%%
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], np.array(results[names[i]]['portlandite']) -
            results[names[i]]['portlandite'][0], label = label[i])
plt.xlabel('Time')
plt.ylabel(r'Change in Portlandite ($\cdot 10^{-12}$ mol/l)')
plt.legend()
plt.show


plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], np.array(results[names[i]]['calcite']), label = label[i])
plt.xlabel('Time')
plt.ylabel(r'Change in Portlandite ($\cdot 10^{-12}$ mol/l)')
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
#%% POROSITY, DE, PH IN CALCITE
#print(np.sort(results[names[1]].keys()))
#print('Porosity in node (1,1)')
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

#%% SAVE TXT
text = ''
for i in range(1,6):
    text += eval('text'+str(i))
np.save(fpath +fname, text)

#%%
mm_CH = 74
mm_CC = 100
mass0 = results[names[0]]['portlandite'][0]*mm_CH + results[names[0]]['calcite'][0]*mm_CC
mass_f = results[names[0]]['portlandite'][-1]*mm_CH + results[names[0]]['calcite'][-1]*mm_CC

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], 
             ((np.array(results[names[i]]['portlandite'])*mm_CH + \
              np.array(results[names[i]]['calcite'])*mm_CC)- \
              results[names[i]]['portlandite'][0]*mm_CH)*1e-15,
             label = label[i], ls = linetype[i])
plt.title('Degree of carbonation')
plt.xlabel('Time (h)')
plt.ylabel('Mass increase (g)')
plt.legend()
plt.show

#%% ADDITIONAL

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][0:10], results[names[i]]['Ca (1, 1)'][0:10],
             ls=linetype[i], label = label[i])
plt.title('Ca (1, 1)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 
#%%
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['C (1, 0)'],
             ls=linetype[i], label = label[i])
plt.title('C boundary')
plt.xlabel('Time (s)')
plt.legend()
plt.show()
#%%
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['vol_CH (1, 7)'],
             ls=linetype[i], label = label[i])
plt.title('CH ')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

#%% CC
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['vol_CC (1, 7)'],
             ls=linetype[i], label = label[i])
plt.title('CC ')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 
#%%
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][0:50], results[names[i]]['vol_CC (1, 2)'][0:50],
             ls=linetype[i], label = label[i])
plt.title('CC (1, 2)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 
#%%
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][5:100], results[names[i]]['pH (1, 0)'][5:100],
             ls=linetype[i], label = label[i])
#plt.ylim(11,12.2)
plt.title('pH (1, 0)')
plt.xlabel('Time (s)')
plt.legend()
plt.show() 
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][5:100], results[names[i]]['pH (1, 1)'][5:100],
             ls=linetype[i], label = label[i])
plt.ylim(11,12.2)
plt.title('pH (1, 2)')
plt.xlabel('Time (s)')
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
#%% Ca
comp = ['Ca (1, %s)'%i for i in range(1,9)]
title = ['Ca in %s'%i for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
    plt.yscale("log")
    plt.title(label[i])
    plt.legend()
    plt.show() 
#%% CC
comp = ['vol_CC (1, %s)'%i for i in range(1,9)]
title = ['CC in %s'%i for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
    #plt.yscale("log")
    plt.title(label[i])
    plt.legend()
    plt.show()  
    #%% CH
comp = ['vol_CH (1, %s)'%i for i in range(1,9)]
title = ['CH in %s'%i for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
    plt.yscale("log")
    plt.title(label[i])
    plt.legend()
    plt.show()  
#%%
      
plt.figure(figsize = (8,4))
for i in range(0, len(names)): 
    pH = np.load(root_dir+'\\results\\output\\06_pco2\\' + names[i] + '\\' + 'pH.npy')
    plt.plot(pH[1,1:-2], label=label[i])
plt.xlabel(r'Length ($\mu$m)')
plt.ylabel('pH')
plt.legend()
plt.show()
    