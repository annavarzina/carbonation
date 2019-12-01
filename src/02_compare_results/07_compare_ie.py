# -*- coding: utf-8 -*-
'''
Compare the results for different IE
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
fname = 'ie'
fpath = root_dir+'\\results\\output\\07_internal_energy\\compare\\'
fn.make_output_dir(fpath)
#names = np.array([ '01_ie01_p05', '02_ie05_p05', '03_ie1_p05'])
names = np.array([ '01_ie01_p005_Db09', '02_ie05_p005_Db09', '03_ie1_p005_Db09'])

label = np.array(['ie 0.1','ie 0.5', 'ie 1'])
linetype = np.array(['-', '--', '-.',':'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\07_internal_energy\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

#%% SCALE
    
scale = 50
for i in range(0, len(names)):
    temp = np.array(results[names[i]]['time'])
    temp *= scale
    results[names[i]]['time']= temp.tolist()
#%%
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

#%% DISSOLUTION RATE

titles = ['Dissolution rate', 'Precipitation rate' ]
comp =  ['portlandite', 'calcite']
suffix = ['_CH_rate', '_CC_rate' ]
r = 500
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'][r:-1], 
                 cf.get_rate(results[names[i]][comp[k]],
                             results[names[i]]['time'][2] - results[names[i]]['time'][1])[r:-1],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.ylabel('Rate (mol/s)')
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
    plt.plot(results[names[i]]['time'][0:500], results[names[i]]['Ca (1, 1)'][0:500],
             ls=linetype[i], label = label[i])
plt.title('Ca (1, 1)')
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
'''
#%% CH DISSOLUTION 
#label = np.array(['0.03%', '10%', '1%', '0.1%', '0.01%'])) 
cf.plot_results(results, names, 'time', 'portlandite', label, linetype,
             'Portlandite', 'Time (s)', 'Portlandite', fpath, fname, '_CH',
             fsize = (8,4))


#%% CC PRECIPITATION
cf.plot_results(results, names, 'time', 'calcite', label, linetype,
             'Calcite', 'Time (s)', 'Calcite', fpath, fname, '_CC',
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
             'Ca in liquid', 'Time (s)', 'Ca*1e-12 (mol)', fpath, fname, '_Ca',
             fsize = (8,4))

#%% LIQUID C
cf.plot_results(results, names, 'time', 'C', label, linetype,
             'C in liquid', 'Time (s)', 'C*1e-12 (mol)', fpath, fname, '_C',
             fsize = (8,4))

#%% AVG PH
cf.plot_results(results, names, 'time', 'pH', label, linetype,
             'Average pH', 'Time (s)', 'pH', fpath, fname, '_pH',
             fsize = (8,4))

#%% INPUT C
cf.plot_results(results, names, 'time', 'C (1, 0)', label, linetype,
             'Input C', 'Time (s)', 'C (mol/l)', fpath, fname, '_inC',
             fsize = (8,4))

#%% POROSITY
cf.plot_results(results, names, 'time', 'avg_poros', label, linetype,
             'Porosity', 'Time (s)', 'Porosity', fpath, fname, '_poros',
             fsize = (8,4))

#%% POROSITY DIFFERENCE
pt = []
dt = results[names[1]]['time'][2] - results[names[1]]['time'][1] 
t1 = 5 * scale
#print(get_porosity_val(results[names[1]], t1, dt))
text2 = ''
for i in range(0, len(names)): 
    nn = names[i]
    p = cf.get_porosity_val(results[nn], t1, dt)
    pt.append(p)
    text2 += ' \nPorosity at time ' + str(t1) + ' is ' + \
              str(p) + ' for ' + str(label[i]) + ' CO2'
plt.figure(figsize=(8,4))
plt.plot(label, pt)
plt.title('Porosity at time %s' %t1)
plt.xlabel('IE')
plt.ylabel('Porosity')
plt.savefig(fpath + fname + '_poros_profile')
plt.show

#%% DISSOLUTION RATE
#dCH[0] = 0
t2 =10 *scale
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][t2:-1], 
             cf.get_rate(results[names[i]]['portlandite'][t2:-1], dt),
             label = label[i], ls = linetype[i])
plt.title('Dissolution rate')
plt.xlabel('Time (s)')
plt.ylabel('Rate (mol/s)')
plt.legend()
plt.savefig(fpath + fname + '_CH_rate')
plt.show

#%% PRECIPITATION RATE
t2 = 10 *scale
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][t2:-1], 
             cf.get_rate(results[names[i]]['calcite'][t2:-1], dt),
             label = label[i], ls = linetype[i])
plt.title('Precipitation rate')
plt.xlabel('Time (s)')
plt.ylabel('Rate (mol/s)')
plt.legend()
plt.savefig(fpath + fname + '_CC_rate')
plt.show
#%%
pt = []
#dt = results[names[1]]['time'][2] - results[names[1]]['time'][1] 
t1 = 15 * scale
for i in range(0, len(names)): 
    nn = names[i]
    p = cf.get_val_at_time(results[nn], 'poros (1, 1)', t1, dt)
    pt.append(p)
plt.figure(figsize=(8,4))
plt.plot(label, pt)
plt.title('Porosity in a single point at time %s' %t1)
plt.xlabel('IE')
plt.ylabel('Porosity')
plt.savefig(fpath + fname + '_poros_11')
plt.show


#%% POROSITY, DE, PH IN CALCITE
#print(np.sort(results[names[1]].keys()))
#print('Porosity in node (1,1)')
text3 = ''
for i in range(0, len(names)): 
    text3 += ' \nPorosity in a node (1,1) ' + \
             str(results[names[i]]['poros (1, 1)'][-1]) + \
             '. for ' + label[i] + ' '+ fname
text4 = ''
for i in range(0, len(names)): 
    text4 += ' \npH in a node (1,1) ' + \
             str(results[names[i]]['pH (1, 1)'][-1]) + \
             '. for ' + label[i] + ' '+fname
text5 = ''
for i in range(0, len(names)): 
    text5 += ' \nDe in a node (1,1) ' + \
             str(results[names[i]]['De (1, 1)'][-1]) + \
             '. for ' + label[i] + ' '+fname
             
#%% SI
si = {}
for i in range(0, len(names)): 
    nn = names[i]
    path = root_dir+'\\results\\output\\simulations\\' + nn + '\\'
    si[names[i]] = np.load(path + 'SI' + '.npy')

text6 = ''
for i in range(0, len(names)):   
    text6 += ' \nSI in a node (1,1) ' + \
             str(si[names[i]][(1,1)]) + \
             '. for ' + label[i] + ' '+fname
   
text7 = ''
for i in range(0, len(names)):   
    text7 += ' \nSI in a node (1,10) ' + \
             str(si[names[i]][(1,10)]) + \
             '. for ' + label[i] + ' '+fname
             
#%% SAVE TXT
text = ''
for i in range(1,8):
    text += eval('text'+str(i))
np.save(fpath +fname, text)
'''