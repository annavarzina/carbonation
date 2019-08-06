# -*- coding: utf-8 -*-
'''

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
import cell_type as ct # change the path to cell_type file
import misc_func as fn
#%% SETTINGS
Ts =100.
fname = 'pco2'
fpath = root_dir+'\\results\\output\\simulations\\compare\\'
fn.make_output_dir(fpath)
names = np.array([ '02_pco2_0', '02_pco2_1', 
                  '02_pco2_2', '02_pco2_3', '01_reference'])

label = np.array(['10%', '1%', '0.1%', '0.01%', '0.03%'])
linetype = np.array(['-', '--', '-.', ':', '-'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\simulations\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

#%% CH DISSOLUTION 
#label = np.array(['0.03%', '10%', '1%', '0.1%', '0.01%']))
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['portlandite'], 
             ls=linetype[i], label = label[i])      
plt.title('Portlandite')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_CH')
plt.show() 
#%% CC PRECIPITATION
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['calcite'], 
             ls=linetype[i], label = label[i])      
plt.title('Calcite')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_CC')
plt.show() 

#%% FULL PORTLANDITE TIME DISSOLUTION
def get_CH_dissolution_time(res, T):
    c = np.array(res['portlandite']) <=0
    if np.any(c):
        t = np.array(results[names[1]]['time'])[c]
        return str(np.min(t)) + ' s'
    else:
        return '>' + str(T) + ' s' 
    
text1 = ''
for i in range(0, len(names)): 
    nn = names[i]
    #print('\nPortlandite dissolved in %s for %s' %(get_CH_dissolution_time(results[nn], Ts), nn))
    text1 += ' \nPortlandite dissolved in ' + \
             get_CH_dissolution_time(results[nn], Ts) + '. for ' + label[i] + ' CO2'
#%% LIQUID CA 
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['Ca'], 
             ls=linetype[i], label = label[i])      
plt.title('Ca in liquid')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_Ca')
plt.show()   

#%% LIQUID C
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['C'], 
             ls=linetype[i], label = label[i])      
plt.title('C in liquid')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_C')
plt.show() 

#%% TOTAL PH
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['pH'], 
             ls=linetype[i], label = label[i])      
plt.title('Average pH')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_pH')
plt.show() 

#%% INPUT C

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['C (1, 0)'], 
             ls=linetype[i], label = label[i])      
plt.title('Input C')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_inC')
plt.show() 

#%% POROSITY

plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'], results[names[i]]['avg_poros'], 
             ls=linetype[i], label = label[i])      
plt.title('Porosity')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig(fpath + fname + '_poros')
plt.show() 

#%% POROSITY DIFFERENCE
pt = []
def get_porosity_val(res, t, dt):
    
    c = np.logical_and(np.array(res['time']) >= t, 
                   np.array(res['time']) < t+dt)
    if np.any(c):
        p = np.array(res['avg_poros'])[c]
        return p[0]
    else:
        return 'Error' #raise error
dt = results[names[1]]['time'][2] - results[names[1]]['time'][1] 
t1 = 30
#print(get_porosity_val(results[names[1]], t1, dt))
text2 = ''
for i in range(0, len(names)): 
    nn = names[i]
    p = get_porosity_val(results[nn], t1, dt)
    pt.append(p)
    text2 += ' \nPorosity at time ' + str(t1) + ' is ' + \
              str(p) + ' for ' + str(label[i]) + ' CO2'
plt.figure(figsize=(8,4))
plt.plot(label, pt)
plt.title('Porosity at time %s' %t1)
plt.xlabel('CO2')
plt.ylabel('Porosity')
plt.savefig(fpath + fname + '_poros_profile')
plt.show

#%% DISSOLUTION RATE
def get_rate(res, dt):
    ch = np.array(res)
    dCH = np.zeros(ch.shape,np.float)
    dCH[0:-1] = np.diff(ch)/dt
    return dCH
#dCH[0] = 0
t2 = 50
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][t2:-1], 
             get_rate(results[names[i]]['portlandite'][t2:-1], dt),
             label = label[i], ls = linetype[i])
plt.show

#%% PRECIPITATION RATE
t2 = 50
plt.figure(figsize=(8,4))
for i in range(0, len(names)):
    plt.plot(results[names[i]]['time'][t2:-1], 
             get_rate(results[names[i]]['calcite'][t2:-1], dt),
             label = label[i], ls = linetype[i])
plt.show
#%% SI, POROSITY, DE, PH IN CALCITE