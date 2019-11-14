# -*- coding: utf-8 -*-
'''
Compare the results for different time step
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
Ts =1000.
scale = 50
fname = 'compare_dt_p005'
fpath = root_dir+'\\results\\output\\01_time_step\\'
fn.make_output_dir(fpath)
#label = np.array(['f 1', 'f 2', 'f 4', 'f 8'])
linetype = np.array(['-', '--', '-.', ':', '-', '--', '-.', ':', ])
colors = np.array([])
names = np.array(['01_dt1_p005_subgrid', '02_dt2_p005_subgrid', 
                  '03_dt5_p005_subgrid', '04_dt10_p005_subgrid'])
#names = np.array(['05_dt1_p005_Dborder', '06_dt2_p005_Dborder',
#                  '07_dt5_p005_Dborder', '08_dt10_p005_Dborder'])
ts = np.array([8.334e-06, 1.667e-05, 4.167e-05, 8.334e-05,
               8.334e-06, 1.667e-05, 4.167e-05, 8.334e-05])
results = {}
for nn in names:
    path = root_dir+'\\results\\output\\01_time_step\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')
    
for i in range(0, len(names)):
    temp = np.array(results[names[i]]['time'])
    temp *= scale
    results[names[i]]['time']= temp.tolist()
    
label = names
#%% CH DISSOLUTION 
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
          'Average pH', 'Input C', 'Porosity']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'C (1, 0)', 'avg_poros']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_average ph', '_input_c', '_poros']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, 4):
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
limits = [[-0.1, 0], [0, 0.1]]
rate = {}
for i in range(0, len(names)):
    rate[names[i]] = {}
    for k in range(0, len(comp)):
        rate[names[i]][comp[k]] = cf.get_rate(results[names[i]][comp[k]],
                             results[names[i]]['time'][2] - results[names[i]]['time'][1])

rstart = 0
rend = len(results[names[1]]['time']) - 1      
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'][rstart:rend], 
                 rate[names[i]][comp[k]][rstart:rend],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    #plt.ylim(limits[k])
    plt.xlabel('Time (s)')
    plt.ylabel('Rate (mol/s)')
    plt.legend()
    #plt.savefig(fpath + fname + suffix[k])
    plt.show()
#plt.savefig(fpath + fname + '_CH_rate')
    

#%% INTERPOLATION
from scipy.interpolate import interp1d

t = {}
p = {}
for n in names:
    t[n] = results[n]['time']
    p[n] = results[n]['portlandite']#['avg_poros']
f = interp1d(t[names[3]], p[names[3]], kind = 'cubic', fill_value="extrapolate")
pi = {}
pi = {}
diff = {}
for i in range(0, len(names)):
    pi[names[i]] = f(t[names[i]])
    diff[names[i]] = np.abs(p[names[i]] - pi[names[i]])

#%% ERROR NORM

def l2_norm(v, o):
    '''
    v - vector
    o - original solution vector
    '''
    return np.sqrt(np.sum((v - o)**2)/ np.sum(o**2))

s = []
npnorm = []
l2norm = []
for i in range(0, len(names)):
    #print('Scale %s' %scale[i])
    s.append(ts[i])
    npnorm.append(np.linalg.norm(diff[names[i]], ord = 1))
    l2norm.append(l2_norm(p[names[i]], pi[names[i]]))
#%% PLOT NORM
plt.figure(figsize=(8,4))
plt.plot(s, npnorm)
plt.title('L2 norm')
plt.xlabel('Time step')
plt.ylabel('Numpy L2 norm')
#plt.savefig(fpath + fname + '_np_l2_norm')
plt.show()

plt.figure(figsize=(8,4))
plt.plot(s, l2norm)
plt.title('Relative error')
plt.xlabel('Time step')
plt.ylabel('Relative error')
#plt.savefig(fpath + fname + '_rel_error')
plt.show()

#%% DIFFERENCE IN PERCENT
k = 'portlandite'
er = np.array([])
for nn in names:
    e = results[nn][k][0]- results[nn][k][-2]
    e = e/results[nn][k][0]*100
    er = np.append(er, e)
#%% 
for nn in names:
    ch24 = np.array(results[nn]['portlandite_cells'])==24
    ch23 = np.array(results[nn]['portlandite_cells'])==23
    cross = np.roll(ch24,-1) +ch23
    t = np.array(results[nn]['time'])[~cross]
    
#%%
    
titles = ['Volume CH in (1,1)', 'Volume CH in (1,2)', 'Volume CH in (1,3)']
comp =  ['vol_CH (1, 1)', 'vol_CH (1, 2)', 'vol_CH (1, 3)']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        #indices = np.where(np.logical_and(np.array(results[names[i]]['time'])<=5,
       #                    np.array(results[names[i]]['time'])>=3.1))
        indices = np.where(np.array(results[names[i]]['time'])<=2)
        plt.plot(np.array(results[names[i]]['time'])[indices],
                 np.array(results[names[i]][comp[k]])[indices],
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.show() 