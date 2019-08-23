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
fname = 'dt'
fpath = root_dir+'\\results\\output\\compare_time_step\\'
fn.make_output_dir(fpath)
#names = np.array(['05_mvol_40', '01_reference', '05_mvol_10', '05_mvol_2', '05_mvol_1'])
#label = np.array(['0.331*40', '0.331*20','0.331*10', '0.331*2', '0.331'])
#linetype = np.array(['-', '--', '-.', ':', '-'])
names = np.array(['01_dt1_p05', '01_dt2_p05', '01_dt4_p05', '01_dt8_p05'])#, '01_fixD_005p_D11'])
label = np.array(['f 1', 'f 2', 'f 4', 'f 8'])
linetype = np.array(['-', '--', '-.', ':'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\time_step\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

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
rstart = 500
rend = len(results[names[1]]['time']) - 1

for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(results[names[i]]['time'][rstart:rend], 
                 cf.get_rate(results[names[i]][comp[k]][rstart:rend],
                             results[names[i]]['time'][2] - \
                             results[names[i]]['time'][1]),
                 ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.ylabel('Rate (mol/s)')
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show()
#plt.savefig(fpath + fname + '_CH_rate')

#%% INTERPOLATION
from scipy.interpolate import interp1d

t = {}
p = {}
for n in names:
    t[n] = results[n]['time']
    p[n] = results[n]['portlandite']#['avg_poros']
f = interp1d(t['01_dt8_p05'], p['01_dt8_p05'], kind = 'cubic', fill_value="extrapolate")
pi = {}
diff = {}
for i in range(0, len(names)):
    pi[names[i]] = f(t[names[i]])
    diff[names[i]] = np.abs(p[names[i]] - pi[names[i]])

#%% ERROR NORM
ts = [1., 2, 4, 8.]
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
plt.xlabel('Time step factor')
plt.ylabel('Numpy L2 norm')
plt.savefig(fpath + fname + '_np_l2_norm')
plt.show()

plt.figure(figsize=(8,4))
plt.plot(s, l2norm)
plt.title('Relative error')
plt.xlabel('Time step factor')
plt.ylabel('Relative error')
plt.savefig(fpath + fname + '_rel_error')
plt.show()

