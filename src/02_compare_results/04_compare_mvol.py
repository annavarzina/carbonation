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
Ts =1000.
fname = 'mvol'
fpath = root_dir+'\\results\\output\\simulations\\compare\\'
fn.make_output_dir(fpath)
#names = np.array(['05_mvol_40', '01_reference', '05_mvol_10', '05_mvol_2', '05_mvol_1'])
#label = np.array(['0.331*40', '0.331*20','0.331*10', '0.331*2', '0.331'])
#linetype = np.array(['-', '--', '-.', ':', '-'])
names = np.array(['05_mvol_500', '05_mvol_100', '05_mvol_50', 
                  '05_mvol_10', '05_mvol_5', '05_mvol_1'])
label = np.array(['500', '100','50', '10', '5', '1'])
linetype = np.array(['-', '--', '-.', ':', '-', '--'])

results = {}
for nn in names:
    path = root_dir+'\\results\\output\\simulations\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

#%% SCALE
scale = [500, 100,50, 10, 5, 1]
keys = ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'time', 'C (1, 0)']
sres = {}
for i in range(0, len(names)):
    temp = {}
    a =np.array(results[names[i]]['time']) <= Ts/scale[i]
    s = np.size(results[names[i]]['time'])
    for k in keys:
        if(np.size(results[names[i]][k])==s):
            temp[k] = np.array(results[names[i]][k])[a]
    temp['time'] *= scale[i]
    temp['portlandite'] *=scale[i]
    temp['calcite'] *=scale[i]
    sres[names[i]] = temp
#%% CH DISSOLUTION 
titles = ['Portlandite', 'Calcite', 'Calcium', 'Carbon',
          'Average pH', 'Input C']
comp =  ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'C (1, 0)']
suffix = ['_portlandite', '_calcite', '_calcium', '_carbon',
          '_average ph', '_input_c']
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(sres[names[i]]['time'], sres[names[i]][comp[k]],
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
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        plt.plot(sres[names[i]]['time'], 
                 cf.get_rate(sres[names[i]][comp[k]],
                             sres[names[i]]['time'][2] - sres[names[i]]['time'][1]),
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
    t[n] = sres[n]['time']
    p[n] = sres[n]['portlandite']
f = interp1d(t['05_mvol_1'], p['05_mvol_1'], kind = 'cubic')
pi = {}
diff = {}
for i in range(0, len(names)-1):
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
for i in range(0, len(names)-1):
    #print('Scale %s' %scale[i])
    s.append(scale[i])
    npnorm.append(np.linalg.norm(diff[names[i]], ord = 1))
    l2norm.append(l2_norm(p[names[i]], pi[names[i]]))
#%% PLOT NORM
plt.figure(figsize=(8,4))
plt.plot(s, npnorm)
plt.title('L2 norm')
plt.xlabel('Scale')
plt.ylabel('Numpy L2 norm')
plt.savefig(fpath + fname + '_np_l2_norm')
plt.show()

plt.figure(figsize=(8,4))
plt.plot(s, l2norm)
plt.title('Relative error')
plt.xlabel('Scale')
plt.ylabel('Relative error')
plt.savefig(fpath + fname + '_rel_error')
plt.show()

