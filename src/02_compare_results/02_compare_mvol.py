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
fname = 'compare_mvol'
fpath = root_dir+'\\results\\output\\02_molar_volume\\'
fn.make_output_dir(fpath)

#names = np.array(['01_mvol1_subgrid', '02_mvol5_subgrid', '03_mvol10_subgrid',  '04_mvol50_subgrid',  '05_mvol100_subgrid',  ])
#names = np.array(['01_mvol1_subgrid_constbc', '02_mvol5_subgrid_constbc',  '03_mvol10_subgrid _constbc',  '04_mvol50_subgrid_constbc', '05_mvol100_subgrid_constbc', ])
#label = np.array([ '1', '5', '10', '50', '100'])
#names = np.array(['02_mvol1_fr04_ll5_Ca0', '02_mvol10_fr04_ll5_Ca0', '02_mvol50_fr04_ll5_Ca0','02_mvol100_fr04_ll5_Ca0'])
#names = np.array(['03_mvol1_fr04_ll5', '03_mvol10_fr04_ll5', '03_mvol50_fr04_ll5','03_mvol100_fr04_ll5'])
#names = np.array(['04_mvol1_fr04_ll1', '04_mvol10_fr04_ll1', '04_mvol50_fr04_ll1','04_mvol100_fr04_ll1'])
#names = np.array([ '05_mvol50_fr04_ll5_Di','05_mvol100_fr04_ll5_Di'])
#names = np.array(['06_mvol1_fr04_ll1_a', '06_mvol50_fr04_ll1_a','06_mvol100_fr04_ll1_a'])

names = np.array(['01_mvol1_mlvl_ll5', '01_mvol10_mlvl_ll5', '01_mvol50_mlvl_ll5','01_mvol100_mlvl_ll5'])
label = np.array([ '1', '10','50', '100'])
scale = [1.,10., 50., 100.]

names = np.array(['05_mvol1_fr04_ll5_Di', '05_mvol2_fr04_ll5_Di', '05_mvol10_fr04_ll5_Di',
                  '05_mvol50_fr04_ll5_Di','05_mvol100_fr04_ll5_Di'])
label = np.array([ '1','2', '10','50', '100'])
scale = [1.,2.,10., 50., 100.]


linetype = np.array(['-', '--', '-.', ':', '-'])
results = {}
for nn in names:
    path = root_dir+'\\results\\output\\02_molar_volume\\' + nn + '\\'
    results[nn] = fn.load_obj(path + nn +'_results')

#%% SCALE
#dtime = [350,350,350,350,200,0]
keys = ['portlandite', 'calcite', 'Ca', 'C', 'pH', 'time',
        'C (1, 0)', 'C (1, 1)', 'avg_poros', 
        'portlandite_cells', 'calcite_cells']
sres = {}
for i in range(0, len(names)):
    temp = {}
    a =np.array(results[names[i]]['time']) <= Ts/scale[i]
    s = np.size(results[names[i]]['time'])
    for k in keys:
        if(np.size(results[names[i]][k])==s):
            temp[k] = np.array(results[names[i]][k])[a]
    temp['time'] *= scale[i] 
    #temp['time']+=dtime[i]
    temp['portlandite'] *=scale[i]
    temp['calcite'] *=scale[i]
    sres[names[i]] = temp
#%% MAIN PROPERTIES
r1 = 1000
r2 = 5500
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
    for i in [0,3,4]:
        plt.plot(sres[names[i]]['time'][:r2], sres[names[i]][comp[k]][:r2],
                 ls=linetype[i], label = label[i])
    plt.ylabel(ylabel[k])
    #plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show() 

#%% DISSOLUTION AND PRECIPITATION RATE

titles = ['Dissolution rate', 'Precipitation rate' ]
comp =  ['portlandite', 'calcite']
suffix = ['_CH_rate', '_CC_rate' ]
s = 1
for k in range(0, len(comp)):
    plt.figure(figsize=(8,4))
    for i in range(0, len(names)):
        rate = np.abs(cf.get_rate(sres[names[i]][comp[k]],
                             sres[names[i]]['time'][2] - sres[names[i]]['time'][1],
                             step = s ))
        t = sres[names[i]]['time']
        r = rate
        l = len(sres[names[i]]['time'])- len(rate)
        if (l>0):
            t = sres[names[i]]['time'][l:]
            r = rate
        elif (l< 0):
            t = sres[names[i]]['time']
            r = rate[np.abs(l):] 
           
        plt.plot(t[10:r2], 
             r[10:r2],
             ls=linetype[i], label = label[i])
    plt.title(titles[k])
    plt.xlabel('Time (s)')
    plt.ylabel('Rate (mol/s)')
    plt.yscale("log")
    plt.legend()
    plt.savefig(fpath + fname + suffix[k])
    plt.show()
#plt.savefig(fpath + fname + '_CH_rate')

#%% INTERPOLATION CH
from scipy.interpolate import interp1d
name = names[0]#'02_mvol5_subgrid'
t = {}
ch = {}
cc = {}
for n in names:
    t[n] = sres[n]['time']
    ch[n] = sres[n]['portlandite']
    cc[n] = sres[n]['calcite']
#interpolation function for name[0]
ich = interp1d(t[name], ch[name], kind = 'cubic', fill_value="extrapolate")
icc = interp1d(t[name], cc[name], kind = 'cubic', fill_value="extrapolate")

diff_ch = {}
diff_cc = {}
for i in range(1, len(names)):
    diff_ch[names[i]] = np.abs(ch[names[i]] - ich(t[names[i]]))
    diff_cc[names[i]] = np.abs(cc[names[i]] - icc(t[names[i]]))

#%%
for i in range(1, len(names)):
    print(names[i])
    print(np.mean(diff_ch[names[i]][:r2]))
    print(np.max(diff_ch[names[i]][:r2]))
    print(np.mean(diff_cc[names[i]][:r2]))
    print(np.max(diff_cc[names[i]][:r2]))
    
'''   
for i in range(1, len(names)):
    print(names[i])
    print(np.mean(diff_ch[names[i]]))
    print(np.max(diff_ch[names[i]]))
    print(np.mean(diff_cc[names[i]]))
    print(np.max(diff_cc[names[i]]))
'''
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
for i in range(1, len(names)):
    #print('Scale %s' %scale[i])
    s.append(scale[i])
    npnorm.append(np.linalg.norm(diff_ch[names[i]], ord = 1))
    l2norm.append(l2_norm(ch[names[i]], ich(t[names[i]])))

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

#%% PERCENT ERROR CALCITE
p = -7000
print(sres[names[0]]['time'][p])
k = 'calcite'# 'avg_poros' #'portlandite''calcite'
er = np.array([])
di = np.array([])
for nn in names:
    d_ch = sres[nn][k][0]- sres[nn][k][p]
    er = np.append(er, d_ch)
    di = np.append(di,  results[nn][k][-1])
print(di)
for i in range(0,len(er)):
    e = (er[i]-er[0])/er[0] * 100    
    print('%s, %s ' %(names[i], str(e)))
#%% PERCENT ERROR PORTLANDITE
p = -7000
k = 'portlandite'# 'avg_poros' #'portlandite''calcite'
er = np.array([])
di = np.array([])
for nn in names:
    d_ch = sres[nn][k][0]- sres[nn][k][p]
    er = np.append(er, d_ch)
    di = np.append(di,  (results[nn][k][0]-results[nn][k][-1])/results[nn][k][0])
print(di)
for i in range(0,len(er)):
    e = (er[i]-er[0])/er[0] * 100    
    print('%s, %s ' %(names[i], str(e)))   
#%% AVERAGE PERCENT PORTLANDITE
k = 'portlandite'
print(sres[names[0]]['time'][r2])
#i = 3

dch =  np.abs(sres[names[0]][k][0]- sres[names[0]][k])
for i in range(1, len(names)):  
    if len(dch[:r2]) < len(diff_ch[names[i]][:r2]):
        e = diff_ch[names[i]][r1:r2]/dch[r1-1:r2]*100
    elif len(dch[:r2]) > len(diff_ch[names[i]][:r2]):
        e = diff_ch[names[i]][r1-1:r2]/dch[r1:r2]*100
    else:
        e = diff_ch[names[i]][r1:r2]/dch[r1:r2]*100
    print(names[i])
    print(np.mean(e))
    print(np.max(e))
    print(e[-1])

#%% AVERAGE PERCENT PORTLANDITE
k = 'calcite'
dcc =  np.abs(sres[names[0]][k][0]- sres[names[0]][k])

for i in range(1, len(names)):  
    if len(dcc[:r2]) < len(diff_cc[names[i]][:r2]):
        e = diff_cc[names[i]][r1:r2]/dcc[r1-1:r2]*100
    elif len(dch[:r2]) > len(diff_ch[names[i][:r2]]):
        e = diff_cc[names[i]][r1-1:r2]/dcc[r1:r2]*100
    else:
        e = diff_cc[names[i]][r1:r2]/dcc[r1:r2]*100
    print(names[i])
    print(np.mean(e))
    print(np.max(e))
    print(e[-1])


#%% TIME WHEN TRANSITION HAPPENS
t = []
for nn in names:
    ch24 = np.array(sres[nn]['portlandite_cells'])==24
    ch23 = np.array(sres[nn]['portlandite_cells'])==23
    cross = np.roll(ch24,-1) + ch23
    t.append(np.array(sres[nn]['time'])[~cross][0])
print(t)
    
#%% Ca
'''
comp = ['Ca (1, %s)'%i for i in range(1,9)]
title = ['Ca in %s'%i for i in range(1,9)]

for i in range(0, len(names)):
    plt.figure(figsize=(8,4))    
    for k in range(0, len(comp)):
        plt.plot(results[names[i]]['time'], results[names[i]][comp[k]],
                 ls=linetype[i], label = title[k])
    plt.xlabel('Time (s)')
    plt.title(label[i])
    plt.legend()
    plt.show() 
'''