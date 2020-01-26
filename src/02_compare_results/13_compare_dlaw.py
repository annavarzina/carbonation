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
fpath = root_dir+'\\results\\output\\14_compare_dlaw\\'
fn.make_output_dir(fpath)
names = np.array(['01_mvol100_fr04_ll5_Da', '01_mvol100_fr04_ll5_Dc', '01_mvol100_fr04_ll5_Di'])
label = np.array([ 'Archie','Const', 'Weighted'])


linetype = np.array(['-', '--', '-.', ':', '-'])
results = {}
for nn in names:
    path = root_dir+'\\results\\output\\14_dlaw\\' + nn + '\\'
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

#%% CC PRECIPITATION
cf.plot_results(results, names, 'time', 'calcite', label, linetype,
             'Calcite', 'Time (h)', 'Calcite*1e-12 (mol)', fpath, fname, '_CC',
             fsize = (8,4))
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

#%% De
cf.plot_results(results, names, 'time', 'avg_D_eff', label, linetype,
             'De', 'Time (h)', 'De', fpath, fname, '_De',
             fsize = (8,4))
