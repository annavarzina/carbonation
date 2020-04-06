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
#import plot_func as cf
#%% SETTINGS
Ts =1000.
fname = 'compare_mvol'
fpath = root_dir+'\\results\\output\\02_compare_mvol\\'
fn.make_output_dir(fpath)


names = np.array(['02_mvol-1_leach', '02_mvol-50_leach', '02_mvol-100_leach'])
label = np.array([ '1','50', '100'])
scale = [1.,50., 100.]

#names = np.array(['06_mvol1_fr04_CH_leach', '06_mvol50_fr04_CH_leach', '06_mvol100_fr04_CH_leach'])
#label = np.array([ '1','50', '100'])
#scale = [1.,50., 100.]

linetype = np.array(['-', '--', '-.', ':', '-'])
ch = {}
time = {} #already scaled
dtime = {}
ca = {}
for nn in names:
    path = root_dir+'\\results\\output\\02_compare_mvol\\' + nn + '\\'
    ch[nn] = np.load(path +'CH.npy')
    time[nn] = np.load(path +'time.npy')
    dtime[nn] = np.load(path +'dis_time.npy')
    ca[nn] = np.load(path +'Ca_prof.npy')


#%% CH

plt.figure(figsize=(8,4)) 
for i in range(0, len(names)):
    plt.plot(time[names[i]], 
             ch[names[i]]*scale[i],
             ls=linetype[i], label = label[i])
plt.xlabel('Time (s)')
plt.legend()
plt.show() 

#%% CH difference
for i in range(0, len(names)):
    print("dCH %s" %((ch[names[i]][0]-ch[names[i]][-1])*scale[i]))