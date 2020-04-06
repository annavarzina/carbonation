# -*- coding: utf-8 -*-
'''
Compare the results for different D
'''
#%% MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import misc_func as fn

#%% UPLOAD
#name = '02_ccD12_long_prev'
name = '00_ccD11_prev'
path = root_dir+'\\results\\output\\13_validation\\' + name + '\\'
results = fn.load_obj(path + '00_ccD11' +'_results')
scale = 100
#%%

mm_CH = 74.093
mm_CC = 100.086
alpha = np.array([0.,0.3,0.55,1.05,1.25,1.55])
time = np.array([0.,2.,4.5,10.,12.,20.]) #days


mass0 = 22.977 *1e-3 #[g] initial mass of CH crystal
dmass = alpha*mass0/100 #delta mass increase
dmm  = mm_CC-mm_CH #[g/mol] molar mass difference
xmol = dmass/dmm #[mol] mol of CH transformed to CC

sa  = 29.19e+6
depth_sa = xmol*0.0331*1e+15/sa
mass_sa = dmass/sa
mol_sa = xmol/sa

print(mass_sa)
print(mol_sa)
#%% Carbonation degree
mass0 = results['portlandite'][0]*mm_CH + results['calcite'][0]*mm_CC
mass_f = results['portlandite'][-1]*mm_CH + results['calcite'][-1]*mm_CC
print(mass0*1e-15*scale)
print(mass_f*1e-15*scale)

plt.figure(figsize=(8,4))
d = ((np.array(results['portlandite'])*mm_CH + \
              np.array(results['calcite'])*mm_CC)- \
              results['portlandite'][0]*mm_CH)*1e-15*scale
t = np.array(results['time'])
print(d[-1])
plt.plot(t, d)
plt.xlabel('Time (s)')
plt.xlim(-1,264)
plt.ylabel('Mass increase (g/um2)')
plt.show

#%%

b = np.array([0.8026666666666116,
 2.4080000000003117,
 7.224000000006197,
 14.19000000001471,
 23.99399999994149,
 36.20599999982626,
 50.8833333330211,
 68.08333333303288,
 87.66266666701618,
 109.6213333344118,
 134.76200000191304,
 162.62600000283803,
 193.81533333720674,
 227.38400000498777,
 263.4466666715817])


x = np.arange(0, len(b))
plt.figure(figsize=(8,4))
plt.plot(b, x)
plt.ylabel("Dissolved length (um)")
plt.xlabel("Time")
plt.xlim(-1,264)
plt.legend()
plt.show()

#%% Depth -> mass increase
dmm  = mm_CC-mm_CH #[g/mol] molar mass difference
mi = x/0.0331*1e-15*dmm
nt = b[-2]/time[-2]*time
plt.figure(figsize=(8,4))
plt.plot(b, mi)
plt.plot(nt, mass_sa)
plt.ylabel("Mass increase (g/um2)")
plt.xlabel("Time")
plt.xlim(-1,264)
plt.legend()
plt.show()
