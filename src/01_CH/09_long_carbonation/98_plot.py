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
#name = '00_ccD11_prev'
name = '04_ccD13_Di'
path = root_dir+'\\results\\output\\13_validation\\' + name + '\\'
results = fn.load_obj(path + '04_ccD13_Di' +'_results')
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

#%%

b = np.array([0.070945216,
0.243092284,
0.58738642,
1.103827624,
1.792415895,
2.653151235,
3.686033642,
4.891063118,
6.268239661,
7.817563273,
9.539033952,
11.4326517,
13.49841651,
15.7363284,
18.14638735,
])


x = np.arange(0, len(b))
plt.figure(figsize=(8,4))
plt.plot(b, x)
plt.plot(time, x)
plt.ylabel("Thickness (um)")
plt.xlabel("Time (d)")
plt.xlim(-1,20)
plt.legend()
plt.show()
