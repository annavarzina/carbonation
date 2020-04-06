# -*- coding: utf-8 -*-
#=======================================================================================
from __future__ import division,print_function

import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import yantra
import numpy as np
from numpy import exp
from scipy.special import erfc 
import matplotlib.pylab as plt

#%% domain generation
l = 1e-5
domain=yantra.Domain2D((0,0),(l,l),l/6,grid_type='midway')
domain.nodetype[3,3]=-2

#%%set Multilevel advection diffusion model
Dm= 1e-11*(domain.nodetype==-2)+\
    1e-9*(domain.nodetype==-1)
poros =1.0*(domain.nodetype==-2)+\
       1.0*(domain.nodetype==-1)
c = 1.*(domain.nodetype==-2)+\
    0.*(domain.nodetype==-1) 
#domain parameters
domain_params={}
domain_params['D0']=Dm
domain_params['app_tort']=1.
domain_params['poros']=poros
domain_params['c']=c
#boundary parameters
bc_params = {}
bc_params['top']=['flux',0.0]
bc_params['left']=['flux',0.0]
bc_params['right']=['flux',0.0]
bc_params['bottom']=['flux',0.0]
#solver parameters
solver_params={}
#create instance
ade=yantra.MultilevelAdvectionDiffusion(domain,domain_params,bc_params,solver_params)
#%% run model
'''
ts = 1
for i in range(0,ts):
    ade.advance()
    print(ade.c[3,3])
'''
ade.run(time = 0.1)   

#%%
print('===')
print(ade.dt)
print(ade.time)
#print(ade.De)
#print(ade.c)
#print(ade.c[3,3])

plt.imshow(ade.c)
plt.colorbar()
plt.show()