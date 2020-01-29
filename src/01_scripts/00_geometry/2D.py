# -*- coding: utf-8 -*-
'''
2D geometry
'''

#%% PYTHON MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import yantra
import cell_type as ct # change the path to cell_type file

np.random.seed(0)
#import phrqc
#%% PROBLEM DEFINITION
__doc__= """ 
Default example. Explains possible values
1D carbonation of portlandite/cement.
"""
#problem type
m = 'CH' #or 'CSH' #TODO case for cement

#%% RANDOM GEOMETRY
lx = 31*1.0e-6
ly = 31*1.0e-6
dx = 1.0e-6
rad = 6*dx
#wall_len_y = wall_len_x 
domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(ly, lx), 
                         dx=dx, 
                         grid_type='nodal')
ra = np.random.random(domain.nodetype.shape)
ra = ra>=0.6

domain.nodetype[ra] = ct.Type.MULTILEVEL
domain.nodetype[:,0] = ct.Type.LIQUID
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()

#%% CUBES
lx = 31*1.0e-6
ly = 14*1.0e-6 #32
dx = 1.0e-6
#wall_len_y = wall_len_x 
domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
for i in range(1,6):
    for j in range(1,6):
        domain.draw_rect(center = (6*i*dx-3*dx,6*j*dx-2*dx), lengths=(2*dx, 2*dx), idx=ct.Type.MULTILEVEL)

domain.nodetype[:,0] = ct.Type.LIQUID
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.colorbar()
plt.show()
