# -*- coding: utf-8 -*-
'''
Import csv file with coordinates of the rough surface
'''
#%% 
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
np.set_printoptions(precision=5, threshold=np.inf)
import yantra
import cell_type as ct # change the path to cell_type file

#%%
df = pd.read_csv('input/rough_surface.csv')
data = df.values
x = data[:, 1]
y = data[:, 2]
h = data[:, 3]
surface = df.drop(df.columns[[0]], axis=1)
surface_np = np.array(surface['h'])

#%% Geometry
lx = np.max(x)
ly = np.max(y)
surface_rs = surface_np.reshape([int(lx+1), int(ly+1)])
surface_rs = surface_rs - np.min(surface_rs) + 5
lz = int(np.max(surface_rs)) + 20
dx = 1

domain = yantra.domain.Domain3D(corner=(0, 0, 0), 
                              lengths=(lx, ly, lz), 
                              dx=dx, 
                              grid_type='nodal')

# place solid cells
for i in range(0, int(lz+1)):
    is_solid = (i <=surface_rs[:,:])
    domain.nodetype[i] = ct.Type.SOLID * is_solid
    
#%%
domain_2d = yantra.domain.Domain2D(corner=(0, 0), 
                                  lengths=(lx, ly), 
                                  dx=dx, 
                                  grid_type='nodal')
Z = 12#int(lz/dx/3)
Y = 20
domain_2d.nodetype = domain.nodetype[:, Y,:]#[Z, :, :]

plt.figure(figsize=(5,5))
plt.imshow(domain_2d.nodetype[1:-1, 1:-1]) 
plt.show()

#%%
#np.save('output/suface3D.npy',domain.nodetype)