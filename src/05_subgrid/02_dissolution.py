# -*- coding: utf-8 -*-
'''
Example with precipitation everywhere
Fixed PCO2 at the boundary
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
import time
import yantra
import cell_type as ct # change the path to cell_type file
import misc_func as fn
import rt2
#import phrqc
#%% PROBLEM DEFINITION
__doc__= """ 
Reference:
    Poros = 0.05
    PCO2 = 3.4
    IE = 0.5
    Archies relation for diffusivity
"""
#problem type
m = 'CH' #or 'CSH'

ll = 2
l = 5 +ll
dx = 1.0e-3
lx = l*dx
ly = 2*dx

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, ll+1:l+ll] = ct.Type.MULTILEVEL

domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()
#%%  VALUES
nn='01_test'
slabels = np.zeros(domain.nodetype.shape)
slabels =  100001*(domain.nodetype!= -5) + 100002*(domain.nodetype== -5)
pqty = 0.5*(domain.nodetype==-5)
porosity = 0.5*(domain.nodetype==-5) + 1.0*(domain.nodetype!=-5)
D = 1e-9*(domain.nodetype== -1)+1e-09*(domain.nodetype!= -1) #1e-15
#domain params

domain_params={}
domain_params['D0']=D
domain_params['database']='cemdata07.dat'
domain_params['phrqc_input_file']='02_mlvl_portlandite.phrq'
domain_params['solution_labels']=slabels
domain_params['eq_names']=['portlandite']
domain_params['solid_phases']={'portlandite':{'type':'diffusive','mvol':1,'c':pqty}}
domain_params['voxel_vol']=1
domain_params['poros']=porosity
#solver parameters
solver_params={}
#solver_params['tfactbased'] = True
#solver_params['tfact'] = 1./6./8

solver_params['collision_model']='trt'
solver_params['magic_para']=1.0/4.0
solver_params['cphi_fact']=1.0/3.0
          
bc_params = {'solution_labels':{'left':100003}, 
            'top':['flux', 0.0],
            'bottom':['flux', 0.0],
            'left':['flux', 0.0],
            'right':['flux', 0.0],}
            

#%% INITIATE THE SOLVER
rt= rt2.DissolutionRT('MultilevelAdvectionDiffusion',  domain, 
                          domain_params, bc_params, solver_params) 
rt.phrqc.phrqc_smart_run_tol = 1e-6
#%% PARAMETERS
#plist =  [(1,2), (1,3), (1,4), (1,5), (1,6), (1,7), (1,8), (1,9), (1,10)]
plist =  [(1,n) for n in np.arange(0, l)]
pavglist = ['avg_poros', 'pH', 'avg_D_eff', 'sum_vol', 'precipitation',
            'dissolution', 'portlandite_cells', 'calcite_cells', 'dt'] 
#'delta_ch', 'delta_cc', 'precipitation','dissolution', 'portlandite_cells', 
#'calcite_cells', 'active_cells','dt', 'pH', 'avg_poros',  'avg_D_eff', 'sum_vol'
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)

#%% TIME SETTINGS
itr = 0 
j = 0
nitr = 1000
Ts = 1000*3600
rt_port = []
rt_time = []
#%% RUN SOLVER
while  itr <= nitr:# rt.time <=Ts: #
    rt.advance()      
    rt_port.append(np.sum(rt.solid.portlandite.c))
    rt_time.append(rt.time/3600)
    itr += 1
    
#%% SIMULATION TIME
#'''
print('Ca %s' %str(np.array(rt.fluid.Ca.c[1,:])))
print('Ca +ss/theta %s' %str(np.array(rt.fluid.Ca.c[1,:]) + np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.fluid.Ca.poros[1,:])))
print('H +ss %s' %str(np.array(rt.fluid.H.c[1,:]) + np.array(rt.fluid.H._ss[1,:])/np.array(rt.fluid.H.poros[1,:])))
print('O +ss %s' %str(np.array(rt.fluid.O.c[1,:]) + np.array(rt.fluid.O._ss[1,:])/np.array(rt.fluid.O.poros[1,:])))
print('CH %s' %str(np.array(rt.solid.portlandite.c[1,:])))
print('dCH %s' %str(np.array(rt.phrqc.dphases['portlandite'][1,:])))
#fn.plot_fields(carb_rt, names={ 'calcite', 'portlandite', 'Ca', 'C'})
#print(rt.phrqc.selected_output())
#'''

#%%
print('time %s hours' %str(rt.time/3600))
plt.figure()
plt.plot(rt_time, rt_port)
plt.show()
