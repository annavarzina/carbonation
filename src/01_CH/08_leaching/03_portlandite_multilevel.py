#!/usr/bin/python
# -*- coding: utf-8 -*-
#=======================================================================================
#This File is part of Yantra: A lattice Boltzmann method based tool for multiscale/
#multiphyics simulations
#=======================================================================================
#
#Copyright (C) 2016-2017  <Author> Ravi A. Patel <Email> ravee.a.patel@gmail.com
#
#This program is free software: you can redistribute it and/or modify it under the
#terms of the GNU General Public License as published by the Free Software 
#Foundation, either version 3 of the License, or any later version.
#This program is distributed in the hope that it will be useful, but WITHOUT ANY 
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
#PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#=======================================================================================
from __future__ import division,print_function
__doc__="""
Portlandite dissolution
"""

#%% import modules
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import yantra
import numpy as np
import matplotlib.pylab as plt
import copy
#%% geometry
ll = 3
l = 5 +ll
lx = l*1.0e-6
ly = 1.0e-6
dx = 1.0e-6
domain = yantra.Domain2D((0,0),(lx,ly),dx, grid_type = 'midway')
domain.nodetype[:, ll+1:l+ll+1] = -5
pqty = 1.*(domain.nodetype>0)
#%% params
slabels = np.zeros(domain.nodetype.shape)
slabels =  100001*(domain.nodetype== -1) + 100002*(domain.nodetype!= -1)
pqty = 1.*(domain.nodetype==-5)
porosity = np.ones(pqty.shape)
D = 1e-9#*(domain.nodetype== -1)+1e-09*(domain.nodetype!= -1) #1e-15
#domain params

domain_params={}
domain_params['D0']=D
domain_params['database']='cemdata07.dat'
#domain_params['phrqc_input_file']='02_mlvl_portlandite.phrq'
domain_params['phrqc_input_file']='portlandite_mlvl.phrq'
domain_params['solution_labels']=slabels
domain_params['eq_names']=['portlandite']
domain_params['solid_phases']={'portlandite':{'type':'diffusive','mvol':1,'c':pqty}}
domain_params['voxel_vol']=1
domain_params['poros']=porosity
#solver parameters
solver_params={}
solver_params['collision_model']='trt'
solver_params['magic_para']= 1./4.
solver_params['cphifact']= 1./3.
#solver_params['phrqc_flags']['smart_run']=True
#solver_params['phrqc_smart_run_tol']=1e-8

bc_params={'top':['flux', 0.0],
           'bottom':['flux', 0.0],
           'left':['flux', 0.0],
           'right':['flux', 0.0]}
rt= yantra.PhrqcReactiveTransport('MultilevelAdvectionDiffusion', domain,
                                  domain_params,bc_params,solver_params)
#%% run model
time=[]
AvgCa =[]
TotCH = []
iters = 20
while  rt.iters < iters:#rt.time<=0.1:#
    rt.fluid.call('advance')
    c=copy.deepcopy(rt.fluid.get_attr('c'))
    #advance phrqc
    ss=rt.phrqc.modify_solution(c,rt.dt,rt.solid.nodetype)
    phaseqty=rt.solid.update(rt.phrqc.dphases)
    if len(phaseqty):
        rt.phrqc.modify_solid_phases(phaseqty)
    rt.fluid.set_attr('ss',ss)
    rt.fluid.set_attr('nodetype',rt.solid.nodetype,component_dict=False)
    
    AvgCa.append(np.sum(rt.fluid.Ca.c)/np.sum(rt.fluid.Ca.nodetype<=0)) 
    TotCH.append(np.sum(rt.solid.portlandite.c))
    time.append(rt.time)
        
#%%plot results
plt.figure()
plt.plot(time,AvgCa)
plt.xlabel('Time [s]')
plt.ylabel('Avg. Ca conc in aqeuous phase [mM]')
plt.show()


plt.figure()
plt.plot(time, TotCH)
plt.xlabel('Time [s]')
plt.ylabel('Total CH mass [mM]')
plt.show()

plt.figure()
plt.imshow(rt.fluid.Ca.c)
plt.colorbar()
plt.title('Ca concentarion [mol/l]')
plt.show()

#%% print       
print('Ca')       
print(rt.fluid.Ca._c[1,:])
print('Ca + ss')       
print(rt.fluid.Ca._c[1,:]+rt.fluid.Ca._ss[1,:])
print('CH')       
print(rt.solid.portlandite.c[1,:])
print('De')       
print(rt.fluid.Ca.De[1,:])
