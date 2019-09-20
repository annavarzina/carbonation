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
Benchmark 10: Plus shaped portlandite dissolution with geometry update. The results  \
of this benchmark are published in Patel et al (2014), Phy & Chem. of Earth.
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
#%%
l = 10
lx = l*1.0e-6
ly = l*1.0e-6
dx = 1.0e-6
domain = yantra.Domain2D((0,0),(lx,ly),dx, grid_type = 'nodal')
#domain.nodetype[3:8,3:8] = -5
domain.nodetype[5,5] = -5
#%%
slabels = np.zeros(domain.nodetype.shape)
slabels =  100001*(domain.nodetype== -1) + 100002*(domain.nodetype!= -1)
pqty = 1.*(domain.nodetype==-5)
porosity = np.ones(pqty.shape)
D = 1e-9*(domain.nodetype== -1)+1e-9*(domain.nodetype!= -1) #1e-15
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
solver_params['collision_model']='trt'
solver_params['tfactbased']=1
solver_params['tfact']= 1./6.
#solver_params['phrqc_flags']['smart_run']=True
#solver_params['phrqc_smart_run_tol']=1e-8
rt= yantra.PhrqcReactiveTransport('MultilevelAdvectionDiffusion', domain,
                                  domain_params,{},solver_params)
#%%run model
time=[]
AvgCa =[]
iters = 1
while rt.iters < iters: #rt.time<=1:#20
    rt.fluid.call('advance')
    c=copy.deepcopy(rt.fluid.get_attr('c'))
    #advance phrqc
    ss=rt.phrqc.modify_solution(c,rt.dt,rt.solid.nodetype)
    phaseqty=rt.solid.update(rt.phrqc.dphases)
    if len(phaseqty):
        rt.phrqc.modify_solid_phases(phaseqty)
    rt.fluid.set_attr('ss',ss)
    rt.fluid.set_attr('nodetype',rt.solid.nodetype,component_dict=False)
    '''
    if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):
        poros=self.solid.porosity()
        app_tort= self.solid.app_tort()
        self.fluid.call('update_transport_params',poros,
                        app_tort,self.auto_time_step)
        self.phrqc.poros=deepcopy(poros)
    '''
    
    AvgCa.append(np.sum(rt.fluid.Ca.c)/np.sum(rt.fluid.Ca.nodetype<=0)) 
    time.append(rt.time)
        
#%%plot results
plt.figure()
plt.plot(time,AvgCa)
plt.xlabel('Time [s]')
plt.ylabel('Avg. Ca conc in aqeuous phase [mM]')
plt.show()

plt.figure()
plt.imshow(rt.fluid.Ca.c)
plt.colorbar()
plt.title('Ca concentarion [mol/l]')
plt.show()

print(rt.fluid.Ca._ss[3:8,3:8] + rt.fluid.Ca._c[3:8,3:8])
print(rt.solid.portlandite.c[5,5])