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

from copy import deepcopy
#%%

from yantra.physics.PhrqcReactiveTransport import PhrqcReactiveTransport 
class DissolutionRT(PhrqcReactiveTransport):
        
    def advance(self):     
        if(self.iters==0):    
            self.prev_port = self.solid.portlandite.c >0
            self.dissolution_time = []
            self.set_volume()
            self.set_porosity()              
            self.nodetype = deepcopy(domain.nodetype)
            
        self.fluid.call('advance')
        
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0): 
            self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid.poros) 
        c=deepcopy( self.fluid.get_attr('c'))
        ss=self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        phaseqty=self.solid.update(self.phrqc.dphases)
        if len(phaseqty):
            self.phrqc.modify_solid_phases(phaseqty)
        self.fluid.set_attr('ss',ss)       
        self.update_nodetype()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)
        self.update_solid_params()     
         
    def set_volume(self):
        self.solid.vol = np.zeros(self.solid.shape)
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            self.solid.vol += val.c * self.solid.mvol[num-1] 
            
    def set_porosity(self):        
        self.solid.poros=1.- self.solid.vol/self.solid.voxel_vol 
        self.solid.app_tort = 1. * self.solid.poros ** (1./3.)
        if np.any(self.solid.poros<=1e-10):
            sys.exit('Negative or zero porosity')   
    
        
    def update_solid_params(self):
        self.set_volume()
        self.set_porosity()
    
    def update_nodetype(self):
        is_port = self.solid.portlandite.c > 0
        if ~(np.all(self.prev_port == is_port)):
            nt = is_port * (-5) +(~is_port) * (-1) 
            self.nodetype = deepcopy(nt)
            self.solid.nodetype = deepcopy(nt)
            self.dissolution_time.append(self.time)
        self.prev_port = is_port
        

#%% geometry

ll = 50
l = 10 +ll
lx = l*1.0e-6
ly = 2.0e-6
dx = 1.0e-6
domain = yantra.Domain2D((0,0),(lx,ly),dx, grid_type = 'nodal')
domain.nodetype[:, ll+1:l+ll+1] = -5

#%% params
slabels =  100003*(domain.nodetype== -1) + 100002*(domain.nodetype!= -1)
mvol = 0.5
porosity = 0.9
pqty =  (1-porosity)/mvol*(domain.nodetype==-5)
porosity = 1.0*(domain.nodetype!=-5) + porosity*(domain.nodetype==-5)
D = 1e-9#*(domain.nodetype== -1)+1e-09*(domain.nodetype!= -1) #1e-15
#domain params

domain_params={}
domain_params['D0']=D
domain_params['database']='cemdata07.dat'
domain_params['phrqc_input_file']='portlandite_mlvl.phrq'
domain_params['solution_labels']=slabels
domain_params['eq_names']=['portlandite']
domain_params['solid_phases']={'portlandite':{'type':'diffusive','mvol':mvol,'c':pqty}}
domain_params['voxel_vol']=1
domain_params['poros']=porosity
#solver parameters
solver_params={}
solver_params['collision_model']='trt'
solver_params['tauref']=1
solver_params['magic_para']= 1./4.
solver_params['cphi_fact']= 1./3.#*3.

bc_params={'solution_labels':{'left':100001}, 
           'top':['flux', 0.0],
           'bottom':['flux', 0.0],
           'left':['conc', 0.0],
           'right':['flux', 0.0]}
rt= DissolutionRT('MultilevelAdvectionDiffusion2', domain,
                                  domain_params,bc_params,solver_params)
rt.phrqc.phrqc_flags['smart_run'] = True
#%% run model

conc_Ca = 'Ca\n'
conc_H = 'H\n'
conc_O = 'O\n'
pqty_CH = 'CH\n'
porosity = 'Poros\n'

time=[]
TotCa =[]
TotCH = []
Tf=1.
iters = 1000
while  rt.time<=Tf:#rt.iters < iters:#
    rt.advance() 
    TotCa.append(np.sum(rt.fluid.Ca.c*rt.fluid.Ca.poros)) 
    TotCH.append(np.sum(rt.solid.portlandite.c))
    time.append(rt.time)     
    
    conc_Ca += str(rt.iters) + ': ' + str(np.array(rt.fluid.Ca.c[1,:]) + np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:])) +'\n'
    conc_H += str(rt.iters) + ': ' + str(np.array(rt.fluid.H.c[1,:]) + np.array(rt.fluid.H._ss[1,:])/np.array(rt.phrqc.poros[1,:]))+'\n'
    conc_O += str(rt.iters) + ': ' + str(np.array(rt.fluid.O.c[1,:]) + np.array(rt.fluid.O._ss[1,:])/np.array(rt.phrqc.poros[1,:]))+'\n'
    pqty_CH += str(rt.iters) + ': ' + str(np.array(rt.solid.portlandite.c[1,:]))+'\n'
    porosity += str(rt.iters) + ': ' + str(np.array(rt.solid.poros[1,:]))+'\n'
        
#%%plot results
'''
plt.figure()
plt.plot(time, TotCH)
plt.xlabel('Time [s]')
plt.ylabel('Total CH mass [mM]')
plt.show()

plt.figure()
plt.plot(time,(np.array(TotCH)/np.sum(pqty)))
plt.xlabel('Time [s]')
plt.ylabel('Total CH fraction [mM]')
plt.show()
'''
plt.figure()
plt.plot(rt.fluid.Ca.c[1,:])
plt.title('Ca concentarion [mol/l]')
plt.show()

#%% print 
print(rt.dissolution_time)   
'''   
print('Ca')       
print(rt.fluid.Ca._c[1,:])
print('Ca + ss')       
print(rt.fluid.Ca._c[1,:]+rt.fluid.Ca._ss[1,:]/rt.phrqc.selected_output()['poros'][1,:])
print('CH')       
print(rt.solid.portlandite.c[1,:])
print('De')       
print(rt.fluid.Ca.De[1,:])
print('tau')       
print(rt.fluid.Ca.tau[1,:])
print(rt.dt)
'''
'''
with open('conc1.txt', 'w') as file:
    file.write(conc_Ca)
    file.write(conc_O)
    file.write(conc_H)
    file.write(pqty_CH)
    file.close()
'''