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
            self.set_volume()
            self.set_porosity()        
            self.nodetype = deepcopy(domain.nodetype)
            self.update_diffusivity()
            self.fluid.call('_set_relaxation_params')  
            
        #self.fluid.call('advance')
        
        
        
        #Ca self.fluid.Ca.advance()
        self.fluid.Ca.time +=self.dt
        self.fluid.Ca.iters+=1
        print("=======================================")
        print("step %s" %self.fluid.Ca.iters)
        print("=======================================")
        print("Before transport ")
        print("\tf: \n\t N3 %s \n\t N4 %s \n\t N5 %s" %(self.fluid.Ca.f[1,3,:],
                                          self.fluid.Ca.f[1,4,:],
                                          self.fluid.Ca.f[1,5,:]))
        print("\tC:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.c[1,3],
                                          self.fluid.Ca.c[1,4],
                                          self.fluid.Ca.c[1,5]))
        print("\tSS:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca._ss[1,3],
                                          self.fluid.Ca._ss[1,4],
                                          self.fluid.Ca._ss[1,5]))
        print("\tPorosity:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.poros[1,3],
                                          self.fluid.Ca.poros[1,4],
                                          self.fluid.Ca.poros[1,5]))
        print("\tCa+source: \t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.c[1,3]+self.fluid.Ca._ss[1,3]*self.fluid.Ca.poros[1,3],
                                          self.fluid.Ca.c[1,4]+self.fluid.Ca._ss[1,4]/self.fluid.Ca.poros[1,4],
                                          self.fluid.Ca.c[1,5]+self.fluid.Ca._ss[1,5]/self.fluid.Ca.poros[1,5]))
        self.fluid.Ca.collide()
        print("Collide")
        print("\tf: \n\t N3 %s \n\t N4 %s \n\t N5 %s" %(self.fluid.Ca.f[1,3,:],
                                          self.fluid.Ca.f[1,4,:],
                                          self.fluid.Ca.f[1,5,:]))
        print("\tC:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.c[1,3],
                                          self.fluid.Ca.c[1,4],
                                          self.fluid.Ca.c[1,5]))
        print("\tSS:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca._ss[1,3],
                                          self.fluid.Ca._ss[1,4],
                                          self.fluid.Ca._ss[1,5]))
        print("\tPorosity:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.poros[1,3],
                                          self.fluid.Ca.poros[1,4],
                                          self.fluid.Ca.poros[1,5]))
        self.fluid.Ca.stream()
        print("Stream")
        print("\tf: \n\t N3 %s \n\t N4 %s \n\t N5 %s" %(self.fluid.Ca.f[1,3,:],
                                          self.fluid.Ca.f[1,4,:],
                                          self.fluid.Ca.f[1,5,:]))
        print("\tC:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.c[1,3],
                                          self.fluid.Ca.c[1,4],
                                          self.fluid.Ca.c[1,5]))
        print("\tSS:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca._ss[1,3],
                                          self.fluid.Ca._ss[1,4],
                                          self.fluid.Ca._ss[1,5]))
        print("\tPorosity:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.poros[1,3],
                                          self.fluid.Ca.poros[1,4],
                                          self.fluid.Ca.poros[1,5]))
        self.fluid.Ca.apply_bc()
        print("BC")
        print("\tf: \n\t N3 %s \n\t N4 %s \n\t N5 %s" %(self.fluid.Ca.f[1,3,:],
                                          self.fluid.Ca.f[1,4,:],
                                          self.fluid.Ca.f[1,5,:]))
        print("\tC:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.c[1,3],
                                          self.fluid.Ca.c[1,4],
                                          self.fluid.Ca.c[1,5]))
        print("\tSS:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca._ss[1,3],
                                          self.fluid.Ca._ss[1,4],
                                          self.fluid.Ca._ss[1,5]))
        print("\tPorosity:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.poros[1,3],
                                          self.fluid.Ca.poros[1,4],
                                          self.fluid.Ca.poros[1,5]))
        self.fluid.Ca.compute_macro_var()
        print("Macro")
        print("\tf: \n\t N3 %s \n\t N4 %s \n\t N5 %s" %(self.fluid.Ca.f[1,3,:],
                                          self.fluid.Ca.f[1,4,:],
                                          self.fluid.Ca.f[1,5,:]))
        print("\tC:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.c[1,3],
                                          self.fluid.Ca.c[1,4],
                                          self.fluid.Ca.c[1,5]))
        print("\tSS:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca._ss[1,3],
                                          self.fluid.Ca._ss[1,4],
                                          self.fluid.Ca._ss[1,5]))
        print("\tPorosity:\t N3 %s \t N4 %s \t N5 %s" %(self.fluid.Ca.poros[1,3],
                                          self.fluid.Ca.poros[1,4],
                                          self.fluid.Ca.poros[1,5]))
        self.fluid.H.advance()
        self.fluid.O.advance()
        
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):            
            self.update_solid_params()        
            self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid.poros) 
        c=deepcopy( self.fluid.get_attr('c'))
        ss=self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        phaseqty=self.solid.update(self.phrqc.dphases)
        if len(phaseqty):
            self.phrqc.modify_solid_phases(phaseqty)
        self.fluid.set_attr('ss',ss)
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)
         
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
  
    def update_diffusivity(self):
        
        D_CH= 1e-9
        Dref = 1e-9        
        is_port = self.solid.portlandite.c >0
        is_liquid = ~is_port
        De = D_CH *is_port+ Dref * is_liquid 
        Dnew_lb = De/self.solid.poros/self.solid.app_tort
        #print(Dnew_lb)
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)
        

#%% geometry

ll = 3
l = 5 +ll
lx = l*1.0e-6
ly = 2.0e-6
dx = 1.0e-6
domain = yantra.Domain2D((0,0),(lx,ly),dx, grid_type = 'nodal')
domain.nodetype[:, ll+1:l+ll+1] = -5

#%% params
slabels =  100003*(domain.nodetype== -1) + 100002*(domain.nodetype!= -1)
pqty =  0.95*(domain.nodetype==-5)
porosity = 1.0*(domain.nodetype!=-5) + 0.05*(domain.nodetype==-5)
D = 1e-9#*(domain.nodetype== -1)+1e-09*(domain.nodetype!= -1) #1e-15
#domain params

domain_params={}
domain_params['D0']=D
domain_params['database']='cemdata07.dat'
domain_params['phrqc_input_file']='portlandite_mlvl.phrq'
domain_params['solution_labels']=slabels
domain_params['eq_names']=['portlandite']
domain_params['solid_phases']={'portlandite':{'type':'diffusive','mvol':1.,'c':pqty}}
domain_params['voxel_vol']=1
domain_params['poros']=porosity
#solver parameters
solver_params={}
solver_params['collision_model']='trt'
#solver_params['tauref']=1
solver_params['magic_para']= 1./4.
solver_params['cphi_fact']= 1./3.

bc_params={'solution_labels':{'left':100001}, 
           'top':['flux', 0.0],
           'bottom':['flux', 0.0],
           'left':['conc', 0.0],
           'right':['flux', 0.0]}
rt= DissolutionRT('MultilevelAdvectionDiffusion', domain,
                                  domain_params,bc_params,solver_params)
#%% run model
time=[]
AvgCa =[]
TotCH = []
iters = 3
while  rt.iters < iters:#rt.time<=0.1:#
    rt.advance() 
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
print(rt.fluid.Ca._c[1,:]+rt.fluid.Ca._ss[1,:]/rt.phrqc.selected_output()['poros'][1,:])
print('CH')       
print(rt.solid.portlandite.c[1,:])
print('De')       
print(rt.fluid.Ca.De[1,:])
