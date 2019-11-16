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
from yantra.physics._phrqc_wrapper import  Phrqc
from yantra._base import Multicomponent 
from yantra.physics.PhrqcReactiveTransport import PhrqcReactiveTransport 
from yantra.physics.PhrqcReactiveTransport import Solid

from copy import deepcopy
#%%

class DissolutionRT(PhrqcReactiveTransport):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params):
        '''
        Reactive transport class for carbonation
        A Lattice Boltzmann method based reactive transport model \
        where reactions are computed using geochemical solver PHREEQC
        '''
        self.auto_time_step = solver_params.get('auto_time_step',True)
        self.phrqc = CarbonationPhrqc(domain,domain_params,bc_params,solver_params)        
        components = self.phrqc.components
        bc = self.phrqc.boundary_conditions 
        init_c = self.phrqc.initial_conditions         
        for name in components:
            if name not in domain_params:
                domain_params[name]={}
            domain_params[name]['c']=init_c[name]            
        for name in bc:
            for comp in components:
                if comp not in bc_params:
                    bc_params[comp]={}
                bc_params[comp][name]=['c',bc[name][comp]]                
        self.fluid=Multicomponent(eqn,components,domain,domain_params,
                                  bc_params,solver_params)
        self.solid=Solid(domain,domain_params,solver_params)   
        self.ptype = 'CSH' if hasattr(self.fluid, 'Si') else 'CH'
        self.set_volume()
        self.set_porosity()
        self.nodetype = deepcopy(domain.nodetype)
        
    def advance(self):     
        if(self.iters==0):                    
            self.set_volume()
            self.set_porosity()               
            self.nodetype = deepcopy(domain.nodetype)
        #self.update_solid_params()
        #self.update_diffusivity()  
        #self.fluid.call('_set_relaxation_params')  
            
        #self.fluid.call('advance')
        #self.fluid.set_attr('Deref', self.Deref,component_dict=False)  
        
        
        
        self.fluid.Ca.advance()
        '''
        self.fluid.Ca.time +=self.dt
        self.fluid.Ca.iters+=1
        
        print("=======================================")
        print("step %s" %self.fluid.Ca.iters)
        print("=======================================")
        print('tau')       
        print(self.fluid.Ca.tau[1,:])
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
        
        print("\tsum f: \t N3 %s \t N4 %s \t N5 %s" %(np.sum(self.fluid.Ca.f[1,3,:]),
                                          np.sum(self.fluid.Ca.f[1,4,:])/self.fluid.Ca.poros[1,4],
                                          np.sum(self.fluid.Ca.f[1,5,:])/self.fluid.Ca.poros[1,5]))
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
        '''
        self.fluid.H.advance()
        self.fluid.O.advance()   
        self.update_solid_params()      
        #prev_port = self.solid.portlandite.c[1,4] 
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):
            self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid.poros) 
        c=deepcopy( self.fluid.get_attr('c'))
        ss=self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        phaseqty=self.solid.update(self.phrqc.dphases)
        if len(phaseqty):
            self.phrqc.modify_solid_phases(phaseqty)
        #self.solid.portlandite.c[1,4] = prev_port  
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
  
        

class CarbonationPhrqc(Phrqc):
    pass
#%% geometry

ll = 20
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
Dref = 1.*1.e-11
Dhigh = 1e-9
D = 1e-9#*(domain.nodetype== -1)+1e-15*(domain.nodetype!= -1) #1e-15
#D[:, ll+1]=Dref
#domain params
domain_params={}
domain_params['D0']=D
domain_params['Deref']=Dref
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
solver_params['tauref']=1#0.2*Dhigh/Dref+0.5#5.5 for 1e-10
#solver_params['Deref']=Dref
solver_params['magic_para']= 1./4.
solver_params['cphi_fact']= 1./3.
#solver_params['cphi']=1./3.
#solver_params['tfact']= 1./6.
#['tfactbased']= True

bc_params={#'solution_labels':{'left':100001}, 
           'top':['flux', 0.0],
           'bottom':['flux', 0.0],
           'left':['conc', 0.0],
           'right':['flux', 0.0]}
rt= DissolutionRT('MultilevelAdvectionDiffusion2', domain,
                                  domain_params,bc_params,solver_params)
#%% run model

conc_Ca = 'Ca\n'
conc_H = 'H\n'
conc_O = 'O\n'
pqty_CH = 'CH\n'
porosity = 'Poros\n'

time=[]
TotCa =[]
TotCH = []
iters = 500
while  rt.iters < iters:#rt.time<=0.1:#
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
plt.figure()
plt.plot(time,TotCa)
plt.xlabel('Time [s]')
plt.ylabel('Total Ca mass in aqeuous phase [mM]')
plt.show()


plt.figure()
plt.plot(time, TotCH)
plt.xlabel('Time [s]')
plt.ylabel('Total CH mass [mM]')
plt.show()

plt.figure()
plt.plot(rt.fluid.Ca.c[1,:])
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
print('tau')       
print(rt.fluid.Ca.tau[1,:])
print('porosity')       
print(rt.phrqc.selected_output()['poros'][1,:])

'''
with open('conc1.txt', 'w') as file:
    file.write(conc_Ca)
    file.write(conc_O)
    file.write(conc_H)
    file.write(pqty_CH)
    file.close()
'''