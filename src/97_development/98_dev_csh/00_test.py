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
Benchmark 11: C-S-H dissolution benchmark as described in chapter 4 Patel (2016), PhD thesis, UGent
"""

#%% import modules
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import matplotlib.pylab as plt
#%%generate mesh
domain = yantra.Domain2D((0,0),(5e-6,5e-6/10),5e-6/50., grid_type = 'nodal')

#%%define physics
domain_params={}
domain_params['D0']=2.2e-10
domain_params['database']='cemdata07.dat'
domain_params['solution_labels']=100005
domain_params['phrqc_input_file']='benchmark11.phrq'
domain_params['ss_names']={'Tob_jen_ss':['CSHQ_TobH','CSHQ_TobD','CSHQ_JenH','CSHQ_JenD']}
domain_params['solid_phases']={
 'CSHQ_TobH': {'c':0.1041,'mvol':55.30e-3,'type':'diffusive'},
 'CSHQ_TobD': {'c':2.5050,'mvol':47.95e-3,'type':'diffusive'},
 'CSHQ_JenH': {'c':2.1555,'mvol':75.63e-3,'type':'diffusive'},
 'CSHQ_JenD': {'c':3.2623,'mvol':80.58e-3,'type':'diffusive'},}
domain_params['tort_model_params']=[1.,1./3.]  
#bc params
bc_params={}
bc_params['left']=['conc',0]
bc_params['right']=['flux',0]
bc_params['solution_labels']={'left':100001}
#generate physics
solver_params={}
solver_params['phrqc_flags']={}
solver_params['phrqc_flags']['smart_run']=True
solver_params['phrqc_smart_run_tol']=1e-8
rt=yantra.PhrqcReactiveTransport('MultilevelDiffusion',domain,domain_params,bc_params,solver_params)
#%%run model
Ts=1.
while rt.time < Ts:
    rt.advance()
    if rt.iters%10==0:
        print ('Time: %s, Time step: %s'%(rt.time,rt.dt))
#%%plot results
loc_y = int(domain.ny/2)
plt.rc('text', usetex=False)
f, axarr = plt.subplots(2,2, sharex=False,sharey=False)
#Ca aqueous
axarr[0,0].set_title('Ca (aqueous)')
axarr[0,0].plot(np.array(domain.x)*1e6,rt.fluid.Ca.c[loc_y,:]*1e3)
axarr[0,0].set_xlabel(r'Distance (\mu m)')
axarr[0,0].set_ylabel('Concentration (mmol/lit)')
#Ca Si
axarr[0,1].set_title('Si (aqueous)')
axarr[0,1].plot(np.array(domain.x)*1e6,rt.fluid.Si.c[loc_y,:]*1e3)
axarr[0,1].set_xlabel(r'Distance (\mu m)')
axarr[0,1].set_ylabel('Concentration (mmol/lit)')
#Ca solid
Ca_s = (0.8333333*rt.solid.CSHQ_TobD.c[loc_y,:] + 0.6666667*rt.solid.CSHQ_TobH.c[loc_y,:] + 
        1.3333333*rt.solid.CSHQ_JenH.c[loc_y,:] + 1.5*rt.solid.CSHQ_JenD.c[loc_y,:])
axarr[1,0].set_title('Si (solid)')
axarr[1,0].plot(np.array(domain.x)*1e6,Ca_s)
axarr[1,0].set_xlabel(r'Distance (\mu m)')
axarr[1,0].set_ylabel('Concentration (mol/lit)')
#Si solid
Si_s = (0.6666667*rt.solid.CSHQ_TobD.c[loc_y,:] + 1.0*rt.solid.CSHQ_TobH.c[loc_y,:] + 
        1.0*rt.solid.CSHQ_JenH.c[loc_y,:] + 0.6666667*rt.solid.CSHQ_JenD.c[loc_y,:])
axarr[1,1].set_title('Si (solid)')
axarr[1,1].plot(np.array(domain.x)*1e6,Si_s)
axarr[1,1].set_xlabel(r'Distance (\mu m)')
axarr[1,1].set_ylabel('Concentration (mol/lit)')
plt.show()        
        
