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
Benchmark 5: Comparison of different collision models for domain with varying diffusion \
coefficient using Yantra's 2D implementation. The benchmark example refers to example of \
Perko and Patel (2014) PRE paper Fig 5 and results in Fig 6.
"""
#%% import necessray
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
import numpy as np
    
#%% generate domain instance
dx = 0.1*1e-3 #0.1 um
lx = 100*dx
ly = 5*dx
domain = yantra.Domain2D((0,0), (lx,ly),dx, grid_type = 'midway')
x,_ = domain.meshgrid()
D  = 2.2e-10
#%% define physics
#domain params
domain_params={}
domain_params['D']=D
domain_params['c']=0.1 
domain_params['poros']=1 * (x<=20*dx) + 0.05*(x>20*dx)*(x<80*dx) +  1 * (x>=80*dx)

bc_params = {}
bc_params['left']=['c',0.]
bc_params['right']=['flux',0]       
solver_params={}         
solver_params['lattice']='D2Q5'
domain_params['D0']=D
solver_params['collision_model']='trt'
solver_params['magic_para']=1./4.
#solver_params['cphi']=1./3.
solver_params['cphi_fact']=1./3.
solver_params['tfactbased']=True
solver_params['tfact']=1./6.

ade=yantra.MultilevelAdvectionDiffusion(domain,domain_params,bc_params,solver_params)
print('Current dt %s' %str(ade.dt))

#%% run models
tf=24*3600 #1 h
ade.run(time=tf)
print('trt collision model finished. Number of iterations: %s'%ade.iters)
#%% plot data to generate fig 6 of the paper
y = int(ade.ny/2)
plt.figure(figsize=(10,10))
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20
plt.plot(x[y,:],ade.c[y,:],'-.',label=ade.collision_model.upper())
plt.legend(loc=4)
plt.xlabel('Distance (m)')
plt.ylabel(r'concentration (mol/m$^3$)')
plt.title('Time %s h'%(ade.time/(3600)))
plt.show()

#%% Criteria
print('Current dt %s' %str(ade.dt))
print('Von Neumann: dt should be less then %s' %str(0.5*dx*dx/np.max(ade.De)))
#print('Von Neumann: dt=1 should be less then %s' %str(0.5*(dx/ade.convfactors['L'])**2/(np.max(ade.De)/ade.convfactors['D'])))
#print('dt should be less then %s' %str(0.5*dx*dx*np.min(ade.poros)/np.max(ade.De)))