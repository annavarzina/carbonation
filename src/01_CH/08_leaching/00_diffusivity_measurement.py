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
import numpy as np
__doc__="""
Check diffusivity
"""
#%% import necessray
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
    
#%% generate domain instance
length = 200
dx = 1e-6
lx = length*dx
ly = 3*dx
domain = yantra.Domain2D((0,0), (lx, ly), dx, grid_type = 'nodal')

#%% define physics
#Dlow  = 1e-11
#Dhigh = 2000 * Dlow
#D = Dlow * (x<=0.02) + Dhigh*(x>0.02)*(x<0.08) +  Dlow * (x>=0.08)
#domain params
D = 1e-9
domain_params={}
domain_params['D']=D
domain_params['c']=0.0 
bc_params = {}
bc_params['left']=['c',1.]
bc_params['right']= ['flux',0]
bc_params['top']=['flux',0]
bc_params['bottom']=['flux',0]       
solver_params={}         
solver_params['lattice']='D2Q5'
solver_params['collision_model']='srt'
ade_srt=yantra.AdvectionDiffusion(domain,domain_params,bc_params,solver_params)
solver_params['collision_model']='trt'
solver_params['magic_para']=1./4.
solver_params['cphi_fact']=1./3.
ade_trt=yantra.AdvectionDiffusion(domain,domain_params,bc_params,solver_params)
#%% run models
#%% TIME SETTINGS
nitr =1000
tf =  1.0 #seconds
tf = tf + 0.001
step = tf/20.
time_points = np.arange(0, tf+step, step)
#%% RUN SOLVER
itr = 0 
l = 1
j = 0
'''
conc_step_srt = []
while  ade_srt.time <=tf: # itr < nitr: #            
    if ( (ade_srt.time <= time_points[j]) and ((ade_srt.time + ade_srt.dt) > time_points[j]) ):  
        conc_step_srt.append(ade_srt.c[1,l])            
        j +=1
    ade_srt.advance() 
    itr += 1
'''   
j = 0 
conc_step_trt = []
conc_step_trt1 = []
conc_step_trt2 = []
time = []
while   itr < nitr: #            ade_trt.time <=tf: #
          
    ade_trt.advance() 
    conc_step_trt.append(ade_trt.c[1,l])            
    conc_step_trt1.append(ade_trt.c[1,l+1])       
    conc_step_trt2.append(ade_trt.c[1,l+2])   
    time.append(ade_trt.time)
    itr += 1
#%% Compute diffusivcity SRT
''' 
print('================================')
#l = 20
x = l*dx #m
print('Length %s' %x)
cm = 1.0

from scipy.optimize import root
from math import erf
def equation(d, x, t, cm, c):
    return (c - cm*(1.-erf(x/2./np.sqrt(d*t))))

diffusivity1 = []
for i in np.arange(0, len(conc_step_srt)):
    res = root(equation, D, args = (x, time_points[i], cm, conc_step_srt[i]))
    #k = res.x[0]
    diffusivity1.append(res.x[0])
    #print('Time %s' %time_points[i])
    #print('Concentration %s' %conc_step_srt[i])
    #print('Diffusivity %s' %k)
'''
#%% Compute diffusivcity TRT
print('================================')

x = l*dx #m
print('Length %s' %x)
cm = 1.0


from scipy.optimize import root
from math import erf
def equation(d, x, t, cm, c):
    return (c - cm*(1.-erf(x/2./np.sqrt(d*t))))

diffusivity0 = []
diffusivity1 = []
diffusivity2 = []
for i in np.arange(0, len(conc_step_trt)):
    res = root(equation, D, args = (x, time[i], cm, conc_step_trt[i]))
    diffusivity0.append(res.x[0])
x = x+dx
for i in np.arange(0, len(conc_step_trt1)):
    res = root(equation, D, args = (x, time[i], cm, conc_step_trt1[i]))
    diffusivity1.append(res.x[0])
x = x+dx
for i in np.arange(0, len(conc_step_trt2)):
    res = root(equation, D, args = (x, time[i], cm, conc_step_trt2[i]))
    diffusivity2.append(res.x[0])
#%%
plt.figure()
plt.loglog(diffusivity0, label = "1")
plt.loglog(diffusivity1, label = "2")
plt.loglog(diffusivity2, label = "3")
plt.legend()
plt.show()

plt.figure()
plt.plot(diffusivity0, label = "1")
plt.plot(diffusivity1, label = "2")
plt.plot(diffusivity2, label = "3")
plt.legend()
plt.show()

print(diffusivity0[-1])
print(diffusivity1[-1])
print(diffusivity2[-1])
print(diffusivity0[0])
print(diffusivity1[0])
print(diffusivity2[0])
#%% CHECK
'''
t = 10
x = 20e-6
d = 1e-9
cm = 1.0
print(cm*(1.-erf(x/2./np.sqrt(d*t))))
#%%
x = domain.meshgrid()[0]
plt.figure()
plt.plot(x[1,:], ade_trt.c[1,:])
plt.show()
#%% CHECK
t = 100.*3600
x = 2.0e-3
d = 1.17e-10
cm = 800.
print(cm*(1.-erf(x/2./np.sqrt(d*t))))
'''