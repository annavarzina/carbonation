# -*- coding: utf-8 -*-
'''
Example with precipitation everywhere
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
import rt 
#import phrqc
#%% PROBLEM DEFINITION
__doc__= """ 
Carbonation of cement. 
Portlandite only.

"""
#problem type
m = 'CH' #or 'CSH'

#%% GEOMETRY
l = 25
lx = l*1.0e-6
ly = 2.0e-6
dx = 1.0e-6
rad = 6*dx
#wall_len_y = wall_len_x 

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, 5:l] = ct.Type.MULTILEVEL

domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()
#%%  VALUES

init_porosCH = 0.05

mvol_ratio = 3.69/3.31
mvolCH = 20
mvol = [mvolCH, mvolCH*mvol_ratio]

mvol = fn.set_mvols(mvol, ptype = m) #m3/mol
max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = fn.set_init_pqty(mvol, init_porosCH)
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m)          
D = 1.0e-09 # default diffusion coefficient in pure liquid
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort = 1. * porosity ** (1./3.)

settings = {'precipitation': 'mineral', # 'interface'/'all'/'mineral' nodes
            'active': 'all', # 'all'/'smart'/'interface'
            'diffusivity':{'type':'fixed', #'fixed' or 'archie'
                           'D_CC': 9e-12,
                           'D_CH': 1e-12},
            'pcs': {'pcs': True, 
                    'pores': 'block', #'block'/'cylinder'
                    'int_energy': 0.485, # internal energy
                    'pore_size': 0.01*dx, # threshold radius or distance/2
                    'crystal_size': 0.5*dx, # crystal or pore length
                    'pore_density': 20000, #pore density per um3 - only for cylinder type
                    #'threshold': 'poresize', #poresize/porosity or si
                    #'threshold_value': 1.0, 
                    }, 
           'velocity': False, 
           #'bc': 'const',
           'dx': dx 
           }
               
tfact =  1./6.   

nn='example_1def'#'acc10'
path = root_dir+'\\results\\output\\'

#%% PARAMETERS (DOMAIN, BC, SOLVER)
domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir +'\\phreeqc_input\\CH_CC_nat.phrq')#'CH_CC-1percent.phrq'
                                     #input_file = 'CH_CC_-2.phrq')
bc_params={'solution_labels':{'left':100001}, 
           'top':['flux', 0.0],
           'bottom':['flux', 0.0],
           'left':['flux', 0.0],
           'right':['flux', 0.0],}
solver_params = fn.set_solver_params(tfact = tfact)
domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL

#%% INITIATE THE SOLVER
carb_rt= rt.CarbonationRT('MultilevelAdvectionDiffusion',  domain, domain_params, bc_params, solver_params) 
carb_rt.settings = settings
fn.apply_settings(carb_rt)
fn.save_settings(carb_rt.settings, bc_params, solver_params, path, nn)

#%% PARAMETERS
plist =  [(1,2), (1,3), (1,4), (1,5), (1,6), (1,7), (1,8), (1,9), (1,10)]#[(1,n) for n in np.array([1, 2, 3])] #v
pavglist = ['avg_poros', 'avg_D_eff']
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)

#%% TIME SETTINGS
itr = 0 
j = 0
ni = 100
nitr = 20
Ts = 1.11#1.001#1.01
step = 0.1
#time_points = np.arange(0, Ts+step, step)
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step)))
it=time.time()

#%% RUN SOLVER
while carb_rt.time <=Ts: #itr < nitr: # 
    if(False):
        if ( (carb_rt.time <= time_points[j]) and ((carb_rt.time + carb_rt.dt) > time_points[j]) ):  
            print(time_points[j])
            #fn.save_figures_minerals(rt,  max_pqty, time_points[j], path, nn, ptype=m)  
            #save_figures_mols(rt, time_points[j], path, nn, ptype=m) 
            #save_vti(rt, phases, time_points[j], path, nn, m)
            #save_pickle(rt, phases, time_points[j], path, nn)
            if(j>0):
                points = [(1,n) for n in np.arange(1,15)]
                fn.print_points(carb_rt, points, names=['calcite', 'portlandite'])
                print('SI %carb_rt' %carb_rt.phrqc.selected_output()['SI_calcite'][1,:])
                print('C %carb_rt' %carb_rt.fluid.C.c[1,:])
                print('Ca %carb_rt' %carb_rt.fluid.Ca.c[1,:])
            j +=1
        
    carb_rt.advance()    
    results = fn.append_results(carb_rt, results)
    itr += 1
    
#%% SIMULATION TIME
simulation_time = time.time()-it
fn.print_time(simulation_time, carb_rt)
            
#%%  SAVE
fresults  = fn.filter_results(results, path, nn)
#fn.save_obj(fresults, path + str(nn) +'_results')

#%% PLOT 
fn.plot_species(results, names=[])#['calcite']
fn.plot_avg(results, names=['avg_poros', 'avg_D_eff'])
fn.plot_points(results, names=['calcite', 'portlandite', 'poros', 'Ca', 'C'])
fn.plot_fields(carb_rt, names=['calcite', 'Ca', 'poros'],fsize=(15,1))

#%% PRINT
points = [(1,n) for n in np.arange(2,15)]
#fn.print_points(rt, points)