# -*- coding: utf-8 -*-
'''
Sript for both CH and CSH systems
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
domain.nodetype[:, 1:l] = ct.Type.MULTILEVEL

domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()

#%%  VALUES
nn='example_11_2'#'acc10'
path = root_dir+'\\results\\output\\'

phrqc_input = {'c_bc':{'type':'pco2', 'value': 1.0}, 
#phrqc_input = {'c_bc':{'type':'conc', 'value': 0.01}, 
               'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
               'c_liq':{'type':'eq', 'value': 'calcite'},
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
               'ca_liq':{'type':'eq', 'value': 'calcite'}}
phrqc = fn.set_phrqc_input(phrqc_input)            
fn.save_phrqc_input(phrqc,root_dir, nn)   

tfact =  1./6./2

init_porosCH = 0.25
mvol_ratio = 3.69/3.31 #m3/mol
mvolCH = 3.31 #l/mol
mvol = [mvolCH, mvolCH*mvol_ratio]

mvol = fn.set_mvols(mvol, ptype = m) #l/mol
max_pqty = fn.get_max_pqty(mvol) #mol/l
init_conc = fn.set_init_pqty(mvol, init_porosCH) #mol/l
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m)          
D = 1.0e-09 # default diffusion coefficient in pure liquid
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort = 1. * porosity ** (1./3.)

settings = {'precipitation': 'all', # 'interface'/'all'/'mineral' nodes
            'active': 'all', # 'all'/'smart'/'interface'
            'diffusivity':{'type':'archie', #'fixed' or 'archie'
                           'D_CC': 5e-11,
                           'D_CH': 1e-11},
            'pcs': {'pcs': True, 
                    'pores': 'block', #'block'/'cylinder'
                    'int_energy': 0.485, # internal energy
                    'pore_size': 0.005*dx, # threshold radius or distance/2
                    'crystal_size': 0.5*dx, # crystal or pore length
                    'pore_density': 20000, #pore density per um3 - only for cylinder type
                    }, 
           'velocity': False, 
           'bc': phrqc_input['c_bc'],
           'pco2': True,
           'dx': dx 
           }
               

#%% PARAMETERS (DOMAIN, BC, SOLVER)
domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir +'\\phreeqc_input\\' + nn + '.phrq')#'CH_CC-1percent.phrq'
                                     #input_file = 'CH_CC_-2.phrq')
bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = tfact)
domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL
fn.save_settings(settings, bc_params, solver_params, path, nn)

#%% INITIATE THE SOLVER
carb_rt= rt.CarbonationRT('MultilevelAdvectionDiffusion',  domain, domain_params, bc_params, solver_params, settings) 

#%% OUTPUT PARAMETERS
plist = [(1,n) for n in np.arange(0, l)]
pavglist = ['avg_poros', 'pH', 'avg_D_eff', 'sum_vol', 'precipitation',
            'dissolution', 'portlandite_cells', 'calcite_cells'] 
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)

#%% TIME SETTINGS
itr = 0 
j = 0
ni = 100
nitr = 2#20
Ts = 5.01#51#1.001#1.01
step = 1.0
#time_points = np.arange(0, Ts+step, step)
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step)))
it=time.time()

#%% RUN SOLVER
while carb_rt.time <=Ts: #itr < nitr: # 
    if(True):
        if ( (carb_rt.time <= time_points[j]) and ((carb_rt.time + carb_rt.dt) > time_points[j]) ):  
            print(time_points[j])
            fn.save_figures_minerals(carb_rt,  max_pqty, time_points[j], path+'fig\\', nn, ptype=m)  
            fn.save_figures_mols(carb_rt, time_points[j], path+'fig\\', nn, ptype=m, cC = 0.1, cCa = 0.05)
            #fn.save_figures_minerals(rt,  max_pqty, time_points[j], path, nn, ptype=m)  
            #save_figures_mols(rt, time_points[j], path, nn, ptype=m) 
            #save_vti(rt, phases, time_points[j], path, nn, m)
            #save_pickle(rt, phases, time_points[j], path, nn)
            if(j>0):                
                fn.plot_fields(carb_rt, names=['calcite', 'C', 'Ca', 'poros'],fsize=(15,1))
            j +=1
        
    carb_rt.advance()    
    results = fn.append_results(carb_rt, results)
    itr += 1
    
#%% SIMULATION TIME
simulation_time = time.time()-it
fn.print_time(simulation_time, carb_rt)
            
#%%  SAVE
fresults  = fn.filter_results(results, path, nn)
fn.save_obj(fresults, path + str(nn) +'_results')

#%% PLOT 
fn.plot_species(results, names=['calcite', 'portlandite'])#['calcite']
fn.plot_avg(results, names=['avg_poros', 'avg_D_eff'])
fn.plot_points(results, names=['calcite', 'portlandite', 'poros', 'Ca', 'C'],fsize=(20,8))
fn.plot_fields(carb_rt, names=['calcite', 'Ca', 'poros', 'C', 'portlandite'],fsize=(15,1))

#%% PRINT
#points = [(1,n) for n in np.arange(2,15)]
#fn.print_points(rt, points)