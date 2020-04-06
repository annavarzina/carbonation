# -*- coding: utf-8 -*-
'''
Sript for both CH and CSH systems
'''

#%% PYTHON MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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
Carbonation of lime paste
"""
#problem type
m = 'CH' #or 'CSH'

#%% GEOMETRY
lx = 11*1.0e-6
ly = 31*1.0e-6
dx = 1.0e-6
rad = 6*dx
#wall_len_y = wall_len_x 
domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(ly, lx), 
                         dx=dx, 
                         grid_type='nodal')

np.random.seed(0)
ra = np.random.random(domain.nodetype.shape)
ra = ra>=0.5

domain.nodetype[ra] = ct.Type.MULTILEVEL
domain.nodetype[:,0] = ct.Type.LIQUID
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()



#%%  VALUES
scale = 50
init_porosCH = 0.05
nn='01_rand_lime_pco2_34'
#fn.make_output_dir(root_dir+'\\results\\output\\simulations\\')
path = root_dir+'\\results\\output\\12_lime_paste\\' + nn + '\\'
fn.make_output_dir(path)
np.save(path + 'init_geometry', domain.nodetype)
'''
phrqc_input = {'c_bc':{'type':'pco2', 'value': 3.4}, #3.05E-02, 3.74E-02, 4.30E-02
               'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
               'c_liq':{'type':'eq', 'value': 'calcite'},
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
               'ca_liq':{'type':'eq', 'value': 'calcite'}}#calcite
'''
phrqc_input = {'c_bc':{'type':'pco2', 'value': 3.4}, #3.05E-02, 3.74E-02, 4.30E-02
               'c_mlvl':{'type':'conc', 'value': '0'}, 
               'c_liq':{'type':'eq', 'value': 'calcite'},
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
               'ca_liq':{'type':'eq', 'value': 'calcite'}}#calcite
phrqc = fn.set_phrqc_input(phrqc_input)            
fn.save_phrqc_input(phrqc,root_dir, nn)   

tfact =  1./6.*2.
default_porosCH = 0.05

mvol_ratio = 3.69/3.31
mvolCH = 0.0331*scale * (1-init_porosCH) / (1-default_porosCH)
mvol = [mvolCH, mvolCH*mvol_ratio]
mvol = fn.set_mvols(mvol, ptype = m) #m3/mol
max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = fn.set_init_pqty(mvol, init_porosCH)
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m)          
D = 1.0e-09 # default diffusion coefficient in pure liquid
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort_degree = 1./3.
app_tort = 1. * porosity ** (app_tort_degree)


settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            'dissolution':'subgrid', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/
            'diffusivity':{'border': D, ##diffusivity at border
                           'CH': ('const', 1e-15), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CC': ('const', 1e-12), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           }, 
            'pcs_mode': {'pcs': True, #Pore-Size Controlled Solubility concept
                         'pores': 'block', #'block'/'cylinder'
                         'int_energy': 0.1, # internal energy
                         'pore_size': 0.005*dx, # threshold radius or distance/2
                         'crystal_size': 0.5*dx, # crystal or pore length
                         }, 
            'subgrid': {'fraction':0.004}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': app_tort_degree}, 
            'velocity': False, 
            'bc': phrqc_input['c_bc'],
            'dx': dx, 
            'Dref':D
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
carb_rt.phrqc.phrqc_smart_run_tol = 1e-4
#%% OUTPUT PARAMETERS
plist =  []
pavglist = ['avg_poros', 'avg_D_eff']
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)
it=time.time()
#%% TIME SETTINGS
itr = 0 
j = 0
nitr = 500

Ts = 3600. #sec
Ts = Ts/scale + 0.001#1.001#1.01 +
step = max(1,int(Ts/60))
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step)))


N = Ts/carb_rt.dt
N_res = 1e+4
S = max(1,int(N/N_res))


#%% RUN SOLVER
while  carb_rt.time <=Ts: #itr < nitr: #
    if(True):
        if ( (carb_rt.time <= time_points[j]) and ((carb_rt.time + carb_rt.dt) > time_points[j]) ):  
            print('Time %s' %time_points[j])
            if(carb_rt.time >0):
                print('Active cells %s'%carb_rt.phrqc.nactive)
            fn.save_figures_minerals(carb_rt,  max_pqty, time_points[j], path, nn, ptype=m)  
            fn.save_figures_mols(carb_rt, time_points[j], path, nn, ptype=m) 
            if(False):
                print('SI %carb_rt' %carb_rt.phrqc.selected_output()['SI_calcite'][1,:])
                print('C %carb_rt' %carb_rt.fluid.C.c[1,:])
                print('Ca %carb_rt' %carb_rt.fluid.Ca.c[1,:])
            j +=1
    carb_rt.advance() 
    results = fn.append_results(carb_rt, results, step=S)
    itr += 1
    
#%% SIMULATION TIME
simulation_time = time.time()-it
fn.print_time(simulation_time, carb_rt)

#%%  SAVE
fn.save_obj(results, path + str(nn) +'_results')
np.save(path + 'SI', carb_rt.phrqc.selected_output()['SI_calcite'] )
np.save(path + 'pH', carb_rt.phrqc.selected_output()['pH'] )
np.save(path + 'Ca', carb_rt.phrqc.selected_output()['Ca'] )
np.save(path + 'C', carb_rt.phrqc.selected_output()['C'] )
np.save(path + 'De', carb_rt.fluid.Ca.De )
np.save(path + 'poros', carb_rt.fluid.Ca.poros)
#%% PLOT
fn.plot_species(results, names=['calcite', 'portlandite', 'Ca', 'C'],fsize=(6,4))
fn.plot_fields(carb_rt, names={ 'calcite', 'portlandite', 'Ca', 'C'})

#%%PRINT

#print(carb_rt.fluid.Ca.De)
#print(carb_rt.fluid.C.De)
#print(carb_rt.fluid.Ca._c+carb_rt.fluid.Ca._ss)
#print(carb_rt.fluid.C._c+carb_rt.fluid.C._ss)
#print(carb_rt.solid.portlandite.c)
#print(carb_rt.solid.calcite.c)
#print(carb_rt.solid.poros)


#%%
'''
from copy import deepcopy

print(carb_rt.fluid.C._c)
c = deepcopy(carb_rt.fluid.get_attr('_c'))
f = deepcopy(carb_rt.fluid.get_attr('_f'))
f2 = deepcopy(carb_rt.fluid.C._f)
flag = False
comp = 'C'
n = c[comp] < 0
if n.any():
    flag = True
    f[comp][n] = f[comp][n] - f[comp][n]*(f[comp][n]<0)
carb_rt.fluid.set_attr('_f',f)
carb_rt.fluid.set_attr('f',f)
carb_rt.fluid.call('compute_macro_var')
print(carb_rt.fluid.C._c)
'''