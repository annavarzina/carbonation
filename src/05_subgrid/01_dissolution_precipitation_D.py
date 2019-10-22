# -*- coding: utf-8 -*-
'''
Example with precipitation everywhere
Fixed PCO2 at the boundary
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
Reference:
    Poros = 0.05
    PCO2 = 3.4
    IE = 0.5
    Archies relation for diffusivity
"""
#problem type
m = 'CH' #or 'CSH'

ll = 1
l = 5 +ll
lx = l*1.0e-6
ly = 2.0e-6
dx = 1.0e-6

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, ll+1:l+ll] = ct.Type.MULTILEVEL

domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()
#%%  VALUES
nn='01_test_D_p005'
scale = 50
init_porosCH = 0.05
fn.make_output_dir(root_dir+'\\results\\temp\\01_subgrid\\')
path = root_dir+'\\results\\temp\\01_subgrid\\' + nn + '\\'
fn.make_output_dir(path)

phrqc_input = {'c_bc':{'type':'conc', 'value': 1e-2}, #3.05E-02, 3.74E-02, 4.30E-02
               'c_mlvl':{'type':'conc', 'value': '0'}, 
               'c_liq':{'type':'conc', 'value': '0'},
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
               'ca_liq':{'type':'conc', 'value': '0'}}#calcite
phrqc = fn.set_phrqc_input(phrqc_input)            
fn.save_phrqc_input(phrqc,root_dir, nn)   

tfact =  1./6.
mvol_ratio = 3.69/3.31
mvolCH = 0.0331*scale
mvol = [mvolCH, mvolCH*mvol_ratio]

mvol = fn.set_mvols(mvol, ptype = m) #m3/mol
max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = fn.set_init_pqty(mvol, init_porosCH)
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m)          
D = 1.0e-09 # default diffusion coefficient in pure liquid
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort = 1. * porosity ** (1./3.)

settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            #'dissolution':'subgrid',
            'active': 'all', # 'all'/'smart'/'interface'
            'diffusivity':{'type':'cc_archie_ch_kin', #'fixed' or 'archie'
                           'D_CH': 1e-12,
                           'D_CH_Ca': 1e-9},
            'pcs': {'pcs': True, 
                    'pores': 'block', #'block'/'cylinder'
                    'int_energy': 0.5, # internal energy
                    'pore_size': 0.01*dx, # threshold radius or distance/2
                    'crystal_size': 0.5*dx, # crystal or pore length
                    'pore_density': 20000, #pore density per um3 - only for cylinder type
                    }, 
           'velocity': False, 
           'bc': phrqc_input['c_bc'],
           'dx': dx 
           }
               
            
#%% PARAMETERS (DOMAIN, BC, SOLVER)
domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir +'\\phreeqc_input\\' + nn + '.phrq')
bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = tfact)
domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL
fn.save_settings(settings, bc_params, solver_params, path, nn)

#%% INITIATE THE SOLVER
carb_rt= rt.CarbonationRT('MultilevelAdvectionDiffusion',  domain, 
                          domain_params, bc_params, solver_params,
                          settings) 
carb_rt.phrqc.phrqc_smart_run_tol = 1e-6
#%% PARAMETERS
#plist =  [(1,2), (1,3), (1,4), (1,5), (1,6), (1,7), (1,8), (1,9), (1,10)]
plist =  [(1,n) for n in np.arange(1,4)]
pavglist = ['avg_poros', 'pH', 'avg_D_eff', 'sum_vol', 'precipitation',
            'dissolution', 'portlandite_cells', 'calcite_cells', 'dt'] 
#'delta_ch', 'delta_cc', 'precipitation','dissolution', 'portlandite_cells', 
#'calcite_cells', 'active_cells','dt', 'pH', 'avg_poros',  'avg_D_eff', 'sum_vol'
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)


#%% TIME SETTINGS
itr = 0 
j = 0
ni = 200
nitr =500
Ts = 1 #second
Ts = Ts/scale + 0.001#1.001#1.01 +
step = max(int(Ts/36.),1)
#time_points = np.arange(0, Ts+step, step)
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step)))

N = Ts/carb_rt.dt
N_res = 1e+4
S = max(1,int(N/N_res))
#%% RUN SOLVER
while  itr <= nitr: #carb_rt.time <=Ts: #
    
    carb_rt.advance()    
    results = fn.append_results(carb_rt, results, step = S )
    itr += 1
    
fn.save_obj(results, path + str(nn) +'_results')
#%% SIMULATION TIME

print('Ca ss %s' %str(np.array(carb_rt.fluid.Ca._ss[1,:])))
print('Ca +ss %s' %str(np.array(carb_rt.fluid.Ca.c[1,:]) + np.array(carb_rt.fluid.Ca._ss[1,:])/np.array(carb_rt.phrqc.poros[1,:])))
print('Ca +ss %s' %str(np.array(carb_rt.fluid.Ca.c[1,:]) + np.array(carb_rt.fluid.Ca._ss[1,:])))
print('C +ss %s' %str(np.array(carb_rt.fluid.C.c[1,:]) + np.array(carb_rt.fluid.C._ss[1,:])/np.array(carb_rt.phrqc.poros[1,:])))
print('H +ss %s' %str(np.array(carb_rt.fluid.H.c[1,:]) + np.array(carb_rt.fluid.H._ss[1,:])/np.array(carb_rt.phrqc.poros[1,:])))
print('O +ss %s' %str(np.array(carb_rt.fluid.O.c[1,:]) + np.array(carb_rt.fluid.O._ss[1,:])/np.array(carb_rt.phrqc.poros[1,:])))
print('CH %s' %str(np.array(carb_rt.solid.portlandite.c[1,:])))
print('dCH %s' %str(np.array(carb_rt.phrqc.dphases['portlandite'][1,:])))
print('CC %s' %str(np.array(carb_rt.solid.calcite.c[1,:])))
print('dCC %s' %str(np.array(carb_rt.phrqc.dphases['calcite'][1,:])))
print('Vol %s' %str(np.array(carb_rt.solid.vol[1,:])))

#print(carb_rt.phrqc.selected_output())

#%%
fn.plot_points(results,names={ 'calcite', 'portlandite', 'Ca'})#, 'C', 'O', 'H'})
#fn.plot_species(results, names={ 'calcite', 'portlandite', 'Ca', 'C', 'O', 'H'})
#print(results['portlandite'][-1]-results['portlandite'][0])
#print(results['calcite'][-1]-results['calcite'][0])
print(results['portlandite (1, 2)'][-1]-results['portlandite (1, 2)'][0])
print(results['portlandite (1, 3)'][-1]-results['portlandite (1, 3)'][0])
print(results['calcite (1, 2)'][-1]-results['calcite (1, 2)'][0])
print(results['calcite (1, 1)'][-1]-results['calcite (1, 1)'][0])