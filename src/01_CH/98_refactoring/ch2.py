# -*- coding: utf-8 -*-
from __future__ import division  
#%% SYS and OS
import sys,os
root_dir = os.path.dirname(
        os.path.dirname(
                os.path.dirname(
                        os.path.dirname(
                                os.path.abspath(__file__)
                        ) ) ) )
src_dir = os.path.dirname(
        os.path.dirname(
                os.path.dirname(
                        os.path.abspath(__file__)
                        ) ) )
ch_dir = os.path.dirname(
        os.path.dirname(
                os.path.abspath(__file__)
                ) )

sys.path.append(root_dir)
sys.path.append(src_dir)
sys.path.append(ch_dir)

#%% PYTHON MODULES
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import time
import yantra
import cell_type as ct # change the path to cell_type file
import misc_func as fn
import rt_carb_ch as rt
#import phrqc
#%% PROBLEM DEFINITION
__doc__= """ 
Default example. Explains possible values
1D carbonation of portlandite/cement.
"""
#problem type
m = 'CH' 

#%% GEOMETRY
ll = 5#liquid lauer in front of portlandite
l_ch = 6#length of portlandite
lx = (l_ch+ll)*1.0e-6
ly = 2.0e-6
dx = 1.0e-6

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, ll+1: ll+l_ch] = ct.Type.MULTILEVEL
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()

#%%  PHREEQC
nn='01_example_default'
fn.make_output_dir(root_dir+'\\results\\output\\00_examples\\')
path = root_dir+'\\results\\output\\00_examples\\' + nn + '\\'
fn.make_output_dir(path)
phrqc_input = {'c_bc':{'type':'pco2', 'value': 3.4},
               'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
               'c_liq':{'type':'eq', 'value': 'calcite'},
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'},    
               'ca_liq':{'type':'eq', 'value': 'calcite'}}

phrqc = rt.PhreeqcInputCH(phrqc_input)            
phrqc.save_phrqc_input(root_dir, nn)   

#%% VALUES
tfact_default = 1./6./1#*init_porosCH
scale = 50 # scale of molar volume
init_poros = [0.05, 1.0] #initial porosity of portlandite or calcite nodes
D = 1.0e-09 # default diffusion coefficient in pure liquid
app_tort_degree = 1./3.
settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            'dissolution':'subgrid', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/
            'diffusivity':{'border': D, ##diffusivity at border
                           'CH': ('const', 1e-15), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CC': ('inverse', 1e-12), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           }, 
            'pcs_mode': {'pcs': True, #Pore-Size Controlled Solubility concept
                         'pores': 'block', #'block'/'cylinder'
                         'int_energy': 0.5, # internal energy
                         'pore_size': 0.01*dx, # threshold radius or distance/2
                         'crystal_size': 0.5*dx, # crystal or pore length
                         }, 
            'subgrid': {'fraction':0.004}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': app_tort_degree}, 
            'velocity': False, 
            'bc': phrqc_input['c_bc'],
            'dx': dx, 
            'Dref':D}
            
mvol = rt.SettingsCH.set_mvols(scale) #dm3/mol
max_pqty = rt.SettingsCH.get_max_pqty(mvol) #mol/dm3
init_conc = rt.SettingsCH.set_init_pqty(mvol, init_poros, scale)
pqty = rt.SettingsCH.get_pqty(init_conc, domain)
slabels = rt.SettingsCH.set_labels(domain)          
porosity = rt.SettingsCH.get_porosity(domain, pqty, mvol)
app_tort =  rt.SettingsCH.get_app_tort(domain, app_tort_degree, porosity)
            
#%% PARAMETERS (DOMAIN, BC, SOLVER)
domain_params = rt.SettingsCH.set_domain_params(D, mvol, pqty, porosity, 
                    app_tort, slabels, database = 'cemdata07.dat',
                    input_file = root_dir +'\\phreeqc_input\\' + nn + '.phrq')
bc_params = rt.SettingsCH.set_bc_params(bc_slabels = {'left':100001})
solver_params = rt.SettingsCH.set_solver_params(tfact = tfact_default, smart_thres = 1e-8, cphi_fact = 1/3.)
domain = rt.SettingsCH.reset_portlandite_nodes(domain)
rt.SettingsCH.save_settings(settings, bc_params, solver_params, path, nn)

#%% INITIATE THE SOLVER
carb_rt= rt.CarbonationCH('MultilevelAdvectionDiffusion',  domain, 
                          domain_params, bc_params, solver_params, settings) 
#%% PARAMETERS
#results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)
res = rt.ResultsCH(nodes= [(1,n) for n in np.arange(1, 8)])
#%% TIME SETTINGS
nitr =1000
Ts =  1#3600.*3 #seconds
Ts = Ts/scale + 0.001
step = max(int(Ts/10.),1)
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step))) #time_points = np.arange(0, Ts+step, step)
N = Ts/carb_rt.dt
N_res = 1e+4
S = max(1,int(N/N_res))
it=time.time()

#%% RUN SOLVER
itr = 0 
j = 0
while carb_rt.time <=Ts: #itr < nitr: # 
    if(True):
        if ( (carb_rt.time <= time_points[j]) and \
            ((carb_rt.time + carb_rt.dt) > time_points[j]) ): 
            print(time_points[j])
            j +=1
        
    carb_rt.advance()    
    res.append_results(carb_rt, step = S)
    itr += 1
    
#%% SIMULATION TIME
simulation_time = time.time()-it
rt.ResultsCH.print_time(simulation_time, carb_rt)
  
#%% PLOT 

#res.plot_species(names=[])#['calcite']
#res.plot_avg()
#res.plot_nodes(names=[])
#res.plot_fields(carb_rt, names=[],fsize=(15,1))

#%% PRINT
rt.ResultsCH.print_profiles(carb_rt)