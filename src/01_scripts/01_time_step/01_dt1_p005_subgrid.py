# -*- coding: utf-8 -*-
'''
Example with precipitation at interface
Fixed PCO2 at the boundary
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
Default example. Explains possible values
1D carbonation of portlandite/cement.
"""
#problem type
m = 'CH' #or 'CSH' #TODO case for cement

#%% GEOMETRY
ll = 1 #liquid lauer in front of portlandite
l_ch = 25 #length of portlandite
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
#%%  VALUES
scale = 50
nn=os.path.basename(__file__)[:-3]
fn.make_output_dir(root_dir+'\\results\\output\\01_time_step\\')
path = root_dir+'\\results\\output\\01_time_step\\' + nn + '\\'
fn.make_output_dir(path)

phrqc_input = {'c_bc':{'type':'pco2', 'value': 3.4}, #3.05E-02, 3.74E-02, 4.30E-02
               'c_mlvl':{'type':'conc', 'value': '0'}, 
               'c_liq':{'type':'conc', 'value': '0'},
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
               'ca_liq':{'type':'eq', 'value': 'portlandite'}}#calcite
phrqc = fn.set_phrqc_input(phrqc_input)            
fn.save_phrqc_input(phrqc,root_dir, nn)   

scale = 50 # scale of molar volume
init_porosCH = 0.05 #initial porosity of portlandite nodes
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
app_tort_degree = 1./3.
app_tort = 1. * porosity ** app_tort_degree

settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            'dissolution':'subgrid', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/'interface'
            'diffusivity':{'type':'archie', #'archie'/'fixed'/'mixed'
                           #'D_border':D, #diffusivity at border
                           #'D_CC': 9e-12, # fixed diffusivity in calcite node
                           #'D_CH': 1e-12, # fixed diffusivity in portlandite node
                           }, 
            'pcs_mode': {'pcs': True, #Pore-Size Controlled Solubility concept
                         'pores': 'block', #'block'/'cylinder'
                         'int_energy': 0.5, # internal energy
                         'pore_size': 0.01*dx, # threshold radius or distance/2
                         'crystal_size': 0.5*dx, # crystal or pore length
                         'pore_density': 2000, #pore density per um3 - only for cylinder type
                         }, 
            'subgrid': {'fraction':None}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': app_tort_degree}, #TODO
            'velocity': False, 
            'bc': phrqc_input['c_bc'],
            'dx': dx, 
            'Dref':D
            }
 
tfact_default = 1./6.*init_porosCH
tfact_scale = 1.
tfact = tfact_default * tfact_scale
#%% PARAMETERS (DOMAIN, BC, SOLVER)
domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir +'\\phreeqc_input\\' + nn + '.phrq')
bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = None, smart_thres = 1e-8)
domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL
fn.save_settings(settings, bc_params, solver_params, path, nn)

#%% INITIATE THE SOLVER
carb_rt= rt.CarbonationRT('MultilevelAdvectionDiffusion',  domain, 
                          domain_params, bc_params, solver_params,
                          settings) 

#%% PARAMETERS
#plist =  [(1,2), (1,3), (1,4), (1,5), (1,6), (1,7), (1,8), (1,9), (1,10)]
plist =  [(1,n) for n in np.arange(0, 6)]
pavglist = ['avg_poros', 'pH', 'avg_D_eff', 'sum_vol', 'precipitation',
            'dissolution', 'portlandite_cells', 'calcite_cells'] 
#'delta_ch', 'delta_cc', 'precipitation','dissolution', 'portlandite_cells', 
#'calcite_cells', 'active_cells','dt', 'pH', 'avg_poros',  'avg_D_eff', 'sum_vol'
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)

#%% TIME SETTINGS
itr = 0 
j = 0
ni = 100
nitr = 100
Ts = 1000.
Ts = Ts/scale + 0.001#1.001#1.01
step = max(int(Ts/10.),1)
#time_points = np.arange(0, Ts+step, step)
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step)))
it=time.time()

N = Ts/carb_rt.dt
N_res = 1e+4
S = max(1,int(N/N_res))
#%% RUN SOLVER
while  carb_rt.time <=Ts: #itr <= nitr: #
    if(True):
        if ( (carb_rt.time <= time_points[j]) and ((carb_rt.time + carb_rt.dt) > time_points[j]) ):  
            print(time_points[j])
            fn.save_figures_minerals(carb_rt,  max_pqty, time_points[j], path, nn, ptype=m)  
            fn.save_figures_mols(carb_rt, time_points[j], path, nn, ptype=m) 
            j +=1
        
    carb_rt.advance()    
    results = fn.append_results(carb_rt, results, step = S )
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

fn.plot_species(results, names=[])#['calcite']
fn.plot_avg(results, names=['avg_poros', 'avg_D_eff'])
fn.plot_points(results, names=['calcite', 'portlandite', 'poros', 'Ca', 'C'])
fn.plot_fields(carb_rt, names=['calcite', 'portlandite', 'Ca', 'C', 'poros'],fsize=(15,1))

#%% PRINT
print('Ca %s' %str(np.array(carb_rt.fluid.Ca._c[1,:])))
print('Ca +ss %s' %str(np.array(carb_rt.fluid.Ca.c[1,:]) + np.array(carb_rt.fluid.Ca._ss[1,:])/np.array(carb_rt.phrqc.poros[1,:])))
print('C %s' %str(np.array(carb_rt.fluid.C._c[1,:])))
print('C +ss %s' %str(np.array(carb_rt.fluid.C.c[1,:]) + np.array(carb_rt.fluid.C._ss[1,:])/np.array(carb_rt.phrqc.poros[1,:])))
print('H +ss %s' %str(np.array(carb_rt.fluid.H.c[1,:]) + np.array(carb_rt.fluid.H._ss[1,:])/np.array(carb_rt.phrqc.poros[1,:])))
print('O +ss %s' %str(np.array(carb_rt.fluid.O.c[1,:]) + np.array(carb_rt.fluid.O._ss[1,:])/np.array(carb_rt.phrqc.poros[1,:])))
print('CH %s' %str(np.array(carb_rt.solid.portlandite.c[1,:])))
print('dCH %s' %str(np.array(carb_rt.phrqc.dphases['portlandite'][1,:])))
print('CC %s' %str(np.array(carb_rt.solid.calcite.c[1,:])))
print('dCC %s' %str(np.array(carb_rt.phrqc.dphases['calcite'][1,:])))
print('Vol %s' %str(np.array(carb_rt.solid.vol[1,:])))
print('D %s' %str(np.array(carb_rt.fluid.C.De[1,:])))
print('SI %s' %str(np.array(carb_rt.solid.target_SI[1,:])))
print('pH %s' %str(np.array(carb_rt.phrqc.selected_output()['pH'][1,:])))
print('poros %s' %str(np.array(carb_rt.solid.poros[1,:])))
print('phrqc poros %s' %str(np.array(np.array(carb_rt.phrqc.poros[1,:]))))


print('Total CH dissolved %s' %(results['portlandite'][-1]-results['portlandite'][0]))
print('Total CC precipitated %s' %(results['calcite'][-1]-results['calcite'][0]))