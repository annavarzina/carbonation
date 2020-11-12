# -*- coding: utf-8 -*-
'''
Molar volume = 100
Liquid layer = 5 um
Fraction = 0.004
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
import yantra
import cell_type as ct # change the path to cell_type file
import  rt 
import misc_func as fn
import time
#%% PROBLEM DEFINITION
__doc__= """ 
...
"""
#%% GEOMETRY
ll = 4 #liquid lauer in front of portlandite
l_ch = 5 #length of CSH
lx = (l_ch+ll+1)*1.0e-6
ly = 2.0e-6
dx = 1.0e-6

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, ll+1: ll+l_ch+1] = ct.Type.MULTILEVEL
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()

#%%  VALUES
m = 'CSH' #or 'CSH'
nn=os.path.basename(__file__)[:-3]
fn.make_output_dir(root_dir+'\\results\\output_csh\\01_default\\')
path = root_dir+'\\results\\output_csh\\01_default\\' + nn + '\\'
fn.make_output_dir(path)


phrqc_input = {'c_bc':{'type':'pco2', 'value': 3.0},
               #'c_bc':{'type':'conc', 'value': 0.01},
               'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
               'c_liq':{'type':'eq', 'value': 'calcite'},
               'ca_mlvl':{'type':'eq', 'value': 'calcite'},    
               'ca_liq':{'type':'eq', 'value': 'calcite'}}

phrqc = fn.set_phrqc_input(phrqc_input, ptype = m)            
fn.save_phrqc_input(phrqc,root_dir, nn)  
#%%
scale = 50# scale of molar volume
init_porosCH = 0.05 #initial porosity of portlandite nodes
mvol = fn.set_mvols({}, scale, ptype = m) #m3/mol

max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = fn.set_init_pqty(mvol, scale, init_porosCH)
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m)     
D = 1.0e-09 # default diffusion coefficient in pure liquid
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort_degree = 1./3.
app_tort = 1. * porosity ** app_tort_degree

settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            'dissolution':'multilevel', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/ better use all to avoid errors
            'diffusivity':{'border': D, ##diffusivity at border
                           'CH': ('const', 1e-15), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CC': ('inverse', 1e-12), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CSH': ('const', 1e-11), # fixed diffusivity in CSH node 'archie'/'const'/'inverse'
                          }, 
            'pcs_mode': {'pcs': True, #Pore-Size Controlled Solubility concept
                         'pores': 'block', #'block'/'cylinder'
                         'int_energy': 0.1, # internal energy
                         'pore_size': 0.01*dx, # threshold radius or distance/2
                         'crystal_size': 0.5*dx, # crystal or pore length
                         'pore_density': 2000, #pore density per um3 - only for cylinder type
                         }, 
            'subgrid': {'fraction':0.004}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': 1./3.}, #TODO
            'velocity': False, 
            'bc': phrqc_input['c_bc'],
            'dx': dx, 
            'Dref':D
            }
tfact_default = 1./6./1#*init_porosCH
domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir +'\\phreeqc_input\\' + nn + '.phrq')
bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = tfact_default, smart_thres = 1e-8, cphi_fact = 1/3.)
#fn.save_settings(settings, bc_params, solver_params, path, nn)
#csh=yantra.PhrqcReactiveTransport('MultilevelDiffusion',domain,domain_params,bc_params,solver_params)
csh=rt.CarbonationRT('MultilevelDiffusion',domain,domain_params,bc_params,solver_params, settings)
#%% results dict

plist =  [(1,n) for n in np.arange(3, 6)] #points list
pavglist = ['avg_poros', 'pH', 'avg_D_eff', 'sum_vol', 'precipitation', #argument list
            'dissolution', 'portlandite_cells', 'calcite_cells'] 
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)

Ts  = 3600 * 3 #s
Ts = Ts/scale + 0.001
N = Ts/csh.dt
N_res = 1e+4
S = max(1,int(N/N_res))
'''
dt = csh.dt
Ts =  1.#seconds
Ts = Ts/scale + 0.001
N_res = 1e+4
S = max(1,int(N/N_res))
'''
#%% run
#n=1000

it=time.time()
while csh.time <=Ts: #itr < nitr: # 
    csh.advance()    
    results = fn.append_results(csh, results, step = S )
    
simulation_time = time.time()-it
fn.print_time(simulation_time, csh)
#%%
fn.plot_species(results, names=[])#['calcite']
fn.plot_avg(results, names=['avg_poros', 'avg_D_eff'])
fn.plot_fields(csh, names=['calcite'],fsize=(15,1)) #,'Ca','Si'
fn.plot_points(results, names=['calcite', 'poros', 'Ca','pH'])

#%% plot ca/si against density
plt.figure()
plt.plot(results['Ca_Si'], results['csh_density'])
plt.legend()
plt.ylabel('CSH density')
plt.xlabel('C/S')
plt.show()


#%% PRINT
'''
print('Ca %s' %str(np.array(csh.fluid.Ca._c[1,:])))
print('Ca +ss %s' %str(np.array(csh.fluid.Ca.c[1,:]) + np.array(csh.fluid.Ca._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('C %s' %str(np.array(csh.fluid.C._c[1,:])))
print('C +ss %s' %str(np.array(csh.fluid.C.c[1,:]) + np.array(csh.fluid.C._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('H +ss %s' %str(np.array(csh.fluid.H.c[1,:]) + np.array(csh.fluid.H._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('O +ss %s' %str(np.array(csh.fluid.O.c[1,:]) + np.array(csh.fluid.O._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('CH %s' %str(np.array(csh.solid.portlandite.c[1,:])))
print('dCH %s' %str(np.array(csh.phrqc.dphases['portlandite'][1,:])))
print('CC %s' %str(np.array(csh.solid.calcite.c[1,:])))
print('dCC %s' %str(np.array(csh.phrqc.dphases['calcite'][1,:])))
print('Vol %s' %str(np.array(csh.solid.vol[1,:])))
print('D %s' %str(np.array(csh.fluid.C.De[1,:])))
print('SI %s' %str(np.array(csh.solid.target_SI[1,:])))
print('pH %s' %str(np.array(csh.phrqc.selected_output()['pH'][1,:])))
print('poros %s' %str(np.array(csh.solid.poros[1,:])))
print('phrqc poros %s' %str(np.array(np.array(csh.phrqc.poros[1,:]))))
'''
#%%  SAVE
#'''
fn.save_obj(results, path + str(nn) +'_results')
np.save(path + 'SI', csh.phrqc.selected_output()['SI_calcite'] )
np.save(path + 'pH', csh.phrqc.selected_output()['pH'] )
np.save(path + 'Ca', csh.phrqc.selected_output()['Ca'] )
np.save(path + 'C', csh.phrqc.selected_output()['C'] )
np.save(path + 'CH', csh.phrqc.selected_output()['portlandite'] )
np.save(path + 'CC', csh.phrqc.selected_output()['calcite'] )
np.save(path + 'De', csh.fluid.Ca.De )
np.save(path + 'poros', csh.fluid.Ca.poros)