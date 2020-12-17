
# -*- coding: utf-8 -*-
'''
Molar volume = 100
Liquid layer = 5 um
Fraction = 0.004
'''

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

#%% MODULES
import time
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import yantra
import cell_type as ct # change the path to cell_type file
import rt_carb_csh_sio2 as rt #rt_r as rt
import misc_func as fn
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



#%%  VALUES
m = 'CSH' #or 'CSH'
nn=os.path.basename(__file__)[:-3]
fn.make_output_dir(root_dir+'\\results\\output_csh\\01_default\\')
path = root_dir+'\\results\\output_csh\\01_default\\' + nn + '\\'
fn.make_output_dir(path)


phrqc_input = {'c_bc':{'type':'pco2', 'value': 1.0},
               #'c_bc':{'type':'conc', 'value': 0.01},
               'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
               'c_liq':{'type':'eq', 'value': 'calcite'},
               'ca_mlvl':{'type':'eq', 'value': 'calcite'},    
               'ca_liq':{'type':'eq', 'value': 'calcite'}}

phrqc = rt.PhreeqcInputCSHQ(phrqc_input)            
phrqc.save_phrqc_input(root_dir, nn)   


tfact_default = 1./6./1#*init_porosCH
scale = 50 # scale of molar volume
init_poros = [0.05, 1.0, 0.99424327, 0.87988525, 0.83697953, 0.73712387, 1.0] #initial porosity of portlandite or calcite nodes
D = 1.0e-09 # default diffusion coefficient in pure liquid
app_tort_degree = 1./3.
settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            'dissolution':'subgrid', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/
            'diffusivity':{'border': D, ##diffusivity at border
                           'CH': ('const', 1e-15), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CC': ('inverse', 1e-11), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CSH': ('const', 1e-11), # fixed diffusivity in CSH node 'archie'/'const'/'inverse'
                           }, 
            'pcs_mode': {'pcs': True, #Pore-Size Controlled Solubility concept
                         'pores': 'block', #'block'/'cylinder'
                         'int_energy': 0.1, # internal energy
                         'pore_size': 0.01*dx, # threshold radius or distance/2
                         'crystal_size': 0.5*dx, # crystal or pore length
                         }, 
            'subgrid': {'fraction':1.0e-6}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': app_tort_degree}, 
            'velocity': False, 
            'bc': phrqc_input['c_bc'],
            'dx': dx, 
            'Dref':D}
            

mvol = rt.SettingsCSHQ.set_mvols(scale) #dm3/mol
max_pqty = rt.SettingsCSHQ.get_max_pqty(mvol) #mol/dm3
init_conc = rt.SettingsCSHQ.set_init_pqty(mvol, init_poros, scale)
pqty = rt.SettingsCSHQ.get_pqty(init_conc, domain)
slabels = rt.SettingsCSHQ.set_labels(domain)          
porosity = rt.SettingsCSHQ.get_porosity(domain, pqty, mvol)
app_tort =  rt.SettingsCSHQ.get_app_tort(domain, app_tort_degree, porosity)
            
#%% PARAMETERS (DOMAIN, BC, SOLVER)
domain_params = rt.SettingsCSHQ.set_domain_params(D, mvol, pqty, porosity, 
                    app_tort, slabels, database = 'cemdata18.dat',
                    input_file = root_dir +'\\phreeqc_input\\' + nn + '.phrq')
bc_params = rt.SettingsCSHQ.set_bc_params(bc_slabels = {'left':100001})
solver_params = rt.SettingsCSHQ.set_solver_params(tfact = tfact_default, smart_thres = 1e-8, cphi_fact = 1/3.)
#domain = rt.SettingsCSHQ.reset_portlandite_nodes(domain)
#rt.SettingsCSHQ.save_settings(settings, bc_params, solver_params, path, nn)

csh= rt.CarbonationCSHQ('MultilevelDiffusion',  domain, 
                          domain_params, bc_params, solver_params, settings) 
res = rt.ResultsCSHQ(nodes= [(1,n) for n in np.arange(3, 6)])
#%% results dict
nitr = 1000
Ts  = 30# 36 * 3 #s
Ts = Ts/scale + 0.001
N = Ts/csh.dt
N_res = 1e+4
S = max(1,int(N/N_res))
step = max(int(Ts/10.),1)
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step)))
it=time.time()


#%% run
itr = 0 
j = 0
while itr < nitr: # csh.time <=Ts: # 
    if(True):
        if ( (csh.time <= time_points[j]) and \
            ((csh.time + csh.dt) > time_points[j]) ): 
            print(time_points[j])
            j +=1        
    csh.advance()    
    res.append_results(csh, step = S)
    itr += 1
       
simulation_time = time.time()-it
fn.print_time(simulation_time, csh)


#%% plot and print

rt.ResultsCSHQ.print_profiles(csh)
#res.plot_species(names=[])#['calcite']
#res.plot_avg()
#res.plot_nodes(names=[])
#res.plot_fields(csh, names=[],fsize=(15,1))
#%% plot ca/si against density
'''
plt.figure()
plt.plot(results['Ca_Si'], results['csh_density'])
plt.legend()
plt.ylabel('CSH density')
plt.xlabel('C/S')
plt.show()
'''

#%% C-S-H profile
'''
linetype = np.array(['dotted', 'solid', 'dashed', 'dashdot', 'dotted'])
x = range(1, 9)

fig, ax1 = plt.subplots(figsize=(8,4))
ax1.set_xlabel(r'Distance ($\mu$m)',fontsize = 14)
ax1.set_ylabel('C-S-H Q phases',fontsize = 14)
ax1.plot(x, csh.phrqc.selected_output()['CSHQ_TobD'][1,1:-2],  label='CSHQ_TobD', ls=linetype[1])
ax1.plot(x, csh.phrqc.selected_output()['CSHQ_TobH'][1,1:-2],  label='CSHQ_TobH', ls=linetype[1])
ax1.plot(x, csh.phrqc.selected_output()['CSHQ_JenD'][1,1:-2],  label='CSHQ_JenD', ls=linetype[1])
ax1.plot(x, csh.phrqc.selected_output()['CSHQ_JenH'][1,1:-2],  label='CSHQ_JenH', ls=linetype[1])
plt.legend(loc = "center left")
plt.show()
'''