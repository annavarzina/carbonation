
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
ll = 1 #liquid lauer in front of portlandite
l_csh = 1 #length of CSH
lx = (l_csh+ll+1)*1.0e-6
ly = 2.0e-6
dx = 1.0e-6

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, ll+1: ll+l_csh+1] = ct.Type.MULTILEVEL
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID



#%%  VALUES
m = 'CSH' #or 'CSH'
nn=os.path.basename(__file__)[:-3]
fn.make_output_dir(root_dir+'\\results\\output_csh\\01_default\\')
path = root_dir+'\\results\\output_csh\\01_default\\' + nn + '\\'
fn.make_output_dir(path)


phrqc_input = {'c_bc':{'type':'pco2', 'value': 3.41},
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
#init_poros = [0.05, 1.0, 0.99424327, 0.87988525, 0.83697953, 0.73712387, 1.0]
D = 1.0e-09 # default diffusion coefficient in pure liquid
app_tort_degree = 7.43#1./3.
settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            'dissolution':'multilevel', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/
            'diffusivity':{'border': D, ##diffusivity at border
                           'CH': ('const', 1e-15), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CC': ('inverse', 1e-13), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CSH': ('inverse', 1e-12), # fixed diffusivity in CSH node 'archie'/'const'/'inverse'
                           }, 
            'pcs_mode': {'pcs': True, #Pore-Size Controlled Solubility concept
                         'pores': 'block', #'block'/'cylinder'
                         'int_energy': 0.5, # internal energy
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
res = rt.ResultsCSHQ(nodes= [(1,n) for n in np.arange(0, 4)])
#%% results dict
nitr = 5
Ts  = 4.0*3600# 36 * 3 #s
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
while csh.time <=Ts: # itr < nitr: # 
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
#res.plot_species(names=[])
#res.plot_avg()
#res.plot_nodes(names=[])
#res.plot_fields(csh, names=[],fsize=(15,1))
#%% plot Ca/Si against density
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(res.results['ratio_CaSi (1, 2)'], res.results['density_CSHQ'])
#plt.legend(fontsize=18)
plt.ylabel('C-S-H density (g/l)', fontsize=16)
plt.xlabel('Ca/Si', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_equi_casi_density.png', bbox_inches='tight')
#plt.show()
#%% plot Ca/Si against density
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['ratio_CaSi (1, 2)'])
#plt.legend(fontsize=18)
plt.ylabel('Ca/Si', fontsize=16)
plt.xlabel('Time (h)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_equi_casi.png', bbox_inches='tight')
#plt.show()
#%% CSHQ phases
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['CSHQ_JenD'], label = 'JenD', ls = '-')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['CSHQ_JenH'], label = 'JenH', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['CSHQ_TobD'], label = 'TobD', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['CSHQ_TobH'], label = 'TobH', ls = ':')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['calcite'], label = '$C\overline{C}$', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['sio2am'], label = '$SiO_2$', ls = '--')
plt.xlabel('Time (h)', fontsize=16)
plt.ylabel(r'C-S-H ($mol/dm^3$)', fontsize=16)
plt.legend(fontsize=10, loc = 'upper left' )
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_equi_csh_phases_5h.png', bbox_inches='tight')
#%% CSHQ phases 30 min
t = 1000
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])[0:t]*scale/60, res.results['CSHQ_JenD'][0:t], label = 'JenD', ls = '-')
plt.plot(np.array(res.results['time'])[0:t]*scale/60, res.results['CSHQ_JenH'][0:t], label = 'JenH', ls = '--')
plt.plot(np.array(res.results['time'])[0:t]*scale/60, res.results['CSHQ_TobD'][0:t], label = 'TobD', ls = '-.')
plt.plot(np.array(res.results['time'])[0:t]*scale/60, res.results['CSHQ_TobH'][0:t], label = 'TobH', ls = ':')
plt.plot(np.array(res.results['time'])[0:t]*scale/60, res.results['calcite'][0:t], label = r'$C\overline{C}$', ls = '--')
plt.plot(np.array(res.results['time'])[0:t]*scale/60, res.results['sio2am'][0:t], label = r'$SiO_2$', ls = '--')
plt.xlabel('Time (min)', fontsize=16)
plt.ylabel(r'C-S-H ($mol/dm^3$)', fontsize=16)
plt.legend(fontsize=10)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_equi_csh_phases_30min.png', bbox_inches='tight')
#%% pH 
t = 1000
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['pH (1, 1)'], label = 'Water voxel', ls = '-')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['pH (1, 2)'], label = 'C-S-H voxel', ls = '--')
plt.xlabel('Time (h)', fontsize=16)
plt.ylabel('pH', fontsize=16)
plt.legend(fontsize=10)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_equi_ph_change.png', bbox_inches='tight')
#%% CC and SiO2
t = 1000
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['calcite (1, 1)'], label = r'$C\overline{C}$ in water voxel', ls = '-', color = 'tab:blue')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['calcite (1, 2)'], label = r'$C\overline{C}$ C-S-H voxel', ls = '-', color = 'tab:orange')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['sio2am (1, 1)'], label = r'$SiO_2$ in water voxel', ls = '--', color = 'tab:blue')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['sio2am (1, 2)'], label = r'$SiO_2$ C-S-H voxel', ls = '--', color = 'tab:orange')
plt.xlabel('Time (h)', fontsize=16)
plt.ylabel(r'Phase ($mol/dm^3$)', fontsize=16)
plt.legend(fontsize=10, loc = 'upper left')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_equi_cc_sio2_change.png', bbox_inches='tight')
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

#%% SAVE
'''
rt.ResultsCSHQ.save_obj(res.results, path + str(nn) +'_results')
np.save(path + 'SI', csh.phrqc.selected_output()['SI_calcite'][1,:])
np.save(path + 'pH', csh.phrqc.selected_output()['pH'][1,:])
np.save(path + 'Ca', np.array(csh.fluid.Ca.c[1,:]) + \
               np.array(csh.fluid.Ca._ss[1,:])/np.array(csh.phrqc.poros[1,:]))
np.save(path + 'Si', np.array(csh.fluid.Si.c[1,:]) + \
               np.array(csh.fluid.Si._ss[1,:])/np.array(csh.phrqc.poros[1,:]))
np.save(path + 'C', np.array(csh.fluid.C.c[1,:]) + \
               np.array(csh.fluid.C._ss[1,:])/np.array(csh.phrqc.poros[1,:]))
np.save(path + 'CC', csh.solid.calcite.c[1,:] )
np.save(path + 'De', csh.fluid.Ca.De[1,:])
np.save(path + 'poros', csh.fluid.Ca.poros[1,:])

np.save(path + 'CSHQ_TobD', np.array(csh.solid.CSHQ_TobD.c[1,:]))
np.save(path + 'CSHQ_TobH', np.array(csh.solid.CSHQ_TobH.c[1,:]))
np.save(path + 'CSHQ_JenD', np.array(csh.solid.CSHQ_JenD.c[1,:]))
np.save(path + 'CSHQ_JenH', np.array(csh.solid.CSHQ_JenH.c[1,:]))
np.save(path + 'SiO2_am', np.array(csh.solid.sio2am.c[1,:]))
'''