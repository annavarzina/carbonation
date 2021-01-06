
# -*- coding: utf-8 -*-
'''
CO2 3%
PCO2 1.52
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
l_csh = 25 #length of CSH
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


phrqc_input = {'c_bc':{'type':'pco2', 'value': 1.52},
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
            'dissolution':'subgrid', #'multilevel'/'subgrid'
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
res = rt.ResultsCSHQ(nodes= [(1,n) for n in np.arange(3, 9)])
#%% results dict
nitr = 5
Ts  = 10*3600# 36 * 3 #s
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
#res.plot_species(names=[])#['calcite']
#res.plot_avg()
#res.plot_nodes(names=[])
#res.plot_fields(csh, names=[],fsize=(15,1))
#%% plot Ca/Si against density
'''
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(res.results['ratio_CaSi (1, 5)'], res.results['density_CSHQ'])
plt.ylabel('C-S-H density (g/l)', fontsize=16)
plt.xlabel('Ca/Si', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_sg_casi_density.png', bbox_inches='tight')
#plt.show()

#%% plot Ca/Si 
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['ratio_CaSi (1, 5)'], label = 'Voxel 6', ls = '-')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['ratio_CaSi (1, 6)'], label = 'Voxel 7', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['ratio_CaSi (1, 7)'], label = 'Voxel 8', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['ratio_CaSi (1, 8)'], label = 'Voxel 9', ls = ':')
plt.legend(fontsize=10)
plt.ylabel('Ca/Si', fontsize=16)
plt.xlabel('Time (h)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_sg_casi_time.png', bbox_inches='tight')
#plt.show()
#%% plot Ca 
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Ca (1, 3)'], label = 'Voxel 4', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Ca (1, 4)'], label = 'Voxel 5', ls = ':')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Ca (1, 5)'], label = 'Voxel 6', ls = '-')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Ca (1, 6)'], label = 'Voxel 7', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Ca (1, 7)'], label = 'Voxel 8', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Ca (1, 8)'], label = 'Voxel 9', ls = ':')
plt.legend(fontsize=10)
plt.ylabel('Ca (mol/l)', fontsize=16)
plt.xlabel('Time (h)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_sg_Ca_time.png', bbox_inches='tight')
#%% plot C
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['C (1, 3)'], label = 'Voxel 4', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['C (1, 4)'], label = 'Voxel 5', ls = ':')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['C (1, 5)'], label = 'Voxel 6', ls = '-')
plt.legend(fontsize=10)
plt.ylabel('C (mol/l)', fontsize=16)
plt.xlabel('Time (h)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_sg_C_time.png', bbox_inches='tight')
#%% plot Si
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Si (1, 3)'], label = 'Voxel 4', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Si (1, 4)'], label = 'Voxel 5', ls = ':')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Si (1, 5)'], label = 'Voxel 6', ls = '-')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Si (1, 6)'], label = 'Voxel 7', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Si (1, 7)'], label = 'Voxel 8', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['Si (1, 8)'], label = 'Voxel 9', ls = ':')
plt.legend(fontsize=10)
plt.ylabel('Si (mol/l)', fontsize=16)
plt.xlabel('Time (h)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_sg_Si_time.png', bbox_inches='tight')
#%% plot pH
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['pH (1, 3)'], label = 'Voxel 4', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['pH (1, 4)'], label = 'Voxel 5', ls = ':')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['pH (1, 5)'], label = 'Voxel 6', ls = '-')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['pH (1, 6)'], label = 'Voxel 7', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['pH (1, 7)'], label = 'Voxel 8', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['pH (1, 8)'], label = 'Voxel 9', ls = ':')
plt.legend(fontsize=10)
plt.ylabel('pH', fontsize=16)
plt.xlabel('Time (h)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_sg_ph_time.png', bbox_inches='tight')
#%% plot vol_CSHQ
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['vol_CSHQ (1, 5)'], label = 'Voxel 6', ls = '-')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['vol_CSHQ (1, 6)'], label = 'Voxel 7', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['vol_CSHQ (1, 7)'], label = 'Voxel 8', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['vol_CSHQ (1, 8)'], label = 'Voxel 9', ls = ':')
plt.legend(fontsize=10)
plt.ylabel(r'C-S-H volume ($\mu m^3$)', fontsize=16)
plt.xlabel('Time (h)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_sg_vol_CSHQ_time.png', bbox_inches='tight')
#%% CSHQ phases
fig = plt.figure(figsize = (5,3), dpi = 200)
plt.plot(np.array(res.results['time'])*scale/3600, res.results['CSHQ_JenD'], label = 'JenD', ls = '-')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['CSHQ_JenH'], label = 'JenH', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['CSHQ_TobD'], label = 'TobD', ls = '-.')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['CSHQ_TobH'], label = 'TobH', ls = ':')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['calcite'], label = '$C\overline{C}$', ls = '--')
plt.plot(np.array(res.results['time'])*scale/3600, res.results['sio2am'], label = '$SiO_2$', ls = '--')
plt.xlabel('Time (h)', fontsize=16)
plt.ylabel(r'C-S-H $\cdot 10^{-15}$ $(mol)$', fontsize=16)
plt.legend(fontsize=10, loc = 'upper left' )
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig.savefig('07_sg_csh_phases_3h.png', bbox_inches='tight')
#%% Ca profile
fig = plt.figure(figsize = (5,3), dpi = 200) 
plt.imshow(csh.fluid.Ca.c[1:-1,0:-1])     
plt.xlabel(r'Distance ($\mu m$)', fontsize = 12)
plt.title('Ca concentration')
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
clb = plt.colorbar(orientation="horizontal", pad = 0.4)
clb.ax.set_title(r'Colorbar: Ca concentration ($mol/l$)', fontsize = 12) 
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)  
plt.show()
fig.savefig('07_sg_Ca_profile.png', bbox_inches='tight')
#%% C profile
fig = plt.figure(figsize = (5,3), dpi = 200) 
plt.imshow(csh.fluid.C.c[1:-1,0:-1])     
plt.xlabel(r'Distance ($\mu m$)', fontsize = 12)
plt.title('C concentration')
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
clb = plt.colorbar(orientation="horizontal", pad = 0.4)
clb.ax.set_title(r'Colorbar: C concentration ($mol/l$)', fontsize = 12) 
clb.ax.locator_params(nbins=5)
#clb.ax.set_xticklabels(clb.ax.get_xticklabels(), rotation=-45)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)  
plt.show()
fig.savefig('07_sg_C_profile.png', bbox_inches='tight')
#%% Si profile
fig = plt.figure(figsize = (5,3), dpi = 200) 
plt.imshow(csh.fluid.Si.c[1:-1,0:-1])     
plt.title('Si concentration')
plt.xlabel(r'Distance ($\mu m$)', fontsize = 12)
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
clb = plt.colorbar(orientation="horizontal", pad = 0.4)
clb.ax.set_title(r'Colorbar: Si concentration ($mol/l$)', fontsize = 12) 
clb.ax.locator_params(nbins=5)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)  
plt.show()
fig.savefig('07_sg_Si_profile.png', bbox_inches='tight')
#%% CC profile
fig = plt.figure(figsize = (5,3), dpi = 200) 
plt.imshow(csh.solid.calcite.c[1:-1,0:-1])     
plt.xlabel(r'Distance ($\mu m$)', fontsize = 12)
plt.title('Calcite ')
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
clb = plt.colorbar(orientation="horizontal", pad = 0.4)
clb.ax.set_title(r'Colorbar: Calcite amount ($mol/dm^3$)', fontsize = 12) 
#clb.ax.locator_params(nbins=5)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)  
plt.show()
fig.savefig('07_sg_CC_profile.png', bbox_inches='tight')
#%% CSH profile
fig = plt.figure(figsize = (5,3), dpi = 200) 
plt.imshow(csh.solid.csh_vol[1:-1,0:-1])     
plt.xlabel(r'Distance ($\mu m$)', fontsize = 12)
plt.title('C-S-H volume')
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
clb = plt.colorbar(orientation="horizontal", pad = 0.4)
clb.ax.set_title(r'Colorbar: C-S-H volume ($\mu m^3$)', fontsize = 12) 
#clb.ax.locator_params(nbins=5)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)  
plt.show()
fig.savefig('07_sg_cshvol_profile.png', bbox_inches='tight')
#%% Cporos profile
fig = plt.figure(figsize = (5,3), dpi = 200) 
plt.imshow(csh.solid.poros[1:-1,0:-1], vmin = 0, vmax = 1)     
plt.xlabel(r'Distance ($\mu m$)', fontsize = 12)
plt.title(r'Free volume ($\mu m^3$)')
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
clb = plt.colorbar(orientation="horizontal", pad = 0.4)
clb.ax.set_title(r'Colorbar: Free volume ($\mu m^3$)', fontsize = 12) 
#clb.ax.locator_params(nbins=5)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)  
plt.show()
fig.savefig('07_sg_poros_profile.png', bbox_inches='tight')
#%% ph
fig = plt.figure(figsize = (5,3), dpi = 200) 
plt.imshow(csh.phrqc.selected_output()['pH'][1:-1,0:-1])     
plt.xlabel(r'Distance ($\mu m$)', fontsize = 12)
plt.title('pH')
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
clb = plt.colorbar(orientation="horizontal", pad = 0.4)
clb.ax.set_title(r'Colorbar: pH', fontsize = 12) 
#clb.ax.locator_params(nbins=5)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)  
plt.show()
fig.savefig('07_sg_ph_profile.png', bbox_inches='tight')
'''
#%% SAVE

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
