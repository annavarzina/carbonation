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
import rt_leach
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

#%%  PHREEQC
nn='03_example_leaching'
fn.make_output_dir(root_dir+'\\results\\output\\00_examples\\')
path = root_dir+'\\results\\output\\00_examples\\' + nn + '\\'
fn.make_output_dir(path)

phrqc_input = {'ca_bc':{'type':'conc', 'value': '0.01'}, 
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
               'ca_liq':{'type':'conc', 'value': '0'}} # another option ca_liq':{'type':'conc', 'value': '0'} or ca_liq':{'type':'eq', 'value': 'portlandite'}

def set_phrqc_input(p, ptype ='CH'):
    def set_phrqc_bc(ca):
        phrqc_input = [] 
        phrqc_input.append('#boundary_solution')    
        phrqc_input.append('SOLUTION\t100001')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']) + '\n')
        else:
            pass
        phrqc_input.append('EQUILIBRIUM_PHASES\t100001\n')
        return phrqc_input
    def set_phrqc_liquid(ca):
        phrqc_input = [] 
        phrqc_input.append('#solution_liquid')    
        phrqc_input.append('SOLUTION\t100002')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
        phrqc_input.append('portlandite\t0\t0')
        return phrqc_input
    
    def set_phrqc_mlvl(ca):
        phrqc_input = [] 
        phrqc_input.append('#solution_multilevel')    
        phrqc_input.append('SOLUTION\t100003')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100003')
        phrqc_input.append('portlandite\t0\t1')
        return phrqc_input
    
    def set_phrqc_solid():
        phrqc_input = [] 
        phrqc_input.append('#solution_solid')    
        phrqc_input.append('SOLUTION\t100005')
        phrqc_input.append('\t-water\t1\n')
        return phrqc_input

    phrqc_input = [] 
    phrqc_input += set_phrqc_bc(p['ca_bc'])
    phrqc_input += set_phrqc_liquid( p['ca_liq'])
    phrqc_input += set_phrqc_mlvl(p['ca_mlvl'])    
    phrqc_input += set_phrqc_solid()
    return phrqc_input

phrqc = set_phrqc_input(phrqc_input)            
fn.save_phrqc_input(phrqc,root_dir, nn)   

#%% VALUES
scale = 50 # scale of molar volume
init_porosCH = 0.05 #initial porosity of portlandite nodes
mvolCH = 0.0331*scale
mvol = [mvolCH]
max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = fn.set_init_pqty(mvol, init_porosCH)
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m)          
D = 1.0e-09 # default diffusion coefficient in pure liquid
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort_degree = 1./3.
app_tort = 1. * porosity ** app_tort_degree

settings = {'dissolution':'subgrid', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/'interface'
            'diffusivity':{'type':'archie', #'archie'/'fixed'
                           'D_border':D, #diffusivity at border
                           'D_CH': 1e-12, # fixed diffusivity in portlandite node
                           },
            'subgrid': {'fraction':None}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': app_tort_degree}, #TODO
            'dx': dx, 
            'Dref':D
            }
               
#%% PARAMETERS (DOMAIN, BC, SOLVER)
domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir + \
                                     '\\phreeqc_input\\' + nn + '.phrq')#'CH_CC-nat.phrq'
bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = None, smart_thres = 1e-8)# optional values, for time step (if tfact => tfactbased tau)
domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL
fn.save_settings(settings, bc_params, solver_params, path, nn)
#%% INITIATE THE SOLVER
rt= rt_leach.LeachingRT('MultilevelAdvectionDiffusion',  domain, 
                          domain_params, bc_params, solver_params,
                          settings) 

#%% PARAMETERS

plist =  [(1,n) for n in np.arange(1, 6)] #points list
pavglist = ['avg_poros', 'pH', 'avg_D_eff', 'sum_vol', #argument list
            'dissolution', 'portlandite_cells'] 
#'delta_ch', 'delta_cc', 'precipitation','dissolution', 'portlandite_cells', 
#'calcite_cells', 'active_cells','dt', 'pH', 'avg_poros',  'avg_D_eff', 'sum_vol'
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)

#%% TIME SETTINGS
nitr =10
Ts =  0.05 #seconds
Ts = Ts/scale + 0.001
step = max(int(Ts/36.),1)
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step))) #time_points = np.arange(0, Ts+step, step)
N = Ts/rt.dt
N_res = 1e+4
S = max(1,int(N/N_res))
it=time.time()

#%% RUN SOLVER
itr = 0 
j = 0
while itr < nitr: # rt.time <=Ts: #
    if(False):
        if ( (rt.time <= time_points[j]) and \
            ((rt.time + rt.dt) > time_points[j]) ):  
            print(time_points[j])
            fn.save_figures_minerals(rt,  max_pqty, time_points[j], path, nn, ptype=m)  
            fn.save_figures_mols(rt, time_points[j], path, nn, ptype=m) 
            #save_vti(rt, phases, time_points[j], path, nn, m)
            #save_pickle(rt, phases, time_points[j], path, nn)
            if(j>0):
                points = [(1,n) for n in np.arange(1,6)]
                fn.print_points(rt, points, names=['portlandite'])
                print('Ca %s' %rt.fluid.Ca.c[1,:])
            j +=1
        
    rt.advance()    
    #results = fn.append_results(rt, results, step = S )
    itr += 1
    
#%% SIMULATION TIME
simulation_time = time.time()-it
fn.print_time(simulation_time, rt)
  
#%%  SAVE
#fn.save_obj(results, path + str(nn) +'_results')
np.save(path + 'pH', rt.phrqc.selected_output()['pH'] )
np.save(path + 'Ca', rt.phrqc.selected_output()['Ca'] )
np.save(path + 'De', rt.fluid.Ca.De )
np.save(path + 'poros', rt.fluid.Ca.poros)

#%% PLOT 
#fn.plot_species(results, names=[])#['calcite']
#fn.plot_avg(results, names=['avg_poros', 'avg_D_eff'])
#fn.plot_points(results, names=['portlandite', 'poros', 'Ca'])
#fn.plot_fields(rt, names=['Ca', 'poros'],fsize=(15,1))

#%% PRINT
print('Ca %s' %str(np.array(rt.fluid.Ca._c[1,:])))
print('Ca +ss %s' %str(np.array(rt.fluid.Ca.c[1,:]) + np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('H +ss %s' %str(np.array(rt.fluid.H.c[1,:]) + np.array(rt.fluid.H._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('O +ss %s' %str(np.array(rt.fluid.O.c[1,:]) + np.array(rt.fluid.O._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('CH %s' %str(np.array(rt.solid.portlandite.c[1,:])))
print('dCH %s' %str(np.array(rt.phrqc.dphases['portlandite'][1,:])))
print('Vol %s' %str(np.array(rt.solid.vol[1,:])))
print('D %s' %str(np.array(rt.fluid.Ca.De[1,:])))
print('pH %s' %str(np.array(rt.phrqc.selected_output()['pH'][1,:])))
print('poros %s' %str(np.array(rt.solid.poros[1,:])))
print('phrqc poros %s' %str(np.array(np.array(rt.phrqc.poros[1,:]))))


#print('Total CH dissolved %s' %(results['portlandite'][-1]-results['portlandite'][0]))