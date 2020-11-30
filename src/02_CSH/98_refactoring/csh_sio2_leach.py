# -*- coding: utf-8 -*-

#%% PYTHON MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ch_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
sys.path.append(ch_dir)
import matplotlib.pylab as plt
import numpy as np
#np.set_printoptions(precision=3, threshold=np.inf)
import time
import yantra
import cell_type as ct # change the path to cell_type file
import misc_func as fn
import rt_leach_csh_sio2 as rtl
from copy import deepcopy
#import phrqc
#%% PROBLEM DEFINITION

#problem type
m = 'CSH' #or 'CSH' #TODO case for cement

#%% LOOP
f = 0.01 # fraction
ll = 1 #liquid lauer in front of portlandite
l_ch = 10 #length of portlandite
lx = (l_ch+ll)*1.0e-2
ly = 2.0e-2
dx = 1.0e-2

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, (ll+1): ll+l_ch] = ct.Type.MULTILEVEL
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()

#%%  PHREEQC
nn=os.path.basename(__file__)[:-3] 
nn += str(f)
fn.make_output_dir(root_dir+'\\results\\output\\10_subgrid_leaching')
path = root_dir+'\\results\\output\\10_subgrid_leaching\\' + nn + '\\'
fn.make_output_dir(path)
phase_name = m
phrqc_input = {'csh_mol':{'value': '0.051555' },
               'csh':{'name':phase_name, 
                      'stochiometry':{'Ca':1.3333333333,
                                  'Si':1.0,
                                  'H2O':2.1666666667},
                      'log_k':-10.96765}, #z - x - y
               'ca_bc':{'type':'conc', 'value': '0.0'}, 
               'si_bc':{'type':'conc', 'value': '0.0'},
               'ca_mlvl':{'type':'eq', 'value': phase_name}, 
               'si_mlvl':{'type':'eq', 'value': phase_name}, 
               'ca_liq':{'type':'conc', 'value': '0'},
               'si_liq':{'type':'conc', 'value': '0'}               
               } 
 # another option ca_liq':{'type':'conc', 'value': '0'} or ca_liq':{'type':'eq', 'value': 'portlandite'}

phrqc = rtl.PhreeqcInputCSH(phrqc_input)            
phrqc.save_phrqc_input(root_dir, nn)   

#%% VALUES
scale = 100 # scale of molar volume
mvol = [75.63e-3*scale, 27.3e-3 *scale]#m3/mol # CSH Jennite and SiO2
#mvol = [mvolCH]
max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = [5.1555/scale, 0]
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m) 
D_CSH =    1.e-15
D_border = 1.e-9 # 5.e-10#8*1.e-12
D_high = 1.e-9
D = D_high*(domain.nodetype==-1) + D_CSH*(domain.nodetype!=-1) 
D[1,ll+1] = D_border # default diffusion coefficient in pure liquid
#D = D_high
porosity = (1- pqty[0]*mvol[0])*(domain.nodetype == ct.Type.MULTILEVEL) +\
            1*(domain.nodetype == ct.Type.LIQUID )+\
            1*(domain.nodetype == ct.Type.INTERFACE )+\
            1*(domain.nodetype == ct.Type.SOLID )
app_tort_degree = 1.#1./3.
app_tort = 1. * porosity ** app_tort_degree

settings = {'dissolution':'subgrid', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/'interface'
            'diffusivity':{'type':'fixed', #'archie'/'fixed'
                           'D_border': D_border, #diffusivity at border
                           'D_CSH': D_CSH, # fixed diffusivity in portlandite node
                           },
            'subgrid': {'fraction':f,
                        'poros': False}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': app_tort_degree}, 
            'dx': dx, 
            'Dref':D_high
            }
          
bc_params = {'solution_labels':{'left':100001}, 
            'top':['flux', 0.0],
            'bottom':['flux', 0.0],
            'left':['flux', 0.0],
            'right':['flux', 0.0],}               
#%% PARAMETERS (DOMAIN, BC, SOLVER)

input_file = root_dir + '\\phreeqc_input\\' + nn + '.phrq'
dp={} #domain parameters
dp['D0'] = D                     
dp['voxel_vol']=1
dp['mvol']= mvol

dp['poros'] = porosity                 
dp['app_tort'] = app_tort 

dp['solution_labels']=slabels
dp['database']='cemdata07.dat'

if(len(mvol) == 2):
    dp['phrqc_input_file']=input_file#'CH_CC_Ceq.phrq'
    dp['eq_names'] = ['CSH', 'SiO2am']
    dp['solid_phases']={'CSH':  {'type':'diffusive','mvol':mvol[0],'c':pqty[0]},
                        'SiO2am': {'type':'diffusive','mvol':mvol[1],'c':pqty[1]}}
        
#bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = 1/6., smart_thres = 1e-8, cphi_fact = 1/3.)# optional values, for time step (if tfact => tfactbased tau)
solver_params['phrqc_flags']['smart_run']=False
domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL
fn.save_settings(settings, bc_params, solver_params, path, nn)
#%% INITIATE THE SOLVER
rt= rtl.CSH_Leaching('MultilevelAdvectionDiffusion', domain, 
                     dp, bc_params, solver_params, settings) 
#%% RUN SOLVER
it=time.time()
itr = 0 
j = 0
l = ll
prev_nodetype = deepcopy(rt.nodetype)
rt.dissolution_time = []
rt_csh = [] 
rt_csh.append(np.sum(rt.solid.CSH.c))
rt_time = []
rt_time.append(rt.time)
dport = []
dport.append(0)
Ts = 1000 
dl = 2
while rt.iters <= 100: #rt.time <=Ts: # len(rt.dissolution_time)<dl:#
    rt.advance()             
    if (~np.all(rt.nodetype == prev_nodetype)):
        prev_nodetype = deepcopy(rt.nodetype)
        rt_time.append(rt.time*scale)
        rt.dissolution_time.append(rt.time)
        rt_csh.append(np.sum(rt.solid.CSH.c)*scale)
        if(len(rt_csh)<2):
            dport.append(0)
        else:
            dport.append((rt_csh[-1]-rt_csh[-2])/(rt_time[-1]-rt_time[-2]))
    itr += 1 

#%% SIMULATION TIME
print("Fraction %s done" %str(f))
print("Time to dissolve %s" %str(rt.time*scale))
simulation_time = time.time()-it
fn.print_time(simulation_time, rt)
  
#%%  PRINT
print('Ca %s' %str(np.array(rt.fluid.Ca._c[1,:])))
print('Ca +ss %s' %str(np.array(rt.fluid.Ca.c[1,:]) + np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('H +ss %s' %str(np.array(rt.fluid.H.c[1,:]) + np.array(rt.fluid.H._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('O +ss %s' %str(np.array(rt.fluid.O.c[1,:]) + np.array(rt.fluid.O._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('Si +ss %s' %str(np.array(rt.fluid.Si.c[1,:]) + np.array(rt.fluid.Si._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('SiO2 %s' %str(np.array(rt.solid.SiO2am.c[1,:])))
print('Vol %s' %str(np.array(rt.solid.vol[1,:])))
print('D %s' %str(np.array(rt.fluid.Ca.De[1,:])))
print('pH %s' %str(np.array(rt.phrqc.selected_output()['pH'][1,:])))
print('poros %s' %str(np.array(rt.solid.poros[1,:])))
print('phrqc poros %s' %str(np.array(np.array(rt.phrqc.poros[1,:]))))
print('CSH %s' %str(np.array(rt.solid.CSH.c[1,:])))
print('PHRQC CSH %s' %str(np.array(rt.phrqc.selected_output()['CSH'][1,:])))