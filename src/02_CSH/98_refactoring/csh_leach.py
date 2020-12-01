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
import rt_leach_csh as rtl
from copy import deepcopy
#import phrqc
#%% PROBLEM DEFINITION

#problem type
m = 'CSH'

f = 1#0.01 # fraction

ll = 1 #liquid lauer in front of csh
l_csh = 40 #40 #length of csh
lx = (l_csh+ll)*1.0e-6
ly = 2.0e-6
dx = 1.0e-6
domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, (ll+1): ll+l_csh] = ct.Type.MULTILEVEL
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()

#%%  PHREEQC
nn=os.path.basename(__file__)[:-3] 
nn += str(f)
fn.make_output_dir(root_dir+'\\results\\output_csh\\03_leaching')
path = root_dir+'\\results\\output_csh\\03_leaching\\' + nn + '\\'
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
mvol = [75.63e-3*scale]#m3/mol # CSH Jennite and SiO2
#mvol = [mvolCH]
max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = [5.1555/scale]
pqty = init_conc[0] * (domain.nodetype == ct.Type.MULTILEVEL)

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

dp['phrqc_input_file']=input_file#'CH_CC_Ceq.phrq'
dp['eq_names'] = ['CSH']
dp['solid_phases']={'CSH':  {'type':'diffusive','mvol':mvol[0],'c':pqty}}
        
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
rt_csh.append(np.sum(rt.solid.CSH.c)*scale)
rt_time = []
rt_time.append(rt.time)
dcsh = [] #dissolution rate
dcsh.append(0)
Ts = 1000 
dl = 10
nitr = 2
while len(rt.dissolution_time)<dl:# rt.iters <= nitr: #rt.time <=Ts: # 
    rt.advance()             
    if (~np.all(rt.nodetype == prev_nodetype)):
        prev_nodetype = deepcopy(rt.nodetype)
        rt_time.append(rt.time*scale)
        rt.dissolution_time.append(rt.time)
        rt_csh.append(np.sum(rt.solid.CSH.c)*scale)
        if(len(rt_csh)<2):
            dcsh.append(0)
        else:
            dcsh.append((rt_csh[-1]-rt_csh[-2])/(rt_time[-1]-rt_time[-2])) #mol/s/um2
    itr += 1 

#%% SIMULATION TIME
print("Fraction %s done" %str(f))
print("Time to dissolve %s" %str(rt.time*scale))
simulation_time = time.time()-it
fn.print_time(simulation_time, rt)

#%%  SAVE
np.save(path + 'dis_time', rt.dissolution_time )
np.save(path + 'dCH', dcsh)
np.save(path + 'time', rt_time)
np.save(path + 'CH', rt_csh)
np.save(path + 'Ca_prof', rt.fluid.Ca.c[1,:]+ np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:]))
np.save(path + 'Si_prof', rt.fluid.Si.c[1,:]+ np.array(rt.fluid.Si._ss[1,:])/np.array(rt.phrqc.poros[1,:]))
    
  
#%%  PRINT
print('Ca %s' %str(np.array(rt.fluid.Ca._c[1,:])))
print('Ca +ss %s' %str(np.array(rt.fluid.Ca.c[1,:]) + np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('H +ss %s' %str(np.array(rt.fluid.H.c[1,:]) + np.array(rt.fluid.H._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('O +ss %s' %str(np.array(rt.fluid.O.c[1,:]) + np.array(rt.fluid.O._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('Si +ss %s' %str(np.array(rt.fluid.Si.c[1,:]) + np.array(rt.fluid.Si._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
print('Vol %s' %str(np.array(rt.solid.vol[1,:])))
print('D %s' %str(np.array(rt.fluid.Ca.De[1,:])))
print('pH %s' %str(np.array(rt.phrqc.selected_output()['pH'][1,:])))
print('poros %s' %str(np.array(rt.solid.poros[1,:])))
print('phrqc poros %s' %str(np.array(np.array(rt.phrqc.poros[1,:]))))
print('CSH %s' %str(np.array(rt.solid.CSH.c[1,:])))
print('PHRQC CSH %s' %str(np.array(rt.phrqc.selected_output()['CSH'][1,:])))

#%% Rate
#dim = 10**3*10**(-15)*10**8 # convert to mmol/l/cm2/s
dim =  10**(-15)*10**12 #convert to mol/l/m2/s
r1 = dcsh[1]*dim #mol/s
print('R1: %s' %str(r1))
n= 10 #-11
rn = dcsh[n]*dim
print('Rn: %s' %str(rn))
#%%
'''
r1 = []
rn = []
rmean = []
for f in fractions:    
    r1.append(dch[nn + str(f)][1]*scale*dim)
    rn.append(dch[nn + str(f)][-11]*scale*dim)
    rmean.append(np.mean(dch[nn + str(f)][:])*scale*dim)
    
    
plt.figure()
plt.plot(fractions, np.abs(r1), label = "1")
plt.plot(fractions, np.abs(rn), label = "20")
plt.plot(fractions, np.abs(rmean), label = "mean")
plt.ylabel(r"Rate of dissolution $k$ $(mmol / (l\cdot s\cdot cm^2)$")
plt.xlabel(r"Fraction of equilibrated water $\sigma$")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.show()


r = rn#rmean
idx = np.where(np.logical_and(np.abs(r)>=0.39e-5, np.abs(r)<=6.2e-5))[0]
plt.figure(figsize = (6,4))
plt.plot(np.abs(r), fractions)
plt.fill_between(np.abs(r)[idx],fractions[idx], color = "#6f8191", alpha=.5,label = "Johannsen K., et.al. (1999)")
plt.xlabel(r"Rate of dissolution $R$ $(mmol / (l\cdot s\cdot cm^2)$",fontsize=14)
plt.ylabel(r"Mixing parameter $\sigma$ ", fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.yscale("log")
plt.xscale("log")
plt.legend(fontsize=12)
plt.show()
'''