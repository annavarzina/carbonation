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
np.set_printoptions(precision=3, threshold=np.inf)
import time
import yantra
import cell_type as ct # change the path to cell_type file
import misc_func as fn
import rt1
from copy import deepcopy
#import phrqc
#%% PROBLEM DEFINITION
__doc__= """ 
Default example. Explains possible values
1D carbonation of portlandite/cement.
"""
#problem type
m = 'CH' #or 'CSH' #TODO case for cement

#%% GEOMETRY
ll = 2 #liquid lauer in front of portlandite
l_ch = 20 #length of portlandite
lx = (l_ch+ll)*1.0e-3
ly = 2.0e-3
dx = 1.0e-3

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
fn.make_output_dir(root_dir+'\\results\\temp')
path = root_dir+'\\results\\temp\\' + nn + '\\'
fn.make_output_dir(path)

phrqc_input = {'ca_bc':{'type':'conc', 'value': '0.0'}, 
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
               'ca_liq':{'type':'conc', 'value': '0'}} # {'type':'eq', 'value': 'portlandite'}}


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
        phrqc_input.append('portlandite\t0\t0')
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
scale = 100 # scale of molar volume
init_porosCH = 0.05 #initial porosity of portlandite nodes
mvolCH =0.0331*scale
mvol = [mvolCH]
max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = fn.set_init_pqty(mvol, init_porosCH)
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m) 
D_CH =    1.e-15
D_border = 1.e-9 # 5.e-10#8*1.e-12
D_high = 1.e-9
D = D_high*(domain.nodetype==-1) + D_CH*(domain.nodetype!=-1) 
D[1,ll+1] = D_border # default diffusion coefficient in pure liquid
#D = D_high
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort_degree = 1.#1./3.
app_tort = 1. * porosity ** app_tort_degree

settings = {'dissolution':'subgrid', #'multilevel'/'subgrid'
            'active_nodes': 'smart', # 'all'/'smart'/'interface'
            'diffusivity':{'type':'fixed', #'archie'/'fixed'
                           'D_border': D_border, #diffusivity at border
                           'D_CH': D_CH, # fixed diffusivity in portlandite node
                           },
            'subgrid': {'fraction':1.,
                        'poros': False}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': app_tort_degree}, #TODO
            'dx': dx, 
            'Dref':D_high
            }
          
bc_params = {#'solution_labels':{'left':100001}, 
            'top':['flux', 0.0],
            'bottom':['flux', 0.0],
            'left':['flux', 0.0],
            'right':['flux', 0.0],}               
#%% PARAMETERS (DOMAIN, BC, SOLVER)
domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir + \
                                     '\\phreeqc_input\\' + nn + '.phrq')#'CH_CC-nat.phrq'
#bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = None, smart_thres = 1e-8)# optional values, for time step (if tfact => tfactbased tau)
solver_params['phrqc_flags']['smart_run']=False
domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL
fn.save_settings(settings, bc_params, solver_params, path, nn)
#%% INITIATE THE SOLVER
rt= rt1.LeachingRT('MultilevelAdvectionDiffusion',  domain, 
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
nitr =100#2000
Ts =  500*3600. #seconds
Ts = Ts/scale + 0.001
N = Ts/rt.dt
N_res = 1e+4
S = max(1,int(N/N_res))
it=time.time()


step = Ts/20.
time_points = np.arange(0, Ts+step, step)

#%% RUN SOLVER
itr = 0 
j = 0
l = ll
prev_nodetype = deepcopy(rt.nodetype)
rt.dissolution_time = []
rt_port = []
rt_time = []
dport = []
avg_Ca = []
tot_Ca = []
while  rt.time <=Ts: # itr < nitr: #
    rt.advance() 
    if (rt.iters%S == 0):
        avg_Ca.append(np.sum(rt.fluid.Ca.c[1,0:ll])/ll)
        tot_Ca.append(np.sum(rt.fluid.Ca.c[1,0:ll]*ll)*40.078*10**(-15))
        rt_port.append(np.sum(rt.solid.portlandite.c))
        if(len(rt_port)<2):
            dport.append(0)
        else:
            dport.append((rt_port[-1]-rt_port[-2])/rt.dt)
        rt_time.append(rt.time)
    if (~np.all(rt.nodetype == prev_nodetype)):
        prev_nodetype = deepcopy(rt.nodetype)
        rt.dissolution_time.append(rt.time)
    itr += 1

rt_time = np.array(rt_time)*scale    
#%% SIMULATION TIME
simulation_time = time.time()-it
fn.print_time(simulation_time, rt)
  
#%%  SAVE
'''
np.save(path + 'dis_time', rt.dissolution_time )
np.save(path + 'dCH', dport)
np.save(path + 'time', rt_time)
np.save(path + 'CH', rt_port)
np.save(path + 'Ca_prof', rt.fluid.Ca.c[1,:]+ np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:]))
'''
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
print('tau %s' %str(np.array(rt.fluid.Ca.tau[1,:])))


#print('Total CH dissolved %s' %(results['portlandite'][-1]-results['portlandite'][0]))
#%%
print('time %s seconds' %str(rt.time))
plt.figure()
plt.plot(np.array(rt_time)/3600, rt_port)
plt.xlabel("Time")
plt.ylabel("Portlandite")
plt.show()


plt.figure()
plt.plot(np.array(rt_time)/3600, np.abs(dport))
plt.xlabel("Time")
plt.ylabel("Portlandite rate")
plt.show()

plt.figure()
plt.plot(np.array(rt_time)/3600, avg_Ca)
#plt.xscale("log")
plt.xlabel("Time")
plt.ylabel("Concentration Ca (average)")
plt.show()


plt.figure()
plt.plot(np.array(rt_time)/3600, tot_Ca)
plt.ylabel("Ca total mass")
plt.xlabel("Time")
#plt.xscale("log")
plt.show()


plt.figure()
plt.plot( rt.fluid.Ca.c[1,:]+ np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.solid.poros[1,:]))
plt.ylabel("Concentration")
plt.xlabel("X")
plt.show()

#%%
'''
#print(rt.dissolution_time)
t = np.array(rt.dissolution_time)*scale
x = np.arange(1,len(rt.dissolution_time)+1)
plt.figure()
plt.plot(t, x)
plt.ylabel("Dissolved length (um)")
plt.xlabel("Time")
plt.show()

x = x*1e-6
plt.figure()
D = x**2/2/t
plt.plot(t, D)
plt.ylabel("D (m2/s)")
plt.xlabel("Time (s)")
plt.show()
'''

#%%
x = np.array([0, 200.])
y = np.array([0, 0.0175])
amCa = 40.078 #g/mol
def mol2ppm(c, am):
    p = c * am * 1000
    return p #ppm
yppm = mol2ppm(y,amCa)
d = yppm[1]/x[1] # ppm/second/um2
print(d)
#y = 0.015/20*x
