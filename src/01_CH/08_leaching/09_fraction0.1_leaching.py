# -*- coding: utf-8 -*-
'''
Sript for both CH and CSH systems
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
np.set_printoptions(precision=3, threshold=np.inf)
import time
import yantra
import cell_type as ct # change the path to cell_type file
import misc_func as fn
import rt_leach as rt1
#import phrqc
#%% PROBLEM DEFINITION
__doc__= """ 
Default example. Explains possible values
1D carbonation of portlandite/cement.
"""
#problem type
m = 'CH' #or 'CSH' #TODO case for cement

#%% GEOMETRY
ll = 100 #liquid lauer in front of portlandite
l_ch = 3 #length of portlandite
lx = (l_ch+ll)*1.0e-6
ly = 2.0e-6
dx = 1.0e-6

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
fn.make_output_dir(root_dir+'\\results\\output\\00_examples\\')
path = root_dir+'\\results\\output\\00_examples\\' + nn + '\\'
fn.make_output_dir(path)

phrqc_input = {'ca_bc':{'type':'conc', 'value': '0.0'}, 
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
scale = 50 # scale of molar volume
init_porosCH = 0.1 #initial porosity of portlandite nodes
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
            'subgrid': {'fraction':0.1,
                        'poros': False}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': app_tort_degree}, #TODO
            'dx': dx, 
            'Dref':D_high
            }
          
bc_params = {#'solution_labels':{'left':100001}, 
            'top':['flux', 0.0],
            'bottom':['flux', 0.0],
            'left':['open'],#['flux', 0.0],
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
nitr =1000#2000
Ts =  0.001*scale #seconds
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
rt_port = []
rt_time = []
dport = []
conc_step = []
conc_step1 = []
conc_step2 = []
while  itr < nitr: #rt.time <=Ts: # 
    '''       
    if ( (rt.time <= time_points[j]) and ((rt.time + rt.dt) > time_points[j]) ):  
        conc_step.append(rt.fluid.Ca.c[1,l])            
        j +=1
    '''
    rt.advance() 
    conc_step.append(rt.fluid.Ca.c[1,l+1]+ np.array(rt.fluid.Ca._ss[1,l+1])/np.array(rt.phrqc.poros[1,l+1]))  
    conc_step1.append(rt.fluid.Ca.c[1,l-0]+ np.array(rt.fluid.Ca._ss[1,l])/np.array(rt.phrqc.poros[1,l]))
    conc_step2.append(rt.fluid.Ca.c[1,l-1])  
    rt_port.append(np.sum(rt.solid.portlandite.c))
    if(rt.iters<2):
        dport.append(0)
    else:
        dport.append((rt_port[-1]-rt_port[-2])/rt.dt)
    rt_time.append(rt.time)
    itr += 1
    
#%% SIMULATION TIME
simulation_time = time.time()-it
fn.print_time(simulation_time, rt)
  
#%%  SAVE
#fn.save_obj(results, path + str(nn) +'_results')

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
plt.plot(rt_time, rt_port)
plt.show()


plt.figure()
plt.plot(rt_time[2:-1], dport[2:-1])
plt.show()

plt.figure()
plt.plot(rt_time, conc_step)
plt.xscale("log")
plt.show()


plt.figure()
plt.plot(rt_time, conc_step1)
plt.xscale("log")
plt.show()

plt.figure()
plt.plot(rt_time, conc_step2)
plt.xscale("log")
plt.show()

plt.figure()
plt.plot( rt.fluid.Ca.c[1,:]+ np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:]))
plt.show()

#%% DIFFUSION COEFFICIENT
D = 1e-9
c = rt.fluid.Ca.c[1,l] #mol/l
print('Concentration %s' %c)
cm = 0.0194923#rt.fluid.Ca.c[1,-2] #mol/l
t = rt.time #s

from math import erf
from scipy.optimize import root
def equation(d, x, t, cm, c):
    return (c - cm*(1.-erf(x/2./np.sqrt(d*t))))

x = 1 #m
diffusivity0 = []
for i in np.arange(1, len(conc_step)):
    res = root(equation, D/dx**2, args = (x, rt_time[i], cm, conc_step[i]), method='lm')
    #print(res.x[0])
    diffusivity0.append(res.x[0]*dx**2)
#print(diffusivity1)
print('Diffusivity %s' %np.mean(diffusivity0))

x = x+1.
diffusivity1 = []
for i in np.arange(1, len(conc_step)):
    res = root(equation,D /dx**2, args = (x, rt_time[i], cm, conc_step1[i]), method='lm')
    diffusivity1.append(res.x[0]*dx**2)
#print(diffusivity1)
print('Diffusivity %s' %np.mean(diffusivity1))

x = x+1.
diffusivity2 = []
for i in np.arange(1, len(conc_step)):
    res = root(equation, D /dx**2, args = (x, rt_time[i], cm, conc_step2[i]), method='lm')
    #print(res.x[0])
    diffusivity2.append(res.x[0]*dx**2)
#print(diffusivity2)
print('Diffusivity %s' %np.mean(diffusivity2))
x = domain.meshgrid()[0]

'''
a = np.where(np.array(diffusivity2)<=1.e-10)[0]
for i in np.arange(0,len(a)):
    diffusivity2.pop(a[i])
    '''
plt.figure()
plt.loglog(diffusivity0[1:-1], label = '0')
plt.loglog(diffusivity1[1:-1], label = '1')
plt.loglog(diffusivity2[1:-1], label = '2')
plt.legend()
plt.show()

plt.figure()
plt.plot(diffusivity0[500:-1], label = '0')
plt.plot(diffusivity1[500:-1], label = '1')
plt.plot(diffusivity2[500:-1], label = '2')
plt.legend()
plt.show()


'''
print(diffusivity0[-3000])
print(diffusivity1[-3000])
print(diffusivity2[-3000])
'''

print(diffusivity0[-1])
print(diffusivity1[-1])
print(diffusivity2[-1])

#%%

plt.figure()
plt.plot(rt_time[200:-1], np.abs(dport[200:-1]), label = "0.1")
plt.xlabel("Time (s)")
plt.ylabel("dCH/dt (mol/s) * 1e-15")
plt.legend()
plt.show()


plt.figure()
plt.loglog(rt_time[200:-1], np.abs(dport[200:-1]), label = "0.1")
plt.xlabel("Time (s)")
plt.ylabel("dCH/dt (mol/s) * 1e-15")
plt.legend()
plt.show()
'''
plt.figure()
plt.plot(rt_time[10:-1], np.abs(dport[10:-1]), label = "0.05")
plt.plot(rt_time[10:-1], np.abs(dport01[10:-1]), label = "0.1")
plt.plot(rt_time[10:-1], np.abs(dport001[10:-1]), label = "0.01")
plt.xlabel("Time (s)")
plt.ylabel("dCH/dt (mol/s) * 1e-15")
plt.legend()
plt.show()
'''