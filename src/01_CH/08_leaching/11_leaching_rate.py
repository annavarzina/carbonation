# -*- coding: utf-8 -*-
'''

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
from copy import deepcopy
#import phrqc
#%% PROBLEM DEFINITION
__doc__= """ 

"""
#problem type
m = 'CH' #or 'CSH' #TODO case for cement

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
#%% LOOP
fractions = np.array([1.,0.005 ])
dl = 2
for f in fractions:
    #%% GEOMETRY
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
    phrqc_input = {'ca_bc':{'type':'conc', 'value': '0.0'}, 
                   'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
                   'ca_liq':{'type':'conc', 'value': '0'}} # another option ca_liq':{'type':'conc', 'value': '0'} or ca_liq':{'type':'eq', 'value': 'portlandite'}
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
                'subgrid': {'fraction':f,
                            'poros': False}, # fraction of interface cell number or None = porosity
                'app_tort':{'degree': app_tort_degree}, #TODO
                'dx': dx, 
                'Dref':D_high
                }
              
    bc_params = {'solution_labels':{'left':100001}, 
                'top':['flux', 0.0],
                'bottom':['flux', 0.0],
                'left':['flux', 0.0],
                'right':['flux', 0.0],}               
    #%% PARAMETERS (DOMAIN, BC, SOLVER)
    domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                         input_file = root_dir + \
                                         '\\phreeqc_input\\' + nn + '.phrq')#'CH_CC-nat.phrq'
    #bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
    solver_params = fn.set_solver_params(tfact = 1/6., smart_thres = 1e-8)# optional values, for time step (if tfact => tfactbased tau)
    solver_params['phrqc_flags']['smart_run']=False
    domain.nodetype[domain.nodetype == ct.Type.MULTILEVEL_CH] = ct.Type.MULTILEVEL
    fn.save_settings(settings, bc_params, solver_params, path, nn)
    #%% INITIATE THE SOLVER
    rt= rt1.LeachingRT('MultilevelAdvectionDiffusion',  domain, 
                              domain_params, bc_params, solver_params,
                              settings) 
    
    #%% RUN SOLVER
    it=time.time()
    itr = 0 
    j = 0
    l = ll
    prev_nodetype = deepcopy(rt.nodetype)
    rt.dissolution_time = []
    rt_port = [] 
    rt_port.append(np.sum(rt.solid.portlandite.c))
    rt_time = []
    rt_time.append(rt.time)
    dport = []
    dport.append(0)
    while  len(rt.dissolution_time)<dl:#rt.time <=Ts: # itr < nitr: #
        rt.advance()             
        if (~np.all(rt.nodetype == prev_nodetype)):
            prev_nodetype = deepcopy(rt.nodetype)
            rt_time.append(rt.time*scale)
            rt.dissolution_time.append(rt.time)
            rt_port.append(np.sum(rt.solid.portlandite.c)*scale)
            if(len(rt_port)<2):
                dport.append(0)
            else:
                dport.append((rt_port[-1]-rt_port[-2])/(rt_time[-1]-rt_time[-2]))
        itr += 1 

    #%% SIMULATION TIME
    print("Fraction %s done" %str(f))
    print("Time to dissolve %s" %str(rt.time*scale))
    simulation_time = time.time()-it
    fn.print_time(simulation_time, rt)
      
    #%%  SAVE
    np.save(path + 'dis_time', rt.dissolution_time )
    np.save(path + 'dCH', dport)
    np.save(path + 'time', rt_time)
    np.save(path + 'CH', rt_port)
    np.save(path + 'Ca_prof', rt.fluid.Ca.c[1,:]+ np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:]))
    

    #%% PLOT
      
    
    plt.figure()
    plt.plot(rt_time[2:-1], np.abs(dport[2:-1]), "r.")
    plt.xlabel("Time")
    plt.ylabel("Portlandite rate")
    plt.show()
    
    t = np.array(rt.dissolution_time)*scale
    x = np.arange(1,len(rt.dissolution_time)+1)*1e-6
    
    plt.figure()
    D = x**2/2/t
    plt.plot(t, D, "r.")
    plt.ylabel("D (m2/s)")
    plt.xlabel("Time (s)")
    plt.show()
    print("Leaching diffusion %s" %str(D[-1]))
    
    
    '''
    plt.figure()
    plt.plot(t, x)
    plt.ylabel("Dissolved length (um)")
    plt.xlabel("Time")
    plt.show()
    '''


    