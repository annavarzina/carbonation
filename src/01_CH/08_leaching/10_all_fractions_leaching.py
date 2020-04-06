# -*- coding: utf-8 -*-
'''
Compare leaching rate for different fractions
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
fractions = np.array([1., 0.7, 0.5, 0.3, 0.2, 0.1, 0.07, 0.05, 0.03, 
                      0.02, 0.01, 0.007, 0.005, 0.003, 0.002, 0.001 ])
#final_time = np.array([500, 500, 500, 1000, 1000, 1000, 1000, 1500, 1500,
#                       1500, 2000, 2000, 2500, 3000, 3500, 4000])  
#fractions = np.array([0.007, 0.005, 0.003, 0.002, 0.001 ])
dl = 20
#final_time = np.array([10, 20])
#for f, ft in zip(fractions, final_time):
for f in fractions:
    #%% GEOMETRY
    ll = 1 #liquid lauer in front of portlandite
    l_ch = 40 #length of portlandite
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
    nitr =5000#2000
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
    
    #%% PRINT
    '''
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
    '''
    
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

#%% COMPARE    
fractions = np.array([1., 0.7, 0.5, 0.3, 0.2, 0.1,
                      0.07, 0.05, 0.03,0.02, 0.01,
                      0.007, 0.005, 0.003, 0.002, 0.001 ])
time = {}
ch = {}
ca ={}
dch = {}
dissolution = {}
nn="03_subgrid_leaching_depth"#os.path.basename(__file__)[:-3] 
for f in fractions:
    path = root_dir+'\\results\\output\\10_subgrid_leaching\\' + nn + str(f)+'\\'
    dissolution[nn+ str(f)] = np.load(path +'dis_time'+'.npy')
    ch[nn+ str(f)] = np.load(path +'CH'+'.npy')
    ca[nn+ str(f)] = np.load(path +'Ca_prof'+'.npy')
    dch[nn+ str(f)] = np.load(path +'dCH'+'.npy')
    time[nn+ str(f)] = np.load(path +'time'+'.npy')
    
#%% PLOT  
plt.figure()
for f in fractions:    
    t = np.array(dissolution[nn + str(f)])*scale
    x = np.arange(1,len(dissolution[nn + str(f)])+1)
    plt.plot(t, x, label = f)
plt.ylabel("Dissolved length (um)")
plt.xlabel("Time")
plt.legend()
plt.show()

plt.figure()
for f in fractions:  
    t = np.array(dissolution[nn + str(f)])*scale
    x = np.arange(1,len(dissolution[nn + str(f)])+1)*1e-6
    D = (x**2)/2/t
    plt.plot(t, D, label = f)
plt.ylabel("D (m2/s)")
plt.xlabel("Time (s)")
plt.legend()
plt.show()

#%%
d1 = []
dn = []
for f in fractions:
    t = np.array(dissolution[nn + str(f)])*scale
    x = np.arange(1,len(dissolution[nn + str(f)])+1)*1e-6
    D = (x**2)/2/t
    d1.append(D[0])
    dn.append(D[-1])
    
plt.figure()
plt.plot(fractions, d1, label = "1")
#plt.plot(fractions, dn, label = "20")
plt.ylabel("D (m2/s)")
plt.xlabel("Fraction")
plt.yscale("log")
plt.legend()
plt.show()
#%% CH    
plt.figure()
for f in fractions:  
    plt.plot(np.array(time[nn + str(f)]), ch[nn + str(f)]*scale, label = f)
plt.ylabel("CH mol/l")
plt.xlabel("Time")
plt.legend()
plt.show()

#%% CH rate    
plt.figure()
for f in fractions:  
    plt.plot(np.array(time[nn + str(f)])[1:-1], np.abs(dch[nn + str(f)])[1:-1]*scale, label = f)
plt.ylabel("CH rate (*10^-15 mol/s)")
plt.xlabel("Time")
plt.yscale("log")
plt.legend()
plt.show()

plt.figure()
for f in fractions:  
    plt.plot( ca[nn + str(f)], label = f)
plt.ylabel("Ca (mol/l)")
plt.xlabel("X (um)")
plt.legend()
plt.show()
#%% CH rate
r1 = []
rn = []
for f in fractions:    
    r1.append(dch[nn + str(f)][1]*scale)
    rn.append(dch[nn + str(f)][-11]*scale)
    
    
plt.figure()
plt.plot(fractions, np.abs(r1), label = "1")
plt.plot(fractions, np.abs(rn), label = "20")
plt.ylabel("dCH (*10^-15 mol/s)")
plt.xlabel("Fraction")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.show()
#%%

from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
def diff_predict(f, c0, c1, c2):
    d = c0-c1/(f+c2)
    #d = c0*(1-np.exp(-f*c1))
    return d
f = 0.1
t = np.array(dissolution[nn + str(f)])*scale
x = np.arange(1,len(dissolution[nn + str(f)])+1)*1e-6
D = x**2/2/t
c, cov = curve_fit(diff_predict, t, D)
print(c)

diff_opt = diff_predict(t, c[0],c[1],c[2])
print('R2: ', r2_score(diff_opt, D))


plt.plot(t, D, 'r.', label = "Result")
plt.plot(t, diff_opt, 'b-', label = "Fitting")
plt.xlabel('Time')
plt.ylabel(r'Diffusivity ($m^2/s$)')
plt.legend()
plt.show()

#%% Portlandite to ppm
amCH = 74.093 #mol/l
amCa = 40.078 #g/mol
def mol2ppm(c, am):
    p = c * am * 1000
    return p #ppm
#ppm/s/um2
rate = mol2ppm(np.abs(r1), amCH)*3600 *10**6 *10**-15
print(rate)
rate_n = mol2ppm(np.abs(rn), amCH)*3600 *10**6 *10**-15
print(rate_n)
plt.figure()
plt.loglog(fractions, rate, label = "1")
plt.loglog(fractions, rate_n, label = "20")
plt.ylabel("dCH (ppm/h/mm2)")
plt.xlabel("Fraction")
plt.legend()
plt.show()

#%% Ca to ppm 
ca_n = []
for f in fractions:    
    ca_n.append(ca[nn + str(f)][-21])
print(ca_n)

ca_ppm = mol2ppm(np.abs(ca_n), amCa)
print(ca_ppm)
    