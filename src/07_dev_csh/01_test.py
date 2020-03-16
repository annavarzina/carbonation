# -*- coding: utf-8 -*-
'''
Molar volume = 100
Liquid layer = 5 um
Fraction = 0.004
'''

#%% PYTHON MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
sys.path.append(src_dir)
import matplotlib.pylab as plt
from copy import deepcopy
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import yantra
import cell_type as ct # change the path to cell_type file
import dev_rt as rt 
import misc_func as fn
#import phrqc
#%% FUNCTIONS
def set_phrqc_input(p, ptype ='CH'):
    phrqc_input = [] 
    if(ptype=='CSH'):
        phrqc_input.append('PHASES')  
        phrqc_input.append('CSHQ_TobH')  
        phrqc_input.append('\t(CaO)0.66666667(SiO2)1(H2O)1.5 = 0.66666667Ca++ + 1 SiO(OH)3- + 0.33333334OH- -0.16666667 H2O') 
        phrqc_input.append('\tlog_K -6.190832') 
        phrqc_input.append('CSHQ_TobD') 
        phrqc_input.append('\t(CaO)0.8333333333(SiO2)0.6666666667(H2O)1.8333333333 = 0.8333333333 Ca++ + 0.6666666667 SiO(OH)3- + 0.99999999990 OH- + 0.3333333333 H2O') 
        phrqc_input.append('\tlog_K -6.8995533') 
        phrqc_input.append('CSHQ_JenH') 
        phrqc_input.append('\t(CaO)1.3333333333(SiO2)1(H2O)2.1666666667 = 1.3333333333 Ca++ + 1 SiO(OH)3- + 1.6666666667 OH- -0.1666666667 H2O') 
        phrqc_input.append('\tlog_K -10.96765') 
        phrqc_input.append('CSHQ_JenD') 
        phrqc_input.append('\t(CaO)1.5(SiO2)0.6666666667(H2O)2.5 = 1.5 Ca++ + 0.6666666667 SiO(OH)3- + 2.3333333333 OH- + 0.3333333333 H2O') 
        phrqc_input.append('\tlog_K -10.47635') 
        phrqc_input.append('knobs') 
        phrqc_input.append('\t-iterations 8000') 

    # Boundary
    phrqc_input.append('#boundary_solution')    
    phrqc_input.append('SOLUTION\t100001')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    if(p['c_bc']['type'] == 'conc'):
        phrqc_input.append('\tC\t' + str(p['c_bc']['value']) + '\n')
    elif(p['c_bc']['type'] == 'pco2'):
        phrqc_input.append('\tC\t1\tCO2(g)\t-' + str(p['c_bc']['value']) + '\n')
    phrqc_input.append('EQUILIBRIUM_PHASES\t100001\n')
    
    #Liquid
    phrqc_input.append('#solution_liquid')    
    phrqc_input.append('SOLUTION\t100002')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    if(p['c_liq']['type'] == 'conc'):
        phrqc_input.append('\tC\t' + str(p['c_liq']['value']))
    elif(p['c_liq']['type'] == 'eq'):
        phrqc_input.append('\tC\t1\t' + str(p['c_liq']['value']))
    else:
        pass        
    if(p['ca_liq']['type'] == 'conc'):
        phrqc_input.append('\tCa\t' + str(p['ca_liq']['value']))
    elif(p['ca_liq']['type'] == 'eq'):
        phrqc_input.append('\tCa\t1\t' + str(p['ca_liq']['value']))
    else:
        pass        
    phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
    phrqc_input.append('portlandite\t0\t0')
    phrqc_input.append('calcite\t0\t0\n')
    
    #Multilevel
    phrqc_input.append('#solution_multilevel')    
    phrqc_input.append('SOLUTION\t100003')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    if(p['c_mlvl']['type'] == 'conc'):
        phrqc_input.append('\tC\t' + str(p['c_mlvl']['value']))
    elif(p['c_mlvl']['type'] == 'eq'):
        phrqc_input.append('\tC\t1\t' + str(p['c_mlvl']['value']))
    else:
        pass              
    if(p['ca_mlvl']['type'] == 'conc'):
        phrqc_input.append('\tCa\t' + str(p['ca_mlvl']['value']))
    elif(p['ca_mlvl']['type'] == 'eq'):
        phrqc_input.append('\tCa\t1\t' + str(p['ca_mlvl']['value']))
    else:
        pass        
    phrqc_input.append('EQUILIBRIUM_PHASES\t100003')
    phrqc_input.append('portlandite\t0\t1')
    phrqc_input.append('calcite\t0\t0\n')    
    if(ptype=='CSH'):
        phrqc_input.append('#solution_csh_multilevel') 
        phrqc_input.append('SOLUTION\t100004') 
        phrqc_input.append('\t-water\t0.448230266981165') 
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\tpH\t12\tcharge') 
        phrqc_input.append('\tCa\t1.955e-002') 
        phrqc_input.append('\tSi\t3.018e-005') 
        phrqc_input.append('SOLID_SOLUTIONS\t100004') 
        phrqc_input.append('Tob_jen_ss') 
        phrqc_input.append('\t-comp\tCSHQ_TobH\t0.1041') 
        phrqc_input.append('\t-comp\tCSHQ_TobD\t2.5050') 
        phrqc_input.append('\t-comp\tCSHQ_JenH\t2.1555') 
        phrqc_input.append('\t-comp\tCSHQ_JenD\t3.2623') 
        phrqc_input.append('EQUILIBRIUM_PHASES\t100004')
        phrqc_input.append('portlandite\t0\t0')
        phrqc_input.append('calcite\t0\t0\n')
            
    #Solid
    phrqc_input.append('#solution_solid')    
    phrqc_input.append('SOLUTION\t100005')
    phrqc_input.append('\t-water\t1\n')
    return phrqc_input

def set_mvols(mvol, scale, ptype = 'CSH'):
    '''
    mvol is dictionary with fields CH, CC and CSH types (TobH, TobD, JenH, JenD)
    '''
    ch_type = (ptype=='CH')
    csh_type = (ptype=='CSH')
    for k in mvol:
        mvol[k] = mvol[k]*scale
    default = {'CH': 33.10e-3*scale, 'CC': 36.90e-3*scale,
               'CSH_TobH': 55.30e-3*scale, 'CSH_TobD': 47.95e-3*scale,
               'CSH_JenH': 75.63e-3*scale, 'CSH_JenD': 80.58e-3*scale}
    merged = default.copy()   
    merged.update(mvol)
    if ch_type:
        return [merged['CH'], merged['CC']]
    elif csh_type:
        return [merged['CH'], merged['CC'],
                merged['CSH_TobH'], merged['CSH_TobD'],
                merged['CSH_JenH'], merged['CSH_JenD']]
    else:
        print('Define ptype as \'CH\' or \'CSH\'. ')    
        return []
    
def set_init_pqty(mvol, porosCH = 0.1, wc = 0.45): #TODO check for CSH
    maxCH = fn.get_max_pqty(mvol)[0]
    initCH = (1 - porosCH) * maxCH
    initCC = 0.0
    init_conc = [initCH, initCC]
    if (len(mvol) == 6):
        initCSH_TobH = 0.1041#(1-init_porosCSH) *maxCSH_TobH
        initCSH_TobD = 2.5050#(1-init_porosCSH) *maxCSH_TobD
        initCSH_JenH = 2.1555#(1-init_porosCSH) *maxCSH_JenH
        initCSH_JenD = 3.2623#(1-init_porosCSH) *maxCSH_JenD
        # set different init if wc ~= 0.45
        init_conc = [initCH, initCC, 
                     initCSH_TobH, initCSH_TobD, 
                     initCSH_JenH, initCSH_JenD]
    return init_conc


#%% PROBLEM DEFINITION
__doc__= """ 
Reference:
    Lime solution \theta = 0.05
    PCO2 = 3.4
    Molar volume = 100
    Liquid layer = 5 um
    Fraction = 0.004
"""
#problem type

#%% GEOMETRY
ll = 5 #liquid lauer in front of portlandite
l_ch = 25 #length of portlandite
lx = (l_ch+ll)*1.0e-6
ly = 2.0e-6
dx = 1.0e-6

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, :] = ct.Type.MULTILEVEL

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()

#%%  VALUES
m = 'CSH' #or 'CSH'
nn=os.path.basename(__file__)[:-3]
fn.make_output_dir(root_dir+'\\results\\temp\\01_develop\\')
path = root_dir+'\\results\\temp\\01_develop\\' + nn + '\\'
fn.make_output_dir(path)


phrqc_input = {'c_bc':{'type':'pco2', 'value': 3.4},
               'c_mlvl':{'type':'conc', 'value': '0'}, 
               'c_liq':{'type':'conc', 'value': '0'},
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'},    
               'ca_liq':{'type':'conc', 'value': '0'}}

phrqc = set_phrqc_input(phrqc_input, ptype = m)            
fn.save_phrqc_input(phrqc,root_dir, nn)  
#%%
scale = 1 # scale of molar volume
init_porosCH = 0.05 #initial porosity of portlandite nodes
mvol = set_mvols({}, scale, ptype = m) #m3/mol

max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = set_init_pqty(mvol, init_porosCH)
pqty = fn.get_pqty(init_conc, domain)

slabels = fn.set_labels(domain, m)     
D = 1.0e-09 # default diffusion coefficient in pure liquid
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort_degree = 1./3.
app_tort = 1. * porosity ** app_tort_degree

settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            'dissolution':'multilevel', #'multilevel'/'subgrid'
            'active_nodes': 'smart', # 'all'/'smart'/
            'diffusivity':{'border': D, ##diffusivity at border
                           'CH': ('const', 1e-15), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CC': ('inverse', 1e-12), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           }, 
            'pcs_mode': {'pcs': True, #Pore-Size Controlled Solubility concept
                         'pores': 'block', #'block'/'cylinder'
                         'int_energy': 0.1, # internal energy
                         'pore_size': 0.01*dx, # threshold radius or distance/2
                         'crystal_size': 0.5*dx, # crystal or pore length
                         'pore_density': 2000, #pore density per um3 - only for cylinder type
                         }, 
            'subgrid': {'fraction':0.004}, # fraction of interface cell number or None = porosity
            'app_tort':{'degree': 1./3.}, #TODO
            'velocity': False, 
            'bc': phrqc_input['c_bc'],
            'dx': dx, 
            'Dref':D
            }
tfact_default = 1./6./1#*init_porosCH
domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir +'\\phreeqc_input\\' + nn + '.phrq')
bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = tfact_default, smart_thres = 1e-8, cphi_fact = 1/3.)
#fn.save_settings(settings, bc_params, solver_params, path, nn)
#csh=yantra.PhrqcReactiveTransport('MultilevelDiffusion',domain,domain_params,bc_params,solver_params)
csh=rt.CSHCarbonation('MultilevelDiffusion',domain,domain_params,bc_params,solver_params, settings)
#%%run model
n=100
while csh.iters<n:
    csh.advance()    
    
#%%
