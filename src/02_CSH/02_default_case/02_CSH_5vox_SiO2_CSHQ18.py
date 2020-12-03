# -*- coding: utf-8 -*-
'''
Molar volume = 100
Liquid layer = 5 um
Fraction = 0.004
'''

#%% PYTHON MODULES
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(root_dir)
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(src_dir)
ch_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ch_dir)
import matplotlib.pylab as plt
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import yantra
import cell_type as ct # change the path to cell_type file
import rt_carb_csh_sio2 as rt
import misc_func as fn
import time
#%% PROBLEM DEFINITION
__doc__= """ 
...
"""

#%% PHREEQC input

from phrqc_input import PhreeqcInput

class PhreeqcInputCSHQ(PhreeqcInput):
    
    def make_phrqc_input(self):
        self.phrqc_phases()
        self.phrqc_boundary_voxel(self.c['c_bc'])
        self.phrqc_liquid_voxel(self.c['c_liq'], self.c['ca_liq'])
        self.phrqc_multilevel_voxel_CH(self.c['c_mlvl'], self.c['ca_mlvl'])        
        self.phrqc_multilevel_voxel_CSH()
        self.phrqc_solid_voxel()
        return(self.phrqc_input)
        
    def phrqc_phases(self):
        phrqc_input = [] 
        
        phrqc_input.append('PHASES')  
        phrqc_input.append('CSHQ_TobH')  
        phrqc_input.append('\t(CaO)0.6667(SiO2)1(H2O)1.5 + 1.3334H+ = 0.6667Ca+2 + 2.1667H2O + SiO2') 
        phrqc_input.append('\t-Vm\t55') 
        phrqc_input.append('\t-analytical_expression\t-12.519254 0 2163.381583 5.476331 0 0 0') 
        phrqc_input.append('\t-log_K\t8.286642\n') 
		
        phrqc_input.append('CSHQ_TobD') 
        phrqc_input.append('\t((CaO)1.25(SiO2)1(H2O)2.75)0.6667 + 1.66675H+ = 0.833375Ca+2 + 2.6668H2O + 0.6667SiO2') 
        phrqc_input.append('\t-Vm\t48') 
        phrqc_input.append('\t-analytical_expression\t-10.916344 0 3959.367696 4.563888 0 0 0') 
        phrqc_input.append('\t-log_K\t13.655314\n') 
		
        phrqc_input.append('CSHQ_JenH') 
        phrqc_input.append('\t(CaO)1.3333(SiO2)1(H2O)2.1667 + 2.6666H+ = 1.3333Ca+2 + 3.5H2O + SiO2') 
        phrqc_input.append('\t-Vm\t76') 
        phrqc_input.append('\t-analytical_expression\t-17.10944 0 6470.553982 7.107847 0 0 0') 
        phrqc_input.append('\t-log_K\t22.179305\n') 
		
        phrqc_input.append('CSHQ_JenD') 
        phrqc_input.append('\t(CaO)1.5(SiO2)0.6667(H2O)2.5 + 3H+ = 1.5Ca+2 + 4H2O + 0.6667SiO2') 
        phrqc_input.append('\t-Vm\t81') 
        phrqc_input.append('\t-analytical_expression\t-15.591756 0 8609.739692 6.24251 0 0 0') 
        phrqc_input.append('\t-log_K\t28.730362\n') 
		
        phrqc_input.append('SiO2am') 
        phrqc_input.append('\tSiO2 = SiO2')
        phrqc_input.append('\t-Vm\t29') 
        phrqc_input.append('\t-analytical_expression\t0 0 -809.189752 0 0 0 0') 
        phrqc_input.append('\t-log_K\t-2.714066\n') 
		
        phrqc_input.append('calcite') 
        phrqc_input.append('\tCaCO3 = CO3-2 + Ca+2')
        phrqc_input.append('\t-Vm\t36.934') 
        phrqc_input.append('\t-analytical_expression\t130.276347 0 -5689.203921 -48.36444 0 0 0') 
        phrqc_input.append('\t-log_K\t-8.479966\n') 
		
        phrqc_input.append('portlandite') 
        phrqc_input.append('\tCa(OH)2 + 2H+ = Ca+2 + 2H2O')
        phrqc_input.append('\t-Vm\t33.06') 
        phrqc_input.append('\t-analytical_expression\t-11.299363 0 7301.394065 3.883957 0 0 0') 
        phrqc_input.append('\t-log_K\t22.799937\n') 
		
        phrqc_input.append('knobs') 
        phrqc_input.append('\t-iterations 8000\n') 
        
        self.phrqc_input +=  phrqc_input
    
    def phrqc_multilevel_voxel_CH(self, c, ca):
        phrqc_input = [] 
        phrqc_input.append('#solution_multilevel')    
        phrqc_input.append('SOLUTION\t100003')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(c['type'] == 'conc'):
            phrqc_input.append('\tC\t' + str(c['value']))
        elif(c['type'] == 'eq'):
            phrqc_input.append('\tC\t1\t' + str(c['value']))
        else:
            pass              
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100003')
        phrqc_input.append('portlandite\t0\t1')
        phrqc_input.append('calcite\t0\t0\n')
        
        self.phrqc_input +=  phrqc_input
                           
    def phrqc_multilevel_voxel_CSH(self):
        phrqc_input = [] 
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
        phrqc_input.append('\t-comp\tSiO2am\t0.0') 
        phrqc_input.append('EQUILIBRIUM_PHASES\t100004')
        phrqc_input.append('portlandite\t0\t0')
        phrqc_input.append('calcite\t0\t0\n')
        
        self.phrqc_input +=  phrqc_input


#%% GEOMETRY
ll = 4 #liquid lauer in front of portlandite
l_ch = 5 #length of CSH
lx = (l_ch+ll+1)*1.0e-6
ly = 2.0e-6
dx = 1.0e-6

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, ll+1: ll+l_ch+1] = ct.Type.MULTILEVEL
domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()

#%%  VALUES
m = 'CSH' #or 'CSH'
nn=os.path.basename(__file__)[:-3]
fn.make_output_dir(root_dir+'\\results\\output_csh\\01_default\\')
path = root_dir+'\\results\\output_csh\\01_default\\' + nn + '\\'
fn.make_output_dir(path)


phrqc_input = {'c_bc':{'type':'pco2', 'value': 3.0},
               'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
               'c_liq':{'type':'eq', 'value': 'calcite'},
               'ca_mlvl':{'type':'eq', 'value': 'calcite'},    
               'ca_liq':{'type':'eq', 'value': 'calcite'}}

phrqc = PhreeqcInputCSHQ(phrqc_input)            
phrqc.save_phrqc_input(root_dir, nn)   
#%%
scale = 100# scale of molar volume
init_porosCH = 0.05 #initial porosity of portlandite nodes
mvol = [33.10e-3*scale, 36.90e-3*scale,
        55.30e-3*scale, 47.95e-3*scale,
        75.63e-3*scale, 80.58e-3*scale, 
        27.3e-3 *scale]

max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = [(1 - init_porosCH) * max_pqty[0], 0.0,
              0.1041/scale, 2.5050/scale, 
              2.1555/scale, 3.2623/scale, 0.0]

pqty_CH = init_conc[0] * (domain.nodetype == ct.Type.MULTILEVEL_CH)    
pqty_CC = init_conc[1] * np.ones(domain.nodetype.shape)   
pqty_CSHQ_TobH = init_conc[2] * (domain.nodetype == ct.Type.MULTILEVEL) 
pqty_CSHQ_TobD = init_conc[3] * (domain.nodetype == ct.Type.MULTILEVEL) 
pqty_CSHQ_JenH = init_conc[4] * (domain.nodetype == ct.Type.MULTILEVEL) 
pqty_CSHQ_JenD = init_conc[5] * (domain.nodetype == ct.Type.MULTILEVEL) 
pqty_SiO2 = init_conc[6] * np.ones(domain.nodetype.shape)
pqty = [pqty_CH, pqty_CC, 
        pqty_CSHQ_TobH, pqty_CSHQ_TobD, 
        pqty_CSHQ_JenH, pqty_CSHQ_JenD,
        pqty_SiO2]
        

slabels = fn.set_labels(domain, m)     
D = 1.0e-09 # default diffusion coefficient in pure liquid
porosity = fn.get_porosity(domain, pqty, mvol, m)
app_tort_degree = 1./3.
app_tort = 1. * porosity ** app_tort_degree

settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
            'dissolution':'multilevel', #'multilevel'/'subgrid'
            'active_nodes': 'all', # 'all'/'smart'/ better use all to avoid errors
            'diffusivity':{'border': D, ##diffusivity at border
                           'CH': ('const', 1e-15), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CC': ('inverse', 1e-11), # fixed diffusivity in portlandite node 'archie'/'const'/'inverse'
                           'CSH': ('const', 1e-11), # fixed diffusivity in CSH node 'archie'/'const'/'inverse'
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
domain_params['database'] = 'cemdata18.dat'
bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
solver_params = fn.set_solver_params(tfact = tfact_default, smart_thres = 1e-8, cphi_fact = 1/3.)
#fn.save_settings(settings, bc_params, solver_params, path, nn)
#csh=yantra.PhrqcReactiveTransport('MultilevelDiffusion',domain,domain_params,bc_params,solver_params)
csh=rt.CSH_Carbonation('MultilevelDiffusion',domain,domain_params,bc_params,solver_params, settings)
#%% results dict

plist =  [(1,n) for n in np.arange(3, 6)] #points list
pavglist = ['avg_poros', 'pH', 'avg_D_eff', 'sum_vol', 'precipitation', #argument list
            'dissolution', 'portlandite_cells', 'calcite_cells'] 
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)

Ts  = 30# 36 * 3 #s
Ts = Ts/scale + 0.001
N = Ts/csh.dt
N_res = 1e+4
S = max(1,int(N/N_res))
step = max(int(Ts/10.),1)
time_points = np.concatenate((np.arange(0, step, step/10.), np.arange(step, Ts+step, step)))
'''
dt = csh.dt
Ts =  1.#seconds
Ts = Ts/scale + 0.001
N_res = 1e+4
S = max(1,int(N/N_res))
'''
#%% run
#n=1000
#'''
j = 0
it=time.time()
nitr = 30
while csh.time <=Ts: #csh.iters <= nitr: #
    if(True):
        if ( (csh.time <= time_points[j]) and ((csh.time + csh.dt) > time_points[j]) ):  
            print(time_points[j])
            #fn.save_figures_minerals(carb_rt,  max_pqty, time_points[j], path, nn, ptype=m)  
            #fn.save_figures_mols(carb_rt, time_points[j], path, nn, ptype=m) 
            j +=1
    csh.advance()    
    results = fn.append_results(csh, results, step = S )
    
simulation_time = time.time()-it
fn.print_time(simulation_time, csh)
#'''
#%%
'''
fn.plot_species(results, names=[])#['calcite']
fn.plot_avg(results, names=['avg_poros', 'avg_D_eff'])
fn.plot_fields(csh, names=['calcite'],fsize=(15,1)) #,'Ca','Si'
fn.plot_points(results, names=['calcite', 'poros', 'Ca','pH'])
'''
#%% plot ca/si against density
'''
plt.figure()
plt.plot(results['Ca_Si'], results['csh_density'])
plt.legend()
plt.ylabel('CSH density')
plt.xlabel('C/S')
plt.show()
'''

#%% PRINT
#'''
print('Ca +ss %s' %str(np.array(csh.fluid.Ca.c[1,:]) + np.array(csh.fluid.Ca._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('Si +ss %s' %str(np.array(csh.fluid.Si.c[1,:]) + np.array(csh.fluid.Si._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('C +ss %s' %str(np.array(csh.fluid.C.c[1,:]) + np.array(csh.fluid.C._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('H +ss %s' %str(np.array(csh.fluid.H.c[1,:]) + np.array(csh.fluid.H._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('O +ss %s' %str(np.array(csh.fluid.O.c[1,:]) + np.array(csh.fluid.O._ss[1,:])/np.array(csh.phrqc.poros[1,:])))
print('CH %s' %str(np.array(csh.solid.portlandite.c[1,:])))
print('dCH %s' %str(np.array(csh.phrqc.dphases['portlandite'][1,:])))
print('CC %s' %str(np.array(csh.solid.calcite.c[1,:])))
print('dCC %s' %str(np.array(csh.phrqc.dphases['calcite'][1,:])))
print('Vol %s' %str(np.array(csh.solid.vol[1,:])))
print('D %s' %str(np.array(csh.fluid.C.De[1,:])))
print('SI %s' %str(np.array(csh.solid.target_SI[1,:])))
print('pH %s' %str(np.array(csh.phrqc.selected_output()['pH'][1,:])))
print('poros %s' %str(np.array(csh.solid.poros[1,:])))
print('phrqc poros %s' %str(np.array(np.array(csh.phrqc.poros[1,:]))))
print('CSHQ TobD %s' %str(np.array(csh.solid.CSHQ_TobD.c[1,:])))
print('CSHQ TobH %s' %str(np.array(csh.solid.CSHQ_TobH.c[1,:])))
print('CSHQ JenD %s' %str(np.array(csh.solid.CSHQ_JenD.c[1,:])))
print('CSHQ JenH %s' %str(np.array(csh.solid.CSHQ_JenH.c[1,:])))
print('SiO2 %s' %str(np.array(csh.solid.SiO2am.c[1,:])))
#'''
#%%  SAVE
'''
fn.save_obj(results, path + str(nn) +'_results')
np.save(path + 'SI', csh.phrqc.selected_output()['SI_calcite'] )
np.save(path + 'pH', csh.phrqc.selected_output()['pH'] )
np.save(path + 'Ca', csh.phrqc.selected_output()['Ca'] )
np.save(path + 'C', csh.phrqc.selected_output()['C'] )
np.save(path + 'CH', csh.phrqc.selected_output()['portlandite'] )
np.save(path + 'CC', csh.phrqc.selected_output()['calcite'] )
np.save(path + 'De', csh.fluid.Ca.De )
np.save(path + 'poros', csh.fluid.Ca.poros)
'''