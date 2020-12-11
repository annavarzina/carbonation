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
#%%import phrqc

from phrqc_input import PhreeqcInput

class PhreeqcInputCSH(PhreeqcInput):
    def make_phrqc_input(self):
        self.phrqc_phases(self.c['csh'])
        self.phrqc_boundary_voxel_CSH(ca = self.c['ca_bc'], si = self.c['si_bc'])
        self.phrqc_liquid_voxel_CSH(ca = self.c['ca_liq'], si = self.c['si_liq'])
        self.phrqc_multilevel_voxel_CSH(ca = self.c['ca_mlvl'], si = self.c['si_mlvl'])        
        self.phrqc_solid_voxel()
        #print(self.phrqc_input)
        return(self.phrqc_input)
        
    def phrqc_phases(self, csh):        
        phrqc_input = [] 
        s = csh['stochiometry']
        h = s['H+']
        #oh = 2*s['Ca'] - s['Si']
        h2o =  s['H2O'] +  s['H+'] - s['Ca']
        sign1 = '+'
        if h < 0: 
            sign1 = '-'
            h *= -1
        sign2 = '+'
        if h2o < 0: 
            sign2 = '-'
            h2o *= -1
        phrqc_input.append('PHASES')  
        phrqc_input.append(csh['name'])  
        phrqc_input.append('\t(CaO)' + str(s['Ca']) +'(SiO2)'+ str(s['Si']) + \
                           '(H2O)' + str(s['H2O']) + ' ' + sign1 + ' ' + str(h) + 'H+ = ' + \
                           str(s['Ca']) + 'Ca+2 + ' + str(s['Si']) + 'SiO2 '  + sign2 +\
                           ' ' + str(h2o) + ' H2O') 
        #phrqc_input.append('\t-Vm\t' + str(csh['vm']) ) 
        phrqc_input.append('\t-log_K\t' + str(csh['log_k']) + '\n')
        phrqc_input.append('knobs') 
        phrqc_input.append('\t-iterations 8000\n') 
        
        self.phrqc_input +=  phrqc_input
        
    def phrqc_boundary_voxel_CSH(self, ca, si):
        phrqc_input = [] 
        phrqc_input.append('#boundary_solution')    
        phrqc_input.append('SOLUTION\t100001')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
            phrqc_input.append('\tSi\t' + str(si['value']) + '\n')
        else:
            pass
        phrqc_input.append('EQUILIBRIUM_PHASES\t100001\n')
        #phrqc_input.append(self.c['csh']['name'] + '\t0\t0\n')
        self.phrqc_input +=  phrqc_input
        
    def phrqc_liquid_voxel_CSH(self, ca, si):
        phrqc_input = [] 
        phrqc_input.append('#solution_liquid')    
        phrqc_input.append('SOLUTION\t100002')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
            phrqc_input.append('\tSi\t' + str(si['value'])+ '\n')
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
            phrqc_input.append('\tSi\t1\t' + str(si['value'])+ '\n')
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
        phrqc_input.append(self.c['csh']['name'] + '\t0\t0\n')
        self.phrqc_input +=  phrqc_input
        
    
    def phrqc_multilevel_voxel_CSH(self, ca, si):
        phrqc_input = [] 
        phrqc_input.append('#solution_multilevel')    
        phrqc_input.append('SOLUTION\t100004')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
            phrqc_input.append('\tSi\t' + str(si['value'])+ '\n')
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
            phrqc_input.append('\tSi\t1\t' + str(si['value'])+ '\n')
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100004')
        phrqc_input.append(self.c['csh']['name'] + '\t0\t0\n')
        #phrqc_input.append(self.c['csh']['name'] + '\t0\t'+self.c['csh_mol']['value']+'\n')
        self.phrqc_input +=  phrqc_input
        
#%% PROBLEM DEFINITION

#problem type
m = 'CSH'

#fractions = np.array([1., 0.7, 0.5, 0.3, 0.2, 0.1, 0.07, 0.05, 0.03, 
#                      0.02, 0.01, 0.007, 0.005, 0.003, 0.002, 0.001 ])

dl = 5
f = 0.1
#fractions = np.array([f])

ll = 1 #liquid lauer in front of csh
l_csh = 10 #40 #length of csh
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
scale = 100 # scale of molar volume
mvol = [51.7e-3*scale]#m3/mol # CSH Jennite and SiO2
max_pqty = fn.get_max_pqty(mvol) #mol/m3
init_conc = fn.set_init_pqty(mvol, scale = scale, porosCH = 0.45) #TODO replace porosity CSH

nn=os.path.basename(__file__)[:-3] 
nn += '_' + str(f)
fn.make_output_dir(root_dir+'\\results\\output_csh\\03_leaching')
path = root_dir+'\\results\\output_csh\\03_leaching\\' + nn + '\\'
fn.make_output_dir(path)
phase_name = m
phrqc_input = {'csh_mol':{'value': str(init_conc[0]) },
               'csh':{'name':phase_name, 
                      'stochiometry':{'Ca':0.83,
                                  'Si':1.0,
                                  'H2O':2.49,
                                  'H+':1.66},
                      'log_k':11.15}, #-11.15
               'ca_bc':{'type':'conc', 'value': '0.0'}, 
               'si_bc':{'type':'conc', 'value': '0.0'},
               'ca_mlvl':{'type':'eq', 'value': phase_name}, 
               'si_mlvl':{'type':'eq', 'value': phase_name}, 
               'ca_liq':{'type':'conc', 'value': '0'},
               'si_liq':{'type':'conc', 'value': '0'}               
               } 
 # another option ca_liq':{'type':'conc', 'value': '0'} or ca_liq':{'type':'eq', 'value': 'portlandite'}

phrqc = PhreeqcInputCSH(phrqc_input)            
phrqc.save_phrqc_input(root_dir, nn)   

#%% VALUES
pqty = init_conc[0] * (domain.nodetype == ct.Type.MULTILEVEL)

slabels = fn.set_labels(domain, m) 
D_CSH =    1.e-11
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
dp['database']='cemdata18.dat'

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
prev_csh = deepcopy(rt.solid.CSH.c > 0)
rt.dissolution_time = []
rt_csh = [] 
rt_csh.append(np.sum(rt.solid.CSH.c)*scale)
rt_time = []
rt_time.append(rt.time)
dcsh = [] #dissolution rate
dcsh.append(0)
Ts = 1000 
nitr = 2
length = 0
while rt.iters <= nitr: # length < dl: #rt.time <=Ts: #  
    rt.advance()             
    if (np.sum(rt.solid.CSH.c > 0) < np.sum(prev_csh)):
        prev_csh = deepcopy(rt.solid.CSH.c > 0)
        rt_time.append(rt.time*scale)
        rt.dissolution_time.append(rt.time)
        rt_csh.append(np.sum(rt.solid.CSH.c)*scale)
        if(len(rt_csh)<2):
            dcsh.append(0)
        else:
            dcsh.append((rt_csh[-1]-rt_csh[-2])/(rt_time[-1]-rt_time[-2])) #mol/s/um2
        length +=1
        print(length)
    itr += 1 
    
#%% SIMULATION TIME
print("Fraction %s done" %str(f))
print("Time to dissolve %s" %str(rt.time*scale))
simulation_time = time.time()-it
fn.print_time(simulation_time, rt)
#%%  SAVE
'''
dim =  10**(-15)*10**12 #convert to mol/l/m2/s
vol =  0.04 * 1e-3 / 60 #volume
np.save(path + 'dis_time', rt.dissolution_time )
np.save(path + 'dCSH', [np.abs(d*dim) for d in dcsh]) # mol/l/m2/s
np.save(path + 'dCSH2', [np.abs(d*dim*vol) for d in dcsh]) #mol/m2/s
np.save(path + 'logK1', [np.log10(np.abs(d*dim)) for d in dcsh[1:]]) #mol/m2/s
np.save(path + 'logK2', [np.log10(np.abs(d*dim*vol)) for d in dcsh[1:]]) #mol/m2/s
np.save(path + 'time', rt_time)
np.save(path + 'CH', rt_csh)
np.save(path + 'Ca_prof', rt.fluid.Ca.c[1,:]+ np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.phrqc.poros[1,:]))
np.save(path + 'Si_prof', rt.fluid.Si.c[1,:]+ np.array(rt.fluid.Si._ss[1,:])/np.array(rt.phrqc.poros[1,:]))
'''
#%% Rate
#dim = 10**3*10**(-15)*10**8 # convert to mmol/l/cm2/s
dim =  10**(-15)*10**12 #convert to mol/l/m2/s
vol = 0.04 * 1e-3 / 60 #l/s from TB paper
print([d*dim for d in dcsh])
print([d*dim*vol for d in dcsh])
print([np.log10(np.abs(d*dim)) for d in dcsh[1:]])
print([np.log10(np.abs(d*dim*vol)) for d in dcsh[1:]])
#%%  PRINT
'''
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
'''