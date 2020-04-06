'''
Miscellaneous functions for the carbonation solver
'''

import numpy as np
from copy import deepcopy
import sys    
import os
import matplotlib.pylab as plt

import pickle
import json
import cell_type as ct 

from yantra._pyevtk.hl  import pointsToVTK
from yantra._pyevtk.hl  import gridToVTK
from yantra._pyevtk.hl  import imageToVTK
#%% PROBLEM DEFINITION

def set_mvols(mvol = [], scale=50, ptype = 'CSH'):
    #TODO similar scaling for CH and CSH
    ch_type = (ptype=='CH')
    csh_type = (ptype=='CSH')
    
    if ch_type:
        mv  = (not bool(mvol)) #mvol is empty
        if (mv or len(mvol)!=2):
            mvolCH = 33.10e-3
            mvolCC = 36.90e-3
            mvol = [mvolCH, mvolCC]
        return mvol
        #return [merged['CH'], merged['CC']]
    elif csh_type:
        for k in mvol:
            mvol[k] = mvol[k]*scale
        default = {'CH': 33.10e-3*scale, 'CC': 36.90e-3*scale,
                   'CSH_TobH': 55.30e-3*scale, 'CSH_TobD': 47.95e-3*scale,
                   'CSH_JenH': 75.63e-3*scale, 'CSH_JenD': 80.58e-3*scale}
        merged = default.copy()   
        merged.update(mvol)
        return [merged['CH'], merged['CC'],
                merged['CSH_TobH'], merged['CSH_TobD'],
                merged['CSH_JenH'], merged['CSH_JenD']]
    else:
        print('Define ptype as \'CH\' or \'CSH\'. ')    
        return []
		
def get_max_pqty(mvol):
    max_pqty = map(lambda x: 1./x, mvol)
    return max_pqty


def set_init_pqty(mvol, scale=50, porosCH = 0.1, wc = 0.45): 
    maxCH = get_max_pqty(mvol)[0]
    initCH = (1 - porosCH) * maxCH
    initCC = 0.0
    init_conc = [initCH, initCC]
    if (len(mvol) == 6):
        initCSH_TobH = 0.1041/scale#(1-init_porosCSH) *maxCSH_TobH
        initCSH_TobD = 2.5050/scale#(1-init_porosCSH) *maxCSH_TobD
        initCSH_JenH = 2.1555/scale#(1-init_porosCSH) *maxCSH_JenH
        initCSH_JenD = 3.2623/scale#(1-init_porosCSH) *maxCSH_JenD
        # set different init if wc ~= 0.45
        init_conc = [initCH, initCC, 
                     initCSH_TobH, initCSH_TobD, 
                     initCSH_JenH, initCSH_JenD]
    return init_conc
    
def get_pqty(init_pq, domain):
    pqty = []
    if (len(init_pq) == 2): #CH type
        pqty_CH = init_pq[0] * (domain.nodetype == ct.Type.MULTILEVEL) + \
            init_pq[0] * (domain.nodetype == ct.Type.MULTILEVEL_CH)
        pqty_CC = init_pq[1] * np.ones(domain.nodetype.shape)
        pqty = [pqty_CH, pqty_CC]
    elif (len(init_pq) == 6):  #CSH type
        pqty_CH = init_pq[0] * (domain.nodetype == ct.Type.MULTILEVEL_CH)    
        pqty_CC = init_pq[1] * np.ones(domain.nodetype.shape)
        pqty_CSHQ_TobH = init_pq[2] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty_CSHQ_TobD = init_pq[3] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty_CSHQ_JenH = init_pq[4] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty_CSHQ_JenD = init_pq[5] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty = [pqty_CH, pqty_CC, 
                pqty_CSHQ_TobH, pqty_CSHQ_TobD, 
                pqty_CSHQ_JenH, pqty_CSHQ_JenD]    
    return pqty


def set_labels(domain, ptype = 'CSH'):
    ch_type = (ptype=='CH')
    csh_type = (ptype=='CSH')
    slabels = np.zeros(domain.nodetype.shape)
    if ch_type:
        slabels = 100002 * (domain.nodetype == ct.Type.LIQUID) + \
          100002 * (domain.nodetype == ct.Type.INTERFACE) + \
          100003 * (domain.nodetype == ct.Type.MULTILEVEL)  + \
          100005 * (domain.nodetype == ct.Type.SOLID)
    elif csh_type:
        slabels = 100002 * (domain.nodetype == ct.Type.LIQUID) + \
          100002 * (domain.nodetype == ct.Type.INTERFACE) + \
          100003 * (domain.nodetype == ct.Type.MULTILEVEL_CH) + \
          100004 * (domain.nodetype == ct.Type.MULTILEVEL)  + \
          100005 * (domain.nodetype == ct.Type.SOLID)
    return slabels

def get_porosity(domain, pqty, mvol, ptype = 'CSH'):
    ch_type = (ptype=='CH')
    csh_type = (ptype=='CSH')
    poros = np.zeros(domain.nodetype.shape)
    if ch_type:
        poros = (1- pqty[0]*mvol[0])*(domain.nodetype == ct.Type.MULTILEVEL) + \
            1*(domain.nodetype != ct.Type.MULTILEVEL)
    elif csh_type:
        poros = (1- pqty[0]*mvol[0])*(domain.nodetype == ct.Type.MULTILEVEL_CH) + \
            (1- pqty[2]*mvol[2] - pqty[3]*mvol[3] - pqty[4]*mvol[4] - \
            pqty[5]*mvol[5])*(domain.nodetype == ct.Type.MULTILEVEL) +\
            1*(domain.nodetype == ct.Type.LIQUID )+\
            1*(domain.nodetype == ct.Type.INTERFACE )+\
            1*(domain.nodetype == ct.Type.SOLID )
    return poros

def set_domain_params(D, mvol, pqty, poros, app_tort, slabels, input_file = 'CH_CC.phrq' ):
    dp={}
    dp['D0'] = D                     
    dp['voxel_vol']=1
    dp['mvol']= mvol
    
    dp['poros'] = poros                 
    dp['app_tort'] = app_tort 
    
    dp['solution_labels']=slabels
    dp['database']='cemdata07.dat'
    
    if(len(mvol) == 1):
        dp['phrqc_input_file']=input_file#'CH_CC_Ceq.phrq'
        dp['eq_names'] = ['portlandite']
        dp['solid_phases']={'portlandite':{'type':'diffusive','mvol':mvol[0],'c':pqty[0]}}
    
    if(len(mvol) == 2):
        dp['phrqc_input_file']=input_file#'CH_CC_Ceq.phrq'
        dp['eq_names'] = ['portlandite', 'calcite']
        dp['solid_phases']={'portlandite':{'type':'diffusive','mvol':mvol[0],'c':pqty[0]},
                            'calcite':    {'type':'diffusive','mvol':mvol[1],'c':pqty[1]}}
    if(len(mvol) == 6):
        dp['phrqc_input_file']=input_file#'CH_CC_Ceq.phrq'
        dp['eq_names'] = ['portlandite', 'calcite']
        dp['ss_names']={'Tob_jen_ss':['CSHQ_TobH','CSHQ_TobD','CSHQ_JenH','CSHQ_JenD']}
        dp['solid_phases']={'portlandite':{'c':pqty[0],'mvol':mvol[0] ,'type':'diffusive'},
                            'calcite':    {'c':pqty[1],'mvol':mvol[1],'type':'diffusive'},
                            'CSHQ_TobH':  {'c':pqty[2],'mvol':mvol[2],'type':'diffusive'},
                            'CSHQ_TobD':  {'c':pqty[3],'mvol':mvol[3],'type':'diffusive'},
                            'CSHQ_JenH':  {'c':pqty[4],'mvol':mvol[4],'type':'diffusive'},
                            'CSHQ_JenD':  {'c':pqty[5],'mvol':mvol[5],'type':'diffusive'},}    
    return dp


def set_solver_params(tfact = None, smart_thres = 1e-8, cphi_fact = 1./3., cphi = 0):
    sp={}
    sp['collision_model']= 'trt' #'diff_vel' #
    sp['magic_para']=1.0/4.0
    sp['cphi_fact']= cphi_fact
    if cphi>0:        
        sp['cphi']= cphi
    
    sp['phrqc_flags'] = {}
    sp['phrqc_flags']['smart_run']=True
    sp['phrqc_smart_run_tol']=smart_thres
    if(tfact):
        sp['tfactbased']=1
        sp['tfact']= tfact
    return sp

def set_bc_params(bc_slabels):
    bcp = {'solution_labels':bc_slabels, 
           'top':['flux', 0.0],
           'bottom':['flux', 0.0],
           'left':['flux', 0.0],
           'right':['flux', 0.0],}
    return bcp

def set_phrqc_input(p, ptype ='CH'):
	# TODO add initial Si concentrations
    # TODO CSH parameters
    '''
    Example:
    p = {'c_bc':{'type':'conc', 'value': 0.01}, 
         'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
         'c_liq':{'type':'eq', 'value': 'calcite'},
         'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
         'ca_liq':{'type':'eq', 'value': 'calcite'} }
    '''
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
    phrqc_input += set_phrqc_bc(p['c_bc'], ptype)
    phrqc_input += set_phrqc_liquid(p['c_liq'], p['ca_liq'], ptype)
    phrqc_input += set_phrqc_mlvl(p['c_mlvl'], p['ca_mlvl'], ptype)    

    phrqc_input += set_phrqc_solid()
    return phrqc_input
    
def set_phrqc_bc(c, ptype):
    phrqc_input = [] 
    phrqc_input.append('#boundary_solution')    
    phrqc_input.append('SOLUTION\t100001')
    phrqc_input.append('\t-units\tmol/kgw')
    phrqc_input.append('\t-water\t1')
    phrqc_input.append('\tpH\t7\tcharge')
    if(c['type'] == 'conc'):
        phrqc_input.append('\tC\t' + str(c['value']) + '\n')
    elif(c['type'] == 'pco2'):
        phrqc_input.append('\tC\t1\tCO2(g)\t-' + str(c['value']) + '\n')
    phrqc_input.append('EQUILIBRIUM_PHASES\t100001\n')
    return phrqc_input
    
def set_phrqc_liquid(c, ca, ptype):
    phrqc_input = [] 
    phrqc_input.append('#solution_liquid')    
    phrqc_input.append('SOLUTION\t100002')
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
    phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
    phrqc_input.append('portlandite\t0\t0')
    phrqc_input.append('calcite\t0\t0\n')
    return phrqc_input

def set_phrqc_mlvl(c, ca, ptype):
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
    return phrqc_input

def set_phrqc_solid():
    phrqc_input = [] 
    phrqc_input.append('#solution_solid')    
    phrqc_input.append('SOLUTION\t100005')
    phrqc_input.append('\t-water\t1\n')
    return phrqc_input

def save_phrqc_input(phrqc,root_dir, name):
    with open(root_dir +'\\phreeqc_input\\' + name + '.phrq', 'w') as f:
        for item in phrqc:
            f.write("%s\n" % item)
#%% PARAMETERS
   
def get_sum_mineral_volume(rt):
    '''
    Total mineral volume
    '''
    #vol = np.sum(rt.solid.vol)
    not_boundary = (rt.phrqc.boundcells !=1)
    vol = np.mean(rt.solid.vol[not_boundary])
    return(vol)
    
def get_average_poros(rt):
    '''
    Average poroisty in the domain
    '''
    #poros = np.mean(rt.solid.poros)
    not_boundary = (rt.phrqc.boundcells !=1)
    poros = np.mean(rt.solid.poros[not_boundary])
    return(poros)
    
def get_average_D(rt):
    #TODO change
    '''
    Avarage effective diffusivity calculated by Archie's formula
    '''
    De = rt.fluid.H.De
    not_boundary = (rt.phrqc.boundcells !=1)  
    return np.mean(De[not_boundary])

def get_delta_portlandite(rt):
    '''
    Differentiation of Portladite
    '''
    s = np.sum(rt.phrqc.dphases['portlandite'])#*getattr(rt.solid, 'portlandite').mvol)
    return(s)

def get_delta_calcite(rt):
    '''
    Differentiation of Calcite
    '''
    s = np.sum(rt.phrqc.dphases['calcite'])#*getattr(rt.solid, 'calcite').mvol)
    return(s)
    
def get_dissolution(rt):
    '''
    N cells that are dissolving
    '''
    s = np.sum(rt.phrqc.dphases['portlandite']<0)
    return(s)
    
def get_precipitation(rt):
    '''
    N cells that are precipitating
    '''
    s = np.sum(rt.phrqc.dphases['calcite']>0)
    return(s)
    
    
def get_portlandite_cells(rt):
    '''
    N cells containing portlandite
    '''
    s = np.sum(rt.solid.portlandite.c>0)
    return(s)
    
def get_calcite_cells(rt):
    '''
    N cells containing calcite
    '''
    s = np.sum(rt.solid.calcite.c>0)
    return(s)
    
    
def get_average_pH(rt):
    '''
    N cells containing calcite
    '''
    nx = rt.fluid.Ca.nx -1
    ny = rt.fluid.Ca.ny -1
    s = np.mean(rt.phrqc.selected_output()['pH'][1:ny,1:nx])
    return(s)
    
def get_co2_uptake(rt):
    return(np.sum(rt.fluid.C._ss[:,0]))

def get_active(rt):
    '''
    N active cells in phreeqc
    '''
    a = rt.phrqc.nactive
    return(a)
    
def get_dt(rt):
    return(rt.dt)


def get_sum_csh(rt): 
    #TODO CSH mass
    return(np.sum(rt.solid.csh))

def get_Ca_solid(rt):
    ca = np.sum(0.8333333*rt.solid.CSHQ_TobD.c[:,:] + 0.6666667*rt.solid.CSHQ_TobH.c[:,:] + 
        1.3333333*rt.solid.CSHQ_JenH.c[:,:] + 1.5*rt.solid.CSHQ_JenD.c[:,:])
    return ca
def get_Si_solid(rt):
    si = np.sum(0.6666667*rt.solid.CSHQ_TobD.c[:,:] + 1.0*rt.solid.CSHQ_TobH.c[:,:] + 
        1.0*rt.solid.CSHQ_JenH.c[:,:] + 0.6666667*rt.solid.CSHQ_JenD.c[:,:])
    return si
def get_Ca_Si(rt):
    r = get_Ca_solid(rt)/get_Si_solid(rt)
    return r

def get_csh_density(rt):
    m_ca = 56.0774 #g/mol
    m_si = 60.08 #g/mol
    m_h2o = 18.01528 #g/mol
    h2o= np.sum(1.8333333333*rt.solid.CSHQ_TobD.c[:,:] + 1.5*rt.solid.CSHQ_TobH.c[:,:] + 
        2.1666666667*rt.solid.CSHQ_JenH.c[:,:] + 2.5*rt.solid.CSHQ_JenD.c[:,:])
    d = h2o*m_h2o + get_Ca_solid(rt)*m_ca + get_Si_solid(rt)*m_si
    return d
    

#%% LISTS OF PARAMETERS
def init_results(pavg=True, pavg_list=[], points=[], ptype='CSH'):
    results = {}     
    params = [] # always save these parameters
    if ptype =='CH':
        params += ['portlandite', 'calcite',
                  'Ca','C','O','H']
    elif ptype =='CSH':
        params += ['CSHQ_TobD', 'CSHQ_JenD', 'CSHQ_JenH', 'CSHQ_TobH', 'portlandite', 'calcite',
                  'Ca','C','O','H','Si']
        params += ['csh','Ca_solid','Si_solid','Ca_Si','csh_density']
    results={name: [] for name in params}       
    results['params'] = params
    if pavg: #average parameters        
        if not pavg_list:           
            pavg_list += ['sum_vol', 'avg_poros', 'avg_D_eff', 
                   'delta_ch', 'delta_cc', 'precipitation','dissolution', 
                   'portlandite_cells', 'calcite_cells', 'active_cells',
                   'dt', 'pH']  #'avg_aperture'
        results.update({name: [] for name in pavg_list})
    results['pavg_list'] = pavg_list
    if points: # points
        pointparamslist = []
        pointparamslist += params
        pointparamslist += ['vol', 'poros', 'pH', 'De', 'vol_CH', 'vol_CC'] 
        if ptype =='CSH': 
            pointparamslist += ['csh', 'vol_CSH']
        l = [[m+ ' '+n for n in [str(n) for n in points]] for m in [str(m) for m in pointparamslist]]
        par_points = [item for sublist in l for item in sublist] #sum(l, [])
        results.update({name: [] for name in par_points}) 
        results['pointparamslist'] = pointparamslist
    results['points'] = points    
    results['time'] = []   
    return results

def append_results(rt, results, step = 1e+2):
    if (rt.iters%step == 0):
        ptype = rt.ptype
        results['time'].append(rt.time)
        for num, phase in enumerate(rt.solid.diffusive_phase_list, start=1):
            results[phase].append(np.sum(rt.solid._diffusive_phaseqty[num-1]))        
        for num, comp in enumerate(rt.fluid.components, start=1):
            results[comp].append(np.sum(getattr(rt.fluid, comp)._c*getattr(rt.fluid, comp).poros))                
        if (ptype == 'CSH'):
            results['csh'].append(get_sum_csh(rt))
            results['Ca_solid'].append(get_Ca_solid(rt))
            results['Si_solid'].append(get_Si_solid(rt))
            results['Ca_Si'].append(get_Ca_Si(rt))
            results['csh_density'].append(get_csh_density(rt))
        # average
        favgall = {'sum_vol':get_sum_mineral_volume,
                'avg_D_eff':get_average_D,
                'avg_poros':get_average_poros,
                'precipitation':get_precipitation,
                'dissolution':get_dissolution,
                'calcite_cells':get_calcite_cells,
                'portlandite_cells':get_portlandite_cells,
                'active_cells':get_active,
                'dt':get_dt,
                'pH':get_average_pH,
                'delta_ch': get_delta_portlandite,
                'delta_cc': get_delta_calcite,
                'co2_uptake': get_co2_uptake}
        favg = {k: favgall[k] for k in results['pavg_list']}
        for key, value in favg.iteritems():
            if key in ['delta_ch', 'delta_cc']:
                if (rt.iters ==0): 
                    results[key].append(0) 
                else:
                    results[key].append(value(rt))
            else:
                results[key].append(value(rt)) 
        # points
        if results['points']: # points
            for p in results['points']:
                results['portlandite'+' ' + str(p)].append(rt.solid.portlandite.c[p])
                results['calcite'+' ' + str(p)].append(rt.solid.calcite.c[p])
                results['Ca'+' ' + str(p)].append(rt.fluid.Ca._c[p]+rt.fluid.Ca._ss[p]/rt.fluid.Ca.poros[p])
                results['C'+' ' + str(p)].append(rt.fluid.C._c[p]+rt.fluid.C._ss[p]/rt.fluid.C.poros[p])
                results['H'+' ' + str(p)].append(rt.fluid.H._c[p]+rt.fluid.H._ss[p]/rt.fluid.H.poros[p])
                results['O'+' ' + str(p)].append(rt.fluid.O._c[p]+rt.fluid.O._ss[p]/rt.fluid.O.poros[p])
                results['poros'+' ' + str(p)].append(rt.solid.poros[p])
                results['vol'+' ' + str(p)].append(rt.solid.vol[p])
                results['De'+' ' + str(p)].append(rt.fluid.H.De[p])
                results['pH'+' ' + str(p)].append(rt.phrqc.selected_output()['pH'][p])
                results['vol_CH'+' ' + str(p)].append(rt.solid.vol_ch[p])
                results['vol_CC'+' ' + str(p)].append(rt.solid.vol_cc[p])
            
            if( ptype=='CSH'):
                results['Si'+' ' + str(p)].append(rt.fluid.Si._c[p]+rt.fluid.Si._ss[p])
                results['CSHQ_TobD'+' ' + str(p)].append(rt.solid.CSHQ_TobD.c[p])
                results['CSHQ_JenD'+' ' + str(p)].append(rt.solid.CSHQ_JenD.c[p])
                results['CSHQ_TobH'+' ' + str(p)].append(rt.solid.CSHQ_TobH.c[p])
                results['CSHQ_JenH'+' ' + str(p)].append(rt.solid.CSHQ_JenH.c[p])
                results['csh'+' ' + str(p)].append(rt.solid.csh[p])
                #results['vol_CSH'+' ' + str(p)].append(rt.solid.vol_csh[p])
                
    return(results)


def filter_results(results, path, name, length = 1e+4):
    l = len(results['time'])
    r = l/length
    filtered_results ={}
    if (r <= 1):
        filtered_results = results
    else:
        filt = int(r)
        filtered_results = {k: v[0::filt] for k, v in results.items()}
    #save_obj(filtered_results, path + str(name) +'_results')
    return filtered_results

#%% SETTINGS
       
def save_settings(settings, bc_params, solver_params, path, name):
    def write_txt(settings, bc_params, solver_params, path, name):
        with open(path + name +'_settings.txt', 'w') as file:
            file.write(json.dumps(settings))
            file.write('\n\n')
            file.write(json.dumps(bc_params))
            file.write('\n\n')
            file.write(json.dumps(solver_params))
        file.close()
    try:
        write_txt(settings, bc_params, solver_params, path, name)
        #os.mkdir(dirName)   
    except IOError:
        os.mkdir(path)
        write_txt(settings, bc_params, solver_params, path, name)
    
#%% PICKLE

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        f.close()

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
def save_pickle(rt, t, path, name):
    save_obj(rt.phases, path + str(name) +'_nodetype_' + str(t))  
    save_obj(rt.solid._poros, path + str(name) +'_porosity_' + str(t)) 

#%% VTI & VTS
    
def save_vti(rt, phases, t, path, name, ptype = 'CSH'):
    
    nx = rt.fluid.Ca.nx -1
    ny = rt.fluid.Ca.ny -1     
    filename = path + str(name) +'_all_' + str(t) 
    
    p = deepcopy(phases[1:ny,1:nx,np.newaxis])
    cCa = rt.fluid.Ca.c[1:ny,1:nx,np.newaxis]
    cH = rt.fluid.H.c[1:ny,1:nx,np.newaxis]
    cO = rt.fluid.O.c[1:ny,1:nx,np.newaxis]
    cC = rt.fluid.C.c[1:ny,1:nx,np.newaxis]
    cCC = rt.solid.calcite.c[1:ny,1:nx,np.newaxis]
    cCH = rt.solid.portlandite.c[1:ny,1:nx,np.newaxis]
    porosity = rt.solid._poros[1:ny,1:nx,np.newaxis]
    
    if (ptype == 'CH'):
        imageToVTK(filename, cellData = {"phases" : p,
                                          "Ca":cCa,
                                          "C":cC,
                                          "O":cO,
                                          "H":cH,
                                          "Calcite":cCC,
                                          "Portlandite":cCH,
                                          "porosity":porosity,
                                          })
    if (ptype == 'CSH'):
        cSi = rt.fluid.Si.c[1:ny,1:nx,np.newaxis]
        cCSHQ_JenD = rt.solid.CSHQ_JenD.c[1:ny,1:nx,np.newaxis]
        cCSHQ_JenH = rt.solid.CSHQ_JenH.c[1:ny,1:nx,np.newaxis]
        cCSHQ_TobD = rt.solid.CSHQ_TobD.c[1:ny,1:nx,np.newaxis]
        cCSHQ_TobH = rt.solid.CSHQ_TobH.c[1:ny,1:nx,np.newaxis]
        imageToVTK(filename, cellData = {"phases" : p,
                                          "Ca":cCa,
                                          "Si":cSi,
                                          "C":cC,
                                          "O":cO,
                                          "H":cH,
                                          "Calcite":cCC,
                                          "Portlandite":cCH,
                                          "CSHQ_JenD":cCSHQ_JenD,
                                          "CSHQ_JenH":cCSHQ_JenH,
                                          "CSHQ_TobD":cCSHQ_TobD,
                                          "CSHQ_TobH":cCSHQ_TobH,
                                          "porosity":porosity,
                                          })
            
def save_vts(domain,phases, t, path, name, ptype = 'CSH'):
    filename = path + str(name) +'_nodetype_' + str(t) #"./fname"
    x,y = domain.meshgrid()
    x = x[:,:,np.newaxis]
    y = y[:,:,np.newaxis]
    z = np.zeros(x.shape)
    p = deepcopy(phases[:,:,np.newaxis])
    pointsToVTK(filename, x, y, z, data = {"phases" : p})
    gridToVTK(filename, x, y, z, pointData = {"phases" : p})
    

#%% PRINT FUNCTIONS
def print_time(st, rt):
    
    print('==========================')
    hours = int(st/3600)
    minutes = int((st/3600 - hours)*60)
    seconds = int(((st/3600 - hours)*60 - minutes)*60)
    print ('time :%s' %(rt.time))
    #print ('real time taken for simulation:%s seconds' %(simulation_time))
    print ('iterations :%s' %(rt.iters))
    print ('real time taken for simulation:%s hours %s min %s sec' %(hours, minutes, seconds))
    print('==========================') 
   
def print_points(rt, points, names=[]):    
    title = get_titles()
    
    fields = {'portlandite': rt.solid.portlandite.c,
              'calcite':     rt.solid.calcite.c,
              'Ca': rt.fluid.Ca.c,
              'C':  rt.fluid.C.c,
              'H':  rt.fluid.H.c,
              'O':  rt.fluid.O.c,
              'phases': rt.solid.phases,
              'poros': rt.solid._poros,
              'vol_CH': rt.solid.vol_ch,
              'vol_CC': rt.solid.vol_cc,
              'target_si':rt.phrqc._target_SI,
              'radius':rt.solid.pore_radius,
              'pore_amount':rt.solid.pore_amount}
    if rt.ptype == 'CSH':
        pass
        #fields.update({'csh':get_csh_conc(rt),'Si':rt.fluid.Si.c}) #'CSHQ_TobD', 'CSHQ_JenD', 'CSHQ_JenH', 'CSHQ_TobH'

    if not names:
        names = ['portlandite', 'calcite', 'Ca', 'C', 'H', 'phases', 'poros',
                 'vol_CH', 'vol_CC', 'target_si', 'radius', 'pore_amount']
        if rt.ptype == 'CSH':
            names +=['csh', 'Si'] #'CSHQ_TobD', 'CSHQ_JenD', 'CSHQ_JenH', 'CSHQ_TobH'
    
    for k in names: 
        print("%s: %s" % (title[k], fields[k][[i[0] for i in points], [i[1] for i in points]]))
              
        
#%% PLOT FUNCTIONS
def get_titles():
    title ={}
    title = {'C':   'Carbon',
             'Ca':  'Calcium',
             'Si':  'Silicium',
             'O':   'Oxygen',
             'H':   'Hydrogen',
             'De':  'Effective diffusivity',
             'CSHQ_JenH':'Jennite (H)',
             'CSHQ_JenD':'Jennite (D)',
             'CSHQ_TobH':'Tobermorite (H)',
             'CSHQ_TobD':'Tobermorite (D)',
             'portlandite':'Portlandite',
             'calcite': 'Calcite',
             'csh':'CSH (*1e-12) [mol]',
             'active_cells':'Number of active nodes (PHREEQC)',
             'avg_aperture':'Average aperture',
             'avg_poros':'Average porosity',
             'avg_D_eff':'Average effective diffusivity',
             'delta_cc':'Calcite differentiation',
             'delta_ch':'Portlandite differentiation',
             'dissolution':'Number of dissolution nodes',
             'precipitation':'Number of precipitation nodes',
             'calcite_cells':'Number of Calcite nodes',
             'portlandite_cells':'Number of Portlandite nodes',
             'sum_vol':'Total mineral volume',
             'vol':'Total mineral volume',
             'dt':'Time step ',
             'pH': 'Average pH ',
             'poros': 'Porosity',
             'vol_CH': 'Portlandite volume',
             'vol_CC': 'Calcite volume',
             'vol_CSH': 'CSH volume',
             'phases': 'Cement phases',
             'target_si': 'Target saturation index',
             'radius': 'Pore radius',
             'pore_amount': 'Amount of pores',             
             'Ca_solid':'solid Ca',
             'Si_solid' :'solid Si',
             'Ca_Si': 'Ca/Si',
             'csh_density': 'Density(CSH)'
             }
    return(title)
     
def get_ylabs():
    ylab = {}
    ylab = {'C':'C (*1e-12) [mol]',
            'Ca':'Ca (*1e-12) [mol]',
            'Si':'Si (*1e-12) [mol]',
            'O':'O (*1e-12) [mol]',
            'H':'H (*1e-12) [mol]',
            'De':  'D_eff [m^2/s]',
            'CSHQ_JenH':'CSH JenH (*1e-12) [mol]',
            'CSHQ_JenD':'CSH JenD (*1e-12) [mol]',
            'CSHQ_TobH':' CSH TobH (*1e-12) [mol]',
            'CSHQ_TobD':'CSH TobD (*1e-12) [mol]',
            'portlandite': 'CH (*1e-12) [mol]',
            'calcite': 'CC (*1e-12) [mol]',
            'csh':'CSH (*1e-12) [mol]',
            'active_cells':'Nodes [-]',
            'avg_aperture':'Aperture [um]',
            'avg_poros':'Porosity [-]',
            'avg_D_eff':'D [m^2/s]',
            'delta_cc':'dCC [mol]',
            'delta_ch':'dCH [mol]',
            'dissolution':'Nodes [-]',
            'precipitation':'Nodes [-]',
            'calcite_cells':'Nodes [-]',
            'portlandite_cells':'Nodes [-]',
            'dt':'dt [s]',
            'pH': 'pH [-]',
            'sum_vol':'Volume [um^3]',
            'vol':'Volume [um^3]',
            'poros': 'Porosity [-]',
            'vol_CH': 'Portlandite volume [um^3]',
            'vol_CC': 'Calcite volume [um^3]',
            'vol_CSH': 'CSH volume [um^3]',
            'phases': 'Phases [-]',
            'target_si': 'Target SI [-]',
            'radius': 'Radius [m]',
            'pore_amount': 'Pores',
            'Ca_solid':'solid Ca [mol]',
            'Si_solid' :'solid Si [mol]',
            'Ca_Si': 'Ca/Si [-]',
            'csh_density': 'Density(CSH) [g/l]'
            }        
    return(ylab)
    
def plot_species(results, names=[], fsize=(8,4)):
    ylab = get_ylabs()
    title = get_titles()
    if not names:
        names = results['params']
    for k in names:        
        plt.figure(figsize=fsize)
        plt.plot(results['time'], results[k])
        plt.legend()
        plt.title(title[k])
        plt.ylabel(ylab[k])
        plt.xlabel('Time [s]')
        plt.show()
        
def plot_avg(results, names=[], fsize=(8,4)):
    ylab = get_ylabs()
    title = get_titles()    
    if not names:
        names = results['pavg_list']
    for k in names:        
        plt.figure(figsize=fsize)
        plt.plot(results['time'], results[k])
        plt.legend()
        plt.title(title[k])
        plt.ylabel(ylab[k])
        plt.xlabel('Time [s]')
        plt.show()
        
def plot_points(results, names=[], fsize=(8,4)):
    ylab = get_ylabs()
    title = get_titles()
    if not results['points']:
        pass    
    else:
        if not names:
            names = results['pointparamslist']
        for k in names:
            plt.figure(figsize=fsize)        
            for p in results['points']:  
                plt.plot(results['time'], results[k+' '+str(p)], label = str(p))
            plt.legend()
            plt.title(title[k] + ' in points ' + str(results['points']))
            plt.ylabel(ylab[k])
            plt.xlabel('Time [s]')
            plt.show()
   
def plot_fields(rt, names={}, fsize = (8,4)):
    #nl = {name: limit}
    
    def plot(field, title, ylab, limit=[], size=fsize):
        plt.figure() 
        if not limit:
            plt.imshow(field)
        else:
            plt.imshow(field, vmin = limit[0], vmax = limit[1])
        
        clb = plt.colorbar()
        clb.ax.get_yaxis().labelpad = 15
        clb.ax.set_ylabel(ylab, rotation=270)
        
        plt.title(title)
        plt.show()
    
    fields = {'portlandite': rt.solid.portlandite.c,
              'calcite':     rt.solid.calcite.c,
              'Ca': rt.fluid.Ca.c,
              'C':  rt.fluid.C.c,
              'H':  rt.fluid.H.c,
              'O':  rt.fluid.O.c,
              'phases': rt.solid.phases,
              'poros': rt.solid.poros}
    if rt.ptype == 'CSH':
        fields.update({'csh':rt.solid.csh,
                       'Si':rt.fluid.Si.c}) #'CSHQ_TobD', 'CSHQ_JenD', 'CSHQ_JenH', 'CSHQ_TobH'
        
        
    if not names:
        names = ['portlandite', 'calcite', 'Ca', 'C', 'H', 'phases', 'poros']
        if rt.ptype == 'CSH':
            names +=['csh', 'Si'] #'CSHQ_TobD', 'CSHQ_JenD', 'CSHQ_JenH', 'CSHQ_TobH'
    
    ylab = get_ylabs()
    title = get_titles()
    for k in names: 
        plot(fields[k], title[k], ylab[k])
    
    

#%% SAVE FIGURES
def  make_output_dir(path):
    try:
        os.mkdir(path)
    except WindowsError:
        pass
    
def save_figures_minerals(rt, max_pqty, t, path, name, ptype = 'CSH', fsize = (12,8)): #time_points[j]
    

    nx = rt.fluid.Ca.nx -1
    ny = rt.fluid.Ca.ny -1
    
    f = plt.figure(figsize = fsize)	
    plt.imshow(rt.solid.portlandite.c[1:ny,1:nx], vmin = 0, vmax = max_pqty[0]) 
    plt.title('Ca(OH)2' + ' time='+ str(t) )
    plt.colorbar()        
    fname = path + name +'_CH_' + str(t) + '.png'
    plt.savefig(fname)
    plt.close(f)  
    
    f = plt.figure(figsize=fsize)
    plt.imshow(rt.solid.calcite.c[1:ny,1:nx], vmin = 0, vmax = max_pqty[1]) 
    plt.title('CaCO3' + ' time='+ str(t) )
    plt.colorbar()        
    fname = path + name +'_CC_' + str(t) + '.png'
    plt.savefig(fname)
    plt.close(f)   
    
    f = plt.figure(figsize=fsize)
    plt.imshow(rt.solid.poros[1:ny,1:nx], vmin = 0, vmax = 1) 
    plt.title('Porosity' + ' time='+ str(t) )
    plt.colorbar()        
    fname = path + name +'_porosity_' + str(t) + '.png'
    plt.savefig(fname)
    plt.close(f)  
    
    f = plt.figure(figsize=fsize)
    plt.imshow(rt.solid.phases[1:ny,1:nx], vmin =-15, vmax = 2) 
    plt.title('Phases' + ' time='+ str(t) )
    plt.colorbar()        
    fname = path + name +'_phases_' + str(t) + '.png'
    plt.savefig(fname)
    plt.close(f) 
    
    
    if (ptype == 'CSH'):
        #csh = get_csh_conc(rt)
        csh = get_sum_csh(rt)
        m = np.sum(max_pqty[2:6])
        
        f = plt.figure(figsize=fsize)
        plt.imshow(csh[1:ny,1:nx], vmin =0, vmax = m) 
        plt.title('CSH' + ' time='+ str(t) )
        plt.colorbar()        
        fname = path + name +'_CSH_' + str(t) + '.png'
        plt.savefig(fname)
        plt.close(f)

def save_figures_mols(rt, t, path, name, ptype = 'CSH', fsize = (12,8)): #time_points[j]
    
    nx = rt.fluid.Ca.nx -1
    ny = rt.fluid.Ca.ny -1
    
    f = plt.figure(figsize=fsize)
    plt.imshow(rt.fluid.Ca.c[1:ny,1:nx]) 
    plt.title('Ca' + ' time='+ str(t) )
    plt.colorbar()        
    fname = path + name +'_Ca_' + str(t) + '.png'
    plt.savefig(fname)
    plt.close(f)   
    
    f = plt.figure(figsize=fsize)
    plt.imshow(rt.fluid.C.c[1:ny,1:nx]) 
    plt.title('CO2' + ' time='+ str(t) )
    plt.colorbar()        
    fname = path + name +'_C_' + str(t) + '.png'
    plt.savefig(fname)
    plt.close(f)  
    
    if (ptype == 'CSH'):
        
        f = plt.figure(figsize=fsize)
        plt.imshow(rt.fluid.Si.c[1:ny,1:nx]) 
        plt.title('Si' + ' time='+ str(t) )
        plt.colorbar()        
        fname = path + name +'_Si_' + str(t) + '.png'
        plt.savefig(fname)
        plt.close(f) 
        
#%% DEPRECATED
        
def get_params_list(ptype = 'CSH'):
    '''
    Use function init_results()
    '''
    params = list()
    def next_params(par, x):
        return (par+ [x] )
    diffusive_phase_list = ['portlandite', 'calcite']
    components = ['Ca','C','O','H']
    if (ptype == 'CSH'):
        diffusive_phase_list = ['CSHQ_TobD', 'CSHQ_JenD', 'CSHQ_JenH', 'CSHQ_TobH', 'portlandite', 'calcite']
        components = ['Ca','C','O','H', 'Si']
    
    params =reduce(next_params, diffusive_phase_list, ['time'])
    params =reduce(next_params, components, params)
    #params =reduce(next_params, rt.solid.diffusive_phase_list, ['time'])
    #params =reduce(next_params, rt.phrqc.components, params)
    params += ['avg_aperture']
    params += ['sum_vol']
    params += ['avg_poros']
    params += ['avg_D_eff']
    params += ['delta_ch']
    params += ['delta_cc']
    params += ['precipitation']
    params += ['dissolution']
    params += ['portlandite_cells']
    params += ['calcite_cells']
    params += ['active_cells']
    params += ['dt']
    params += ['pH']
    if (ptype == 'CSH'):
        params += ['csh']
    return(params)

def append_params(results, rt, D):
    '''
    Use function append_results()
    '''
    #dx = rt.fluid.H.dx
    ptype = rt.ptype
    results['time'].append(rt.time)
    for num, phase in enumerate(rt.solid.diffusive_phase_list, start=1):
        results[phase].append(np.sum(rt.solid._diffusive_phaseqty[num-1]))#* \
                    #8getattr(rt.solid, phase).mvol))  
        
    for num, comp in enumerate(rt.fluid.components, start=1):
        results[comp].append(np.sum(getattr(rt.fluid, comp)._c*getattr(rt.fluid, comp).poros))
        
    results['avg_aperture'].append(get_average_aperture(rt))
    results['sum_vol'].append(get_sum_mineral_volume(rt))#*dx*dx
    results['avg_poros'].append(get_average_poros(rt))
    results['avg_D_eff'].append(get_average_Archie_D(rt))
    results['precipitation'].append(get_precipitation(rt))
    results['dissolution'].append(get_dissolution(rt))
    results['calcite_cells'].append(get_calcite_cells(rt))
    results['portlandite_cells'].append(get_portlandite_cells(rt))
    results['active_cells'].append(get_active(rt))
    results['dt'].append(get_dt(rt))
    results['pH'].append(get_average_pH(rt))
    
    if (rt.iters >0):        
        results['delta_ch'].append(get_delta_portlandite(rt))
        results['delta_cc'].append(get_delta_calcite(rt))
    else:
        results['delta_ch'].append(0)
        results['delta_cc'].append(0)
        
    if (ptype == 'CSH'):
        results['csh'].append(get_sum_csh(rt))#(get_csh_conc(rt))
        
    return(results)
    