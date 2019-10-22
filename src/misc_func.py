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

#TODO move functions to func
def set_mvols(mvol = [], ptype = 'CSH'):
    '''
    Molar volumes:
        1. CH; 2. CC; 3. TobH; 4. TobD; 5. JenH; 6. JenD
    '''
    ch_type = (ptype=='CH')
    csh_type = (ptype=='CSH')
    mvol = mvol
    mv  = (not bool(mvol)) #mvol is empty
    if ch_type:
        if (mv or len(mvol)!=2):
            mvolCH = 33.10e-3
            mvolCC = 36.90e-3
            mvol = [mvolCH, mvolCC]
    elif csh_type:
        if (mv or len(mvol)!=6):
            mvolCH = 33.10e-3
            mvolCC = 36.90e-3
            mvolCSH_TobH = 55.30e-3
            mvolCSH_TobD = 47.95e-3
            mvolCSH_JenH = 75.63e-3
            mvolCSH_JenD = 80.58e-3
            mvol = [mvolCH, mvolCC, 
                    mvolCSH_TobH, mvolCSH_TobD,
                    mvolCSH_JenH, mvolCSH_JenD]
    else:
        print('Define ptype as \'CH\' or \'CSH\'. ')    
    return mvol

def get_max_pqty(mvol):
    max_pqty = map(lambda x: 1./x, mvol)
    return max_pqty


def set_init_pqty(mvol, porosCH = 0.1, wc = 0.45):
    maxCH = get_max_pqty(mvol)[0]
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
    
    if(len(mvol) == 2):
        dp['phrqc_input_file']=input_file#'CH_CC_Ceq.phrq'
        dp['eq_names'] = ['portlandite', 'calcite']
        dp['solid_phases']={'portlandite':{'type':'diffusive','mvol':mvol[0],'c':pqty[0]},
                            'calcite':    {'type':'diffusive','mvol':mvol[1],'c':pqty[1]}}
    if(len(mvol) == 6):
        dp['phrqc_input_file']='CH_CSH_CC.phrq'#'CH_CC_Ceq.phrq'
        dp['eq_names'] = ['portlandite', 'calcite']
        dp['ss_names']={'Tob_jen_ss':['CSHQ_TobH','CSHQ_TobD','CSHQ_JenH','CSHQ_JenD']}
        dp['solid_phases']={'portlandite':{'c':pqty[0],'mvol':mvol[0] ,'type':'diffusive'},
                            'calcite':    {'c':pqty[1],'mvol':mvol[1],'type':'diffusive'},
                            'CSHQ_TobH':  {'c':pqty[2],'mvol':mvol[2],'type':'diffusive'},
                            'CSHQ_TobD':  {'c':pqty[3],'mvol':mvol[3],'type':'diffusive'},
                            'CSHQ_JenH':  {'c':pqty[4],'mvol':mvol[4],'type':'diffusive'},
                            'CSHQ_JenD':  {'c':pqty[5],'mvol':mvol[5],'type':'diffusive'},}    
    return dp


def set_solver_params(tfact = 1./6.*2):
    sp={}
    sp['collision_model']= 'trt' #'diff_vel' #
    sp['magic_para']=1.0/4.0
    sp['cphi_fact']=1.0/3.0
    
    sp['phrqc_flags'] = {}
    sp['phrqc_flags']['smart_run']=True
    sp['phrqc_smart_run_tol']=1e-8
    
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
    '''
    Example:
    p = {'c_bc':{'type':'conc', 'value': 0.01}, 
         'c_mlvl':{'type':'eq', 'value': 'calcite'}, 
         'c_liq':{'type':'eq', 'value': 'calcite'},
         'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
         'ca_liq':{'type':'eq', 'value': 'calcite'} }
    '''
    phrqc_input = [] 
    phrqc_input += set_phrqc_bc(p['c_bc'])
    phrqc_input += set_phrqc_liquid(p['c_liq'], p['ca_liq'])
    phrqc_input += set_phrqc_mlvl(p['c_mlvl'], p['ca_mlvl'])    
    if(ptype=='CSH'):
        pass 
        #TODO CSH cells
    phrqc_input += set_phrqc_solid()
    return phrqc_input
    
def set_phrqc_bc(c ):
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
    
def set_phrqc_liquid(c, ca):
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
        phrqc_input.append('\tC\t' + str(ca['value']))
    elif(ca['type'] == 'eq'):
        phrqc_input.append('\tCa\t1\t' + str(ca['value']))
    else:
        pass        
    phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
    phrqc_input.append('portlandite\t0\t0')
    phrqc_input.append('calcite\t0\t0\n')
    return phrqc_input

def set_phrqc_mlvl(c, ca):
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
        phrqc_input.append('\tC\t' + str(ca['value']))
    elif(ca['type'] == 'eq'):
        phrqc_input.append('\tCa\t1\t' + str(ca['value']))
    else:
        pass        
    phrqc_input.append('EQUILIBRIUM_PHASES\t100003')
    phrqc_input.append('portlandite\t0\t1')
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
    

def get_average_aperture(rt):
    '''
    Aperture = average between top cells and mineral conteining cells
    '''
    phases = rt.solid.phases
    #ptype = rt.ptype
    nx = rt.fluid.Ca.nx -1
    
    mineral_cells = np.where((phases<=-5) + (phases >0))
    widths = list()
    for i in np.arange(1,nx):
        width = 0
        mineral = np.array([])
        if (mineral_cells[0].size != 0):
            # idx cells containing calcite (CC)
            mineral = np.any(mineral_cells[1] ==i)    
        if (not mineral):
            #no CC and no CH
            width = rt.fluid.Ca.ny
        else:
            # only CC
            width = np.min(mineral_cells[0][np.where(mineral_cells[1] ==i)])
        widths.append(width-1)
    return np.mean(widths)

   
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
    
def get_average_Archie_D(rt):
    '''
    Avarage effective diffusivity calculated by Archie's formula
    '''
    D0 = rt.fluid.H.D0
    Dnew = rt.solid.poros * rt.solid.app_tort * D0    
    not_boundary = (rt.phrqc.boundcells !=1)  
    return np.mean(Dnew[not_boundary])

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

def get_active(rt):
    '''
    N active cells in phreeqc
    '''
    a = rt.phrqc.nactive
    return(a)
    
def get_dt(rt):
    return(rt.dt)

def get_sum_csh(rt): 
    return(np.sum(rt.solid.csh))
    

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
        params += ['csh']
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
            results[comp].append(np.sum((getattr(rt.fluid, comp)._c+getattr(rt.fluid, comp)._ss)/getattr(rt.fluid, comp).poros))        
        if (ptype == 'CSH'):
            results['csh'].append(get_sum_csh(rt))
        # average
        favgall = {'avg_aperture': get_average_aperture,
                'sum_vol':get_sum_mineral_volume,
                'avg_D_eff':get_average_Archie_D,
                'avg_poros':get_average_poros,
                'precipitation':get_precipitation,
                'dissolution':get_dissolution,
                'calcite_cells':get_calcite_cells,
                'portlandite_cells':get_portlandite_cells,
                'active_cells':get_active,
                'dt':get_dt,
                'pH':get_average_pH,
                'delta_ch': get_delta_portlandite,
                'delta_cc': get_delta_calcite}
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
                results['vol_CSH'+' ' + str(p)].append(rt.solid.vol_csh[p])
                
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
    print ('time :%s' %(rt.time))
    #print ('real time taken for simulation:%s seconds' %(simulation_time))
    print ('iterations :%s' %(rt.iters))
    print ('real time taken for simulation:%s hours %s min' %(hours,minutes))
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
             'pore_amount': 'Amount of pores'
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
            'pore_amount': 'Pores'
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
            plt.title(title[k] + ' for points ' + str(results['points']))
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
        fields.update({'csh':rt.get_csh_conc(),
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
    