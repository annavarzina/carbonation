# -*- coding: utf-8 -*-
'''
Example with precipitation everywhere
Fixed PCO2 at the boundary
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
from copy import deepcopy
np.set_printoptions(precision=5, threshold=np.inf)
import yantra
import cell_type as ct # change the path to cell_type file
import misc_func as fn

from yantra._base import Multicomponent 
from yantra.physics.PhrqcReactiveTransport import PhrqcReactiveTransport 
from yantra.physics.PhrqcReactiveTransport import Solid

from yantra.physics._phrqc_wrapper import  Phrqc


#%% Class

class DissolutionRT(PhrqcReactiveTransport):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        '''
        Reactive transport class for carbonation
        A Lattice Boltzmann method based reactive transport model \
        where reactions are computed using geochemical solver PHREEQC
        '''
        self.auto_time_step = solver_params.get('auto_time_step',True)
        #self.phrqc = Phrqc(domain,domain_params,bc_params,solver_params)  
        self.phrqc = CarbonationPhrqc(domain,domain_params,bc_params,solver_params)   
        components = self.phrqc.components
        bc = self.phrqc.boundary_conditions 
        init_c = self.phrqc.initial_conditions         
        for name in components:
            if name not in domain_params:
                domain_params[name]={}
            domain_params[name]['c']=init_c[name]
            
        for name in bc:
            for comp in components:
                if comp not in bc_params:
                    bc_params[comp]={}
                bc_params[comp][name]=['c',bc[name][comp]]                
        self.fluid=Multicomponent(eqn,components,domain,domain_params,
                                  bc_params,solver_params)
        self.solid=Solid(domain,domain_params,solver_params)        
        self.dx = self.fluid.H.dx  #pcs settings      
        self.ptype = 'CSH' if hasattr(self.fluid, 'Si') else 'CH'
        self.set_volume()
        self.set_porosity()        
        self.nodetype = deepcopy(domain.nodetype)
        self.apply_settings(settings)  
        self.set_phrqc()   
        self.solid.phases = self.update_phases()
        self.update_nodetype()
        self.phrqc.mvol = self.solid.mvol
        self.phrqc.phaselist = self.solid.phaselist
    def advance(self):
        #print('\n===========================================\n')
        #print(self.iters) 
        #if(self.settings['diffusivity']['type']=='fixed'): 
        #    self.update_diffusivity()
        #self.fluid.call('_set_relaxation_params') 
        #print(self.solid.poros[1,:]) 
        self.fluid.call('advance')
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):
            self.update_solid_params()
            self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid.poros) 
        c=deepcopy( self.fluid.get_attr('c'))
        #advance phrqc
        ss=self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        phaseqty=self.solid.update(self.phrqc.dphases)
        if len(phaseqty):
            self.phrqc.modify_solid_phases(phaseqty)
        self.fluid.set_attr('ss',ss)
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)
        #self.update_solid_params()  
        '''  
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1])
        self.phrqc.modify_solid_phases(phaseqty) 
        pqty=self.solid.update(self.phrqc.dphases)    
        self.fluid.set_attr('ss',ss)
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)
        self.update_solid_params()        
        self.solid.phases = self.update_phases()
        '''
        '''
        print('Ca %s' %str(np.array(self.fluid.Ca.c[1,:]) ))
        print('Ca +ss %s' %str(np.array(self.fluid.Ca.c[1,:]) +\
                               np.array(ss['Ca'][1,:]/self.fluid.Ca.convfactors['ss'])/ \
                               np.array(self.fluid.Ca.poros[1,:])))
        print('O %s' %str(np.array(self.fluid.O.c[1,:]) ))
        print('O +ss %s' %str(np.array(self.fluid.O.c[1,:]) + \
                              np.array(ss['O'][1,:]/self.fluid.O.convfactors['ss'])/ \
                              np.array(self.fluid.O.poros[1,:])))
        print('H %s' %str(np.array(self.fluid.H.c[1,:]) ))
        print('H +ss %s' %str(np.array(self.fluid.H.c[1,:]) + \
                              np.array(ss['H'][1,:]//self.fluid.H.convfactors['ss'])/ \
                              np.array(self.fluid.H.poros[1,:])))
        print('CH %s' %str(np.array(self.solid.portlandite.c[1,:])))
        print('dCH %s' %str(np.array(self.phrqc.dphases['portlandite'][1,:])))
        print('poros %s' %str(self.solid.poros))
        '''
    #%% SETTINGS
    def set_volume(self):
        self.solid.vol = np.zeros(self.solid.shape)
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            self.solid.vol += val.c * self.solid.mvol[num-1] 
            
    def set_porosity(self):        
        self.solid.poros=1.- self.solid.vol/self.solid.voxel_vol 
        self.solid.app_tort = 1. * self.solid.poros ** (1./3.)
        if np.any(self.solid.poros<=1e-10):
            sys.exit('Negative or zero porosity')
        
    def set_boundary_cells(self):
        nodes = np.zeros(np.shape(self.nodetype))
        nodes[0,:]  = ct.Type.SOLID
        nodes[-1,:] = ct.Type.SOLID
        nodes[:,0]  = ct.Type.SOLID
        nodes[:,-1] = ct.Type.SOLID
        return nodes        
     
    def set_bc(self, settings):
        ptype = self.ptype
        for key in self.phrqc.boundary_solution_labels:  
            self.fluid.Ca.bc[key] = ['flux', 0.0]
            self.fluid.Ca._bc[key+'bc'] = 'flux'            
            if (settings['bc']['type'] == 'pco2'):                
                self.fluid.C.bc[key] = ['flux', 0.0]
                self.fluid.C._bc[key+'bc'] = 'flux'
                self.fluid.H.bc[key] = ['flux', 0.0]
                self.fluid.H._bc[key+'bc'] = 'flux'
                self.fluid.O.bc[key] = ['flux', 0.0]
                self.fluid.O._bc[key+'bc'] = 'flux'
            if(ptype == 'CSH'):
                self.fluid.Si.bc[key] = ['flux', 0.0]
                self.fluid.Si._bc[key+'bc'] = 'flux'
    
    def set_phrqc(self):
        self.phrqc.boundcells = self.set_boundary_cells()
        self.phrqc.init_port = self.solid.portlandite.c>0 
        self.phrqc.nodetype = deepcopy(self.nodetype) 
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty) 
        
    def apply_settings(self, settings):
        self.settings = settings
        self.phrqc.active = settings['active']
        if(settings['active'] == 'all'):
            self.phrqc.phrqc_flags['smart_run'] = False
            self.phrqc.phrqc_flags['only_interface'] = False
            self.phrqc.nodetype = deepcopy(self.nodetype)
        elif(settings['active'] == 'interface'):        
            self.phrqc.phrqc_flags['smart_run'] = False
            self.phrqc.phrqc_flags['only_interface'] = True
        elif(settings['active'] == 'smart'):
            self.phrqc.phrqc_flags['smart_run'] = True
            self.phrqc.phrqc_flags['only_interface'] = False
        else:
            sys.exit()
                
    #%% UPDATES
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.solid.prev_vol = deepcopy(self.solid.vol)
        self.set_volume()
        self.set_porosity()
        
    def update_phases(self, thres = 1.0e-3):
        ch = self.solid.portlandite.c * self.solid.portlandite.mvol
        tot = ch 
        is_cl = (self.nodetype == 1) #clinker
        is_liq = (tot < thres)* ((self.nodetype ==-2)+(self.nodetype ==-1))
        phases = is_cl*2 + is_liq*(-1)
        if(self.ptype=='CH'):
            is_ch =  ~is_liq * ~is_cl 
            phases += is_ch*(-10)
        return phases
    
    def update_nodetype(self):
        '''
        find neighbous of the multilevel cells
        '''
        prev_nodetype = self.nodetype
        is_port = (self.solid.portlandite.c>0)
        is_solid = self.solid.nodetype == ct.Type.SOLID
        is_critical = np.zeros(np.shape(is_port))
        is_interface = np.zeros(np.shape(is_port))
        #calc_c = self.phrqc.solid_phase_conc['calcite']
        if self.ptype == 'CH':
            is_liquid = (~is_port)
            if(self.iters>1):
                is_critical = (self.solid.pore_size <= self.solid.threshold_pore_size)& (~is_port)
                is_liquid =  (~is_critical)&(~is_port)&(~is_solid)#&((prev_nodetype==-1)|(prev_nodetype==-2))
                is_interface = (~is_critical)&(~is_port)&(~is_solid)&(~is_liquid)
            self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
                ct.Type.MULTILEVEL * is_critical + \
                ct.Type.MULTILEVEL * is_port +\
                ct.Type.INTERFACE * is_interface + \
                ct.Type.SOLID * is_solid
        
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)       
        border, interface = self.get_b_i(self.solid.nodetype)
        #self.solid.border = border
        #self.solid.interface = interface
        self.nodetype = self.solid.nodetype
        
    def get_b_i(self, nodetype):
        is_port = (self.solid.portlandite.c>0)
        is_mineral = is_port | (nodetype==ct.Type.SOLID)
        val = -16
        temp = deepcopy(nodetype)
        temp = temp*(~is_mineral) + val*is_mineral
        border = self.get_border(temp, val, is_port)
        interface = self.get_interfaces(border, temp, val)
        return border, interface
        
    def get_border(self, nt, val, is_port):
        rolled = np.roll(nt, -1, axis = 1)/4+np.roll(nt, 1, axis = 1)/4 +\
              np.roll(nt, -1, axis = 0)/4+np.roll(nt, 1, axis = 0)/4
        is_border = (rolled!=val)&is_port
        return(is_border)
        
    def get_interfaces(self, is_border, nt, val):        
        down = is_border & (np.roll(nt, 1, axis = 0) != val)
        up =  is_border & (np.roll(nt, -1, axis = 0) != val)
        left =  is_border & (np.roll(nt, -1, axis = 1) != val)
        right =  is_border & (np.roll(nt, 1, axis = 1) != val)       
        
        interfaces = {'left': left,
                      'right': right,
                      'up': up,
                      'down': down,
                      'sum' : 1*down + 1*up + 1*left + 1*right}
        
        return interfaces
    def update_diffusivity(self):
        D_CH = self.settings['diffusivity']['D_CH']
        Dref = self.fluid.H.Deref 
        is_port = self.solid.portlandite.c >0
        is_liquid = ~is_port
        De = D_CH *is_port+ Dref * is_liquid 
        Dnew_lb = De/self.solid.poros/self.solid.app_tort
        #print(Dnew_lb)
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)

        
    def correct(self):
        c = deepcopy(self.fluid.get_attr('_c'))
        f = deepcopy(self.fluid.get_attr('_f'))
        flag = False
        for comp in self.fluid.components:
            n = c[comp] < 0
            if n.any():
                flag = True            
                f[comp][n] = f[comp][n] - f[comp][n]*(f[comp][n]<0)
        if flag: 
            self.fluid.set_attr('_f',f)
            self.fluid.set_attr('f',f)
            self.fluid.call('compute_macro_var')
            
class CarbonationPhrqc(Phrqc):
    def active_nodes(self,c,nodetype):
        active =np.ones(self.array_shape).flatten()
        #is_solid = (nodetype > 0)
        prev_c = self.component_conc
        smart_inactive = np.zeros(self.array_shape)
        inactive=0
        if self.iters >1 and self.phrqc_flags['smart_run']:
            for name in self.active_components:
                diff =np.abs( c[name]*(c[name]-prev_c[name])/(c[name]+1e-30)**2)
                smart_inactive +=(1*(diff<self.phrqc_smart_run_tol))
        inactive+=1*(smart_inactive.flatten(order='c') >0)# + 1*is_solid.flatten(order='c')
        if self.phrqc_flags['only_interface']:
            inactive += (1*(nodetype!=-2)-1*(nodetype==-2)).flatten()
        if self.phrqc_flags['only_fluid']: 
            inactive += (1*(nodetype>0)-1*(nodetype<=0)).flatten()
        #tot_solid_phase_conc = self.add_dict(self.flatten_dict(self.solid_phase_conc)) #TODO test w/o solid
        #inactive += -1*(tot_solid_phase_conc>0)
        active -= 1*(inactive>0)
        self.nactive = np.sum(active)
        return active

    
    def modify_eq(self, phaseqty):#, c, nodetype):
        '''
        modifies the phaseqty in phreeqc
        updates saturation index (SI) for calcite
        
        Parameters
        ----------
        phaseqty: dict
            dictionary containing an ndarray of new quantities of equilibrium phases
        '''
        #active_nodes = self.active_nodes(c,nodetype)
        phaseqty = self.flatten_dict(phaseqty)
        modifystr = []
        is_boundary = (self.boundcells ==1).flatten(order='C')
        if (self.iters==0):
            for i,cell in enumerate(range(self.startcell,self.stopcell+1,1)):
                #if active_nodes[i]:
                    modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)
                    if(is_boundary[cell-1]):  
                        self.modify_bc(modifystr)
                    else:
                        for key in phaseqty.keys():
                            modifystr.append("\t -component %s" %(key))
                            modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1]))   
                            modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
        else:
            is_liquid = (self.nodetype.flatten(order='C') == -1) 
            is_mineral = (self.init_port>0).flatten(order='C')  
            si = self._target_SI.flatten(order='C')  
            for i,cell in enumerate(range(self.startcell,self.stopcell+1,1)):
                #if active_nodes[i]:
                    modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)                     
                    if(is_boundary[cell-1]):  
                        self.modify_bc(modifystr)
                    else:
                        for key in phaseqty.keys():
                            modifystr.append("\t -component %s" %(key))
                            modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1]))   
                            if key == 'portlandite':   
                                modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
                            if (key == 'calcite'):
                                if (self.pcs == False):
                                    modifystr.append("\t\t%s\t%s" %('-precipitate_only', 1))
                                elif(self.pcs == True and self.precipitation == 'interface' and self.active == 'interface'):
                                    modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))
                                    modifystr.append("\t\t%s\t%s" %('-precipitate_only', 1))
                                elif(self.pcs == True and self.precipitation == 'interface' and self.active != 'interface'):                            
                                    if (is_liquid[cell-1]):                                   
                                        modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
                                    else:
                                        modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))
                                        modifystr.append("\t\t%s\t%s" %('-precipitate_only', 1))
                                elif(self.pcs == True and self.precipitation == 'all'):
                                    modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))
                                    modifystr.append("\t\t%s\t%s" %('-precipitate_only', 1))                
                                elif(self.pcs == True and self.precipitation == 'mineral'): 
                                    if (is_mineral[cell-1]): 
                                        modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))
                                        modifystr.append("\t\t%s\t%s" %('-precipitate_only', 1))                             
                                    else:   
                                        modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
                                else: 
                                    pass  
                                if (self.ppt == True):
                                    if (phaseqty['portlandite'][cell-1]<=0):# and phaseqty['calcite'][cell-1]<=0 ):
                                        modifystr.append("\t -component\tCO2(g)") 
                                        modifystr.append("\t\tsi\t-%.20e" %self.pinput['value'])
        modifystr.append("end") 
        modifystr ='\n'.join(modifystr)
        #print(modifystr)
        '''
        import json
        with open('modify_eq.txt', 'w') as file:
            file.write(json.dumps(modifystr))
            file.close()
        '''    
        self.IPhreeqc.RunString(modifystr)
        
    def modify_bc(self, modifystr):
        modifystr.append("\t -component calcite")
        modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
        modifystr.append("\t -component portlandite")
        modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
     
    def modify_solution(self,c,dt,nodetype):
        """
        inconc/prevconc flatten dictonary of conc subset.
        nodetype/solid_phase_qty:matrix or flattened either is possible
        to be run by current iphreeqc solver.

        input:
        -----
        inconc = dictionary  containing each component as key and ncells \
        values of conc for each key
        current_time=current time in sec
        nodetype= node type matrix

        output:
        -------
        updates currentselecout
        """
        self.dt = dt
        active_nodes = self.active_nodes(c,nodetype)
        c_trans=deepcopy(c)
        moles = self.flatten_dict(self.to_moles(c))
        modify_str=[]
        runcells=[]
        runcell_str=[]
        for i,cell in enumerate(range(self.startcell,self.stopcell+1,1)):
            if active_nodes[i]:
                runcells.append(str(cell))
                modify_str.append('SOLUTION_MODIFY %i' % cell)
                if 'H' in moles:
                    c = moles['H'][i]
                    if c <=0: c= 1e-30
                    modify_str.append('\t%s\t%.20e' % ('total_h', c))
                if 'O' in moles:
                    c = moles['O'][i]
                    if c <=0: c= 1e-30
                    modify_str.append('\t%s\t%.20e' % ('total_o', c))
                modify_str.append('\t-totals')
                for name,val in moles.iteritems():
                    if (name != 'H') and (name!='O'):
                        c = val[i]
                        if c <=0: c= 1e-30
                        modify_str.append('\t\t%s\t%.20e' % (name, c))
        modify_str.append('end')
        modify_str ='\n'.join(modify_str)
        runcell_str.append('RUN_CELLS')
        runcell_str.append('\t-cells %s'%'\n\t\t'.join(runcells))
        runcell_str.append('\t-start_time %s'%self.time)
        runcell_str.append('\t-time_step %s'%self.dt)
        runcell_str.append('end')
        runcell_str='\n'.join(runcell_str)
        #print(modify_str)
        self.IPhreeqc.AccumulateLine(modify_str)
        self.IPhreeqc.AccumulateLine(runcell_str)
        self.IPhreeqc.RunAccumulated()
        '''
        import json
        with open('modity_solution.txt', 'w') as file:
            file.write(json.dumps(modify_str))
            file.write('\n\n')
            file.write(json.dumps(runcell_str))
            file.close()
        '''
#        c_prev= deepcopy(self.component_conc)
        self.phrqc_flags['update_output'] = True
        self.selected_output(merge_with_previous=True,toarray=True,reshape=True)
        self.poros = self.selected_output()['poros']
        c_current = self.component_conc
        ss={}
        #vol = np.zeros(self.selected_output()['poros'].shape)
        #for i, p in enumerate(self.phaselist):            
        #    vol += self.selected_output()[p] * self.mvol[i]
        #self.poros = 1- vol
        for name in self.components:
            ss[name] = (c_current[name]-c_trans[name])/self.dt
            ss[name] = ss[name] *(active_nodes.reshape(self.array_shape)>0)
            ss[name] *= self.poros
        #print(self.poros[1,:])
        self.time+=self.dt
        self.iters+=1
        return ss
#%% PROBLEM DEFINITION
__doc__= """ 
Reference:
    Poros = 0.05
    PCO2 = 3.4
    IE = 0.5
    Archies relation for diffusivity
"""
#problem type
m = 'CH' #or 'CSH'

ll = 2
l = 5 +ll
dx = 1.0e-6
lx = l*dx
ly = 2*dx

domain = yantra.Domain2D(corner=(0, 0), 
                         lengths=(lx, ly), 
                         dx=dx, 
                         grid_type='nodal')
domain.nodetype[:, ll+1:l+ll] = ct.Type.MULTILEVEL

domain.nodetype[0,:] = ct.Type.SOLID
domain.nodetype[-1,:] = ct.Type.SOLID
domain.nodetype[:,-1] = ct.Type.SOLID

plt.figure(figsize=(5,5))
plt.imshow(domain.nodetype) 
plt.show()
#%%  VALUES
nn='01_test'
slabels = np.zeros(domain.nodetype.shape)
slabels =  100001*(domain.nodetype!= -5) + 100002*(domain.nodetype== -5)
pqty = 0.95*(domain.nodetype==-5)
porosity = 0.05*(domain.nodetype==-5) + 1.0*(domain.nodetype!=-5)
#porosity[1,3] = 0.1
#pqty[1,3] = 0.9
D = 1e-9
#domain params

domain_params={}
domain_params['D0']=D
domain_params['database']='cemdata07.dat'
domain_params['phrqc_input_file']='02_mlvl_portlandite.phrq'
domain_params['solution_labels']=slabels
domain_params['eq_names']=['portlandite']
domain_params['solid_phases']={'portlandite':{'type':'diffusive','mvol':1,'c':pqty}}
domain_params['voxel_vol']=1
domain_params['poros']=porosity
#solver parameters
solver_params={}
solver_params['tauref'] = 0.75#0.55
#solver_params['tfactbased'] = True
#solver_params['tfact'] = 1./6./8

solver_params['collision_model']='trt'
solver_params['magic_para']=1.0/4.00
solver_params['cphi_fact']=1.0/3.0
          
bc_params = {'solution_labels':{'left':100003}, 
            'top':['flux', 0.0],
            'bottom':['flux', 0.0],
            'left':['flux', 0.0],
            'right':['flux', 0.0],}
            

settings = {'active': 'all', # 'all'/'smart'/'interface'
            'diffusivity':{'type':'fixed', #'fixed' or 'archie'cc_archie_ch_kin
                           'D_CH': 1e-9},
           'dx': dx 
           }
            
#%% INITIATE THE SOLVER
rt= DissolutionRT('MultilevelAdvectionDiffusion',  domain, 
                          domain_params, bc_params, solver_params,
                          settings) 
rt.phrqc.phrqc_smart_run_tol = 1e-6
rt.Dref = D
#%% PARAMETERS
#plist =  [(1,2), (1,3), (1,4), (1,5), (1,6), (1,7), (1,8), (1,9), (1,10)]
plist =  [(1,n) for n in np.arange(0, l)]
pavglist = ['avg_poros', 'pH', 'avg_D_eff', 'sum_vol', 'precipitation',
            'dissolution', 'portlandite_cells', 'calcite_cells', 'dt'] 
#'delta_ch', 'delta_cc', 'precipitation','dissolution', 'portlandite_cells', 
#'calcite_cells', 'active_cells','dt', 'pH', 'avg_poros',  'avg_D_eff', 'sum_vol'
results = fn.init_results(pavg=True, pavg_list=pavglist, points=plist, ptype=m)


#%% TIME SETTINGS
itr = 0 
j = 0
nitr = 50
Ts = 0.01#000*3600
rt_port = []
rt_time = []
#%% RUN SOLVER
#print(rt.fluid.Ca.tau)
print(rt.fluid.Ca.De)
print(rt.fluid.Ca.tau)
while  rt.time <=Ts:#itr <= nitr:  #
    #print(rt.iters)
    rt.advance()      
    rt_port.append(np.sum(rt.solid.portlandite.c))
    rt_time.append(rt.time)
    #print('Ca %s' %str(np.array(rt.fluid.Ca.c[1,:])))
    #print('Ca ss%s' %str(np.array(rt.fluid.Ca._ss[1,:])))
    #print('CH %s' %str(np.array(rt.solid.portlandite.c[1,:])))
    itr += 1
    #print(rt.fluid.Ca.tau)
    
#%% SIMULATION TIME
#'''
print('Ca %s' %str(np.array(rt.fluid.Ca.c[1,:])))
print('Ca +ss/theta %s' %str(np.array(rt.fluid.Ca.c[1,:]) + np.array(rt.fluid.Ca._ss[1,:])/np.array(rt.fluid.Ca.poros[1,:])))
print('H +ss %s' %str(np.array(rt.fluid.H.c[1,:]) + np.array(rt.fluid.H._ss[1,:])/np.array(rt.fluid.H.poros[1,:])))
print('O +ss %s' %str(np.array(rt.fluid.O.c[1,:]) + np.array(rt.fluid.O._ss[1,:])/np.array(rt.fluid.O.poros[1,:])))
print('CH %s' %str(np.array(rt.solid.portlandite.c[1,:])))
print('dCH %s' %str(np.array(rt.phrqc.dphases['portlandite'][1,:])))
#fn.plot_fields(carb_rt, names={ 'calcite', 'portlandite', 'Ca', 'C'})
#print(rt.phrqc.selected_output())
#'''

#%%
print('time %s sec' %str(rt.time))
plt.figure()
plt.plot(rt_time, rt_port)
plt.show()

print(rt.fluid.Ca.De)
print(rt.fluid.Ca.tau)