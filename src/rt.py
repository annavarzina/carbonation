# -*- coding: utf-8 -*-
import sys
import numpy as np
from copy import deepcopy
import yantra
from yantra._base import Multicomponent 
from yantra.physics.PhrqcReactiveTransport import PhrqcReactiveTransport 
from yantra.physics.PhrqcReactiveTransport import Solid

import cell_type as ct # change the path to cell_type file
import defaults as df
import phrqc 

class CarbonationRT(PhrqcReactiveTransport):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        '''
        Reactive transport class for carbonation
        A Lattice Boltzmann method based reactive transport model \
        where reactions are computed using geochemical solver PHREEQC
        '''
        self.auto_time_step = solver_params.get('auto_time_step',True)
        self.phrqc = phrqc.CarbonationPhrqc(domain,domain_params,bc_params,solver_params)        
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
        self.ptype = 'CSH' if hasattr(self.fluid, 'Si') else 'CH'
        self.set_volume()
        self.set_porosity()  
        self.set_app_tort() #TODO degree      
        self.nodetype = deepcopy(domain.nodetype)  
        self.apply_settings(settings)  
        self.set_phrqc()        
        self.set_bc(settings)
        self.solid.phases = self.update_phases()
        self.update_nodetype()
        
    
    def advance(self):
        self.update_border_and_interface(self.nodetype)
        self.update_target_SI()
        self.update_diffusivity() 
        self.fluid.call('advance')         
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):
            self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid.poros)      
        self.phrqc.is_calc = self.solid.calcite.c>0  
        self.fluid.call('_set_relaxation_params')        
        if(self.phrqc.precipitation == 'interface'):
            if(self.phrqc.active != 'interface'):
                self.phrqc.nodetype = deepcopy(self.nodetype)
        
        c=deepcopy( self.fluid.get_attr('c'))
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1])
        self.phrqc.modify_solid_phases(phaseqty) 
        ss = self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        pqty=self.solid.update(self.phrqc.dphases)   
        if(self.settings['dissolution']=='subgrid'):
            ss=self.update_border_solution(c,ss)
        ss = self.update_no_flux(ss)
        self.fluid.set_attr('ss',ss)
        
        self.update_solid_params()  
        self.solid.phases = self.update_phases()    
        if(self.settings['velocity']):
            self.update_velocity()
        self.update_nodetype()
        
    #%% SETTINGS
    def set_volume(self):
        self.solid.vol = np.zeros(self.solid.shape)
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            self.solid.vol += val.c * self.solid.mvol[num-1] 
            
    def set_porosity(self):        
        self.solid.poros=1.- self.solid.vol/self.solid.voxel_vol
        if np.any(self.solid.poros<=1e-10):
            sys.exit('Negative or zero porosity')
        
    def set_app_tort(self, d = 1./3.):  
        self.solid.app_tort = 1. * self.solid.poros ** d
        
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
        self.phrqc.is_calc = self.solid.calcite.c>0
        self.phrqc._target_SI = np.zeros(self.solid.shape)  
        self.phrqc.nodetype = deepcopy(self.nodetype) 
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty) 
        
    def apply_settings(self, settings):
        self.settings = settings
        self.Dref = settings['Dref']
        self.dx =  settings['dx']
        self.phrqc.pinput = settings['bc']        
        self.phrqc.active = settings['active_nodes']
        self.phrqc.precipitation = settings['precipitation']
        self.phrqc.pcs = settings['pcs_mode']['pcs']
        if(settings['active_nodes'] == 'all'):
            self.phrqc.phrqc_flags['smart_run'] = False
            self.phrqc.phrqc_flags['only_interface'] = False
            self.phrqc.nodetype = deepcopy(self.nodetype)
        elif(settings['active_nodes'] == 'interface'):        
            self.phrqc.phrqc_flags['smart_run'] = False
            self.phrqc.phrqc_flags['only_interface'] = True
        elif(settings['active_nodes'] == 'smart'):
            self.phrqc.phrqc_flags['smart_run'] = True
            self.phrqc.phrqc_flags['only_interface'] = False
        else:
            sys.exit()
        if ('pco2' in settings) and (settings['bc']['type']=='pco2'):            
            self.phrqc.ppt = settings['pco2']
        else:         
            self.phrqc.ppt = False
                
    #%% UPDATES
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.solid.prev_vol = deepcopy(self.solid.vol)
        self.set_volume()
        self.set_porosity()
        self.set_app_tort() #TODO degree
        if(self.settings['velocity'] == True):
            self.solid.dvol = self.solid.vol-self.solid.prev_vol
        
    def update_phases(self, thres = 1.0e-3):
        ch = self.solid.portlandite.c * self.solid.portlandite.mvol
        cc = self.solid.calcite.c * self.solid.calcite.mvol
        tot = cc + ch 
        is_cl = (self.nodetype == 1) #clinker
        is_liq = (tot < thres)* ((self.nodetype ==-2)+(self.nodetype ==-1))
        phases = is_cl*2 + is_liq*(-1)
        if(self.ptype=='CH'):
            is_cc = (cc==np.maximum.reduce([cc,ch])) * ~is_liq * ~is_cl
            is_ch = (ch==np.maximum.reduce([cc,ch])) * ~is_liq * ~is_cl 
            phases += is_ch*(-10) + is_cc*(-15)
        if(self.ptype=='CSH'):
            csh = self.solid.CSHQ_TobH.c * self.solid.CSHQ_TobH.mvol + \
                  self.solid.CSHQ_TobD.c * self.solid.CSHQ_TobH.mvol + \
                  self.solid.CSHQ_JenH.c * self.solid.CSHQ_TobH.mvol + \
                  self.solid.CSHQ_JenD.c * self.solid.CSHQ_TobH.mvol              
            tot += csh            
            is_cc = (cc==np.maximum.reduce([cc,ch,csh]))   * ~is_liq * ~is_cl
            is_ch = (ch==np.maximum.reduce([cc,ch,csh]))   * ~is_liq * ~is_cl 
            is_csh = (csh==np.maximum.reduce([cc,ch,csh])) * ~is_liq * ~is_cl
            phases += is_csh*(-5) + is_ch*(-10) + is_cc*(-15)
        return phases
    
    def update_nodetype(self):
        '''
        find neighbous of the multilevel cells
        '''
        prev_nodetype = self.nodetype
        is_port = (self.solid.portlandite.c>0)
        is_calc = (self.solid.calcite.c>0)
        is_solid = self.solid.nodetype == ct.Type.SOLID
        is_critical = np.zeros(np.shape(is_port))
        is_interface = np.zeros(np.shape(is_port))
        #calc_c = self.phrqc.solid_phase_conc['calcite']
        if self.ptype == 'CH':
            is_liquid = (~is_port)
            if(self.iters>1):
                is_critical = (self.solid.pore_size <= self.solid.threshold_pore_size) & is_calc & (~is_port)
                is_liquid =  (~is_critical)&(~is_port)&(~is_solid)&(~is_calc)#&((prev_nodetype==-1)|(prev_nodetype==-2))
                is_interface = (~is_critical)&(~is_port)&(~is_solid)&(~is_liquid)
            self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
                ct.Type.MULTILEVEL * is_critical + \
                ct.Type.MULTILEVEL * is_port +\
                ct.Type.INTERFACE * is_interface + \
                ct.Type.SOLID * is_solid
        if self.ptype == 'CSH':
            is_csh= (self.solid.CSHQ_JenD.c>0) | (self.solid.CSHQ_JenH.c>0) | \
                    (self.solid.CSHQ_TobD.c>0) |(self.solid.CSHQ_TobH.c>0)
        
            is_liquid =  (~is_port) & (~is_csh)
            if(self.iters>1):
                is_critical = (self.solid.target_SI >= self.settings['si_params']['threshold_SI']) & is_calc& (~is_port) & (~is_csh)
                is_liquid =  (~is_critical)&(~is_port) & (~is_csh)
                
            self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
                ct.Type.MULTILEVEL * is_port  + \
                ct.Type.MULTILEVEL * is_csh  + \
                ct.Type.MULTILEVEL * is_critical +\
                ct.Type.SOLID * is_solid
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)                   
        self.solid.prev_calc_c = deepcopy(self.phrqc.solid_phase_conc['calcite'])
        self.nodetype = self.solid.nodetype
        #self.phrqc.nodetype = self.solid.nodetype

    def update_diffusivity(self):
        if(self.settings['diffusivity']['type']=='fixed'): 
            self.update_diffusivity_fixed()
        elif(self.settings['diffusivity']['type']=='mixed'): 
            self.update_diffusivity_mixed()
        elif(self.settings['diffusivity']['type']=='archie'):
            pass
        else:
            print('Incorrect diffusivity type. Archie\'s relationship is used.')
            
    def update_diffusivity_fixed(self):
        Dref = self.Dref
        D_border = self.Dref
        D_CC = None
        D_CH = None
        if('D_border' in self.settings['diffusivity']):
            D_border = self.settings['diffusivity']['D_border']
        if('D_CC' in self.settings['diffusivity']):
            D_CC = self.settings['diffusivity']['D_CC']
        if('D_CH' in self.settings['diffusivity']):
            D_CH = self.settings['diffusivity']['D_CH']
        
        is_border = self.solid.border
        is_port = (self.solid.portlandite.c >0) & (~is_border)
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
        is_calc = np.logical_and(is_calc,~is_border)
        is_liquid = np.logical_and(~is_port, ~is_calc)
        is_liquid = np.logical_and(is_liquid, ~is_border)
        
        Dnew_lb = Dref*np.ones(np.shape(self.nodetype))
        if(D_CC is None and D_CH is None): 
            De = D_border*is_border + Dref*is_liquid 
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
            Dnew_lb = Dref*is_port + Dref*is_calc + Dnew_lb*np.logical_or(~is_port, ~is_calc)
        elif(D_CC is None and D_CH is not None):
            De = D_CH*is_port+ D_border*is_border + Dref*is_liquid 
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
            Dnew_lb = Dref*is_calc + Dnew_lb*(~is_calc)
        elif(D_CC is not None and D_CH is None):            
            De = D_CC*is_calc+ D_border*is_border + Dref*is_liquid 
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
            Dnew_lb = Dref*is_port + Dnew_lb*(~is_port)
        else:
            De = D_CH*is_port + D_CC*is_calc + D_border*is_border + Dref*is_liquid 
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
            
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        #self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)  
        
    def update_diffusivity_mixed(self):
        #TODO check correctness
        Dref = self.Dref
        D_border = self.Dref
        D_CC = None
        D_CH = None
        if('D_border' in self.settings['diffusivity']):
            D_border = self.settings['diffusivity']['D_border']
        if('D_CC' in self.settings['diffusivity']):
            D_CC = self.settings['diffusivity']['D_CC']
        if('D_CH' in self.settings['diffusivity']):
            D_CH = self.settings['diffusivity']['D_CH']
        
        is_border = self.solid.border
        is_port = (self.solid.portlandite.c >0) & (~is_border)
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
        is_calc = np.logical_and(is_calc,~is_border)
        is_liquid = np.logical_and(~is_port, ~is_calc)
        is_liquid = np.logical_and(is_liquid, ~is_border)
        
        mCH = self.solid.portlandite.c * self.solid.portlandite.mvol
        mCC = self.solid.calcite.c * self.solid.calcite.mvol
        Dnew_lb = Dref*np.ones(np.shape(self.nodetype))
        if(D_CC is None and D_CH is None): 
            De = D_border*is_border + Dref*is_liquid 
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
            Dnew_lb = Dref*is_port + Dref*is_calc + Dnew_lb*np.logical_or(~is_port, ~is_calc)
        elif(D_CC is None and D_CH is not None):
            D_mixed = np.nan_to_num(1./((1-mCH)/Dref/self.solid.poros/self.solid.app_tort + mCH/D_CH), Dref)
            De = D_mixed *np.logical_or(~is_liquid, ~is_border)+ D_border*is_border+ Dref * is_liquid
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
        elif(D_CC is not None and D_CH is None):            
            D_mixed = np.nan_to_num(1./((1-mCC)/Dref/self.solid.poros/self.solid.app_tort + mCC/D_CC), Dref)
            De = D_mixed *np.logical_or(~is_liquid, ~is_border)+ D_border*is_border+ Dref * is_liquid
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
        else:
            D_mixed = np.nan_to_num(1./((1-mCC-mCH)/Dref/self.solid.poros/self.solid.app_tort + mCC/D_CC+ mCH/D_CH), Dref)
            De = D_mixed *np.logical_or(~is_liquid, ~is_border)+ D_border*is_border+ Dref * is_liquid
            Dnew_lb = De/self.solid.poros/self.solid.app_tort                    
        
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        #self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)  
        
    def update_velocity(self):
        is_solid = self.solid.nodetype >=1
        not_solid = ~is_solid
        neighbors = np.roll(not_solid, shift = 1, axis= 1).astype(int) + \
            np.roll(not_solid, shift = -1, axis= 1).astype(int) +\
            np.roll(not_solid, shift = 1, axis= 0).astype(int) +\
            np.roll(not_solid, shift = -1, axis= 0).astype(int)
        self.fluid.u = self.solid.dvol/(neighbors + 1*(neighbors==0))*not_solid*self.solid.poros    
        ux = (np.roll(self.fluid.u,1,1) - np.roll(self.fluid.u,-1,1))* \
            not_solid
        uy = (np.roll(self.fluid.u,1,0) - np.roll(self.fluid.u,-1,0))* \
            not_solid
        velocity = np.array((uy,ux))
        velocity={}
        for name in self.fluid.components:
            velocity[name] = np.array((ux,uy))
        self.fluid.set_attr('u',velocity)
    
    def update_target_SI(self):
        if (self.ptype == 'CH'):  
            self.solid.vol_ch = self.volume_CH()
            self.solid.vol_cc = self.volume_CC()          
            if (self.settings['pcs_mode']['pores']=='cylinder'):
                #TODO check this case
                if (self.iters<=1):
                    self.solid.threshold_pore_size = self.settings['pcs_mode']['crystal_size']
                    self.solid.pore_length = self.settings['pcs_mode']['crystal_size']  
                self.solid.pore_amount = self.pore_amount()  
                self.solid.free_vol_cc = (self.solid.voxel_vol-self.solid.vol_ch)* self.dx**3#self.free_volume()  
                self.solid.pore_size = self.pore_size()
                self.solid.pore_volume_cc = self.pore_volume_CC() 
            elif(self.settings['pcs_mode']['pores']=='block'):                
                if (self.iters<1):      
                    self.solid.threshold_crystal = self.settings['pcs_mode']['crystal_size'] #default crystal size
                    self.solid.threshold_pore_size = self.settings['pcs_mode']['pore_size']
                    self.solid.threshold_distance = self.solid.threshold_pore_size*2
                    self.solid.block_size = self.block_size()
                self.solid.free_vol_cc = (self.solid.voxel_vol-self.solid.vol_ch)* self.dx**3#self.solid._poros * self.dx**3#self.get_pore_volume(self)
                self.solid.ncrystals = self.block_amount()                 
                self.solid.crystal_size = self.crystal_size()
                self.solid.distance = self.distance()
                self.solid.pore_size =  self.pore_size()
            else:
                pass        
            self.solid.target_SI = self.target_SI()
            self.phrqc._target_SI = self.solid.target_SI
        elif (self.ptype == 'CSH'): 
            pass     
    
    def update_no_flux(self, ss):        
        ss['Ca'] = ss['Ca']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1) 
        return(ss)
        
    def update_border_and_interface(self, nodetype):
        is_port = (self.solid.portlandite.c>0)
        is_mineral = is_port | (nodetype==ct.Type.SOLID)
        val = -16
        temp = deepcopy(nodetype)
        temp = temp*(~is_mineral) + val*is_mineral
        border = self.get_border(temp, val, is_port)
        interface = self.get_interfaces(border, temp, val)
        self.solid.border = border
        self.solid.interface = interface
    
    def get_border(self, nt, val, is_port):
        rolled = np.roll(nt, -1, axis = 1)/4+np.roll(nt, 1, axis = 1)/4 +\
              np.roll(nt, -1, axis = 0)/4+np.roll(nt, 1, axis = 0)/4
        is_border = (rolled!=val)&is_port
        return(is_border)
        
    def get_interfaces(self, is_border, nt, val):        
        down = is_border & (np.roll(nt, 1, axis = 0) != val)
        up =  is_border & (np.roll(nt, -1, axis = 0) != val)
        left =  is_border & (np.roll(nt, 1, axis = 1) != val)
        right =  is_border & (np.roll(nt, -1, axis = 1) != val)       
        
        interfaces = {'left': left,
                      'right': right,
                      'up': up,
                      'down': down,
                      'sum' : 1*down + 1*up + 1*left + 1*right}        
        return interfaces
    
    def update_border_solution(self,c,ss):
        phrqc_poros = self.phrqc.selected_output()['poros']
        fraction = self.settings['subgrid']['fraction']
        result = {}
        by = np.where(self.solid.border)[0]
        bx = np.where(self.solid.border)[1]
        df = np.where(self.solid.border.flatten())[0]
        lx = self.nodetype.shape[1]
        for i in np.arange(0, np.sum(self.solid.border)):
            if fraction is None:
                fraction = 1-phrqc_poros[by[i], bx[i]]
                #fraction = fraction - phrqc_poros[by[i], bx[i]]
            if (self.solid.interface['down'][by[i], bx[i]]):
                cell_i = df[i]+1-lx
                cell_m = df[i]+1
                result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                      self.solid.portlandite.c[by[i], bx[i]],
                                      self.solid.calcite.c[by[i], bx[i]], 
                                      self.solid.calcite.c[by[i]-1, bx[i]],
                                      fraction)
                ssnew={}
                for name in self.phrqc.components:
                    ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i]-1, bx[i]])/self.dt
                    ssnew[name] *= phrqc_poros[by[i]-1, bx[i]]
                    ss[name][by[i]-1, bx[i]] = ssnew[name]
            if (self.solid.interface['up'][by[i], bx[i]]):
                cell_i = df[i]+1+lx
                cell_m = df[i]+1
                result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                      self.solid.portlandite.c[by[i], bx[i]],
                                      self.solid.calcite.c[by[i], bx[i]], 
                                      self.solid.calcite.c[by[i]+1, bx[i]],
                                      fraction)#self.solid.poros[1,1])
                ssnew={}
                for name in self.phrqc.components:
                    ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i]+1, bx[i]])/self.dt
                    ssnew[name] *= phrqc_poros[by[i]+1, bx[i]]
                    ss[name][by[i]+1, bx[i]] = ssnew[name]
            if (self.solid.interface['left'][by[i], bx[i]]):
                cell_i = df[i]
                cell_m = df[i]+1
                result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                      self.solid.portlandite.c[by[i], bx[i]],
                                      self.solid.calcite.c[by[i], bx[i]], 
                                      self.solid.calcite.c[by[i], bx[i]-1],
                                      fraction)
                ssnew={}
                for name in self.phrqc.components:
                    ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i], bx[i]-1])/self.dt
                    ssnew[name] *= phrqc_poros[by[i], bx[i]-1]
                    ss[name][by[i], bx[i]-1] = ssnew[name]
               
            if (self.solid.interface['right'][by[i], bx[i]]):
                cell_i = df[i]+2
                cell_m = df[i] +1
                result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                      self.solid.portlandite.c[by[i], bx[i]],
                                      self.solid.calcite.c[by[i], bx[i]], 
                                      self.solid.calcite.c[by[i], bx[i]+1],
                                      fraction)
                ssnew={}
                for name in self.phrqc.components:
                    ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i], bx[i]+1])/self.dt
                    ssnew[name] *= phrqc_poros[by[i], bx[i]+1]
                    ss[name][by[i], bx[i]+1] = ssnew[name]
        
        for i in np.arange(0, np.sum(self.solid.border)):
            if (self.solid.interface['down'][by[i], bx[i]]):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]-lx+1)]['portlandite_m']
                #self.solid.calcite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]-lx+1)]['calcite_m']
                #self.solid.calcite.c[by[i]-1, bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]-lx+1)]['calcite_i']
            if (self.solid.interface['up'][by[i], bx[i]]):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]+lx+1)]['portlandite_m']
                #self.solid.calcite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]+lx+1)]['calcite_m']
                #self.solid.calcite.c[by[i]+1, bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]+lx+1)]['calcite_i']
            if (self.solid.interface['left'][by[i], bx[i]]):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i])]['portlandite_m']  
                self.solid.calcite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i])]['calcite_m']  
                self.solid.calcite.c[by[i], bx[i]-1] = result[str(df[i]+1) + ' ' +str(df[i])]['calcite_i']    
            if (self.solid.interface['right'][by[i], bx[i]]):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]+2)]['portlandite_m'] 
                #self.solid.calcite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]+2)]['calcite_m']  
                #self.solid.calcite.c[by[i], bx[i]+1] = result[str(df[i]+1) + ' ' +str(df[i]+2)]['calcite_i']    
        
        return ss
        
    def update_neighbour_solution(self, result, n_int, n_ch, m_ch, m_cc, i_cc, fraction=1):
        ncell = 123456
        modify_str = []
        modify_str.append("EQUILIBRIUM_PHASES %i" % ncell)
        modify_str.append("Portlandite 0 %.20e dissolve only" %(m_ch))   
        modify_str.append("Calcite 0 %.20e dissolve only" %(m_cc))   # m_cc
        modify_str.append("END") 
        modify_str.append('MIX %i' % ncell)         
        modify_str.append('%i %.20e' %( n_int, max(fraction,0.))) #modify_str.append('%i 1' %n_int)  
        modify_str.append('SAVE solution %i' % ncell)  
        modify_str.append("END") 
        modify_str.append('USE solution %i' % ncell)  
        modify_str.append('USE equilibrium_phase %i' % ncell)
        modify_str.append('SAVE solution %i' % ncell)   
        modify_str.append("END") 
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str) 
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        port = output[2][10] 
        calc = output[2][12] 
        modify_str = [] 
        modify_str.append("EQUILIBRIUM_PHASES %i" %n_ch)
        modify_str.append("Portlandite 0 %.20e dissolve only" %(port))
        modify_str.append("Calcite 0 %.20e precipitate only" %(calc)) #
        modify_str.append("END") 
        modify_str.append('USE equilibrium_phase %i' %n_int)      
        modify_str.append('MIX %i' % ncell)      
        modify_str.append('%i 1' % ncell)   
        modify_str.append('%i %.20e' %( n_int, min(1.-fraction,1.))) #modify_str.append('%i 0' %n_int)    
        modify_str.append('SAVE solution %i' %(n_int))  
        modify_str.append('SAVE equilibrium_phase %i' %n_int)
        modify_str.append('SAVE equilibrium_phase %i' %(n_ch))  
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str)  
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        comp = {}
        comp['portlandite_m'] = port
        comp['calcite_m'] = calc
        comp['calcite_i'] = output[1][12]
        comp['C'] = output[1][6]
        comp['Ca'] = output[1][7]
        comp['H'] = (output[1][8] - self.phrqc.H_norm)
        comp['O'] = (output[1][9] - self.phrqc.O_norm)
        result[str(n_ch) + ' ' + str(n_int)] = comp        
        return(result)     
    
    #%% VOLUMES
    def volume_CH(self):
        CH_vol = self.solid.portlandite.c * self.solid.portlandite.mvol
        return CH_vol
    
    def volume_CC(self):
        CC_vol = self.solid.calcite.c * self.solid.calcite.mvol
        return CC_vol
    
    def volume_CSH(self):
        CSH_vol = self.solid.CSHQ_JenD.c * self.solid.CSHQ_JenD.mvol +\
                  self.solid.CSHQ_JenH.c * self.solid.CSHQ_JenH.mvol +\
                  self.solid.CSHQ_TobD.c * self.solid.CSHQ_TobD.mvol +\
                  self.solid.CSHQ_TobH.c * self.solid.CSHQ_TobH.mvol 
        return CSH_vol

    def pore_volume_CC(self):
        v = 0
        if(self.settings['pcs_mode']['pores'] == 'cylinder'):
            v = self.solid.pore_amount * np.pi * self.solid.pore_size**2 * \
                self.solid.pore_length #*2
        else:
            v = self.free_volume()
        return v
       
    def free_volume(self):
        v = self.solid.voxel_vol - self.solid.vol_ch
        return v
    
    #%% PORE PROPERTIES    
    def pore_amount(self, ptype='CH'):
        '''
        Return matrix of pore amounts per each cell
        '''  
        vol = self.free_volume()
        n = self.settings['pcs_mode']['pore_density'] * vol #* np.ones(np.shape(pore_vol))
        n[n<1] = 1
        return n   
    
    @staticmethod
    def get_cylinder_radius(pore_vol, pore_num, pore_length, dx):
        '''
        Formula for cylinder pore radius
        Pore volume is total free volume or porosity
        '''
        r = (pore_vol *dx**3 / (pore_num * np.pi * pore_length))**(1./2.) 
        return r    
    
    def pore_size(self):
        '''
        Return matrix of pore radiuses per each cell
        '''
        pore_size = np.ones(np.shape(self.solid.shape)) 
        if(self.settings['pcs_mode']['pores'] == 'cylinder'):
            pore_size = self.get_cylinder_radius(self.solid.poros, 
                                                 self.solid.pore_amount, #self.settings['pcs']['pore_density']
                                                 self.solid.pore_length, 
                                                 self.dx)
            #pore_size[pore_size>self.solid.threshold_pore_size ] = self.solid.threshold_pore_size
        else:
            pore_size=self.solid.distance/2.
        return pore_size
    
    def threshold_pore_size(self):
        si = self.settings['pcs_mode']['threshold_SI']
        omega = 10.**si
        ps = (-1) * self.settings['si_params']['mvol'] *\
            np.cos(np.pi*self.settings['si_params']['angle'])*2*\
            self.settings['pcs_mode']['iene'] /\
            self.settings['pcs_mode']['R'] /\
            self.settings['pcs_mode']['T'] / np.log(omega)
        return ps
    
    def block_size(self):
        return self.solid.threshold_crystal+self.solid.threshold_distance
    
    def block_amount(self): 
        v_cp = (self.solid.threshold_crystal+self.solid.threshold_distance)**3
        N = self.solid.free_vol_cc/v_cp
        return N
    
    def crystal_size(self):
        return (self.solid.vol_cc/self.solid.ncrystals)**(1./3.) * self.dx 
    
    def distance(self):   
        return(self.solid.block_size - self.solid.crystal_size)
        
    def block_porosity(self):
        v_cp = ((self.solid.crystal_size+self.solid.distance)/self.dx )**3
        v_c = (self.solid.crystal_size+self.dx )**3
        p = 1 - v_c/v_cp
        return p
    
    #%% SATURATION INDEX    
    def target_SI(self):        
        pore_size = self.solid.pore_size
        si = self.saturation_index(pore_size)    
        is_port = (self.solid.vol_ch >0)
        not_critical = (pore_size >= self.solid.threshold_pore_size) #rt.solid.radius_max 
        si[np.logical_and(is_port, not_critical)]= 0
        return si   
    
    @staticmethod
    def saturation_ratio(r, int_energy=0.094, angle=np.pi, mvol=3.69e-5, gas_c=8.314, temp_kelv=298.3):
        '''
        Saturation ratio for pore-size controlled solubility
        '''
        omega = np.exp(-1*mvol*np.cos(angle)* 2 * int_energy / (gas_c * temp_kelv * r ))
        return omega
    
    def saturation_index(self, size):   
        '''
        Saturation index 
        '''
        params = self.settings['pcs_mode']
        omega = self.saturation_ratio(size, 
                                      params.get('int_energy'), 
                                      df.ANGLE, 
                                      df.MVOL_CC, 
                                      df.GAS_CONSTANT, 
                                      df.TEMPERATURE )
        si = np.log10(omega)
        return si
    
    #%% PROPERTIES
    def csh_conc(self):    
        csh = self.solid.CSHQ_TobH.c+ self.solid.CSHQ_TobD.c +\
                self.solid.CSHQ_JenH.c + self.solid.CSHQ_JenD.c
        return csh 
                