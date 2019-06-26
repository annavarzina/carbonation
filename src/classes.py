# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 16:30:02 2019

@author: avarzina
"""

import numpy as np
from copy import deepcopy
import yantra
from yantra._base import Multicomponent 
from yantra.physics._phrqc_wrapper import  Phrqc
from yantra.physics.PhrqcReactiveTransport import PhrqcReactiveTransport 
from yantra.physics.PhrqcReactiveTransport import Solid

import modules.add_ons.cell_type as ct # change the path to cell_type file
import func as fn

class MyPhrqcReactiveTransport(PhrqcReactiveTransport):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params):
        '''
        A Lattice Boltzmann method based reactive transport model where reactions are computed using \
        geochemical solver phreeqc
        '''
        self.auto_time_step = solver_params.get('auto_time_step',True)
        self.phrqc = MyPhrqc(domain,domain_params,bc_params,solver_params)
        
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
        self.dx = self.fluid.H.dx
        
        self.ptype = 'CSH' if hasattr(self.fluid, 'Si') else 'CH'
        self.nodetype = deepcopy(domain.nodetype)
                
        self.solid._voxel_vol = 1
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        self.solid._vol = np.zeros(self.solid.nodetype.shape)
        
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            self.solid._vol += val.c * self.solid.mvol[num-1] # 
        
        self.solid._poros=1.- self.solid._vol/self.solid._voxel_vol 
        self.solid._app_tort = 1. * self.solid._poros ** (1./3.)
        self.phrqc.boundcells = fn.boundary_cells(self)
        self.phrqc.init_port = self.solid.portlandite.c>0
        self.phrqc.is_calc = self.solid.calcite.c>0
        self.phrqc._target_SI = np.zeros(self.solid.shape)
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty)
            
        self.update_solid_params()
        self.set_bc()
        self.solid.phases = self.update_phases()
        self.update_nodetype()
        
    
    def advance(self):
        self.update_target_SI()
        if(self.settings['velocity']):
            self.update_velocity()
        if(self.settings['diffusivity']['type']=='fixed'): 
            self.update_diffusivity()
        self.fluid.call('advance')
        self.phrqc.is_calc = self.solid.calcite.c>0
        
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):
            poros = self.solid._poros
            app_tort = self.solid._app_tort
            self.fluid.call('update_transport_params',poros,
                            app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(poros)
        
        self.fluid.call('_set_relaxation_params')
        
        c=deepcopy( self.fluid.get_attr('c'))
        #advance phrqc        
        
        phaseqty=self.solid.phaseqty_for_phrqc()		
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1])
        self.phrqc.modify_solid_phases(phaseqty) 
		
        ss=self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        pqty=self.solid.update(self.phrqc.dphases)
            
        self.fluid.set_attr('ss',ss)
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)
        if hasattr(self.phrqc, 'pm'):
            if (self.phrqc.pm == 'interface'):
                self.phrqc.nodetype = deepcopy(self.solid.nodetype)
        self.update_solid_params() # or after if?
        self.solid.phases = self.update_phases()
        self.update_nodetype()

    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.solid._prev_vol = deepcopy(self.solid._vol)
        self.solid._voxel_vol = 1
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        self.solid._vol = np.zeros(self.solid.nodetype.shape)
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            self.solid._vol += val.c * self.solid.mvol[num-1] 
        self.solid._poros= 1.- self.solid._vol/self.solid._voxel_vol 
        #print(self.solid._poros)
        self.solid._app_tort = 1. * self.solid._poros ** (1./3.)
        self.solid._dvol = self.solid._vol-self.solid._prev_vol
    
    def set_bc(self):
        ptype = self.ptype
        for key in self.phrqc.boundary_solution_labels:
            self.fluid.Ca.bc[key] = ['flux', 0.0]
            self.fluid.Ca._bc[key+'bc'] = 'flux'
            if(ptype == 'CSH'):
                self.fluid.Si.bc[key] = ['flux', 0.0]
                self.fluid.Si._bc[key+'bc'] = 'flux'
                
    def update_phases(self, thres = 1.0e-3):
        ptype = self.ptype
        ch = self.solid.portlandite.c * self.solid.portlandite.mvol
        cc = self.solid.calcite.c * self.solid.calcite.mvol
        tot = cc + ch 
        is_cl = (self.nodetype == 1) #clinker
        is_liq = (tot < thres)* ((self.nodetype ==-2)+(self.nodetype ==-1))
        phases = is_cl*2 + is_liq*(-1)
        if(ptype=='CH'):
            is_cc = (cc==np.maximum.reduce([cc,ch])) * ~is_liq * ~is_cl
            is_ch = (ch==np.maximum.reduce([cc,ch])) * ~is_liq * ~is_cl 
            phases += is_ch*(-10) + is_cc*(-15)
        if(ptype=='CSH'):
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
        #calc_c = self.phrqc.solid_phase_conc['calcite']
        if self.ptype == 'CH':
            is_liquid =  (~is_port)
            if(self.iters>1):
                is_critical = (self.solid.pore_radius <= self.solid.threshold_radius) & is_calc
                #is_critical = (self.solid.target_SI >= self.settings['si_params']['threshold_SI']) & is_calc
                is_liquid =  (~is_critical)&(~is_port)&(~is_solid)
            not_critical = np.logical_not(is_critical)
            self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
                ct.Type.INTERFACE * np.logical_and(is_port, not_critical) + \
                ct.Type.MULTILEVEL * is_critical +\
                ct.Type.SOLID * is_solid #ct.Type.INTERFACE * ((is_calc) &(~is_critical)) + 
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

    def update_velocity(self):
        is_solid = self.solid.nodetype >=1
        not_solid = ~is_solid
        #not_solid = np.ones(self.solid.shape)
        '''
        is_ch = self.solid.portlandite.c >= 0.05
        not_liquid = np.logical_or(is_ch,is_solid)
        ch_neighbors = not_liquid.astype(int) + \
            np.roll(not_liquid, shift = 1, axis= 1).astype(int) + \
            np.roll(not_liquid, shift = -1, axis= 1).astype(int) +\
            np.roll(not_liquid, shift = 1, axis= 0).astype(int) +\
            np.roll(not_liquid, shift = -1, axis= 0).astype(int)
            
        #not_solid =  np.logical_and(~is_solid, ~(ch_neighbors==5))
        '''
        neighbors = np.roll(not_solid, shift = 1, axis= 1).astype(int) + \
            np.roll(not_solid, shift = -1, axis= 1).astype(int) +\
            np.roll(not_solid, shift = 1, axis= 0).astype(int) +\
            np.roll(not_solid, shift = -1, axis= 0).astype(int)
        self.fluid.u = self.solid._dvol/(neighbors + 1*(neighbors==0))*not_solid/self.solid._poros
        
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
            if (self.settings['pores']=='cylinder'):
                if (self.iters<=1):
                    self.solid.threshold_radius = fn.get_thres_radius(self) /self.dx
                    self.solid.radius_max = self.settings['si_params']['radius'] /self.dx # in lb units
                    self.solid.pore_length = self.settings['si_params']['L'] /self.dx # in lb units 
                self.solid.vol_ch = fn.get_CH_volume(self)
                self.solid.vol_cc = fn.get_CC_volume(self)          
                self.solid.pore_amount = fn.get_pore_amount(self)  
                self.solid.free_vol = fn.get_free_volume(self)  
                self.solid.tot_free_vol = self.solid._poros#fn.get_pore_volume(self)
                self.solid.pore_radius = fn.get_radius(self)
                self.solid.cc_pore_volume = fn.get_calcite_pore_volume(self)            
                self.solid.taget_SI = fn.get_target_SI(self)
                self.phrqc._target_SI = self.solid.taget_SI
            elif(self.settings['pores']=='blocks'):                
                if (self.iters<1):      
                    self.solid.threshold_crystal = self.settings['si_params']['threshold_crystal'] #default crystal size
                    self.solid.threshold_distance = self.settings['si_params']['threshold_distance']
                    self.solid.threshold_radius = self.solid.threshold_distance/2
                    self.solid.block_size = fn.get_block_size(self.solid.threshold_distance,
                                            self.solid.threshold_crystal)
                self.solid.vol_ch = fn.get_CH_volume(self)   
                self.solid.vol_cc = fn.get_CC_volume(self) 
                self.solid.tot_free_vol = self.solid._poros * self.dx**3#fn.get_pore_volume(self)
                self.solid.ncrystals = fn.get_block_amount(self.solid.tot_free_vol,
                                            self.solid.threshold_distance, 
                                            self.solid.threshold_crystal) 
                
                self.solid.crystal_size = fn.get_crystal_size(self.solid.vol_cc,
                                        N = self.solid.ncrystals,
                                        dx = self.dx)
                self.solid.distance = fn.get_distance(self.solid.crystal_size,
                                                            self.solid.block_size)
                self.solid.pore_radius =  self.solid.distance/2.
                self.solid.taget_SI = fn.get_target_SI2(self)
                self.phrqc._target_SI = self.solid.taget_SI
            else:
                pass
        elif (self.ptype == 'CSH'): 
            pass        
      
    def update_diffusivity(self):
        D_CC = self.settings['diffusivity']['calcite']
        D_CH = self.settings['diffusivity']['portlandite']
        D0 = self.fluid.Ca.D0  
        is_port = self.solid.portlandite.c >0
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
        is_liquid = np.logical_and(~is_port, ~is_calc)
        De = D_CC *is_calc + D_CH *is_port+ D0 * is_liquid 
        Dnew_lb = De/self.solid._poros/self.solid._app_tort
        #print(Dnew_lb)
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)
    
class MyPhrqc(Phrqc):
    def active_nodes(self,c,nodetype):
        active =np.ones(self.array_shape).flatten()
        #is_boundary = np.zeros(self.array_shape)
        is_boundary = (self.boundcells ==1)
        prev_c = self.component_conc
        smart_inactive = np.zeros(self.array_shape)
        inactive=0
        if self.iters >1 and self.phrqc_flags['smart_run']:
            for name in self.active_components:
                diff =np.abs( c[name]*(c[name]-prev_c[name])/(c[name]+1e-30)**2)
                smart_inactive +=(1*(diff<self.phrqc_smart_run_tol))
        
        inactive+=1*(smart_inactive.flatten(order='c') >0) + 1*is_boundary.flatten(order='c')
        
        if self.phrqc_flags['only_interface']:
            inactive += (1*(nodetype!=-2)-1*(nodetype==-2)+1*is_boundary).flatten()
        #if self.phrqc_flags['only_mlvl']:
        #    inactive += (1*(nodetype==-1)-1*(nodetype!=-1)).flatten()
        if self.phrqc_flags['only_fluid']: 
            inactive += (1*(nodetype>0)-1*(nodetype<=0)+1*is_boundary).flatten()
        tot_solid_phase_conc = self.add_dict(self.flatten_dict(self.solid_phase_conc))
        inactive += -1*(tot_solid_phase_conc>0)
        active -= 1*(inactive>0)
        self.nactive = np.sum(active)
        return active

    
    def modify_eq(self, phaseqty):
        '''
        modifies the phaseqty in phreeqc
        updates saturation index (SI) for calcite
        
        Parameters
        ----------
        phaseqty: dict
            dictionary containing an ndarray of new quantities of equilibrium phases
        '''
        #print("!")
        phaseqty = self.flatten_dict(phaseqty)
        modifystr = []
        pm = 'none'
        #is_boundary = (self.boundcells ==1).flatten(order='C') 
        is_liquid = np.ones(np.shape(self.poros)).flatten(order='C') 
        if hasattr(self, 'pm'):      
            pm = self.pm
            is_liquid = (self.nodetype.flatten(order='C') == -1) 
        if hasattr(self, 'pc'):   
            pass
        else:
            si = self._target_SI.flatten(order='C')      
            for cell in range(self.startcell,self.stopcell+1,1):
                modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)
                for key in phaseqty.keys():
                    modifystr.append("\t -component %s" %(key))
                    modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1]))   
                    if key == 'portlandite':   
                        modifystr.append("\t\t%s\t%.20e" %('-dissolve_only', 1))
                    if (key == 'calcite'):
                        modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))
                        modifystr.append("\t\t%s\t%.20e" %('-precipitate_only', 1))
                        
        modifystr.append("end") 
        modifystr ='\n'.join(modifystr)
        self.IPhreeqc.RunString(modifystr)