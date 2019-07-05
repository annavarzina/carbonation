# -*- coding: utf-8 -*-
import numpy as np
from copy import deepcopy
import yantra
from yantra._base import Multicomponent 
from yantra.physics.PhrqcReactiveTransport import PhrqcReactiveTransport 
from yantra.physics.PhrqcReactiveTransport import Solid

import cell_type as ct # change the path to cell_type file
import phrqc 


class CarbonationRT(PhrqcReactiveTransport):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params):
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
        
        self.dx = self.fluid.H.dx        
        self.ptype = 'CSH' if hasattr(self.fluid, 'Si') else 'CH'
        self.set_volume()
        self.set_porosity()
        
        self.nodetype = deepcopy(domain.nodetype)         
        self.phrqc.boundcells = self.set_boundary_cells()
        self.phrqc.init_port = self.solid.portlandite.c>0
        self.phrqc.is_calc = self.solid.calcite.c>0
        self.phrqc._target_SI = np.zeros(self.solid.shape)        
        
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty)       
        
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
            self.fluid.call('update_transport_params',self.solid._poros,
                            self.solid._app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid._poros)        
        self.fluid.call('_set_relaxation_params')  
        
        c=deepcopy( self.fluid.get_attr('c'))
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

    def set_volume(self):
        self.solid._vol = np.zeros(self.solid.shape)
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            self.solid._vol += val.c * self.solid.mvol[num-1] 
            
    def set_porosity(self):        
        self.solid._poros=1.- self.solid._vol/self.solid.voxel_vol 
        self.solid._app_tort = 1. * self.solid._poros ** (1./3.)
        
    def set_boundary_cells(self):
        nodes = np.zeros(np.shape(self.nodetype))
        nodes[0,:]  = ct.Type.SOLID
        nodes[-1,:] = ct.Type.SOLID
        nodes[:,0]  = ct.Type.SOLID
        nodes[:,-1] = ct.Type.SOLID
        return nodes    
    
    
    def set_bc(self):
        ptype = self.ptype
        for key in self.phrqc.boundary_solution_labels:
            self.fluid.Ca.bc[key] = ['flux', 0.0]
            self.fluid.Ca._bc[key+'bc'] = 'flux'
            if(ptype == 'CSH'):
                self.fluid.Si.bc[key] = ['flux', 0.0]
                self.fluid.Si._bc[key+'bc'] = 'flux'
            
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.solid._prev_vol = deepcopy(self.solid._vol)
        self.set_volume()
        self.set_porosity()
        self.solid._dvol = self.solid._vol-self.solid._prev_vol
        
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
                    self.solid.threshold_radius = self.get_thres_radius() /self.dx
                    self.solid.radius_max = self.settings['si_params']['radius'] /self.dx # in lb units
                    self.solid.pore_length = self.settings['si_params']['L'] /self.dx # in lb units 
                self.solid.vol_ch = self.get_CH_volume()
                self.solid.vol_cc = self.get_CC_volume()          
                self.solid.pore_amount = self.get_pore_amount()  
                self.solid.free_vol = self.get_free_volume()  
                self.solid.tot_free_vol = self.solid._poros#self.get_pore_volume()
                self.solid.pore_radius = self.get_radius()
                self.solid.cc_pore_volume = self.get_calcite_pore_volume()            
                self.solid.taget_SI = self.get_target_SI()
                self.phrqc._target_SI = self.solid.taget_SI
            elif(self.settings['pores']=='blocks'):                
                if (self.iters<1):      
                    self.solid.threshold_crystal = self.settings['si_params']['threshold_crystal'] #default crystal size
                    self.solid.threshold_distance = self.settings['si_params']['threshold_distance']
                    self.solid.threshold_radius = self.solid.threshold_distance/2
                    self.solid.block_size = self.get_block_size(self.solid.threshold_distance,
                                            self.solid.threshold_crystal)
                self.solid.vol_ch = self.get_CH_volume()   
                self.solid.vol_cc = self.get_CC_volume() 
                self.solid.tot_free_vol = (self.solid.voxel_vol-self.solid.vol_ch)* self.dx**3#self.solid._poros * self.dx**3#self.get_pore_volume(self)
                self.solid.ncrystals = self.get_block_amount(self.solid.tot_free_vol,
                                            self.solid.threshold_distance, 
                                            self.solid.threshold_crystal) 
                
                self.solid.crystal_size = self.get_crystal_size(self.solid.vol_cc,
                                        N = self.solid.ncrystals,
                                        dx = self.dx)
                self.solid.distance = self.get_distance(self.solid.crystal_size,
                                                            self.solid.block_size)
                self.solid.pore_radius =  self.solid.distance/2.
                self.solid.taget_SI = self.get_target_SI2() #TODO
                self.phrqc._target_SI = self.solid.taget_SI
            else:
                pass
        elif (self.ptype == 'CSH'): 
            pass        
      
    #%% VOLUMES
    def get_CH_volume(self):
        CH_vol = self.solid.portlandite.c * self.solid.portlandite.mvol
        return CH_vol
    
    def get_CSH_volume(self):
        CSH_vol = self.solid.CSHQ_JenD.c * self.solid.CSHQ_JenD.mvol +\
                  self.solid.CSHQ_JenH.c * self.solid.CSHQ_JenH.mvol +\
                  self.solid.CSHQ_TobD.c * self.solid.CSHQ_TobD.mvol +\
                  self.solid.CSHQ_TobH.c * self.solid.CSHQ_TobH.mvol 
        return CSH_vol
    
    def get_CC_volume(self):
        CC_vol = self.solid.calcite.c * self.solid.calcite.mvol
        return CC_vol
    
    def get_pore_volume(self):
        pore_vol = self.solid._poros
        return pore_vol
    
    def get_calcite_pore_volume(self):
        v = self.solid.pore_amount * np.pi * self.solid.pore_radius**2 * \
            self.solid.pore_length #*2
        return v
    
    def get_free_volume(self): #volume free from mineral and inner pores
        v = self.solid._voxel_vol - self.solid.vol_ch - self.solid.vol_cc - \
                    ((self.solid.radius_max**2) *np.pi * self.solid.pore_length* \
                     self.solid.pore_amount)
        return v
    
    
    def get_free_volume2(self): #TODO rename
        v = self.solid._voxel_vol - self.solid.vol_ch - self.solid.vol_cc
        return v
    #%% PORE RADIUS
    
    def get_pore_amount(free_vol, ptype='CH'):
        #TODO
        '''
        Return matrix of pore amounts per each cell
        '''  
        cc = free_vol.solid.vol_cc   
        n = free_vol.settings['si_params']['N'] * cc #* np.ones(np.shape(pore_vol))
        n[n<1] = 1
        return n    
    #TODO static
    def get_cylinder_radius(pore_vol, pore_num, pore_length, dx):
        '''
        Formula for cylinder pore radius
        Pore volume is total free volume or porosity
        '''
        r = (pore_vol *dx**3 / (pore_num * np.pi * pore_length))**(1./2.) 
        return r
        
    def get_radius(self):
        '''
        Return matrix of pore radiuses per each cell
        '''
        pore_vol = self.solid._poros
        pore_num = self.solid.pore_amount
        pore_length = self.solid.pore_length
        rad_max = self.solid.radius_max 
        rad = np.ones(np.shape(pore_vol)) 
        rad = self.get_cylinder_radius(pore_vol, pore_num, pore_length, 1)
        rad[rad>rad_max] = rad_max
        return rad
    
    def get_thres_radius(self):
        si = self.settings['si_params']['threshold_SI']
        omega = 10.**si
        r = (-1) * self.settings['si_params']['mvol'] *\
            np.cos(np.pi*self.settings['si_params']['angle'])*2*\
            self.settings['si_params']['iene'] /\
            self.settings['si_params']['R'] /\
            self.settings['si_params']['T'] / np.log(omega)
        return r
    
    def get_block_size(self, distance = 1e-8, crystal_size = 1e-7):
        return crystal_size+distance
    
    def get_block_amount(self, free_vol, distance = 1e-8, crystal_size = 1e-7): #TODO
        v_cp = (crystal_size+distance)**3
        N = free_vol/v_cp
        return N
    
    def get_crystal_size(self, vol_cc, N, dx): #TODO
        s = (vol_cc/N)**(1./3.)
        return s * dx
    
    def get_distance(self, crystal_size, block_size): #TODO
        return(block_size - crystal_size)
        
    def get_block_porosity(self, distance = 1e-8, crystal_size = 1e-7, dx=1e-6):
        v_cp = ((crystal_size+distance)/dx)**3
        v_c = (crystal_size/dx)**3
        p = 1 - v_c/v_cp
        return p
    #%% SATURATION INDEX
    
    def get_target_SI(self):
        rad = self.solid.pore_radius
        '''
        rad_max= rt.solid.radius_max 
        is_pcs = (rad < rad_max)
        si = np.zeros(np.shape(rad)) #*rt.settings['si_params']['critical_SI']
        si[is_pcs] = get_PCS_SI(rt, rt.solid.pore_radius[is_pcs]* rt.dx)
        '''
        si = self.get_PCS_SI(rad* self.dx)    
        is_port = self.solid.vol_ch >0
        not_critical = rad >= self.solid.treshold_radius #rt.solid.radius_max 
        si[np.logical_and(is_port, not_critical)]= 0
        return si   
    
    
    def get_target_SI2(self):
        rad = self.solid.pore_radius
        si = self.get_PCS_SI( rad)    
        is_port = (self.solid.vol_ch >0)
        not_critical = (rad >= self.solid.threshold_radius) #rt.solid.radius_max 
        si[np.logical_and(is_port, not_critical)]= 0
        return si   
    
    def get_sat_ratio_PCS(self, r, angle=np.pi, mvol=3.69e-5, int_energy=0.094, gas_c=8.314, temp_kelv=298.3):
        '''
        Calculates saturation ratio 
        '''
        omega = np.exp(-1*mvol*np.cos(angle)* 2 * int_energy / (gas_c * temp_kelv * r ))
        return omega
    
    def get_PCS_SI(self, pore_rad):   
        '''
        Calculates saturation index 
        '''
        si_params = self.settings['si_params']
        angle = si_params.get('angle') * np.pi 
        mvol = si_params.get('mvol') #calcite molar volume
        int_energy = si_params.get('iene')
        gas_c = si_params.get('R')
        temp_kelv = si_params.get('T')
        
        omega = self.get_sat_ratio_PCS(pore_rad, angle, mvol, int_energy, gas_c, temp_kelv )
        si = np.log10(omega)
        return si
    
    def get_csh_conc(self):    
        csh = self.solid.CSHQ_TobH.c+ self.solid.CSHQ_TobD.c +\
                self.solid.CSHQ_JenH.c + self.solid.CSHQ_JenD.c
        return csh 

    def set_feq(self):
        def get_feq(c,poros):
            shape = np.shape(c) #(ly, lx)    
            f = np.zeros((shape[0], shape[1], 5))    
            f[:,:,0] = c*poros
            f[:,:,1] = 0
            f[:,:,2] = 0
            f[:,:,3] = 0
            f[:,:,4] = 0
            return f
        c = self.fluid.get_attr('c')
        poros = self.solid._poros
        f= self.fluid.get_attr('f')
        for key, value in f.iteritems():
            f[key] = get_feq(c[key], poros)
        self.fluid.set_attr('f',f)