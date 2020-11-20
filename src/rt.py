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
        if(eqn == 'MultilevelAdvectionDiffusion'): 
            eqn = 'MultilevelAdvectionDiffusion2'
        self.auto_time_step = solver_params.get('auto_time_step',True)
        self.phrqc = phrqc.CarbonationPhrqc(domain,domain_params,
                                            bc_params,solver_params)        
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
        self.nodetype = deepcopy(domain.nodetype)  
        self.apply_settings(settings)  
        self.update_volume()
        self.update_porosity()  
        self.update_app_tort() 
        self.update_init_phrqc()        
        self.update_bc(settings)
        self.update_phases()
        self.update_nodetype()
        self.update_target_SI() 
        
    
    def advance(self):
        self.update_diffusivity() 
        self.fluid.call('advance') 
        self.update_target_SI() 
        self.update_source()       
        self.update_phases()    
        self.update_velocity()
        self.update_nodetype()
        
    #%% UPDATES
    def update_volume(self):
        vol = np.zeros(self.solid.shape)
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            vol += val.c * self.solid.mvol[num-1] 
        self.solid.vol = vol
            
    def update_porosity(self):        
        self.solid.poros=1.- self.solid.vol/self.solid.voxel_vol
        if np.any(self.solid.poros<=1e-10):
            sys.exit('Negative or zero porosity')
        
    def update_app_tort(self, ):  
        self.solid.app_tort = 1. * self.solid.poros ** self.solid.app_tort_degree
        
    def update_boundary_cells(self):
        nodes = np.zeros(np.shape(self.nodetype))
        nodes[0,:]  = ct.Type.SOLID
        nodes[-1,:] = ct.Type.SOLID
        nodes[:,0]  = ct.Type.SOLID
        nodes[:,-1] = ct.Type.SOLID
        return nodes        
     
    def update_bc(self, settings):
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
    
    def update_init_phrqc(self):
        self.phrqc.boundcells = self.update_boundary_cells()
        if 'portlandite' in self.solid.diffusive_phase_list:
            self.phrqc.init_port = self.solid.portlandite.c>0
        if self.ptype == 'CSH':
            self.phrqc.init_csh = np.logical_or(np.logical_or(self.solid.CSHQ_TobD.c>0, self.solid.CSHQ_TobH.c>0),
                                                np.logical_or(self.solid.CSHQ_JenD.c>0, self.solid.CSHQ_JenH.c>0))
        if 'calcite' in self.solid.diffusive_phase_list:
            self.phrqc.is_calc = self.solid.calcite.c>0
        self.phrqc._target_SI = np.zeros(self.solid.shape)  
        self.phrqc.nodetype = deepcopy(self.nodetype) 
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty) 
                
            
    def update_phrqc(self):
        self.phrqc.poros=deepcopy(self.solid.poros)      
        self.phrqc.is_calc = self.solid.calcite.c>0   #TODO check if necessary - move to update phrqc function
        if(self.phrqc.precipitation == 'interface'):
            if(self.phrqc.active != 'interface'):
                self.phrqc.nodetype = deepcopy(self.nodetype)
        
    def apply_settings(self, settings):
        self.settings = settings
        self.Dref = settings['Dref']
        self.dx =  settings['dx']
        self.solid.app_tort_degree = settings['app_tort']['degree']
        self.phrqc.pinput = settings['bc']        
        self.phrqc.active = settings['active_nodes']
        self.phrqc.precipitation = settings['precipitation']
        self.phrqc.pcs = settings['pcs_mode']['pcs']
        self.phrqc.phrqc_flags['only_mineral'] = False
        self.phrqc.phrqc_flags['only_interface'] = False
        self.phrqc.phrqc_flags['only_fluid'] = False
        self.phrqc.phrqc_flags['smart_run'] = False
        if(settings['active_nodes'] == 'all'):
            self.phrqc.nodetype = deepcopy(self.nodetype)
        elif(settings['active_nodes'] == 'interface'):    
            self.phrqc.phrqc_flags['only_interface'] = True
        elif(settings['active_nodes'] == 'smart'):
            self.phrqc.phrqc_flags['smart_run'] = True
        elif(settings['active_nodes'] == 'mineral'):
            self.phrqc.phrqc_flags['only_mineral'] = True
        elif(settings['active_nodes'] == 'fluid'):
            self.phrqc.phrqc_flags['only_fluid'] = True
        else:
            sys.exit()
        if ('pco2' in settings) and (settings['bc']['type']=='pco2'):            
            self.phrqc.ppt = settings['pco2']
        else:         
            self.phrqc.ppt = False
            
            
    def update_source(self):        
        
        c=deepcopy( self.fluid.get_attr('c'))
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1]) 
        if(self.settings['dissolution']=='subgrid'):
            phaseqty['portlandite'][np.where(self.solid.border)]=0
            #phaseqty['portlandite'] = np.zeros(phaseqty['portlandite'].shape)
        self.phrqc.modify_solid_phases(phaseqty)        
        ss = self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)  
        if (self.settings['dissolution']=='multilevel'): 
            pqty=self.solid.update(self.phrqc.dphases)  
        else:
            if self.iters >1 :pqty=self.solid.update(self.phrqc.dphases)  
            self.update_solid_params() 
            ss=self.update_border_solution(c,ss)
                  
        self.update_solid_params()  
        self.update_target_SI() 
        self.update_phrqc()  
        ss = self.update_no_flux(ss)
        self.fluid.set_attr('ss',ss)  
                
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.solid.prev_vol = deepcopy(self.solid.vol)
        self.update_volume()
        self.update_porosity()
        self.update_app_tort()
        if(self.settings['velocity'] == True):
            self.solid.dvol = self.solid.vol-self.solid.prev_vol
        self.fluid.call('update_transport_params',self.solid.poros, #move to solid params function
                            self.solid.app_tort,self.auto_time_step)
        
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
            self.solid.csh=csh
        self.solid.phases = phases
        #self.c = []
        self.solid.prev_calc = deepcopy(self.solid.calcite.c > 1e-4)
        self.solid.prev_port = deepcopy(self.solid.portlandite.c > 0)
    
    def update_nodetype(self):
        '''
        find neighbous of the multilevel cells
        '''
        prev_nodetype = deepcopy(self.nodetype)
        is_port = (self.solid.portlandite.c>0)
        is_calc = (self.solid.calcite.c>0)
        #is_clinker = self.solid.nodetype == ct.Type.CLINKER
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
                    (self.solid.CSHQ_TobD.c>0) |(self.solid.CSHQ_TobH.c>0) &(~is_port) 
            is_liquid = (~is_port) & (~is_csh)
            if(self.iters>1):
                is_critical = (self.solid.pore_size <= self.solid.threshold_pore_size) & is_calc & (~is_port)
                is_liquid =  (~is_critical)&(~is_port)&(~is_solid)&(~is_calc)&(~is_csh)#&((prev_nodetype==-1)|(prev_nodetype==-2))
                is_interface = (~is_critical)&(~is_port)&(~is_solid)&(~is_liquid)&(~is_csh)
            self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
                ct.Type.MULTILEVEL * is_critical + \
                ct.Type.MULTILEVEL * is_port +\
                ct.Type.MULTILEVEL * is_csh +\
                ct.Type.INTERFACE * is_interface + \
                ct.Type.SOLID * is_solid
        
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        self.solid.prev_calc_c = deepcopy(self.phrqc.solid_phase_conc['calcite'])
        self.nodetype = self.solid.nodetype
        self.update_border()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)  
        
    def update_diffusivity(self):
        Dref = self.Dref        
        if self.ptype == 'CH':
            cc = self.settings['diffusivity']['CC']
            ch = self.settings['diffusivity']['CH']
            #D_border = self.Dref
            #if('border' in self.settings['diffusivity']):
            #    D_border = self.settings['diffusivity']['border']
            
            is_border = self.solid.border
            is_port = (self.solid.portlandite.c >0) & (~is_border)
            is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
            #is_calc = np.logical_and(is_calc,~is_border)
            is_liquid = np.logical_and(~is_port, ~is_calc)
            is_liquid = np.logical_and(is_liquid, ~is_border)
            
            D_archie = Dref * self.solid.poros * self.solid.app_tort 
            D_CC = D_archie
            D_CH = D_archie
            if(cc[0] == 'const'):
                D_CC = cc[1]*np.ones(D_CC.shape)
            elif(cc[0] == 'inverse'):
                mineral = self.solid.calcite.c * self.solid.calcite.mvol
                #D_CC =np.nan_to_num(1./((1-mineral)/Dref/self.solid.poros/self.solid.app_tort + mineral/cc[1]), Dref)
                D_CC =np.nan_to_num(1./((1-mineral)/Dref + mineral/cc[1]), Dref)
            
            if(ch[0] == 'const'):
                D_CH = ch[1]*np.ones(D_CH.shape)
            elif(ch[0] == 'inverse'):            
                mineral = self.solid.vol
                D_CH =np.nan_to_num(1./((1-mineral)/Dref/self.solid.poros/self.solid.app_tort + mineral/ch[1]), Dref)
            
            De = D_CH*is_port + D_CC*is_calc + D_CC*is_border + Dref*is_liquid 
            De = De/self.solid.poros/self.solid.app_tort 
                
            self.fluid.set_attr('D0',De,component_dict=False)
            self.fluid.set_attr('Deref',np.max(De),component_dict=False)
            self.fluid.call("_set_relaxation_params")           
        if self.ptype == 'CSH':
            cc = self.settings['diffusivity']['CC']
            ch = self.settings['diffusivity']['CH']
            csh = self.settings['diffusivity']['CSH']
            #D_border = self.Dref
            #if('border' in self.settings['diffusivity']):
            #    D_border = self.settings['diffusivity']['border']
            is_border = self.solid.border
            is_port = (self.solid.portlandite.c >0) & (~is_border)
            is_csh = (self.solid.CSHQ_JenD.c > 0) | (self.solid.CSHQ_JenH.c > 0) | \
                (self.solid.CSHQ_TobD.c > 0) |(self.solid.CSHQ_TobH.c > 0) 
            is_csh =  np.logical_and(is_csh, ~is_port)
            is_csh =  np.logical_and(is_csh, ~is_border)
            is_calc = np.logical_and(self.solid.calcite.c >0, ~is_port)
            is_calc = np.logical_and(is_calc, ~is_csh)
            is_calc = np.logical_and(is_calc,~is_border)
            is_liquid = np.logical_and(~is_port, ~is_calc)
            is_liquid = np.logical_and(is_liquid, ~is_csh)
            is_liquid = np.logical_and(is_liquid, ~is_border)
            
            D_archie = Dref * self.solid.poros * self.solid.app_tort 
            D_CC = D_archie
            D_CH = D_archie
            D_CSH = D_archie
            if(cc[0] == 'const'):
                D_CC = cc[1]*np.ones(D_CC.shape)
            elif(cc[0] == 'inverse'):
                mineral = self.solid.calcite.c * self.solid.calcite.mvol
                #D_CC =np.nan_to_num(1./((1-mineral)/Dref/self.solid.poros/self.solid.app_tort + mineral/cc[1]), Dref)
                D_CC =np.nan_to_num(1./((1-mineral)/Dref + mineral/cc[1]), Dref)
            
            if(ch[0] == 'const'):
                D_CH = ch[1]*np.ones(D_CH.shape)
            elif(ch[0] == 'inverse'):            
                mineral = self.solid.vol
                D_CH =np.nan_to_num(1./((1-mineral)/Dref/self.solid.poros/self.solid.app_tort + \
                                        mineral/ch[1]), Dref)           
            
            if(csh[0] == 'const'):
                D_CSH = csh[1]*np.ones(D_CSH.shape)
            
            De = D_CH*is_port + D_CSH*is_csh + \
                D_CC*is_calc + D_CC*is_border + Dref*is_liquid 
                
            De = De/self.solid.poros/self.solid.app_tort 
            #print(De)
            self.fluid.set_attr('D0',De,component_dict=False)
            #self.fluid.set_attr('Deref',np.max(De),component_dict=False)
            self.fluid.set_attr('Deref',Dref,component_dict=False)
            self.fluid.call("_set_relaxation_params")           
        
        
    def update_velocity(self):
        '''
        Optional function that defines velocity change due to calcite precipitation.
        TODO: optmize
        TODO: C-S-H phase
        '''
        if(self.settings['velocity']):
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
        else:
            pass
        
    
    def update_target_SI(self):
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
    
    def update_no_flux(self, ss):        
        ss['Ca'] = ss['Ca']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1)
        if (self.ptype =='CSH'):
            ss['Si'] = ss['Si']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1) 
        return(ss)
        
    def update_border(self):
        if self.ptype == 'CH':
            is_port = (self.solid.portlandite.c>0)
            is_mineral = is_port | (self.solid.nodetype==ct.Type.SOLID)
            val = -16
            temp = deepcopy(self.solid.nodetype)
            temp = temp*(~is_mineral) + val*is_mineral        
            rolled = np.roll(temp, -1, axis = 1)/4+np.roll(temp, 1, axis = 1)/4 +\
                  np.roll(temp, -1, axis = 0)/4+np.roll(temp, 1, axis = 0)/4
            border = (rolled!=val)&is_port
            self.solid.border = border
        
        if self.ptype == 'CSH':
            is_port = (self.solid.portlandite.c>0)
            is_csh = (self.solid.CSHQ_JenD.c > 0) | (self.solid.CSHQ_JenH.c > 0) | \
                (self.solid.CSHQ_TobD.c > 0) |(self.solid.CSHQ_TobH.c > 0) 
            is_mineral = is_port | is_csh | (self.solid.nodetype==ct.Type.SOLID)
            val = -16
            temp = deepcopy(self.solid.nodetype)
            temp = temp*(~is_mineral) + val*is_mineral        
            rolled = np.roll(temp, -1, axis = 1)/4+np.roll(temp, 1, axis = 1)/4 +\
                  np.roll(temp, -1, axis = 0)/4+np.roll(temp, 1, axis = 0)/4
            border = (rolled!=val)&(is_port|is_csh)
            self.solid.border = border
    
    def update_border_solution(self,c,ss):
        #poros = self.solid.poros#
        phrqc_poros = self.solid.poros#self.phrqc.selected_output()['poros'] #
        fraction = self.settings['subgrid']['fraction']
        result = {}
        by = np.where(self.solid.border)[0]
        bx = np.where(self.solid.border)[1]
        bf = np.where(self.solid.border.flatten())[0]
        for i in np.arange(0, np.sum(self.solid.border)):          
            if fraction >0:    
                cell_m = bf[i]+1                
                result  = self.update_equilibrium(result, cell_m,  
                                          self.solid.portlandite.c[by[i], bx[i]],
                                          phrqc_poros[by[i], bx[i]],
                                          fraction)
                ssnew={}
                for name in self.phrqc.components:
                    ssnew[name] = (result[str(cell_m)][name]-c[name][by[i], bx[i]])/self.dt
                    ssnew[name] *= phrqc_poros[by[i], bx[i]]
                    ss[name][by[i], bx[i]] = ssnew[name]
        
        for i in np.arange(0, np.sum(self.solid.border)):
            self.solid.portlandite.c[by[i], bx[i]] = result[str(bf[i]+1)]['portlandite_m']
            #self.solid.calcite.c[by[i], bx[i]] = result[str(bf[i]+1)]['calcite_m']
        return ss
    
    def update_equilibrium(self, result, n_ch, m_ch, porosity, fraction=1):
        ncell = 123456
        fract = fraction/porosity
        if fract>1: fract =1
        #print(fract)
        modify_str = []
        modify_str.append("EQUILIBRIUM_PHASES %i" % ncell)
        modify_str.append("Portlandite 0 %.20e dissolve only" %(m_ch))   
        modify_str.append("END") 
        modify_str.append('MIX %i' % ncell)         
        modify_str.append('%i %.20e' %( n_ch, fract)) #modify_str.append('%i 1' %n_int)  
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
        #print(modify_str)
        #print(output)
        
        port = output[2][10] 
        modify_str = [] 
        modify_str.append("EQUILIBRIUM_PHASES %i" %n_ch)
        modify_str.append("Portlandite 0 %.20e" %(0)) # dissolve only
        modify_str.append("END") 
        modify_str.append('USE equilibrium_phase %i' %n_ch)      
        modify_str.append('MIX %i' % ncell)      
        modify_str.append('%i 1' % ncell)   
        modify_str.append('%i %.20e' %(n_ch, 1.-fract)) #modify_str.append('%i 0' %n_int)   
        modify_str.append('SAVE solution %i' %(n_ch))  
        modify_str.append('SAVE equilibrium_phase %i' %(n_ch))  
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str)  
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        #print(modify_str)
        #print(output)
        comp = {}        
        comp['portlandite_m'] = port
        comp['C'] = output[1][6]
        comp['Ca'] = output[1][7]
        comp['H'] = (output[1][8] - self.phrqc.H_norm)
        comp['O'] = (output[1][9] - self.phrqc.O_norm)
        result[str(n_ch)] = comp        
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
        #is_port = (self.solid.vol_ch >0)
        not_critical = (pore_size >= self.solid.threshold_pore_size) #rt.solid.radius_max 
        #si[np.logical_or(is_port, not_critical)]= 0
        si[not_critical]= 0
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
                