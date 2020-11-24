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
#%% PORE-SIZE CONTROLLED SOLUBILITY     
class PCS():
    def __init__(self, settings, solid, dx):
        self.dx = dx
        self.pore_type = settings['pcs_mode']['pores']
        self.internal_energy = settings['pcs_mode']['int_energy']
        if (self.pore_type == 'block'):
            self.threshold_crystal = settings['pcs_mode']['crystal_size']
            self.threshold_pore_size = settings['pcs_mode']['pore_size']
            self.threshold_distance = self.threshold_pore_size*2            
            self.block_size = self.get_block_size()
            self.free_vol_cc = (solid.voxel_vol -
                                solid.dissolving_mineral_vol)* dx**3
            self.ncrystals = self.get_block_amount(solid)   
            self.crystal_size = self.get_crystal_size(solid)
            self.distance = self.get_distance()
            self.pore_size = self.get_pore_size_block()
        if (self.pore_type == 'cylinder'):
            self.threshold_pore_size = settings['pcs_mode']['pore_size']
            self.pore_length = settings['pcs_mode']['crystal_size']  
            self.pore_density = settings['pcs_mode']['pore_density']
            self.solid.pore_amount = self.get_pore_amount()  
            self.free_vol_cc = (solid.voxel_vol-
                                solid.dissolving_mineral_vol)*dx**3
            self.pore_size = self.pore_size_cylinder()
            self.pore_volume_cc = self.pore_volume_CC() 
    
    def update_target_SI(self, solid):
        if(self.pore_type=='block'): 
            # Pores are considered as a space between calcite crystals   
			self.free_vol_cc = (solid.voxel_vol-
                                solid.dissolving_mineral_vol)*self.dx**3
			self.ncrystals = self.get_block_amount(solid)                 
			self.crystal_size = self.get_crystal_size(solid)
			self.distance = self.get_distance()
			self.pore_size =  self.get_pore_size_block()
        elif (self.pore_type=='cylinder'):
            # Pores are considered as cylinders 
			self.pore_amount = self.pore_amount()  
			self.free_vol_cc = (solid.voxel_vol-
                                solid.dissolving_mineral_vol)*self.dx**3
			self.solid.pore_size = self.pore_size_cylinder()
			self.solid.pore_volume_cc = self.pore_volume_CC(solid) 
        else:
			pass  
        target_SI = self.target_SI()
        return(target_SI)
    
    def target_SI(self):    
        '''
        Target saturation index 
        '''    
        pore_size = self.pore_size
        si = self.saturation_index(pore_size)  
        not_critical = (pore_size >= self.threshold_pore_size) 
        si[not_critical]= 0
        return si  
    
    def saturation_index(self, size):   
        '''
        Target saturation index 
        '''
        sr = self.saturation_ratio(size,  self.internal_energy, 
                                      df.ANGLE,  df.MVOL_CC,  
                                      df.GAS_CONSTANT, df.TEMPERATURE)
        si = np.log10(sr)
        return si
    
    @staticmethod
    def saturation_ratio(radius, int_energy=0.094, angle=np.pi, 
                         mvol=3.69e-5, gas_c=8.314, temp_kelv=298.3):
        '''
        Saturation ratio for pore-size controlled solubility
        '''
        sr = np.exp(-1*mvol*np.cos(angle)* 2 * 
                       int_energy / (gas_c * temp_kelv * radius ))
        return sr
    
    # block    
    def get_block_size(self):
        return(self.threshold_crystal + self.threshold_distance)
        
    def get_pore_size_block(self): 
        return(self.distance/2)
        
    def get_block_amount(self, solid): 
        v_cp = (self.threshold_crystal+self.threshold_distance)**3
        N = self.free_vol_cc/v_cp
        return N
    
    def get_crystal_size(self, solid):
        cs = (solid.calcite.vol/self.ncrystals)**(1./3.) * self.dx 
        return(cs)
    
    def get_distance(self):   
        return(self.block_size - self.crystal_size)
        
    def get_block_porosity(self):
        v_cp = ((self.crystal_size+self.distance)/self.dx )**3
        v_c = (self.crystal_size+self.dx )**3
        p = 1 - v_c/v_cp
        return p
       
    # cylinder        
    def pore_size_cylinder(self, solid): 
        pore_size = self.get_cylinder_radius(solid.poros, 
                                             self.pore_amount, #self.settings['pcs']['pore_density']
                                             self.pore_length, 
                                             self.dx)
        return(pore_size)
    
    @staticmethod
    def get_cylinder_radius(pore_vol, pore_num, pore_length, dx):
        '''
        Formula for cylinder pore radius
        Pore volume is total free volume or porosity
        '''
        r = (pore_vol *dx**3 / (pore_num * np.pi * pore_length))**(1./2.) 
        return(r) 
    
    def pore_volume_CC(self, solid):
        v = 0
        v = self.pore_amount * np.pi * self.pore_size**2 * \
                self.pore_length #*2
        return(v)

    def pore_amount(self, solid):
        '''
        Return matrix of pore amounts per each cell
        '''  
        vol =(solid.voxel_vol - solid.dissolving_mineral_vol)#*self.dx**3
        n = self.pore_density * vol 
        n[n<1] = 1
        return(n)
            
#%% Carbonation 
class CarbonationRT(PhrqcReactiveTransport):
                
    #%% INITIALIZATION      
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        '''
        Reactive transport class for carbonation
        A Lattice Boltzmann method based reactive transport model \
        where reactions are computed using geochemical solver PHREEQC
        '''
        if(eqn == 'MultilevelAdvectionDiffusion'): 
            eqn = 'MultilevelAdvectionDiffusion2' # replace with corrected solver
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
        self.nodetype = deepcopy(domain.nodetype)  
        self.apply_settings(settings)  
        
        self.update_volume()
        self.update_porosity()  
        self.update_app_tort() 
        self.update_init_phrqc()        
        self.update_bc(settings)
        self.update_phases()
        self.update_nodetype()
        
        self.solid.pcs = PCS(settings, self.solid, self.dx)
        self.update_target_SI() 
        
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
            self.phrqc.nodetype = self.nodetype 
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
            
    def update_init_phrqc(self):
        pass
    
    def update_boundary_cells(self):
        nodes = np.zeros(np.shape(self.nodetype))
        nodes[0,:]  = ct.Type.SOLID
        nodes[-1,:] = ct.Type.SOLID
        nodes[:,0]  = ct.Type.SOLID
        nodes[:,-1] = ct.Type.SOLID
        return nodes        
     
    def update_bc(self, settings):
        pass
    #%% LOOP STEP
    
    def advance(self):
        self.update_diffusivity() 
        self.fluid.call('advance') 
        self.update_source()       
        #self.update_phases()    
        #self.update_velocity()
        self.update_nodetype()
        
        
    def update_diffusivity(self):
        pass        
      
    def update_source(self):   
        pass
    
    def update_nodetype(self):
        pass  
     
    def update_phrqc(self):
        self.phrqc.poros=self.solid.poros 
        self.phrqc.is_calc = self.solid.calcite.c>0   
        if(self.phrqc.precipitation == 'interface'):
            if(self.phrqc.active != 'interface'):
                self.phrqc.nodetype = self.nodetype
            
    def update_no_flux(self, ss):        
        return(ss)
        
    def update_border(self):
        pass

    #%% KINETICS    
    def update_border_solution(self,c,ss):
        return(ss)
    
    def update_equilibrium(self, result, n_ch, m_ch, porosity, fraction=1):
        return(result) 
        
    #%% SOLID PARAMETERS
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.solid.prev_vol = deepcopy(self.solid.vol)
        self.update_volume()
        self.update_porosity()
        self.update_app_tort()
        self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
    def update_volume(self):
        pass
            
    def update_porosity(self):        
        self.solid.poros=1.- self.solid.vol/self.solid.voxel_vol
        if np.any(self.solid.poros<=1e-10):
            sys.exit('Negative or zero porosity')
        
    def update_app_tort(self, ):  
        self.solid.app_tort = 1. * self.solid.poros ** self.solid.app_tort_degree
                 
    def update_target_SI(self):
        self.solid.target_SI = self.solid.pcs.update_target_SI(self.solid)
        self.phrqc._target_SI = self.solid.target_SI 
                  
    #%% OPTIONAL OR DEPRECATED
    def update_velocity(self):
        '''
        Optional function that defines velocity change due to calcite precipitation.
        !Need to optmize and add C-S-H phase
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
        
        
    def update_phases(self, thres = 1.0e-3):
        '''
        Optional function to plot dominated phase.
        !Need to optmize and add C-S-H phase
        '''
        pass
        '''
        ch = self.solid.portlandite.c * self.solid.portlandite.mvol
        cc = self.solid.calcite.c * self.solid.calcite.mvol
        tot = cc + ch 
        is_cl = (self.nodetype == 1) #clinker
        is_liq = (tot < thres)* ((self.nodetype ==-2)+(self.nodetype ==-1))
        phases = is_cl*2 + is_liq*(-1)
        
        is_cc = (cc==np.maximum.reduce([cc,ch])) * ~is_liq * ~is_cl
        is_ch = (ch==np.maximum.reduce([cc,ch])) * ~is_liq * ~is_cl 
        phases += is_ch*(-10) + is_cc*(-15)
        self.solid.phases = phases
        '''
        '''
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
        '''
        