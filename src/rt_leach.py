import sys
import numpy as np
from copy import deepcopy
import yantra
from yantra._base import Multicomponent 
from yantra.physics.PhrqcReactiveTransport import Solid
from yantra.physics.PhrqcReactiveTransport import PhrqcReactiveTransport 
from yantra.physics._phrqc_wrapper import  Phrqc
#from rt import CarbonationRT
import cell_type as ct # change the path to cell_type file
import defaults as df

class LeachingRT(PhrqcReactiveTransport):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        '''
        Reactive transport class for carbonation
        A Lattice Boltzmann method based reactive transport model \
        where reactions are computed using geochemical solver PHREEQC
        ''' 
        if(settings['diffusivity']['type']=='fixed'):
            print("newMlvl")
            eqn = 'MultilevelAdvectionDiffusion2'
            #solver_params['tauref'] = 1
            #domain_params['Deref'] = settings['diffusivity']['D_border']
            #solver_params['tauref'] = 0.5*np.max(settings['Dref'])/domain_params['Deref']+0.5#5.5 for 1e-10
        self.auto_time_step = solver_params.get('auto_time_step',True)
        self.phrqc = Phrqc(domain,domain_params,bc_params,solver_params)    
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
        self.app_tort_degree = settings['app_tort']['degree'] 
        self.update_solid_params()    
        self.nodetype = deepcopy(domain.nodetype) 
        self.apply_settings(settings)  
        self.update_nodetype()
        self.update_border_and_interface(self.nodetype)
        #self.fluid.call("_set_relaxation_params")
        self.solid.prev_border = deepcopy(self.solid.border)
        self.Dref = settings['Dref']
    
    def advance(self):
        #print(self.iters)
        
        self.update_border_and_interface(self.nodetype) 
        if(~np.all(self.solid.border == self.solid.prev_border)):
            print("update")            
            self.update_diffusivity() 
            self.solid.prev_border = deepcopy(self.solid.border)
        self.fluid.call('advance')    
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):                  
            self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid.poros)
        self.update_source()
        
    def update_source(self):
        pass
            
    def modify_eq(self, phaseqty):
        pass
    
    def update_diffusivity(self):
        pass
        
    def update_nodetype(self):
        pass

    def update_no_flux(self, ss): 
        pass
        
    def update_border_and_interface(self, nodetype):
        pass
    
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
        
    def set_app_tort(self):  
        d = self.app_tort_degree
        self.solid.app_tort = 1. * self.solid.poros ** d
        
    def set_boundary_cells(self):
        nodes = np.zeros(np.shape(self.nodetype))
        nodes[0,:]  = ct.Type.SOLID
        nodes[-1,:] = ct.Type.SOLID
        nodes[:,0]  = ct.Type.SOLID
        nodes[:,-1] = ct.Type.SOLID
        return nodes        
     
    def apply_settings(self, settings):
        self.settings = settings
        self.dx =  settings['dx']
                
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.set_volume()
        self.set_porosity()
        self.set_app_tort() 
    
    
    def get_border(self, nt, val, is_mineral):
        rolled = np.roll(nt, -1, axis = 1)/4+np.roll(nt, 1, axis = 1)/4 +\
              np.roll(nt, -1, axis = 0)/4+np.roll(nt, 1, axis = 0)/4
        is_border = (rolled!=val)&is_mineral
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
    