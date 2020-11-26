# -*- coding: utf-8 -*-
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
            eqn = 'MultilevelAdvectionDiffusion2'
        self.auto_time_step = solver_params.get('auto_time_step',True)
        self.phrqc = Phrqc(domain,domain_params,bc_params,solver_params)    
        components = self.phrqc.components
        bc = self.phrqc.boundary_conditions 
        #print(self.phrqc.nactive)
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
        self.app_tort_degree = settings['app_tort']['degree']
        self.update_solid_params()    
        self.nodetype = deepcopy(domain.nodetype) 
        self.apply_settings(settings)  
        self.update_nodetype()
        self.update_border_and_interface(self.nodetype)
        self.solid.prev_border = deepcopy(self.solid.border)
        self.Dref = settings['Dref']
    
    def advance(self):
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
        #self.fluid.call("_set_relaxation_params")
        
    
    def update_source(self):        
        c = self.fluid.get_attr('c')
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1])
        ss = self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        if self.iters >1 :pqty=self.solid.update(self.phrqc.dphases)  
        ss = self.update_equilibrium(c,ss)
        self.set_volume()
        self.set_porosity()  
        self.set_app_tort()        
        self.update_nodetype()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)      
        self.fluid.set_attr('ss',ss) 
        
    def modify_eq(self, phaseqty):
        """
        modifies the phaseqty in phreeqc
        
        Parameters
        ----------
        phaseqty: dict
            dictionary containing an ndarray of new quantities of equilibrium phases
        """
        phaseqty = self.phrqc.flatten_dict(phaseqty)
        modifystr = []
        bord = np.where(self.solid.border.flatten())[0] + 1
        for cell in range(self.phrqc.startcell,self.phrqc.stopcell+1,1):
            if (cell not in bord):
                modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)
                for key in phaseqty.keys():
                    modifystr.append("\t -component %s" %(key))
                    modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1]))
                    #modifystr.append("\t\t%s\t%.20e" %('-si', 2))
        modifystr.append("end")       
        modifystr ='\n'.join(modifystr)
        #print(modifystr)
        self.phrqc.IPhreeqc.RunString(modifystr)
        
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
                
    #%% UPDATES
    
    def update_diffusivity(self):
        pass
        
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.set_volume()
        self.set_porosity()
        self.set_app_tort() 
    
    def update_nodetype(self):
        # override
        pass

    def update_no_flux(self, ss):
        # override        
        return(ss)
        
    def update_border_and_interface(self, nodetype):
        # override
        pass
    
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
    
    def update_border_solution(self,c,ss):
        phrqc_poros = self.phrqc.selected_output()['poros']
        fr = np.zeros(phrqc_poros.shape)
        fraction = self.settings['subgrid']['fraction']
        result = {}
        by = np.where(self.solid.border)[0]
        bx = np.where(self.solid.border)[1]
        bf = np.where(self.solid.border.flatten())[0]
        lx = self.nodetype.shape[1]
        for i in np.arange(0, np.sum(self.solid.border)):
            if fraction is None:
                fraction = 1-phrqc_poros[by[i], bx[i]]
            else:
                if (self.settings['subgrid']['poros']):
                    fraction = fraction - phrqc_poros[by[i], bx[i]]
                        #print(fraction)
            if fraction >0:    
                fr[by[i], bx[i]] = fraction
                if (self.solid.interface['down'][by[i], bx[i]]):
                    cell_i = bf[i]+1-lx
                    cell_m = bf[i]+1
                    result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                          self.solid.portlandite.c[by[i], bx[i]],
                                          fraction)
                    ssnew={}
                    for name in self.phrqc.components:
                        ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i]-1, bx[i]])/self.dt
                        ssnew[name] *= phrqc_poros[by[i]-1, bx[i]]
                        ss[name][by[i]-1, bx[i]] = ssnew[name]
                if (self.solid.interface['up'][by[i], bx[i]]):
                    cell_i = bf[i]+1+lx
                    cell_m = bf[i]+1
                    result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                          self.solid.portlandite.c[by[i], bx[i]],
                                          fraction)#self.solid.poros[1,1])
                    ssnew={}
                    for name in self.phrqc.components:
                        ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i]+1, bx[i]])/self.dt
                        ssnew[name] *= phrqc_poros[by[i]+1, bx[i]]
                        ss[name][by[i]+1, bx[i]] = ssnew[name]
                if (self.solid.interface['left'][by[i], bx[i]]):
                    cell_i = bf[i]
                    cell_m = bf[i]+1
                    result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                          self.solid.portlandite.c[by[i], bx[i]],
                                          fraction)
                    ssnew={}
                    for name in self.phrqc.components:
                        ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i], bx[i]-1])/self.dt
                        ssnew[name] *= phrqc_poros[by[i], bx[i]-1]
                        ss[name][by[i], bx[i]-1] = ssnew[name]
                   
                if (self.solid.interface['right'][by[i], bx[i]]):
                    cell_i = bf[i]+2
                    cell_m = bf[i] +1
                    result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                          self.solid.portlandite.c[by[i], bx[i]],
                                          fraction)
                    ssnew={}
                    for name in self.phrqc.components:
                        ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i], bx[i]+1])/self.dt
                        ssnew[name] *= phrqc_poros[by[i], bx[i]+1]
                        ss[name][by[i], bx[i]+1] = ssnew[name]
        
        for i in np.arange(0, np.sum(self.solid.border)):
            if (self.solid.interface['down'][by[i], bx[i]] and fr[by[i], bx[i]]>0):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(bf[i]+1) + ' ' +str(bf[i]-lx+1)]['portlandite_m']
            if (self.solid.interface['up'][by[i], bx[i]] and fr[by[i], bx[i]]>0):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(bf[i]+1) + ' ' +str(bf[i]+lx+1)]['portlandite_m']
            if (self.solid.interface['left'][by[i], bx[i]] and fr[by[i], bx[i]]>0):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(bf[i]+1) + ' ' +str(bf[i])]['portlandite_m']  
            if (self.solid.interface['right'][by[i], bx[i]] and fr[by[i], bx[i]]>0):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(bf[i]+1) + ' ' +str(bf[i]+2)]['portlandite_m'] 
       
        return ss
        
    
    def update_equilibrium(self,c,ss):
        # override
        return ss
    
    
    def update_pqrqc_equilibrium(self, result, n_ch, m_ch, porosity, fraction=1):
        # override
        return(result)     
    
    #%% PROPERIES
    
       
    '''
    def volume_CH(self):
        CH_vol = self.solid.portlandite.c * self.solid.portlandite.mvol
        return CH_vol
    def free_volume(self):
        v = self.solid.voxel_vol - self.solid.vol_ch
        return v
    def volume_CSH(self):
        CSH_vol = self.solid.CSHQ_JenD.c * self.solid.CSHQ_JenD.mvol +\
                  self.solid.CSHQ_JenH.c * self.solid.CSHQ_JenH.mvol +\
                  self.solid.CSHQ_TobD.c * self.solid.CSHQ_TobD.mvol +\
                  self.solid.CSHQ_TobH.c * self.solid.CSHQ_TobH.mvol 
        return CSH_vol

    
    def csh_conc(self):    
        csh = self.solid.CSHQ_TobH.c+ self.solid.CSHQ_TobD.c +\
                self.solid.CSHQ_JenH.c + self.solid.CSHQ_JenD.c
        return csh 
    '''