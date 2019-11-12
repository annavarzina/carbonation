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
        self.ptype = 'CSH' if hasattr(self.fluid, 'Si') else 'CH'
        self.update_solid_params()    
        self.nodetype = deepcopy(domain.nodetype) 
        self.apply_settings(settings)  
        self.update_nodetype()
        self.update_border_and_interface(self.nodetype)
        self.fluid.call('_set_relaxation_params')  
        
    
    def advance(self):
        #print(self.iters)
        self.update_border_and_interface(self.nodetype)
        self.update_diffusivity() 
        self.fluid.call('_set_relaxation_params')        
        self.fluid.call('advance')    
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):
            self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid.poros)
        
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
        self.set_volume()
        self.set_porosity()        
        self.update_nodetype()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)      
        self.fluid.set_attr('ss',ss) 
        
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
     
    def apply_settings(self, settings):
        self.settings = settings
        self.Dref = settings['Dref']
        self.dx =  settings['dx']
                
    #%% UPDATES
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.set_volume()
        self.set_porosity()
        self.set_app_tort() 
    
    def update_nodetype(self):
        prev_nodetype = self.nodetype
        is_port = (self.solid.portlandite.c>0)
        is_solid = self.solid.nodetype == ct.Type.SOLID
        is_interface = np.zeros(np.shape(is_port))
        #calc_c = self.phrqc.solid_phase_conc['calcite']
        if self.ptype == 'CH':
            is_liquid = (~is_port)
            if(self.iters>1):
                is_liquid =  (~is_port)&(~is_solid)
            self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
                ct.Type.MULTILEVEL * is_port +\
                ct.Type.INTERFACE * is_interface + \
                ct.Type.SOLID * is_solid
        if self.ptype == 'CSH':
            pass
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)      
        self.nodetype = self.solid.nodetype

    def update_diffusivity(self):
        if(self.settings['diffusivity']['type']=='fixed'): 
            self.update_diffusivity_fixed()
        elif(self.settings['diffusivity']['type']=='archie'):
            pass
        else:
            print('Incorrect diffusivity type. Archie\'s relationship is used.')
            
    def update_diffusivity_fixed(self):
        Dref = self.Dref
        D_border = self.Dref
        D_CH = None
        if('D_border' in self.settings['diffusivity']):
            D_border = self.settings['diffusivity']['D_border']
        if('D_CH' in self.settings['diffusivity']):
            D_CH = self.settings['diffusivity']['D_CH']
        
        is_border = self.solid.border
        is_port = (self.solid.portlandite.c >0) & (~is_border)
        is_liquid = np.logical_and(~is_port, ~is_border)
        
        Dnew_lb = Dref*np.ones(np.shape(self.nodetype))
        if(D_CH is None): 
            De = D_border*is_border + Dref*is_liquid 
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
            Dnew_lb = Dref*is_port + Dnew_lb*(~is_port)
        elif(D_CH is not None):
            De = D_CH*is_port+ D_border*is_border + Dref*is_liquid 
            Dnew_lb = De/self.solid.poros/self.solid.app_tort
            
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        #self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)  
        
        
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
                if(fraction <= 1e-15):
                    fraction = 1e-6
                print(fraction)
            if (self.solid.interface['down'][by[i], bx[i]]):
                cell_i = df[i]+1-lx
                cell_m = df[i]+1
                result  = self.update_neighbour_solution(result,cell_i, cell_m,  
                                      self.solid.portlandite.c[by[i], bx[i]],
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
                                      fraction)
                ssnew={}
                for name in self.phrqc.components:
                    ssnew[name] = (result[str(cell_m) + ' ' +str(cell_i)][name]-c[name][by[i], bx[i]+1])/self.dt
                    ssnew[name] *= phrqc_poros[by[i], bx[i]+1]
                    ss[name][by[i], bx[i]+1] = ssnew[name]
        
        for i in np.arange(0, np.sum(self.solid.border)):
            if (self.solid.interface['down'][by[i], bx[i]]):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]-lx+1)]['portlandite_m']
            if (self.solid.interface['up'][by[i], bx[i]]):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]+lx+1)]['portlandite_m']
            if (self.solid.interface['left'][by[i], bx[i]]):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i])]['portlandite_m']  
            if (self.solid.interface['right'][by[i], bx[i]]):
                self.solid.portlandite.c[by[i], bx[i]] = result[str(df[i]+1) + ' ' +str(df[i]+2)]['portlandite_m'] 
       
        return ss
        
    def update_neighbour_solution(self, result, n_int, n_ch, m_ch, fraction=1):
        ncell = 123456
        modify_str = []
        modify_str.append("EQUILIBRIUM_PHASES %i" % ncell)
        modify_str.append("Portlandite 0 %.20e dissolve only" %(m_ch))   
        modify_str.append("END") 
        modify_str.append('MIX %i' % ncell)         
        modify_str.append('%i %.20e' %( n_int, fraction)) #modify_str.append('%i 1' %n_int)  
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
        
        port = output[2][9] 
        modify_str = [] 
        modify_str.append("EQUILIBRIUM_PHASES %i" %n_ch)
        modify_str.append("Portlandite 0 %.20e" %(port)) # dissolve only
        modify_str.append("END") 
        modify_str.append('USE equilibrium_phase %i' %n_int)      
        modify_str.append('MIX %i' % ncell)      
        modify_str.append('%i 1' % ncell)   
        modify_str.append('%i %.20e' %(n_int, 1.-fraction)) #modify_str.append('%i 0' %n_int)   
        modify_str.append('SAVE solution %i' %(n_int))  
        modify_str.append('SAVE equilibrium_phase %i' %n_int)
        modify_str.append('SAVE equilibrium_phase %i' %(n_ch))  
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str)  
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        #print(modify_str)
        #print(output)
        comp = {}        
        comp['portlandite_m'] = port
        comp['Ca'] = output[1][6]
        comp['H'] = (output[1][7] - self.phrqc.H_norm)
        comp['O'] = (output[1][8] - self.phrqc.O_norm)
        result[str(n_ch) + ' ' + str(n_int)] = comp        
        return(result)     
    
    #%% PROPERIES
    def volume_CH(self):
        CH_vol = self.solid.portlandite.c * self.solid.portlandite.mvol
        return CH_vol
    
    def volume_CSH(self):
        CSH_vol = self.solid.CSHQ_JenD.c * self.solid.CSHQ_JenD.mvol +\
                  self.solid.CSHQ_JenH.c * self.solid.CSHQ_JenH.mvol +\
                  self.solid.CSHQ_TobD.c * self.solid.CSHQ_TobD.mvol +\
                  self.solid.CSHQ_TobH.c * self.solid.CSHQ_TobH.mvol 
        return CSH_vol

       
    def free_volume(self):
        v = self.solid.voxel_vol - self.solid.vol_ch
        return v
    
    def csh_conc(self):    
        csh = self.solid.CSHQ_TobH.c+ self.solid.CSHQ_TobD.c +\
                self.solid.CSHQ_JenH.c + self.solid.CSHQ_JenD.c
        return csh 
                