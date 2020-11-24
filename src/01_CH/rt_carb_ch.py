# -*- coding: utf-8 -*-

import sys, os
src_dir = os.path.dirname(__file__)
sys.path.append(src_dir)
import numpy as np
from copy import deepcopy
import yantra
import cell_type as ct # change the path to cell_type file
from rt import CarbonationRT


class CH_Carbonation(CarbonationRT):  
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        super(CH_Carbonation, self).__init__(eqn,domain,domain_params,bc_params,solver_params, settings)
        self.ptype ='CH'
     
    def update_bc(self, settings):
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
    
    def update_init_phrqc(self):
        self.phrqc.boundcells = self.update_boundary_cells()
        if 'portlandite' in self.solid.diffusive_phase_list:
            self.phrqc.init_port = self.solid.portlandite.c>0
        if 'calcite' in self.solid.diffusive_phase_list:
            self.phrqc.is_calc = self.solid.calcite.c>0
        self.phrqc._target_SI = np.zeros(self.solid.shape)  
        self.phrqc.nodetype = self.nodetype 
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty) 
            
            
    def update_source(self):                
        c=self.fluid.get_attr('c')
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = self.solid.diffusive_phase_list
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
        
    def update_nodetype(self):
        '''
        find neighbous of the multilevel cells
        '''
        prev_nodetype = deepcopy(self.nodetype) 
        
        self.solid.is_port = (self.solid.portlandite.c>0)
        self.solid.is_calc = (self.solid.calcite.c>0)
        self.solid.is_solid = (self.solid.nodetype == ct.Type.SOLID)
        self.solid.is_critical = (self.solid.nodetype == ct.Type.SOLID)
        is_port = self.solid.is_port
        is_calc = self.solid.is_calc
        #is_clinker = self.solid.nodetype == ct.Type.CLINKER
        is_solid = self.solid.is_solid
        is_critical = np.zeros(np.shape(is_port))
        is_interface = np.zeros(np.shape(is_port))
        #calc_c = self.phrqc.solid_phase_conc['calcite']
        is_liquid = (~is_port)
        if(self.iters>1):
            is_critical = (self.solid.pcs.pore_size <= 
                           self.solid.pcs.threshold_pore_size) & is_calc & (~is_port)
            is_liquid =  (~is_critical)&(~is_port)&(~is_solid)&(~is_calc)#&((prev_nodetype==-1)|(prev_nodetype==-2))
            is_interface = (~is_critical)&(~is_port)&(~is_solid)&(~is_liquid)
        self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
            ct.Type.MULTILEVEL * is_critical + \
            ct.Type.MULTILEVEL * is_port +\
            ct.Type.INTERFACE * is_interface + \
            ct.Type.SOLID * is_solid
        self.solid.is_liquid = is_liquid
        self.solid.is_interface = is_interface
        self.solid.is_critical = is_critical 
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        #self.solid.prev_calc_c = deepcopy(self.phrqc.solid_phase_conc['calcite']) 
        self.nodetype = self.solid.nodetype
        self.update_border()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)  
        self.solid.prev_nodetype = prev_nodetype
        
        
    def update_diffusivity(self):
        Dref = self.Dref        
        cc = self.settings['diffusivity']['CC']
        ch = self.settings['diffusivity']['CH']
        
        is_border = self.solid.border
        is_port = (self.solid.portlandite.c >0) & (~is_border)
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
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
        
    def update_no_flux(self, ss):        
        ss['Ca'] = ss['Ca']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1) 
        return(ss)
        
    def update_border(self):
        is_port = (self.solid.portlandite.c>0)
        is_mineral = is_port | (self.solid.nodetype==ct.Type.SOLID)
        val = -16
        temp = deepcopy(self.solid.nodetype) 
        temp = temp*(~is_mineral) + val*is_mineral        
        rolled = np.roll(temp, -1, axis = 1)/4+np.roll(temp, 1, axis = 1)/4 +\
              np.roll(temp, -1, axis = 0)/4+np.roll(temp, 1, axis = 0)/4
        border = (rolled!=val)&is_port
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
        
    def update_volume(self):
        vol = np.zeros(self.solid.shape)
        phase_list = self.solid.diffusive_phase_list
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            vol += val.c * self.solid.mvol[num-1] 
        self.solid.vol = vol
        self.solid.portlandite.vol = self.volume_CH()
        self.solid.calcite.vol = self.volume_CC() 
        self.solid.dissolving_mineral_vol = self.solid.portlandite.vol
        
    def volume_CH(self):
        CH_vol = self.solid.portlandite.c * self.solid.portlandite.mvol
        return CH_vol
    
    def volume_CC(self):
        CC_vol = self.solid.calcite.c * self.solid.calcite.mvol
        return CC_vol
       
    def free_volume(self):
        v = self.solid.voxel_vol - self.solid.vol_ch
        return v

    
    def csh_conc(self):    
        csh = self.solid.CSHQ_TobH.c+ self.solid.CSHQ_TobD.c +\
                self.solid.CSHQ_JenH.c + self.solid.CSHQ_JenD.c
        return csh 
                