# -*- coding: utf-8 -*-

import sys, os
src_dir = os.path.dirname(__file__)
sys.path.append(src_dir)
import numpy as np
from copy import deepcopy
import yantra
import cell_type as ct # change the path to cell_type file
import defaults as df
from rt import CarbonationRT
from phrqc_input import PhreeqcInput

class PhreeqcInputCSHQ(PhreeqcInput):
    
    def make_phrqc_input(self):
        self.phrqc_CSHQ_phases()
        self.phrqc_boundary_voxel(self.c['c_bc'])
        self.phrqc_liquid_voxel(self.c['c_liq'], self.c['ca_liq'])
        self.phrqc_multilevel_voxel_CH(self.c['c_mlvl'], self.c['ca_mlvl'])        
        self.phrqc_multilevel_voxel_CSH()
        self.phrqc_solid_voxel()
        return(self.phrqc_input)
        
    def phrqc_CSHQ_phases(self):
        phrqc_input = [] 
        
        phrqc_input.append('PHASES')  
        phrqc_input.append('CSHQ_TobH')  
        phrqc_input.append('\t(CaO)0.66666667(SiO2)1(H2O)1.5 = 0.66666667Ca++ + 1 SiO(OH)3- + 0.33333334OH- -0.16666667 H2O') 
        phrqc_input.append('\tlog_K -6.190832') 
        phrqc_input.append('CSHQ_TobD') 
        phrqc_input.append('\t(CaO)0.8333333333(SiO2)0.6666666667(H2O)1.8333333333 = 0.8333333333 Ca++ + 0.6666666667 SiO(OH)3- + 0.99999999990 OH- + 0.3333333333 H2O') 
        phrqc_input.append('\tlog_K -6.8995533') 
        phrqc_input.append('CSHQ_JenH') 
        phrqc_input.append('\t(CaO)1.3333333333(SiO2)1(H2O)2.1666666667 = 1.3333333333 Ca++ + 1 SiO(OH)3- + 1.6666666667 OH- -0.1666666667 H2O') 
        phrqc_input.append('\tlog_K -10.96765') 
        phrqc_input.append('CSHQ_JenD') 
        phrqc_input.append('\t(CaO)1.5(SiO2)0.6666666667(H2O)2.5 = 1.5 Ca++ + 0.6666666667 SiO(OH)3- + 2.3333333333 OH- + 0.3333333333 H2O') 
        phrqc_input.append('\tlog_K -10.47635') 
        phrqc_input.append('knobs') 
        phrqc_input.append('\t-iterations 8000\n') 
        
        self.phrqc_input +=  phrqc_input
    
    def phrqc_multilevel_voxel_CH(self, c, ca):
        phrqc_input = [] 
        phrqc_input.append('#solution_multilevel')    
        phrqc_input.append('SOLUTION\t100003')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(c['type'] == 'conc'):
            phrqc_input.append('\tC\t' + str(c['value']))
        elif(c['type'] == 'eq'):
            phrqc_input.append('\tC\t1\t' + str(c['value']))
        else:
            pass              
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100003')
        phrqc_input.append('portlandite\t0\t1')
        phrqc_input.append('calcite\t0\t0\n')
        
        self.phrqc_input +=  phrqc_input
                           
    def phrqc_multilevel_voxel_CSH(self):
        phrqc_input = [] 
        phrqc_input.append('#solution_csh_multilevel') 
        phrqc_input.append('SOLUTION\t100004') 
        phrqc_input.append('\t-water\t0.448230266981165') 
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\tpH\t12\tcharge') 
        phrqc_input.append('\tCa\t1.955e-002') 
        phrqc_input.append('\tSi\t3.018e-005') 
        
        phrqc_input.append('SOLID_SOLUTIONS\t100004') 
        phrqc_input.append('Tob_jen_ss') 
        phrqc_input.append('\t-comp\tCSHQ_TobH\t0.1041') 
        phrqc_input.append('\t-comp\tCSHQ_TobD\t2.5050') 
        phrqc_input.append('\t-comp\tCSHQ_JenH\t2.1555') 
        phrqc_input.append('\t-comp\tCSHQ_JenD\t3.2623') 
        phrqc_input.append('EQUILIBRIUM_PHASES\t100004')
        phrqc_input.append('portlandite\t0\t0')
        phrqc_input.append('calcite\t0\t0\n')
        
        self.phrqc_input +=  phrqc_input

class CSH_Carbonation(CarbonationRT):  
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        super(CSH_Carbonation, self).__init__(eqn,domain,domain_params,bc_params,solver_params, settings)
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
                self.fluid.Si.bc[key] = ['flux', 0.0]
                self.fluid.Si._bc[key+'bc'] = 'flux'
                
    def update_init_phrqc(self):
        self.phrqc.boundcells = self.update_boundary_cells()
        if 'portlandite' in self.solid.diffusive_phase_list:
            self.phrqc.init_port = self.solid.portlandite.c>0
        self.phrqc.init_csh = np.logical_or(np.logical_or(self.solid.CSHQ_TobD.c>0, self.solid.CSHQ_TobH.c>0),
                                                np.logical_or(self.solid.CSHQ_JenD.c>0, self.solid.CSHQ_JenH.c>0))
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
                
        
    def update_phases(self, thres = 1.0e-3):
        pass
        
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
        
        is_solid = self.solid.is_solid
        is_critical = np.zeros(np.shape(is_port))
        is_interface = np.zeros(np.shape(is_port))
        
        is_csh= (self.solid.CSHQ_JenD.c>0) | (self.solid.CSHQ_JenH.c>0) | \
                (self.solid.CSHQ_TobD.c>0) |(self.solid.CSHQ_TobH.c>0) &(~is_port) 
        is_liquid = (~is_port) & (~is_csh)
        if(self.iters>1):
            is_critical = (self.solid.pcs.pore_size <= 
                           self.solid.pcs.threshold_pore_size) & is_calc & (~is_port)
            is_liquid =  (~is_critical)&(~is_port)&(~is_solid)&(~is_calc)&(~is_csh)#&((prev_nodetype==-1)|(prev_nodetype==-2))
            is_interface = (~is_critical)&(~is_port)&(~is_solid)&(~is_liquid)&(~is_csh)
        self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
            ct.Type.MULTILEVEL * is_critical + \
            ct.Type.MULTILEVEL * is_port +\
            ct.Type.MULTILEVEL * is_csh +\
            ct.Type.INTERFACE * is_interface + \
            ct.Type.SOLID * is_solid
        self.solid.is_liquid = is_liquid
        self.solid.is_interface = is_interface
        self.solid.is_critical = is_critical
        self.solid.is_csh = is_csh
        
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
           
    def update_no_flux(self, ss):        
        ss['Ca'] = ss['Ca']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1)
        ss['Si'] = ss['Si']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1) 
        return(ss)
        
    def update_border(self):
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
        return ss
    
    def update_equilibrium(self, result, n_ch, m_ch, porosity, fraction=1):
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
        self.solid.csh_vol = self.volume_CSH()
        self.solid.dissolving_mineral_vol = self.solid.portlandite.vol + \
            self.solid.csh_vol
        
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

       
    def free_volume(self):
        v = self.solid.voxel_vol - self.solid.vol_ch
        return v

    
    def csh_conc(self):    
        csh = self.solid.CSHQ_TobH.c+ self.solid.CSHQ_TobD.c +\
                self.solid.CSHQ_JenH.c + self.solid.CSHQ_JenD.c
        return csh 
                