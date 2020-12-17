# -*- coding: utf-8 -*-

import sys, os
src_dir = os.path.dirname(__file__)
sys.path.append(src_dir)
import numpy as np
import matplotlib.pylab as plt
from copy import deepcopy
import yantra
from sklearn.linear_model import LinearRegression
import cell_type as ct # change the path to cell_type file
import defaults as df
from rt_carb import CarbonationRT
from phrqc_input import PhreeqcInput
from settings import Settings, Results
from yantra._pyevtk.hl  import imageToVTK

class PhreeqcInputCSHQ(PhreeqcInput):
    #TODO check that for subgrid there is no mineral
    def make_phrqc_input(self):
        self.phrqc_CSHQ_phases_cemdata18()
        self.phrqc_boundary_voxel(self.c['c_bc'])
        self.phrqc_liquid_voxel(self.c['c_liq'], self.c['ca_liq'])
        self.phrqc_multilevel_voxel_CH(self.c['c_mlvl'], self.c['ca_mlvl'])        
        self.phrqc_multilevel_voxel_CSH()
        self.phrqc_solid_voxel()
        return(self.phrqc_input)
        
    def phrqc_CSHQ_phases_cemdata18(self):
        phrqc_input = [] 
        
        phrqc_input.append('PHASES')  
        phrqc_input.append('calcite')  
        phrqc_input.append('\tCaCO3 = CO3-2 + Ca+2') 
        phrqc_input.append('\t-log_K\t-8.479966')  
        phrqc_input.append('\t-analytical_expression\t130.276347\t0\t-5689.203921\t-48.36444\t0\t0\t0')  
        phrqc_input.append('sio2am')  
        phrqc_input.append('\tSiO2 = SiO2') 
        phrqc_input.append('\t-log_K\t-2.714066')  
        phrqc_input.append('\t-analytical_expression\t0\t0\t-809.189752\t0\t0\t0\t0') 
        phrqc_input.append('CSHQ_TobH')  
        phrqc_input.append('\t(CaO)0.6667(SiO2)1(H2O)1.5 + 1.3334H+ = 0.6667Ca+2 + 2.1667H2O + SiO2') 
        phrqc_input.append('\t-log_K\t8.286642') 
        phrqc_input.append('\t-analytical_expression\t-12.519254 0 2163.381583 5.476331 0 0 0') 
        phrqc_input.append('CSHQ_TobD') 
        phrqc_input.append('\t((CaO)1.25(SiO2)1(H2O)2.75)0.6667 + 1.66675H+ = 0.833375Ca+2 + 2.6668H2O + 0.6667SiO2') 
        phrqc_input.append('\t-log_K\t13.655314') 
        phrqc_input.append('\t-analytical_expression\t-10.916344 0 3959.367696 4.563888 0 0 0') 
        phrqc_input.append('CSHQ_JenH') 
        phrqc_input.append('\t(CaO)1.3333(SiO2)1(H2O)2.1667 + 2.6666H+ = 1.3333Ca+2 + 3.5H2O + SiO2') 
        phrqc_input.append('\t-log_K\t22.179305') 
        phrqc_input.append('\t-analytical_expression\t-17.10944 0 6470.553982 7.107847 0 0 0') 
        phrqc_input.append('CSHQ_JenD') 
        phrqc_input.append('\t(CaO)1.5(SiO2)0.6667(H2O)2.5 + 3H+ = 1.5Ca+2 + 4H2O + 0.6667SiO2') 
        phrqc_input.append('\t-log_K\t28.730362')  
        phrqc_input.append('\t-analytical_expression\t-15.591756 0 8609.739692 6.24251 0 0 0') 
        
        phrqc_input.append('knobs') 
        phrqc_input.append('\t-iterations 8000\n') 
        self.phrqc_input +=  phrqc_input
    
    def phrqc_CSHQ_phases_cemdata07(self):
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
        phrqc_input.append('sio2am') 
        phrqc_input.append('\tSiO2 + 1OH- + 1H2O = SiO(OH)3-') 
        phrqc_input.append('\tlog_K -1.475988') 
        phrqc_input.append('\t-analytical_expression\t-2.14181238156\t0\t664.05554528339\t0.5620295123\t0') 
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
        phrqc_input.append('\tCa\t1.955e-002')  #0?
        phrqc_input.append('\tSi\t3.018e-005') 
        
        phrqc_input.append('SOLID_SOLUTIONS\t100004') 
        phrqc_input.append('Tob_jen_ss') 
        phrqc_input.append('\t-comp\tCSHQ_TobH\t0')#0.1041') 
        phrqc_input.append('\t-comp\tCSHQ_TobD\t0')#2.5050') 
        phrqc_input.append('\t-comp\tCSHQ_JenH\t0')#2.1555') 
        phrqc_input.append('\t-comp\tCSHQ_JenD\t0')#3.2623') 
        phrqc_input.append('\t-comp\tsio2am\t0.0') 
        phrqc_input.append('EQUILIBRIUM_PHASES\t100004')
        phrqc_input.append('portlandite\t0\t0')
        phrqc_input.append('calcite\t0\t0\n')
        
        self.phrqc_input +=  phrqc_input

class CarbonationCSHQ(CarbonationRT):  
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        super(CarbonationCSHQ, self).__init__(eqn,domain,domain_params,bc_params,solver_params, settings)
        self.ptype ='CSH'
        self.fraction_model = self.model_fraction()
        
    
    def model_fraction(self):
        casi = np.array([1.67, 1.65, 1.60, 1.55, 1.50, 1.45, 1.40, 1.35, 1.30, 1.25,
                  1.20, 1.15, 1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.83])
        d = np.array([-8.40, -8.47, -8.53, -8.53, -8.56, -8.67, -8.89, -9.20, -9.58,
              -9.98, -10.36, -10.67, -10.88, -10.97, -10.95, -10.86, -10.79, 
              -10.87, -10.99])
        s = np.array([865, 880, 925,  970, 1015, 1060, 1105, 1150, 1195, 1240, 1285, 1330,
              1375, 1420, 1465, 1510, 1555, 1600, 1615])
        sigma = (10**d) * s
        
        X = casi.reshape(-1,1) # put your dates in here
        y = sigma # put your kwh in here
        model = LinearRegression()
        model.fit(X, y)       
        return(model)        

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
            phaseqty['CSHQ_TobH'][np.where(self.solid.csh_vol>0)]=0 #border
            phaseqty['CSHQ_TobD'][np.where(self.solid.csh_vol)]=0
            phaseqty['CSHQ_JenH'][np.where(self.solid.csh_vol)]=0
            phaseqty['CSHQ_JenD'][np.where(self.solid.csh_vol)]=0
            #phaseqty['portlandite'][np.where(self.solid.border)]=0
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

       
    def casi_ratio(self):
        ca = 0.8333333/3.3333333333*self.solid.CSHQ_TobD.c[:,:] + \
             0.6666667/3.16666667*self.solid.CSHQ_TobH.c[:,:] + \
             1.3333333/4.5*self.solid.CSHQ_JenH.c[:,:] + \
             1.5/4.666666666699999*self.solid.CSHQ_JenD.c[:,:]
        si = 0.6666667/3.3333333333*self.solid.CSHQ_TobD.c[:,:] + \
             1.0/3.16666667*self.solid.CSHQ_TobH.c[:,:] + \
             1.0/4.5*self.solid.CSHQ_JenH.c[:,:] + \
             0.6666667/4.666666666699999*self.solid.CSHQ_JenD.c[:,:]
        r = ca/si     
        self.solid.casi = r[self.solid.csh_vol>0]
        return(r[self.solid.csh_vol>0])
            
    def predict_fraction(self):
        model = self.fraction_model
        y_predict = np.zeros(np.shape(self.casi_ratio()))
        for i in range(0, len(self.casi_ratio())):
            casi = self.casi_ratio()[i]
            #print(casi)
            if(casi>1.67):
                y_predict[i] = 5.0e-06
            elif(casi<0.83):
                y_predict[i] = 1.0e-11
            else:
                X_predict = np.array([casi]).reshape(-1,1)  # put the dates of which you want to predict kwh here
                y_predict[i] = model.predict(X_predict)
        return(y_predict)
        
    def update_border_solution(self,c,ss): #TODO rename csh
        phrqc_poros = self.solid.poros#self.phrqc.selected_output()['poros'] #
        #fraction = self.settings['subgrid']['fraction']
        #TODO fraction array!
        fraction = self.predict_fraction()
        #print(self.solid.casi)
        #print(fraction)
        result = {}
        by = np.where(self.solid.csh_vol)[0]#np.where(self.solid.border)[0]
        bx = np.where(self.solid.csh_vol)[1]#np.where(self.solid.border)[1]
        bf = np.where(self.solid.csh_vol.flatten())[0]#np.where(self.solid.border.flatten())[0]
        #for i in np.arange(0, np.sum(self.solid.border)):          
        for i in np.arange(0, np.sum(self.solid.csh_vol>0)): 
            if fraction[i] >0:    
                cell_m = bf[i]+1                
                #print(cell_m)
                m_csh = [self.solid.CSHQ_TobH.c[by[i], bx[i]],
                         self.solid.CSHQ_TobD.c[by[i], bx[i]],
                         self.solid.CSHQ_JenH.c[by[i], bx[i]],
                         self.solid.CSHQ_JenD.c[by[i], bx[i]],
                         self.solid.sio2am.c[by[i], bx[i]]]
                result  = self.update_equilibrium(result, cell_m,  
                                          m_csh,
                                          phrqc_poros[by[i], bx[i]],
                                          fraction[i])
                ssnew={}
                for name in self.phrqc.components:
                    ssnew[name] = (result[str(cell_m)][name]-c[name][by[i], bx[i]])/self.dt
                    ssnew[name] *= phrqc_poros[by[i], bx[i]]
                    ss[name][by[i], bx[i]] = ssnew[name]
        
        #for i in np.arange(0, np.sum(self.solid.border)):
        for i in np.arange(0, np.sum(self.solid.csh_vol>0)):
            self.solid.CSHQ_TobH.c[by[i], bx[i]] = result[str(bf[i]+1)]['CSHQ_TobH_m']
            self.solid.CSHQ_TobD.c[by[i], bx[i]] = result[str(bf[i]+1)]['CSHQ_TobD_m']
            self.solid.CSHQ_JenH.c[by[i], bx[i]] = result[str(bf[i]+1)]['CSHQ_JenH_m']
            self.solid.CSHQ_JenD.c[by[i], bx[i]] = result[str(bf[i]+1)]['CSHQ_JenD_m']
            self.solid.sio2am.c[by[i], bx[i]] = result[str(bf[i]+1)]['sio2am_m']
            #self.solid.calcite.c[by[i], bx[i]] = result[str(bf[i]+1)]['calcite_m']
        return ss
    
    def update_equilibrium(self, result, n_csh, m_csh, porosity, fraction=1):
        ncell = 123456
        fract = fraction/porosity
        if fract>1: fract =1
        #print(fract)
        modify_str = []
        modify_str.append("SOLID_SOLUTIONS %i" % ncell)
        modify_str.append('Tob_jen_ss') #TODO dissolve only?
        modify_str.append("-comp\tCSHQ_TobH\t%.20e" %(m_csh[0]))
        modify_str.append("-comp\tCSHQ_TobD\t%.20e" %(m_csh[1]))
        modify_str.append("-comp\tCSHQ_JenH\t%.20e" %(m_csh[2]))
        modify_str.append("-comp\tCSHQ_JenD\t%.20e" %(m_csh[3]))
        modify_str.append("-comp\tsio2am\t%.20e" %(m_csh[4]))
        modify_str.append("END") 
        modify_str.append('MIX %i' % ncell)         
        modify_str.append('%i %.20e' %( n_csh, fract)) #modify_str.append('%i 1' %n_int)  
        modify_str.append('SAVE SOLUTION %i' % ncell)  
        modify_str.append("END") 
        modify_str.append('USE SOLUTION %i' % ncell)  
        modify_str.append('USE SOLID_SOLUTIONS %i' % ncell)
        modify_str.append('SAVE SOLUTION %i' % ncell)   
        modify_str.append("END") 
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str) 
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        #print(modify_str)
        #print(output)
        
        tobh = output[2][15]
        tobd = output[2][17]
        jenh = output[2][19]
        jend = output[2][21] 
        sio2 = output[2][23]
        modify_str = [] 
        modify_str.append("SOLID_SOLUTIONS %i" %n_csh)
        modify_str.append('Tob_jen_ss')
        modify_str.append("-comp\tCSHQ_TobH\t%.20e" %(0))
        modify_str.append("-comp\tCSHQ_TobD\t%.20e" %(0))
        modify_str.append("-comp\tCSHQ_JenH\t%.20e" %(0))
        modify_str.append("-comp\tCSHQ_JenD\t%.20e" %(0))
        modify_str.append("-comp\tsio2am\t%.20e" %(0))
        modify_str.append("END") 
        modify_str.append('USE SOLID_SOLUTIONS %i' %n_csh)      
        modify_str.append('MIX %i' % ncell)      
        modify_str.append('%i 1' % ncell)   
        modify_str.append('%i %.20e' %(n_csh, 1.-fract)) #modify_str.append('%i 0' %n_int)   
        modify_str.append('SAVE SOLUTION %i' %(n_csh))  
        modify_str.append('SAVE SOLID_SOLUTIONS %i' %(n_csh))  
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str)  
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        #print(modify_str)
        #print(output)
        comp = {}        
        comp['CSHQ_TobH_m'] = tobh
        comp['CSHQ_TobD_m'] = tobd
        comp['CSHQ_JenH_m'] = jenh
        comp['CSHQ_JenD_m'] = jend
        comp['sio2am_m'] = sio2
        comp['C'] = output[1][6] 
        comp['Ca'] = output[1][7]
        comp['Si'] = output[1][8]
        comp['H'] = (output[1][9] - self.phrqc.H_norm)
        comp['O'] = (output[1][10] - self.phrqc.O_norm)
        result[str(n_csh)] = comp        
        #print(comp)
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
                  self.solid.CSHQ_TobH.c * self.solid.CSHQ_TobH.mvol + \
                  self.solid.sio2am.c * self.solid.sio2am.mvol
        return CSH_vol

       
    def free_volume(self):
        v = self.solid.voxel_vol - self.solid.portlandite.vol - \
            self.solid.csh_vol - self.solid.calcite.vol
        return v

    
    def csh_conc(self):    
        csh = self.solid.CSHQ_TobH.c+ self.solid.CSHQ_TobD.c +\
                self.solid.CSHQ_JenH.c + self.solid.CSHQ_JenD.c +\
                self.solid.sio2am.c
        return csh 
                
#%% Additional funcrions              
class SettingsCSHQ(Settings):
    
    @staticmethod
    def set_mvols(scale=50):
        '''
        idx 0 - portlandite,
        idx 1 - calcite
        '''
        mvolCH = df.MOLAR_VOLUME['CH']*scale
        mvolCC = df.MOLAR_VOLUME['CC']*scale
        mvolCSH_TobH = df.MOLAR_VOLUME['CSH_TobH']*scale
        mvolCSH_TobD = df.MOLAR_VOLUME['CSH_TobD']*scale
        mvolCSH_JenH = df.MOLAR_VOLUME['CSH_JenH']*scale
        mvolCSH_JenD = df.MOLAR_VOLUME['CSH_JenD']*scale
        mvolSiO2 = df.MOLAR_VOLUME['SiO2']*scale
        mvol =  [mvolCH, mvolCC, mvolCSH_TobH, mvolCSH_TobD,
                 mvolCSH_JenH, mvolCSH_JenD, mvolSiO2]
        return(mvol)
    
    @staticmethod
    def set_init_pqty(mvol, poros, scale=50): 
        max_pqty = Settings.get_max_pqty(mvol)
        init_conc = [0.0]*len(mvol)
        for i in range(len(mvol)):
            init_conc[i] = (1 - poros[i]) * max_pqty[i]
        return init_conc
        
    @staticmethod
    def get_pqty(init_pqty, domain):
        pqty_CH = init_pqty[0] * (domain.nodetype == ct.Type.MULTILEVEL_CH)
        pqty_CC = init_pqty[1] * np.ones(domain.nodetype.shape)   
        pqty_CSHQ_TobH = init_pqty[2] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty_CSHQ_TobD = init_pqty[3] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty_CSHQ_JenH = init_pqty[4] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty_CSHQ_JenD = init_pqty[5] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty_SiO2 = init_pqty[6] * (domain.nodetype == ct.Type.MULTILEVEL) 
        pqty = [pqty_CH, pqty_CC, 
                pqty_CSHQ_TobH, pqty_CSHQ_TobD, 
                pqty_CSHQ_JenH, pqty_CSHQ_JenD,
                pqty_SiO2]    
        return(pqty)
        
    @staticmethod
    def set_labels(domain):
        slabels = np.zeros(domain.nodetype.shape)
        slabels = 100002 * (domain.nodetype == ct.Type.LIQUID) + \
          100002 * (domain.nodetype == ct.Type.INTERFACE) + \
          100003 * (domain.nodetype == ct.Type.MULTILEVEL_CH) + \
          100004 * (domain.nodetype == ct.Type.MULTILEVEL)  + \
          100005 * (domain.nodetype == ct.Type.SOLID)
        return(slabels)
    
    @staticmethod
    def get_porosity(domain, pqty, mvol):
        #print(pqty)
        #print(mvol)
        poros = np.zeros(domain.nodetype.shape)
        poros = (1- pqty[2]*mvol[2] - pqty[3]*mvol[3] - pqty[4]*mvol[4] - \
            pqty[5]*mvol[5]-pqty[6]*mvol[6])*(domain.nodetype == ct.Type.MULTILEVEL) +\
            (1 - pqty[0]*mvol[0])*(domain.nodetype == ct.Type.MULTILEVEL_CH) +\
            1*(domain.nodetype == ct.Type.LIQUID )+\
            1*(domain.nodetype == ct.Type.INTERFACE )+\
            1*(domain.nodetype == ct.Type.SOLID )
        return(poros)
    
    @staticmethod
    def set_domain_params(D, mvol, pqty, poros, app_tort, slabels, 
                          input_file, database = 'cemdata18.dat' ):
        dp={}
        
        dp['D0'] = D                     
        dp['voxel_vol']=1
        dp['mvol']= mvol        
        dp['poros'] = poros                 
        dp['app_tort'] = app_tort         
        dp['solution_labels'] =slabels
        dp['database'] = database
        
        dp['phrqc_input_file']=input_file#'CH_CC_Ceq.phrq'
        dp['eq_names'] = ['portlandite', 'calcite']
        dp['ss_names']={'Tob_jen_ss':['CSHQ_TobH','CSHQ_TobD','CSHQ_JenH','CSHQ_JenD', 'sio2am']}
        dp['solid_phases']={'portlandite': {'c':pqty[0],'mvol':mvol[0], 'type':'diffusive'},
                            'calcite':     {'c':pqty[1],'mvol':mvol[1], 'type':'diffusive'},
                            'CSHQ_TobH':   {'c':pqty[2],'mvol':mvol[2],'type':'diffusive'},
                            'CSHQ_TobD':   {'c':pqty[3],'mvol':mvol[3],'type':'diffusive'},
                            'CSHQ_JenH':   {'c':pqty[4],'mvol':mvol[4],'type':'diffusive'},
                            'CSHQ_JenD':   {'c':pqty[5],'mvol':mvol[5],'type':'diffusive'},
                            'sio2am':   {'c':pqty[6],'mvol':mvol[6],'type':'diffusive'},}    
        return dp

#%% Results

class ResultsCSHQ(Results):
    
    def init_results(self, nodes=[]):
        '''
        pavg_list - list of average values to save 
        nodes - list of parameters per node to save 
        '''
        results = {}     
        params = [] 
        params += ['portlandite', 'calcite', 'sio2am',
                   'CSHQ_TobD', 'CSHQ_JenD', 'CSHQ_JenH', 'CSHQ_TobH',                   
                   'Ca', 'Si', 'C', 'O', 'H'] # always save these parameters
        params += ['solid_CSHQ', 'solid_Ca', 'solid_Si', 
                   'ratio_CaSi','density_CSHQ']
        results={name: [] for name in params}       
        results['params'] = params
        pavg_list = ['tot_vol', 'mean_D_eff', 'mean_poros', 
                          'dissolution_CH', 'nodes_CH', 'active_nodes', 
                          'mean_pH', 'dCH', 'dt', 'nodes_CSH',
                          'precipitation_CC', 'nodes_CC', 'dCC']  
        results.update({name: [] for name in pavg_list})
        results['pavg_list'] = pavg_list
        if nodes: 
            pointparamslist = []
            pointparamslist += params
            pointparamslist += ['vol', 'poros', 'pH', 'De', 'vol_CH', 'vol_CC',
                                'vol_CSHQ', 'solid_CSHQ'] 
            l = [[m+ ' '+n for n in [str(n) for n in nodes]] for m in [str(m) for m in pointparamslist]]
            par_points = [item for sublist in l for item in sublist] #sum(l, [])
            results.update({name: [] for name in par_points}) 
            results['pointparamslist'] = pointparamslist
        results['nodes'] = nodes    
        results['time'] = []   
        self.results = results
    
    def append_results(self, rtmodel, step = 1e+2):
        if (rtmodel.iters%step == 0):
            self.results['time'].append(rtmodel.time)
            for num, phase in enumerate(rtmodel.solid.diffusive_phase_list, start=1):
                self.results[phase].append(np.sum(rtmodel.solid._diffusive_phaseqty[num-1]))        
            for num, comp in enumerate(rtmodel.fluid.components, start=1):
                self.results[comp].append(np.mean(getattr(rtmodel.fluid, comp)._c*getattr(rtmodel.fluid, comp).poros))  
            #cshq
            self.results['solid_CSHQ'].append(self.__class__.CSHQ_solid(rtmodel))
            self.results['solid_Ca'].append(self.__class__.Ca_solid(rtmodel))
            self.results['solid_Si'].append(self.__class__.Si_solid(rtmodel))
            self.results['ratio_CaSi'].append(self.__class__.CaSi_ratio(rtmodel))
            self.results['density_CSHQ'].append(self.__class__.CSHQ_density(rtmodel))
            #general
            self.results['tot_vol'].append(self.__class__.total_mineral_volume(rtmodel))
            self.results['mean_D_eff'].append(self.__class__.mean_effective_D(rtmodel))
            self.results['mean_poros'].append(self.__class__.mean_porosity(rtmodel))
            self.results['nodes_CH'].append(self.__class__.CH_nodes(rtmodel))
            self.results['nodes_CC'].append(self.__class__.CC_nodes(rtmodel))
            self.results['nodes_CSH'].append(self.__class__.CSH_nodes(rtmodel))
            self.results['active_nodes'].append(self.__class__.active_nodes(rtmodel))
            self.results['dt'].append(self.__class__.get_dt(rtmodel))
            self.results['mean_pH'].append(self.__class__.mean_pH(rtmodel))
            self.results['dissolution_CH'].append(self.__class__.CH_dissolution(rtmodel))
            self.results['precipitation_CC'].append(self.__class__.CC_precipitation(rtmodel))
            if (rtmodel.iters ==0): 
                self.results['dCH'].append(0) 
                self.results['dCC'].append(0) 
            else:
                self.results['dCH'].append(self.__class__.delta_CH(rtmodel))
                self.results['dCC'].append(self.__class__.delta_CC(rtmodel))
            
            # nodes
            if self.results['nodes']: # points
                for p in self.results['nodes']:
                    self.results['portlandite'+' ' + str(p)].append(rtmodel.solid.portlandite.c[p])
                    self.results['Ca'+' ' + str(p)].append(rtmodel.fluid.Ca._c[p]+\
                                 rtmodel.fluid.Ca._ss[p]/rtmodel.fluid.Ca.poros[p])
                    self.results['H'+' ' + str(p)].append(rtmodel.fluid.H._c[p]+\
                                 rtmodel.fluid.H._ss[p]/rtmodel.fluid.H.poros[p])
                    self.results['O'+' ' + str(p)].append(rtmodel.fluid.O._c[p]+\
                                 rtmodel.fluid.O._ss[p]/rtmodel.fluid.O.poros[p])
                    self.results['poros'+' ' + str(p)].append(rtmodel.solid.poros[p])
                    self.results['vol'+' ' + str(p)].append(rtmodel.solid.vol[p])
                    self.results['De'+' ' + str(p)].append(rtmodel.fluid.H.De[p])
                    self.results['vol_CH'+' ' + str(p)].append(rtmodel.solid.portlandite.vol[p])
                    self.results['pH'+' ' + str(p)].append(rtmodel.phrqc.selected_output()['pH'][p])
                    self.results['calcite'+' ' + str(p)].append(rtmodel.solid.calcite.c[p])
                    self.results['C'+' ' + str(p)].append(rtmodel.fluid.C._c[p]+\
                                 rtmodel.fluid.C._ss[p]/rtmodel.fluid.C.poros[p])
                    self.results['vol_CC'+' ' + str(p)].append(rtmodel.solid.calcite.vol[p])
                    self.results['solid_CSHQ'+' ' + str(p)].append(rtmodel.solid.CSHQ_TobD.c[p] + \
                                 rtmodel.solid.CSHQ_TobH.c[p] + rtmodel.solid.CSHQ_JenD.c[p] + \
                                 rtmodel.solid.CSHQ_JenH.c[p] )
                    self.results['vol_CSHQ'+' ' + str(p)].append(rtmodel.solid.csh_vol[p])
                    
    @staticmethod
    def delta_CC(rt):
        '''
        Differentiation of Portladite
        '''
        s = np.sum(rt.phrqc.dphases['calcite'])#*getattr(rt.solid, 'calcite').mvol)
        return(s)
    
    @staticmethod
    def CC_precipitation(rt):
        '''
        N cells that are dissolving
        '''
        s = np.sum(rt.phrqc.dphases['calcite']>0)
        return(s)
         
    @staticmethod
    def CC_nodes(rt):
        '''
        N cells containing portlandite
        '''
        s = np.sum(rt.solid.calcite.c>0)
        return(s)
        
    @staticmethod
    def CSH_nodes(rt):
        '''
        N cells containing portlandite
        '''
        s = np.sum(rt.solid.csh_vol>0)
        return(s)
        
    @staticmethod
    def CSHQ_solid(rt): #TODO check
        csh = np.sum(rt.solid.CSHQ_TobD.c + rt.solid.CSHQ_TobH.c + 
                     rt.solid.CSHQ_JenH.c + 1.5*rt.solid.CSHQ_JenD.c)
        return(csh)
        
    @staticmethod
    def Ca_solid(rt): #TODO check
        ca = np.sum(0.8333333*rt.solid.CSHQ_TobD.c[:,:] + 0.6666667*rt.solid.CSHQ_TobH.c[:,:] + 
                    1.3333333*rt.solid.CSHQ_JenH.c[:,:] + 1.5*rt.solid.CSHQ_JenD.c[:,:])
        return(ca)
        
    @staticmethod
    def Si_solid(rt): #TODO check
        si = np.sum(0.6666667*rt.solid.CSHQ_TobD.c[:,:] + 1.0*rt.solid.CSHQ_TobH.c[:,:] + 
                    1.0*rt.solid.CSHQ_JenH.c[:,:] + 0.6666667*rt.solid.CSHQ_JenD.c[:,:])
        return (si)
    
    @staticmethod
    def CaSi_ratio(rt):
        r = ResultsCSHQ.Ca_solid(rt)/ResultsCSHQ.Si_solid(rt)
        return(r)

    @staticmethod
    def CSHQ_density(rt):
        m_ca = 56.0774 #g/mol
        m_si = 60.08 #g/mol
        m_h2o = 18.01528 #g/mol
        h2o= np.sum(1.8333333333*rt.solid.CSHQ_TobD.c[:,:] + 1.5*rt.solid.CSHQ_TobH.c[:,:] + 
            2.1666666667*rt.solid.CSHQ_JenH.c[:,:] + 2.5*rt.solid.CSHQ_JenD.c[:,:])
        d = h2o*m_h2o + ResultsCSHQ.Ca_solid(rt)*m_ca + ResultsCSHQ.Si_solid(rt)*m_si
        return(d)
    #%% PRINT and PLOT
    
    @staticmethod
    def print_profiles(rt):    
        Results.print_profiles(rt)
        print('C +ss %s' %str(np.array(rt.fluid.C.c[1,:]) + \
               np.array(rt.fluid.C._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
        print('Si +ss %s' %str(np.array(rt.fluid.Si.c[1,:]) + \
               np.array(rt.fluid.Si._ss[1,:])/np.array(rt.phrqc.poros[1,:])))
        print('CC %s' %str(np.array(rt.solid.calcite.c[1,:])))
        print('dCC %s' %str(np.array(rt.phrqc.dphases['calcite'][1,:])))
        print('CSHQ TobD %s' %str(np.array(rt.solid.CSHQ_TobD.c[1,:])))
        #print('PHRQC CSHQ TobD %s' %str(np.array(rt.phrqc.selected_output()['CSHQ_TobD'][1,:])))
        print('CSHQ TobH %s' %str(np.array(rt.solid.CSHQ_TobH.c[1,:])))
        #print('PHRQC CSHQ TobH %s' %str(np.array(rt.phrqc.selected_output()['CSHQ_TobH'][1,:])))
        print('CSHQ JenD %s' %str(np.array(rt.solid.CSHQ_JenD.c[1,:])))
        print('CSHQ JenH %s' %str(np.array(rt.solid.CSHQ_JenH.c[1,:])))
        print('SiO2 am %s' %str(np.array(rt.solid.sio2am.c[1,:])))
        
    def get_titles(self):
        title = {'Ca': 'Calcium',
                 'C': 'Carbon',
                 'O': 'Oxygen',
                 'H': 'Hydrogen',
                 'portlandite':'Portlandite',
                 'calcite': 'Calcite',
                 'De':  'Effective diffusivity',
                 'poros': 'Porosity',
                 'vol_CH': 'Portlandite volume',
                 'vol_CC': 'Calcite volume',
                 'vol_CSHQ': 'CSH volume',
                 'solid_CSHQ': 'CSH',
                 'vol':'Mineral volume',
                 'dt':'Time step',
                 'pH': 'pH',
                 'active_nodes':'Number of active nodes (PHREEQC)',
                 'mean_poros':'Average porosity',
                 'mean_D_eff':'Average effective diffusivity',
                 'mean_pH': 'Average pH ',
                 'dCH': 'Portlandite differentiation',
                 'dissolution_CH':'Number of dissolving nodes',
                 'nodes_CH':'Number of Portlandite nodes',
                 'tot_vol':'Total mineral volume',
                 'dCC': 'Calcite differentiation',
                 'precipitation_CC':'Number of precipitating nodes',
                 'nodes_CC':'Number of Calcite nodes',
                 'Si': 'Silicium',
                 'CSHQ_TobH': 'Tobermorite H',
                 'CSHQ_TobD': 'Tobermorite D',
                 'CSHQ_JenH': 'Jennite H',
                 'CSHQ_JenD': 'Jennite D',
                 'sio2am': 'SiO2',
                 'nodes_CSH':'Number of CSH nodes',
                 'solid_Ca': 'Ca in CSH',
                 'solid_Si': 'Si in CSH',
                 'ratio_CaSi': 'Ca/Si',
                 'density_CSHQ': 'CSH density',                 
                 }
        return(title)
         
    def get_ylabs(self):
        ylab = {'Ca':'Ca (*1e-12) [mol]',
                'C':'C (*1e-12) [mol]',
                'O':'O (*1e-12) [mol]',
                'H':'H (*1e-12) [mol]',
                'portlandite': 'CH (*1e-12) [mol]',
                'calcite': 'CC (*1e-12) [mol]',
                'De':  'D_eff [m^2/s]',
                'poros': 'Porosity [-]',
                'vol_CH': 'Portlandite volume [um^3]',
                'vol_CC': 'Calcite volume [um^3]',
                'vol_CSHQ': 'CSH volume [um^3]',
                'solid_CSHQ': 'CSH (*1e-12) [mol]',
                'vol':'Mineral volume [um^3]',
                'dt':'dt [s]',
                'pH': 'pH [-]',                
                'active_nodes':'Nodes [-]',
                'mean_poros':'Average porosity [-]',
                'mean_D_eff':'Average D_eff [m^2/s]',
                'mean_pH': 'Average  pH [-]',              
                'dCH':'dCH [mol]',
                'dissolution_CH':'Nodes [-]',
                'nodes_CH':'Nodes [-]',
                'tot_vol': 'Total mineral volume [um^3]',
                'dCC':'dCC [mol]',
                'precipitation_CC':'Nodes [-]',
                'nodes_CC':'Nodes [-]',
                 'Si': 'Si (*1e-12) [mol]',
                 'CSHQ_TobH': 'TobH (*1e-12) [mol]',
                 'CSHQ_TobD': 'TobD (*1e-12) [mol]',
                 'CSHQ_JenH': 'JenH (*1e-12) [mol]',
                 'CSHQ_JenD': 'JenD (*1e-12) [mol]',
                 'sio2am': 'SiO2 (*1e-12) [mol]',
                 'nodes_CSH':'Nodes [-]',
                 'solid_Ca': 'Ca in CSH (*1e-12) [mol]',
                 'solid_Si': 'Si in CSH (*1e-12) [mol]',
                 'ratio_CaSi': 'Ca/Si [-]',
                 'density_CSHQ': 'CSH density [g/l]',      
                }
        return(ylab)
    
    def plot_fields(self, rt, names={}, fsize = (8,4)):
        #nl = {name: limit}
        
        def plot(field, title, ylab, limit=[], size=fsize):
            plt.figure() 
            if not limit:
                plt.imshow(field)
            else:
                plt.imshow(field, vmin = limit[0], vmax = limit[1])
            
            clb = plt.colorbar()
            clb.ax.get_yaxis().labelpad = 15
            clb.ax.set_ylabel(ylab, rotation=270)
            
            plt.title(title)
            plt.show()
        
        fields = {'portlandite': rt.solid.portlandite.c[1:-1,1:-1],
                  'calcite': rt.solid.calcite.c[1:-1,1:-1],
                  'Ca': rt.fluid.Ca.c[1:-1,0:-1],
                  'C': rt.fluid.C.c[1:-1,0:-1],
                  'H':  rt.fluid.H.c[1:-1,0:-1],
                  'O':  rt.fluid.O.c[1:-1,0:-1],
                  'poros': rt.solid.poros[1:-1,0:-1]}
            
        if not names:
            names = fields.keys()
        
        ylab = self.get_ylabs()
        title = self.get_titles()
        for k in names: 
            plot(fields[k], title[k], ylab[k])
        
    #%% SAVE
    def save_vti(rt, t, path, name):
        
        nx = rt.fluid.Ca.nx -1
        ny = rt.fluid.Ca.ny -1     
        filename = path + str(name) +'_all_' + str(t) 
        
        cCa = rt.fluid.Ca.c[1:ny,1:nx,np.newaxis]
        cSi = rt.fluid.Si.c[1:ny,1:nx,np.newaxis]
        cH = rt.fluid.H.c[1:ny,1:nx,np.newaxis]
        cO = rt.fluid.O.c[1:ny,1:nx,np.newaxis]
        cC = rt.fluid.C.c[1:ny,1:nx,np.newaxis]
        cCC = rt.solid.calcite.c[1:ny,1:nx,np.newaxis]
        cCH = rt.solid.portlandite.c[1:ny,1:nx,np.newaxis]
        cSiO2 = rt.solid.sio2am.c[1:ny,1:nx,np.newaxis]
        cTobH = rt.solid.CSHQ_TobH.c[1:ny,1:nx,np.newaxis]
        cTobD = rt.solid.CSHQ_TobD.c[1:ny,1:nx,np.newaxis]
        cJenH = rt.solid.CSHQ_JenH.c[1:ny,1:nx,np.newaxis]
        cJenD = rt.solid.CSHQ_JenD.c[1:ny,1:nx,np.newaxis]
        porosity = rt.solid._poros[1:ny,1:nx,np.newaxis]
        
        imageToVTK(filename, cellData = {"Ca":cCa,
                                         "Si":cSi,
                                         "C":cC,
                                         "O":cO,
                                         "H":cH,
                                         "Calcite":cCC,
                                         "Portlandite":cCH,
                                         "SiO2":cSiO2,
                                         "TobH":cTobH,
                                         "TobD":cTobD,
                                         "JenH":cJenH,
                                         "JenD":cJenD,
                                         "porosity":porosity,
                                         })