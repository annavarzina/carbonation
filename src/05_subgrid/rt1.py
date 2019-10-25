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
        self.dx = self.fluid.H.dx  #pcs settings      
        self.ptype = 'CSH' if hasattr(self.fluid, 'Si') else 'CH'
        self.set_volume()
        self.set_porosity()        
        self.nodetype = deepcopy(domain.nodetype)  
        self.apply_settings(settings)  
        c=deepcopy( self.fluid.get_attr('c'))
        self.set_phrqc(c, self.nodetype)        
        self.set_bc(settings)
        self.solid.phases = self.update_phases()
        self.update_nodetype()
        
    
    def advance(self):
        #print(self.iters)
        
        self.update_target_SI()
        if(self.settings['diffusivity']['type']=='fixed'): 
            self.update_diffusivity()
        elif(self.settings['diffusivity']['type']=='cc_archie'): 
            self.update_diffusivity_ch()
        elif(self.settings['diffusivity']['type']=='cc_archie_ch_kin'): 
            self.update_diffusivity_ch_kin()
        elif(self.settings['diffusivity']['type']=='mixed'): 
            self.update_diffusivity_mixed()
        self.fluid.call('advance')  
        self.phrqc.is_calc = self.solid.calcite.c>0        
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):
            self.fluid.call('update_transport_params',self.solid.poros,
                            self.solid.app_tort,self.auto_time_step)
            self.phrqc.poros=deepcopy(self.solid.poros)      
        #self.correct()                
        self.fluid.call('_set_relaxation_params')  
        if(self.phrqc.precipitation == 'interface' and self.phrqc.active != 'interface'):
            self.phrqc.nodetype = deepcopy(self.nodetype) 
            
        
        
        c=deepcopy( self.fluid.get_attr('c'))
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1])
        self.phrqc.modify_solid_phases(phaseqty, c, self.solid.nodetype) 
        ss=self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        new_c = self.update_neighbour_solution(9,10, self.solid.portlandite.c[1,2], self.solid.calcite.c[1,2], self.solid.poros[1,1])
        #print('new c %s' %str(new_c))
        #print(ss)
        #self.solid.portlandite.c[1,2] =new_c['portlandite']
        self.solid.calcite.c[1,1] =new_c['calcite']
        '''
        print('Ca %s' %str(np.array(self.fluid.Ca.c[1,:]) ))
        print('CH %s' %str(np.array(self.solid.portlandite.c[1,:])))
        print('CC %s' %str(np.array(self.solid.calcite.c[1,:])))
        print('dCH %s' %str(np.array(self.phrqc.dphases['portlandite'][1,:])))
        print('dCC %s' %str(np.array(self.phrqc.dphases['calcite'][1,:])))
        '''
        self.update_solid_params() # or after if?    
        #self.fluid.call('update_transport_params',self.solid.poros,
        #                    self.solid.app_tort,self.auto_time_step)
        #self.phrqc.poros=deepcopy(self.solid.poros)      
        ss = self.modify_all(9,10, new_c, c, self.solid.nodetype,phaseqty)
        #print(ss)
        '''
        print('Ca %s' %str(np.array(self.fluid.Ca.c[1,:]) ))
        print('C %s' %str(np.array(self.fluid.C.c[1,:]) ))
        print('Ca ss %s' %str(np.array(self.fluid.Ca._ss[1,:])))
        print('C ss %s' %str(np.array(self.fluid.C._ss[1,:])))
        '''
        
        #print(self.phrqc.selected_output()['portlandite'])
        ss['Ca'] = ss['Ca']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1) 
        pqty=self.solid.update(self.phrqc.dphases)  
        self.fluid.set_attr('ss',ss)
        '''
        print('Ca %s' %str(np.array(self.fluid.Ca.c[1,:])))
        print('Ca %s' %str(np.array(self.fluid.Ca.c[1,:])+np.array(self.fluid.Ca._ss[1,:])/np.array(self.solid.poros[1,:]) ))
        print('CH %s' %str(np.array(self.solid.portlandite.c[1,:])))
        print('CC %s' %str(np.array(self.solid.calcite.c[1,:])))
        print('dCH %s' %str(np.array(self.phrqc.dphases['portlandite'][1,:])))
        print('dCC %s' %str(np.array(self.phrqc.dphases['calcite'][1,:])))
        '''
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)
        self.update_solid_params() # or after if?        
        if(self.settings['velocity']):
            self.update_velocity()
        self.solid.phases = self.update_phases()
        if(self.phrqc.precipitation == 'interface'):
            self.update_nodetype()
        
    def update_neighbour_solution(self, n_int, n_ch, p_ch, p_cc, poros_i):
        n = 123456
        modify_str = []
        modify_str.append("EQUILIBRIUM_PHASES %i" % n)
        modify_str.append("Portlandite 0 %.20e dissolve only" %(p_ch))   
        modify_str.append("Calcite 0 %.20e" %(p_cc))   # precipitate only
        modify_str.append("END") 
        modify_str.append('MIX %i' % n) 
        modify_str.append('%i 1' % n_int)  
        modify_str.append('SAVE solution %i' % n)  
        modify_str.append("END") 
        modify_str.append('USE solution %i' % n)  
        modify_str.append('USE equilibrium_phase %i' % n)
        modify_str.append('SAVE solution %i' % n)   
        modify_str.append("END") 
        modify_str.append("END") 
        modify_str ='\n'.join(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str) 
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        #print(modify_str)
        #print(output)
        port = output[2][10]
        calc = output[2][12]
        modify_str = [] 
        modify_str.append("EQUILIBRIUM_PHASES %i" % n_ch)
        modify_str.append("Portlandite 0 %.20e" %(port))   
        modify_str.append("Calcite 0 %.20e" %(calc)) # precipitate only
        modify_str.append("END") 
        modify_str.append('USE equilibrium_phase %i' % n_int)        
        modify_str.append('MIX %i' % n)      
        modify_str.append('%i 1' % n)     
        #modify_str.append('%i %.20e' %( n, (poros_i) ))
        modify_str.append('%i 0' % n_int)  
        modify_str.append('SAVE solution %i' % n_int)  
        modify_str.append('SAVE equilibrium_phase %i' % n_int) 
        modify_str.append("END") 
               
        modify_str ='\n'.join(modify_str)
        #print(modify_str)
        self.phrqc.IPhreeqc.RunString(modify_str)  
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        #print(output)
        
        comp = {}
        comp['portlandite'] = port
        comp['calcite'] = output[1][12]
        comp['C'] = output[1][6]
        comp['Ca'] = output[1][7]
        comp['H'] = output[1][8]
        comp['O'] = output[1][9]
        
        return(comp)
    def modify_all(self, n_i, n_m, new_c, c, nodetype, phaseqty):
        '''
        phaseqty = self.phrqc.flatten_dict(phaseqty)
        modifystr = []
        for cell in range(self.phrqc.startcell,self.phrqc.stopcell+1,1):
            modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)
            if cell == n_m:                
                modifystr.append("\t -component portlandite")
                modifystr.append("\t\t%s\t%.20e" %('-moles', new_c['portlandite']))
                modifystr.append("\t -component calcite")
                modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty['calcite'][cell-1]))
            else:
                for key in phaseqty.keys():
                    modifystr.append("\t -component %s" %(key))
                    modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1]))
                    #modifystr.append("\t\t%s\t%.20e" %('-si', 2))
        modifystr.append("end")       
        modifystr ='\n'.join(modifystr)
        self.phrqc.IPhreeqc.RunString(modifystr)
        #print(modifystr)
        '''
        
        moles_Ca = new_c['Ca']#*self.phrqc.poros.flatten(order='C')[n_i]
        moles_C = new_c['C']#*self.phrqc.poros.flatten(order='C')[n_i]
        moles_H = new_c['H']#(self.solid.vol.flatten(order='C')[n_i])#-self.phrqc.H_norm
        moles_O = new_c['O']#(self.solid.vol.flatten(order='C')[n_i])#-self.phrqc.O_norm
        
        active_nodes = self.phrqc.active_nodes(c,nodetype)
        c_trans=deepcopy(c)
        moles = self.phrqc.flatten_dict(self.phrqc.to_moles(c))
        modify_str=[]
        runcells=[]
        runcell_str=[]
        for i,cell in enumerate(range(self.phrqc.startcell,self.phrqc.stopcell+1,1)):
            if active_nodes[i]:
                if cell == n_i:
                    runcells.append(str(cell))
                    modify_str.append('SOLUTION_MODIFY %i' % n_i)
                    if moles_H <=0: moles_H= 1e-30
                    modify_str.append('\t%s\t%.20e' % ('total_h', moles_H))
                    if moles_O <=0: moles_O= 1e-30
                    modify_str.append('\t%s\t%.20e' % ('total_o', moles_O))
                    modify_str.append('\t-totals')
                    if moles_Ca <=0: moles_Ca= 1e-30
                    modify_str.append('\t\t%s\t%.20e' % ('Ca', moles_Ca))
                    modify_str.append('end')   
                    if moles_C <=0: moles_C= 1e-30
                    modify_str.append('\t\t%s\t%.20e' % ('C', moles_C))
                    modify_str.append('end')   
                else:#elif cell == n_m:
                    
                    runcells.append(str(cell))
                    modify_str.append('SOLUTION_MODIFY %i' % cell)
                    if 'H' in moles:
                        c = moles['H'][i]
                        if c <=0: c= 1e-30
                        modify_str.append('\t%s\t%.20e' % ('total_h', c))
                    if 'O' in moles:
                        c = moles['O'][i]
                        if c <=0: c= 1e-30
                        modify_str.append('\t%s\t%.20e' % ('total_o', c))
                    modify_str.append('\t-totals')
                    for name,val in moles.iteritems():
                        if (name != 'H') and (name!='O'):
                            c = val[i]
                            if c <=0: c= 1e-30
                            modify_str.append('\t\t%s\t%.20e' % (name, c))
        modify_str ='\n'.join(modify_str)
        runcell_str.append('RUN_CELLS')
        runcell_str.append('\t-cells %s'%'\n\t\t'.join(runcells))
        runcell_str.append('\t-start_time %s'%self.time)
        runcell_str.append('\t-time_step %s'%self.dt)
        runcell_str.append('end')
        runcell_str='\n'.join(runcell_str)
        '''
        import json
        with open('modity_solution_new.txt', 'w') as file:
            file.write(json.dumps(modify_str))
            file.write('\n\n')
            file.write(json.dumps(runcell_str))
            file.close()
        '''
        #print(modify_str)
        self.phrqc.IPhreeqc.AccumulateLine(modify_str)
        self.phrqc.IPhreeqc.AccumulateLine(runcell_str)
        self.phrqc.IPhreeqc.RunAccumulated()
        
        output=self.phrqc.IPhreeqc.GetSelectedOutputArray()
        selected_output={}
        if len(output) > 0:
            header = output[0]
            for head in header:
                 selected_output[head] = []
            for row in output[1:]:
                for (i, head) in enumerate(header):
                    results = row[i]
                    if head == 'H':
                        results -= self.phrqc.H_norm
                    elif head == 'O':
                        results -= self.phrqc.O_norm
                    selected_output[head].append(results)
        merged_output = self.phrqc.list_dict(self.phrqc.flatten_dict(self.phrqc.prev_selected_output))
        for name,val in selected_output.iteritems():
            for i,cell in enumerate(selected_output[u'soln']): 
                merged_output[name][cell-1]=val[i]             
        selected_output = merged_output
        selected_output = self.phrqc.ndarray_dict(selected_output)
        selected_output = self.phrqc.reshape_dict(selected_output,self.phrqc.array_shape)
        self.phrqc._selected_output = selected_output
        self.phrqc.phrqc_flags['update_output'] = False
        self.phrqc.poros = self.phrqc.selected_output()['poros']#self.solid.poros#
        c_current = self.phrqc.component_conc
        ss={}
        #print('c_trans %s' %str(c_trans))
        #print('c_curr %s'%str(c_current))
        for name in self.phrqc.components:
            ss[name] = (c_current[name]-c_trans[name])/self.phrqc.dt
            ss[name] = ss[name] *(active_nodes.reshape(self.phrqc.array_shape)>0)
            ss[name] *= self.phrqc.poros
        '''
        print('ss %s'%str(ss))
        print('dt %s'%str(self.phrqc.dt))
        print('T %s'%str(self.fluid.Ca.convfactors['T']))
        print('poros s %s' %str(self.solid.poros))
        print('poros p %s' %str(self.phrqc.poros))
        print('poros so %s' %str(self.phrqc.selected_output()['poros']))
        '''
        
        return(ss) 
        
    #%% SETTINGS
    def set_volume(self):
        self.solid.vol = np.zeros(self.solid.shape)
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            val = getattr(self.solid, phase)
            self.solid.vol += val.c * self.solid.mvol[num-1] 
            
    def set_porosity(self):        
        self.solid.poros=1.- self.solid.vol/self.solid.voxel_vol 
        self.solid.app_tort = 1. * self.solid.poros ** (1./3.)
        if np.any(self.solid.poros<=1e-10):
            sys.exit('Negative or zero porosity')
        
    def set_boundary_cells(self):
        nodes = np.zeros(np.shape(self.nodetype))
        nodes[0,:]  = ct.Type.SOLID
        nodes[-1,:] = ct.Type.SOLID
        nodes[:,0]  = ct.Type.SOLID
        nodes[:,-1] = ct.Type.SOLID
        return nodes        
     
    def set_bc(self, settings):
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
    
    def set_phrqc(self, c, nodetype):
        self.phrqc.boundcells = self.set_boundary_cells()
        self.phrqc.init_port = self.solid.portlandite.c>0
        self.phrqc.is_calc = self.solid.calcite.c>0
        self.phrqc._target_SI = np.zeros(self.solid.shape)   
        #self.phrqc.active = 'all'
        #self.phrqc.precipitation = 'all'
        #self.phrqc.pcs = False
        self.phrqc.nodetype = deepcopy(self.nodetype) 
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty, c, nodetype) 
        
    def apply_settings(self, settings):
        self.settings = settings
        self.phrqc.pinput = settings['bc']
        
        self.phrqc.active = settings['active']
        self.phrqc.precipitation = settings['precipitation']
        self.phrqc.pcs = settings['pcs']['pcs']
        if(settings['active'] == 'all'):
            self.phrqc.phrqc_flags['smart_run'] = False
            self.phrqc.phrqc_flags['only_interface'] = False
            self.phrqc.nodetype = deepcopy(self.nodetype)
        elif(settings['active'] == 'interface'):        
            self.phrqc.phrqc_flags['smart_run'] = False
            self.phrqc.phrqc_flags['only_interface'] = True
        elif(settings['active'] == 'smart'):
            self.phrqc.phrqc_flags['smart_run'] = True
            self.phrqc.phrqc_flags['only_interface'] = False
        else:
            sys.exit()
        if ('pco2' in settings) and (settings['bc']['type']=='pco2'):            
            self.phrqc.ppt = settings['pco2']
        else:         
            self.phrqc.ppt = False
                
    #%% UPDATES
    def update_solid_params(self):
        '''
        Calculates porosity, apparent tortuosity,
        volume and change in volume last time step
        '''
        self.solid.prev_vol = deepcopy(self.solid.vol)
        self.set_volume()
        self.set_porosity()
        if(self.settings['velocity'] == True):
            self.solid.dvol = self.solid.vol-self.solid.prev_vol
        
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
            #is_liquid =  (~is_port)
            #if(self.iters>1):
            #    is_critical = (self.solid.pore_size <= self.solid.threshold_pore_size) & is_calc
            #    #is_critical = (self.solid.target_SI >= self.settings['si_params']['threshold_SI']) & is_calc
            #    is_liquid =  (~is_critical)&(~is_port)&(~is_solid)&(prev_nodetype==-1)
            #not_critical = np.logical_not(is_critical)
            #is_interface = (not_critical) & (is_port|(prev_nodetype==-2)| (prev_nodetype==-5))
            #self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
            #    ct.Type.INTERFACE * is_interface + \
            #    ct.Type.MULTILEVEL * is_critical +\
            #    ct.Type.SOLID * is_solid #ct.Type.INTERFACE * ((is_calc) &(~is_critical)) + 
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
        D_CC = self.settings['diffusivity']['D_CC']
        D_CH = self.settings['diffusivity']['D_CH']
        Dref = self.fluid.H.Deref 
        is_port = self.solid.portlandite.c >0
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
        is_liquid = np.logical_and(~is_port, ~is_calc)
        De = D_CC *is_calc + D_CH *is_port+ Dref * is_liquid 
        Dnew_lb = De/self.solid.poros/self.solid.app_tort
        #print(Dnew_lb)
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)

    def update_diffusivity_ch(self):
        D_CH = self.settings['diffusivity']['D_CH']
        Dref = self.fluid.H.Deref 
        is_port = self.solid.portlandite.c >0
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
        is_liquid = np.logical_and(~is_port, ~is_calc)
        De = D_CH *is_port+ Dref * is_liquid 
        Dnew_lb = De/self.solid.poros/self.solid.app_tort
        Dnew_lb = Dref*is_calc + Dnew_lb*(~is_calc)
        #print(Dnew_lb)
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)  
        
    def update_diffusivity_ch_kin(self):
        D_CH_Ca = self.settings['diffusivity']['D_CH_Ca']
        D_CH = self.settings['diffusivity']['D_CH']
        Dref = self.fluid.H.Deref 
        is_port = self.solid.portlandite.c >0
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
        is_liquid = np.logical_and(~is_port, ~is_calc)
        De = D_CH *is_port+ Dref * is_liquid 
        Dnew_lb = De/self.solid.poros/self.solid.app_tort
        Dnew_lb = Dref*is_calc + Dnew_lb*(~is_calc)
        #print(Dnew_lb)
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)  
        
        De = D_CH_Ca *is_port+ Dref * is_liquid 
        Dnew_lb = De/self.solid.poros/self.solid.app_tort
        Dnew_lb = Dref*is_calc + Dnew_lb*(~is_calc)
        
        self.fluid.Ca.D0 = Dnew_lb
        self.fluid.Ca.Deref=np.max(Dnew_lb)
        self.fluid.Ca.Dr=Dnew_lb
        
    def update_diffusivity_mixed(self):
        D_CH_Ca = self.settings['diffusivity']['D_CH_Ca']
        D_CH = self.settings['diffusivity']['D_CH']
        Dref = self.fluid.H.Deref
        
        is_port = self.solid.portlandite.c >0
        is_calc = np.logical_and(self.solid.calcite.c >0,~is_port)
        is_liquid = np.logical_and(~is_port, ~is_calc)
        
        mCH = self.solid.portlandite.c * self.solid.portlandite.mvol
        D_mixed = np.nan_to_num(1./((1-mCH)/Dref/self.solid.poros/self.solid.app_tort + mCH/D_CH), Dref)
        De = D_mixed *(~is_liquid)+ Dref * is_liquid
        Dnew_lb = De/self.solid.poros/self.solid.app_tort
        #print(Dnew_lb)
        self.fluid.set_attr('D0',Dnew_lb,component_dict=False)
        self.fluid.set_attr('Deref',np.max(Dnew_lb),component_dict=False)
        self.fluid.set_attr('Dr',Dnew_lb,component_dict=False)  
        
        D_mixed = np.nan_to_num(1./((1-mCH)/Dref/self.solid.poros/self.solid.app_tort + mCH/D_CH_Ca), Dref)
        De = D_mixed *(~is_liquid)+ Dref * is_liquid
        Dnew_lb = De/self.solid.poros/self.solid.app_tort
        
        self.fluid.Ca.D0 = Dnew_lb
        self.fluid.Ca.Deref=np.max(Dnew_lb)
        self.fluid.Ca.Dr=Dnew_lb
        
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
    
    def update_target_SI(self):
        if (self.ptype == 'CH'):  
            self.solid.vol_ch = self.volume_CH()
            self.solid.vol_cc = self.volume_CC()          
            if (self.settings['pcs']['pores']=='cylinder'):
                #TODO check this case
                if (self.iters<=1):
                    self.solid.threshold_pore_size = self.settings['pcs']['crystal_size']
                    self.solid.pore_length = self.settings['pcs']['crystal_size']  
                self.solid.pore_amount = self.pore_amount()  
                self.solid.free_vol_cc = (self.solid.voxel_vol-self.solid.vol_ch)* self.dx**3#self.free_volume()  
                self.solid.pore_size = self.pore_size()
                self.solid.pore_volume_cc = self.pore_volume_CC() 
            elif(self.settings['pcs']['pores']=='block'):                
                if (self.iters<1):      
                    self.solid.threshold_crystal = self.settings['pcs']['crystal_size'] #default crystal size
                    self.solid.threshold_pore_size = self.settings['pcs']['pore_size']
                    self.solid.threshold_distance = self.solid.threshold_pore_size*2
                    self.solid.block_size = self.block_size()
                self.solid.free_vol_cc = (self.solid.voxel_vol-self.solid.vol_ch)* self.dx**3#self.solid._poros * self.dx**3#self.get_pore_volume(self)
                self.solid.ncrystals = self.block_amount()                 
                self.solid.crystal_size = self.crystal_size()
                self.solid.distance = self.distance()
                self.solid.pore_size =  self.pore_size()
            else:
                pass        
            self.solid.taget_SI = self.target_SI()
            self.phrqc._target_SI = self.solid.taget_SI
        elif (self.ptype == 'CSH'): 
            pass        
    def correct(self):
        c = deepcopy(self.fluid.get_attr('_c'))
        f = deepcopy(self.fluid.get_attr('_f'))
        flag = False
        for comp in self.fluid.components:
            n = c[comp] < 0
            if n.any():
                flag = True            
                f[comp][n] = f[comp][n] - f[comp][n]*(f[comp][n]<0)
        if flag: 
            self.fluid.set_attr('_f',f)
            self.fluid.set_attr('f',f)
            self.fluid.call('compute_macro_var')      
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
        if(self.settings['pcs']['pores'] == 'cylinder'):
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
        n = self.settings['pcs']['pore_density'] * vol #* np.ones(np.shape(pore_vol))
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
        if(self.settings['pcs']['pores'] == 'cylinder'):
            pore_size = self.get_cylinder_radius(self.solid.poros, 
                                                 self.solid.pore_amount, #self.settings['pcs']['pore_density']
                                                 self.solid.pore_length, 
                                                 self.dx)
            #pore_size[pore_size>self.solid.threshold_pore_size ] = self.solid.threshold_pore_size
        else:
            pore_size=self.solid.distance/2.
        return pore_size
    
    def threshold_pore_size(self):
        si = self.settings['pcs']['threshold_SI']
        omega = 10.**si
        ps = (-1) * self.settings['si_params']['mvol'] *\
            np.cos(np.pi*self.settings['si_params']['angle'])*2*\
            self.settings['pcs']['iene'] /\
            self.settings['pcs']['R'] /\
            self.settings['pcs']['T'] / np.log(omega)
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
        is_port = (self.solid.vol_ch >0)
        not_critical = (pore_size >= self.solid.threshold_pore_size) #rt.solid.radius_max 
        si[np.logical_and(is_port, not_critical)]= 0
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
        params = self.settings['pcs']
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
    #%% MISCELANEOUS
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
        poros = self.solid.poros
        f= self.fluid.get_attr('f')
        for key, value in f.iteritems():
            f[key] = get_feq(c[key], poros)
        self.fluid.set_attr('f',f)
        
        
 
            