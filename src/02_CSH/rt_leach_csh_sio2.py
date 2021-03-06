# -*- coding: utf-8 -*-
import sys, os
src_dir = os.path.dirname(__file__)
sys.path.append(src_dir)

import numpy as np
from copy import deepcopy
import yantra
#from rt import CarbonationRT
import cell_type as ct # change the path to cell_type file
from rt_leach import LeachingRT
from phrqc_input import PhreeqcInput

class PhreeqcInputCSH(PhreeqcInput):
    def make_phrqc_input(self):
        self.phrqc_phases(self.c['csh'])
        self.phrqc_boundary_voxel_CSH(ca = self.c['ca_bc'], si = self.c['si_bc'])
        self.phrqc_liquid_voxel_CSH(ca = self.c['ca_liq'], si = self.c['si_liq'])
        self.phrqc_multilevel_voxel_CSH(ca = self.c['ca_mlvl'], si = self.c['si_mlvl'])        
        self.phrqc_solid_voxel()
        #print(self.phrqc_input)
        return(self.phrqc_input)
        
    def phrqc_phases(self, csh):        
        phrqc_input = [] 
        s = csh['stochiometry']
        oh = 2*s['Ca'] - s['Si']
        h2o =  s['H2O'] - s['Ca']- s['Si']
        sign = '+'
        if h2o < 0: 
            sign = '-'
            h2o *= -1
        phrqc_input.append('PHASES')  
        phrqc_input.append(csh['name'])  
        phrqc_input.append('\t(CaO)' + str(s['Ca']) +'(SiO2)'+ str(s['Si']) + \
                           '(H2O)' + str(s['H2O']) + ' = ' + str(s['Ca']) + \
                           'Ca++ + ' + str(s['Si']) + ' SiO(OH)3- + ' + \
                           str(oh) + 'OH- ' + sign + ' ' + str(h2o) + ' H2O') 
        phrqc_input.append('\tlog_K ' + str(csh['log_k']) )
        phrqc_input.append('knobs') 
        phrqc_input.append('\t-iterations 8000\n') 
        
        self.phrqc_input +=  phrqc_input
        
    def phrqc_boundary_voxel_CSH(self, ca, si):
        phrqc_input = [] 
        phrqc_input.append('#boundary_solution')    
        phrqc_input.append('SOLUTION\t100001')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
            phrqc_input.append('\tSi\t' + str(si['value']) + '\n')
        else:
            pass
        phrqc_input.append('EQUILIBRIUM_PHASES\t100001\n')
        #phrqc_input.append(self.c['csh']['name'] + '\t0\t0\n')
        self.phrqc_input +=  phrqc_input
        
    def phrqc_liquid_voxel_CSH(self, ca, si):
        phrqc_input = [] 
        phrqc_input.append('#solution_liquid')    
        phrqc_input.append('SOLUTION\t100002')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
            phrqc_input.append('\tSi\t' + str(si['value'])+ '\n')
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
            phrqc_input.append('\tSi\t1\t' + str(si['value'])+ '\n')
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
        phrqc_input.append('SiO2am\t0\t0')
        phrqc_input.append(self.c['csh']['name'] + '\t0\t0\n')
        self.phrqc_input +=  phrqc_input
        
    
    def phrqc_multilevel_voxel_CSH(self, ca, si):
        phrqc_input = [] 
        phrqc_input.append('#solution_multilevel')    
        phrqc_input.append('SOLUTION\t100004')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
            phrqc_input.append('\tSi\t' + str(si['value'])+ '\n')
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
            phrqc_input.append('\tSi\t1\t' + str(si['value'])+ '\n')
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100004')
        phrqc_input.append('SiO2am\t0\t0')
        phrqc_input.append(self.c['csh']['name'] + '\t0\t1\n')
        #phrqc_input.append(self.c['csh']['name'] + '\t0\t'+self.c['csh_mol']['value']+'\n')
        self.phrqc_input +=  phrqc_input

class CSH_Leaching(LeachingRT):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        super(CSH_Leaching, self).__init__(eqn,domain,domain_params,
             bc_params,solver_params, settings)
        self.ptype ='CSH'
    
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
        
    
    def update_init_phrqc(self):
        #self.phrqc.boundcells = self.update_boundary_cells()
        self.phrqc.init_csh = (self.solid.CSH.c>0)
        self.phrqc.nodetype = self.solid.nodetype 
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty) 
            
    def update_source(self):
        c = self.fluid.get_attr('c')
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1])
        ss = self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        if self.iters >1: pqty=self.solid.update(self.phrqc.dphases)  
        #ss = self.update_equilibrium(c,ss)
        self.set_volume()
        self.set_porosity()  
        self.set_app_tort()        
        self.update_nodetype()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)      
        self.fluid.set_attr('ss',ss) 
    
                        
    
    def update_diffusivity(self):
        Dref = self.Dref
        D_border = self.Dref
        if('D_border' in self.settings['diffusivity']):
            D_border = self.settings['diffusivity']['D_border']
        if('D_CSH' in self.settings['diffusivity']):
            D_CSH = self.settings['diffusivity']['D_CSH']
        
        is_border = self.solid.border
        is_csh = (self.solid.CSH.c >0) & (~is_border)
        is_liquid = np.logical_and(~is_csh, ~is_border)
        
        De = D_CSH*is_csh + D_border*is_border + Dref*is_liquid 
        
        self.fluid.set_attr('D0',De,component_dict=False)
        self.fluid.set_attr('Deref',np.max(De),component_dict=False)
        self.fluid.call("_set_relaxation_params")
        
    
    def update_nodetype(self):
        prev_nodetype = self.nodetype
        is_csh = (self.solid.CSH.c>0)
        is_solid = self.solid.nodetype == ct.Type.SOLID
        is_interface = np.zeros(np.shape(is_csh))
        
        is_liquid = (~is_csh)
        if(self.iters>1):
            is_liquid =  (~is_csh)&(~is_solid)
        self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
            ct.Type.MULTILEVEL * is_csh +\
            ct.Type.INTERFACE * is_interface + \
            ct.Type.SOLID * is_solid
        
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)      
        self.nodetype = self.solid.nodetype

    def update_no_flux(self, ss):        
        ss['Ca'] = ss['Ca']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1) 
        ss['Si'] = ss['Si']*(self.phrqc.boundcells==0) + 0* (self.phrqc.boundcells==1)
        return(ss)
        
    def update_border_and_interface(self, nodetype):
        is_csh = (self.solid.CSH.c>0)
        is_mineral = is_csh | (nodetype==ct.Type.SOLID)
        val = -16
        temp = deepcopy(nodetype)
        temp = temp*(~is_mineral) + val*is_mineral
        border = self.get_border(temp, val, is_csh)
        interface = self.get_interfaces(border, temp, val)
        self.solid.border = border
        self.solid.interface = interface
    
    def get_border(self, nt, val, is_port):
        rolled = np.roll(nt, -1, axis = 1)/4+np.roll(nt, 1, axis = 1)/4 +\
              np.roll(nt, -1, axis = 0)/4+np.roll(nt, 1, axis = 0)/4
        is_border = (rolled!=val)&is_port
        return(is_border)
        
    
    def update_equilibrium(self,c,ss):
        phrqc_poros = self.phrqc.selected_output()['poros']
        fr = np.zeros(phrqc_poros.shape)
        fraction = self.settings['subgrid']['fraction']
        result = {}
        by = np.where(self.solid.border)[0]
        bx = np.where(self.solid.border)[1]
        bf = np.where(self.solid.border.flatten())[0]
        for i in np.arange(0, np.sum(self.solid.border)):            
            if fraction >0:    
                fr[by[i], bx[i]] = fraction
                cell_m = bf[i]+1
                result  = self.update_phrqc_equilibrium(result, cell_m,  
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
        
        return ss
    
    
    def update_phrqc_equilibrium(self, result, n_ch, m_ch, porosity, fraction=1):
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
        
        port = output[2][9] 
        modify_str = [] 
        modify_str.append("EQUILIBRIUM_PHASES %i" %n_ch)
        modify_str.append("Portlandite 0 %.20e" %(0)) # dissolve only
        modify_str.append("END") 
        modify_str.append('USE equilibrium_phase %i' %n_ch)      
        modify_str.append('MIX %i' % ncell)      
        modify_str.append('%i 1' % ncell)   
        modify_str.append('%i %.20e' %(n_ch, 1.-fract)) #modify_str.append('%i 0' %n_int)   
        modify_str.append('SAVE solution %i' %(n_ch))  
        #modify_str.append('SAVE equilibrium_phase %i' %(n_ch))  
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
        result[str(n_ch)] = comp        
        return(result)     
    
    #%% PROPERIES
    def volume_CH(self):
        CH_vol = self.solid.portlandite.c * self.solid.portlandite.mvol
        return CH_vol