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

class PhreeqcInputCH(PhreeqcInput):
    def make_phrqc_input(self):
        self.phrqc_boundary_voxel_Ca(self.c['ca_bc'])
        self.phrqc_liquid_voxel_Ca(self.c['ca_liq'])
        self.phrqc_multilevel_voxel_Ca(self.c['ca_mlvl'])        
        self.phrqc_solid_voxel()
        #print(self.phrqc_input)
        return(self.phrqc_input)
        
    def phrqc_boundary_voxel_Ca(self, ca):
        phrqc_input = [] 
        phrqc_input.append('#boundary_solution')    
        phrqc_input.append('SOLUTION\t100001')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']) + '\n')
        else:
            pass
        phrqc_input.append('EQUILIBRIUM_PHASES\t100001\n')
        self.phrqc_input +=  phrqc_input
        
    def phrqc_liquid_voxel_Ca(self, ca):
        phrqc_input = [] 
        phrqc_input.append('#solution_liquid')    
        phrqc_input.append('SOLUTION\t100002')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100002')
        phrqc_input.append('portlandite\t0\t0')
        self.phrqc_input +=  phrqc_input
        
    
    def phrqc_multilevel_voxel_Ca(self, ca):
        phrqc_input = [] 
        phrqc_input.append('#solution_multilevel')    
        phrqc_input.append('SOLUTION\t100003')
        phrqc_input.append('\t-units\tmol/kgw')
        phrqc_input.append('\t-water\t1')
        phrqc_input.append('\tpH\t7\tcharge')
        if(ca['type'] == 'conc'):
            phrqc_input.append('\tCa\t' + str(ca['value']))
        elif(ca['type'] == 'eq'):
            phrqc_input.append('\tCa\t1\t' + str(ca['value']))
        else:
            pass        
        phrqc_input.append('EQUILIBRIUM_PHASES\t100003')
        phrqc_input.append('portlandite\t0\t0')
        self.phrqc_input +=  phrqc_input

class CH_Leaching(LeachingRT):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        super(CH_Leaching, self).__init__(eqn,domain,domain_params,
             bc_params,solver_params, settings)
        self.ptype ='CH'
    
    def update_source(self):
        c=deepcopy( self.fluid.get_attr('c'))
        phaseqty=self.solid.phaseqty_for_phrqc()
        phase_list = deepcopy(self.solid.diffusive_phase_list)
        for num, phase in enumerate(phase_list, start=1):
            phaseqty[phase] = deepcopy(self.solid._diffusive_phaseqty[num-1])
        #if self.iters ==0:
        #    self.phrqc.modify_solid_phases(phaseqty)
        #self.modify_eq(phaseqty)  
        ss = self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)
        if self.iters >1 :pqty=self.solid.update(self.phrqc.dphases)  
        ss = self.update_portlandite_eq(c,ss)
        #if(self.settings['dissolution']=='subgrid'):
            #ss=self.update_border_solution(c,ss) 
        self.set_volume()
        self.set_porosity()  
        self.set_app_tort()        
        self.update_nodetype()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)      
        self.fluid.set_attr('ss',ss) 
        #self.fluid.call("_set_relaxation_params")
                   
    def update_diffusivity(self):
        Dref = self.Dref
        D_border = self.Dref
        if('D_border' in self.settings['diffusivity']):
            D_border = self.settings['diffusivity']['D_border']
        if('D_CH' in self.settings['diffusivity']):
            D_CH = self.settings['diffusivity']['D_CH']
        
        is_border = self.solid.border
        is_port = (self.solid.portlandite.c >0) & (~is_border)
        is_liquid = np.logical_and(~is_port, ~is_border)
        
        De = D_CH*is_port + D_border*is_border + Dref*is_liquid 
        
        self.fluid.set_attr('D0',De,component_dict=False)
        self.fluid.set_attr('Deref',np.max(De),component_dict=False)
        self.fluid.call("_set_relaxation_params")
        
    
    def update_nodetype(self):
        prev_nodetype = self.nodetype
        is_port = (self.solid.portlandite.c>0)
        is_solid = self.solid.nodetype == ct.Type.SOLID
        is_interface = np.zeros(np.shape(is_port))
        is_liquid = (~is_port)
        if(self.iters>1):
            is_liquid =  (~is_port)&(~is_solid)
        self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
            ct.Type.MULTILEVEL * is_port +\
            ct.Type.INTERFACE * is_interface + \
            ct.Type.SOLID * is_solid
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)      
        self.nodetype = self.solid.nodetype

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
    
    def volume_CH(self):
        CH_vol = self.solid.portlandite.c * self.solid.portlandite.mvol
        return CH_vol
    
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
        
    
    def update_portlandite_eq(self,c,ss):
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
    