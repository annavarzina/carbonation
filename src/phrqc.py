# -*- coding: utf-8 -*-
import numpy as np
from copy import deepcopy
from yantra.physics._phrqc_wrapper import  Phrqc

class CarbonationPhrqc(Phrqc):
    def active_nodes(self,c,nodetype):
        active =np.ones(self.array_shape).flatten()
        prev_c = self.component_conc
        smart_inactive = np.zeros(self.array_shape)
        inactive=0
        if self.iters >1 and self.phrqc_flags['smart_run']:
            for name in self.active_components:
                diff =np.abs( c[name]*(c[name]-prev_c[name])/(c[name]+1e-30)**2)
                smart_inactive +=(1*(diff<self.phrqc_smart_run_tol))
            inactive+=1*(smart_inactive.flatten(order='c') >0) #+ 1*is_solid.flatten(order='c') #is_solid = (nodetype > 0)
        if self.phrqc_flags['only_interface']:
            inactive += (1*(nodetype!=-2)-1*(nodetype==-2)).flatten()
            tot_solid_phase_conc = self.add_dict(self.flatten_dict(self.solid_phase_conc)) 
            inactive += -1*(tot_solid_phase_conc>0)
        if self.phrqc_flags['only_fluid']: 
            inactive += (1*(nodetype>0)-1*(nodetype<=0)).flatten()
        if self.phrqc_flags['only_mineral']: 
            tot_solid_phase_conc = self.add_dict(self.flatten_dict(self.solid_phase_conc))
            inactive += -1*(tot_solid_phase_conc>0)
        active -= 1*(inactive>0)
        self.nactive = np.sum(active)
        return active

    
    def modify_eq(self, phaseqty):#, c, nodetype):
        '''
        modifies the phaseqty in phreeqc
        updates saturation index (SI) for calcite
        
        Parameters
        ----------
        phaseqty: dict
            dictionary containing an ndarray of new quantities of equilibrium phases
        '''
        #active_nodes = self.active_nodes(c,nodetype)
        phaseqty = self.flatten_dict(phaseqty)
        modifystr = []
        is_boundary = (self.boundcells ==1).flatten(order='C')
        if (self.iters==0):
            for i,cell in enumerate(range(self.startcell,self.stopcell+1,1)):
                modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)
                if(is_boundary[cell-1]):  
                    self.modify_bc(modifystr)
                else:
                    for key in phaseqty.keys():
                        modifystr.append("\t -component %s" %(key))
                        modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1])) 
        else:
            is_liquid = (self.nodetype.flatten(order='C') == -1) 
            is_mineral = (self.init_port>0).flatten(order='C')  
            si = self._target_SI.flatten(order='C')  
            for i,cell in enumerate(range(self.startcell,self.stopcell+1,1)):
                modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)                     
                if(is_boundary[cell-1]):  
                    self.modify_bc(modifystr)
                else:
                    for key in phaseqty.keys():
                        modifystr.append("\t -component %s" %(key))
                        modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1]))  
                        if (key == 'calcite' and self.pcs == True):
                            if(self.precipitation == 'interface'):                            
                                if (is_liquid[cell-1]):                                   
                                    modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
                                else:
                                    modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1])) 
                                    modifystr.append("\t\t%s\t%s" %('-dissolve_only', 0))
                            elif(self.precipitation == 'all'):
                                modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))              
                            elif(self.precipitation == 'mineral'): 
                                if (is_mineral[cell-1]): 
                                    modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))
                                    modifystr.append("\t\t%s\t%s" %('-dissolve_only', 0))                            
                                else:   
                                    modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
                            else: 
                                pass  
                            if (self.ppt == True):
                                if (phaseqty['portlandite'][cell-1]<=0):# and phaseqty['calcite'][cell-1]<=0 ):
                                    modifystr.append("\t -component\tCO2(g)") 
                                    modifystr.append("\t\tsi\t-%.20e" %self.pinput['value'])
        modifystr.append("end") 
        modifystr ='\n'.join(modifystr)
        #print(modifystr)
        '''
        import json
        with open('modify_eq.txt', 'w') as file:
            file.write(json.dumps(modifystr))
            file.close()
        '''    
        self.IPhreeqc.RunString(modifystr)
        
    def modify_bc(self, modifystr):
        modifystr.append("\t -component calcite")
        modifystr.append("\t\t%s\t%.20e" %('-moles', 0))  
        modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))
        modifystr.append("\t -component portlandite")
        modifystr.append("\t\t%s\t%.20e" %('-moles', 0))  
        modifystr.append("\t\t%s\t%s" %('-dissolve_only', 1))       
        modifystr.append("\t -component\tO2(g)") 
        modifystr.append("\t\tsi\t-3")
        '''
        modifystr.append("\t -component\tH2(g)") 
        modifystr.append("\t\tsi\t-3")
        '''
        
        if(self.pinput['type']=='pco2'):
            modifystr.append("\t -component\tCO2(g)") 
            modifystr.append("\t\tsi\t-%.20e" %self.pinput['value'])
    
         
    def modify_solid_phases(self,phaseqty):#, c, nodetype):
        """
        modifies solid phases to the qunatities given as input in phaseqty
        
        Parameters
        ----------
        phaseqty: dict
            dictonary of ndarray giving quantities of solid phases
        """
        #first sort the phases
        for name,val in phaseqty.iteritems():
            self._selected_output[name] = val
        eqphases, ssphases, kinphases = self.sort_phases(phaseqty)
        if len(eqphases)>0:
            self.modify_eq(eqphases)#, c, nodetype)
        if len(ssphases)>0:
            self.modify_ss(ssphases)            
        if len(kinphases)>0:
            self.modify_kin(kinphases)     