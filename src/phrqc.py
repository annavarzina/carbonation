# -*- coding: utf-8 -*-
import numpy as np
from yantra.physics._phrqc_wrapper import  Phrqc

class CarbonationPhrqc(Phrqc):
    def active_nodes(self,c,nodetype):
        active =np.ones(self.array_shape).flatten()
        #is_boundary = np.zeros(self.array_shape)
        is_boundary = (self.boundcells ==1)
        prev_c = self.component_conc
        smart_inactive = np.zeros(self.array_shape)
        inactive=0
        if self.iters >1 and self.phrqc_flags['smart_run']:
            for name in self.active_components:
                diff =np.abs( c[name]*(c[name]-prev_c[name])/(c[name]+1e-30)**2)
                smart_inactive +=(1*(diff<self.phrqc_smart_run_tol))
        
        inactive+=1*(smart_inactive.flatten(order='c') >0) + 1*is_boundary.flatten(order='c')
        
        if self.phrqc_flags['only_interface']:
            inactive += (1*(nodetype!=-2)-1*(nodetype==-2)+1*is_boundary).flatten()
        #if self.phrqc_flags['only_mlvl']:
        #    inactive += (1*(nodetype==-1)-1*(nodetype!=-1)).flatten()
        if self.phrqc_flags['only_fluid']: 
            inactive += (1*(nodetype>0)-1*(nodetype<=0)+1*is_boundary).flatten()
        tot_solid_phase_conc = self.add_dict(self.flatten_dict(self.solid_phase_conc))
        inactive += -1*(tot_solid_phase_conc>0)
        active -= 1*(inactive>0)
        self.nactive = np.sum(active)
        return active

    
    def modify_eq(self, phaseqty):
        '''
        modifies the phaseqty in phreeqc
        updates saturation index (SI) for calcite
        
        Parameters
        ----------
        phaseqty: dict
            dictionary containing an ndarray of new quantities of equilibrium phases
        '''
        #print("!")
        phaseqty = self.flatten_dict(phaseqty)
        modifystr = []
        #is_boundary = (self.boundcells ==1).flatten(order='C') 
        is_liquid = np.ones(np.shape(self.poros)).flatten(order='C') 
        if self.precipitation == 'interface':    
            if self.active == 'interface':  
                si = self._target_SI.flatten(order='C')      
                for cell in range(self.startcell,self.stopcell+1,1):
                    modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)
                    for key in phaseqty.keys():
                        modifystr.append("\t -component %s" %(key))
                        modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1]))   
                        if key == 'portlandite':   
                            modifystr.append("\t\t%s\t%.20e" %('-dissolve_only', 1))
                        if (key == 'calcite'):
                            modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))
                            modifystr.append("\t\t%s\t%.20e" %('-precipitate_only', 1))
            elif self.active == 'all' or 'smart':  
                is_liquid = (self.nodetype.flatten(order='C') == -1)   
                si = self._target_SI.flatten(order='C')      
                for cell in range(self.startcell,self.stopcell+1,1):
                    modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)
                    for key in phaseqty.keys():
                        modifystr.append("\t -component %s" %(key))
                        modifystr.append("\t\t%s\t%.20e" %('-moles', phaseqty[key][cell-1]))   
                        if key == 'portlandite':   
                            modifystr.append("\t\t%s\t%.20e" %('-dissolve_only', 1))
                        if (key == 'calcite'):
                            if (is_liquid[cell-1]):                              
                                #modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))     
                                modifystr.append("\t\t%s\t%.20e" %('-dissolve_only', 1))
                            else:
                                modifystr.append("\t\t%s\t%.20e" %('-si', si[cell-1]))
                                modifystr.append("\t\t%s\t%.20e" %('-precipitate_only', 1))
            else:
                pass          
                        
        modifystr.append("end") 
        modifystr ='\n'.join(modifystr)
        self.IPhreeqc.RunString(modifystr)