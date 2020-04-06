#!/usr/bin/python
# -*- coding: utf-8 -*-
#=======================================================================================
#This File is part of Yantra: A lattice Boltzmann method based tool for multiscale/
#multiphyics simulations
#=======================================================================================
#
#Copyright (C) 2016-2017  <Author> Ravi A. Patel <Email> ravee.a.patel@gmail.com
#
#This program is free software: you can redistribute it and/or modify it under the
#terms of the GNU General Public License as published by the Free Software 
#Foundation, either version 3 of the License, or any later version.
#This program is distributed in the hope that it will be useful, but WITHOUT ANY 
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
#PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#=======================================================================================
from __future__ import division,print_function
from copy import deepcopy
import numpy as np
import yantra
from yantra._base import Multicomponent 
from yantra.physics.PhrqcReactiveTransport import Solid
from rt import CarbonationRT
import dev_phrqc as phrqc
import cell_type as ct # change the path to cell_type file

class CSHCarbonation(CarbonationRT):
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params, settings):
        if(eqn == 'MultilevelAdvectionDiffusion'): 
            eqn = 'MultilevelAdvectionDiffusion2'
        self.auto_time_step = solver_params.get('auto_time_step',True)
        self.phrqc = phrqc.CarbonationPhrqc(domain,domain_params,
                                            bc_params,solver_params)        
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
        self.nodetype = deepcopy(domain.nodetype)  
        self.apply_settings(settings)  
        self.update_volume()
        self.update_porosity()  
        self.update_app_tort() 
        self.update_init_phrqc()        
        self.update_bc(settings)
        self.update_phases()
        self.solid.transition_time = []
        self.solid.trans_port = []
        self.solid.trans_calc = []
        self.update_nodetype()
        self.update_target_SI() 
            

        
    def update_init_phrqc(self):
        self.phrqc.boundcells = self.update_boundary_cells()
        if 'portlandite' in self.solid.diffusive_phase_list:
            self.phrqc.init_port = self.solid.portlandite.c>0
        if self.ptype == 'CSH':
            self.phrqc.init_csh = np.logical_or(np.logical_or(self.solid.CSHQ_TobD.c>0, self.solid.CSHQ_TobH.c>0),
                                                np.logical_or(self.solid.CSHQ_JenD.c>0, self.solid.CSHQ_JenH.c>0))
        if 'calcite' in self.solid.diffusive_phase_list:
            self.phrqc.is_calc = self.solid.calcite.c>0
        self.phrqc._target_SI = np.zeros(self.solid.shape)  
        self.phrqc.nodetype = deepcopy(self.nodetype) 
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty) 
                
    def update_nodetype(self):
        prev_nodetype = deepcopy(self.nodetype)
        is_port = (self.solid.portlandite.c>0)
        is_calc = (self.solid.calcite.c>0)
        #is_clinker = self.solid.nodetype == ct.Type.CLINKER
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
        if self.ptype == 'CSH':
            is_csh= (self.solid.CSHQ_JenD.c>0) | (self.solid.CSHQ_JenH.c>0) | \
                    (self.solid.CSHQ_TobD.c>0) |(self.solid.CSHQ_TobH.c>0) &(~is_port) 
            is_liquid = (~is_port) & (~is_csh)
            if(self.iters>1):
                is_critical = (self.solid.pore_size <= self.solid.threshold_pore_size) & is_calc & (~is_port)
                is_liquid =  (~is_critical)&(~is_port)&(~is_solid)&(~is_calc)&(~is_csh)#&((prev_nodetype==-1)|(prev_nodetype==-2))
                is_interface = (~is_critical)&(~is_port)&(~is_solid)&(~is_liquid)&(~is_csh)
            self.solid.nodetype = ct.Type.LIQUID * is_liquid + \
                ct.Type.MULTILEVEL * is_critical + \
                ct.Type.MULTILEVEL * is_port +\
                ct.Type.MULTILEVEL * is_csh +\
                ct.Type.INTERFACE * is_interface + \
                ct.Type.SOLID * is_solid
        
        yantra._solvers.update2d.reassign_mlvl(self.solid.nodetype) 
        self.solid.nodetype[prev_nodetype == ct.Type.SOLID] = ct.Type.SOLID
        yantra._solvers.update2d.reassign_mlvl_solid(self.solid.nodetype) 
        self.solid.prev_calc_c = deepcopy(self.phrqc.solid_phase_conc['calcite'])
        self.nodetype = self.solid.nodetype
        self.update_border()
        self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)  
        if (~np.array_equal(self.nodetype, prev_nodetype)):
            self.solid.transition_time.append(self.time)
        if (~np.array_equal(self.solid.prev_calc, (self.solid.calcite.c > 1e-4))):
            self.solid.trans_calc.append(self.time)
        if (~np.array_equal(self.solid.prev_port, is_port)):
            self.solid.trans_port.append(self.time)
            
    def update_target_SI(self):
        #if (self.ptype == 'CH'):  
            self.solid.vol_ch = self.volume_CH()
            self.solid.vol_cc = self.volume_CC()          
            if (self.settings['pcs_mode']['pores']=='cylinder'):
                #TODO check this case
                if (self.iters<=1):
                    self.solid.threshold_pore_size = self.settings['pcs_mode']['crystal_size']
                    self.solid.pore_length = self.settings['pcs_mode']['crystal_size']  
                self.solid.pore_amount = self.pore_amount()  
                self.solid.free_vol_cc = (self.solid.voxel_vol-self.solid.vol_ch)* self.dx**3#self.free_volume()  
                self.solid.pore_size = self.pore_size()
                self.solid.pore_volume_cc = self.pore_volume_CC() 
            elif(self.settings['pcs_mode']['pores']=='block'):                
                if (self.iters<1):      
                    self.solid.threshold_crystal = self.settings['pcs_mode']['crystal_size'] #default crystal size
                    self.solid.threshold_pore_size = self.settings['pcs_mode']['pore_size']
                    self.solid.threshold_distance = self.solid.threshold_pore_size*2
                    self.solid.block_size = self.block_size()
                self.solid.free_vol_cc = (self.solid.voxel_vol-self.solid.vol_ch)* self.dx**3#self.solid._poros * self.dx**3#self.get_pore_volume(self)
                self.solid.ncrystals = self.block_amount()                 
                self.solid.crystal_size = self.crystal_size()
                self.solid.distance = self.distance()
                self.solid.pore_size =  self.pore_size()
            else:
                pass        
            self.solid.target_SI = self.target_SI()
            self.phrqc._target_SI = self.solid.target_SI
        #elif (self.ptype == 'CSH'): 
        #    pass     