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
#PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#=======================================================================================
from __future__ import division, print_function
from yantra import _solvers
from yantra._base import LBMeta,Fmethod,Variable
from copy import deepcopy
from AdvectionDiffusion import AdvectionDiffusion
import numpy as np
from yantra import __version__,__license__
__name__ = 'Diffusion'
__author__ = 'Ravi Patel'   
__email__ = 'ravee.a.patel@gmail.com'

class MultilevelAdvectionDiffusion2(AdvectionDiffusion):
    """
    Solves advection diffusion equation written as
    ..math::
        \partial_{t} \phi c = -\vec{nabla}  \cdot \vec{j}+ ss
        \vec{j}=  \phi c\vec{u} - \phi \zeta_{a} D_{0}\vec{\nabla}c 
    where,
    c = concentration in the aqeuous phase ..math::[N^{1}L^{-3}]
    u = velocity of the fluid (not darcy velocity) ..math::[L^{1}T^{-2}]
    j = flux ..math::[N^{1}L^{-2}T^{-1}]
    $\phi$ = porosity of the media
    $\zeta_{a}$ = apparent tortousity of the media [-]
    $D_{0}$ = diffusion coefficient in pore water..math::[L^{2}T^{-1}]
    ss = Source/sink term ..math::[N^{1}L^{-3}T^{-1}]
    The term \phi \zeta_{a} D_{0} is referred to as effective diffusivity (D_e)
    
    The bounce back condition in streaming step imposes following condition for all nodetype > 0
    \vec{\nabla}c  \cdot \hat{n} = 0 
    """
    __metaclass__= LBMeta
    _signature = 'yantra.physics.MultilevelAdvectionDiffusion.MultilevelAdvectionDiffusion'
    _vars=deepcopy(AdvectionDiffusion._vars)
    _vars['D0']= _vars.pop('D')
    _vars['D0'].type = 'parameter'
    _vars['Deref']=_vars.pop('Dref')
    _vars['collision_model'].default='TRT'
    _vars['poros']=Variable('scalar',1.,{},True,False,'domain_params','porosity')
    _vars['cphi_fact']=Variable('parameter',1./3.5,{},False,False,'solver_params','cphi fact for TRT multilevel formulation')
    _vars['cphi']=Variable('parameter',1./3.5,{},False,False,'solver_params','cphi_fact*max(poros)')
    _vars['app_tort']= Variable('scalar',1.,{},True,False,'domain_params','apparent tortousity')
    
    def _set_relaxation_params(self):
        """
        sets relaxation parameter based on effective diffusion
        """
        self.tau = self.diff2tau(D=self._De)
       
    def _get_conv_factors(self,*args):
        """
        computes conversion factors for different variables
        """
        if len(args)==3:
            domain,domain_params,solver_params = args
        klist=['D0','poros','app_tort','cphi_fact','tauref','tfactbased','tfact']
        vals=[] 
        if len(args)==3:
            for k in klist:
                if self._vars[k].belongs_to=='domain_params':
                   vals.append(domain_params.get(k,self._vars[k].default)) 
                elif self._vars[k].belongs_to=='solver_params':
                   vals.append(solver_params.get(k,self._vars[k].default))                       
                dx = getattr(domain,'dx')
        else:
            for k in klist:
                vals.append(getattr(self,k))
            dx = self.dx
        D0,poros,app_tort,cphi_fact,tauref,tfactbased,tfact=vals
        De=(D0*poros*app_tort)
        #print(De)
        #De[De==0]=np.NaN
        Deref= np.nanmax(De)
        cphi =  cphi_fact *np.nanmin(poros)
        if len(args)==3:         
            Deref= solver_params.get('Deref',Deref)
            #cphi= solver_params.get('cphi',cphi)
            print('Deref %s' %Deref)
        if len(args)==3:
            self._vars['cphi'].default = cphi
            self._vars['Deref'].default= Deref
        else:
            self.cphi = cphi
            self.Deref= Deref
        convfactors={}
        if tfactbased:
            dt= tfact * dx**2/Deref
            #print('Deref %s' %Deref)
            #print('tfact %s' %tfact)
        else:
            Dereflb=self.tau2diff(tau = tauref,cphi=cphi)
            D0= Deref/Dereflb
            #print('Dereflb %s' %Dereflb)
            #print('Deref %s' %Deref)
            #print('D0 %s' %D0)
            dt=dx**2/D0
        base = {}
        base['L']= dx
        base['N']=base['L']**3
        base['T']= dt
        convfactors = deepcopy(base)
        for k,v in self._vars.iteritems():
            if v.isphyvar:
                cfact = 1
                for p,v in v.dimension.iteritems():
                    cfact*=base.get(p,1)**v  
                convfactors[k]=cfact
        self.convfactors= convfactors
    
    def update_transport_params(self,poros,app_tort,auto_time_step=True):
        self.poros = poros
        self.app_tort = app_tort
        if auto_time_step:
            old_convfactors = deepcopy(self.convfactors)
            self._get_conv_factors()
            #reset values in LB units
            for key,val in self._vars.iteritems():
                if val.isphyvar:
                    try:
                        v=getattr(self,'_'+key)
                        v*=(old_convfactors[key]/self.convfactors[key])
                        setattr(self,'_'+key,v)          
                    except AttributeError:
                        pass
    @property
    def De(self):
        """
        effective diffusivity in physical units
        """
        #D = self.D0 * self.app_tort * self.poros
        D = self.D0 * np.ones(self.poros.shape)
        return D
    
    @property
    def _De(self):
        """
        effective diffusivity in physical units
        """
        #D = self._D0 * self.app_tort * self.poros
        D = self._D0 * np.ones(self.poros.shape)
        return D


    def tau2diff(self,**kwargs):
        """
        Gets  diffusion coefficient(in LB units) from tau
        
        Parameters
        ----------
        tau: int, float or ndarray
            relaxation parameter
        
        Returns
        -------
        D:    int, float or ndarray
            Diffusion coefficient in LB units
        """
        key_args=['tau','cphi']
        val =[]
        for k in key_args:
            try:
                val.append(kwargs.get(k,getattr(self,k)))
            except AttributeError:
                 val.append(kwargs.get(k,0))               
        return val[1]*(val[0]-0.5)
        
    def diff2tau(self,**kwargs):
        """
        gets tau from diffusion coefficient (in LB units)

        Parameters
        ----------
        De: int,float or ndarray
            effective diffusion in LB units

        Returns
        -------
        tau:    int, float or ndarray
            relaxation parameter             
        """
        key_args=['De','cphi']
        val =[]
        for k in key_args:
            try:
                if k=='De': k ='_'+k
                val.append(kwargs.get(k,getattr(self,k)))
            except AttributeError:
                 val.append(kwargs.get(k,0))
        return val[0]/val[1]+0.5

    
    @classmethod
    def _construct_fort_solver(cls,inst,solver_params):
        """
        constructs fortran solver based on given inputs
        """
        d = solver_params.get('d',inst._vars['d'].default)
        lattice = solver_params.get('lattice',inst._vars['lattice'].default)
        collision_model = solver_params.get('collision_model',inst._vars['collision_model'].default)
        lattice = lattice.upper()
        collision_model = collision_model.lower()
        if int(lattice[1]) == 2 and d == 2:
                try:
                    assert lattice == 'D2Q5'
                except AssertionError:
                    raise ValueError("Only D2Q5 lattice implementation exists for 2D")
        elif int(lattice[1]) == 3 and d == 2:
            inst._vars['d'].default=d=3
            try:
                assert lattice == 'D3Q7'
                inst._vars['q'].default=7
            except AssertionError:
                raise ValueError("Only D3Q7 lattice implementation exists for 3D")
        elif int(lattice[1]) == 3 and d == 3:
            try:
                assert lattice == 'D3Q7'
                inst._vars['q'].default=7
            except AssertionError:
                raise ValueError("Only D3Q7 lattice implementation exists for 3D")
        elif int(lattice[1]) == 2 and d == 3:
            inst._vars['lattice'].default=lattice='D3Q7'
            inst._vars['q'].default=7
        if d==2:
            ade = getattr(_solvers,'multilevel_ade2d')
        elif d==3:
            ade = getattr(_solvers,'multilevel_ade3d')   
        #compute_macro_var
        inst.compute_macro_var = Fmethod(inst,ade.compute_macro_var,
                                        ['f', '_c', '_flux', '_u','poros','nodetype','tau'])
										
        if collision_model == 'trt':
            inst.compute_edf = Fmethod(inst,ade.compute_edf,
                                       ['_c', '_u','cphi','poros','nodetype'],inplace_update=False)
            inst.collide = Fmethod(inst,ade.collide,
                                   ['f', '_c',  '_u','cphi','poros','nodetype',
                                    'tau','magic_para','_ss'])    
        if collision_model == 'diff_vel':
            inst.compute_edf = Fmethod(inst,ade.compute_edf_diff_vel,
                                       ['_c', '_u','nodetype'],inplace_update=False)
            inst.collide = Fmethod(inst,ade.collide_diff_vel,
                                   ['f', '_c', 'Dr',  '_u','u_adv','nodetype', 'poros',
                                    'tau','dc','dporos','_ss']) 
        if collision_model == 'srt':
            raise ValueError("SRT collision model is unavailable")
        '''
        #compute edf
        inst.compute_edf = Fmethod(inst,ade.compute_edf,
                                        ['_c', '_u','cphi','poros','nodetype'],inplace_update=False)

        #collison
        if collision_model == 'trt':
             inst.collide = Fmethod(inst,ade.collide,
                                   ['f', '_c',  '_u','cphi','poros','nodetype',
                                    'tau','magic_para','_ss'])      
        else:
            raise ValueError("Only TRT collision model is available")
		'''
		
        #stream
        inst.stream = Fmethod(inst,ade.stream_and_bounce,['f', 'nodetype'])      
        #apply_bc
        inst.apply_bc = Fmethod(inst,ade.boundary_conditions,
        ['f', '_u', 'nodetype', 'poros','cphi','tau', 'interp','grid_type'],['_bc'])
        
