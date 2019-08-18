# -*- coding: utf-8 -*-
"""
Mass Balance class
@author: avarzina
"""

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy 
import cell_type as ct 

def get_cell_volumes_2D(dx, nx, ny):
    default_cell_volume = dx ** 2
    volume = np.zeros((ny+1, nx+1))
    # inner cells
    volume[1:ny, 1:nx] = 1. * default_cell_volume 
    # face cells - 4 faces
    volume[0, 1:nx] = 1./2. * default_cell_volume 
    volume[ny, 1:nx] = 1./2. * default_cell_volume 
    volume[1:ny, 0] = 1./2. * default_cell_volume 
    volume[1:ny, nx] = 1./2. * default_cell_volume 
    # corner cells - 4 corners
    volume[0, 0] = 1./4. * default_cell_volume 
    volume[0, nx] = 1./4. * default_cell_volume
    volume[ny, 0] = 1./4. * default_cell_volume 
    volume[ny, nx] = 1./4. * default_cell_volume 
    return volume      


def get_cell_volumes_3D(dx, nx, ny, nz):
    default_cell_volume = dx ** 3
    volume = np.zeros((nz+1, ny+1, nx+1))
    # inner cells
    volume[1:nz, 1:ny, 1:nx] = 1. * default_cell_volume 
    # edge cells - 6 walls
    volume[0, 1:ny, 1:nx] = 1./2. * default_cell_volume  
    volume[nz, 1:ny, 1:nx] = 1./2. * default_cell_volume  
    volume[1:nz, 0, 1:nx] = 1./2. * default_cell_volume  
    volume[1:nz, ny, 1:nx] = 1./2. * default_cell_volume  
    volume[1:nz, 1:ny, 0] = 1./2. * default_cell_volume  
    volume[1:nz, 1:ny, nx] = 1./2. * default_cell_volume  
    # face cells - 12 faces
    volume[1:nz, 0, 0] = 1./4. * default_cell_volume  
    volume[1:nz, 0, nx] = 1./4. * default_cell_volume  
    volume[1:nz, ny, 0] = 1./4. * default_cell_volume  
    volume[1:nz, ny, nx] = 1./4. * default_cell_volume 
    volume[0, 1:ny, 0] = 1./4. * default_cell_volume  
    volume[0, 1:ny, nx] = 1./4. * default_cell_volume  
    volume[nz, 1:ny, 0] = 1./4. * default_cell_volume  
    volume[nz, 1:ny, nx] = 1./4. * default_cell_volume  
    volume[0, 0, 1:nx] = 1./4. * default_cell_volume  
    volume[0, ny, 1:nx] = 1./4. * default_cell_volume  
    volume[nz, 0, 1:nx] = 1./4. * default_cell_volume  
    volume[nz, ny, 1:nx] = 1./4. * default_cell_volume   
    # corner cells - 8 corners
    volume[0, 0, 0] = 1./8. * default_cell_volume  
    volume[0, 0, nx] = 1./8. * default_cell_volume  
    volume[0, ny, 0] = 1./8. * default_cell_volume  
    volume[0, ny, nx] = 1./8. * default_cell_volume 
    volume[nz, 0, 0] = 1./8. * default_cell_volume  
    volume[nz, 0, nx] = 1./8. * default_cell_volume  
    volume[nz, ny, 0] = 1./8. * default_cell_volume  
    volume[nz, ny, nx] = 1./8. * default_cell_volume 
    return volume

def get_mass_2D(field, cellvol, poros, nx, ny):
    mass = 0
    cv = cellvol  
    f = field               
    # mass in inner cells + face cells (4) + corner cells (4)
    mass = np.sum(f[1:ny, 1:nx] * cv[1:ny, 1:nx] * poros[1:ny, 1:nx]) + \
           np.sum(f[0, 1:nx]  * cv[0, 1:nx] * poros[0, 1:nx]) + \
           np.sum(f[ny, 1:nx] * cv[ny, 1:nx] * poros[ny, 1:nx]) + \
           np.sum(f[1:ny, 0]  * cv[1:ny, 0] * poros[1:ny, 0]) + \
           np.sum(f[1:ny, nx] * cv[1:ny, nx] * poros[1:ny, nx]) + \
           f[0, 0]  * cv[0, 0] * poros[0, 0]+ \
           f[ny, 0] * cv[ny, 0] * poros[ny, 0]+ \
           f[0, nx] * cv[0, nx]  * poros[0, nx]+ \
           f[ny,nx] * cv[ny, nx]  * poros[ny, nx]
    return mass

def get_mass_3D(field, cellvol, poros, nx, ny, nz):
    mass = 0
    cv = cellvol  
    f = field           
    # mass in inner cells + edge cells (6) + face cells (12) + corner cells (8)
    mass = np.sum(f[1:nz, 1:ny, 1:nx] * cv[1:nz, 1:ny, 1:nx] * poros[1:nz, 1:ny, 1:nx]) + \
           np.sum(f[0, 1:ny, 1:nx]  * cv[0, 1:ny, 1:nx] * poros[0, 1:ny, 1:nx]) + \
           np.sum(f[nz, 1:ny, 1:nx] * cv[nz, 1:ny, 1:nx] * poros[nz, 1:ny, 1:nx]) + \
           np.sum(f[1:nz, 0, 1:nx]  * cv[1:nz, 0, 1:nx] * poros[1:nz, 0, 1:nx]) + \
           np.sum(f[1:nz, ny, 1:nx] * cv[1:nz, ny, 1:nx] * poros[1:nz, ny, 1:nx]) + \
           np.sum(f[1:nz, 1:ny, 0]  * cv[1:nz, 1:ny, 0] * poros[1:nz, 1:ny, 0]) + \
           np.sum(f[1:nz, 1:ny, nx] * cv[1:nz, 1:ny, nx] * poros[1:nz, 1:ny, nx]) + \
           np.sum(f[1:nz, 0, 0]   * cv[1:nz, 0, 0] * poros[1:nz, 0, 0]) + \
           np.sum(f[1:nz, 0, nx]  * cv[1:nz, 0, nx] * poros[1:nz, 0, nx]) + \
           np.sum(f[1:nz, ny, 0]  * cv[1:nz, ny, 0] * poros[1:nz, ny, 0]) + \
           np.sum(f[1:nz, ny, nx] * cv[1:nz, ny, nx] * poros[1:nz, ny, nx]) + \
           np.sum(f[0, 1:ny, 0]   * cv[0, 1:ny, 0] * poros[0, 1:ny, 0]) + \
           np.sum(f[0, 1:ny, nx]  * cv[0, 1:ny, nx] * poros[0, 1:ny, nx]) + \
           np.sum(f[nz, 1:ny, 0]  * cv[nz, 1:ny, 0] * poros[nz, 1:ny, 0]) + \
           np.sum(f[nz, 1:ny, nx] * cv[nz, 1:ny, nx] * poros[nz, 1:ny, nx]) + \
           np.sum(f[0, 0, 1:nx]   * cv[0, 0, 1:nx] * poros[0, 0, 1:nx]) + \
           np.sum(f[0, ny, 1:nx]  * cv[0, ny, 1:nx] * poros[0, ny, 1:nx]) + \
           np.sum(f[nz, 0, 1:nx]  * cv[nz, 0, 1:nx] * poros[nz, 0, 1:nx]) + \
           np.sum(f[nz, ny, 1:nx] * cv[nz, ny, 1:nx] * poros[nz, ny, 1:nx]) + \
           f[0, 0, 0]   * cv[0, 0, 0] * poros[0, 0, 0]+ \
           f[0, 0, nx]  * cv[0, 0, nx] * poros[0, 0, nx] + \
           f[0, ny, 0]  * cv[0, ny, 0] * poros[0, ny, 0] + \
           f[0, ny, nx] * cv[0, ny, nx] * poros[0, ny, nx] + \
           f[nz, 0, 0]  * cv[nz, 0, 0] * poros[nz, 0, 0] + \
           f[nz, 0, nx] * cv[nz, 0, nx] * poros[nz, 0, nx] + \
           f[nz, ny, 0] * cv[nz, ny, 0] * poros[nz, ny, 0] + \
           f[nz, ny,nx] * cv[nz, ny, nx] * poros[nz, ny, nx] 
    return mass

def get_boundary_mass_2D(field, cellvol, poros, nx, ny):
    mass = 0
    f = field
    cv = cellvol
    #left
    mass -= np.sum(f[0][1:ny,0] * cv[1:ny, 0] * poros[1:ny, 0]) * 2. + \
                   f[0][ny, 0]  * cv[ny, 0] * poros[ny, 0] * 2.  + \
                   f[0][0, 0]   * cv[0, 0] * poros[0, 0] * 2. 
    # right
    mass += np.sum(f[0][1:ny,nx] * cv[1:ny, nx] * poros[1:ny, nx])  * 2. +  \
                   f[0][0, nx]   * cv[0, nx] * poros[0, nx] * 2. + \
                   f[0][ny, nx]  * cv[ny, nx] * poros[ny, nx] * 2.
    # top
    mass += np.sum(f[1][0, 1:nx] * cv[0, 1:nx] * poros[0, 1:nx]) * 2. + \
                   f[1][0, 0]   * cv[0, 0] * poros[0, 0] * 2. + \
                   f[1][0, nx]  * cv[0, nx] * poros[0, nx] * 2.
    # bottom
    mass -= np.sum(f[1][ny, 1:nx]* cv[ny, 1:nx] * poros[ny, 1:nx]) * 2. + \
                   f[1][ny, nx]  * cv[ny, nx] * poros[ny, nx] * 2.+ \
                   f[1][ny, 0]   * cv[ny, 0] * poros[ny, 0] * 2. 
    #mass = mass * lbs.tfact 
    #mass *= lbs.dx
    return mass

def get_boundary_mass_3D(field, cellvol, poros, nx, ny, nz):
    mass = 0
    f = field
    cv = cellvol
    # left
    mass -= np.sum(f[0][1:nz, 1:ny, 0] * cv[1:nz, 1:ny, 0] * poros[1:nz, 1:ny, 0]) * 2. + \
         np.sum(f[0][1:nz, 0, 0] * cv[1:nz, 0, 0] * poros[1:nz, 0, 0])  * 2.  + \
         np.sum(f[0][1:nz, ny, 0]* cv[1:nz, ny, 0]* poros[1:nz, ny, 0]) * 2. +\
         np.sum(f[0][0, 1:ny, 0] * cv[0, 1:ny, 0] * poros[0, 1:ny, 0])  * 2.  + \
         np.sum(f[0][nz, 1:ny, 0]* cv[nz, 1:ny, 0]* poros[nz, 1:ny, 0]) * 2. +\
                f[0][0, 0, 0]  * cv[0, 0, 0]  * poros[0, 0, 0]  * 2.  + \
                f[0][0, ny, 0] * cv[0, ny, 0] * poros[0, ny, 0] * 2. +\
                f[0][nz, 0, 0] * cv[nz, 0, 0] * poros[nz, 0, 0] * 2.  + \
                f[0][nz, ny, 0]* cv[nz, ny, 0]* poros[nz, ny, 0]* 2. 
                 
    # right
    mass += np.sum(f[0][1:nz, 1:ny, nx] * cv[1:nz, 1:ny, nx] * poros[1:nz, 1:ny, nx]) * 2. +  \
          np.sum(f[0][1:nz, 0, nx]  * cv[1:nz, 0, nx] * poros[1:nz, 0, nx]) * 2.  + \
          np.sum(f[0][1:nz, ny, nx] * cv[1:nz, ny, nx] * poros[1:nz, ny, nx])* 2. +\
          np.sum(f[0][0, 1:ny, nx]  * cv[0, 1:ny, nx] * poros[0, 1:ny, nx]) * 2.  + \
          np.sum(f[0][nz, 1:ny, nx] * cv[nz, 1:ny, nx] * poros[nz, 1:ny, nx])* 2. +\
                 f[0][0, 0, nx]   * cv[0, 0, nx] * poros[0, 0, nx]  * 2.  + \
                 f[0][0, ny, nx]  * cv[0, ny, nx] * poros[0, ny, nx]  * 2. +\
                 f[0][nz, 0, nx]  * cv[nz, 0, nx] * poros[nz, 0, nx]  * 2.  + \
                 f[0][nz, ny, nx] * cv[nz, ny, nx]  * poros[nz, ny, nx]* 2. 
    # top
    mass += np.sum(f[1][1:nz, 0, 1:nx] * cv[1:nz, 0, 1:nx] * poros[1:nz, 0, 1:nx]) * 2. + \
          np.sum(f[1][1:nz, 0, 0]  * cv[1:nz, 0, 0] * poros[1:nz, 0, 0]) * 2.  + \
          np.sum(f[1][1:nz, 0, nx] * cv[1:nz, 0, nx] * poros[1:nz, 0, nx]) * 2. +\
          np.sum(f[1][0, 0, 1:nx]  * cv[0, 0, 1:nx] * poros[0, 0, 1:nx])  * 2.  + \
          np.sum(f[1][nz, 0, 1:nx] * cv[nz, 0, 1:nx] * poros[nz, 0, 1:nx]) * 2. +\
                 f[1][0, 0, 0]  * cv[0, 0, 0] * poros[0, 0, 0]    * 2.  + \
                 f[1][0, 0, nx]  * cv[0, 0, nx] * poros[0, 0, nx]  * 2. +\
                 f[1][nz, 0, 0]  * cv[nz, 0, 0] * poros[nz, 0, 0]  * 2.  + \
                 f[1][nz, 0, nx] * cv[nz, 0, nx] * poros[nz, 0, nx] * 2. 
    # bottom
    mass -= np.sum(f[1][1:nz, ny, 1:nx] * cv[1:nz, ny, 1:nx] * poros[1:nz, ny, 1:nx]) * 2. + \
          np.sum(f[1][1:nz, ny, 0]  * cv[1:nz, ny, 0] * poros[1:nz, ny, nx]) * 2.  + \
          np.sum(f[1][1:nz, ny, nx] * cv[1:nz, ny, nx] * poros[1:nz, ny, nx]) * 2. +\
          np.sum(f[1][0, ny, 1:nx]  * cv[0, ny, 1:nx] * poros[0, ny, 1:nx])  * 2.  + \
          np.sum(f[1][nz, ny, 1:nx] * cv[nz, ny, 1:nx] * poros[nz, ny, 1:nx]) * 2. +\
                 f[1][0, ny, 0]  * cv[0, ny, 0] * poros[0, ny, 0]    * 2.  + \
                 f[1][0, ny, nx]  * cv[0, ny, nx] * poros[0, ny, nx]  * 2. +\
                 f[1][nz, ny, 0]  * cv[nz, ny, 0] * poros[nz, ny, 0]  * 2.  + \
                 f[1][nz, ny, nx] * cv[nz, ny, nx] * poros[nz, ny, nx] * 2. 
    
    # front 
    mass += np.sum(f[2][0, 1:ny, 1:nx] * cv[0, 1:ny, 1:nx] * poros[0, 1:ny, 1:nx]) * 2. + \
          np.sum(f[2][0, 1:ny, 0]  * cv[0, 1:ny, 0] * poros[0, 1:ny, 0]) * 2.  + \
          np.sum(f[2][0, 1:ny, nx] * cv[0, 1:ny, nx] * poros[0, 1:ny, nx]) * 2. +\
          np.sum(f[2][0, 0, 1:nx]  * cv[0, 0, 1:nx] * poros[0, 0, 1:nx])  * 2.  + \
          np.sum(f[2][0, ny, 1:nx] * cv[0, ny, 1:nx] * poros[0, ny, 1:nx]) * 2. +\
                 f[2][0, 0, 0]   * cv[0, 0, 0]  * poros[0, 0, 0]   * 2.  + \
                 f[2][0, 0, nx]  * cv[0, 0, nx]  * poros[0, 0, nx] * 2. +\
                 f[2][0, ny, 0]  * cv[0, ny, 0]  * poros[0, ny, 0] * 2.  + \
                 f[2][0, ny, nx] * cv[0, ny, nx]  * poros[0, ny, nx]* 2. 
    # back
    
    mass -= np.sum(f[2][nz, 1:ny, 1:nx] * cv[nz, 1:ny, 1:nx] * poros[nz, 1:ny, 1:nx]) * 2. + \
          np.sum(f[2][nz, 1:ny, 0]  * cv[nz, 1:ny, 0] * poros[nz, 1:ny, 0]) * 2.  + \
          np.sum(f[2][nz, 1:ny, nx] * cv[nz, 1:ny, nx] * poros[nz, 1:ny, nx]) * 2. +\
          np.sum(f[2][nz, 0, 1:nx]  * cv[nz, 0, 1:nx] * poros[nz, 0, 1:nx])  * 2.  + \
          np.sum(f[2][nz, ny, 1:nx] * cv[nz, ny, 1:nx] * poros[nz, ny, 1:nx]) * 2. +\
                 f[2][nz, 0, 0]   * cv[nz, 0, 0]  * poros[nz, 0, 0]   * 2.  + \
                 f[2][nz, 0, nx]  * cv[nz, 0, nx]  * poros[nz, 0, nx] * 2. +\
                 f[2][nz, ny, 0]  * cv[nz, ny, 0]  * poros[nz, ny, 0] * 2.  + \
                 f[2][nz, ny, nx] * cv[nz, ny, nx]  * poros[nz, ny, nx]* 2. 
    #mass = mass  * lbs.tfact 
    return mass


class Mass_MCRT(object): #multicomponent
    def __init__(self, mcrt, domain):
        
        #assert (mcrt.grid_type == 'nodal'), "Mass balance for \'midway\' was not implemented yet"
        #assert (lbs.d==2 or lbs.d ==3), "Dimension should be set to 2 or 3" 
        
        phase_list = deepcopy(mcrt.solid.diffusive_phase_list)
        component_list = mcrt.fluid.components
        
        self.itr = 0
        self.iterations = []
        self.time = []
        self.set_domain_params(domain)
        self.cellvol = self.get_cell_volume() # lattice cell volume           
        for phase in phase_list:
            p = Phase(self.cellvol)
            setattr(self, phase, p)
        for component in component_list:
            c = Component()
            setattr(self, component, c)
        self.add_mass(mcrt)
        
    def add_mass(self, mcrt ):  
        """
        Adds mass values, time and iteration to class variables
        """
        self.phase_list = deepcopy(mcrt.solid.diffusive_phase_list)
        self.component_list = mcrt.fluid.components
        
        self.iterations.append(self.itr)
        self.itr += 1
        self.time.append(mcrt.time)       
        self.poros = mcrt.solid.poros
        for phase in self.phase_list:
            getattr(self, phase).get_mass(self, getattr(mcrt.solid, phase))
        for component in self.component_list:  
            getattr(self, component).get_mass(self, getattr(mcrt.fluid, component))
     
    def set_domain_params(self, domain):
        self.d = domain.d
        self.dx = domain.dx
        self.nx = domain.nx
        self.ny = domain.ny
        if domain.d == 3:
            self.nz = domain.nz
            
    def get_cell_volume(self):
        
        if self.d == 2:
            return get_cell_volumes_2D(self.dx, self.nx, self.ny)
        elif self.d ==3:
            return get_cell_volumes_3D(self.dx, self.nx, self.ny, self.nz)
        else:
            raise(ValueError)
        
    def get_cumflux(self):        
        for component in self.component_list:
            getattr(self, component).get_cumflux_mass(getattr(self, component).flux)
    
    def get_cumss(self):
        for component in self.component_list:
            getattr(self, component).get_cumss_mass(getattr(self, component).sourcesink)
    
    def get_sum_total(self):        
        for component in self.component_list:
            getattr(self, component).get_total_mass(getattr(self, component))
        
    def is_conserved(self):
        err = np.abs(self.total[self.itr-1]-self.total[0])
        is_mb = (err < 1e-8) 
        #TODO add dependency on dx
        if (is_mb):
            print("======\n Mass balance is OKAY\n======\n")
        else:
            print("======\n Mass balance is NOT OKAY. Error = %s \n=====\n" %err)
        return is_mb
        
class Phase(object):
    def __init__(self, cellvol):
        self.solid = []
        
    def get_mass(self, mass, phase):
        self.solid.append(self.get_solid_mass(mass, phase))    
    
    def get_solid_mass(self, mass, phase):
        field = getattr(phase, 'c') 
        m = 0
        if mass.d == 2:
            m = get_mass_2D(field, mass.cellvol, mass.poros, mass.nx-1, mass.ny-1)     
        elif mass.d ==3:
            m = get_mass_2D(field, mass.cellvol, mass.poros, mass.nx-1, mass.ny-1)    
        else:
           raise(ValueError) 
        return m    
        

class Component(object):
    def __init__(self):
        self.liquid = []
        self.flux = []
        self.sourcesink = []
    
    def get_mass(self, mass, comp):
        self.liquid.append(self.get_liquid_mass(mass, comp))    
        self.flux.append(self.get_flux_mass(mass, comp))    
        self.sourcesink.append(self.get_ss_mass(mass, comp))   
        
    def get_liquid_mass(self, mass, comp): 
        c = getattr(comp, 'c')
        m = 0
        if mass.d == 2:
            m = get_mass_2D(c, mass.cellvol, mass.poros, mass.nx-1, mass.ny-1)  
        #TODO 3D case
        return m
    
    def get_flux_mass(self, mass, comp):  
        flux = getattr(comp, 'flux')
        if mass.d == 2:
            m = get_boundary_mass_2D(flux, mass.cellvol, mass.poros, 
                                     mass.nx-1, mass.ny-1) 
        if mass.d == 3:
            m = get_boundary_mass_3D(flux, mass.cellvol, mass.poros, 
                                     mass.nx-1, mass.ny-1, mass.nz-1) 
            #np.mean(comp.D)#/ lbs.tfact * 10
        #TODO Find out why tfact
        return m
    
    def get_ss_mass(self, mass, comp): 
        m = 0
        if (hasattr(comp, '_ss')):
            ss = getattr(comp, '_ss')
            if mass.d == 2:
                m = get_mass_2D(ss, mass.cellvol, mass.poros, 
                                         mass.nx-1, mass.ny-1) 
            if mass.d == 3:
                m = get_mass_3D(ss, mass.cellvol, mass.poros, 
                                         mass.nx-1, mass.ny-1, mass.nz-1) 
        return m
    
    def get_cumflux_mass(self): 
        pass
        '''
        n = self.itr - 1
        mass = 0
        if(n==0): 
            mass = deepcopy(self.flux[n])
        else:
            mass = deepcopy(self.cumflux[n-1]) + deepcopy(self.flux[n])
        return mass
        '''
    
    def get_cumss_mass(self):
        pass
        '''
        n = self.itr - 1
        mass = 0
        if(n==0): 
            mass = deepcopy(self.sourcesink[n])
        else:
            mass = deepcopy(self.cumss[n-1]) + deepcopy(self.sourcesink[n])
        return mass
        '''
    
    def get_total_mass(self):  
        pass
        '''
        n = self.itr - 1        
        mass = deepcopy(self.liquid[n]) + \
               deepcopy(self.solid[n]) + \
               deepcopy(self.cumflux[n])
        return mass
        '''
    
    
#self.cumflux.append(self.cumflux_mass())  
#self.cumss.append(self.cum_ss_mass()) 
#self.total.append(self.total_mass())
#self.totalvol.append(self.total_volume(lbs))    
'''        
def get_mass(self, lbs, field):
    if(lbs.grid_type == 'nodal'): 
        mass = 0
        cv = self.cellvol 
        # 2D case
        if (lbs.d == 2):
            
            nx = lbs.nx - 1
            ny = lbs.ny - 1 
            poros = np.ones([lbs.ny,lbs.nx])
            if hasattr(lbs, 'poros'):
                poros = lbs.poros    
            mass = get_mass_2D(field, cv, poros, nx, ny)
        # 3D case
        if (lbs.d == 3):
            nx = lbs.nx - 1
            ny = lbs.ny - 1
            nz = lbs.nz - 1
            poros = np.ones([lbs.nz,lbs.ny,lbs.nx])
            if hasattr(lbs, 'poros'):
                poros = lbs.poros   
            mass = get_mass_3D(field, cv, poros, nx, ny, nz)    
        return mass
    if(lbs.grid_type == 'midway'):
        raise NotImplementedError

def get_boundary_mass(self, lbs, field):
    if(lbs.grid_type == 'nodal'): 
        mass = 0
        f = field# np array
        cv = self.cellvol
        # 2D case
        if (lbs.d == 2):
            nx = lbs.nx - 1
            ny = lbs.ny - 1
            poros = np.ones([lbs.ny,lbs.nx])
            if hasattr(lbs, 'poros'):
                poros = lbs.poros                        
            mass = get_boundary_mass_2D(field, cv, poros, nx, ny)
            
        if (lbs.d == 3):
            nx = lbs.nx - 1
            ny = lbs.ny - 1
            nz = lbs.nz - 1   
            poros = np.ones([lbs.nz,lbs.ny,lbs.nx])
            if hasattr(lbs, 'poros'):
                poros = lbs.poros  
            mass = get_boundary_mass_3D(field, cv, poros, nx, ny, nz)
        return mass
        if(lbs.grid_type == 'midway'):
            pass
        
@staticmethod
def cell_volume(lbs):
    """
    Calculates and returns volume of every grid cell on lattice.
    Called once during Mass class initialization.
    At the moment implemented only for nodal grid type. 
    Midway grid cell raises 'NotImplementedError'.
    """
    if(lbs.grid_type == 'nodal'):  
        # 2D case
        if (lbs.d == 2):
            dx = lbs.dx
            nx = lbs.nx - 1
            ny = lbs.ny - 1
            volume = get_cell_volumes_2D(dx, nx, ny)
        # 3D case
        if (lbs.d == 3):
            dx = lbs.dx
            nx = lbs.nx - 1
            ny = lbs.ny - 1
            nz = lbs.nz - 1
            default_cell_volume = dx ** 3
            volume = get_cell_volumes_3D(dx, nx, ny, nz)
        return volume
            
    if(lbs.grid_type == 'midway'):
        raise NotImplementedError

def get_error(self):
    err = np.abs(self.total[0] - self.total[self.itr - 1])
    return err
'''


   
