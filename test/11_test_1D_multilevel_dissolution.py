# -*- coding: utf-8 -*-
'''
Test the case with pores treated as blocks
CO2 partial pressure is fixed at the boundary
'''
#python -m unittest discover -s test run all tests
from __future__ import division  #using floating everywhere
import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import yantra
import src.cell_type as ct # change the path to cell_type file
import src.misc_func as fn
import src.rt  as rt
import unittest

class TestSum(unittest.TestCase):
    def test_run(self):
        """
        Test initialization and few iterations
        """
        #%%
        m = 'CH' #or 'CSH'
        ll = 1 #liquid lauer in front of portlandite
        l_ch = 25 #length of portlandite
        lx = (l_ch+ll)*1.0e-6
        ly = 2.0e-6
        dx = 1.0e-6
        
        domain = yantra.Domain2D(corner=(0, 0), 
                                 lengths=(lx, ly), 
                                 dx=dx, 
                                 grid_type='nodal')
        domain.nodetype[:, ll+1: ll+l_ch] = ct.Type.MULTILEVEL
        domain.nodetype[0,:] = ct.Type.SOLID
        domain.nodetype[-1,:] = ct.Type.SOLID
        domain.nodetype[:,-1] = ct.Type.SOLID
        
        #%%
        nn='10_unit_test'
        phrqc_input = {'c_bc':{'type':'conc', 'value': 0.01}, # another option 'c_bc':{'type':'pco2', 'value': 3.4}
               'c_mlvl':{'type':'conc', 'value': '0'},  # another option c_mlvl':{'type':'conc', 'value': '0'}
               'c_liq':{'type':'conc', 'value': '0'}, # another option c_liq':{'type':'conc', 'value': '0'}
               'ca_mlvl':{'type':'eq', 'value': 'portlandite'}, 
               'ca_liq':{'type':'eq', 'value': 'portlandite'}} # another option ca_liq':{'type':'conc', 'value': '0'} or ca_liq':{'type':'eq', 'value': 'portlandite'}

        phrqc = fn.set_phrqc_input(phrqc_input)            
        fn.save_phrqc_input(phrqc,root_dir, nn)   
        
        scale = 50 # scale of molar volume
        init_porosCH = 0.05 #initial porosity of portlandite nodes
        mvol_ratio = 3.69/3.31
        mvolCH = 0.0331*scale
        mvol = [mvolCH, mvolCH*mvol_ratio]
        mvol = fn.set_mvols(mvol, ptype = m) #m3/mol
        #max_pqty = fn.get_max_pqty(mvol) #mol/m3
        init_conc = fn.set_init_pqty(mvol, init_porosCH)
        pqty = fn.get_pqty(init_conc, domain)
        
        slabels = fn.set_labels(domain, m)          
        D = 1.0e-09 # default diffusion coefficient in pure liquid
        porosity = fn.get_porosity(domain, pqty, mvol, m)
        app_tort_degree = 1./3.
        app_tort = 1. * porosity ** app_tort_degree
        
        settings = {'precipitation': 'interface', # 'interface'/'all'/'mineral' nodes
                    'dissolution':'multilevel', #'multilevel'/'subgrid'
                    'active_nodes': 'all', # 'all'/'smart'/'interface'
                    'diffusivity':{'type':'archie', #'archie'/'fixed'/'mixed'
                                   }, 
                    'pcs_mode': {'pcs': True, #Pore-Size Controlled Solubility concept
                                 'pores': 'block', #'block'/'cylinder'
                                 'int_energy': 0.5, # internal energy
                                 'pore_size': 0.01*dx, # threshold radius or distance/2
                                 'crystal_size': 0.5*dx, # crystal or pore length
                                 'pore_density': 2000, #pore density per um3 - only for cylinder type
                                 }, 
                    'subgrid': {'fraction':None}, # fraction of interface cell number or None = porosity
                    'app_tort':{'degree': app_tort_degree}, 
                    'velocity': False, 
                    'bc': phrqc_input['c_bc'],
                    'dx': dx, 
                    'Dref':D
                    }
         
        domain_params = fn.set_domain_params(D, mvol, pqty, porosity, app_tort, slabels,
                                     input_file = root_dir + \
                                     '\\phreeqc_input\\' + nn + '.phrq')#'CH_CC-nat.phrq'
        bc_params = fn.set_bc_params(bc_slabels = {'left':100001})
        solver_params = fn.set_solver_params(tfact = None, smart_thres = 1e-8)# optional values, for time step (if tfact => tfactbased tau)

        #%%   
        carb_rt= rt.CarbonationRT('MultilevelAdvectionDiffusion',  domain, 
                          domain_params, bc_params, solver_params,
                          settings)
        
        itr = 0 
        nitr = 10
        while itr < nitr:     
            carb_rt.advance()    
            itr += 1
        #%%
        self.assertEqual(carb_rt.iters, nitr)

        #%%

if __name__ == '__main__':
    unittest.main()
