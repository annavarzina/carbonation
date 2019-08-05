# -*- coding: utf-8 -*-

import sys,os
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)
import numpy as np
np.set_printoptions(precision=5, threshold=np.inf)
import unittest

import src.misc_func as fn

class TestSum(unittest.TestCase):
    def test_molar_volume(self):
        f_input = [1,1]
        self.assertEqual(fn.set_mvols(mvol=f_input, ptype = 'CH'), f_input)
        self.assertNotEqual(fn.set_mvols(mvol=f_input), f_input)
        f_input = [1,1,1,1,1,1]
        self.assertNotEqual(fn.set_mvols(mvol=f_input, ptype = 'CH'), f_input)
        self.assertEqual(fn.set_mvols(mvol=f_input), f_input)
        self.assertEqual(fn.set_mvols(ptype = 'CH'), [33.10e-3, 36.90e-3])
        self.assertEqual(fn.set_mvols(ptype = 'CSH'), [33.10e-3, 36.90e-3,
                         55.30e-3, 47.95e-3, 75.63e-3, 80.58e-3])
    
    def test_max_pqty(self):
        self.assertEqual(fn.get_max_pqty([2]), [0.5])
        self.assertEqual(fn.get_max_pqty([1,2]), [1,0.5])        
        
    def test_init_pqty(self):  
        self.assertEqual(fn.set_init_pqty([1]), [0.9, 0.0])
        self.assertEqual(fn.set_init_pqty([1, 2], porosCH = 0.2), [0.8, 0.0])
        
    def test_results_dict(self):
        res= fn.init_results()
        self.assertEqual(res['calcite'],[])
        self.assertTrue('sum_vol' in res.keys())
        self.assertFalse('pointparamslist' in res.keys())
        
if __name__ == '__main__':
    unittest.main()
