'''
Created on Jun 4, 2015

@author: Francois
'''
import unittest

from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent

from openmdao.lib.drivers.api import *
from openmdao.main.api import Assembly
import numpy as np
from AerodynamicsOptimization import AerodynamicOptimization

#from pyopt_driver.pyopt_driver import pyOptDriver
class TestDLLM_MDAO_opt(unittest.TestCase):
    
    def test_optim(self):
        N = 10
        opt_problem = AerodynamicOptimization(N)
        parameters_dict = opt_problem.driver.get_parameters()
        out_dvid_tuple, F_list_names_tuple = opt_problem.DLLMOpenMDAO.list_deriv_vars()
        self.assertEqual(len(parameters_dict.keys()), len(out_dvid_tuple))
        
        Initial_Lift =  606570.049598
        Initial_drag = 15751.3017003
        Initial_LoD = Initial_Lift/Initial_drag
        opt_problem.llw_comp.execute()
        self.assertAlmostEqual(Initial_Lift, opt_problem.llw_comp.Lift, places =3)
        self.assertAlmostEqual(Initial_drag, opt_problem.llw_comp.Drag, places =3)
        self.assertAlmostEqual(Initial_LoD, opt_problem.llw_comp.LoD, places =3)        
                
        print "Running optimization"
        opt_problem.run()
        
        self.assertAlmostEqual(Initial_Lift, opt_problem.llw_comp.Lift, places =3)
        self.assertGreater(Initial_drag, opt_problem.llw_comp.Drag)
        self.assertGreater( opt_problem.llw_comp.LoD,Initial_LoD)
        return

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLM_MDAO_opt)
    unittest.TextTestRunner(verbosity=2).run(suite)

