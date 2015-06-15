'''
Created on Jun 4, 2015

@author: Francois
'''

import unittest
from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent
from DLLM.DLLMKernel.DLLMPost import DLLMPost

class TestDLLM_wrapper(unittest.TestCase):
    
     
    def test_init_component(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        default_lift = 0.00
        self.assertAlmostEqual(DLLMOpenMDAO.Lift, default_lift, places=6)
  
    def test_change_OC(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        default_altitude = 10000.
        self.assertAlmostEqual(default_altitude, DLLMOpenMDAO.OC.get_altitude(), places=6)
        new_altitude = 11000.
        DLLMOpenMDAO.OC.set_altitude(11000.)
        self.assertAlmostEqual(new_altitude, DLLMOpenMDAO.OC.get_altitude(), places=6)
 
    def test_change_dv(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        new_span = 48.0
        DLLMOpenMDAO.wing_param.set_value('span', new_span)
        DLLMOpenMDAO.wing_param.update()
        self.assertAlmostEqual(new_span, DLLMOpenMDAO.wing_param.get_span())
 
    def test_exec(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        Target_Lift =  606570.049598     
        DLLMOpenMDAO.DLLM.set_target_Lift(Target_Lift)
        DLLMOpenMDAO.execute()
#         print "Final lift:",DLLMOpenMDAO.Lift
        self.assertAlmostEqual(Target_Lift, DLLMOpenMDAO.Lift,places=3)
         
    def test_list_deriv_vars(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        Target_Lift =  606570.049598     
        DLLMOpenMDAO.DLLM.set_target_Lift(Target_Lift)
        DLLMOpenMDAO.execute()
        self.assertAlmostEqual(Target_Lift, DLLMOpenMDAO.Lift,places=3)
        out_dvid_tuple, F_list_names_tuple = DLLMOpenMDAO.list_deriv_vars()
        self.assertTrue(type(out_dvid_tuple)==tuple)
        self.assertTrue(type(F_list_names_tuple)==tuple)
        F_list_names_ref = tuple(DLLMPost.DEF_F_LIST_NAMES)
        self.assertTupleEqual(F_list_names_ref, F_list_names_tuple)

    def test_jacobian(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        Target_Lift =  606570.049598     
        DLLMOpenMDAO.DLLM.set_target_Lift(Target_Lift)
        DLLMOpenMDAO.execute()
        out_dvid_tuple, F_list_names_tuple = DLLMOpenMDAO.list_deriv_vars()
        n_F = len(F_list_names_tuple)
        Jacobian = DLLMOpenMDAO.provideJ()
        n_dv = len(DLLMOpenMDAO.wing_param.get_dv_id_list())
        self.assertEqual(n_dv, Jacobian.shape[0])
        self.assertEqual(n_F, Jacobian.shape[1])
        DLLMOpenMDAO.check_gradient()
        


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLM_wrapper)
    unittest.TextTestRunner(verbosity=2).run(suite)
