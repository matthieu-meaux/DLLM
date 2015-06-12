'''
Created on Jun 4, 2015

@author: Francois
'''

import unittest
from AerodynamicsOptimization import AerodynamicOptimization
from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent
   
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
#     def test_init_optim(self):
#         N = 5
#         opt_problem = AerodynamicOptimization(N)
#         assert('Cd' in dir(opt_problem.llw_comp))
        #

 #     DLLMOpenMDAO.OC.set_altitude(190000.)
#     DLLMOpenMDAO.DLLM.set_target_Lift(606570.049598)
#     DLLMOpenMDAO.run()
# 
#     print "Final point :"
#     print "    -lift :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Lift()
#     print "    -total drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag()
#     print "    -induced drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag_Induced()
#     print "    -wave drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag_Pressure()
#     print "    -wave drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag_Wave()
#     print "    -friction drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag_Friction()



if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLM_wrapper)
    unittest.TextTestRunner(verbosity=2).run(suite)
# if __name__ == "__main__":
#     import time
#     tt = time.time()
#     N = 10
#     opt_problem = AerodynamicOptimization(N)
#     opt_problem.llw_comp.execute()
#     print "Initial lift =", opt_problem.llw_comp.Lift
#     print "Initial drag =", opt_problem.llw_comp.Drag
#     print "Initial lift over drag =", opt_problem.llw_comp.LoD
#     print "Running optimization"
#     opt_problem.run()
#     print "*** Elapsed time: ", time.time() - tt, "seconds ***"
# 
#     print "Number of function evaluations =", opt_problem.llw_comp.exec_count
#     print "Number of gradient evaluations =", opt_problem.llw_comp.derivative_exec_count
#     print "Final lift =", opt_problem.llw_comp.Lift
#     print "Final drag =", opt_problem.llw_comp.Drag
#     print "Final lift over drag =", opt_problem.llw_comp.LoD
#     for i_dv_info in opt_problem.DLLMOpenMDAO.wing_param.get_dv_info_list():
#         print 5 * " ", i_dv_info
