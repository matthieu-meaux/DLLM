'''
Created on Jun 4, 2015

@author: Francois
'''
import unittest
from AerodynamicsOptimization import AerodynamicOptimization


class TestDLLM_MDAO_opt(unittest.TestCase):

    def test_optim(self):
        N = 10
        Target_Lift = 606570.049598
        opt_problem = AerodynamicOptimization(N=N, Target_Lift=Target_Lift)
        Initial_drag = 15805.247796847883
        Initial_LoD = Target_Lift / Initial_drag
        print "Running optimization"
        opt_problem.run()
        self.assertAlmostEqual(Target_Lift,
                               opt_problem.llw_comp.Lift,
                               places=3)
        self.assertGreater(Initial_drag, opt_problem.llw_comp.Drag)
        self.assertGreater(opt_problem.llw_comp.LoD, Initial_LoD)
        return

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLM_MDAO_opt)
    unittest.TextTestRunner(verbosity=2).run(suite)
