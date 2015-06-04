'''
Created on Jun 4, 2015

@author: Francois
'''

from openmdao.main.api import Assembly
from openmdao.lib.drivers.api import SLSQPdriver
from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent


# class StructuralOptimization(Assembly):
# 
#     """Unconstrained optimization of the DSM Component."""
# 
#     def __init__(self, N):
#         self.N = N
#         super(StructuralOptimization, self).__init__()
# 
#     def configure(self):
# 
#         # Create Optimizer instance
#         self.add('driver', SLSQPdriver())
# 
#         self.add('dsm_comp', DSM_component(self.N))
# 
#         # Iteration Hierarchy
#         self.driver.workflow.add('dsm_comp')
# 
#         # SLSQP Flags
#         self.driver.iprint = 3
#         self.driver.accuracy = 1e-3
# 
#         # Objective
#         self.driver.add_objective('dsm_comp.mass')
# 
#         # Design Variables
#         self.driver.add_parameter('dsm_comp.h', low=0.5, high=3.)
#         self.driver.add_parameter('dsm_comp.w', low=1., high=3.)


class AerodynamicOptimization(Assembly):

    """Constrained optimization of the LLW Component."""

    def __init__(self, N):
        self.N = N
        super(AerodynamicOptimization, self).__init__()

    def configure(self):

        # Create Optimizer instance
        self.add('driver', SLSQPdriver())
        # Create Paraboloid component instances
        self.add('llw_comp', DLLMOpenMDAOComponent(self.N))

        # Iteration Hierarchy
        self.driver.workflow.add('llw_comp')

        # SLSQP Flags
        self.driver.iprint = 3
        self.driver.accuracy = 1e-3

        # Objective
        self.driver.add_objective('llw_comp.Cd')
        self.driver.add_constraint('llw_comp.Cl=0.5')

        # Design Variables
        self.driver.add_parameter('llw_comp.AoA', low=0., high=0.3)
        self.driver.add_parameter(target='llw_comp.twist', low=-0.5, high=0.5)
        self.driver.add_parameter(target='llw_comp.thick', low=0.05, high=0.3)

if __name__ == "__main__":
    import time
    tt = time.time()
    N = 40
    opt_problem_a = AerodynamicOptimization(N)
    opt_problem_a.run()
    print "*** Elapsed time: ", time.time() - tt, "seconds ***"
    print "Minimum structure found at h=" + str(opt_problem_a.dsm_comp.h)
    print "w=" + str(opt_problem_a.dsm_comp.w)
    print "Minimum aero found at twist=" + str(opt_problem_a.llw_comp.twist)
    print "thick=" + str(opt_problem_a.llw_comp.thick)

#    print "\n"
#     print "Minimum found at (%f, %f, %f, %f)" % (opt_problem.llw_comp.thick, \
#                                          opt_problem.llw_comp.h1,\
#                                          opt_problem.llw_comp.h2,\
#                                      opt_problem.llw_comp.Cl)
