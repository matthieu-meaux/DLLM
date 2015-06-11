'''
Created on Jun 4, 2015

@author: Francois
'''

from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent

from openmdao.lib.drivers.api import *
from openmdao.main.api import Assembly
import numpy as np


class AerodynamicOptimization(Assembly):

    """Constrained optimization of the LLW Component."""

    def __init__(self, N):
        self.N = N
        super(AerodynamicOptimization, self).__init__()

    def configure(self):
        print "Configuration openmdao optimizater"
        # Create Optimizer instance
        self.add('driver', SLSQPdriver())
        # Create Paraboloid component instances
        workflow_name = 'llw_comp'
        self.DLLMOpenMDAO = DLLMOpenMDAOComponent(self.N, verbose=0)
        self.DLLMOpenMDAO.DLLM.set_target_Lift(606570.049598)
        # DLLMOpenMDAO.check_gradient()
        self.add(workflow_name, self.DLLMOpenMDAO)

        # Iteration Hierarchy
        self.driver.workflow.add(workflow_name)

        # SLSQP Flags
        self.driver.iprint = 0
        self.driver.accuracy = 1e-6
        self.driver.maxiter = 40

        # Objective

        self.driver.add_objective('%s.Drag' % workflow_name)
        # self.driver.

        # Design Variables
        wing_param = self.DLLMOpenMDAO.wing_param
        print "Adding desing variables to openMDAO problem :"
        # print 5*" "+"Variable name"+20*' '+'Lower bound'+20*" "+'Upper bound'
        rtwist_lb_array = np.zeros(self.N)
        rtwist_ub_array = np.zeros(self.N)
        twist_dv_added = False
        for i_dv_info in wing_param.get_dv_info_list():
            if i_dv_info[0].startswith('rtwist'):
                index_twist = int(i_dv_info[0].replace('rtwist', ''))
                rtwist_lb_array[index_twist] = i_dv_info[1]
                rtwist_ub_array[index_twist] = i_dv_info[3]
        for i_dv_info in wing_param.get_dv_info_list():
            if not i_dv_info[0].startswith('rtwist'):
                self.driver.add_parameter(
                    '%s.%s' %
                    (workflow_name,
                     i_dv_info[0]),
                    low=i_dv_info[1],
                    high=i_dv_info[3])
            else:
                if not twist_dv_added:
                    self.driver.add_parameter(
                        '%s.rtwist' %
                        workflow_name,
                        low=rtwist_lb_array,
                        high=rtwist_ub_array)
                    twist_dv_added = True
            # print 5*" "+i_dv_info[0]+  10*' '+   str(i_dv_info[1]) +  10*' '+   str(i_dv_info[3])
#

if __name__ == "__main__":
    import time
    tt = time.time()
    N = 10
    opt_problem = AerodynamicOptimization(N)
    print "Running optimization"
    opt_problem.run()
    print "*** Elapsed time: ", time.time() - tt, "seconds ***"

    print "Number of function evaluations =", opt_problem.llw_comp.exec_count
    print "Number of gradient evaluations =", opt_problem.llw_comp.derivative_exec_count
    print "Final lift =", opt_problem.llw_comp.Lift
    print "Final drag =", opt_problem.llw_comp.Drag
    print "Final lift over drag =", opt_problem.llw_comp.LoD
    for i_dv_info in opt_problem.DLLMOpenMDAO.wing_param.get_dv_info_list():
        print 5 * " ", i_dv_info
