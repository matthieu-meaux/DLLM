'''
Created on Jun 4, 2015

@author: Francois
'''

from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent

from openmdao.lib.drivers.api import *
from openmdao.main.api import Assembly
import numpy as np

#from pyopt_driver.pyopt_driver import pyOptDriver
class AerodynamicOptimization(Assembly):

    """Constrained optimization of the LLW Component."""

    def __init__(self, N):
        self.N = N
        super(AerodynamicOptimization, self).__init__()

    def configure(self):
        print "Configuration openmdao optimizer"
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
        rtwist_lb_array = -12.*np.ones(self.N)
        rtwist_ub_array = 12.*np.ones(self.N)
        self.driver.add_parameter('%s.rtwist' %workflow_name,
            low=rtwist_lb_array,high=rtwist_ub_array)
        self.driver.add_parameter('%s.span'%workflow_name,10.,50.,34.)
        self.driver.add_parameter('%s.sweep'%workflow_name,0.,40.,34.)
        self.driver.add_parameter('%s.break_percent'%workflow_name,20.,40.,33.)
        self.driver.add_parameter('%s.root_chord'%workflow_name,5.,7.,6.1)
        self.driver.add_parameter('%s.break_chord'%workflow_name,3.,5.,4.6)
        self.driver.add_parameter('%s.tip_chord'%workflow_name,1.,2.,1.5)
        self.driver.add_parameter('%s.root_height'%workflow_name,1.,1.5,1.28)
        self.driver.add_parameter('%s.break_height'%workflow_name,0.8,1.2,0.97)
        self.driver.add_parameter('%s.tip_height'%workflow_name,0.2,0.5,0.33)        
        
        
        #

if __name__ == "__main__":
    import time
    tt = time.time()
    N = 10
    opt_problem = AerodynamicOptimization(N)
    opt_problem.llw_comp.execute()
    print "Initial lift =", opt_problem.llw_comp.Lift
    print "Initial drag =", opt_problem.llw_comp.Drag
    print "Initial lift over drag =", opt_problem.llw_comp.LoD
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
