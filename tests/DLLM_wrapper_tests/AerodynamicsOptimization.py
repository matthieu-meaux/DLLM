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
        self.driver.maxiter = 400

        # Objective

        self.driver.add_objective('%s.Drag' % workflow_name)

        # Design Variables
        wing_param = self.DLLMOpenMDAO.wing_param
        print "Adding desing variables to openMDAO problem "
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
