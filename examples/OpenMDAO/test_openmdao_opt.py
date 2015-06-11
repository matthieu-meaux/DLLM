'''
Created on Jun 4, 2015

@author: Francois
'''

from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent

from openmdao.lib.drivers.api import *
from openmdao.main.api import Assembly




class AerodynamicOptimization(Assembly):

    """Constrained optimization of the LLW Component."""

    def __init__(self, N):
        self.N = N
        super(AerodynamicOptimization, self).__init__()

    def configure(self):

        # Create Optimizer instance
        self.add('driver', SLSQPdriver())
        # Create Paraboloid component instances
        workflow_name = 'llw_comp'
        DLLMOpenMDAO = DLLMOpenMDAOComponent(self.N,verbose = 0)
        DLLMOpenMDAO.run()
        Lift_init = DLLMOpenMDAO.Lift
        Drag_init = DLLMOpenMDAO.Drag
        print "Initial lift =",Lift_init
        print "Initial drag =",Drag_init
        print "Initial LoD=",Lift_init/Drag_init
        self.add(workflow_name, DLLMOpenMDAO)

        # Iteration Hierarchy
        self.driver.workflow.add(workflow_name)

        # SLSQP Flags
        self.driver.iprint = 0
        self.driver.accuracy = 1e-6
        self.driver.max_fun = 100

        # Objective
        
        self.driver.add_objective('%s.Drag'%workflow_name)
        #self.driver.
        self.driver.add_constraint('llw_comp.Lift='+str(Lift_init))

        # Design Variables
        wing_param = DLLMOpenMDAO.wing_param
        print "Adding desing variables to openMDAO problem :"
        print 5*" "+"Variable name"+20*' '+'Lower bound'+20*" "+'Upper bound'
        for i, i_dv_info in enumerate(wing_param.get_dv_info_list()) :
            print 5*" "+i_dv_info[0]+  10*' '+   str(i_dv_info[2]) +  10*' '+   str(i_dv_info[3])       
            self.driver.add_parameter('%s.%s' %(workflow_name,i_dv_info[0]), low=i_dv_info[2], high=i_dv_info[3])
#

if __name__ == "__main__":
    import time
    tt = time.time()
    N = 10
    opt_problem = AerodynamicOptimization(N)
    
    opt_problem.run()
    print "*** Elapsed time: ", time.time() - tt, "seconds ***"
    
    print "Final lift =",opt_problem.llw_comp.Lift
    print "Final drag =",opt_problem.llw_comp.Drag
    print "Final lift over drag =",opt_problem.llw_comp.LoD
    
