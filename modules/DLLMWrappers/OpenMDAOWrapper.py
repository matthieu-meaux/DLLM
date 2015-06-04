"""
@author Francois Gallard
"""

# pylint: disable-msg=E0611,F0401
from openmdao.main.api import Component
from openmdao.main.datatypes.api import Float, Array
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition
import numpy


class DLLMOpenMDAOComponent(Component):
    # set up interface to the framework
    # pylint: disable-msg=E1101
    Cd = Float(iotype='out', desc='Cd')
    Cl = Float(iotype='out', desc='Cl')
    Cm = Float(iotype='out', desc='Cm')
    AoA = Float(desc='AoA', default_value=0., iotype="in")
    thick = Array([], desc='Thickesses', iotype="in")
    twist = Array([], desc='Twists', iotype="in")

    def __init__(self, N):
        self.N = N
        self.thick = numpy.zeros(N)
        self.twist = numpy.zeros(N)
        super(DLLMOpenMDAOComponent, self).__init__()

    def execute(self):
        OC=OperatingCondition('cond1',atmospheric_model='ISA')
        OC.set_Mach(0.8)
        OC.set_AoA(3.5)
        OC.set_altitude(10000.)
        OC.set_T0_deg(15.)
        OC.set_P0(101325.)
        OC.set_humidity(0.)
        OC.compute_atmosphere()
        
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        wing_param.build_wing()
        wing_param.set_value('span',34.1)
        wing_param.set_value('sweep',34.)
        wing_param.set_value('break_percent',33.)
        wing_param.set_value('root_chord',6.1)
        wing_param.set_value('break_chord',4.6)
        wing_param.set_value('tip_chord',1.5)
        wing_param.set_value('root_height',1.28)
        wing_param.set_value('break_height',0.97)
        wing_param.set_value('tip_height',0.33)
        wing_param.convert_to_design_variable('span',(10.,50.))
        wing_param.convert_to_design_variable('sweep',(0.,40.))
        wing_param.convert_to_design_variable('break_percent',(20.,40.))
        wing_param.convert_to_design_variable('root_chord',(5.,7.))
        wing_param.convert_to_design_variable('break_chord',(3.,5.))
        wing_param.convert_to_design_variable('tip_chord',(1.,2.))
        wing_param.convert_to_design_variable('root_height',(1.,1.5))
        wing_param.convert_to_design_variable('break_height',(0.8,1.2))
        wing_param.convert_to_design_variable('tip_height',(0.2,0.5))
        wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        wing_param.build_airfoils_from_ref()
        wing_param.update()
        
        print wing_param
        
        DLLM = DLLMSolver('test',wing_param,OC)
        DLLM.run_direct()
        DLLM.run_post()
        DLLM.run_adjoint()

        print "DLLM flist=", DLLM.get_F_list()
        
#         
#         we = Wing_evaluator(self.N)
#         data = we.build_data()
#         we.add_dv('AoA', lbnd=0., value=self.AoA, ubnd=0.3)
#         for i, th in enumerate(self.thick):
#             we.add_dv('thick' + str(i), lbnd=0.01, value=th, ubnd=0.3)
#         for i, tw in enumerate(self.twist):
#             we.add_dv('twist' + str(i), lbnd=-0.3, value=tw, ubnd=0.3)
# 
#         # Run aero simu
#         we.run_aero(data, alpha=0., Mach=0.3)
# 
#         # Retrive functions and gradients
#         self.Cd = data['y']['Cd']
#         self.Cl = data['y']['Cl']
#         self.Cm = data['y']['Cm']
#         self.J = numpy.array(
#             [data['dy']['Cd'], data['dy']['Cl'], data['dy']['Cm']])

    def list_deriv_vars(self):
        """specified the inputs and outputs where derivatives are defined"""
        return ('AoA', 'thick', 'twist'), ('Cd', 'Cl', 'Cm',)

    def provideJ(self):
        """Calculate the Jacobian"""
        print "LLW jacobian matrix = " + str(self.J)
        return self.J
