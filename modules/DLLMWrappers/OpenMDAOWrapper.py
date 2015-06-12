"""
@author Francois Gallard
"""

# pylint: disable-msg=E0611,F0401
from openmdao.main.api import Component
from openmdao.main.datatypes.api import Float, Array
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMTargetLift import DLLMTargetLift
from MDOTools.OC.operating_condition import OperatingCondition
import numpy as np


class DLLMOpenMDAOComponent(Component):
    # set up interface to the framework
    # pylint: disable-msg=E1101
    # Outputs of lifting line problem
    Lift = Float(iotype='out', desc='Lift')
    Drag = Float(iotype='out', desc='Drag')
    Drag_Pressure = Float(iotype='out', desc='Drag_Pressure')
    Drag_Induced = Float(iotype='out', desc='Drag_Induced')
    Drag_Wave = Float(iotype='out', desc='Drag_Wave')
    Drag_Friction = Float(iotype='out', desc='Drag_Friction')
    Cd = Float(iotype='out', desc='Cd')
    Cdp = Float(iotype='out', desc='Cdp')
    Cdi = Float(iotype='out', desc='Cdi')
    Cdw = Float(iotype='out', desc='Cdw')
    Cdf = Float(iotype='out', desc='Cdf')
    Cl = Float(iotype='out', desc='Cl')
    LoD = Float(iotype='out', desc='LoD')
    Sref = Float(iotype='out', desc='Sref')

    # Design variables of lifting line problem
    rtwist = Array([], desc='rtwist', iotype="in")
    span = Float(desc='span', default_value=34.1, iotype="in")
    sweep = Float(desc='sweep', default_value=34., iotype="in")
    break_percent = Float(desc='break_percent', default_value=33., iotype="in")
    root_chord = Float(desc='root_chord', default_value=6.1, iotype="in")
    break_chord = Float(desc='break_chord', default_value=4.6, iotype="in")
    tip_chord = Float(desc='tip_chord', default_value=1.5, iotype="in")
    root_height = Float(desc='root_height', default_value=1.28, iotype="in")
    break_height = Float(desc='break_height', default_value=0.97, iotype="in")
    tip_height = Float(desc='tip_height', default_value=0.33, iotype="in")
    # Operating conditions variables
    Mach = Float(iotype='in', default_value=0.7, desc='Mach')
    altitude = Float(iotype='in', default_value=10000., desc='altitude')
    T0 = Float(iotype='in', default_value=15., desc='T0')
    P0 = Float(iotype='in', default_value=101325., desc='P0')
    humidity = Float(iotype='in', default_value=0.0, desc='humidity')

    def __init__(self, N, verbose=0):
        self.N = N
        self.OC = self.__set_OC(OC_name="LLW_OC")
        self.rtwist = np.zeros(N)
        super(DLLMOpenMDAOComponent, self).__init__()
        self.__display_wing_param = True
        self.__verbose = verbose
        self.__set_wing_param()

        self.DLLM = DLLMTargetLift('test', self.wing_param,
                                   self.OC, verbose=self.__verbose)
        self.DLLM.run_direct()
        self.DLLM.run_post()

    def __set_wing_param(self, wing_param_name='test_param'):

        self.__set_wing_param_values(wing_param_name=wing_param_name)

        self.__set_wing_param_bounds()
        self.wing_param.build_linear_airfoil(
            self.OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        self.wing_param.build_airfoils_from_ref()
        self.wing_param.update()

        if self.__display_wing_param:
            self.__display_wing_param = False
            print self.wing_param

    def __set_wing_param_values(self, wing_param_name='test_param'):
        self.wing_param = Wing_param(wing_param_name,
            geom_type='Broken', n_sect=self.N * 2)
        self.wing_param.build_wing()
        for param in Wing_param.DISCRETE_ATTRIBUTES_LIST:
            self.wing_param.set_value(param, eval('self.%s' % param))

    def __set_wing_param_bounds(self):
        for i in xrange(self.N):
            self.wing_param.convert_to_design_variable(
#                'rtwist%s' % i, (-float('inf'), float('inf')))
                'rtwist%s' % i, (-float('inf'), float('inf')))
        for param in Wing_param.DISCRETE_ATTRIBUTES_LIST:
            self.wing_param.convert_to_design_variable(
                param, (-float('inf'), float('inf')))
        return

    def __set_OC(self, OC_name='OC1'):
        """ Set default operating conditions"""
        OC = OperatingCondition(OC_name, atmospheric_model='ISA')
        OC.set_Mach(self.Mach)
        OC.set_altitude(self.altitude)
        OC.set_T0_deg(self.T0)
        OC.set_P0(self.P0)
        OC.set_humidity(self.humidity)
        OC.compute_atmosphere()
        return OC

    def execute(self):

        for dv_id in self.wing_param.get_dv_id_list():
            if dv_id.startswith('rtwist'):
                i_twist = int(dv_id.replace('rtwist', ''))
                self.wing_param.set_value(dv_id, self.rtwist[i_twist])
            else:
                exec("self.wing_param.set_value('%s',float(self.%s))" %
                     (dv_id, dv_id))

        self.wing_param.update()
        self.DLLM.run_direct()
        self.DLLM.run_post()
        output = self.DLLM.get_F_list()
        for f, f_name in zip(output, self.DLLM.get_F_list_names()):
            exec("self." + f_name + "=" + str(f))

    def list_deriv_vars(self):
        """specified the inputs and outputs where derivatives are defined"""

        out_dvid = []
        for dv_id in self.wing_param.get_dv_id_list():
            if dv_id.startswith('rtwist0'):
                out_dvid.append('rtwist')
            elif not dv_id.startswith('rtwist'):
                out_dvid.append(dv_id)
        return tuple(out_dvid), tuple(self.DLLM.get_F_list_names())

    def provideJ(self):
        """Calculate the Jacobian"""
        return self.DLLM.get_dpF_list_dpchi()

if __name__ == "__main__":
    import time
    tt = time.time()
    N = 20
    #options_dict = {'Mach':0.80001,'span':34.100001,'span_bounds':(10.0001,40.00001)}
    #DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=1,Mach=0.800,span=34.100001,span_bounds=(10.0001,40.00001))
    DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=1)
    DLLMOpenMDAO.OC.set_altitude(190000.)
    DLLMOpenMDAO.DLLM.set_target_Lift(606570.049598)
    DLLMOpenMDAO.run()

    print "Final point :"
    print "    -lift :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Lift()
    print "    -total drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag()
    print "    -induced drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag_Induced()
    print "    -wave drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag_Pressure()
    print "    -wave drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag_Wave()
    print "    -friction drag :", DLLMOpenMDAO.DLLM.get_DLLMPost().get_Drag_Friction()
#     print DLLMOpenMDAO.provideJ()
#     print DLLMOpenMDAO.list_deriv_vars()
#     print DLLMOpenMDAO.Cd

#    print "\n"
