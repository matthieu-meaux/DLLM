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

    def __init__(self, N, verbose=0,**kwargs):
        self.N = N
        options = kwargs
        
        self.rtwist = np.zeros(N)
        super(DLLMOpenMDAOComponent, self).__init__()
        self.__display_wing_param = True
        self.__verbose = verbose
        #self.__workflow_name = workflow_name

        self.__init_OC()
        self.__init_wing_param()

        self.DLLM = DLLMTargetLift(
            'test',
            self.wing_param,
            self.OC,
            verbose=self.__verbose)
        self.DLLM.run_direct()
        self.DLLM.run_post()

    def __init_wing_param(self,wing_param_name = 'test_param', n_sect = 20):
        wing_param = Wing_param(wing_param_name, geom_type='Broken', n_sect=n_sect)
        wing_param.build_wing()
        wing_param.set_value('span', 34.1)
        wing_param.set_value('sweep', 34.)
        wing_param.set_value('break_percent', 33.)
        wing_param.set_value('root_chord', 6.1)
        wing_param.set_value('break_chord', 4.6)
        wing_param.set_value('tip_chord', 1.5)
        wing_param.set_value('root_height', 1.28)
        wing_param.set_value('break_height', 0.97)
        wing_param.set_value('tip_height', 0.33)
        for i in xrange(self.N):
            wing_param.convert_to_design_variable(
                'rtwist%s' %
                i, (-10.0, 10.0))
        wing_param.convert_to_design_variable('span', (10., 50.))
        wing_param.convert_to_design_variable('sweep', (0., 40.))
        wing_param.convert_to_design_variable('break_percent', (20., 40.))
        wing_param.convert_to_design_variable('root_chord', (5., 7.))
        wing_param.convert_to_design_variable('break_chord', (3., 5.))
        wing_param.convert_to_design_variable('tip_chord', (1., 2.))
        wing_param.convert_to_design_variable('root_height', (1., 1.5))
        wing_param.convert_to_design_variable('break_height', (0.8, 1.2))
        wing_param.convert_to_design_variable('tip_height', (0.2, 0.5))
        wing_param.build_linear_airfoil(
            self.OC,
            AoA0=-2.,
            Cm0=-0.1,
            set_as_ref=True)
        wing_param.build_airfoils_from_ref()
        wing_param.update()
        self.wing_param = wing_param

        if self.__display_wing_param:
            self.__display_wing_param = False
            print wing_param

    def __init_OC(self, OC_name = 'cond1',default = True, **kwargs):
        """Initialize operating condtions. Some default values exist
        but they can be overwritten using appropriate key word defined by 
        keywordname = keywordvalue
        example :
        __init_OC(self, OC_name = 'cond1',default = True, Mach = 0.75)
        """
        OC = self.__set_default_OC(OC_name = OC_name)
        if not default :
            for name, value in  kwargs.items():
                if name == 'Mach' :
                    OC.set_Mach(float(value))
                elif name == 'AoA' :
                    OC.set_AoA(float(value))
                elif name == 'altitude' :
                    OC.set_altitude(float(value))
                elif name == 'humidity' :
                    OC.set_humidity(float(value))
                elif name == 'P0' :
                    OC.set_P0(float(value))
                elif name == 'T0_deg' :
                    OC.set_T0_deg(float(value))
                    
        OC.compute_atmosphere()
        self.OC = OC
    
    def __set_default_OC(self, OC_name):
        """ Set default operating conditions"""
        OC = OperatingCondition(OC_name, atmospheric_model='ISA')
        OC.set_Mach(0.8)
        OC.set_AoA(3.5)
        OC.set_altitude(10000.)
        OC.set_T0_deg(15.)
        OC.set_P0(101325.)
        OC.set_humidity(0.)
        OC.compute_atmosphere()
        return OC

    def execute(self):

        for dv_id in self.wing_param.get_dv_id_list():
            if dv_id.startswith('rtwist'):
                i_twist = int(dv_id.replace('rtwist', ''))
                self.wing_param.set_value(dv_id, self.rtwist[i_twist])
                # print dv_id,self.rtwist[i_twist]
            else:
                #                    self.wing_param.set_value(dv_id,eval('self.%s'%dv_id))
                exec(
                    "self.wing_param.set_value('%s',self.%s)" %
                    (dv_id, dv_id))
                # print dv_id,eval('self.%s'%dv_id)

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
        # print "LLW jacobian matrix = " + str(self.J)
        return self.DLLM.get_dpF_list_dpchi()

if __name__ == "__main__":
    import time
    tt = time.time()
    N = 10
    DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=1)
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
