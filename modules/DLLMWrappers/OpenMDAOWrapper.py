"""
@author Francois Gallard
"""

# pylint: disable-msg=E0611,F0401
from openmdao.main.api import Component
from openmdao.main.datatypes.api import Float, Array
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition
import numpy as np


class DLLMOpenMDAOComponent(Component):
    # set up interface to the framework
    # pylint: disable-msg=E1101
    Lift            = Float(iotype='out', desc='Lift')
    Drag            = Float(iotype='out', desc='Drag')
    Drag_Pressure   = Float(iotype='out', desc='Drag_Pressure')
    Drag_Induced    = Float(iotype='out', desc='Drag_Induced')
    Drag_Wave       = Float(iotype='out', desc='Drag_Wave')
    Drag_Friction   = Float(iotype='out', desc='Drag_Friction')
    Cd              = Float(iotype='out', desc='Cd')
    Cdp             = Float(iotype='out', desc='Cdp')
    Cdi             = Float(iotype='out', desc='Cdi')
    Cdw             = Float(iotype='out', desc='Cdw')
    Cdf             = Float(iotype='out', desc='Cdf')
    Cl              = Float(iotype='out', desc='Cl')
    LoD             = Float(iotype='out', desc='LoD')
    Sref            = Float(iotype='out', desc='Sref')
    rtwist0         = Float(desc='rtwist0', default_value=0., iotype="in")
    rtwist1         = Float(desc='rtwist1', default_value=0., iotype="in")
    rtwist2         = Float(desc='rtwist2', default_value=0., iotype="in")
    rtwist3         = Float(desc='rtwist3', default_value=0., iotype="in")
    rtwist4         = Float(desc='rtwist4', default_value=0., iotype="in")
    rtwist5         = Float(desc='rtwist5', default_value=0., iotype="in")
    rtwist6         = Float(desc='rtwist6', default_value=0., iotype="in")
    rtwist7         = Float(desc='rtwist7', default_value=0., iotype="in")
    rtwist8         = Float(desc='rtwist8', default_value=0., iotype="in")
    rtwist9         = Float(desc='rtwist9', default_value=0., iotype="in")
    span            = Float(desc='span', default_value=34.1, iotype="in")
    sweep           = Float(desc='sweep', default_value=34., iotype="in")
    break_percent   = Float(desc='break_percent', default_value=33., iotype="in")
    root_chord      = Float(desc='root_chord', default_value=6.1, iotype="in")
    break_chord     = Float(desc='break_chord', default_value=4.6, iotype="in")
    tip_chord       = Float(desc='tip_chord', default_value=1.5, iotype="in")
    root_height     = Float(desc='root_height', default_value=1.28, iotype="in")
    break_height    = Float(desc='break_height', default_value=0.97, iotype="in")
    tip_height      = Float(desc='tip_height', default_value=0.33, iotype="in")

    def __init__(self, N, verbose = 0):
        self.N = N
        self.thick = np.zeros(N)
        self.twist = np.zeros(N)
        super(DLLMOpenMDAOComponent, self).__init__()
        self.__display_wing_param = True
        self.__verbose = verbose
        #self.__workflow_name = workflow_name

        OC=OperatingCondition('cond1',atmospheric_model='ISA')
        OC.set_Mach(0.8)
        OC.set_AoA(3.5)
        OC.set_altitude(10000.)
        OC.set_T0_deg(15.)
        OC.set_P0(101325.)
        OC.set_humidity(0.)
        OC.compute_atmosphere()
        self.OC = OC
        
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
        for i in xrange(self.N) :
            wing_param.set_value('rtwist%s'%i,0.0)
            wing_param.convert_to_design_variable('rtwist%s'%i,(-10.0, 10.0))
        wing_param.convert_to_design_variable('span',(10.,50.))
        wing_param.convert_to_design_variable('sweep',(0.,40.))
        wing_param.convert_to_design_variable('break_percent',(20.,40.))
        wing_param.convert_to_design_variable('root_chord',(5.,7.))
        wing_param.convert_to_design_variable('break_chord',(3.,5.))
        wing_param.convert_to_design_variable('tip_chord',(1.,2.))
        wing_param.convert_to_design_variable('root_height',(1.,1.5))
        wing_param.convert_to_design_variable('break_height',(0.8,1.2))
        wing_param.convert_to_design_variable('tip_height',(0.2,0.5))
        wing_param.build_linear_airfoil(self.OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        wing_param.build_airfoils_from_ref()
        wing_param.update()
        self.wing_param = wing_param
        
        if self.__display_wing_param :
            self.__display_wing_param = False
            print wing_param

        self.DLLM = DLLMSolver('test',self.wing_param,self.OC,verbose = self.__verbose)


    def execute(self):
        for dv_id in self.wing_param.get_dv_id_list():
            #print eval('self.rtwist%s'%i)
            self.wing_param.set_value(dv_id,eval('self.%s'%dv_id))
        self.wing_param.update()
    
        
        self.DLLM.run_direct()
        self.DLLM.run_post()
        output = self.DLLM.get_F_list()
        self.Lift=645580.617955
        for f,f_name in zip(output, self.DLLM.get_F_list_names()) :
             #eval("self.%s = %s" %(f_name,f))
             exec("self."+f_name+"="+str(f))
             
            
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
#        return ('thick','twist','sweep'),('Cl','Cd','Cm')
        return tuple(self.wing_param.get_dv_id_list()), tuple(self.DLLM.get_F_list_names())

    def provideJ(self):
        """Calculate the Jacobian"""
        #print "LLW jacobian matrix = " + str(self.J)
        self.DLLM.run_adjoint()
        return self.DLLM.get_dpF_list_dpchi()
    
if __name__ == "__main__":
    import time
    tt = time.time()
    N = 10
    DLLMOpenMDAO = DLLMOpenMDAOComponent(N,verbose = 1)
    DLLMOpenMDAO.run()
#     print DLLMOpenMDAO.provideJ()
#     print DLLMOpenMDAO.list_deriv_vars()
#     print DLLMOpenMDAO.Cd

#    print "\n"
