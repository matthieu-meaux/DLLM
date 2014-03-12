from OpenDACE.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMTargetCl import DLLMTargetCl
from MDOTools.OC.operating_condition import OperatingCondition
import numpy
import sys

OC=OperatingCondition('cond1',atmospheric_model='simple')
OC.set_Mach(0.3)
OC.set_AoA(0.)
OC.set_altitude(3000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
wing_param.build_wing()
wing_param.set_value('test_param.span',40.)
wing_param.set_value('test_param.sweep',0.)
wing_param.set_value('test_param.break_percent',33.)
wing_param.set_value('test_param.root_chord',1.0)
wing_param.set_value('test_param.break_chord',1.0)
wing_param.set_value('test_param.tip_chord',1.0)
wing_param.set_value('test_param.root_height',0.15)
wing_param.set_value('test_param.break_height',0.15)
wing_param.set_value('test_param.tip_height',0.15)
# wing_param.convert_to_design_variable('test_param.span',10.,50.)
# wing_param.convert_to_design_variable('test_param.sweep',0.,40.)
# wing_param.convert_to_design_variable('test_param.break_percent',20.,40.)
# wing_param.convert_to_design_variable('test_param.root_chord',5.,7.)
# wing_param.convert_to_design_variable('test_param.break_chord',3.,5.)
# wing_param.convert_to_design_variable('test_param.tip_chord',1.,2.)
# wing_param.convert_to_design_variable('test_param.root_height',1.,1.5)
# wing_param.convert_to_design_variable('test_param.break_height',0.8,1.2)
# wing_param.convert_to_design_variable('test_param.tip_height',0.2,0.5)
wing_param.build_linear_airfoil(OC, AoA0=0., Cm0=0.0, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

print wing_param

print 'AR=',wing_param.get_AR()

DLLM = DLLMTargetCl('Rectsweep0',wing_param,OC)
DLLM.set_target_Cl(0.5)
DLLM.run_direct()
DLLM.run_post()



