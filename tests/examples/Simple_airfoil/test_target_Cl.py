# Imports
import numpy
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMTargetCl import DLLMTargetCl
from MDOTools.OC.operating_condition import OperatingCondition

OC=OperatingCondition('cond1',atmospheric_model='simple')
OC.set_Mach(0.8)
OC.set_AoA(3.5)
OC.set_altitude(10000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

wing_param=Wing_param('test_param',geom_type='Broken',n_sect=40)
wing_param.build_wing()
wing_param.set_value('test_param.span',34.1)
wing_param.set_value('test_param.sweep',34.)
wing_param.set_value('test_param.break_percent',33.)
wing_param.set_value('test_param.root_chord',6.1)
wing_param.set_value('test_param.break_chord',4.6)
wing_param.set_value('test_param.tip_chord',1.5)
wing_param.set_value('test_param.root_height',1.28)
wing_param.set_value('test_param.break_height',0.97)
wing_param.set_value('test_param.tip_height',0.33)
wing_param.convert_to_design_variable('test_param.span',10.,50.)
wing_param.convert_to_design_variable('test_param.sweep',0.,40.)
wing_param.convert_to_design_variable('test_param.break_percent',20.,40.)
wing_param.convert_to_design_variable('test_param.root_chord',5.,7.)
wing_param.convert_to_design_variable('test_param.break_chord',3.,5.)
wing_param.convert_to_design_variable('test_param.tip_chord',1.,2.)
wing_param.convert_to_design_variable('test_param.root_height',1.,1.5)
wing_param.convert_to_design_variable('test_param.break_height',0.8,1.2)
wing_param.convert_to_design_variable('test_param.tip_height',0.2,0.5)
wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

print wing_param

DLLM = DLLMTargetCl('TCl',wing_param,OC)
#DLLM.set_method('scipy')
DLLM.set_target_Cl(0.5)
DLLM.run_direct()
DLLM.run_post()
DLLM.run_adjoint()

dF_list_dchi=DLLM.get_dF_list_dchi()

# for i,name in enumerate(DLLM.get_F_list_names()):
#     print 'gradients for '+name+' = ',dF_list_dchi[i]

