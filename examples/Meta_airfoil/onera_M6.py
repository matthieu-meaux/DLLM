# Imports
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition

OC=OperatingCondition('cond1')
OC.set_Mach(0.3)
OC.set_AoA(4.)
OC.set_altitude(1000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
wing_param.build_wing()
wing_param.set_value('test_param.span',2.392)  #34.1
wing_param.set_value('test_param.sweep',26.7)
wing_param.set_value('test_param.break_percent',33.) #33.
wing_param.set_value('test_param.root_chord',0.806)  #6.1
wing_param.set_value('test_param.break_chord',0.689)
wing_param.set_value('test_param.tip_chord',0.451) #1.5
wing_param.set_value('test_param.root_height',0.0782)
wing_param.set_value('test_param.break_height',0.0668)
wing_param.set_value('test_param.tip_height',0.0438)
wing_param.convert_to_design_variable('test_param.span',(10.,50.))
wing_param.convert_to_design_variable('test_param.sweep',(0.,40.))
wing_param.convert_to_design_variable('test_param.break_percent',(20.,40.))
wing_param.convert_to_design_variable('test_param.root_chord',(5.,7.))
wing_param.convert_to_design_variable('test_param.break_chord',(3.,5.))
wing_param.convert_to_design_variable('test_param.tip_chord',(1.,2.))
wing_param.convert_to_design_variable('test_param.root_height',(0.7,1.))
wing_param.convert_to_design_variable('test_param.break_height',(0.45,0.8))
wing_param.convert_to_design_variable('test_param.tip_height',(0.10,0.26))
wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
wing_param.build_meta_airfoil(OC, '../ONERA_D.xml', relative_thickness=.0, camber=0., Sref=1., Lref=1., sweep=.0, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

print wing_param

DLLM = DLLMSolver('M6wing',wing_param,OC)
DLLM.run_direct()
DLLM.run_post()
#DLLM.run_adjoint()
