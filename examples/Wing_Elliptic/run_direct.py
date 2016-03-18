from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition
import numpy
import string

OC=OperatingCondition('cond1')
#OC.set_Mach(0.8)
#OC.set_Mach(0.6)
#OC.set_AoA(3.5)
AoA_list = [float(xx) for xx in range(0, 10)]
#AoA_list = [9.]
#Mach_list = [float(xx)/10. for xx in range(3, 9)]
Mach_list = [0.2]#,0.6,0.8]

OC.set_altitude(3000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)

wing_param=Wing_param('test_param',geom_type='Elliptic',n_sect=50)
wing_param.set_distrib_type('cos_law')
wing_param.build_wing()
wing_param.set_value('span',40.)
wing_param.set_value('root_chord',4.)
wing_param.set_value('root_height',0.1)
wing_param.set_value('tip_height',0.01)
wing_param.build_linear_airfoil(OC, AoA0=0., set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()
wing_param.plot()
print wing_param

DLLM = DLLMSolver('Elliptic',wing_param,OC)
DLLM.run_direct()
DLLM.run_post()

output = DLLM.get_F_list()
F_list_names = DLLM.get_F_list_names()   

