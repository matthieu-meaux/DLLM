# Imports
from MDOTools.OC.operating_condition import OperatingCondition
from DLLM.DLLMGeom.CTA_param import CTA_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver

OC=OperatingCondition('cond1', atmospheric_model='ISA')
OC.set_Mach(0.79)
OC.set_AoA(3.5)
OC.set_altitude_feet(35000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

y_def_list=[0.,1.520, 3.800, 6.080, 9.120, 12.160, 15.200, 18.240, 19.000]
file_def_list=['../sections_data/section0.dat','../sections_data/section1.dat','../sections_data/section2.dat','../sections_data/section3.dat','../sections_data/section4.dat','../sections_data/section5.dat','../sections_data/section6.dat','../sections_data/section7.dat','../sections_data/section8.dat']

wing_param=CTA_param('broken_wing',n_sect=40, grad_active=False)
wing_param.import_BC_from_file('CTA_parameters.par')
wing_param.build_ref_CTA_airfoil(OC, y_def_list=y_def_list, file_def_list=file_def_list, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()
wing_param.plot()

wing_param.test_airfoils()