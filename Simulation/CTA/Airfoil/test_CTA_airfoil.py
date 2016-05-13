import numpy as np
from DLLM.polarManager.RefCTAAirfoil import RefCTAAirfoil
from MDOTools.OC.operating_condition import OperatingCondition

#Mach = 0.79
#Cl = 0.21
#Altitude 35000 ft 

OC=OperatingCondition('Cond_RefCTA',atmospheric_model='ISA')
OC.set_Mach(0.79)
OC.set_AoA(3.5)
OC.set_altitude_feet(35000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

print OC

RefCTA_AF = RefCTAAirfoil(OC)
RefCTA_AF.set_y_def_list([0.,1.520, 3.800, 6.080, 9.120, 12.160, 15.200, 18.240, 19.000])
RefCTA_AF.set_file_def_list(['../sections_data/section0.dat','../sections_data/section1.dat','../sections_data/section2.dat','../sections_data/section3.dat','../sections_data/section4.dat','../sections_data/section5.dat','../sections_data/section6.dat','../sections_data/section7.dat','../sections_data/section8.dat'])

RefCTA_AF.init_interpolators()
RefCTA_AF.set_y_pos(0.0)
AoA = -3.*np.pi/180.
RefCTA_AF.compute(AoA, 0.79)

CopyAF = RefCTA_AF.get_scaled_copy(1., 1.)