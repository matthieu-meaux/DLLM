# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
#  DLLM (non-linear Differentiated Lifting Line Model, open source software)
# 
#  Copyright (C) 2013-2015 Airbus Group SAS
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
#  https://github.com/matthieu-meaux/DLLM.git
#
# @author : Matthieu MEAUX
#
# Imports
from MDOTools.OC.operating_condition import OperatingCondition
from DLLM.DLLMGeom.wing_broken import Wing_Broken
from DLLM.DLLMKernel.DLLMTargetCl import DLLMTargetCl

OC=OperatingCondition('cond1', atmospheric_model='ISA')
OC.set_Mach(0.8)
OC.set_AoA(3.0)
OC.set_altitude(10000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

wing_param=Wing_Broken('broken_wing',n_sect=20)
wing_param.import_BC_from_file('input_parameters.par')
wing_param.build_linear_airfoil(OC, AoA0=0.0, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()
wing_param.plot()

DLLM = DLLMTargetCl('TCl',wing_param,OC)
DLLM.set_target_Cl(0.5)
DLLM.run_direct()
DLLM.run_post()
DLLM.run_adjoint()

print 'Cl = ',DLLM.get_DLLMPost().Cl,' Target = ',0.5
print 'AoA = ',OC.get_AoA()

dF_list_dchi=DLLM.get_dF_list_dchi()