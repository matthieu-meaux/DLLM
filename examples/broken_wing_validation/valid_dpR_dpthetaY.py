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
from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMGeom.wing_broken import Wing_Broken
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition
import numpy

OC=OperatingCondition('cond1')
OC.set_Mach(0.8)
OC.set_AoA(3.5)
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

print wing_param

thetaY0=wing_param.get_thetaY()
print 'thetaY0 shape',thetaY0.shape
print 'thetaY0=',thetaY0

DLLM = DLLMSolver('Simple',wing_param,OC)
DLLM.run_direct()
iAoA=DLLM.get_iAoA()

def f(x):
    wing_param.set_thetaY(x)
    func=DLLM.comp_R(iAoA)
    return func

def df(x):
    wing_param.set_thetaY(x)
    func_grad=DLLM.comp_dpR_dpthetaY()
    return func_grad

val_grad=FDValidGrad(2,f,df,fd_step=1.e-9)
ok,df_fd,df=val_grad.compare(thetaY0,treshold=1.e-6,return_all=True)


print '\n****************************************************'
if ok:
    print 'dpR_dpthetaY is valid.'
else:
    print '!!!! dpR_dpthetaY is NOT valid !!!!'
print '****************************************************'
