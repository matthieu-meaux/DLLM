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
from DLLM.DLLMGeom.wing_broken import Wing_Broken
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
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

x0=wing_param.get_dv_array()
print 'dv array shape',x0.shape
print 'dv_array=',x0

DLLM = DLLMSolver('Simple',wing_param,OC)
DLLM.run_direct()
iAoA=DLLM.get_iAoA()
AoA0=OC.get_AoA_rad()

# Split validation in 2 steps since the norm depends on the function value scale...
def f(x):
    OC.set_AoA_rad(x[0])
    DLLM.comp_R(iAoA)
    DLLM.set_direct_computed()
    DLLM.run_post()
    func=DLLM.get_F_list()
    return func

def df(x):
    OC.set_AoA_rad(x[0])
    DLLM.comp_R(iAoA)
    DLLM.set_direct_computed()
    DLLM.run_post()
    func_grad=DLLM.get_dpF_list_dpAoA()
    N=len(func_grad)
    np_func_grad=numpy.zeros((N,1))
    np_func_grad[:,0]=func_grad[:]
    return np_func_grad

val_grad=FDValidGrad(2,f,df,fd_step=1.e-8)
ok,df_fd,df=val_grad.compare([AoA0],treshold=1.e-6,split_out=True,return_all=True)

print '\n****************************************************'
if ok:
    print 'dpF_list_dpAoA is valid.'
else:
    print '!!!! dpF_list_dpAoA is NOT valid !!!!'
print '\n****************************************************'
