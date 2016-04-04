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
from DLLM.DLLMGeom.wing_elliptic import Wing_Elliptic
from DLLM.DLLMKernel.DLLMTargetLift import DLLMTargetLift
from MDOTools.OC.operating_condition import OperatingCondition
import numpy as np
import sys

OC=OperatingCondition('cond1')
OC.set_Mach(0.8)
OC.set_AoA(3.5)
OC.set_altitude(10000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

wing_param=Wing_Elliptic('elliptic_wing',n_sect=20)
wing_param.import_BC_from_file('input_parameters.par')
wing_param.build_linear_airfoil(OC, AoA0=0.0, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

print wing_param

x0=wing_param.get_dv_array()

DLLM = DLLMTargetLift('Simple',wing_param,OC)
DLLM.set_target_Lift(769200.)
DLLM.run_direct()
DLLM.run_post()
Ref_F = DLLM.get_F_list()

def f(x):
    wing_param.update_from_x_list(x)
    DLLM = DLLMTargetLift('Simple',wing_param,OC)
    DLLM.set_target_Lift(769200.)
    DLLM.run_direct()
    DLLM.run_post()
    func=DLLM.get_F_list()/Ref_F
    return func

def df(x):
    wing_param.update_from_x_list(x)
    DLLM = DLLMTargetLift('Simple',wing_param,OC)
    DLLM.set_target_Lift(769200.)
    DLLM.run_direct()
    DLLM.run_post()
    DLLM.run_adjoint()
    func_grad=np.array(DLLM.get_dF_list_dchi())
    N = func_grad.shape[0]
    ndv = func_grad.shape[1]
    out_grad = np.zeros((N,ndv))
    for i in xrange(N):
        out_grad[i,:] = func_grad[i,:]/Ref_F[i]
    return out_grad

val_grad=FDValidGrad(2,f,df,fd_step=1.e-8)
ok,df_fd,df=val_grad.compare(x0,treshold=1.e-5,split_out=True,return_all=True)

for j in xrange(len(df[:,0])):
    #print 'index =',j,DLLM.get_F_list_names()[j]
    fid=open('gradient_file'+str(j)+'.dat','w')
    for i in xrange(len(x0)):
        fid.write(str(i)+' '+str(df_fd[j,i])+' '+str(df[j,i])+'\n')
    fid.close()

print 'target lift  = ',DLLM.get_DLLMPost().Lift, 'target =',769200.

print '\n****************************************************'
if ok:
    print 'DLLMTargetLift gradients are valid.'
else:
    print 'DLLMTargetLift gradients are not valid!'
print '****************************************************'
