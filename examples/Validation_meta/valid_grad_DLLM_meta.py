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
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
from MDOTools.OC.operating_condition import OperatingCondition
import numpy

OC=OperatingCondition('cond1')
OC.set_Mach(0.6) #.7
OC.set_AoA(6.0) #3.
OC.set_altitude(10000.) #5000
OC.set_T0_deg(20.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
wing_param.build_wing()
wing_param.set_value('span',34.1)
wing_param.set_value('sweep',32.) #32.
wing_param.set_value('break_percent',33.)
wing_param.set_value('root_chord',5.4)#6.1
wing_param.set_value('break_chord',4.6)#4.6
wing_param.set_value('tip_chord',1.2)#1.5
wing_param.set_value('root_height',.98)#1.28
wing_param.set_value('break_height',0.7)#0.97
wing_param.set_value('tip_height',0.18)#0.33
wing_param.convert_to_design_variable('span',(10.,50.))
wing_param.convert_to_design_variable('sweep',(0.,40.))
wing_param.convert_to_design_variable('break_percent',(20.,40.))
wing_param.convert_to_design_variable('root_chord',(5.,7.))
wing_param.convert_to_design_variable('break_chord',(3.,5.))
wing_param.convert_to_design_variable('tip_chord',(1.,2.))
wing_param.convert_to_design_variable('root_height',(0.7,1.5))
wing_param.convert_to_design_variable('break_height',(0.45,1.2))
wing_param.convert_to_design_variable('tip_height',(0.1,0.5))
#wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
wing_param.build_meta_airfoil(OC, '../MetaModelFixed.xml', relative_thickness=.12, camber=0., Sref=1., Lref=1., sweep=.0, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

print wing_param

x0=wing_param.get_dv_array()
print 'dv array shape',x0.shape
print 'dv_array=',x0

def f(x):
    wing_param.update_from_x_list(x)
    DLLM = DLLMSolver('Meta',wing_param,OC)
    DLLM.run_direct()
    #DLLM.run_post(func_list=['Cl'])
    DLLM.run_post()
    func=DLLM.get_F_list()
    return func

def df(x):
    wing_param.update_from_x_list(x)
    DLLM = DLLMSolver('Meta',wing_param,OC)
    DLLM.run_direct()
    #DLLM.run_post(func_list=['Cl'])
    DLLM.run_post()
    DLLM.run_adjoint()
    func_grad=numpy.array(DLLM.get_dF_list_dchi())
    return func_grad

val_grad=FDValidGrad(2,f,df,fd_step=1.e-8)
ok,df_fd,df=val_grad.compare(x0,treshold=1.e-6, split_out=True,return_all=True)

for j in xrange(len(df[:,0])):
    fid=open('gradient_file'+str(j)+'.dat','w')
    for i in xrange(len(x0)):
        fid.write(str(i)+' '+str(df_fd[j,i])+' '+str(df[j,i])+'\n')
    fid.close()

print '\n****************************************************'
if ok:
    print 'Gradients are valid.'
else:
    print 'Gradients are not valid!'
print '****************************************************'
