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
from MDOTools.OC.operating_condition import OperatingCondition

OC=OperatingCondition('cond1')
OC.set_Mach(0.7)
OC.set_AoA(3.0)
OC.set_altitude(5000.)
OC.set_T0_deg(20.)
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
x_id=wing_param.get_dv_id_list()
print 'dv array shape',x0.shape
print 'dv_array=',x0

def f1(x):
    wing_param.update_from_x_list(x)
    func=wing_param.get_twist()
    return func

def df1(x):
    wing_param.update_from_x_list(x)
    func_grad=wing_param.get_twist_grad()
    return func_grad

val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-8)
ok1,df_fd1,df1=val_grad1.compare(x0,treshold=1.e-6,return_all=True)

def f2(x):
    wing_param.update_from_x_list(x)
    func=wing_param.get_chords()
    return func

def df2(x):
    wing_param.update_from_x_list(x)
    func_grad=wing_param.get_chords_grad()
    return func_grad

val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-8)
ok2,df_fd2,df2=val_grad2.compare(x0,treshold=1.e-6,return_all=True)

def f3(x):
    wing_param.update_from_x_list(x)
    func=wing_param.get_rel_thicks()
    return func

def df3(x):
    wing_param.update_from_x_list(x)
    func_grad=wing_param.get_rel_thicks_grad()
    return func_grad

val_grad3=FDValidGrad(2,f3,df3,fd_step=1.e-8)
ok3,df_fd3,df3=val_grad3.compare(x0,treshold=1.e-6,return_all=True)

def f4(x):
    wing_param.update_from_x_list(x)
    func=wing_param.get_eta()
    return func

def df4(x):
    wing_param.update_from_x_list(x)
    func_grad=wing_param.get_eta_grad()
    return func_grad

val_grad4=FDValidGrad(2,f4,df4,fd_step=1.e-8)
ok4,df_fd4,df4=val_grad4.compare(x0,treshold=1.e-6,return_all=True)

def f5(x):
    wing_param.update_from_x_list(x)
    func=wing_param.get_XYZ()
    return func

def df5(x):
    wing_param.update_from_x_list(x)
    func_grad=wing_param.get_XYZ_grad()
    return func_grad

val_grad5=FDValidGrad(2,f5,df5,fd_step=1.e-8)
ok5,df_fd5,df5=val_grad5.compare(x0,treshold=1.e-6,return_all=True)

print '\n****************************************************'
if ok1:
    print '  Twist: Gradients are valid.'
else:
    print '  Twist: Gradients are not valid !!!'
print '****************************************************'

print '\n****************************************************'
if ok2:
    print '  Chords: Gradients are valid.'
else:
    print '  Chords: Gradients are not valid !!!'
print '****************************************************'

print '\n****************************************************'
if ok3:
    print '  Rel_thicks: Gradients are valid.'
else:
    print '  Rel_thicks: Gradients are not valid !!!'
print '****************************************************'

print '\n****************************************************'
if ok4:
    print '  Eta: Gradients are valid.'
else:
    print '  Eta: Gradients are not valid !!!'
print '****************************************************'

print '\n****************************************************'
if ok5:
    print '  XYZ: Gradients are valid.'
else:
    print '  XYZ: Gradients are not valid !!!'
print '****************************************************'
