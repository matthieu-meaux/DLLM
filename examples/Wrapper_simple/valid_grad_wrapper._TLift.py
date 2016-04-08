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
from DLLM.DLLMEval.DLLMWrapper import DLLMWrapper

config_dict={}
# Operating condition configuration
config_dict['cond1.OC.Mach']=0.8
config_dict['cond1.OC.AoA']=3.5
config_dict['cond1.OC.altitude']=10000.

# Parameterisation configuration
config_dict['cond1.param.geom_type']='Broken'
config_dict['cond1.param.n_sect']=20
config_dict['cond1.param.BCfilename']='input_parameters.par'
config_dict['cond1.param.airfoil.type']='simple'
config_dict['cond1.param.airfoil.AoA0']=-2.
config_dict['cond1.param.airfoil.Cm0']=-0.1

# DLLM configuration
config_dict['cond1.DLLM.type']='TargetLift'
config_dict['cond1.DLLM.method']='inhouse'
config_dict['cond1.DLLM.F_list_names']=['Drag']
config_dict['cond1.DLLM.relax_factor']=0.99
config_dict['cond1.DLLM.stop_residual']=1e-9
config_dict['cond1.DLLM.max_iterations']=100
config_dict['cond1.DLLM.gamma_file_name']='gamma.dat'
config_dict['cond1.DLLM.target_Lift']=769200.

DLLMcond1=DLLMWrapper('cond1', verbose=0)
DLLMcond1.configure(config_dict)
DLLMcond1.set_out_format('numpy')
DLLMcond1.set_grad_format('numpy')

x0=DLLMcond1.get_x0()
print 'dv array shape',x0.shape
print 'dv_array=',x0

val_grad=FDValidGrad(2,DLLMcond1.run,DLLMcond1.run_grad,fd_step=1.e-8)
ok,df_fd,df=val_grad.compare(x0,treshold=2.e-6,split_out=True,return_all=True)

for j in xrange(len(df[:,0])):
    fid=open('gradient_file'+str(j)+'.dat','w')
    for i in xrange(len(x0)):
        fid.write(str(i)+' '+str(df_fd[j,i])+' '+str(df[j,i])+'\n')
    fid.close()

# list_names= DLLMcond1.get_DLLM_solver().get_F_list_names()
# print DLLMcond1.get_DLLM_solver().get_F_list()
# 
# print 'Lift=',DLLMcond1.get_F_value('Lift'), 'target =',769200.
print '\n****************************************************'
if ok:
    print 'DLLMWrapper TargetLift gradients are valid.'
else:
    print 'DLLMWrapper TargetLift gradients are not valid!'
print '****************************************************'
