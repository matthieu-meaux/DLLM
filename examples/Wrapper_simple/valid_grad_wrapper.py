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
#  http://github.com/TBD
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
config_dict['cond1.param.airfoil.type']='simple'
config_dict['cond1.param.airfoil.AoA0']=-2.
config_dict['cond1.param.airfoil.Cm0']=-0.1
# config_dict['cond1.param.airfoil.type']='meta'
# config_dict['cond1.param.airfoil.surrogate_model']='../MetaModelCleaning.xml'
config_dict['cond1.param.desc.span.type']='DesignVariable'
config_dict['cond1.param.desc.span.value']=34.1
config_dict['cond1.param.desc.span.bounds']=(10.,50.)
config_dict['cond1.param.desc.sweep.type']='DesignVariable'
config_dict['cond1.param.desc.sweep.value']=34.
config_dict['cond1.param.desc.sweep.bounds']=(0.,40.)
config_dict['cond1.param.desc.break_percent.type']='Variable'
config_dict['cond1.param.desc.break_percent.value']=33.
#config_dict['cond1.param.desc.break_percent.bounds']=(20.,40.)
config_dict['cond1.param.desc.root_chord.type']='DesignVariable'
config_dict['cond1.param.desc.root_chord.value']=6.1
config_dict['cond1.param.desc.root_chord.bounds']=(5.,7.)
config_dict['cond1.param.desc.break_chord.type']='DesignVariable'
config_dict['cond1.param.desc.break_chord.value']=4.6
config_dict['cond1.param.desc.break_chord.bounds']=(3.,5.)
config_dict['cond1.param.desc.tip_chord.type']='DesignVariable'
config_dict['cond1.param.desc.tip_chord.value']=1.5
config_dict['cond1.param.desc.tip_chord.bounds']=(1.,2.)
config_dict['cond1.param.desc.root_height.type']='DesignVariable'
config_dict['cond1.param.desc.root_height.value']=1.28
config_dict['cond1.param.desc.root_height.bounds']=(1.,1.5)
config_dict['cond1.param.desc.break_height.type']='DesignVariable'
config_dict['cond1.param.desc.break_height.value']=0.97
config_dict['cond1.param.desc.break_height.bounds']=(0.8,1.2)
config_dict['cond1.param.desc.tip_height.type']='DesignVariable'
config_dict['cond1.param.desc.tip_height.value']=0.33
config_dict['cond1.param.desc.tip_height.bounds']=(0.2,0.5)

# DLLM configuration
config_dict['cond1.DLLM.type']='Solver'
config_dict['cond1.DLLM.method']='inhouse'
config_dict['cond1.DLLM.relax_factor']=0.99
config_dict['cond1.DLLM.stop_residual']=1e-9
config_dict['cond1.DLLM.max_iterations']=100
config_dict['cond1.DLLM.gamma_file_name']='gamma.dat'
#config_dict['cond1.DLLM.F_list_names']=['Lift','Drag','Drag_Pressure','Drag_Friction','Cl', 'Cd', 'Cdp', 'Cdf', 'LoD']
#config_dict['cond1.DLLM.F_list_names']=['Lift','Drag','Drag_Pressure','Drag_Friction','LoD']
#config_dict['cond1.DLLM.target_Cl']=0.5
#config_dict['cond1.DLLM.target_Lift']=769200.

DLLMcond1=DLLMWrapper('cond1')
DLLMcond1.configure(config_dict)
DLLMcond1.set_out_format('numpy')
DLLMcond1.set_grad_format('numpy')

x0=DLLMcond1.get_x0()
print 'dv array shape',x0.shape
print 'dv_array=',x0

val_grad=FDValidGrad(2,DLLMcond1.run,DLLMcond1.run_grad,fd_step=1.e-8)
ok,df_fd,df=val_grad.compare(x0,treshold=1.e-6,split_out=True,return_all=True)

for j in xrange(len(df[:,0])):
    fid=open('gradient_file'+str(j)+'.dat','w')
    for i in xrange(len(x0)):
        fid.write(str(i)+' '+str(df_fd[j,i])+' '+str(df[j,i])+'\n')
    fid.close()

print '\n****************************************************'
if ok:
    print 'DLLMWrapper gradients are valid.'
else:
    print 'DLLMWrapper gradients are not valid!'
print '****************************************************'
