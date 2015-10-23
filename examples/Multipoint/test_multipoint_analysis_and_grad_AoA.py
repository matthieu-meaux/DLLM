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
#  @author : Matthieu MEAUX
#
from DLLM.DLLMEval.DLLMMP import DLLMMP
import os
from glob import glob

config_dict={}
config_dict['Case.nb_conditions']=3
config_dict['Case.condition_name']='cond'
config_dict['Case.AoA_id_list']=['AoA1','AoA2','AoA3']
# cond1 Operating condition information
config_dict['Case.cond1.OC.Mach']=0.8
config_dict['Case.cond1.OC.AoA']=3.5
config_dict['Case.cond1.OC.altitude']=10000.
# cond2 Operating condition information
config_dict['Case.cond2.OC.Mach']=0.6
config_dict['Case.cond2.OC.AoA']=4.5
config_dict['Case.cond2.OC.altitude']=5000.
# cond3 Operating condition information
config_dict['Case.cond3.OC.Mach']=0.4
config_dict['Case.cond3.OC.AoA']=6.
config_dict['Case.cond3.OC.altitude']=1000.

# Common parameterisation configuration
config_dict['Case.param.geom_type']='Broken'
config_dict['Case.param.n_sect']=20
config_dict['Case.param.airfoil.type']='simple'
config_dict['Case.param.airfoil.AoA0']=-2.
config_dict['Case.param.airfoil.Cm0']=-0.1
# config_dict['Case.param.airfoil.type']='meta'
# config_dict['Case.param.airfoil.surrogate_model']='../MetaModelCleaning.xml'
config_dict['Case.param.desc.span.type']='DesignVariable'
config_dict['Case.param.desc.span.value']=34.1
config_dict['Case.param.desc.span.bounds']=(10.,50.)
config_dict['Case.param.desc.sweep.type']='DesignVariable'
config_dict['Case.param.desc.sweep.value']=34.
config_dict['Case.param.desc.sweep.bounds']=(0.,40.)
config_dict['Case.param.desc.break_percent.type']='Variable'
config_dict['Case.param.desc.break_percent.value']=33.
#config_dict['cond1.param.desc.break_percent.bounds']=(20.,40.)
config_dict['Case.param.desc.root_chord.type']='DesignVariable'
config_dict['Case.param.desc.root_chord.value']=6.1
config_dict['Case.param.desc.root_chord.bounds']=(5.,7.)
config_dict['Case.param.desc.break_chord.type']='DesignVariable'
config_dict['Case.param.desc.break_chord.value']=4.6
config_dict['Case.param.desc.break_chord.bounds']=(3.,5.)
config_dict['Case.param.desc.tip_chord.type']='DesignVariable'
config_dict['Case.param.desc.tip_chord.value']=1.5
config_dict['Case.param.desc.tip_chord.bounds']=(1.,2.)
config_dict['Case.param.desc.root_height.type']='DesignVariable'
config_dict['Case.param.desc.root_height.value']=1.28
config_dict['Case.param.desc.root_height.bounds']=(1.,1.5)
config_dict['Case.param.desc.break_height.type']='DesignVariable'
config_dict['Case.param.desc.break_height.value']=0.97
config_dict['Case.param.desc.break_height.bounds']=(0.8,1.2)
config_dict['Case.param.desc.tip_height.type']='DesignVariable'
config_dict['Case.param.desc.tip_height.value']=0.33
config_dict['Case.param.desc.tip_height.bounds']=(0.2,0.5)

config_dict['Case.DLLM.type']='Solver'
config_dict['Case.DLLM.method']='inhouse'
config_dict['Case.DLLM.relax_factor']=0.99
config_dict['Case.DLLM.stop_residual']=1e-9
config_dict['Case.DLLM.max_iterations']=100
config_dict['Case.DLLM.F_list_names']=['Lift']
# AoA design variable must not be used with TargetCl or TargetLift
config_dict['Case.param.desc.AoA1.type']='DesignVariable'
config_dict['Case.param.desc.AoA1.value']=0.
config_dict['Case.param.desc.AoA1.bounds']=(-20.,+20.)
config_dict['Case.param.desc.AoA2.type']='DesignVariable'
config_dict['Case.param.desc.AoA2.value']=0.
config_dict['Case.param.desc.AoA2.bounds']=(-20.,+20.)
config_dict['Case.param.desc.AoA3.type']='DesignVariable'
config_dict['Case.param.desc.AoA3.value']=0.
config_dict['Case.param.desc.AoA3.bounds']=(-20.,+20.)

list_log=glob('*.log')
for log in list_log:
    os.remove(log)

MP=DLLMMP('Case')
MP.configure(config_dict)
F_list,F_list_grad=MP.analysis_and_grad()
for i,grad in enumerate(F_list_grad):
    print 'cond'+str(i+1)+' Lift grad = ',grad


# # Parameterisation configuration

# 
# # DLLM configuration
# config_dict['cond1.DLLM.type']='Solver'
# config_dict['cond1.DLLM.method']='inhouse'
# config_dict['cond1.DLLM.relax_factor']=0.99
# config_dict['cond1.DLLM.stop_residual']=1e-9
# config_dict['cond1.DLLM.max_iterations']=100
# config_dict['cond1.DLLM.gamma_file_name']='gamma.dat'
# #config_dict['cond1.DLLM.F_list_names']=['Lift','Drag','Drag_Pressure','Drag_Friction','Cl', 'Cd', 'Cdp', 'Cdf', 'LoD']
# config_dict['cond1.DLLM.F_list_names']=['Lift','Drag','Drag_Pressure','Drag_Friction','LoD']
# #config_dict['cond1.DLLM.target_Cl']=0.5
# #config_dict['cond1.DLLM.target_Lift']=769200.
# 
# DLLMcond1=DLLMWrapper('cond1')
# DLLMcond1.configure(config_dict)
# DLLMcond1.analysis_and_grad()

