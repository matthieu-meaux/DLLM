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
from DLLM.DLLMEval.DLLMMP import DLLMMP
import os
from glob import glob

config_dict={}
config_dict['Case.nb_conditions']=3
config_dict['Case.condition_name']='cond'
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

# Parameterisation configuration
config_dict['Case.param.geom_type']='Broken'
config_dict['Case.param.n_sect']=20
config_dict['Case.param.BCfilename']='input_parameters.par'
config_dict['Case.param.airfoil.type']='simple'
config_dict['Case.param.airfoil.AoA0']=-2.
config_dict['Case.param.airfoil.Cm0']=-0.1

# DLLM configuration
config_dict['Case.DLLM.type']='Solver'
config_dict['Case.DLLM.method']='inhouse'
config_dict['Case.DLLM.relax_factor']=0.99
config_dict['Case.DLLM.stop_residual']=1e-9
config_dict['Case.DLLM.max_iterations']=100
config_dict['Case.cond1.DLLM.gamma_file_name']='cond1_gamma.dat'
config_dict['Case.cond2.DLLM.gamma_file_name']='cond2_gamma.dat'
config_dict['Case.cond3.DLLM.gamma_file_name']='cond3_gamma.dat'

list_log=glob('*.log')
for log in list_log:
    os.remove(log)

MP=DLLMMP('Case')
MP.configure(config_dict)
MP.analysis()
