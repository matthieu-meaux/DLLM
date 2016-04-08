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
# @author : Matthieu Meaux

import unittest
import numpy as np
from numpy import zeros, array

from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMEval.DLLMMP import DLLMMP

class TestDLLMMP(unittest.TestCase):
    
    def __get_base_config_dict(self):
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
                
        return config_dict
        
    def test_DLLM_MP_valid_grad(self):
        config_dict = self.__get_base_config_dict()
         
        MP=DLLMMP('Case')
        MP.configure(config_dict)
        MP.set_out_format('numpy')
        MP.set_grad_format('numpy')
   
        x0=MP.get_x0()
         
        Ref_list=MP.analysis()
         
        def f(x):
            func=MP.run(x)
            return func/Ref_list
           
        def df(x):
            func_grad=MP.run_grad(x)
            N = func_grad.shape[0]
            ndv = func_grad.shape[1]
            out_grad = np.zeros((N,ndv))
            for i in xrange(N):
                out_grad[i,:] = func_grad[i,:]/Ref_list[i]
            return out_grad
       
        val_grad=FDValidGrad(2,f,df,fd_step=1.e-8)
        ok,df_fd,df=val_grad.compare(x0,treshold=1.e-5,split_out=True,return_all=True)
        assert(ok)
        
    def test_DLLM_MP_AoA_valid_grad(self):
        config_dict = self.__get_base_config_dict()
        config_dict['Case.AoA_id_list']=['AoA1','AoA2','AoA3']
        config_dict['Case.param.BCfilename']='input_parameters_AoA.par'
        
        MP=DLLMMP('Case')
        MP.configure(config_dict)
        MP.set_out_format('numpy')
        MP.set_grad_format('numpy')
   
        x0=MP.get_x0()
         
        Ref_list=MP.analysis()
         
        def f(x):
            func=MP.run(x)
            return func/Ref_list
           
        def df(x):
            func_grad=MP.run_grad(x)
            N = func_grad.shape[0]
            ndv = func_grad.shape[1]
            out_grad = np.zeros((N,ndv))
            for i in xrange(N):
                out_grad[i,:] = func_grad[i,:]/Ref_list[i]
            return out_grad
       
        val_grad=FDValidGrad(2,f,df,fd_step=1.e-8)
        ok,df_fd,df=val_grad.compare(x0,treshold=1.e-5,split_out=True,return_all=True)
        assert(ok)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLMMP)
    unittest.TextTestRunner(verbosity=2).run(suite)