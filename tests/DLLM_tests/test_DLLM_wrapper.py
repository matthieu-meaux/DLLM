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
from DLLM.DLLMEval.DLLMWrapper import DLLMWrapper

class TestDLLMWrapper(unittest.TestCase):
    
    def __get_base_config_dict(self):
        config_dict={}
        # Operating condition configuration
        config_dict['test.OC.Mach']=0.8
        config_dict['test.OC.AoA']=3.5
        config_dict['cond1.OC.altitude']=10000.
        
        # Parameterisation configuration
        config_dict['test.param.geom_type']='Broken'
        config_dict['test.param.n_sect']=20
        config_dict['test.param.BCfilename']='input_parameters.par'
        config_dict['test.param.airfoil.type']='simple'
        config_dict['test.param.airfoil.AoA0']=-2.
        
        # DLLM configuration
        config_dict['test.DLLM.type']='Solver'
        config_dict['test.DLLM.method']='inhouse'
        config_dict['test.DLLM.relax_factor']=0.99
        config_dict['test.DLLM.stop_residual']=1e-9
        config_dict['test.DLLM.max_iterations']=100
        config_dict['test.DLLM.gamma_file_name']='gamma.dat'
        
        return config_dict
        
    def test_DLLM_wrapper_valid_grad(self):
        config_dict = self.__get_base_config_dict()
         
        DLLMWrap = DLLMWrapper('test', verbose=0)
        DLLMWrap.configure(config_dict)
        DLLMWrap.set_out_format('numpy')
        DLLMWrap.set_grad_format('numpy')
   
        x0=DLLMWrap.get_x0()
         
        DLLMWrap.analysis()
         
        Ref_list=DLLMWrap.get_F_list()
         
        def f(x):
            DLLMWrap.run(x)
            func=DLLMWrap.get_F_list()
            return func/Ref_list
           
        def df(x):
            DLLMWrap.run_grad(x)
            func_grad=np.array(DLLMWrap.get_F_list_grad())
            N = func_grad.shape[0]
            ndv = func_grad.shape[1]
            out_grad = np.zeros((N,ndv))
            for i in xrange(N):
                out_grad[i,:] = func_grad[i,:]/Ref_list[i]
            return out_grad
       
        val_grad=FDValidGrad(2,f,df,fd_step=1.e-8)
        ok,df_fd,df=val_grad.compare(x0,treshold=1.e-5,split_out=True,return_all=True)
        assert(ok)
        
    def test_DLLM_wrapper_json_valid_grad(self):
        DLLMWrap = DLLMWrapper('test', verbose=0)
        DLLMWrap.config_from_file('conf_file.json')
        DLLMWrap.set_out_format('numpy')
        DLLMWrap.set_grad_format('numpy')
   
        x0=DLLMWrap.get_x0()
         
        DLLMWrap.analysis()
         
        Ref_list=DLLMWrap.get_F_list()
         
        def f(x):
            DLLMWrap.run(x)
            func=DLLMWrap.get_F_list()
            return func/Ref_list
           
        def df(x):
            DLLMWrap.run_grad(x)
            func_grad=np.array(DLLMWrap.get_F_list_grad())
            N = func_grad.shape[0]
            ndv = func_grad.shape[1]
            out_grad = np.zeros((N,ndv))
            for i in xrange(N):
                out_grad[i,:] = func_grad[i,:]/Ref_list[i]
            return out_grad
       
        val_grad=FDValidGrad(2,f,df,fd_step=1.e-8)
        ok,df_fd,df=val_grad.compare(x0,treshold=1.e-5,split_out=True,return_all=True)
        assert(ok)
        
    def test_DLLM_wrapper_TCl(self):
        config_dict = self.__get_base_config_dict()
        config_dict['test.DLLM.type']='TargetCl'
        config_dict['test.DLLM.target_Cl']=0.5
        
        DLLMWrap = DLLMWrapper('test', verbose=0)
        DLLMWrap.configure(config_dict)
        DLLMWrap.set_out_format('numpy')
        DLLMWrap.set_grad_format('numpy')
        
        DLLMWrap.analysis()
        Cl=DLLMWrap.get_F_value('Cl')
        assert((Cl-0.5)<1.e-8)
        
    def test_DLLM_wrapper_TCl_valid_grad(self):
        config_dict = self.__get_base_config_dict()
        config_dict['test.DLLM.type']='TargetCl'
        config_dict['test.DLLM.target_Cl']=0.5
        
        DLLMWrap = DLLMWrapper('test', verbose=0)
        DLLMWrap.configure(config_dict)
        DLLMWrap.set_out_format('numpy')
        DLLMWrap.set_grad_format('numpy')
        
        x0=DLLMWrap.get_x0()
         
        DLLMWrap.analysis()
         
        Ref_list=DLLMWrap.get_F_list()
         
        def f(x):
            DLLMWrap.run(x)
            func=DLLMWrap.get_F_list()
            return func/Ref_list
           
        def df(x):
            DLLMWrap.run_grad(x)
            func_grad=np.array(DLLMWrap.get_F_list_grad())
            N = func_grad.shape[0]
            ndv = func_grad.shape[1]
            out_grad = np.zeros((N,ndv))
            for i in xrange(N):
                out_grad[i,:] = func_grad[i,:]/Ref_list[i]
            return out_grad
       
        val_grad=FDValidGrad(2,f,df,fd_step=1.e-8)
        ok,df_fd,df=val_grad.compare(x0,treshold=1.e-5,split_out=True,return_all=True)
        assert(ok)
        
    def test_DLLM_wrapper_TLift(self):
        config_dict = self.__get_base_config_dict()
        config_dict['test.DLLM.type']='TargetLift'
        config_dict['test.DLLM.target_Lift']=769200.
        
        DLLMWrap = DLLMWrapper('test', verbose=0)
        DLLMWrap.configure(config_dict)
        DLLMWrap.set_out_format('numpy')
        DLLMWrap.set_grad_format('numpy')
        
        DLLMWrap.analysis()
        Lift=DLLMWrap.get_F_value('Lift')
        assert((Lift-769200.)<1.e-2)
        
    def test_DLLM_wrapper_TLift_valid_grad(self):
        config_dict = self.__get_base_config_dict()
        config_dict['test.DLLM.type']='TargetLift'
        config_dict['test.DLLM.target_Lift']=769200.
        
        DLLMWrap = DLLMWrapper('test', verbose=0)
        DLLMWrap.configure(config_dict)
        DLLMWrap.set_out_format('numpy')
        DLLMWrap.set_grad_format('numpy')
        
        x0=DLLMWrap.get_x0()
         
        DLLMWrap.analysis()
         
        Ref_list=DLLMWrap.get_F_list()
         
        def f(x):
            DLLMWrap.run(x)
            func=DLLMWrap.get_F_list()
            return func/Ref_list
           
        def df(x):
            DLLMWrap.run_grad(x)
            func_grad=np.array(DLLMWrap.get_F_list_grad())
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
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLMWrapper)
    unittest.TextTestRunner(verbosity=2).run(suite)