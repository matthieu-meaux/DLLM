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

from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMGeom.wing_broken import Wing_Broken
from DLLM.DLLMKernel.DLLMTargetCl import DLLMTargetCl
from DLLM.DLLMKernel.DLLMTargetLift import DLLMTargetLift
from MDOTools.OC.operating_condition import OperatingCondition

class TestDLLMSimpleTClTLift(unittest.TestCase):
    
    def __init_wing_param(self):
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
        
        return OC,wing_param
    
    def test_DLLM_valid_TCl(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMTargetCl('Simple',wing_param,OC)
#         F_list=DLLM.get_F_list()
        F_list_names=DLLM.get_F_list_names()
        DLLM.set_target_Cl(0.5)
        DLLM.run_direct()
        DLLM.run_post()
        F_list=DLLM.get_F_list()
        F_list_names=DLLM.get_F_list_names()
        cl_index = F_list_names.index('Cl')
        Cl=F_list[cl_index]
        assert(abs(Cl-0.5)<1.e-10)
          
    def test_DLLM_valid_grad_TCl(self):
        OC,wing_param = self.__init_wing_param()
        x0=wing_param.get_dv_array()
        
        DLLM = DLLMTargetCl('Simple',wing_param,OC)
        DLLM.set_target_Cl(0.5)
        DLLM.run_direct()
        DLLM.run_post()
        Ref_list=DLLM.get_F_list()
        
        def f1(x):
            wing_param.update_from_x_list(x)
            DLLM = DLLMTargetCl('Simple',wing_param,OC)
            DLLM.set_target_Cl(0.5)
            DLLM.run_direct()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func/Ref_list
           
        def df1(x):
            wing_param.update_from_x_list(x)
            DLLM = DLLMTargetCl('Simple',wing_param,OC)
            DLLM.set_target_Cl(0.5)
            DLLM.run_direct()
            DLLM.run_post()
            DLLM.run_adjoint()
            func_grad=np.array(DLLM.get_dF_list_dchi())
            N = func_grad.shape[0]
            ndv = func_grad.shape[1]
            out_grad = np.zeros((N,ndv))
            for i in xrange(N):
                out_grad[i,:] = func_grad[i,:]/Ref_list[i]
            return out_grad
           
        val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-8)
        ok1,df_fd1,df1=val_grad1.compare(x0,treshold=1.e-5,split_out=True,return_all=True)
        assert(ok1)
        
    def test_DLLM_valid_TLift(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMTargetLift('Simple',wing_param,OC)
#         F_list=DLLM.get_F_list()
        F_list_names=DLLM.get_F_list_names()
        DLLM.set_target_Lift(769200.)
        DLLM.run_direct()
        DLLM.run_post()
        F_list=DLLM.get_F_list()
        F_list_names=DLLM.get_F_list_names()
        print F_list_names
        lift_index = F_list_names.index('Lift')
        Lift=F_list[lift_index]
        assert(abs(Lift-769200.)<1.e-2)
        
    def test_DLLM_valid_grad_TLift(self):
        OC,wing_param = self.__init_wing_param()
        x0=wing_param.get_dv_array()
        
        DLLM = DLLMTargetLift('Simple',wing_param,OC)
        DLLM.set_target_Lift(769200.)
        DLLM.run_direct()
        DLLM.run_post()
        Ref_list=DLLM.get_F_list()
  
        def f2(x):
            wing_param.update_from_x_list(x)
            DLLM = DLLMTargetLift('Simple',wing_param,OC)
            DLLM.set_target_Lift(769200.)
            DLLM.run_direct()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func/Ref_list
          
        def df2(x):
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
                out_grad[i,:] = func_grad[i,:]/Ref_list[i]
            return out_grad
      
        val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-8)
        ok2,df_fd2,df2=val_grad2.compare(x0,treshold=1.e-5,split_out=True,return_all=True)
        assert(ok2)
        
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLMSimpleTClTLift)
    unittest.TextTestRunner(verbosity=2).run(suite)