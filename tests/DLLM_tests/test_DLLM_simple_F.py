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
from numpy import zeros, array

from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMGeom.wing_broken import Wing_Broken
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition

class TestDLLMSimpleF(unittest.TestCase):
    
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
        
    def test_DLLM_valid_dpF_list_dpW(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver('test',wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA0=DLLM.get_iAoA()
        def f1(x):
            DLLM.comp_R(x)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func
        
        def df1(x):
            DLLM.comp_R(x)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func_grad=DLLM.get_dpF_list_dpW()
            return func_grad
        
        val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-8)
        ok1,df_fd1,df1=val_grad1.compare(iAoA0,treshold=1.e-6,return_all=True)
        assert(ok1)
        
    def test_DLLM_valid_dpF_list_dpchi(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver('test',wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA=DLLM.get_iAoA()
        x0=wing_param.get_dv_array()
        def f2(x):
            wing_param.update_from_x_list(x)
            DLLM.set_geom(wing_param)
            DLLM.comp_R(iAoA)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func
        
        def df2(x):
            wing_param.update_from_x_list(x)
            DLLM.set_geom(wing_param)
            DLLM.comp_R(iAoA)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func_grad=DLLM.get_dpF_list_dpchi()
            return func_grad
        
        val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-8)
        ok2,df_fd2,df2=val_grad2.compare(x0,treshold=1.e-6,split_out=True,return_all=True)
        assert(ok2)
        
    def test_DLLM_valid_dpF_list_dpAoA(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver('test',wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA=DLLM.get_iAoA()
        AoA0=OC.get_AoA_rad()
        
        def f3(x):
            OC.set_AoA_rad(x[0])
            DLLM.comp_R(iAoA)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func
        
        def df3(x):
            OC.set_AoA_rad(x[0])
            DLLM.comp_R(iAoA)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func_grad=DLLM.get_dpF_list_dpAoA()
            N=len(func_grad)
            np_func_grad=zeros((N,1))
            np_func_grad[:,0]=func_grad[:]
            return np_func_grad
        
        val_grad3=FDValidGrad(2,f3,df3,fd_step=1.e-8)
        ok3,df_fd3,df3=val_grad3.compare([AoA0],treshold=1.e-6,split_out=True,return_all=True)
        assert(ok3)
        
    def test_DLLM_valid_dF_list_dchi(self):
        OC,wing_param = self.__init_wing_param()
        x0=wing_param.get_dv_array()
        def f4(x):
            wing_param.update_from_x_list(x)
            DLLM = DLLMSolver('test',wing_param,OC)
            DLLM.run_direct()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func
        
        def df4(x):
            wing_param.update_from_x_list(x)
            DLLM = DLLMSolver('test',wing_param,OC)
            DLLM.run_direct()
            DLLM.run_post()
            DLLM.run_adjoint()
            func_grad=array(DLLM.get_dF_list_dchi())
            return func_grad

        val_grad4=FDValidGrad(2,f4,df4,fd_step=1.e-8)
        ok4,df_fd4,df4=val_grad4.compare(x0,treshold=1.e-6,split_out=True,return_all=True)
        assert(ok4)
        
           
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLMSimpleF)
    unittest.TextTestRunner(verbosity=2).run(suite)