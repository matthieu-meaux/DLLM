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
# @author : Louis Blanchard

import unittest
from numpy import zeros

from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition
from numpy.linalg import norm
import matplotlib.pyplot as plt
import numpy as np

class TestDLLMSimpleLoads(unittest.TestCase):

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init_wing_param(self):
        OC=OperatingCondition('cond1')
        OC.set_Mach(0.8)
        OC.set_AoA(3.5)
        OC.set_altitude(10000.)
        OC.set_T0_deg(15.)
        OC.set_P0(101325.)
        OC.set_humidity(0.)
        OC.compute_atmosphere()
        
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        wing_param.build_wing()
        wing_param.set_value('span',34.1)
        wing_param.set_value('sweep',34.)
        wing_param.set_value('break_percent',33.)
        wing_param.set_value('root_chord',6.1)
        wing_param.set_value('break_chord',4.6)
        wing_param.set_value('tip_chord',1.5)
        wing_param.set_value('root_height',1.28)
        wing_param.set_value('break_height',0.97)
        wing_param.set_value('tip_height',0.33)
        wing_param.convert_to_design_variable('span',(10.,50.))
        wing_param.convert_to_design_variable('sweep',(0.,40.))
        wing_param.convert_to_design_variable('break_percent',(20.,40.))
        wing_param.convert_to_design_variable('root_chord',(5.,7.))
        wing_param.convert_to_design_variable('break_chord',(3.,5.))
        wing_param.convert_to_design_variable('tip_chord',(1.,2.))
        wing_param.convert_to_design_variable('root_height',(1.,1.5))
        wing_param.convert_to_design_variable('break_height',(0.8,1.2))
        wing_param.convert_to_design_variable('tip_height',(0.2,0.5))
        wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        wing_param.build_airfoils_from_ref()
        wing_param.update()
        
        return OC,wing_param

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    def plot_valid_FD(self,X, name1,df1,df1_fd,name2,df2,df2_fd, norm1, norm2, name):

        M1_per_error = 100*abs(df1_fd-df1)/abs(df1_fd)
        M2_per_error = 100*abs(df2_fd-df2)/abs(df2_fd)
        M1_per_error = M1_per_error[~np.isnan(M1_per_error)].T
        M2_per_error = M2_per_error[~np.isnan(M2_per_error)].T
        fig =plt.figure(figsize=(15,10))
        nb=10e1
        x1=np.linspace(1,nb,num= len(M1_per_error))
        x2=np.linspace(1,nb,num= len(M2_per_error))
        plt.semilogy(x1,M1_per_error,'ro-',markersize=5)
        plt.semilogy(x2,M2_per_error,'bo-',markersize=5)
        plt.grid()
        plt.ylabel("y")
        plt.legend(['percentage error : '+name1,'percentage error : '+name2 ], loc='best')
        txt=" Jacobian matrix validation by FD : norm error : "+name1+" = %.3e, "+name2+" = %.3e"
        plt.title(txt % (norm1,norm2) )
        plt.savefig("DLLM_Jacobian_validation_FD_"+name+".png")


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def print_valid_FD(self,X,name1,df1,df1_fd, name2,df2,df2_fd):
        diff1  = abs(df1_fd-df1).sum()
        diff2  = abs(df2_fd-df2).sum()
        norm1  = norm(df1_fd-df1)/norm(df1_fd)
        norm2  = norm(df2_fd-df2)/norm(df2_fd)
        print  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print  "            Finite Difference Validation for Jacobian matrix : "
        print  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print  " X = ",X
        print  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print  name1,": diff = %.5e  norm = %.5e" % (diff1,norm1)
        print  name2,": diff = %.5e  norm = %.5e" % (diff2,norm2)
        print  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        return norm1, norm2

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def test_DLLM_valid_dpLoads_distrib_dpiAoA(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver('test',wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA0=DLLM.get_iAoA()
        Post=DLLM.get_DLLMPost()

        def f1(x):
            R=DLLM.comp_R(x)
            func=Post.comp_Lift_distrib()
            return func
        
        def df1(x):
            R=DLLM.comp_R(x)
            func_grad=Post.comp_dpLift_distrib_dpiAoA()
            return func_grad

        def f2(x):
            R=DLLM.comp_R(x)
            func=Post.comp_Drag_distrib()
            return func
        
        def df2(x):
            R=DLLM.comp_R(x)
            func_grad=Post.comp_dpDrag_distrib_dpiAoA()
            return func_grad
        
        val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-8)
        ok1,df_fd1,df1=val_grad1.compare(iAoA0,treshold=1.e-6,return_all=True)
       	val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-8)
        ok2,df_fd2,df2=val_grad2.compare(iAoA0,treshold=1.e-6,return_all=True)
    	norm1, norm2 = self.print_valid_FD(iAoA0,'dpLift_distrib_dpiAoA',df1,df_fd1,'dpDrag_distrib_dpiAoA',df2,df_fd2)
    	self.plot_valid_FD(iAoA0,'dpLift_distrib_dpiAoA',df1,df_fd1,'dpDrag_distrib_dpiAoA',df2,df_fd2, norm1, norm2, 'Loads_distrib_dpiAoA')
    	if (ok1==True) and (ok2==True): 
            ok = True
    	else : 	
            ok = False
        assert(ok)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    def test_DLLM_valid_dpLoads_distrib_dpAoA(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver('test',wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA0=DLLM.get_iAoA()
        AoA0=OC.get_AoA_rad()
        Post=DLLM.get_DLLMPost()

        def f1(x):
            OC.set_AoA_rad(x[0])
            R=DLLM.comp_R(iAoA0)
            func=Post.comp_Lift_distrib()
            return func
        
        def df1(x):
            OC.set_AoA_rad(x[0])
            R=DLLM.comp_R(iAoA0)
            func_grad=Post.comp_dpLift_distrib_dpAoA()
            N=len(func_grad)
            np_func_grad=np.zeros((N,1))
            np_func_grad[:,0]=func_grad[:]
            return np_func_grad

     	def f2(x):
            OC.set_AoA_rad(x[0])
            R=DLLM.comp_R(iAoA0)
            func=Post.comp_Drag_distrib()
            return func
        
        def df2(x):
            OC.set_AoA_rad(x[0])
            R=DLLM.comp_R(iAoA0)
            func_grad=Post.comp_dpDrag_distrib_dpAoA()
            N=len(func_grad)
            np_func_grad=np.zeros((N,1))
            np_func_grad[:,0]=func_grad[:]
            return np_func_grad
        
        val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-8)
        ok1,df_fd1,df1=val_grad1.compare([AoA0],treshold=1.e-6,return_all=True)
        val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-8)
        ok2,df_fd2,df2=val_grad2.compare([AoA0],treshold=1.e-6,return_all=True)
        norm1, norm2 = self.print_valid_FD([AoA0],'dpLift_distrib_dpAoA',df1,df_fd1,'dpDrag_distrib_dpAoA',df2,df_fd2)
        self.plot_valid_FD([AoA0],'dpLift_distrib_dpAoA',df1,df_fd1,'dpDrag_distrib_dpAoA',df2,df_fd2, norm1, norm2, 'Loads_distrib_dpAoA')
        if (ok1==True) and (ok2==True): 
            ok = True
        else : 	
            ok = False
        assert(ok)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    def test_DLLM_valid_dpLoads_distrib_dpchi(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver('test',wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA=DLLM.get_iAoA()
        x0=wing_param.get_dv_array()
    	Post=DLLM.get_DLLMPost()

        def f1(x):
            wing_param.update_from_x_list(x)
            DLLM.set_geom(wing_param)
            DLLM.comp_R(iAoA)
            func=Post.comp_Lift_distrib()
            return func
        
        def df1(x):
            wing_param.update_from_x_list(x)
            DLLM.set_geom(wing_param)
            DLLM.comp_R(iAoA)
            func_grad=Post.comp_dpLift_distrib_dpchi()
            return func_grad
        
        def f2(x):
            wing_param.update_from_x_list(x)
            DLLM.set_geom(wing_param)
            DLLM.comp_R(iAoA)
            func=Post.comp_Drag_distrib()
            return func
        
        def df2(x):
            wing_param.update_from_x_list(x)
            DLLM.set_geom(wing_param)
            DLLM.comp_R(iAoA)
            func_grad=Post.comp_dpDrag_distrib_dpchi()
            return func_grad

        val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-8)
        ok1,df_fd1,df1=val_grad1.compare(x0,treshold=1.e-6,return_all=True)
        val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-8)
        ok2,df_fd2,df2=val_grad2.compare(x0,treshold=1.e-6,return_all=True)
    	norm1, norm2 = self.print_valid_FD(x0,'dpLift_distrib_dpchi',df1,df_fd1,'dpDrag_distrib_dpchi',df2,df_fd2)
    	self.plot_valid_FD(x0,'dpLift_distrib_dpchi',df1,df_fd1,'dpDrag_distrib_dpchi',df2,df_fd2, norm1, norm2, 'Loads_distrib_dpchi')
    	if (ok1==True) and (ok2==True):
            ok = True
    	else : 
            ok = False
    	assert(ok)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    def test_DLLM_valid_dpLoads_distrib_dpthetaY(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver('test',wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA=DLLM.get_iAoA()
        thetaY0=wing_param.get_thetaY()
    	Post=DLLM.get_DLLMPost()

        def f1(x):
            wing_param.set_thetaY(x)
            DLLM.comp_R(iAoA)
	    func=Post.comp_Lift_distrib()
            return func
        
        def df1(x):
            wing_param.set_thetaY(x)
            func_grad=Post.comp_dpLift_distrib_dpthetaY()
            return func_grad

        def f2(x):
            wing_param.set_thetaY(x)
            DLLM.comp_R(iAoA)
	    func=Post.comp_Drag_distrib()
            return func
        
        def df2(x):
            wing_param.set_thetaY(x)
            func_grad=Post.comp_dpDrag_distrib_dpthetaY()
            return func_grad
        
        val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-8)
        ok1,df_fd1,df1=val_grad1.compare(thetaY0,treshold=1.e-6,return_all=True)
    	val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-8)
        ok2,df_fd2,df2=val_grad2.compare(thetaY0,treshold=1.e-6,return_all=True)
    	norm1, norm2 = self.print_valid_FD(thetaY0,'dpLift_distrib_dpthetaY',df1,df_fd1,'dpDrag_distrib_dpthetaY',df2,df_fd2)
    	self.plot_valid_FD(thetaY0,'dpLift_distrib_dpthetaY',df1,df_fd1,'dpDrag_distrib_dpthetaY',df2,df_fd2, norm1, norm2, 'Loads_distrib_dpthetaY')
    	if (ok1==True) and (ok2==True):
            ok = True
    	else : 
            ok = False
        assert(ok)
           
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLMSimpleLoads)
    unittest.TextTestRunner(verbosity=2).run(suite)
