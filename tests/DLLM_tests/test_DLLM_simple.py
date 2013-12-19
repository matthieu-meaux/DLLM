import unittest
from numpy import zeros, array

from OpenDACE.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition

class TestDLLMSimple(unittest.TestCase):
    
    def __init_wing_param(self):
        OC=OperatingCondition('cond1',atmospheric_model='simple')
        OC.set_Mach(0.8)
        OC.set_AoA(3.5)
        OC.set_altitude(10000.)
        OC.set_T0_deg(15.)
        OC.set_P0(101325.)
        OC.set_humidity(0.)
        OC.compute_atmosphere()
        
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=40)
        wing_param.build_wing()
        wing_param.set_value('test_param.span',34.1)
        wing_param.set_value('test_param.sweep',34.)
        wing_param.set_value('test_param.break_percent',33.)
        wing_param.set_value('test_param.root_chord',6.1)
        wing_param.set_value('test_param.break_chord',4.6)
        wing_param.set_value('test_param.tip_chord',1.5)
        wing_param.set_value('test_param.root_height',1.28)
        wing_param.set_value('test_param.break_height',0.97)
        wing_param.set_value('test_param.tip_height',0.33)
        wing_param.convert_to_design_variable('test_param.span',10.,50.)
        wing_param.convert_to_design_variable('test_param.sweep',0.,40.)
        wing_param.convert_to_design_variable('test_param.break_percent',20.,40.)
        wing_param.convert_to_design_variable('test_param.root_chord',5.,7.)
        wing_param.convert_to_design_variable('test_param.break_chord',3.,5.)
        wing_param.convert_to_design_variable('test_param.tip_chord',1.,2.)
        wing_param.convert_to_design_variable('test_param.root_height',1.,1.5)
        wing_param.convert_to_design_variable('test_param.break_height',0.8,1.2)
        wing_param.convert_to_design_variable('test_param.tip_height',0.2,0.5)
        wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        wing_param.build_airfoils_from_ref()
        wing_param.update()
        
        return OC,wing_param

    def test_DLLM_instantiation(self):
        """
        test class instantiation
        """
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        assert(DLLM is not None)
        
    def test_DLLM_run_direct(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        try:
            print ''
            DLLM.run_direct()
            ok=True
        except:
            ok=False
        assert(ok)
        
    def test_DLLM_run_direct_post(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        try:
            print ''
            DLLM.run_direct()
            DLLM.run_post()
            ok=True
        except:
            ok=False
        assert(ok)
    
    def test_DLLM_run_direct_post_adjoint(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        try:
            print ''
            DLLM.run_direct()
            DLLM.run_post()
            DLLM.run_adjoint()
            ok=True
        except:
            ok=False
        assert(ok)
        
    def test_DLLM_valid_dpR_dpiAoA(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA0=DLLM.get_iAoA()
        def f1(x):
            func=DLLM.comp_R(x)
            return func
        
        def df1(x):
            func_grad=DLLM.comp_dpR_dpiAoA(x)
            return func_grad
        
        val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-3)
        ok1,df_fd1,df1=val_grad1.compare(iAoA0,treshold=1.e-2,return_all=True)
        assert(ok1)
        
    def test_DLLM_valid_dpR_dpchi(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA=DLLM.get_iAoA()
        x0=wing_param.get_dv_array()
        def f2(x):
            wing_param.update_from_x_list(x)
            DLLM.set_wing_param(wing_param)
            func=DLLM.comp_R(iAoA)
            return func
        
        def df2(x):
            wing_param.update_from_x_list(x)
            DLLM.set_wing_param(wing_param)
            func=DLLM.comp_R(iAoA)
            func_grad=DLLM.comp_dpR_dpchi()
            return func_grad
        
        val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-3)
        ok2,df_fd2,df2=val_grad2.compare(x0,treshold=1.e-2,return_all=True)
        assert(ok2)
        
    def test_DLLM_valid_dpR_dpthetaY(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA=DLLM.get_iAoA()
        thetaY0=wing_param.get_thetaY()
        def f3(x):
            wing_param.set_thetaY(x)
            func=DLLM.comp_R(iAoA)
            return func
        
        def df3(x):
            wing_param.set_thetaY(x)
            func_grad=DLLM.comp_dpR_dpthetaY()
            return func_grad
        
        val_grad3=FDValidGrad(2,f3,df3,fd_step=1.e-3)
        ok3,df_fd3,df3=val_grad3.compare(thetaY0,treshold=1.e-2,return_all=True)
        assert(ok3)
        
    def test_DLLM_valid_dpF_list_dpiAoA(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA0=DLLM.get_iAoA()
        def f4(x):
            DLLM.comp_R(x)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func
        
        def df4(x):
            DLLM.comp_R(x)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func_grad=DLLM.get_dpF_list_dpiAoA()
            return func_grad
        
        val_grad4=FDValidGrad(2,f4,df4,fd_step=1.e-3)
        ok4,df_fd4,df4=val_grad4.compare(iAoA0,treshold=1.e-2,return_all=True)
        assert(ok4)
        
    def test_DLLM_valid_dpF_list_dpchi(self):
        OC,wing_param = self.__init_wing_param()
        DLLM = DLLMSolver(wing_param,OC)
        print ''
        DLLM.run_direct()
        iAoA=DLLM.get_iAoA()
        x0=wing_param.get_dv_array()
        def f5(x):
            wing_param.update_from_x_list(x)
            DLLM.set_wing_param(wing_param)
            DLLM.comp_R(iAoA)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func
        
        def df5(x):
            wing_param.update_from_x_list(x)
            DLLM.set_wing_param(wing_param)
            DLLM.comp_R(iAoA)
            DLLM.set_direct_computed()
            DLLM.run_post()
            func_grad=DLLM.get_dpF_list_dpchi()
            return func_grad
        
        val_grad5=FDValidGrad(2,f5,df5,fd_step=1.e-3)
        ok5,df_fd5,df5=val_grad5.compare(x0,treshold=1.e-2,return_all=True)
        assert(ok5)
        
    def test_DLLM_valid_adjoint_gradients(self):
        OC,wing_param = self.__init_wing_param()
        x0=wing_param.get_dv_array()
        def f6(x):
            wing_param.update_from_x_list(x)
            DLLM = DLLMSolver(wing_param,OC)
            DLLM.run_direct()
            DLLM.run_post()
            func=DLLM.get_F_list()
            return func
        
        def df6(x):
            wing_param.update_from_x_list(x)
            DLLM = DLLMSolver(wing_param,OC)
            DLLM.run_direct()
            DLLM.run_post()
            DLLM.run_adjoint()
            func_grad=array(DLLM.get_dF_list_dchi())
            return func_grad

        val_grad6=FDValidGrad(2,f6,df6,fd_step=1.e-3)
        ok6,df_fd6,df6=val_grad6.compare(x0,treshold=1.e-2,return_all=True)
        assert(ok6)
                
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLMSimple)
    unittest.TextTestRunner(verbosity=2).run(suite)