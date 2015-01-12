import unittest
from numpy import zeros, array

from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMGeom.wing_param import Wing_param
from MDOTools.OC.operating_condition import OperatingCondition

class TestWingParam(unittest.TestCase):
    
    def __init_OC(self):
        OC=OperatingCondition('cond1')
        OC.set_Mach(0.7)
        OC.set_AoA(3.0)
        OC.set_altitude(5000.)
        OC.set_T0_deg(20.)
        OC.set_P0(101325.)
        OC.set_humidity(0.)
        OC.compute_atmosphere()
        return OC
    
    def __init_wing_param(self):
        OC=OperatingCondition('cond1')
        OC.set_Mach(0.7)
        OC.set_AoA(3.0)
        OC.set_altitude(5000.)
        OC.set_T0_deg(20.)
        OC.set_P0(101325.)
        OC.set_humidity(0.)
        OC.compute_atmosphere()
        
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        wing_param.build_wing()
        wing_param.set_value('test_param.span',34.1)
        wing_param.set_value('test_param.sweep',32.)
        wing_param.set_value('test_param.break_percent',33.)
        wing_param.set_value('test_param.root_chord',6.1)
        wing_param.set_value('test_param.break_chord',4.6)
        wing_param.set_value('test_param.tip_chord',1.5)
        wing_param.set_value('test_param.root_height',1.28)
        wing_param.set_value('test_param.break_height',0.97)
        wing_param.set_value('test_param.tip_height',0.33)
        wing_param.convert_to_design_variable('test_param.span',(10.,50.))
        wing_param.convert_to_design_variable('test_param.sweep',(0.,40.))
        wing_param.convert_to_design_variable('test_param.break_percent',(20.,40.))
        wing_param.convert_to_design_variable('test_param.root_chord',(5.,7.))
        wing_param.convert_to_design_variable('test_param.break_chord',(3.,5.))
        wing_param.convert_to_design_variable('test_param.tip_chord',(1.,2.))
        wing_param.convert_to_design_variable('test_param.root_height',(1.,1.5))
        wing_param.convert_to_design_variable('test_param.break_height',(0.8,1.2))
        wing_param.convert_to_design_variable('test_param.tip_height',(0.2,0.5))
        wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        wing_param.build_airfoils_from_ref()
        wing_param.update()
        
        return wing_param

    def test_Wing_param_instantiation(self):
        """
        test class instantiation
        """
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        assert(wing_param is not None)
        
    def test_Wing_param_build_wing(self):
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        try:
            wing_param.build_wing()
            ok=True
        except:
            ok=False
        assert(ok)
        
    def test_Wing_param_set_value(self):
        """
        test set value for Wing_param
        """
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        wing_param.build_wing()
        try:
            wing_param.set_value('test_param.span',34.1)
            wing_param.set_value('test_param.sweep',32.)
            wing_param.set_value('test_param.break_percent',33.)
            wing_param.set_value('test_param.root_chord',6.1)
            wing_param.set_value('test_param.break_chord',4.6)
            wing_param.set_value('test_param.tip_chord',1.5)
            wing_param.set_value('test_param.root_height',1.28)
            wing_param.set_value('test_param.break_height',0.97)
            wing_param.set_value('test_param.tip_height',0.33)
            ok=True
        except:
            ok=False
        assert(ok)
        
    def test_Wing_param_convert(self):
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        wing_param.build_wing()
        wing_param.set_value('test_param.span',34.1)
        wing_param.set_value('test_param.sweep',32.)
        wing_param.set_value('test_param.break_percent',33.)
        wing_param.set_value('test_param.root_chord',6.1)
        wing_param.set_value('test_param.break_chord',4.6)
        wing_param.set_value('test_param.tip_chord',1.5)
        wing_param.set_value('test_param.root_height',1.28)
        wing_param.set_value('test_param.break_height',0.97)
        wing_param.set_value('test_param.tip_height',0.33)
        try:
            wing_param.convert_to_design_variable('test_param.span',(10.,50.))
            wing_param.convert_to_design_variable('test_param.sweep',(0.,40.))
            wing_param.convert_to_design_variable('test_param.break_percent',(20.,40.))
            wing_param.convert_to_design_variable('test_param.root_chord',(5.,7.))
            wing_param.convert_to_design_variable('test_param.break_chord',(3.,5.))
            wing_param.convert_to_design_variable('test_param.tip_chord',(1.,2.))
            wing_param.convert_to_design_variable('test_param.root_height',(1.,1.5))
            wing_param.convert_to_design_variable('test_param.break_height',(0.8,1.2))
            wing_param.convert_to_design_variable('test_param.tip_height',(0.2,0.5))
            ok=True
        except:
            ok=False
        assert(ok)
        
    def test_Wing_param_linear_airfoil(self):
        OC=self.__init_OC()
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        wing_param.build_wing()
        wing_param.set_value('test_param.span',34.1)
        wing_param.set_value('test_param.sweep',32.)
        wing_param.set_value('test_param.break_percent',33.)
        wing_param.set_value('test_param.root_chord',6.1)
        wing_param.set_value('test_param.break_chord',4.6)
        wing_param.set_value('test_param.tip_chord',1.5)
        wing_param.set_value('test_param.root_height',1.28)
        wing_param.set_value('test_param.break_height',0.97)
        wing_param.set_value('test_param.tip_height',0.33)
        wing_param.convert_to_design_variable('test_param.span',(10.,50.))
        wing_param.convert_to_design_variable('test_param.sweep',(0.,40.))
        wing_param.convert_to_design_variable('test_param.break_percent',(20.,40.))
        wing_param.convert_to_design_variable('test_param.root_chord',(5.,7.))
        wing_param.convert_to_design_variable('test_param.break_chord',(3.,5.))
        wing_param.convert_to_design_variable('test_param.tip_chord',(1.,2.))
        wing_param.convert_to_design_variable('test_param.root_height',(1.,1.5))
        wing_param.convert_to_design_variable('test_param.break_height',(0.8,1.2))
        wing_param.convert_to_design_variable('test_param.tip_height',(0.2,0.5))
        try:
            wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
            wing_param.build_airfoils_from_ref()
            ok=True
        except:
            ok=False
        assert(ok)
        
    def test_Wing_param_update(self):
        OC=self.__init_OC()
        wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
        wing_param.build_wing()
        wing_param.set_value('test_param.span',34.1)
        wing_param.set_value('test_param.sweep',32.)
        wing_param.set_value('test_param.break_percent',33.)
        wing_param.set_value('test_param.root_chord',6.1)
        wing_param.set_value('test_param.break_chord',4.6)
        wing_param.set_value('test_param.tip_chord',1.5)
        wing_param.set_value('test_param.root_height',1.28)
        wing_param.set_value('test_param.break_height',0.97)
        wing_param.set_value('test_param.tip_height',0.33)
        wing_param.convert_to_design_variable('test_param.span',(10.,50.))
        wing_param.convert_to_design_variable('test_param.sweep',(0.,40.))
        wing_param.convert_to_design_variable('test_param.break_percent',(20.,40.))
        wing_param.convert_to_design_variable('test_param.root_chord',(5.,7.))
        wing_param.convert_to_design_variable('test_param.break_chord',(3.,5.))
        wing_param.convert_to_design_variable('test_param.tip_chord',(1.,2.))
        wing_param.convert_to_design_variable('test_param.root_height',(1.,1.5))
        wing_param.convert_to_design_variable('test_param.break_height',(0.8,1.2))
        wing_param.convert_to_design_variable('test_param.tip_height',(0.2,0.5))
        wing_param.build_linear_airfoil(OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        wing_param.build_airfoils_from_ref()
        try:
            wing_param.update()
            ok=True
        except:
            ok=False
        assert(ok)
        
    def test_Wing_param_valid_grad_twist(self):
        wing_param=self.__init_wing_param()
        x0=wing_param.get_dv_array()
        def f1(x):
            wing_param.update_from_x_list(x)
            func=wing_param.get_twist()
            return func
        
        def df1(x):
            wing_param.update_from_x_list(x)
            func_grad=wing_param.get_twist_grad()
            return func_grad
        
        val_grad1=FDValidGrad(2,f1,df1,fd_step=1.e-8)
        ok1,df_fd1,df1=val_grad1.compare(x0,treshold=1.e-6,return_all=True)
        assert(ok1)
        
    def test_Wing_param_valid_grad_chords(self):
        wing_param=self.__init_wing_param()
        x0=wing_param.get_dv_array()
        def f2(x):
            wing_param.update_from_x_list(x)
            func=wing_param.get_chords()
            return func
        
        def df2(x):
            wing_param.update_from_x_list(x)
            func_grad=wing_param.get_chords_grad()
            return func_grad
        
        val_grad2=FDValidGrad(2,f2,df2,fd_step=1.e-8)
        ok2,df_fd2,df2=val_grad2.compare(x0,treshold=1.e-6,return_all=True)
        assert(ok2)
        
    def test_Wing_param_valid_grad_rel_thicks(self):
        wing_param=self.__init_wing_param()
        x0=wing_param.get_dv_array()
        def f3(x):
            wing_param.update_from_x_list(x)
            func=wing_param.get_rel_thicks()
            return func
        
        def df3(x):
            wing_param.update_from_x_list(x)
            func_grad=wing_param.get_rel_thicks_grad()
            return func_grad
        
        val_grad3=FDValidGrad(2,f3,df3,fd_step=1.e-8)
        ok3,df_fd3,df3=val_grad3.compare(x0,treshold=1.e-6,return_all=True)
        assert(ok3)
        
    def test_Wing_param_valid_grad_eta(self):
        wing_param=self.__init_wing_param()
        x0=wing_param.get_dv_array()
        def f4(x):
            wing_param.update_from_x_list(x)
            func=wing_param.get_eta()
            return func
        
        def df4(x):
            wing_param.update_from_x_list(x)
            func_grad=wing_param.get_eta_grad()
            return func_grad
        
        val_grad4=FDValidGrad(2,f4,df4,fd_step=1.e-8)
        ok4,df_fd4,df4=val_grad4.compare(x0,treshold=1.e-6,return_all=True)
        assert(ok4)
        
    def test_Wing_param_valid_grad_XYZ(self):
        wing_param=self.__init_wing_param()
        x0=wing_param.get_dv_array()
        def f5(x):
            wing_param.update_from_x_list(x)
            func=wing_param.get_XYZ()
            return func
        
        def df5(x):
            wing_param.update_from_x_list(x)
            func_grad=wing_param.get_XYZ_grad()
            return func_grad
        
        val_grad5=FDValidGrad(2,f5,df5,fd_step=1.e-8)
        ok5,df_fd5,df5=val_grad5.compare(x0,treshold=1.e-6,return_all=True)
        assert(ok5)
        
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestWingParam)
    unittest.TextTestRunner(verbosity=2).run(suite)