'''
Created on Jun 4, 2015

@author: Francois
'''

from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent
from DLLM.DLLMKernel.DLLMPost import DLLMPost
from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
from MDOTools.OC.operating_condition import OperatingCondition
import unittest
import numpy as np

class TestDLLM_wrapper(unittest.TestCase):
    
      
    def test_init_component(self):
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N=5, verbose=0)
        default_lift = 0.00
        self.assertEqual(DLLMOpenMDAO.Lift, default_lift)
   
    def test_change_OC(self):
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N = 5, verbose=0)
        self.assertEqual(0., DLLMOpenMDAO.OC.get_altitude())
        DLLMOpenMDAO.altitude=0.
        DLLMOpenMDAO.Mach = 0.
        DLLMOpenMDAO.T0 = OperatingCondition.T0+10
        DLLMOpenMDAO.P0 = OperatingCondition.P0+10
        
        DLLMOpenMDAO.execute()
        self.assertEqual(0., DLLMOpenMDAO.OC.get_altitude())
        self.assertEqual(0., DLLMOpenMDAO.OC.get_Mach())
        self.assertEqual(OperatingCondition.T0+10, DLLMOpenMDAO.OC.get_T0())
        self.assertEqual(OperatingCondition.P0+10, DLLMOpenMDAO.OC.get_P0())
  
    def test_change_dv(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        new_span = 48.0
        DLLMOpenMDAO.span= new_span
        DLLMOpenMDAO.execute()
        self.assertEqual(new_span, DLLMOpenMDAO.wing_param.get_span())
  
    def test_exec(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        Target_Lift =  606570.049598     
        DLLMOpenMDAO.DLLM.set_target_Lift(Target_Lift)
        DLLMOpenMDAO.execute()
#         print "Final lift:",DLLMOpenMDAO.Lift
        self.assertAlmostEqual(Target_Lift, DLLMOpenMDAO.Lift,places=3)
          
    def test_list_deriv_vars(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        Target_Lift =  606570.049598     
        DLLMOpenMDAO.DLLM.set_target_Lift(Target_Lift)
        DLLMOpenMDAO.execute()
        self.assertAlmostEqual(Target_Lift, DLLMOpenMDAO.Lift,places=3)
        out_dvid_tuple, F_list_names_tuple = DLLMOpenMDAO.list_deriv_vars()
        self.assertTrue(type(out_dvid_tuple)==tuple)
        self.assertTrue(type(F_list_names_tuple)==tuple)
        F_list_names_ref = tuple(DLLMPost.DEF_F_LIST_NAMES)
        self.assertTupleEqual(F_list_names_ref, F_list_names_tuple)

    def test_jacobian(self):
        N = 5
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N, verbose=0)
        Target_Lift =  606570.049598     
        DLLMOpenMDAO.DLLM.set_target_Lift(Target_Lift)
        DLLMOpenMDAO.execute()
        out_dvid_tuple, F_list_names_tuple = DLLMOpenMDAO.list_deriv_vars()
        n_F = len(F_list_names_tuple)
        Jacobian = DLLMOpenMDAO.provideJ()
        n_dv = len(DLLMOpenMDAO.wing_param.get_dv_id_list())
        self.assertEqual(n_dv, Jacobian.shape[0])
        self.assertEqual(n_F, Jacobian.shape[1])
        
        def update_dllm_from_x(x):
            DLLMOpenMDAO.rtwist = x[:5]
            DLLMOpenMDAO.span = x[5]
            DLLMOpenMDAO.sweep = x[6]
            DLLMOpenMDAO.break_percent = x[7]
            DLLMOpenMDAO.root_chord = x[8]
            DLLMOpenMDAO.break_chord = x[9]
            DLLMOpenMDAO.tip_chord = x[10]
            DLLMOpenMDAO.root_height = x[11]
            DLLMOpenMDAO.break_height = x[12]
            DLLMOpenMDAO.tip_height=x[13]
            assert(len(x)==14)
        
        def d_obj(x):
            update_dllm_from_x(x)
            DLLMOpenMDAO.execute()
            return DLLMOpenMDAO.provideJ()
            
        def obj(x):
            update_dllm_from_x(x)
            DLLMOpenMDAO.execute()
            return DLLMOpenMDAO.DLLM.get_F_list()
        
        treshold=1e-4
        val_grad=FDValidGrad(1,obj,d_obj,fd_step=1.e-7)
        x0=DLLMOpenMDAO.wing_param.get_dv_array()
        ok,df_fd,df=val_grad.compare(x0,treshold=treshold,split_out=True,return_all=True, iprint =False)
        varnames, funcnames=DLLMOpenMDAO.list_deriv_vars()
        varnames=np.concatenate((['tw0','tw1','tw2','tw3','tw4'],np.array(varnames[1:])))
        funcnames=np.array(funcnames)
        k=0
        ok = True
        for fd, gr in zip(df_fd,df):
            err=np.where(((fd-gr)/fd>treshold) & (fd>1e-5))[0]
            if err.size>0:
                #Compute absolute error when gradient or fd is lower than threshold
                inf_tol=err[np.where((np.abs(gr[err])<treshold) | (np.abs(fd[err])<treshold))[0]]
                #Scale gradient by function value for the lift issue in particular
                f_val=getattr(DLLMOpenMDAO,funcnames[k])
                actually_ok=np.where(((fd[inf_tol]-gr[inf_tol])<treshold*f_val))[0]
                real_err=np.delete(err,actually_ok)
                if real_err.size >0:
                    print "Errror in d",funcnames[k]," /d ",varnames[real_err]
                    print "df analytic = ",gr[real_err]
                    print "df finite d = ",fd[real_err]
                    ok=False
            k=k+1
        assert(ok)
        #This is not that reliable... threshold is hard coded at 1e-5, only relative error is computed...
        #out = DLLMOpenMDAO.check_gradient(fd_form='central', fd_step=1e-5, fd_step_type='relative')
        #assert(len(out[-1])==0)
        


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLM_wrapper)
    unittest.TextTestRunner(verbosity=2).run(suite)
