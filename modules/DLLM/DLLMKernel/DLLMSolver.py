# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard
# @author: Matthieu MEAUX (for refactoring)
 
# - Local imports -
from DLLM.DLLMKernel.DLLMMesh import DLLMMesh
from DLLM.DLLMKernel.DLLMDirect import DLLMDirect
from DLLM.DLLMKernel.DLLMPost import DLLMPost
from DLLM.DLLMKernel.DLLMAdjoint import DLLMAdjoint

class DLLMSolver:
    ERROR_MSG='ERROR in DLLMSolver.'
    def __init__(self, wing_param, OC):
        '''
        Constructor for wings based on lifting line theory
        @param wing_geom : the wing geometry
        @param Lref : Reference length for the moments computation
        @param relaxFactor : relaxation factor for induced angles computation
        @param stopCriteria : the stop criteria for the iterative method for computing the wing circulation.
        '''
        self.__wing_param = wing_param
        self.__OC         = OC
        
        self.__Lref       = 0.
        self.__Sref       = 0.
        self.__Lref_grad  = None
        self.__Sref_grad  = None
        
        self.__DLLMMesh   = DLLMMesh(self)
        self.__DLLMDirect = DLLMDirect(self)
        self.__DLLMPost   = DLLMPost(self)
        self.__DLLMAdjoint= DLLMAdjoint(self)
        
        self.set_OC(OC)
    
    #-- Accessors
    def get_wing_param(self):
        return self.__wing_param
    
    def get_airfoils(self):
        return self.get_wing_param().get_linked_airfoils()
    
    def get_OC(self):
        return self.__OC
    
    def get_Lref(self):
        return self.__Lref
    
    def get_Lref_grad(self):
        return self.__Lref_grad
    
    def get_Sref(self):
        return self.__Sref
    
    def get_Sref_grad(self):
        return self.__Sref_grad
    
    def get_DLLMMesh(self):
        return self.__DLLMMesh
    
    def get_DLLMDirect(self):
        return self.__DLLMDirect
    
    def get_DLLMPost(self):
        return self.__DLLMPost
    
    def get_DLLMAdjoint(self):
        return self.__DLLMAdjoint
    
    #-- DLLMMesh accessors
    def get_K(self):
        return self.__DLLMMesh.get_K()
    
    def get_dK_dchi(self):
        return self.__DLLMMesh.get_dK_dchi()
    
    #-- DLLMDirect accessors
    def set_direct_computed(self, bool=True):
        self.__DLLMDirect.set_computed(bool)
    
    def is_direct_computed(self):
        return self.__DLLMDirect.is_computed()
    
    def get_convergence_history(self):
        return self.__DLLMDirect.get_convergence_history()
    
    def get_localAoA(self):
        return self.__DLLMDirect.get_localAoA()
    
    def get_iAoA(self):
        return self.__DLLMDirect.get_iAoA()
    
    def get_R(self):
        return self.__DLLMDirect.get_R()
    
    def get_dpR_dpW(self):
        return self.__DLLMDirect.get_dpR_dpiAoA()
    
    def get_dpR_dpchi(self):
        return self.__DLLMDirect.get_dpR_dpchi()
    
    def get_dplocalAoA_dpiAoA(self):
        return self.__DLLMDirect.get_dplocalAoA_dpiAoA()
    
    def get_dplocalAoA_dpAoA(self):
        return self.__DLLMDirect.get_dplocalAoA_dpAoA()
    
    def get_dplocalAoA_dpchi(self):
        return self.__DLLMDirect.get_dplocalAoA_dpchi()
    
    def comp_R(self, iAoA):
        return self.__DLLMDirect.comp_R(iAoA)
    
    def comp_dpR_dpiAoA(self, iAoA):
        return self.__DLLMDirect.comp_dpR_dpiAoA(iAoA)
    
    def comp_dpR_dpchi(self):
        return self.__DLLMDirect.comp_dpR_dpchi()
    
    def comp_dpR_dpthetaY(self):
        return self.__DLLMDirect.comp_dpR_dpthetaY()
    
    #-- DLLMPost accessors
    def is_post_computed(self):
        return self.__DLLMPost.is_computed()
    
    def get_F_list_names(self):
        return self.__DLLMPost.get_F_list_names()
    
    def get_F_list(self):
        return self.__DLLMPost.get_F_list()
    
    def get_dpF_list_dpW(self):
        return self.__DLLMPost.get_dpF_list_dpiAoA()
    
    def get_dpF_list_dpchi(self):
        return self.__DLLMPost.get_dpF_list_dpchi()
    
    #-- DLLMAdjoint accessors
    def get_adjoint_list(self):
        return self.__DLLMAdjoint.get_adjoint_list()
    
    def get_adjoint_convergence_correction_list(self):
        return self.__DLLMAdjoint.get_adjoint_convergence_correction_list()
    
    def get_dF_list_dchi(self):
        return self.__DLLMAdjoint.get_dF_list_dchi()
    
    #-- Setters
    def __reinit_modules(self):
        self.__DLLMDirect.set_computed(False)
        self.__DLLMPost.set_computed(False)
        self.__DLLMAdjoint.set_computed(False)
        
    def set_OC(self, OC):
        self.__OC = OC
        self.__reinit_modules()
        
    def set_wing_param(self, wing_param):
        self.__wing_param  = wing_param
        self.__DLLMMesh.recompute()
        self.__reinit_modules()
        
    def set_airfoil(self, airfoils):
        self.__airfoils = airfoils
        self.__DLLMMesh.recompute()
        self.__reinit_modules()
        
    def set_Lref(self, Lref):
        self.__Lref = Lref
        
    def set_Lref_grad(self, Lref_grad):
        self.__Lref_grad = Lref_grad
        
    def set_Sref(self, Sref):
        self.__Sref = Sref
        
    def set_Sref_grad(self, Sref_grad):
        self.__Sref_grad = Sref_grad
    
    #-- DLLMDirect setters
    def set_relax_factor(self, relax_factor):
        self.__DLLMDirect.set_relax_factor(relax_factor)
        
    def set_stop_residual(self, residual):
        self.__DLLMDirect.set_stop_residual(residual)
    
    def set_max_iterations(self, max_it):
        self.__DLLMDirect.set_max_iterations(max_it)
        
    def set_gamma_file_name(self, gamma_f_name):
        self.__DLLMDirect.set_gamma_file_name(gamma_f_name)
        
    #-- Run methods
    def run_direct(self):
        self.__DLLMDirect.run()
        
    def run_post(self,F_list_names=None):
        ERROR_MSG=self.ERROR_MSG+'run_post: '
        if self.is_direct_computed():
            self.__DLLMPost.run(F_list_names=F_list_names)
        else:
            print ERROR_MSG+'Cannot run post-processing if solution is not computed'
            
    def run_adjoint(self):
        ERROR_MSG=self.ERROR_MSG+'run_adjoint: '
        if self.is_post_computed():
            self.__DLLMAdjoint.run()
        else:
            print ERROR_MSG+'Cannot run adjoint if post-processing is not computed'       
    
