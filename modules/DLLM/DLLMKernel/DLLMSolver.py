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
# @author : Francois Gallard
# @author : Matthieu MEAUX
#

# - Local imports -
from DLLM.DLLMKernel.DLLMMesh import DLLMMesh
from DLLM.DLLMKernel.DLLMDirect import DLLMDirect
from DLLM.DLLMKernel.DLLMPost import DLLMPost
from DLLM.DLLMKernel.DLLMAdjoint import DLLMAdjoint

class DLLMSolver:
    ERROR_MSG='ERROR in DLLMSolver.'
    def __init__(self, tag, geom, OC, verbose = 0, grad_active=True):
        '''
        Constructor for wings based on lifting line theory
        @param wing_geom : the wing geometry
        @param Lref : Reference length for the moments computation
        @param relaxFactor : relaxation factor for induced angles computation
        @param stopCriteria : the stop criteria for the iterative method for computing the wing circulation.
        '''
        self.__tag         = tag

        self.__verbose     = verbose
        self.__grad_active = grad_active 

        
        self.__geom        = geom
        self.__OC          = OC
        
        self.__Lref        = 0.
        self.__Sref        = 0.
        self.__Lref_grad   = None
        self.__Sref_grad   = None
        
        self.__DLLMMesh    = DLLMMesh(self, verbose = self.__verbose)
        self.__DLLMDirect  = DLLMDirect(self, verbose = self.__verbose)
        self.__DLLMPost    = DLLMPost(self, verbose = self.__verbose)
        if grad_active:
            self.__DLLMAdjoint = DLLMAdjoint(self, verbose = self.__verbose)

        self.set_OC(OC)
        
    
    #-- Accessors
    def get_tag(self):
        return self.__tag
    
    def get_grad_active(self):
        return self.__grad_active
    
    def get_geom(self):
        return self.__geom
    
    def get_airfoils(self):
        return self.get_geom().get_linked_airfoils()
    
    def get_OC(self):
        return self.__OC
    
    def get_Lref(self):
        return self.__geom.get_Lref()
    
    def get_Lref_grad(self):
        return self.__geom.get_Lref_grad()
    
    def get_Sref(self):
        return self.__geom.get_Sref()
    
    def get_Sref_grad(self):
        return self.__geom.get_Sref_grad()
    
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

    def get_dplocalAoA_dpthetaY(self):
        return self.__DLLMDirect.get_dplocalAoA_dpthetaY()
    
    def comp_R(self, iAoA):
        return self.__DLLMDirect.comp_R(iAoA)
    
    def comp_dpR_dpiAoA(self, iAoA):
        return self.__DLLMDirect.comp_dpR_dpiAoA(iAoA)
    
    def comp_dpR_dpchi(self):
        return self.__DLLMDirect.comp_dpR_dpchi()
    
    def comp_dpR_dpthetaY(self):
        return self.__DLLMDirect.comp_dpR_dpthetaY()
    
    def comp_dpR_dpAoA(self):
        return self.__DLLMDirect.comp_dpR_dpAoA()
    
    #-- DLLMPost accessors
    def is_post_computed(self):
        return self.__DLLMPost.is_computed()
    
    def get_F_list_names(self):
        return self.__DLLMPost.get_F_list_names()
    
    def get_F_list(self):
        return self.__DLLMPost.get_F_list()
    
    def get_dpF_list_dpW(self):
        return self.__DLLMPost.get_dpF_list_dpiAoA()
    
    def get_dpF_list_dpAoA(self):
        return self.__DLLMPost.get_dpF_list_dpAoA()
    
    def get_dpF_list_dpchi(self):
        return self.__DLLMPost.get_dpF_list_dpchi()

    def get_dpF_list_dpthetaY(self):
        return self.__DLLMPost.get_dpF_list_dpthetaY()

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
        if self.get_grad_active():
            self.__DLLMAdjoint.set_computed(False)
        
    def set_OC(self, OC):
        self.__OC = OC
        self.__reinit_modules()
        
    def set_geom(self, geom):
        self.__geom  = geom
        self.__DLLMMesh.recompute()
        self.__reinit_modules()
        
    def set_F_list_names(self, F_list_names):
        self.__DLLMPost.set_F_list_names(F_list_names)
        
    #-- DLLMDirect setters
    def set_relax_factor(self, relax_factor):
        self.__DLLMDirect.set_relax_factor(relax_factor)
        
    def set_stop_residual(self, residual):
        self.__DLLMDirect.set_stop_residual(residual)
    
    def set_max_iterations(self, max_it):
        self.__DLLMDirect.set_max_iterations(max_it)
        
    def set_method(self, method):
        self.__DLLMDirect.set_method(method)
        
    def set_gamma_file_name(self, gamma_f_name):
        self.__DLLMDirect.set_gamma_file_name(gamma_f_name)
        
    #-- Run methods
    def run_direct(self):
        self.__DLLMDirect.run()
        
    def run_post(self, F_list_names=None):
        ERROR_MSG=self.ERROR_MSG+'run_post: '
        if self.is_direct_computed():
            self.__DLLMPost.run(F_list_names=F_list_names)
        else:
            print ERROR_MSG+'Cannot run post-processing if solution is not computed'
            
    def run_adjoint(self):
        ERROR_MSG=self.ERROR_MSG+'run_adjoint: '
        if self.get_grad_active():
            if self.is_post_computed():
                self.__DLLMAdjoint.run()
            else:
                print ERROR_MSG+'Cannot run adjoint if post-processing is not computed'
        else:
            print ERROR_MSG+'Cannot run adjoint if gradient is not active'
            
    #-- Export methods
    def plot(self):
        if self.__DLLMDirect.is_computed():
            self.__DLLMDirect.plot()
        if self.__DLLMPost.is_computed():
            self.__DLLMPost.plot()
    
    def export_F_list(self, filename=None):
        self.__DLLMPost.export_F_list(filename=filename)
        
    def export_dF_list_dchi(self, filename=None):
        self.__DLLMAdjoint.export_dF_list_dchi(filename=filename)
