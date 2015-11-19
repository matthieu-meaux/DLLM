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
# @author : Matthieu MEAUX
#

import numpy
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.Solvers.newton_raphson_problem import NewtonRaphsonProblem

class DLLMTargetCl(DLLMSolver):
    ERROR_MSG='ERROR in DLLMTargetCl.'
    def __init__(self, tag,  geom, OC, verbose = 0):
        self.__verbose = verbose
        DLLMSolver.__init__(self, tag, geom, OC, verbose = self.__verbose)
        
        self.__N            = self.get_geom().get_n_sect()
        self.__ndv          = self.get_geom().get_ndv()
        
        self.__R_TCl         = numpy.zeros(self.__N+1)
        self.__dpR_TCl_dpW   = numpy.zeros((self.__N+1,self.__N+1))
        self.__dpR_TCl_dpchi = numpy.zeros((self.__N+1,self.__ndv))
        
        self.__dpF_list_dpW = None 
        
        self.__target_Cl    = 0.5

        # initialize the Newton-Raphson problem
        self.__NRPb = None 
        self.__init_Newton_Raphson()
        
    # -- Accessors
    def get_target_Cl(self):
        return self.__target_Cl
    
    def set_target_Cl(self, target):
        self.__target_Cl   = target
        
    #-- Newton-Raphson related methods
    def __init_Newton_Raphson(self):
        W0= numpy.zeros(self.__N+1)
        self.__NRPb = NewtonRaphsonProblem(W0, self.comp_R_TCl, self.comp_dpR_TCl_dpW, verbose = self.__verbose)
        self.__NRPb.set_relax_factor(0.99)
        self.__NRPb.set_stop_residual(1.e-9)
        self.__NRPb.set_max_iterations(100)
        
    def set_relax_factor(self, relax_factor):
        self.__NRPb.set_relax_factor(relax_factor)
    
    def set_stop_residual(self, residual):
        self.__NRPb.set_stop_residual(residual)
        
    def set_max_iterations(self, max_it):
        self.__NRPb.set_max_iterations(max_it)
        
    def set_method(self, method):
        self.__NRPb.set_method(method)
        
    #-- Residual Related method
    def comp_R_TCl(self, W):
        DLLMDirect = self.get_DLLMDirect()
        DLLMPost   = self.get_DLLMPost()
        OC = self.get_OC()
        
        iAoA = W[0:self.__N]
        AoA  = W[self.__N]
        OC.set_AoA_rad(AoA)
        
        R    = DLLMDirect.comp_R(iAoA)
        Cl   = DLLMPost.comp_Cl()

        self.__R_TCl[:self.__N] = R[:]
        self.__R_TCl[self.__N]  = Cl - self.__target_Cl
        
        return self.__R_TCl
    
    def comp_dpR_TCl_dpW(self, W):
        DLLMDirect = self.get_DLLMDirect()
        DLLMPost   = self.get_DLLMPost()
        OC =  self.get_OC()
        
        iAoA = W[0:self.__N]
        AoA  = W[self.__N]
        OC.set_AoA_rad(AoA)
        
        dpR_dpiAoA  = DLLMDirect.comp_dpR_dpiAoA(iAoA)
        dpCl_dpiAoA = DLLMPost.comp_dpCl_dpiAoA()
        
        dpR_dpAoA   = DLLMDirect.comp_dpR_dpAoA()
        dpCl_dpAoA  = DLLMPost.dpCl_dpAoA()        
        
        self.__dpR_TCl_dpW[0:self.__N,0:self.__N] = dpR_dpiAoA[:,:]
        self.__dpR_TCl_dpW[self.__N,0:self.__N]   = dpCl_dpiAoA[:]
        
        self.__dpR_TCl_dpW[0:self.__N,self.__N]   = dpR_dpAoA[:]
        self.__dpR_TCl_dpW[self.__N,self.__N]     = dpCl_dpAoA
        
        return self.__dpR_TCl_dpW
    
    def get_R(self):
        return self.__R_TCl
    
    def get_dpR_dpW(self):
        return self.__dpR_TCl_dpW
    
    def get_dpR_dpchi(self):
        DLLMDirect = self.get_DLLMDirect()
        DLLMPost   = self.get_DLLMPost()
        dpR_dpchi  = DLLMDirect.get_dpR_dpchi()
        dpCl_dpchi = DLLMPost.dpCl_dpchi()
        
        self.__dpR_TCl_dpchi[:self.__N,:] = dpR_dpchi[:,:]
        self.__dpR_TCl_dpchi[self.__N,:]  = dpCl_dpchi[:]
        
        return self.__dpR_TCl_dpchi
    
    def get_dpF_list_dpW(self):
        DLLMPost   = self.get_DLLMPost()
        
        dpF_list_dpiAoA = DLLMPost.get_dpF_list_dpiAoA()
        dpF_list_dpAoA  = DLLMPost.get_dpF_list_dpAoA()
        
        dim = len(dpF_list_dpiAoA)
        
        self.__dpF_list_dpW = numpy.zeros((dim,self.__N+1))
        
        for i in xrange(dim):
            self.__dpF_list_dpW[i,:self.__N] = dpF_list_dpiAoA[i,:]
            self.__dpF_list_dpW[i,self.__N]  = dpF_list_dpAoA[i]
        
        return self.__dpF_list_dpW
    
    
    def run_direct(self):
        DLLMDirect = self.get_DLLMDirect()
#        W0= numpy.zeros(self.__N+1)
#        self.__NRPb.valid_jacobian(W0, iprint=True)
        self.__NRPb.solve()
        DLLMDirect.set_computed(True)
        DLLMDirect.write_gamma_to_file()
        DLLMDirect.comp_dpR_dpchi()
            