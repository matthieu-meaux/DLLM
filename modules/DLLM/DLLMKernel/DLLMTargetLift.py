# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus group innovation
# @version: 1.0
# @author: Matthieu MEAUX
import numpy
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.Solvers.newton_raphson_problem import NewtonRaphsonProblem

class DLLMTargetLift(DLLMSolver):
    ERROR_MSG='ERROR in DLLMTargetLift.'
    def __init__(self, tag, wing_param, OC):
        DLLMSolver.__init__(self, tag, wing_param, OC)
        
        self.__N            = self.get_wing_param().get_n_sect()
        self.__ndv          = self.get_wing_param().get_ndv()
        
        self.__R_TL         = numpy.zeros(self.__N+1)
        self.__dpR_TL_dpW   = numpy.zeros((self.__N+1,self.__N+1))
        self.__dpR_TL_dpchi = numpy.zeros((self.__N+1,self.__ndv))
        
        self.__dpF_list_dpW = None 
        
        self.__target_Lift    = 10000.

        # initialize the Newton-Raphson problem
        self.__NRPb = None 
        self.__init_Newton_Raphson()
        
    # -- Accessors
    def get_target_Lift(self):
        return self.__target_Lift
    
    def set_target_Lift(self, target):
        self.__target_Lift   = target
        
    #-- Newton-Raphson related methods
    def __init_Newton_Raphson(self):
        W0= numpy.zeros(self.__N+1)
        self.__NRPb = NewtonRaphsonProblem(W0, self.comp_R_TL, self.comp_dpR_TL_dpW)
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
    def comp_R_TL(self, W):
        DLLMDirect = self.get_DLLMDirect()
        DLLMPost   = self.get_DLLMPost()
        OC = self.get_OC()
        
        iAoA = W[0:self.__N]
        AoA  = W[self.__N]
        OC.set_AoA_rad(AoA)
        
        R    = DLLMDirect.comp_R(iAoA)
        Lift = DLLMPost.comp_Lift()

        self.__R_TL[:self.__N] = R[:]
        self.__R_TL[self.__N]  = Lift - self.__target_Lift
        
        return self.__R_TL
    
    def comp_dpR_TL_dpW(self, W):
        DLLMDirect = self.get_DLLMDirect()
        DLLMPost   = self.get_DLLMPost()
        OC =  self.get_OC()
        
        iAoA = W[0:self.__N]
        AoA  = W[self.__N]
        OC.set_AoA_rad(AoA)
        
        dpR_dpiAoA    = DLLMDirect.comp_dpR_dpiAoA(iAoA)
        dpLift_dpiAoA = DLLMPost.comp_dpLift_dpiAoA()
        
        dpR_dpAoA   = DLLMDirect.comp_dpR_dpAoA()
        dpLift_dpAoA  = DLLMPost.dpLift_dpAoA()        
        
        self.__dpR_TL_dpW[0:self.__N,0:self.__N] = dpR_dpiAoA[:,:]
        self.__dpR_TL_dpW[self.__N,0:self.__N]   = dpLift_dpiAoA[:]
        
        self.__dpR_TL_dpW[0:self.__N,self.__N]   = dpR_dpAoA[:]
        self.__dpR_TL_dpW[self.__N,self.__N]     = dpLift_dpAoA
        
        return self.__dpR_TL_dpW
    
    def get_R(self):
        return self.__R_TL
    
    def get_dpR_dpW(self):
        return self.__dpR_TL_dpW
    
    def get_dpR_dpchi(self):
        DLLMDirect = self.get_DLLMDirect()
        DLLMPost   = self.get_DLLMPost()
        dpR_dpchi    = DLLMDirect.get_dpR_dpchi()
        dpLift_dpchi = DLLMPost.dpLift_dpchi()
        
        self.__dpR_TL_dpchi[:self.__N,:] = dpR_dpchi[:,:]
        self.__dpR_TL_dpchi[self.__N,:]  = dpLift_dpchi[:]
        
        return self.__dpR_TL_dpchi
    
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
        
        
        

        
        
        
        
        

        