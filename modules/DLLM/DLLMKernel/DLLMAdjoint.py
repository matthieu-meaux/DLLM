# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: Francois Gallard
# @author: Matthieu MEAUX (for refactoring)

from numpy import dot, zeros
from numpy.linalg import norm, solve

from MDOTools.Solvers.adjoint_problem import AdjointProblem

class DLLMAdjoint:
    def __init__(self, LLW):
        """
        Adjoint module for the lifting line model
        """
        self.__LLW = LLW
        self.__computed = False

        self.__adj_list           = None
        self.__adj_conv_corr_list = None
        self.__dJ_dchi_list       = None
        
    #-- accessors
    def get_adjoint_list(self):
        return self.__adj_list
    
    def get_adjoint_convergence_correction_list(self):
        return self.__adj_conv_corr_list
    
    def get_dJ_dchi_list(self):
        return self.__dJ_dchi_list
        
    #-- computed related methods
    def is_computed(self):
        return self.__computed
    
    def set_computed(self, bool=True):
        self.__computed = bool
    
    #-- Residual accessors
    def get_R(self):
        return self.__LLW.get_R()
    
    def get_DR_DiAoA(self):
        return self.__LLW.get_DR_DiAoA()
    
    def get_DR_Dchi(self):
        return self.__LLW.get_DR_Dchi()
    
    #-- Post-processing accessors
    def get_func_list(self):
        return self.__LLW.get_func_list()
    
    def get_dfunc_diAoA(self):
        return self.__LLW.get_dfunc_diAoA()
    
    def get_dpJ_dpchi(self):
        return self.__LLW.get_dpJ_dpchi()
    
    #-- Run adjoint problem for each cost function
    def run(self):
        func_list = self.get_func_list()
        
        self.__adj_list           = []
        self.__adj_conv_corr_list = []
        self.__dJ_dchi_list       = []
        
        R        = self.get_R()
        dpRdpW   = self.get_DR_DiAoA()
        dpRdpchi = self.get_DR_Dchi()
        
        dpFdpW_list   = self.get_dfunc_diAoA()
        dpFdpchi_list = self.get_dpJ_dpchi()
        
        for i, func in enumerate(func_list):
           print 'Run adjoint problem for func = '+str(func)
           dpFdpW   = dpFdpW_list[i]
           dpFdpchi = dpFdpchi_list[i]
           
           AdjPb = AdjointProblem(R, dpRdpW, dpRdpchi, dpFdpW, dpFdpchi)
           AdjPb.solve()
           
           self.__adj_list.append(AdjPb.get_adjoint_state())
           self.__adj_conv_corr_list.append(AdjPb.get_convergence_correction())
           self.__dJ_dchi_list.append(AdjPb.get_dFdchi())
           
           print '  - Convergence adjoint correction for '+str(func)+' = '+str(self.__adj_conv_corr_list[i])
           