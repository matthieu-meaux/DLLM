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
        self.__dF_list_dchi       = None
        
    #-- accessors
    def get_adjoint_list(self):
        return self.__adj_list
    
    def get_adjoint_convergence_correction_list(self):
        return self.__adj_conv_corr_list
    
    def get_dF_list_dchi(self):
        return self.__dF_list_dchi
        
    #-- computed related methods
    def is_computed(self):
        return self.__computed
    
    def set_computed(self, bool=True):
        self.__computed = bool
    
    #-- Residual accessors
    def get_R(self):
        return self.__LLW.get_R()
    
    def get_dpR_dpW(self):
        return self.__LLW.get_dpR_dpW()
    
    def get_dpR_dpchi(self):
        return self.__LLW.get_dpR_dpchi()
    
    #-- Post-processing accessors
    def get_F_list_names(self):
        return self.__LLW.get_F_list_names()
    
    def get_dpF_list_dpW(self):
        return self.__LLW.get_dpF_list_dpW()
    
    def get_dpF_list_dpchi(self):
        return self.__LLW.get_dpF_list_dpchi()
    
    #-- Run adjoint problem for each cost function
    def run(self):
        F_list_names = self.get_F_list_names()
        
        self.__adj_list           = []
        self.__adj_conv_corr_list = []
        self.__dF_list_dchi       = []
        
        R        = self.get_R()
        dpRdpW   = self.get_dpR_dpW()
        dpRdpchi = self.get_dpR_dpchi()
        
        dpFdpW_list   = self.get_dpF_list_dpW()
        dpFdpchi_list = self.get_dpF_list_dpchi()
        
        for i, F_name in enumerate(F_list_names):
           print 'Run adjoint problem for func = '+str(F_name)
           dpFdpW   = dpFdpW_list[i]
           dpFdpchi = dpFdpchi_list[i]
           
           AdjPb = AdjointProblem(R, dpRdpW, dpRdpchi, dpFdpW, dpFdpchi)
           AdjPb.solve()
           
           self.__adj_list.append(AdjPb.get_adjoint_state())
           self.__adj_conv_corr_list.append(AdjPb.get_convergence_correction())
           self.__dF_list_dchi.append(AdjPb.get_dFdchi())
           
           print '  - Convergence adjoint correction for '+str(F_name)+' = '+str(self.__adj_conv_corr_list[i])
           