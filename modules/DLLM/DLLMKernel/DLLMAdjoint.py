# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: Francois Gallard
# @author: Matthieu MEAUX (for refactoring)

from numpy import dot, zeros
from numpy.linalg import norm, solve

class DLLMAdjoint:
    def __init__(self, LLW):
        """
        Adjoint module for the lifting line model
        """
        self.__LLW = LLW
        self.__computed = False
        self.__N        = self.get_wing_param().get_n_sect()
        self.__ndv      = self.get_wing_param().get_ndv()
        self.__adj_list = None
        self.__adj_corr = None
        self.__dJ_dchi  = None
        
    #-- accessors
    def get_adj_list(self):
        return self.__adj_list
    
    def get_adj_corr(self):
        return self.__adj_corr
    
    def get_dJ_dchi(self):
        return self.__dJ_dchi
        
    #-- computed related methods
    def is_computed(self):
        return self.__computed
    
    def set_computed(self, bool=True):
        self.__computed = bool
        
    #-- Accessors
    def get_wing_param(self):
        return self.__LLW.get_wing_param()
    
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
    
    # Adjoint methods
    def adjoint(self,i):
        '''
        Computes the adjoint of the problem
        @param dJ_diAoA : sensibility of the objective function to the state vector : induced aoa
        '''
        DR_DiAoA    = self.get_DR_DiAoA()
        dfunc_diAoA = self.get_dfunc_diAoA()
        dJ_diAoA    = dfunc_diAoA[i]
        adj = solve(DR_DiAoA.T,-dJ_diAoA.T)
        return adj
    
    def __adjoint_correction(self,i):
        '''
        Computes the adjoint correction
        @param adj : adjoint vector
        '''
        adj=self.__adj_list[i]
        corr = dot(adj.T,self.get_R())
        return corr
    
    def run(self):
        func_list = self.get_func_list()
        n_func    = len(func_list)
        self.__adj_list  = zeros((n_func,self.__N))
        self.__adj_corr  = zeros(n_func)
        self.__dJ_dchi   = zeros((n_func,self.__ndv))

        for i, func in enumerate(func_list):
            print 'Run adjoint for func = '+str(func)
            self.__adj_list[i,:] = self.adjoint(i)[:]
            self.__adj_corr[i]   = self.__adjoint_correction(i)
            print '  - Convergence adjoint correction for '+str(func)+' = '+str(self.__adj_corr[i])
            self.__dJ_dchi[i,:]  = self.__DJ_Dchi(i)[:]
        print self.__adj_list
        print self.__adj_corr
        print self.__dJ_dchi
        
    def __DJ_Dchi(self,i):
        dR_dchi   = self.get_DR_Dchi()
        dpJ_dpchi = self.get_dpJ_dpchi()[i]
        adjoint   = self.__adj_list[i]
        DJ_Dchi   = dpJ_dpchi + dot(adjoint.T,dR_dchi)
        return DJ_Dchi
   