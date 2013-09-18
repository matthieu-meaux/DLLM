# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: Francois Gallard
# @author: Matthieu MEAUX (for refactoring)

from numpy import dot
from numpy.linalg import norm, solve

class DLLMAdjoint:
    def __init__(self, LLW):
        """
        Adjoint module for the lifting line model
        """
        self.__LLW = LLW
        self.__computed = False
        self.__N        = self.get_wing_geom().get_n_elements()
        self.__adj_list = None
        self.__adj_corr = None
        self.__dJ_dTwist = None
        self.__dJ_dAoA   = None
        self.__dJ_dThick = None
        
    #-- computed related methods
    def is_computed(self):
        return self.__computed
    
    def set_computed(self, bool=True):
        self.__computed = bool
        
    #-- Accessors
    def get_wing_geom(self):
        return self.__LLW.get_wing_geom()
    
    #-- Residual accessors
    def get_R(self):
        return self.__LLW.get_R()
    
    def get_DR_DiAoA(self):
        return self.__LLW.get_DR_DiAoA()
    
    def get_DR_DTwist(self):
        return self.__LLW.get_DR_DTwist()
    
    def get_DR_DAoA(self):
        return self.__LLW.get_DR_DAoA()
    
    def get_DR_DThickness(self):
        return self.__LLW.get_DR_DThickness()
    
    #-- Post-processing accessors
    def get_func_list(self):
        return self.__LLW.get_func_list()
    
    def get_dfunc_diAoA(self):
        return self.__LLW.get_dfunc_diAoA()
    
    def get_dpJ_dpTwist(self):
        return self.__LLW.get_dpJ_dpTwist()
    
    def get_dpJ_dpAoA(self):
        return self.__LLW.get_dpJ_dpAoA()
    
    def get_dpJ_dpThickness(self):
        return self.__LLW.get_dpJ_dpThickness()
    
    # Adjoint methods
    def adjoint(self,dJ_diAoA):
        '''
        Computes the adjoint of the problem
        @param dJ_diAoA : sensibility of the objective function to the state vector : induced aoa
        '''
        DR_DiAoA=self.get_DR_DiAoA()
        adj = solve(DR_DiAoA.T,-dJ_diAoA.T)
        return adj
    
    def __adjoint_correction(self,adj):
        '''
        Computes the adjoint correction
        @param adj : adjoint vector
        '''
        corr = dot(adj.T,self.get_R())
        return corr
    
    def run(self):
        n_func=len(self.get_func_list())
        self.__adj_list  = [None]*n_func
        self.__adj_corr  = [None]*n_func
        self.__dJ_dTwist = [None]*n_func
        self.__dJ_dAoA   = [None]*n_func
        self.__dJ_dThick = [None]*n_func
        dpJ_dpTwist = self.get_dpJ_dpTwist()
        dpJ_dpAoA   = self.get_dpJ_dpAoA()
        dpJ_dpThick = self.get_dpJ_dpThickness()
        for i, func in enumerate(self.get_func_list()):
            print 'Run adjoint for func = '+str(func)
            self.__adj_list[i]  = self.adjoint(self.get_dfunc_diAoA()[i])
            self.__adj_corr[i]  = self.__adjoint_correction(self.__adj_list[i])
            print '  - Convergence adjoint correction for '+str(func)+' = '+str(self.__adj_corr[i])
            self.__dJ_dTwist[i] = self.__DJ_DTwist(dpJ_dpTwist[i], self.__adj_list[i])
            self.__dJ_dAoA[i]   = self.__DJ_DAoA(dpJ_dpAoA[i], self.__adj_list[i])
            self.__dJ_dThick[i] = self.__DJ_DThick(dpJ_dpThick[i], self.__adj_list[i])
            
    def __DJ_DTwist(self,dpJ_dpTwist,adjoint):
        dR_dTwist=self.get_DR_DTwist()
        DJ_DTwist=dpJ_dpTwist+dot(adjoint.T,dR_dTwist)
        return DJ_DTwist
    
    def __DJ_DAoA(self,dpJ_dpAoA,adjoint):
        dR_dAoA=self.get_DR_DAoA()
        DJ_DAoA=dpJ_dpAoA+dot(adjoint.T,dR_dAoA)
        return DJ_DAoA

    def __DJ_DThick(self,dpJ_dpThickness,adjoint):
        dR_dThick=self.get_DR_DThickness()
        DJ_DThick=dpJ_dpThickness+dot(adjoint.T,dR_dThick)
        return DJ_DThick
    