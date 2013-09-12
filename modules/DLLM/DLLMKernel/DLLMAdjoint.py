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
        
    #-- computed related methods
    def is_computed(self):
        return self.__computed
    
    def set_computed(self, bool=True):
        self.__computed = bool
        
    #-- Accessors
    def get_wing_geom(self):
        return self.__LLW.get_wing_geom()
    
    def get_R(self):
        return self.__LLW.get_R()
    
    def get_DR_DiAoA(self):
        return self.__LLW.get_DR_DiAoA()
        
    def get_func_list(self):
        return self.__LLW.get_func_list()
    
    def get_dfunc_diAoA(self):
        return self.__LLW.get_dfunc_diAoA()
    
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
        self.__adj_list = [None]*len(self.get_func_list())
        for i, func in enumerate(self.get_func_list()):
            print 'Run adjoint for func = '+str(func)
            self.__adj_list[i]=self.adjoint(self.get_dfunc_diAoA()[i])
            corr=self.__adjoint_correction(self.__adj_list[i])
            print '  - Convergence adjoint correction for '+str(func)+' = '+str(corr)
            
    def DJ_DTwist(self,dJ_dTwist,adjoint,alpha,Mach):
        '''
        Builds the gradient of J with the adjoint field
        @param dJ_dTwist : partial derivative of the objective function to the twist
        @param adjoint : the adjoint field
        @param alpha : the wing angle of attack
        '''
        self.__Iterate(alpha,Mach)
        dR=self.dRDTwist(self.__iAoA,alpha,Mach)
        
        return dJ_dTwist+dot(adjoint.T,dR)
    
    def DJ_DAoA(self,dJ_dAoA,adjoint,alpha,Mach):
        '''
        Builds the gradient of J with the adjoint field
        @param dJ_dTwist : partial derivative of the objective function to the twist
        @param adjoint : the adjoint field
        @param alpha : the wing angle of attack
        '''
        self.__Iterate(alpha,Mach)
        dR=self.dRDAoA(self.__iAoA,alpha,Mach)
        
        return dJ_dAoA+dot(adjoint.T,dR)

    def DJ_DThick(self,dJ_dThickness,adjoint,alpha,Mach):
        '''
        Builds the gradient of J with the adjoint field
        @param dJ_dThickness : partial derivative of the objective function to the twist
        @param adjoint : the adjoint field
        @param alpha : the wing angle of attack
        '''
        self.__Iterate(alpha,Mach)
        dR=self.dRDThickness(self.__iAoA,alpha,Mach)
        
        return dJ_dThickness+dot(adjoint.T,dR)
    