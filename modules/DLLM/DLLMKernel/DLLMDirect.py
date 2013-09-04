# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: Francois Gallard
# @author: Matthieu MEAUX (for refactoring)

import numpy
from numpy import array, transpose,outer, ones, zeros, copy, divide, diag, dot
from numpy.linalg import norm, solve

class DLLMDirect:
    """
    Direct solver for the lifting line wing model
    """
    STOP_CRITERIA_RESIDUAL='Residual decrease'
    STOP_CRITERIA_N_ITER='Niterations'
    
    def __init__(self, LLW):
        self.__LLW = LLW
        self.__K   = self.__LLW.get_K()
        self.__N   = self.get_wing_geom().get_n_elements()
        
        # initialize numerical parameters
        self.__init_numerical_parameters()
        
        # initialize local variables
        self.__init_local_variables()
        
    #-- Accessors
    def get_wing_geom(self):
        return self.__LLW.get_wing_geom()
    
    def get_iAoA(self):
        return self.__iAoA
    
    def get_OC(self):
        return self.__LLW.get_OC()
    
    #-- Numerical parameters related methods
    def __init_numerical_parameters(self):
        self.__relaxFactor=0.99
        self.__stopCriteria=1e-6
        self.__stop_criteria_type=self.STOP_CRITERIA_RESIDUAL
        self.set_stop_criteria()
        
    def set_relax_factor(self, relax_factor):
        self.__relaxFactor = relax_factor
        
    def set_stop_criteria(self,residual=None,n_it=None):
        if n_it is not None:
            if type(n_it) != type(1):
                raise Exception, "n_it stop criteria must be an integer"
            self.__stop_criteria_type=self.STOP_CRITERIA_N_ITER
            self.__stopCriteria=n_it
        else :
            self.__stop_criteria_type=self.STOP_CRITERIA_RESIDUAL
            if residual is not None:
                if type(residual) != type(1.1):
                    raise Exception, "residual stop criteria must be a float"
                self.__stopCriteria=residual
            else:
                self.__stopCriteria=1e-6
    
    #-- Computation related methods
    def __init_local_variables(self):
        #Initializing local variables for lifting line computations
        self.__R0       = None
        self.__residual = None
        self.__residuals_hist=[]
        self.__localAoA=zeros([self.__N])
        self.__DlocalAoA_DIAoA=-diag(ones([self.__N]))
        self.__DlocalAoA_DTwist=-diag(ones([self.__N]))
        
        self.__DlocalAoA_DAoA=ones(self.__N)
        self.__gamma=zeros(self.__N)
        self.__Dgamma_DlocalAoA=zeros([self.__N,self.__N])
        self.__Dgamma_Dthickness=zeros([self.__N,self.__N])
        self.__dGamma=zeros([self.__N+1])
        self.__DdGammaDy_DGamma=zeros([self.__N+1,self.__N])
        self.__iAoA=zeros([self.__N])
        self.__iAoANew=zeros([self.__N])
        self.__DiAoA_DdGamma=zeros([self.__N,self.__N+1])
        self.__DR_DiAoA=zeros([self.__N])
        self.__R=zeros([self.__N])
        
    def run(self):
        self.__init_iterations()
        
        if self.__stop_criteria_type == self.STOP_CRITERIA_RESIDUAL:
            while(self.__residual>self.__stopCriteria):
                self.__sub_iteration()
        else:
            for i in xrange(self.__stopCriteria):
                self.__sub_iteration() 
        
        
    

        
    