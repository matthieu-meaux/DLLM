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
    
    def get_airfoils(self):
        return self.__LLW.get_airfoils()
    
    def get_iAoA(self):
        return self.__iAoA
    
    def get_OC(self):
        return self.__LLW.get_OC()
    
    def get_convergence_history(self):
        """
        Accessor to the last computation convergence history as a list of residuals normalized by the first iteration residual.
        """
        return self.__residuals_hist
    
    def get_DLocalAOa_DTwist(self):
        return self.__DlocalAoA_DTwist
 
    def get_DLocalAOa_DAoA(self):
        return self.__DlocalAoA_DAoA
    
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
    def run(self):
        print 'Start direct solver newton iteration...'
        self.__init_iterations()
        
        if self.__stop_criteria_type == self.STOP_CRITERIA_RESIDUAL:
            while(self.__residual>self.__stopCriteria):
                self.__iteration()
        else:
            for i in xrange(self.__stopCriteria):
                self.__iteration()
        print 'Done.'
                
    def __init_local_variables(self):
        # Initializing local variables for lifting line computations
        # Residual variables
        self.__R0       = None
        self.__R        = zeros([self.__N])
        self.__DR_DiAoA = zeros([self.__N])
        self.__residual = None
        self.__residuals_hist=[]
        
        
        # Angle of attack variables
        self.__localAoA = zeros([self.__N])
        self.__DlocalAoA_DiAoA  = -diag(ones([self.__N])) # This is a constant matrix
        self.__DlocalAoA_DTwist = diag(ones([self.__N]))  # This is a constant matrix
        self.__DlocalAoA_DAoA   = ones(self.__N)          # This is a constant matrix
        
        # Induced angle of attack variables
        self.__iAoA    = zeros([self.__N])
        self.__iAoANew = zeros([self.__N])
        self.__DiAoA_DdGamma=zeros([self.__N,self.__N+1])
        self.__DiAoA_DiAoA = None
        
        # Circulation variables
        self.__gamma  = zeros(self.__N)
        self.__Dgamma_DiAoA = None
        self.__DdGammaDy_DiAoA = None
        self.__Dgamma_DlocalAoA=zeros([self.__N,self.__N])
        self.__Dgamma_Dthickness=zeros([self.__N,self.__N])
        self.__dGamma=zeros([self.__N+1])
        self.__DdGammaDy_DGamma=zeros([self.__N+1,self.__N])
        
    def __init_iterations(self):
        self.__iAoA=zeros([self.__N])
        self.__Residual()
        self.__R0       = norm(self.__R)
        self.__residual = 1.0
        self.__residuals_hist = []
        
    def __iteration(self):
        self.__Residual()
        self.__newton_iteration()
        
    def __newton_iteration(self):
        #Newton zero search method
        self.__iAoA-=self.__relaxFactor*solve(self.__DR_DiAoA,self.__R)
        
        #Compute stop criteria
        self.__residual=norm(self.__R)/self.__R0
        self.__residuals_hist.append(self.__residual)
        print "||R||/||R0||= "+str(self.__residual)
                
    #-- Residual related methods
    def __Residual(self):
        '''
        Computes the residual and its derivatives
        '''
        self.__compute_localAoA()
        self.__compute_gamma()
        self.__compute_dGamma()
        self.__compute_iAoA()
        
        self.__R=self.__iAoA-self.__iAoANew
        self.__DR_DiAoA=numpy.diag(ones([self.__N]))-self.__DiAoA_DiAoA
    
    def __compute_localAoA(self):
        '''
        Computes the local angle of attack = AoA + twist - induced downwash angle
        '''
        AoA = self.get_OC().get_AoA_rad()
        
        # Why this formula ? twist increases the local airfoil angle of attack normally...
        #self.__localAoA=alpha-iaOa-self.get_wing_geom().get_twist()
        self.__localAoA = AoA + self.get_wing_geom().get_twist() - self.__iAoA 
        
        for i in range(self.__N):
            if self.__localAoA[i] > numpy.pi/2. or self.__localAoA[i] < -numpy.pi/2.:
                raise Exception, "Local angle of attack out of bounds [-pi/2, pi/2]"
            
        #print 'local AoA = ',self.__localAoA
        
    def __compute_gamma(self):
        '''
        Update the circulation
        '''
        Mach = self.get_OC().get_Mach()
        for i in range(self.__N):
            self.__gamma[i] = self.get_airfoils()[i].gamma(self.__localAoA[i],Mach)
            self.__Dgamma_DlocalAoA[i,i]  = self.get_airfoils()[i].DGammaDAoA(self.__localAoA[i],Mach)
            self.__Dgamma_Dthickness[i,i] = self.get_airfoils()[i].DGammaDThickness(self.__localAoA[i],Mach)
        
        self.__Dgamma_DiAoA = dot(self.__Dgamma_DlocalAoA,self.__DlocalAoA_DiAoA)
        
    def __compute_dGamma(self):
        '''
        Computes the circulation y derivation
        '''
        self.__dGamma[0:self.__N]   = self.__gamma
        self.__dGamma[self.__N]     = 0.0
        self.__dGamma[1:self.__N+1]-= self.__gamma
        
        #Differentiation
        self.__DdGammaDy_DGamma[0:self.__N,:]   = diag(ones([self.__N]))
        self.__DdGammaDy_DGamma[self.__N,:]     = 0.0
        self.__DdGammaDy_DGamma[1:self.__N+1,:]-= diag(ones([self.__N]))
        
        self.__DdGammaDy_DiAoA = dot(self.__DdGammaDy_DGamma,self.__Dgamma_DiAoA)
        
    def __compute_iAoA(self):
        '''
        Computes the induced angle on an airfoil for a given circulation on the wing.
        '''
        self.__iAoANew       = dot(self.__K,self.__dGamma)
        self.__DiAoA_DdGamma = dot(self.__K,diag(ones(self.__N+1)))
        self.__DiAoA_DiAoA   = dot(self.__DiAoA_DdGamma,self.__DdGammaDy_DiAoA)
            