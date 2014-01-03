# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: Francois Gallard
# @author: Matthieu MEAUX (for refactoring)

import numpy
from numpy import array, transpose,outer, ones, zeros, copy, divide, diag, dot
from numpy.linalg import norm, solve

from MDOTools.Solvers.newton_raphson_problem import NewtonRaphsonProblem

class DLLMDirect:
    """
    Direct solver for the lifting line wing model
    """
    def __init__(self, LLW):
        self.__LLW          = LLW
        self.__K            = self.__LLW.get_K()
        self.__dK_dchi      = self.__LLW.get_dK_dchi()
        self.__N            = self.get_wing_param().get_n_sect()
        self.__ndv          = self.get_wing_param().get_ndv()
        self.__computed     = False
        self.__gamma_f_name = None
        
        # initialize local variables
        self.__init_local_variables()
        
        # initialize the Newton-Raphson problem
        self.__NRPb = None 
        self.__init_Newton_Raphson()
        
    #-- Accessors
    def get_wing_param(self):
        return self.__LLW.get_wing_param()
    
    def get_airfoils(self):
        return self.__LLW.get_airfoils()
    
    def get_OC(self):
        return self.__LLW.get_OC()
    
    def get_iAoA(self):
        return self.__iAoA
    
    def get_localAoA(self):
        return self.__localAoA
    
    def get_R(self):
        return self.__R
    
    def get_dpR_dpiAoA(self):
        return self.__dpR_dpiAoA
    
    def get_dpR_dpchi(self):
        return self.__dpR_dpchi
    
    def get_dplocalAoA_dpiAoA(self):
        return self.__dplocalAoA_dpiAoA
    
    def get_dplocalAoA_dpchi(self):
        return self.__dplocalAoA_dpchi
    
    def get_convergence_history(self):
        """
        Accessor to the last computation convergence history as a list of residuals normalized by the first iteration residual.
        """
        return self.__residuals_hist
    
    def is_computed(self):
        return self.__computed
    
    #-- Setters
    def set_computed(self, bool=True):
        self.__computed = bool
        
    def set_gamma_file_name(self, gamma_f_name):
        self.__gamma_f_name = gamma_f_name
    
    #-- Newton-Raphson related methods
    def __init_Newton_Raphson(self):
        iAoA0= zeros(self.__N)
        self.__NRPb = NewtonRaphsonProblem(iAoA0, self.comp_R, self.comp_dpR_dpiAoA)
        self.__NRPb.set_relax_factor(0.99)
        self.__NRPb.set_stop_residual(1.e-9)
        self.__NRPb.set_max_iterations(100)
        
    def set_relax_factor(self, relax_factor):
        self.__NRPb.set_relax_factor(relax_factor)
    
    def set_stop_residual(self, residual):
        self.__NRPb.set_stop_residual(residual)
        
    def set_max_iterations(self, max_it):
        self.__NRPb.set_max_iterations(max_it)

    #-- Computation related methods
    def run(self):
        self.__NRPb.solve()
        self.set_computed(True)
        self.__write_gamma_to_file()
        self.comp_dpR_dpchi()
                
    def __init_local_variables(self):
        # Initializing local variables for lifting line computations
        # Residual variables
        self.__R            = zeros([self.__N])
        self.__dpR_dpiAoA   = None
        self.__dpR_dpchi    = None
        self.__dpR_dpthetaY = None
        
        # Angle of attack variables
        self.__localAoA            = zeros([self.__N])
        self.__dplocalAoA_dpiAoA   = -diag(ones([self.__N])) # This is a constant matrix
        self.__dplocalAoA_dpthetaY = diag(ones([self.__N]))  # This is a constant matrix
        self.__dplocalAoA_dpchi    = zeros([self.__N,self.__ndv])
        
        # Induced angle of attack variables
        self.__iAoA    = None
        self.__iAoANew = None
        self.__dpiAoAnew_dpchi  = zeros([self.__N,self.__ndv])   
        self.__dpiAoAnew_dpiAoA = None
        
        # Circulation variables
        self.__gamma  = zeros(self.__N)
        self.__dpgamma_dpiAoA   = None
        self.__dpgamma_dpthetaY = None
        self.__dpgamma_dplocalAoA = zeros([self.__N,self.__N])
        self.__dpgamma_dpchi      = zeros([self.__N,self.__ndv])
                
    #-- Residual related methods
    def comp_R(self, iAoA):
        self.__iAoA = iAoA
        self.__compute_localAoA()
        self.__compute_gamma()
        self.__compute_iAoAnew()
        
        self.__R = self.__iAoA - self.__iAoANew
        
        return self.__R
    
    def comp_dpR_dpiAoA(self, iAoA):
        R=self.comp_R(iAoA)
        
        # dplocalAoA_dpiAoA is a constant matrix, no need to compute it (self.__dplocalAoA_dpiAoA)
        self.__compute_dpgamma_dpiAoA()
        self.__compute_dpiAoAnew_dpiAoA()
        
        self.__dpR_dpiAoA=numpy.diag(ones([self.__N]))-self.__dpiAoAnew_dpiAoA
        
        return self.__dpR_dpiAoA
    
    def comp_dpR_dpchi(self):
        self.__compute_dplocalAoA_dpchi()
        self.__compute_dpgamma_dpchi()
        self.__compute_dpiAoAnew_dpchi()
        self.__dpR_dpchi = -dot(self.__K,self.__dpgamma_dpchi) - self.__dpiAoAnew_dpchi
        
        return self.__dpR_dpchi
    
    def comp_dpR_dpthetaY(self):
        self.__compute_dpgamma_dpthetaY()    
        self.__dpR_dpthetaY = -dot(self.__K,self.__dpgamma_dpthetaY)
        
        return self.__dpR_dpthetaY
    
    def __compute_dpgamma_dpthetaY(self):
         Mach = self.get_OC().get_Mach()
         for i in xrange(self.__N):
             self.__dpgamma_dplocalAoA[i,i]  = self.get_airfoils()[i].dpgamma_dpAoA(self.__localAoA[i],Mach)
             
         self.__dpgamma_dpthetaY = dot(self.__dpgamma_dplocalAoA,self.__dplocalAoA_dpthetaY)
    
    def __compute_dplocalAoA_dpchi(self):
        twist_grad  = self.get_wing_param().get_twist_grad()
        AoA_grad    = self.get_wing_param().get_AoA_grad()
        if AoA_grad is None:
            for i in xrange(self.__N):
                self.__dplocalAoA_dpchi[i,:] =twist_grad[i,:] # + dAoAdksi ?? 
        else:
            for i in xrange(self.__N):
                self.__dplocalAoA_dpchi[i,:] = AoA_grad[:] + twist_grad[i,:] 
            
    def __compute_dpgamma_dpchi(self):
        Mach = self.get_OC().get_Mach()
        for i in xrange(self.__N):
            self.__dpgamma_dpchi[i,:] = self.get_airfoils()[i].dpgamma_dpchi(self.__localAoA[i],Mach)
        
        self.__dpgamma_dpchi = self.__dpgamma_dpchi + dot(self.__dpgamma_dplocalAoA,self.__dplocalAoA_dpchi) 
            
    def __compute_dpiAoAnew_dpchi(self):
        for n in xrange(self.__ndv):
            self.__dpiAoAnew_dpchi[:,n] = dot(self.__dK_dchi[:,:,n],self.__gamma)

    def __compute_localAoA(self):
        Thetay = self.get_wing_param().get_thetaY()
        twist  = self.get_wing_param().get_twist()
        AoA    = self.get_wing_param().get_AoA()
        if AoA is None:
           AoA    = self.get_OC().get_AoA_rad()
        else:
           self.get_OC().set_AoA(AoA*180./numpy.pi)
         
         
        # Why this formula ? twist increases the local airfoil angle of attack normally...
        #self.__localAoA=alpha-iaOa-self.get_wing_geom().get_twist()
        self.__localAoA = AoA + twist - self.__iAoA + Thetay
        
        for i in xrange(self.__N):
            if self.__localAoA[i] > numpy.pi/2. or self.__localAoA[i] < -numpy.pi/2.:
                raise Exception, "Local angle of attack out of bounds [-pi/2, pi/2]"

    def __compute_gamma(self):
        """
        Update the circulation
        """
        Mach = self.get_OC().get_Mach()
        for i in xrange(self.__N):
            self.__gamma[i] = self.get_airfoils()[i].gamma(self.__localAoA[i],Mach)
            
    def __compute_dpgamma_dpiAoA(self):
         Mach = self.get_OC().get_Mach()
         for i in xrange(self.__N):
             self.__dpgamma_dplocalAoA[i,i]  = self.get_airfoils()[i].dpgamma_dpAoA(self.__localAoA[i],Mach)
             
         self.__dpgamma_dpiAoA = dot(self.__dpgamma_dplocalAoA,self.__dplocalAoA_dpiAoA)

    def __compute_iAoAnew(self):
        '''
        Computes the induced angle on an airfoil for a given circulation on the wing.
        '''
        self.__iAoANew        = dot(self.__K,self.__gamma)
        
    def __compute_dpiAoAnew_dpiAoA(self):
        """
        Computes the derivative dpiAoAnew_dpiAoA
        """
        self.__dpiAoAnew_dpiAoA = dot(self.__K,self.__dpgamma_dpiAoA)
        
    def __write_gamma_to_file(self):
        '''
        Writes the circulation repartition in a file
        '''
        if self.__gamma_f_name is None:
            gamma_f_name='gamma.dat'
        else:
            gamma_f_name=self.__gamma_f_name
            
        fid=open(gamma_f_name,'w')
        line="#Slice\t%24s"%"Circulation"+"\n"
        fid.write(line)
        i=0
        for i in range(len(self.__gamma)):
            line=str(i)+"\t%24.16e"%self.__gamma[i]+"\n"
            fid.write(line)
        fid.close()
