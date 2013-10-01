# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.1
# @author: François Gallard

# - Local imports -
from differentiatedAeroShape import DifferentiatedAeroShape
from math import sqrt

class Airfoil(DifferentiatedAeroShape):
    '''
    Airfoil class for lifting line computations. 
    Supports the computation of the circulation and the sensibility of moments to AoA.
    '''
    
    def __init__(self, Sref, Lref, rel_thick=0.0):
        '''
        Constructor
        @param Sref : reference surface
        @param Lref : reference length
        '''
        DifferentiatedAeroShape.__init__(self,Sref,Lref)
        self.__rel_thick = rel_thick
        self.__rel_thick_grad = None 
        
    def set_rel_thick(self, rel_thick):
        self.__rel_thick = rel_thick
        
    def set_rel_thick_grad(self, rel_thick_grad):
        self.__rel_thick_grad = rel_thick_grad
        
    def get_rel_thick(self):
        return self.__rel_thick
    
    def get_rel_thick_grad(self):
        return self.__rel_thick_grad
    
    def gamma(self,alpha,Mach=0.0):
        return 0.5*self.get_Lref()*self.Cl(alpha,Mach=Mach)
    
    def dCl_dchi(self,alpha,Mach=0.0):
        return 0.0
    
    def dCd_dchi(self,alpha,Mach=0.0):
        return 0.0    
 
    def dCm_dchi(self,alpha,Mach=0.0):
        return 0.0  
     
    def DGammaDAoA(self,alpha,Mach):
        return 0.5*self.get_Lref()*self.ClAlpha(alpha,Mach)
    
    def DGammaDchi(self, alpha, Mach):
        dgammadchi =  0.5*self.get_Lref_grad()*self.Cl(alpha,Mach=Mach) + 0.5*self.get_Lref()*self.dCl_dchi(alpha,Mach)
        return dgammadchi
    
    def get_scaled_copy(self,Sref,Lref,rel_thick=0.0):
        return Airfoil(Sref,Lref,rel_thick=rel_thick)
    

