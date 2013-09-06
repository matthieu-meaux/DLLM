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
    
    def __init__(self, Sref, Lref, relative_thickness=0.0):
        '''
        Constructor
        @param Sref : reference surface
        @param Lref : reference length
        '''
        DifferentiatedAeroShape.__init__(self,Sref,Lref)
        self.__relative_thickness=relative_thickness
        
    def set_relative_thickness(self,relative_thickness):
        """
        Setter for the relative_thickness
        """
        self.__relative_thickness=relative_thickness
        
    def get_relative_thickness(self):
        """
        Accessor for the relative_thickness
        """
        return self.__relative_thickness
    
    def gamma(self,alpha,Mach=0.0):
        '''
        Computes the circulation due to lift at aOa alpha
        @param alpha: angle of Attack 
        '''
        return 0.5*self.get_Lref()*self.Cl(alpha,Mach=Mach)
    
    def dCl_dthickness(self,alpha,Mach=0.0):
        """
        Derivative of Cl with respect to relative thickness
        """
        return 0.0
    
    def dCd_dthickness(self,alpha,Mach=0.0):
        """
        Derivative of Cd with respect to relative thickness
        """
        return 0.0    
 
    def dCm_dthickness(self,alpha,Mach=0.0):
        """
        Derivative of Cm with respect to relative thickness
        """
        return 0.0  
     
    def DGammaDAoA(self,alpha,Mach):
        '''
        Computes the sensibility of the circulation to the angle of attack
        @param alpha: angle of Attack 
        '''
        return 0.5*self.get_Lref()*self.ClAlpha(alpha,Mach)

    def DGammaDThickness(self,alpha,Mach):
        """
        Sensibility of Gamma to the relative thickness of the airfoil.
        """
        return 0.5*self.get_Lref()*self.dCl_dthickness(alpha,Mach)
    
    def get_scaled_copy(self,Sref,Lref,relative_thickness=0.0):
        """
        Instances a copy of this airfoil with different Lref and Sref
        """
        return Airfoil(Sref,Lref,relative_thickness=relative_thickness)
    

