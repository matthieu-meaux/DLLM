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
    Supports the computation of the circulation and the sensibility of moments to aOa.
    '''
    
    def __init__(self,Sref,Lref,relative_thickness=0.0):
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
        return 0.5*self.getLref()*self.Cl(alpha,Mach=Mach)
    
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
     
    def dGammaDAoA(self,alpha,Mach=0.0):
        '''
        Computes the sensibility of the circulation to the angle of attack
        @param alpha: angle of Attack 
        '''
        return 0.5*self.getLref()*self.ClAlpha(alpha,Mach=Mach)
    
    def get_scaled_copy(self,Sref,Lref,relative_thickness=0.0):
        """
        Instances a copy of this airfoil with different Lref and Sref
        """
        return Airfoil(Sref,Lref,relative_thickness=relative_thickness)
    
    def Prandtl_correction(self,Mach):
        """
        Computes the Prandls correction to manage compressibility effects for 0.3<Mach<0.7.
        @param Mach  : Mach number
        @param type Mach : Float
        @return : the Prandtl Glauert correction
        @return type : Float
        """
        if Mach<0.:
            raise Exception, "Mach number should be positive."
        elif Mach<0.3:
            return 1.
        elif Mach<0.7:
            return 1./sqrt(abs(1.-Mach**2))
        else:
            raise Exception, "Mach numbers higher than 0.7 are not handled by this correction due to shocks."
    
    def dGammaDThickness(self,alpha,Mach=0.0):
        """
        Sensibility of Gamma to the relative thickness of the airfoil.
        """
        return 0.5*self.getLref()*self.dCl_dthickness(alpha,Mach)
