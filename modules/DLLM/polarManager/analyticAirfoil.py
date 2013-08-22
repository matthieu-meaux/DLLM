# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0.1
# @author: François Gallard

# - Local imports -
from airfoil import Airfoil
from math import pi

class AnalyticAirfoil(Airfoil):
    """
    An analytic airfoil based on linear theory
    """
    THICKNESS_CORRECTION=0.7698
    
    def __init__(self,Sref, Lref, alpha0,Cd0,ClAlpha=2.*pi,Cm0=0.0,relative_thickness=0.0):
        '''
        Constructor for airfoils
        @param Sref : reference surface
        @param Lref : reference length
        @param alpha0:angle of attack of null lift
        @param Cd0 : drag
        @param ClAlpha: lift gradient
        @param Cm0 : pitch moment
        @param relative_thickness : relative thickness for addition of a correction of ClAlpha
        '''
        if relative_thickness<0. or relative_thickness>0.25:
            raise Exception, "Relative thickness must be >0. and <0.25"
        
        Airfoil.__init__(self,Sref,Lref,relative_thickness)
        self.__alpha0=alpha0
        self.__Cd0=Cd0
        self.__ClAlpha=ClAlpha
        self.__Cm0=Cm0
        
    def Cl(self,alpha,Mach=0.0):
        '''
        Lift coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        return self.ClAlpha(alpha,Mach)*(alpha-self.__alpha0)
    
    def dCl_dthickness(self,alpha,Mach=0.0):
        """
        Derivative of Cl with respect to relative thickness
        """
        return self.__d_ClAlpha_d_thickness(alpha,Mach)*(alpha-self.__alpha0)
    
    def __d_ClAlpha_d_thickness(self,alpha,Mach=0.0):
        """
        Derivative of ClAlpha with respect to relative thickness
        """        
        prandlt_corr=self.Prandtl_correction(Mach)
        dthick_corr=self.THICKNESS_CORRECTION
        
        return self.__ClAlpha*dthick_corr*prandlt_corr

    def Cd(self,alpha,Mach=0.0):
        '''
        Drag coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        return self.__Cd0
    
    def Cm(self,alpha,Mach=0.0):
        '''
        Pitch moment coefficient
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        return self.__Cm0
    
    def ClAlpha(self,alpha,Mach=0.0):
        '''
        Sensibility of Cl to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        
        thick_corr=1.+self.THICKNESS_CORRECTION*self.get_relative_thickness()
        prandlt_corr=self.Prandtl_correction(Mach)
        
        return self.__ClAlpha*thick_corr*prandlt_corr
    
    def CdAlpha(self,alpha,Mach=0.0):
        '''
        Sensibility of Cd to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        return 0.0
    
    def CmAlpha(self,alpha,Mach=0.0):
        '''
        Sensibility of Cm to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        return 0.0
    
    def get_scaled_copy(self,Sref,Lref,relative_thickness=0.0):
        """
        Instances a copy of this airfoil with different Lref and Sref
        @param Sref : reference surface
        @param Lref : reference length
        """
        return AnalyticAirfoil(Sref,Lref,self.__alpha0,self.__Cd0,self.__ClAlpha,self.__Cm0,relative_thickness=relative_thickness)
