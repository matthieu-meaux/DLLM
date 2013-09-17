# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0.2
# @author: François Gallard

# - Local imports -
import numpy
from airfoil import Airfoil

class AnalyticAirfoil(Airfoil):
    """
    An analytic airfoil based on linear theory
    """
    THICKNESS_CORRECTION=0.7698
    
    def __init__(self, AoA0=0., Cd0=0., Cm0=0.0, relative_thickness=0.0, Sref=1., Lref=1.):
        '''
        Constructor for airfoils
        @param AoA0:angle of attack of null lift
        @param Cd0 : drag
        @param Cm0 : pitch moment
        @param relative_thickness : relative thickness for addition of a correction of ClAlpha
        @param Sref : reference surface
        @param Lref : reference length
        '''
        if relative_thickness<0. or relative_thickness>0.25:
            raise Exception, "Relative thickness must be >0. and <0.25"
        
        Airfoil.__init__(self,Sref,Lref,relative_thickness)
        self.__AoA0=AoA0
        self.__Cd0=Cd0
        self.__Cm0=Cm0
        
        self.__ClAlpha_base=2.*numpy.pi
        
    def Cl(self,alpha,Mach=0.0):
        '''
        Lift coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        return self.ClAlpha(alpha,Mach)*(alpha-self.__AoA0)
    
    def dCl_dthickness(self,alpha,Mach=0.0):
        """
        Derivative of Cl with respect to relative thickness
        """
        return self.__d_ClAlpha_d_thickness(alpha,Mach)*(alpha-self.__AoA0)
    
    def __d_ClAlpha_d_thickness(self,alpha,Mach=0.0):
        """
        Derivative of ClAlpha with respect to relative thickness
        """        
        prandlt_corr=self.Prandtl_correction(Mach)
        dthick_corr=self.THICKNESS_CORRECTION
        
        return self.__ClAlpha_base*dthick_corr*prandlt_corr

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
    
    def ClAlpha(self,alpha,Mach):
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
        
        return self.__ClAlpha_base*thick_corr*prandlt_corr
    
    def CdAlpha(self,alpha,Mach):
        '''
        Sensibility of Cd to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        return 0.0
    
    def CmAlpha(self,alpha,Mach):
        '''
        Sensibility of Cm to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        '''
        return 0.0
    
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
            return 1./numpy.sqrt(abs(1.-Mach**2))
        elif Mach<1.:
            return 1./numpy.sqrt(abs(1.-Mach**2))
            print "WARNING: Mach numbers higher than 0.7 are not handled by this correction due to compressibility effects."
        else:
            raise Exception, 'ERROR: Analytic airfoil with Prandtl correction cannot be used for Mach > 1.'
    
    def get_scaled_copy(self,relative_thickness=None, Sref=None,Lref=None):
        """
        Instances a copy of this airfoil with different Lref and Sref
        @param Sref : reference surface
        @param Lref : reference length
        """
        if relative_thickness is None:
            relative_thickness = self.get_relative_thickness()
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        return AnalyticAirfoil(self.__AoA0,self.__Cd0, self.__Cm0,relative_thickness=relative_thickness, Sref=Sref, Lref=Lref)
