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
    
    def __init__(self, AoA0=0., Cd0=0., Cm0=0.0, rel_thick=0.0, Sref=1., Lref=1.):
        '''
        Constructor for airfoils
        @param AoA0:angle of attack of null lift
        @param Cd0 : drag
        @param Cm0 : pitch moment
        @param relative_thickness : relative thickness for addition of a correction of ClAlpha
        @param Sref : reference surface
        @param Lref : reference length
        '''
        if rel_thick<0. or rel_thick>0.25:
            raise Exception, "Relative thickness must be >0. and <0.25"
        
        Airfoil.__init__(self,Sref,Lref,rel_thick)
        self.__AoA0=AoA0
        self.__Cd0=Cd0
        self.__Cm0=Cm0
        
        self.__ClAlpha_base=2.*numpy.pi
        
    def Cl(self,alpha,Mach=0.0):
        return self.ClAlpha(alpha,Mach)*(alpha-self.__AoA0)
    
    def dCl_dchi(self,alpha,Mach=0.0):
        """
        Derivative of Cl with respect to relative thickness
        """
        return self.__d_ClAlpha_d_chi(alpha,Mach)*(alpha-self.__AoA0)
    
    def __d_ClAlpha_d_chi(self, alpha,Mach=0.0):
        prandlt_corr=self.Prandtl_correction(Mach)
        dthick_corr=self.THICKNESS_CORRECTION*self.get_rel_thick_grad()
        dClalpha_dchi=self.__ClAlpha_base*dthick_corr*prandlt_corr
        return dClalpha_dchi
    
    def Cd(self,alpha,Mach=0.0):
        return self.__Cd0
    
    def Cm(self,alpha,Mach=0.0):
        return self.__Cm0
    
    def ClAlpha(self,alpha,Mach):
        thick_corr=1.+self.THICKNESS_CORRECTION*self.get_rel_thick()
        prandlt_corr=self.Prandtl_correction(Mach)
        
        return self.__ClAlpha_base*thick_corr*prandlt_corr
    
    def CdAlpha(self,alpha,Mach):
        return 0.0
    
    def CmAlpha(self,alpha,Mach):
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
    
    def get_scaled_copy(self, Sref=None,Lref=None, rel_thick=None):
        """
        Instances a copy of this airfoil with different Lref and Sref
        @param Sref : reference surface
        @param Lref : reference length
        """
        if rel_thick is None:
            rel_thick = self.get_rel_thick()
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        return AnalyticAirfoil(self.__AoA0,self.__Cd0, self.__Cm0,rel_thick=rel_thick, Sref=Sref, Lref=Lref)
