# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0.2
# @author: François Gallard

# - Local imports -
from numpy import pi,cos, sin, sqrt
from airfoil import Airfoil

class AnalyticAirfoil(Airfoil):
    """
    An analytic airfoil based on linear theory
    """
    THICKNESS_CORRECTION=0.7698
    CLALPHA_BASE=2.*pi
    
    def __init__(self, AoA0=0., Cd0=0., Cm0=0.0, Sref=1., Lref=1., rel_thick=0.0, sweep=0.0):
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
        
        Airfoil.__init__(self,Sref,Lref,rel_thick,sweep)
        self.__AoA0=AoA0
        self.__Cd0=Cd0
        self.__Cm0=Cm0
        
    #-- Cl related methods       
    def Cl(self,alpha,Mach=0.0):
        Cl=self.ClAlpha(alpha,Mach)*(alpha-self.__AoA0)
        return Cl
    
    def ClAlpha(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*cos(sweep)
        prandlt_corr=self.__Prandtl_corr(Mach_normal)
        thick_corr=self.__thick_corr()
        sweep_corr=cos(sweep)**2
        
        Clalpha=self.CLALPHA_BASE*prandlt_corr*thick_corr*sweep_corr
        return Clalpha
    
    def dCl_dchi(self,alpha,Mach=0.0):
        return self.__dClAlpha_dchi(alpha,Mach)*(alpha-self.__AoA0)
    
    def __dClAlpha_dchi(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*cos(sweep)
        prandlt_corr=self.__Prandtl_corr(Mach_normal)
        thick_corr=self.__thick_corr()
        sweep_corr=cos(sweep)**2
        
        dsweep=self.get_sweep_grad()
        dprandlt_corr_dMn  = self.__dPrandtl_corr_dMach(Mach_normal)
        dprandlt_corr      = -dprandlt_corr_dMn*Mach*sin(sweep)*dsweep
        dthick_corr        = self.__dthick_corr_dchi()
        dsweep_corr        = -2.*sin(sweep)*cos(sweep)*dsweep

        dClalpha_dchi = self.CLALPHA_BASE*dprandlt_corr*thick_corr*sweep_corr \
                      + self.CLALPHA_BASE*prandlt_corr*dthick_corr*sweep_corr \
                      + self.CLALPHA_BASE*prandlt_corr*thick_corr*dsweep_corr
                      
        return dClalpha_dchi
    
    def __Prandtl_corr(self,Mach):
        if Mach<0.:
            raise Exception, "Mach number should be positive."
        elif Mach<0.3:
            corr=1.
        elif Mach<1.:
            corr=1./sqrt(abs(1.-Mach**2))
        else:
            raise Exception, 'ERROR: Analytic airfoil with Prandtl correction cannot be used for Mach >= 1.'
        
        return corr
        
    def __dPrandtl_corr_dMach(self, Mach):
        if Mach<0.:
            raise Exception, "Mach number should be positive."
        elif Mach<0.3:
            dcorr=0.
        elif Mach<1.:
            dcorr=Mach/(1-Mach**2)**(1.5)
        else:
            raise Exception, 'ERROR: Analytic airfoil with Prandtl correction cannot be used for Mach >= 1.'
        return dcorr
    
    def __thick_corr(self):
        thick_corr = 1. + self.THICKNESS_CORRECTION*self.get_rel_thick()
        return thick_corr
    
    def __dthick_corr_dchi(self):
        dthick_corr = self.THICKNESS_CORRECTION*self.get_rel_thick_grad()
        return dthick_corr
    
    #-- Cd related methods
    def Cd(self,alpha,Mach=0.0):
        return self.__Cd0
    
    def CdAlpha(self,alpha,Mach):
        return 0.0
    
    #-- Cm related methods
    def Cm(self,alpha,Mach=0.0):
        return self.__Cm0
    
    def CmAlpha(self,alpha,Mach):
        return 0.0
    
    def get_scaled_copy(self, Sref=None,Lref=None, rel_thick=None, sweep=None):
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        if rel_thick is None:
            rel_thick = self.get_rel_thick()
        if sweep is None:
            sweep = self.get_sweep()
        return AnalyticAirfoil(self.__AoA0,self.__Cd0, self.__Cm0, Sref=Sref, Lref=Lref, rel_thick=rel_thick, sweep=sweep)
