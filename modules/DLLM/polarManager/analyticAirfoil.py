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
    
    def __init__(self, OC, AoA0=0., Cm0=0.0, Sref=1., Lref=1., rel_thick=0.0, sweep=0.0, Ka=0.95):
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
        
        Airfoil.__init__(self,OC, Sref,Lref,rel_thick,sweep)
        self.__AoA0 = AoA0
        self.__Cm0  = Cm0
        self.__Ka   = Ka
        
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
    
    def dCl_dchi(self, alpha, Mach):
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
            dcorr=Mach/(1.-Mach**2)**(1.5)
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
        Cdw = self.__Cd_wave(alpha, Mach)
        Cdf = self.__Cd_friction(alpha, Mach)
        Cd  = Cdw+ Cdf
        return Cd
    
    def CdAlpha(self,alpha,Mach):
        dCdw = self.__dCd_wave_dAoA(alpha, Mach)
        dCdf = self.__dCd_friction_dAoA(alpha, Mach)
        dCd = dCdw+dCdf
        return dCd
    
    def dCd_dchi(self, alpha, Mach):
        dCdw = self.__dCd_wave_dchi(alpha, Mach)
        dCdf = self.__dCd_friction_dchi(alpha, Mach)
        dCd  = dCdw + dCdf
        return dCd
    
    #-- Cdp related methods
    def Cdp(self,alpha,Mach=0.0):
        Cdw = self.__Cd_wave(alpha, Mach)
        Cd  = Cdw
        return Cd
    
    def CdpAlpha(self,alpha,Mach):
        dCdw = self.__dCd_wave_dAoA(alpha, Mach)
        dCd = dCdw
        return dCd
    
    def dCdp_dchi(self, alpha, Mach):
        dCdw = self.__dCd_wave_dchi(alpha, Mach)
        dCd  = dCdw
        return dCd
    
    #-- Cdf related methods
    def Cdf(self,alpha,Mach=0.0):
        Cdf = self.__Cd_friction(alpha, Mach)
        Cd  = Cdf
        return Cd
    
    def CdfAlpha(self,alpha,Mach):
        dCdf = self.__dCd_friction_dAoA(alpha, Mach)
        dCd = dCdf
        return dCd
    
    def dCdf_dchi(self, alpha, Mach):
        dCdf = self.__dCd_friction_dchi(alpha, Mach)
        dCd  = dCdf
        return dCd
    
        
    def __Cd_wave(self, alpha, Mach):
        sweep  = self.get_sweep()
        toc    = self.get_rel_thick()
        Cl     = self.Cl(alpha, Mach)
        
        Mdd    = self.__Ka/cos(sweep) - toc/(cos(sweep))**2 - Cl/(10.*cos(sweep)**3)
        Mcrit  = Mdd - (0.1/80)**(1./3.)
        if Mach < Mcrit:
            Cdw = 0.0
        else:
            Cdw = 20*(Mach-Mcrit)**4
        return Cdw
    
    def __dCd_wave_dAoA(self, alpha, Mach):
        sweep  = self.get_sweep()
        toc    = self.get_rel_thick()
        Cl     = self.Cl(alpha, Mach)
        dCl    = self.ClAlpha(alpha, Mach)
        
        Mdd    = self.__Ka/cos(sweep) - toc/cos(sweep)**2 - Cl/(10.*cos(sweep)**3)
        dMdd   = -dCl/(10.*cos(sweep)**3)
        Mcrit  = Mdd - (0.1/80)**(1./3.)
        dMcrit = dMdd
        
        if Mach < Mcrit:
            dCdw = 0.0
        else:
            dCdw = -80.*(Mach-Mcrit)**3*dMcrit
        return dCdw
    
    def __dCd_wave_dchi(self, alpha, Mach):
        sweep  = self.get_sweep()
        dsweep = self.get_sweep_grad()
        toc    = self.get_rel_thick()
        dtoc   = self.get_rel_thick_grad()
        Cl     = self.Cl(alpha, Mach)
        dCl    = self.dCl_dchi(alpha, Mach)
        
        Mdd    = self.__Ka/cos(sweep) - toc/cos(sweep)**2 - Cl/(10.*cos(sweep)**3)
        
        dMdd   = self.__Ka*sin(sweep)*dsweep/(cos(sweep))**2 \
               - dtoc/(cos(sweep))**2 - toc*(2.*sin(sweep)*dsweep)/cos(sweep)**3 \
               - dCl/(10.*cos(sweep)**3) - Cl*(3.*sin(sweep)*dsweep)/(10.*cos(sweep)**4)
               
        Mcrit  = Mdd - (0.1/80)**(1./3.)
        dMcrit = dMdd
        
        if Mach < Mcrit:
            dCdw = 0.0
        else:
            dCdw = -80.*(Mach-Mcrit)**3*dMcrit
            
        return dCdw
    
    def __Cd_friction(self, alpha, Mach):
        # Re= V_corr*Lref/nu=Mach*cos(sweep)*c*Lref/nu = Mach*cos(sweep)*c*Lref/nu
        OC = self.get_OC()
        sweep  = self.get_sweep()
        c  = OC.get_c()
        nu = OC.get_nu()
        L  = self.get_Lref()
        Re = Mach*cos(sweep)*c*L/nu
        toc         = self.get_rel_thick()
        thick_coeff = 1.+2.1*toc # a rough guess since 1.21 for 10% relative thickness
        if Re < 1.e-12:
            # Prevent division by 0. Drag is null at zero Re number anyway
            Cdf = 0.
        elif Re < 1.e5:
            # Laminar flow
            Cdf=1.328/sqrt(Re)*thick_coeff
        else:
            # Turbulent flow
            Cdf=0.074*Re**(-0.2)*thick_coeff
            
        return Cdf
    
    def __dCd_friction_dAoA(self, alpha, Mach):
        return 0.
    
    def __dCd_friction_dchi(self, alpha, Mach):
        OC  = self.get_OC()
        sweep  = self.get_sweep()
        dsweep = self.get_sweep_grad()
        c   = OC.get_c()
        nu  = OC.get_nu()
        L   = self.get_Lref()
        dL  = self.get_Lref_grad()
        Re = Mach*cos(sweep)*c*L/nu
        dRe = Mach*cos(sweep)*c*dL/nu-Mach*sin(sweep)*c*L/nu*dsweep
        toc    = self.get_rel_thick()
        dtoc   = self.get_rel_thick_grad()
        thick_coeff = 1.+2.1*toc
        dthick_coeff = 2.1*dtoc
        
        if Re < 1.e-6:
            # Prevent division by 0. Drag is null at zero Re number anyway
            dCdf = 0.
        elif Re < 1.e5:
            # Laminar flow
            dCdf=-0.664/(Re**1.5)*dRe*thick_coeff+1.328/sqrt(Re)*dthick_coeff
        else:
            # Turbulent flow
            dCdf=-0.0148*Re**(-1.2)*dRe*thick_coeff+0.074*Re**(-0.2)*dthick_coeff
        return dCdf
    
    #-- Cm related methods
    def Cm(self,alpha,Mach=0.0):
        return self.__Cm0
    
    def CmAlpha(self,alpha,Mach):
        return 0.0
    
    def get_scaled_copy(self, OC=None, Sref=None,Lref=None, rel_thick=None, sweep=None):
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        if rel_thick is None:
            rel_thick = self.get_rel_thick()
        if sweep is None:
            sweep = self.get_sweep()
        if OC is None:
            OC = self.get_OC()
        return AnalyticAirfoil(OC, self.__AoA0, self.__Cm0, Sref=Sref, Lref=Lref, rel_thick=rel_thick, sweep=sweep)
