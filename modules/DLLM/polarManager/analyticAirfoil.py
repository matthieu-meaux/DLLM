# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
#  DLLM (non-linear Differentiated Lifting Line Model, open source software)
# 
#  Copyright (C) 2013-2015 Airbus Group SAS
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
#  https://github.com/matthieu-meaux/DLLM.git
#
# @author : Francois Gallard
# @author : Matthieu Meaux
#
# - Local imports -
import numpy as np
from airfoil import Airfoil

class AnalyticAirfoil(Airfoil):
    """
    An analytic airfoil based on linear theory
    """
    THICKNESS_CORRECTION = 0.7698
    
    def __init__(self, OC, AoA0=0., Sref=1., Lref=1., rel_thick=0.0, sweep=0.0, pcop=0.25, Ka=0.95):
        '''
        Constructor for airfoils
        @param AoA0:angle of attack of null lift
        @param Cd0 : drag
        @param Cm0 : pitch moment
        @param relative_thickness : relative thickness for addition of a correction of ClAlpha
        @param Sref : reference surface
        @param Lref : reference length
        '''
        if rel_thick<0. or rel_thick>0.30:
            raise Exception, "Relative thickness must be >0. and <0.30"
        
        Airfoil.__init__(self, OC, Sref, Lref)
        
        self.__AoA0_deg  = None
        self.__AoA0      = None
        self.set_AoA0(AoA0)
        
        self.__pcop      = None
        self.set_pcop(pcop)
        
        self.__Ka        = None
        self.set_Ka(Ka)

        self.__rel_thick      = None
        self.__rel_thick_grad = None
        self.set_rel_thick(rel_thick)
        
        self.__sweep          = None
        self.__sweep_grad     = None
        self.set_sweep(sweep)
        
    #-- Setters
    def set_AoA0(self, AoA0_deg):
        self.__AoA0_deg = AoA0_deg
        self.__AoA0     = AoA0_deg*np.pi/180.
        
    def set_pcop(self, pcop):
        self.__pcop = pcop
        
    def set_Ka(self, Ka):
        self.__Ka = Ka
        
    def set_rel_thick(self, rel_thick):
        self.__rel_thick = rel_thick
        
    def set_rel_thick_grad(self, rel_thick_grad):
        self.__rel_thick_grad = rel_thick_grad
        
    def set_sweep(self, sweep):
        self.__sweep = sweep
        
    def set_sweep_grad(self, sweep_grad):
        self.__sweep_grad = sweep_grad
    
    #-- Accessors
    def get_rel_thick(self):
        return self.__rel_thick
    
    def get_rel_thick_grad(self):
        return self.__rel_thick_grad
    
    def get_sweep(self):
        return self.__sweep
    
    def get_sweep_grad(self):
        return self.__sweep_grad
    
    #-- Cl related methods
    def Cl(self, AoA, Mach):
        dCl_dAoA = self.dCl_dAoA(AoA, Mach)
        Cl = dCl_dAoA*(AoA-self.__AoA0)
        return Cl
    
    def dCl_dAoA(self, AoA, Mach):
        sweep = self.get_sweep()
        toc   = self.get_rel_thick()
        
        Mach_normal = Mach*np.cos(sweep)
        # Base slope with Prandtl correction
        if Mach_normal<0.9:
            base_slope=2.*np.pi/np.sqrt(1.-Mach_normal**2)
        elif 0.9<=Mach_normal<=1.1:
            s_sub = 2.*np.pi/np.sqrt(abs(1.-0.9**2))
            s_sup = 4./np.sqrt(abs(1.1**2-1.))
            fact = (1.1-Mach_normal)/0.2
            base_slope = s_sub*fact+s_sup*(1.-fact)
        else:
            base_slope = 4./np.sqrt(Mach_normal**2-1.)
            
        # thickness correction
        thick_corr = 1. + self.THICKNESS_CORRECTION*toc/np.cos(sweep)
        
        # sweep correction
        sweep_corr = np.cos(sweep)**2
        
        dCl_dAoA = base_slope*thick_corr*sweep_corr
        
        return dCl_dAoA
    
    def __d2Cl_dAoAdchi(self, AoA, Mach):
        sweep  = self.get_sweep()
        dsweep = self.get_sweep_grad()
        
        toc    = self.get_rel_thick()
        dtoc   = self.get_rel_thick_grad()
        
        Mach_normal  = Mach*np.cos(sweep)
        dMach_normal = -Mach*np.sin(sweep)*dsweep
        
        # Base slope with Prandtl correction
        if Mach_normal<0.9:
            base_slope = 2.*np.pi/np.sqrt(1.-Mach_normal**2)
            dbase_slope = 2.*np.pi*(Mach_normal*dMach_normal)/((1.-Mach_normal**2)*np.sqrt(1.-Mach_normal**2))
        elif 0.9<=Mach_normal<=1.1:
            s_sub = 2.*np.pi/np.sqrt(abs(1.-0.9**2))
            s_sup = 4./np.sqrt(abs(1.1**2-1.))
            fact  = (1.1-Mach_normal)/0.2
            dfact = -dMach_normal/0.2
            base_slope  = s_sub*fact+s_sup*(1.-fact)
            dbase_slope = s_sub*dfact+s_sup*(1.-dfact)
        else:
            base_slope  = 4./np.sqrt(Mach_normal**2-1.)
            dbase_slope = -4.*(Mach_normal*dMach_normal)/((1.-Mach_normal**2)*np.sqrt(1.-Mach_normal**2))
            
        # thickness correction
        thick_corr  = 1. + self.THICKNESS_CORRECTION*toc/np.cos(sweep)
        dthick_corr = self.THICKNESS_CORRECTION*(dtoc*np.cos(sweep)+toc*np.sin(sweep)*dsweep)/np.cos(sweep)**2
        
        # sweep correction
        sweep_corr  = np.cos(sweep)**2
        dsweep_corr = -2.*np.sin(sweep)*np.cos(sweep)*dsweep
        
        
        d2Cl_dAoAdchi = dbase_slope*thick_corr*sweep_corr\
                      + base_slope*dthick_corr*sweep_corr\
                      + base_slope*thick_corr*dsweep_corr
                      
        return d2Cl_dAoAdchi
                      
    def dCl_dchi(self, AoA, Mach):
        d2Cl_dAoAdchi = self.__d2Cl_dAoAdchi(AoA, Mach)
        dCl_dchi      =  d2Cl_dAoAdchi*(AoA-self.__AoA0)
        return dCl_dchi
    
    #-- Cdw related methods
    def Cdw(self, AoA, Mach):
        sweep  = self.get_sweep()
        toc    = self.get_rel_thick()
        Cl     = self.Cl(AoA, Mach)
        
        Mdd    = self.__Ka/np.cos(sweep) - toc/(np.cos(sweep))**2 - Cl/(10.*np.cos(sweep)**3)
        Mcrit  = Mdd - (0.1/80)**(1./3.)
        if Mach < Mcrit:
            Cdw = 0.0
        else:
            Cdw = 20.*(Mach-Mcrit)**4
        return Cdw
    
    def dCdw_dAoA(self, AoA, Mach):
        sweep  = self.get_sweep()
        toc    = self.get_rel_thick()
        Cl     = self.Cl(AoA, Mach)
        dCl    = self.dCl_dAoA(AoA, Mach)
        
        Mdd    = self.__Ka/np.cos(sweep) - toc/np.cos(sweep)**2 - Cl/(10.*np.cos(sweep)**3)
        dMdd   = -dCl/(10.*np.cos(sweep)**3)
        Mcrit  = Mdd - (0.1/80)**(1./3.)
        dMcrit = dMdd
        
        if Mach < Mcrit:
            dCdw_dAoA = 0.0
        else:
            dCdw_dAoA = -80.*(Mach-Mcrit)**3*dMcrit
        return dCdw_dAoA
    
    def dCdw_dchi(self, AoA, Mach):
        sweep  = self.get_sweep()
        dsweep = self.get_sweep_grad()
        
        toc    = self.get_rel_thick()
        dtoc   = self.get_rel_thick_grad()
        
        Cl     = self.Cl(AoA, Mach)
        dCl    = self.dCl_dchi(AoA, Mach)
        
        Mdd    = self.__Ka/np.cos(sweep) - toc/np.cos(sweep)**2 - Cl/(10.*np.cos(sweep)**3)
        
        dMdd   = self.__Ka*np.sin(sweep)*dsweep/(np.cos(sweep))**2 \
               - dtoc/(np.cos(sweep))**2 - toc*(2.*np.sin(sweep)*dsweep)/np.cos(sweep)**3 \
               - dCl/(10.*np.cos(sweep)**3) - Cl*(3.*np.sin(sweep)*dsweep)/(10.*np.cos(sweep)**4)
               
        Mcrit  = Mdd - (0.1/80)**(1./3.)
        dMcrit = dMdd
        
        if Mach < Mcrit:
            dCdw_dchi = 0.0
        else:
            dCdw_dchi = -80.*(Mach-Mcrit)**3*dMcrit
            
        return dCdw_dchi
    
    #-- Cdvp related methods
    def Cdvp(self, AoA, Mach):
        sweep = self.get_sweep()
        Cl  = self.Cl(AoA, Mach)
        Cdf = self.Cdf(AoA, Mach)
        toc = self.get_rel_thick()/np.cos(sweep)
        
        Cdvp_min = 60.*toc**4*Cdf
        
        Cdvp = Cdvp_min + 2.*(Cdvp_min + Cdf)*Cl**2
        return Cdvp
    
    def dCdvp_dAoA(self, AoA, Mach):
        sweep = self.get_sweep()
        Cl  = self.Cl(AoA, Mach)
        Cdf = self.Cdf(AoA, Mach)
        toc = self.get_rel_thick()/np.cos(sweep)
        
        dCl = self.dCl_dAoA(AoA, Mach)
        
        Cdvp_min = 60.*toc**4*Cdf
        
        dCdvp_dAoA = 4.*(Cdvp_min + Cdf)*Cl*dCl
        return dCdvp_dAoA
    
    def dCdvp_dchi(self, AoA, Mach):
        sweep  = self.get_sweep()
        dsweep = self.get_sweep_grad()
        
        Cl  = self.Cl(AoA, Mach)
        dCl = self.dCl_dchi(AoA, Mach)
        
        Cdf  = self.Cdf(AoA, Mach)
        dCdf = self.dCdf_dchi(AoA, Mach)
        
        toc = self.get_rel_thick()/np.cos(sweep)
        dtoc   = (self.get_rel_thick_grad()*np.cos(sweep)+self.get_rel_thick()*np.sin(sweep)*dsweep)/np.cos(sweep)**2
        
        Cdvp_min  = 60.*toc**4*Cdf
        dCdvp_min = 60.*(4.*dtoc*toc**3*Cdf+toc**4*dCdf)
        
        dCdvp_dchi = dCdvp_min + 2.*(dCdvp_min+dCdf)*Cl**2+4.*(Cdvp_min+Cdf)*Cl*dCl
        
        return dCdvp_dchi
    
    #-- Cdf related methods
    def Cdf(self, AoA, Mach):
        OC = self.get_OC()
        sweep  = self.get_sweep()
        c  = OC.get_c()
        nu = OC.get_nu()
        L  = self.get_Lref()
        Re = Mach*np.cos(sweep)*c*L/nu
        toc         = self.get_rel_thick()/np.cos(sweep)
        thick_coeff = 1.+2.1*toc # a rough guess since 1.21 for 10% relative thickness
        if Re < 1.e-12:
            # Prevent division by 0. Drag is null at zero Re number anyway
            Cdf = 0.
        elif Re < 1.e5:
            # Laminar flow
            Cdf=1.328/np.sqrt(Re)*thick_coeff
        else:
            # Turbulent flow
            Cdf=0.074*Re**(-0.2)*thick_coeff
            
        return Cdf
    
    def dCdf_dAoA(self, AoA, Mach):
        dCdf_dAoA = 0.
        return dCdf_dAoA
    
    def dCdf_dchi(self, AoA, Mach):
        OC  = self.get_OC()
        
        sweep  = self.get_sweep()
        dsweep = self.get_sweep_grad()
        
        c   = OC.get_c()
        nu  = OC.get_nu()
        
        L   = self.get_Lref()
        dL  = self.get_Lref_grad()
        
        Re = Mach*np.cos(sweep)*c*L/nu
        dRe = Mach*np.cos(sweep)*c*dL/nu-Mach*np.sin(sweep)*c*L/nu*dsweep
        
        toc    = self.get_rel_thick()/np.cos(sweep)
        dtoc   = (self.get_rel_thick_grad()*np.cos(sweep)+self.get_rel_thick()*np.sin(sweep)*dsweep)/np.cos(sweep)**2
        
        thick_coeff = 1.+2.1*toc
        dthick_coeff = 2.1*dtoc
        
        if Re < 1.e-12:
            # Prevent division by 0. Drag is null at zero Re number anyway
            dCdf_dchi = 0.
        elif Re < 1.e5:
            # Laminar flow
            dCdf_dchi=-0.664/(Re**1.5)*dRe*thick_coeff+1.328/np.sqrt(Re)*dthick_coeff
        else:
            # Turbulent flow
            dCdf_dchi=-0.0148*Re**(-1.2)*dRe*thick_coeff+0.074*Re**(-0.2)*dthick_coeff
        return dCdf_dchi

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
        return AnalyticAirfoil(OC, self.__AoA0_deg, Sref=Sref, Lref=Lref, rel_thick=rel_thick, sweep=sweep, Ka=self.__Ka)

