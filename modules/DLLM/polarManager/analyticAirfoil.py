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
    
    def __init__(self, OC, AoA0=0., Sref=1., Lref=1., rel_thick=0.0, pcop=0.25, Ka=0.95, grad_active=True):
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
        
        Airfoil.__init__(self, OC, Lref=Lref, Sref=Sref, grad_active=grad_active)
        
        self.__AoA0_deg  = None
        self.__AoA0      = None
        self.set_AoA0(AoA0)
        
        self.__Ka        = None
        self.set_Ka(Ka)

        self.__rel_thick      = None
        self.__rel_thick_grad = None
        self.set_rel_thick(rel_thick)
        
        self.pcop        = pcop
        
    #-- Setters
    def set_AoA0(self, AoA0_deg):
        self.__AoA0_deg = AoA0_deg
        self.__AoA0     = AoA0_deg*np.pi/180.
        
    def set_Ka(self, Ka):
        self.__Ka = Ka
        
    def set_rel_thick(self, rel_thick):
        self.__rel_thick = rel_thick
        
    def set_rel_thick_grad(self, rel_thick_grad):
        self.__rel_thick_grad = rel_thick_grad
    
    #-- Accessors
    def get_rel_thick(self):
        return self.__rel_thick
    
    def get_rel_thick_grad(self):
        return self.__rel_thick_grad
    
    #-- Methods to compute aero coefficients
    def comp_aero_coeffs(self, AoA, Mach):
        OC = self.get_OC()
        c  = OC.get_c()
        nu = OC.get_nu()
        L     = self.get_Lref()
        sweep = self.get_sweep()
        toc   = self.get_rel_thick()
        toc2  = self.get_rel_thick()/np.cos(sweep)
        Mach_normal = Mach*np.cos(sweep)
        
        if self.is_grad_active():
            dL  = self.get_Lref_grad()
            dsweep = self.get_sweep_grad()
            dtoc   = self.get_rel_thick_grad()
            dtoc2   = (dtoc*np.cos(sweep)+toc*np.sin(sweep)*dsweep)/np.cos(sweep)**2
            dMach_normal = -Mach*np.sin(sweep)*dsweep

        
        #-- compute dCl_dAoA         
        # Base slope with Prandtl correction
        if Mach_normal<0.9:
            base_slope=2.*np.pi/np.sqrt(1.-Mach_normal**2)
            if self.is_grad_active():
                dbase_slope = 2.*np.pi*(Mach_normal*dMach_normal)/((1.-Mach_normal**2)*np.sqrt(1.-Mach_normal**2))
        elif 0.9<=Mach_normal<=1.1:
            s_sub = 2.*np.pi/np.sqrt(abs(1.-0.9**2))
            s_sup = 4./np.sqrt(abs(1.1**2-1.))
            fact = (1.1-Mach_normal)/0.2
            base_slope = s_sub*fact+s_sup*(1.-fact)
            if self.is_grad_active():
                dfact = -dMach_normal/0.2
                dbase_slope = s_sub*dfact+s_sup*(1.-dfact)
        else:
            base_slope = 4./np.sqrt(Mach_normal**2-1.)
            if self.is_grad_active():
                dbase_slope = -4.*(Mach_normal*dMach_normal)/((1.-Mach_normal**2)*np.sqrt(1.-Mach_normal**2))

        # thickness correction
        thick_corr = 1. + self.THICKNESS_CORRECTION*toc/np.cos(sweep)
        if self.is_grad_active():
            dthick_corr = self.THICKNESS_CORRECTION*(dtoc*np.cos(sweep)+toc*np.sin(sweep)*dsweep)/np.cos(sweep)**2
        # sweep correction
        sweep_corr = np.cos(sweep)**2
        if self.is_grad_active():
            dsweep_corr = -2.*np.sin(sweep)*np.cos(sweep)*dsweep
            
        self.dCl_dAoA = base_slope*thick_corr*sweep_corr
        
        if self.is_grad_active():
            d2Cl_dAoAdchi = dbase_slope*thick_corr*sweep_corr\
                          + base_slope*dthick_corr*sweep_corr\
                          + base_slope*thick_corr*dsweep_corr
                          
        #-- Compute Cl
        self.Cl = self.dCl_dAoA*(AoA-self.__AoA0)
        
        if np.isnan(self.Cl):
            print 'dCl_dAoA = ',self.dCl_dAoA
            print 'sweep_corr = ',sweep_corr
            print 'thick_corr = ',thick_corr
            print 'toc = ',toc
        
        #-- Compute dCl_dchi
        if self.is_grad_active():
            self.dCl_dchi = d2Cl_dAoAdchi*(AoA-self.__AoA0)
    
                      
        #-- Compute Cdw
        Mdd    = self.__Ka/np.cos(sweep) - toc/(np.cos(sweep))**2 - self.Cl/(10.*np.cos(sweep)**3)
        Mcrit  = Mdd - (0.1/80)**(1./3.)
        if Mach < Mcrit:
            self.Cdw = 0.0
        else:
            self.Cdw = 20.*(Mach-Mcrit)**4
            
        
        #-- Compute dCdw_dAoA
        dMdd_dAoA   = -self.dCl_dAoA/(10.*np.cos(sweep)**3)
        dMcrit_dAoA =  dMdd_dAoA
        if Mach < Mcrit:
            self.dCdw_dAoA = 0.0
        else:
            self.dCdw_dAoA = -80.*(Mach-Mcrit)**3*dMcrit_dAoA
    
        #-- compute dCdw_dchi
        if self.is_grad_active():
            dMdd   = self.__Ka*np.sin(sweep)*dsweep/(np.cos(sweep))**2 \
               - dtoc/(np.cos(sweep))**2 - toc*(2.*np.sin(sweep)*dsweep)/np.cos(sweep)**3 \
               - self.dCl_dchi/(10.*np.cos(sweep)**3) - self.Cl*(3.*np.sin(sweep)*dsweep)/(10.*np.cos(sweep)**4)
            dMcrit = dMdd
             
            if Mach < Mcrit:
                self.dCdw_dchi = 0.0
            else:
                self.dCdw_dchi = -80.*(Mach-Mcrit)**3*dMcrit
                
        #-- compute Cdf
        Re = Mach*np.cos(sweep)*c*L/nu
        thick_coeff = 1.+2.1*toc2 # a rough guess since 1.21 for 10% relative thickness
        if Re < 1.e-12:
            # Prevent division by 0. Drag is null at zero Re number anyway
            self.Cdf = 0.
        elif Re < 1.e5:
            # Laminar flow
            self.Cdf=1.328/np.sqrt(Re)*thick_coeff
        else:
            # Turbulent flow
            self.Cdf=0.074*Re**(-0.2)*thick_coeff
            
        #-- Compute dCdf_dAoA
        self.dCdf_dAoA = 0.
        
        #-- Compute dCdf_dchi
        if self.is_grad_active():
            dRe = Mach*np.cos(sweep)*c*dL/nu-Mach*np.sin(sweep)*c*L/nu*dsweep
            dthick_coeff = 2.1*dtoc2
            if Re < 1.e-12:
                # Prevent division by 0. Drag is null at zero Re number anyway
                self.dCdf_dchi = 0.
            elif Re < 1.e5:
                # Laminar flow
                self.dCdf_dchi=-0.664/(Re**1.5)*dRe*thick_coeff+1.328/np.sqrt(Re)*dthick_coeff
            else:
                # Turbulent flow
                self.dCdf_dchi=-0.0148*Re**(-1.2)*dRe*thick_coeff+0.074*Re**(-0.2)*dthick_coeff
                
        #-- compute Cdvp
        Cdvp_min = 60.*toc2**4*self.Cdf
        self.Cdvp = Cdvp_min + 2.*(Cdvp_min + self.Cdf)*self.Cl**2
        #-- compute dCdvp_dAoA
        # dCdvp_min_dAoA = 0. since dCdf_dAoA = 0. 
        self.dCdvp_dAoA = 4.*(Cdvp_min + self.Cdf)*self.Cl*self.dCl_dAoA
        
        if self.is_grad_active():
            dCdvp_min= 60.*(4.*dtoc2*toc2**3*self.Cdf+toc2**4*self.dCdf_dchi)
            self.dCdvp_dchi = dCdvp_min + 2.*(dCdvp_min+self.dCdf_dchi)*self.Cl**2+4.*(Cdvp_min + self.Cdf)*self.Cl*self.dCl_dchi
            
    def get_scaled_copy(self, OC=None, Sref=None, Lref=None, rel_thick=None, grad_active=True):
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        if rel_thick is None:
            rel_thick = self.get_rel_thick()
        if OC is None:
            OC = self.get_OC()
        return AnalyticAirfoil(OC, self.__AoA0_deg, Sref=Sref, Lref=Lref, rel_thick=rel_thick, Ka=self.__Ka, grad_active=self.is_grad_active())
