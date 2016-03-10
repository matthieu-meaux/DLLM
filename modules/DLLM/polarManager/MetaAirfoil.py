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
# @author : Andre Mendes
# @author : Regis Lebrun
#

from airfoil import Airfoil
import numpy as np 
from surrogate_coefs import SurrogateCoefs
import math

class MetaAirfoil(Airfoil):
    """
    A meta model based on 2D CFD data
    Provide lift, drag and moment coefficients
    Provide sensitivities of aerodynamic surrogate function w.r.t. input parameters
    """
    POW_COS  = 1.0
    
    def __init__(self, OC, surrogate_model, relative_thickness=.12, camber=0., Sref=1., Lref=1., sweep=0., surrogate_fcs=None, grad_active=True):

        Airfoil.__init__(self, OC, Sref, Lref, relative_thickness, sweep, camber, grad_active=grad_active)
        if not surrogate_fcs:
            self.__coefs = SurrogateCoefs(surrogate_model)
        else:
            self.__coefs = surrogate_fcs
        self.__surrogate_model = surrogate_model
        
    #-- Note: to be updated with new way of working
    def Cl(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal = Mach*np.cos(sweep)**self.POW_COS
        ToC=self.get_rel_thick()/np.cos(sweep)*100.
        Cl = self.__coefs.meta_Cl(ToC, self.get_camber(), alpha, Mach_normal)
        sweep_corr = np.cos(self.get_sweep())**(2.*self.POW_COS)
        
        return Cl*sweep_corr
    
    def ClAlpha(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)**self.POW_COS
        ToC=self.get_rel_thick()/np.cos(sweep)**self.POW_COS*100.
        dCl_dAoA = self.__coefs.grad_Cl(ToC, self.get_camber(), alpha, Mach_normal)[2]
        sweep_corr = np.cos(self.get_sweep())**(2.*self.POW_COS)
        
        return dCl_dAoA*sweep_corr
    
    def Cd(self, alpha, Mach):
        Cdf = self.__Cd_friction(alpha, Mach)
        Cdp = self.Cdp(alpha, Mach)
        
        return Cdf + Cdp
    
    def Cdp(self, alpha, Mach):
        
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)**self.POW_COS
        ToC=self.get_rel_thick()/np.cos(sweep)**self.POW_COS*100.
        Cd = self.__coefs.meta_Cd(ToC, self.get_camber(), alpha, Mach_normal)  
        
        Cds = self.Cd_spurious(alpha, Mach)
        
        Cdp = max([Cd-Cds,  0.])
        
        return Cdp
    
    def CdpAlpha(self, alpha, Mach):
        
        Cdp = self.Cdp(alpha, Mach)
        
        if Cdp == 0.:
            dCdp_dAoA = 0.
            
        else:
            Cdsalpha = self.CdsAlpha(alpha, Mach)
            
            sweep=self.get_sweep()
            Mach_normal= Mach*np.cos(sweep)**self.POW_COS
            ToC=self.get_rel_thick()/np.cos(sweep)**self.POW_COS*100.
            dCdp_dAoA = self.__coefs.grad_Cd(ToC, self.get_camber(), alpha, Mach_normal)[2] - Cdsalpha
        
        return dCdp_dAoA
        
    def CdAlpha(self, alpha, Mach):
        dCdp_dAoA = self.CdpAlpha(alpha, Mach) 
        dCdf_AoA = self.__dCd_friction_dAoA(alpha, Mach)
        
        return dCdp_dAoA+dCdf_AoA
        
    def Cm(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)**self.POW_COS
        ToC=self.get_rel_thick()/np.cos(sweep)**self.POW_COS*100.
        Cm = self.__coefs.meta_Cd(ToC, self.get_camber(), alpha, Mach_normal)  
        sweep_corr = np.cos(self.get_sweep())**(2.*self.POW_COS)
        
        return Cm*sweep_corr
    
    def dCl_dchi(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal = Mach*np.cos(sweep)**self.POW_COS
        ToC=self.get_rel_thick()/np.cos(sweep)**self.POW_COS*100.
        Cl = self.__coefs.meta_Cl(ToC, self.get_camber(), alpha, Mach_normal)
        
        dsweep=self.get_sweep_grad()
        sweep_corr = np.cos(self.get_sweep())**(2.*self.POW_COS)
        dsweep_corr = -2.*self.POW_COS*np.sin(sweep)*dsweep*np.cos(sweep)**(2.*self.POW_COS-1.)
        
        dCl_dToC= self.gradCl(alpha, Mach)[0]
        dToC_dchi = 100.*(self.get_rel_thick_grad()*np.cos(sweep)**self.POW_COS+self.POW_COS*self.get_rel_thick()*np.sin(sweep)*dsweep*np.cos(sweep)**(self.POW_COS-1.))/np.cos(sweep)**(2*self.POW_COS)
        
        dCl_dMach = self.gradCl(alpha, Mach)[3]
        dMach_dchi = -self.POW_COS*Mach*np.sin(sweep)*dsweep*np.cos(sweep)**(self.POW_COS-1.)
        
    
        return (dCl_dToC*dToC_dchi + dCl_dMach*dMach_dchi)*sweep_corr + Cl*dsweep_corr
        
    def dCd_dchi(self, alpha, Mach):
        dCdf_dchi = self.__dCd_friction_dchi(alpha, Mach)
        dCdp_dchi = self.dCdp_dchi(alpha, Mach)
        
        return dCdf_dchi + dCdp_dchi
        
    def dCdp_dchi(self, alpha, Mach):
        
        Cdp = self.Cdp(alpha, Mach)
        
        if Cdp != 0:
            sweep=self.get_sweep()
            Mach_normal = Mach*np.cos(sweep)**self.POW_COS
            ToC=self.get_rel_thick()/np.cos(sweep)**self.POW_COS*100.
            Cd = self.__coefs.meta_Cd(ToC, self.get_camber(), alpha, Mach_normal)
            
            dsweep=self.get_sweep_grad()
            
            dCdp_dToC = self.gradCd(alpha, Mach)[0]
            dToC_dchi = 100.*(self.get_rel_thick_grad()*np.cos(sweep)**self.POW_COS+self.POW_COS*self.get_rel_thick()*np.sin(sweep)*dsweep*np.cos(sweep)**(self.POW_COS-1.))/np.cos(sweep)**(2*self.POW_COS)
            
            dCdp_dMach = self.gradCd(alpha, Mach)[3]
            dMach_dchi = -self.POW_COS*Mach*np.sin(sweep)*dsweep*np.cos(sweep)**(self.POW_COS-1.)
            
            dCdsdchi = self.dCds_dchi(alpha, Mach)
            
            dCdpdchi = dCdp_dToC*dToC_dchi + dCdp_dMach*dMach_dchi - dCdsdchi
        
        else:
            dCdpdchi = 0.
         
        return dCdpdchi
    
    def Cdf(self, alpha, Mach):
        return self.__Cd_friction(alpha, Mach)
        
    def __Cd_friction(self, alpha, Mach):
        # Re= V_corr*Lref_corr/nu=Mach*cos(sweep)*c*Lref/cos(sweep)/nu = Mach*c*Lref/cos(sweep)
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
            Cdf=1.328/math.sqrt(Re)*thick_coeff
        else:
            # Turbulent flow
            Cdf=0.074*Re**(-0.2)*thick_coeff
            
        return Cdf
    
    def CdfAlpha(self, alpha, Mach):
        return self.__dCd_friction_dAoA(alpha, Mach)
    
    def __dCd_friction_dAoA(self, alpha, Mach):
        return 0.
    
    def  dCdf_dchi(self, alpha, Mach):
        return self. __dCd_friction_dchi(alpha, Mach)
    
    def __dCd_friction_dchi(self, alpha, Mach):
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
        
        if Re < 1.e-6:
            # Prevent division by 0. Drag is null at zero Re number anyway
            dCdf = 0.
        elif Re < 1.e5:
            # Laminar flow
            dCdf=-0.664/(Re**1.5)*dRe*thick_coeff+1.328/math.sqrt(Re)*dthick_coeff
        else:
            # Turbulent flow
            dCdf=-0.0148*Re**(-1.2)*dRe*thick_coeff+0.074*Re**(-0.2)*dthick_coeff
        return dCdf
        
    def Cd_spurious(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal = Mach*np.cos(sweep)
        ToC=self.get_rel_thick()/np.cos(sweep)*100.
        Cl = self.__coefs.meta_Cl(ToC, self.get_camber(), alpha, Mach_normal)
        
        k, coef = self.k_coef()
        Cl0 = self.Cl0()
       
        Cds = k + coef*(Cl - Cl0)**2.
        return Cds
        
    def CdsAlpha(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal = Mach*np.cos(sweep)
        ToC=self.get_rel_thick()/np.cos(sweep)*100.
        
        Cl = self.__coefs.meta_Cl(ToC, self.get_camber(), alpha, Mach_normal)
        Cl_alpha = self.__coefs.grad_Cl(ToC, self.get_camber(), alpha, Mach_normal)[2]
        
        k, coef = self.k_coef()
        Cl0 = self.Cl0()
        
        Cds_alpha = coef*2.*(Cl - Cl0)*Cl_alpha
        
        return Cds_alpha
    
    def dCds_dchi(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal = Mach*np.cos(sweep)
        ToC=self.get_rel_thick()/np.cos(sweep)*100.
        
        Cl = self.__coefs.meta_Cl(ToC, self.get_camber(), alpha, Mach_normal)
        #dCldchi = self.dCl_dchi(alpha, Mach)
        
        dCl_dToC = self.gradCl(alpha, Mach)[0]
        dToC_dchi = 100.*(self.get_rel_thick_grad()*np.cos(sweep)+self.get_rel_thick()*np.sin(sweep))/np.cos(sweep)**2
        
        dCl_dMach = self.gradCl(alpha, Mach)[3]
        dsweep = self.get_sweep_grad()
        dMach_dchi = -Mach*np.sin(sweep)*dsweep
        
        dCldchi = dCl_dToC*dToC_dchi + dCl_dMach*dMach_dchi
        
        k, coef = self.k_coef()
        dk_dchi, d_coef_dchi = self.d_k_coef_dchi()
        
        Cl0 = self.Cl0()
        dCl0dchi = self.dCl0_dchi()
        
        dCdsdchi = dk_dchi + d_coef_dchi*(Cl - Cl0)**2. + 2.*coef*(Cl - Cl0)*(dCldchi - dCl0dchi) 
        
        return dCdsdchi
        
    def k_coef(self):
        
        k = 0.0045235
        coef = 0.03132053
        
        return k, coef 
    
    def d_k_coef_dchi(self):
        return 0., 0.
    
    def Cl0(self):
        
        Cl0 = 0.16372357
        
        return Cl0
        
    def dCl0_dchi(self):
        return 0.
    
    def gradCl(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        ToC=self.get_rel_thick()/np.cos(sweep)*100.
        gradCl = self.__coefs.grad_Cl(ToC, self.get_camber(), alpha, Mach_normal)
        return gradCl
        
    def gradCd(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        ToC=self.get_rel_thick()/np.cos(sweep)*100.
        gradCd = self.__coefs.grad_Cd(ToC, self.get_camber(), alpha, Mach_normal)
        return gradCd
        
    def gradCm(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        ToC=self.get_rel_thick()/np.cos(sweep)*100.
        gradCm = self.__coefs.grad_Cm(ToC, self.get_camber(), alpha, Mach_normal)
        return gradCm
    
#     def dCl_dthickness(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCl_dthic = self.__coefs.grad_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[0]*100.
#         
#         return dCl_dthic*sweep_corr
#         
#     def dCl_dcamber(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCl_dcamb = self.__coefs.grad_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[1]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCl_dcamb*sweep_corr
#         
#     def dCl_dMach(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCl_dM = self.__coefs.grad_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[3]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCl_dM*sweep_corr
#         
#     def dCd_dthickness(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCd_dthic = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[0]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCd_dthic*sweep_corr
#         
#     def dCd_dcamber(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCd_dcamb = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[1]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCd_dcamb*sweep_corr
#         
#     def dCd_dalpha(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCd_dAoA = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[2]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCd_dAoA*sweep_corr
#         
#     def dCd_dMach(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCd_dM = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[3]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCd_dM*sweep_corr
        
#     def dCm_dthickness(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCm_dthic = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[0]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCm_dthic*sweep_corr
#         
#     def dCm_dcamber(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCm_dcamb = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[1]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCm_dcamb*sweep_corr
#         
#     def dCm_dalpha(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCm_dAoA = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[2]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCm_dAoA*sweep_corr
#         
#     def dCm_dMach(self):
#         sweep=self.get_sweep()
#         Mach_normal= Mach*np.cos(sweep)
#         dCm_dM = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[3]
#         sweep_corr = np.cos(self.get_sweep())**2
#         
#         return dCm_dM*sweep_corr
        
    def get_scaled_copy(self, OC=None, Sref=None,Lref=None, rel_thick=None, camber=None, sweep=None):
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        if rel_thick is None:
            rel_thick = self.get_rel_thick()
        if camber is None:
            camber = self.get_camber()
        if sweep is None:
            sweep = self.get_sweep()
        if OC is None:
            OC = self.get_OC()
        return MetaAirfoil(OC, self.__surrogate_model, relative_thickness=rel_thick, camber=camber, Sref=Sref, Lref=Lref, sweep=sweep, surrogate_fcs=self.__coefs)
        
        
