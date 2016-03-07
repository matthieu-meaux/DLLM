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

# - Local imports -
class Airfoil:
    '''
    Airfoil class for lifting line computations. 
    Supports the computation of the circulation and all partial derivatives necessary for optimisation
    '''
    
    def __init__(self, OC, Lref=0., Sref=0.):
        '''
        Constructor
        @param OC        : Operating condition object to access fluid properties
        @param Sref      : reference surface
        @param Lref      : reference length
        '''
        self.__OC = OC
        
        self.__Lref      = Lref
        self.__Lref_grad = None
        
        self.__Sref      = Sref
        self.__Sref_grad = None
        
        # Airfoil aerodynamic coefficients and gradient
#         self.__Cl         = None
#         self.__dCl_dAoA   = None
#         self.__dCl_dchi   = None
#         
#         self.__Cdw        = None
#         self.__dCdw_dAoA  = None
#         self.__dCdw_dchi  = None
#         
#         self.__Cdvp       = None
#         self.__dCdvp_dAoA = None
#         self.__dCdvp_dchi = None
        
    
    #-- Setters
    def set_Lref(self, Lref):
        self.__Lref = Lref
    
    def set_Lref_grad(self, Lref_grad):
        self.__Lref_grad = Lref_grad
    
    def set_Sref(self, Sref):
        self.__Sref = Sref
        
    def set_Sref_grad(self, Sref_grad):
        self.__Sref_grad = Sref_grad
        
    #-- Accessors
    def get_Lref(self):
        return self.__Lref
    
    def get_Lref_grad(self):
        return self.__Lref_grad
    
    def get_Sref(self):
        return self.__Sref
    
    def get_Sref_grad(self):
        return self.__Sref_grad
    
    def get_OC(self):
        return self.__OC
    
    #-- Aerodynamic methods - to be overloaded by child classes
    #-- Cl related mehods
    def Cl(self, AoA, Mach):
        return None
    
    def dCl_dAoA(self, AoA, Mach):
        return None
    
    def dCl_dchi(self, AoA, Mach):
        return None
    
    #-- Cdw related methods
    def Cdw(self, AoA, Mach):
        return None
    
    def dCdw_dAoA(self, AoA, Mach):
        return None
    
    def dCdw_dchi(self, AoA, Mach):
        return None
    
    #-- Cdvp related methods
    def Cdvp(self, AoA, Mach):
        return None
    
    def dCdvp_dAoA(self, AoA, Mach):
        return None
    
    def dCdvp_dchi(self, AoA, Mach):
        return None
    
    #-- Cdf related methods
    def Cdf(self, AoA, Mach):
        return None
    
    def dCdf_dAoA(self, AoA, Mach):
        return None
    
    def dCdf_dchi(self, AoA, Mach):
        return None
    
    #-- Cd related methods
    def Cd(self, AoA, Mach):
        Cdw  = self.Cdw(AoA, Mach)
        Cdvp = self.Cdvp(AoA, Mach)
        Cdf  = self.Cdf(AoA, Mach)
        Cd = Cdw + Cdvp + Cdf
        return Cd
        
    def dCd_dAoA(self, AoA, Mach):
        dCdw  = self.dCdw_dAoA(AoA, Mach)
        dCdvp = self.dCdvp_dAoA(AoA, Mach)
        dCdf  = self.dCdf_dAoA(AoA, Mach)
        dCd   = dCdw + dCdvp + dCdf
        return dCd
    
    def dCd_dchi(self, AoA, Mach):
        dCdw  = self.dCdw_dchi(AoA, Mach)
        dCdvp = self.dCdvp_dchi(AoA, Mach)
        dCdf  = self.dCdf_dchi(AoA, Mach)
        dCd   = dCdw + dCdvp + dCdf
        return dCd
    
    #-- Cdp related methods
    def Cdp(self, AoA, Mach):
        Cdw  = self.Cdw(AoA, Mach)
        Cdvp = self.Cdvp(AoA, Mach)
        Cdp = Cdw + Cdvp
        return Cdp
        
    def dCdp_dAoA(self, AoA, Mach):
        dCdw  = self.dCdw_dAoA(AoA, Mach)
        dCdvp = self.dCdvp_dAoA(AoA, Mach)
        dCdp  = dCdw + dCdvp
        return dCdp
    
    def dCdp_dchi(self, AoA, Mach):
        dCdw  = self.dCdw_dchi(AoA, Mach)
        dCdvp = self.dCdvp_dchi(AoA, Mach)
        dCdp  = dCdw + dCdvp
        return dCdp
    
    #-- gamma related methods
    def gamma(self, AoA, Mach):
        L     = self.get_Lref()
        Cl    = self.Cl(AoA, Mach)
        gamma = 0.5*L*Cl
        return gamma
    
    def dgamma_dAoA(self, AoA, Mach):
        L = self.get_Lref()
        dCl_dAoA = self.dCl_dAoA(AoA, Mach)
        dgamma_dAoA = 0.5*L*dCl_dAoA
        return dgamma_dAoA
    
    def dgamma_dchi(self, AoA, Mach):
        L   = self.get_Lref()
        dL  = self.get_Lref_grad()
        Cl  = self.Cl(AoA, Mach)
        dCl = self.dCl_dchi(AoA, Mach)
        dgamma_dchi =  0.5*dL*Cl+ 0.5*L*dCl 
        return dgamma_dchi
    
    #-- scaled copy (useless in this case)
    def get_scaled_copy(self, Sref, Lref):
        OC = self.get_OC()
        return Airfoil(OC, Sref, Lref)  
