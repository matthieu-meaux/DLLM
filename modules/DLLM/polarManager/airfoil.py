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
    
    def __init__(self, OC, Lref=0., Sref=0., grad_active=True):
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
        
        self.__grad_active = grad_active
        
        # Airfoil aerodynamic coefficients and gradient
        self.Cl          = None
        self.dCl_dAoA    = None
        self.dCl_dchi    = None
         
        self.Cdw         = None
        self.dCdw_dAoA   = None
        self.dCdw_dchi   = None
         
        self.Cdvp        = None
        self.dCdvp_dAoA  = None
        self.dCdvp_dchi  = None
        
        self.Cdf         = None
        self.dCdf_dAoA   = None
        self.dCdf_dchi   = None
        
        self.pcop        = None
        self.dpcop_dAoA  = None
        self.dpcop_dchi  = None
        
        self.gamma       = None
        self.dgamma_dAoA = None
        self.dgamma_dchi = None
        
    #-- Setters
    def set_Lref(self, Lref):
        self.__Lref = Lref
    
    def set_Lref_grad(self, Lref_grad):
        self.__Lref_grad = Lref_grad
    
    def set_Sref(self, Sref):
        self.__Sref = Sref
        
    def set_Sref_grad(self, Sref_grad):
        self.__Sref_grad = Sref_grad
        
    def set_grad_active(self, grad_active):
        self.__grad_active = grad_active
        
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
    
    def is_grad_active(self):
        return self.__grad_active
    
    #-- Aerodynamic methods - to be overloaded by child classes
    #-- Compute the aerodynamic coefficients and partial derivatives
    def compute(self, AoA, Mach):
        self.comp_aero_coeffs(AoA, Mach)
        self.comp_gamma_infos()
    
    def comp_aero_coeffs(self, AoA, Mach):
        """
        method to compute aero coefficients to be overloaded in sub classes
        """
        pass
    
    def comp_gamma_infos(self):
        """
        Compute gamma and gamma derivatives
        """
        L = self.get_Lref()
        self.gamma       = 0.5*L*self.Cl
        self.dgamma_dAoA = 0.5*L*self.dCl_dAoA
        if self.is_grad_active():
            dL = self.get_Lref_grad()
            self.dgamma_dchi = 0.5*(dL*self.Cl+L*self.dCl_dchi)
    
    #-- scaled copy (useless in this case)
    def get_scaled_copy(self, Sref, Lref):
        OC = self.get_OC()
        return Airfoil(OC, Sref, Lref, grad_active=self.__grad_active)  
    
    def print_coeffs(self):
        print '\n*** Airfoil aerodynamic oefficients ***'
        print '  Cl   = ', self.Cl , '[-]'
        print '  Cdw  = ', self.Cdw, '[-]'
        print '  Cdvp = ', self.Cdvp, '[-]'
        print '  Cdf  = ', self.Cdf, '[-]'


