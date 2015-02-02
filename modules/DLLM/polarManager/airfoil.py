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
#  http://github.com/TBD
#
# - Local imports -
from differentiatedAeroShape import DifferentiatedAeroShape
from numpy import sqrt,sin,cos

class Airfoil(DifferentiatedAeroShape):
    '''
    Airfoil class for lifting line computations. 
    Supports the computation of the circulation and the sensibility of moments to AoA.
    '''
    
    def __init__(self, OC, Sref, Lref, rel_thick=0.0, sweep=0.0, camber=0.0):
        '''
        Constructor
        @param Sref : reference surface
        @param Lref : reference length
        '''
        DifferentiatedAeroShape.__init__(self,Sref,Lref)
        self.__OC = OC
        
        self.__rel_thick = rel_thick
        self.__rel_thick_grad = 0.
        
        self.__camber      = camber
        self.__camber_grad = 0.
        
        self.__sweep      = sweep
        self.__sweep_grad = 0.
        
    def set_OC(self, OC):
        self.__OC = OC
    
    def get_OC(self):
        return self.__OC
        
    def set_rel_thick(self, rel_thick):
        self.__rel_thick = rel_thick
        
    def set_rel_thick_grad(self, rel_thick_grad):
        self.__rel_thick_grad = rel_thick_grad
        
    def get_rel_thick(self):
        return self.__rel_thick
    
    def get_rel_thick_grad(self):
        return self.__rel_thick_grad
        
    def set_camber(self, camber):
        self.__camber = camber
        
    def set_camber_grad(self, camber_grad):
        self.__camber_grad = camber_grad
        
    def get_camber(self):
        return self.__camber
    
    def get_camber_grad(self):
        return self.__camber_grad
    
    def set_sweep(self, sweep):
        self.__sweep = sweep
        
    def set_sweep_grad(self, sweep_grad):
        self.__sweep_grad = sweep_grad
        
    def get_sweep(self):
        return self.__sweep
    
    def get_sweep_grad(self):
        return self.__sweep_grad
    
    def gamma(self,alpha,Mach=0.0):
        #V= self.get_OC().get_V()
        sweep=self.get_sweep()
        L=self.get_Lref()
        Cl=self.Cl(alpha,Mach=Mach)
        gamma =  0.5*L*Cl
        return gamma
    
    def dCl_dchi(self,alpha,Mach=0.0):
        return 0.0
    
    def dCd_dchi(self,alpha,Mach=0.0):
        return 0.0    
 
    def dCm_dchi(self,alpha,Mach=0.0):
        return 0.0  
     
    def dpgamma_dpAoA(self,alpha,Mach):
        #V= self.get_OC().get_V()
        sweep=self.get_sweep()
        L=self.get_Lref()
        dCl_dAoA=self.ClAlpha(alpha,Mach)
        dpgamma_dpAoA = 0.5*L*dCl_dAoA
        return dpgamma_dpAoA
    
    def dpgamma_dpchi(self, alpha, Mach):
        #V= self.get_OC().get_V()
        L=self.get_Lref()
        dL=self.get_Lref_grad()
        Cl=self.Cl(alpha,Mach=Mach)
        dCl=self.dCl_dchi(alpha,Mach)
        sweep=self.get_sweep()
        dsweep=self.get_sweep_grad()
        dpgamma_dpchi =  0.5*dL*Cl+ 0.5*L*dCl 
        return dpgamma_dpchi
    
    def get_scaled_copy(self,Sref, Lref, rel_thick=0.0, sweep=0.):
        return Airfoil(Sref,Lref,rel_thick=rel_thick, sweep=sweep)
    

