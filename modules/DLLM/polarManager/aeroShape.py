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

# - Local imports -
class AeroShape:

    def __init__(self,Sref,Lref):
        '''
        Constructor
        '''
        self.__Sref = Sref
        self.__Lref = Lref
        self.__Sref_grad = None
        self.__Lred_grad = None
    
    #-- Accessors
    def get_Sref(self):
        return self.__Sref
    
    def get_Lref(self):
        return self.__Lref
    
    def get_Sref_grad(self):
        return self.__Sref_grad
    
    def get_Lref_grad(self):
        return self.__Lref_grad
    
    #-- Setters
    def set_Sref(self,Sref):
        self.__Sref=Sref
        
    def set_Lref(self,Lref):
        self.__Lref=Lref
        
    def set_Sref_grad(self, Sref_grad):
        self.__Sref_grad = Sref_grad
    
    def set_Lref_grad(self, Lref_grad):
        self.__Lref_grad = Lref_grad
    
    #-- Methods
    def Cl(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return 0.0 
    
    def Cd(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return 0.0 
    
    def Cm(self,alpha,beta=0.0,Mach=0.0):
        '''
        Pitch moment coefficient
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return 0.0 
    
    def LoD(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return self.Cl(alpha,beta)/self.Cd(alpha,beta)
    
    def lift(self,alpha,rho,V,c=0.0,beta=0.0):
        '''
        Lift  function
        @param alpha: angle of Attack 
        @param rho air volumic mass
        @param V : true air speed
        @param beta : sideslip angle
        '''
        Mach=V/c
        return 0.5*rho*self.__Sref*V*V*self.Cl(alpha,beta,Mach)
    
    def drag(self,alpha,rho,V,c=0.0,beta=0.0):
        '''
        Drag  function
        @param alpha: angle of Attack 
        @param rho air volumic mass
        @param V : true air speed
        @param beta : sideslip angle
        '''
        Mach=V/c
        return 0.5*rho*self.__Sref*V*V*self.Cd(alpha,beta,Mach)
    
    def pitchMoment(self,alpha,rho,V,c=0.0,beta=0.0):
        '''
        Pitch moment
        @param alpha: angle of Attack 
        @param rho air volumic mass
        @param V : true air speed
        @param c : speed of sound
        @param beta : sideslip angle
        '''
        Mach=V/c
        return 0.5*rho*self.__Lref*self.__Sref*V*V*self.Cm(alpha,beta,Mach)
    
