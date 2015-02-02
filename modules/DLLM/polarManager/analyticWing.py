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
from aeroShape import AeroShape

class AnalyticWing(AeroShape):
    
    def __init__(self,Sref, Lref, alpha0,kInducedDrag,Cd0,ClAlpha):
        '''
        Constructor for wings
        @param alpha0:angle of attack of null lift
        @param kInducedDrag: induced drag factor
        @param Cx0: skin friction drag
        @param CzAlpha: lift gradient
        '''
        AeroShape.__init__(self,Sref,Lref)
        self.__alpha0=alpha0
        self.__kInducedDrag=kInducedDrag
        self.__Cd0=Cd0
        self.__ClAlpha=ClAlpha
    
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
        return self.__ClAlpha*(alpha-self.__alpha0)
    
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
        return self.__Cd0+self.__kInducedDrag*self.Cl(alpha)**2