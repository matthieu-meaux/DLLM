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
#
# - Local imports -
from aeroShape import AeroShape

class DifferentiatedAeroShape(AeroShape):
    
    epsilonFD=1e-9
    
    def __init__(self,Sref,Lref):
        '''
        Constructor
        @param Sref : reference surface
        @param Lref : reference length
        '''
        AeroShape.__init__(self,Sref,Lref)
    
    def ClAlpha(self,alpha,beta=0.0,Mach=0.0):
        '''
        Sensibility of Cl to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return (self.Cl(alpha+DifferentiatedAeroShape.epsilonFD,Mach=Mach)-self.Cl(alpha-DifferentiatedAeroShape.epsilonFD,Mach=Mach))/(2.*DifferentiatedAeroShape.epsilonFD)
    
    def CdAlpha(self,alpha,beta=0.0,Mach=0.0):
        '''
        Sensibility of Cd to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return (self.Cd(alpha+DifferentiatedAeroShape.epsilonFD,Mach=Mach)-self.Cd(alpha-DifferentiatedAeroShape.epsilonFD,Mach=Mach))/(2.*DifferentiatedAeroShape.epsilonFD)
    
    def CmAlpha(self,alpha,beta=0.0,Mach=0.0):
        '''
        Sensibility of Cm to alpha.
        Done by finites differences by default.
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        return (self.Cm(alpha+DifferentiatedAeroShape.epsilonFD,Mach=Mach)-self.Cm(alpha-DifferentiatedAeroShape.epsilonFD,Mach=Mach))/(2.*DifferentiatedAeroShape.epsilonFD)
    

