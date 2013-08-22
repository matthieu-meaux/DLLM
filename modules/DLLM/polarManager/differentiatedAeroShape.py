# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard

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
    

