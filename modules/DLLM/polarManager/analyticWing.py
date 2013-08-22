# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard

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