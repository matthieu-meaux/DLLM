# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard

# - Local imports -
class AeroShape:

    def __init__(self,Sref,Lref):
        '''
        Constructor
        @param Sref : reference surface
        @param Lref : reference length
        '''
        self.__Sref=Sref
        self.__Lref=Lref
    
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
    
    def getSref(self):
        '''
        Accessor for reference surface
        @return :Sref : the reference surface
        '''
        return self.__Sref
    
    def getLref(self):
        '''
        Accessor for reference length
        @return Lref: the reference length
        '''
        return self.__Lref
    
    def setSref(self,Sref):
        '''
        Setter for reference surface
        @param Sref : the reference surface
        '''
        self.__Sref=Sref
        
    def setLref(self,Lref):
        '''
        Setter for reference length
        @param Lref : the reference length
        '''
        self.__Lref=Lref