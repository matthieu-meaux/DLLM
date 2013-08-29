# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.1
# @author: François Gallard

# - Local imports -
from airfoil import Airfoil
from OpenDACE.database.Database import Database

from numpy import array
from spline_interpolator import Spline2D_interpolator

class AirfoilPolar(Airfoil):
    DATABASE_AERO_COEFS_TAG='Cl_Cd_Cm'
    DATABASE_DF_TAG='Df'
    DATABASE_ABSCISSA_TAG='AoA_Thickness'
    
    def __init__(self,Sref,Lref,database=None,relative_thickness=0.15,aero_coefs_interpolator=None):
        '''
        Constructor for a polar of a simple wing.
        @param Sref : reference surface
        @param Lref : reference length        
        @param database : the dataBase in which the discrete values are stored, or None
        '''
        Airfoil.__init__(self,Sref,Lref,relative_thickness)
        
        self.__database=database
        self.__fd_database=Database()
        
        if aero_coefs_interpolator is None:
                print "building new spline 2D Interpolator"
                x,y=self.__database.get_all_in_outs(self.DATABASE_ABSCISSA_TAG,self.DATABASE_AERO_COEFS_TAG)
                y=array(y).T
                x=array(x)
                scale=array([1.,1.])
                self.__cl_interpolator=Spline2D_interpolator(array(x),array(y[0,:]),scale)
                self.__cd_interpolator=Spline2D_interpolator(array(x),array(y[1,:]),scale)
                self.__cm_interpolator=Spline2D_interpolator(array(x),array(y[2,:]),scale)
                print "Done."
        else:
            self.__cl_interpolator=aero_coefs_interpolator[0]
            self.__cd_interpolator=aero_coefs_interpolator[1]
            self.__cm_interpolator=aero_coefs_interpolator[2]
            
    
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
        cl=self.__cl_interpolator.f(array([alpha,self.get_relative_thickness()]))
        return cl

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
        cd=self.__cd_interpolator.f(array([alpha,self.get_relative_thickness()]))
        if cd <0.:
            raise Exception, "Negative drag : at aoa = "+str(alpha)+" thick = "+str(self.get_relative_thickness())
        return cd
    
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
        cm=self.__cm_interpolator.f(array([alpha,self.get_relative_thickness()]))
        return cm
    
    def dCl_dthickness(self,alpha,beta=0.0,Mach=0.0):
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
        dcl=self.__cl_interpolator.df(array([alpha,self.get_relative_thickness()]),comp=1)
        return dcl

    def dCd_dthickness(self,alpha,beta=0.0,Mach=0.0):
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
        dcd=self.__cd_interpolator.df(array([alpha,self.get_relative_thickness()]),comp=1)
        return dcd
    
    def dCm_dthickness(self,alpha,beta=0.0,Mach=0.0):
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
        dcm=self.__cm_interpolator.df(array([alpha,self.get_relative_thickness()]),comp=1)
        return dcm
    
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
        dcl=self.__cl_interpolator.df(array([alpha,self.get_relative_thickness()]),comp=0)
        return dcl

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
        dcd=self.__cd_interpolator.df(array([alpha,self.get_relative_thickness()]),comp=0)
        return dcd
    
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
        dcm=self.__cm_interpolator.df(array([alpha,self.get_relative_thickness()]),comp=0)
        return dcm
    
    def get_database(self):
        '''
        Accessor for the database in which the values are stored
        @return : the database
        '''
        return self.__database
        
    def get_Interpolator(self):
        return (self.__cl_interpolator,self.__cd_interpolator,self.__cm_interpolator)
    
    def get_scaled_copy(self,Sref,Lref,relative_thickness=None):
        """
        Instances a copy of this airfoil with different Lref and Sref
        @param Sref : reference surface
        @param Lref : reference length
        """
        thick=relative_thickness
        if thick==None:
            thick=self.get_relative_thickness()
            
        return AirfoilPolar(Sref=Sref,Lref=Lref,database=self.get_database(),relative_thickness=thick\
                            ,aero_coefs_interpolator=self.get_Interpolator())
    