# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: Francois Gallard
# @author: Matthieu MEAUX (for refactoring)
import numpy
from numpy import array, transpose,outer, ones, zeros, copy, divide, diag, dot
from numpy.linalg import norm, solve

class DLLMMesh:
    """
    Class that deals with geometry for the lifting mine wing solver
    """
    def __init__(self, LLW):
        self.__LLW = LLW
        self.recompute()
    
    #-- Accessors
    def get_airfoils(self):
        return self.__LLW.get_airfoils()
    
    def get_wing_geom(self):
        return self.__LLW.get_wing_geom()
    
    def get_K(self):
        return self.__K
    
    #-- Setters
    def set_Sref(self, Sref):
        self.__LLW.set_Sref(Sref)
        
    def set_Lref(self, Lref):
        self.__LLW.set_Lref(Lref)
    
    #-- Methods
    def recompute(self):
        self.__N   = self.get_wing_geom().get_n_elements()
        
        # compute Lref and Sref from LLW attributes
        self.__set_Lref_Sref()
        
        # Set computational geometry
        self.__K   = None 
        self.__setGeom()
        
    def getAspectRatio(self):
        '''
        Computes the aspect ratio wingspan²/S
        '''
        wingspan = self.get_wing_geom().get_wing_span()
        Sref     = self.get_wing_geom().get_Sref()
        return wingspan**2/Sref
        
    #-- Private methods
    def __set_Lref_Sref(self):
        """
        Compute Lref and Sref from airfoils information
        """
        Lref=0.
        Sref=0.
        for airfoil in self.get_airfoils():
            Sref+=airfoil.get_Sref()
            Lref+=airfoil.get_Lref()
        Lref/=self.__N
        
        self.set_Sref(Sref)
        self.set_Lref(Lref)
        
    def __setGeom(self):
        '''
        Sets the geometry of the wing, builds the stiffness geometry matrix
        '''
        eta=self.get_wing_geom().get_eta()[:,1]
        y=self.get_wing_geom().get_XYZ()[:,1]
        
        YminEta=transpose(outer(ones([self.__N+1]),y))-outer(ones([self.__N]),eta)
        self.__K=divide(ones([self.__N,self.__N+1]),YminEta)
        self.__K/=4.*numpy.pi
        
        
#     Why this method is not in the geometry handling ? If more geometrical parameters, what happens?
#     def set_relative_thickness(self,thickness):
#         """
#         Setter for the height of the airfoils
#         """
#         if type(thickness)==type([]):
#             thick=array(thickness)
#         elif type(thickness)==type(array([0.])):
#             thick=thickness
#         else:
#             raise Exception, "Incorrect type for thickness : "+str(type(thickness))
#         
#         self.get_wing_geom().set_relative_thickness(thick)
#         
#         print "LLW set_relative_thickness thick = "+str(thick)
#         for airfoil, thickness in zip(self.__airfoils,thick):
#             airfoil.set_relative_thickness(thickness)

#     def set_twist(self,twistLaw):
#         '''
#         Sets the twist law of the wing
#         @param twistLaw : the twist law
#         '''
#         if type(twistLaw)==type([]):
#             twist=array(twistLaw)
#         elif type(twistLaw)==type(array([0.])):
#             twist=twistLaw
#         else:
#             raise Exception, "Incorrect type for twistLaw : "+str(type(twistLaw))
#         
#         self.get_wing_geom().set_twist(twist)

        
#     What is the use of these methods ??? estimation of the structual stiffness matrix ? Why there ?
#     def __setGeomWeissinger(self):
#         '''
#         Sets the geometry of the wing, builds the stiffness geometry matrix
#         '''
#         eta=self.get_wing_geom().get_eta()[:,1]
#         self.__K=zeros([self.__N,self.__N+1])
#         
#         x=self.get_wing_geom().get_XYZ()[:,0]
#         y=self.get_wing_geom().get_XYZ()[:,1]
#         
#         for i in range(self.__N):
#             for j in range(self.__N+1):
#                 #Weissinger
#                 self.__K[i,j]=(1.+sqrt(1.+((y[i]-eta[j])/x[i])**2))/(y[i]-eta[j])
#  
#         self.__K/=4.*math.pi        
#         
#     def __setGeomWeissinger2(self):
#         '''
#         Sets the geometry of the wing, builds the stiffness geometry matrix
#         '''
#         eta=self.get_wing_geom().get_eta()[:,1]
#         
#         x=self.get_wing_geom().get_XYZ()[:,0]
#         y=self.get_wing_geom().get_XYZ()[:,1]
#         
#         self.__K=zeros([self.__N,self.__N+1])
#         for i in range(self.__N):
#             for j in range(self.__N+1):
#                 #Weissinger
#                 self.__K[i,j]=(1./(y[i]-eta[j,0])**2)*(1.+(x[i]-eta[j,1])/sqrt((x[i]-eta[j,1])**2+(y[i]-eta[j,0])**2))
#  
#         self.__K/=4.*math.pi
    