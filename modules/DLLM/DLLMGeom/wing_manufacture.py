# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard

# - Local imports -
from DLLM.DLLMGeom.discrete_wing import Discrete_wing
from numpy import zeros, array
from math import pi, sqrt

class Wing_manufacture():
    
    GEOMETRY_POSSIBLE_TYPES=["Rectangular","Elliptic","Broken"]
    
    def __init__(self):
        """
        Constructor
        generates numpy arrays to describe discrete wing stored in discrete_wing
        """
        self.__wing_geometry_type = None
        self.__wingspan           = None
        self.__root_chord         = None
        self.__root_height        = None
        self.__tip_height         = None
        self.__break_chord        = None
        self.__tip_chord          = None
        self.__break_percent      = None
        self.__break_height       = None
        self.__discrete_wing      = None

    # -- Setters
    def set_wing_geometry_type(self, wing_geometry_type):
        """
        @param wing_geometry_type : Rectangular or Elliptic planforms
        @param type wing_geometry_type : String
        """
        if not wing_geometry_type in self.GEOMETRY_POSSIBLE_TYPES:
            raise Exception, "wing_geometry_type :"+str(wing_geometry_type)+" not in possible types : "+str(self.GEOMETRY_POSSIBLE_TYPES)
        self.__wing_geometry_type=wing_geometry_type
        
    def set_wingspan(self, wingspan):
        self.__wingspan = wingspan
        
    def set_root_chord(self, root_chord):
        self.__root_chord = root_chord
    
    def set_root_height(self, root_height):
        self.__root_height = root_height
        
    def set_tip_height(self, tip_height):
        self.__tip_height = tip_height
        
    def set_break_chord(self, break_chord):
        self.__break_chord = break_chord
        
    def set_tip_chord(self, tip_chord):
        self.__tip_chord = tip_chord
    
    def set_break_percent(self, break_percent):
        self.__break_percent = break_percent
        
    def set_break_height(self, break_height):
        self.__break_height = break_height
        
    def build_Discrete_wing(self,N,twist_law=None):
        '''
        Builds an elliptic liftingLineWing
        @param reference_airfoil : The reference airfoil from which others will be built
        @param N : the number of elements for the whole wing
        @param twist_law : the twist law
        '''
        if N%2!=0:
            raise Exception, "The total number of elements in the wing must be odd."
        
        n=int(N/2)
        XYZ=zeros([N,3])
        eta=zeros([N+1,3])
        chords=zeros(N)
        heights=zeros(N)
        wbs = zeros(N)  # wingbox length
        wingspan=self.__wingspan
        twist=self.__cast_twist(n,twist_law)
        for i in range(2*n+1):
            eta[i,1]=self.__compute_eta(i,n)
            
        for i in range(N):
            y=self.__compute_y(i,n)
            Lloc=self.__compute_chord(i,n)
            x=Lloc/2.#3quarters chord
            XYZ[i,:]=array([x,y,0])
            chords[i]=Lloc
            heights[i]=self.__compute_h(i,n)
     
        dw=Discrete_wing(XYZ,eta,chords,heights,twist,wingspan)
        
        self.__discrete_wing = dw
        
        return dw

    def __linear(self, x,y_min,y_max,x_min=0.,x_max=1.):
        return (y_max-y_min)*float(x-x_min)/(float(x_max-x_min))+y_min
    
    def __compute_eta(self,i,n):
        ifl=float(i-n)
        wingspan=self.__wingspan
        return (ifl/float(n))*wingspan/2.
    
    def __compute_y(self,i,n):
        r=float(i+0.5-n)/float(n)
        wingspan=self.__wingspan
        return r*wingspan/2.
        
    def __compute_chord(self,i,n):
        root_chord=self.__root_chord
        
        if self.__wing_geometry_type=="Elliptic":
            r=float(i+0.5-n)/float(n)
            return root_chord*sqrt(1.-r**2)
        
        if self.__wing_geometry_type=="Rectangular":
            return root_chord
        
        if self.__wing_geometry_type=="Broken":
            b_root = self.__root_chord
            b_break = self.__break_chord
            b_tip = self.__tip_chord
            p = self.__break_percent/100.
            
            r=abs(float(i+0.5-n)/float(n))
            
            if r<= p:
                return self.__linear(r,b_root,b_break,0.,p)
            else :
                return self.__linear(r,b_break,b_tip,p,1.)
            
    def __compute_h(self,i,n):
        root_chord=self.__root_chord
        
        if self.__wing_geometry_type=="Broken":
            h_root = self.__root_height
            h_break = self.__break_height
            h_tip = self.__tip_height
            p = self.__break_percent/100.
            
            r=abs(float(i+0.5-n)/float(n))
            
            if r<= p:
                return self.__linear(r,h_root,h_break,0.,p)
            else :
                return self.__linear(r,h_break,h_tip,p,1.)
            
        else:
            h_root = self.__root_height
            h_tip = self.__tip_height
            
            r=abs((i+0.5-n)/float(n))
            
            return self.__linear(r,h_root,h_tip)
        
    def __cast_twist(self,n,twist_law):
        if twist_law is None:
            twist=zeros([2*n])
        else:
            if type(twist_law)==type([]):
                twist=array(twist_law)
            elif type(twist_law)==type(array([0.])):
                twist=twist_law
            else:
                raise Exception, "Incorrect type for twistLaw : "+str(type(twist_law))
            
        return twist
