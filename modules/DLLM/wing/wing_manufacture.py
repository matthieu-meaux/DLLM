# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard

# - Local imports -
from AeroElastAdj.liftingLine.liftingLineWing import LiftingLineWing
from AeroElastAdj.differentiatedStructuralModel.Components.wing_structure import Wing_structure
from AeroElastAdj.wing.discrete_wing import Discrete_wing
from numpy import zeros, array
from math import pi, sqrt

class Wing_manufacture():
    
    GEOMETRY_POSSIBLE_TYPES=["Rectangular","Elliptic","Broken"]
    
    def __init__(self,wing_geometry_type):
        """
        Initializes the type of wing that will be built.
        @param wing_geometry_type : Rectangular or Elliptic planforms
        @param type wing_geometry_type : String
        """
        if not wing_geometry_type in self.GEOMETRY_POSSIBLE_TYPES:
            raise Exception, "wing_geometry_type :"+str(wing_geometry_type)+" not in possible types : "+str(self.GEOMETRY_POSSIBLE_TYPES)

        self.__wing_geometry_type=wing_geometry_type
        
    def build_Lifting_line_wing(self,discrete_wing,reference_airfoil):
        airfoils=[]
        n=discrete_wing.get_n_elements()/2
        thick=discrete_wing.get_relative_thickness()
        Lref=0.
        for i in range(2*n):
            Lloc=discrete_wing.get_chords()[i]
            Sloc=self.__compute_local_Sref(Lloc,n,discrete_wing.get_wingspan())
            airfoils.append(reference_airfoil.get_scaled_copy(Sloc,Lloc,thick[i]))
            Lref+=Lloc
            
        Lref/=2*n
        llw=LiftingLineWing(discrete_wing,Lref,airfoils)
        
        return llw
        
    def build_Wing_structure(self,discrete_wing,material_properties,structure_elements_type):
        sw=Wing_structure(discrete_wing, material_properties,structure_elements_type)
        return sw
        
    def build_Discrete_wing(self,N,geom_parameters,twist_law=None):
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
        wingspan=geom_parameters["wingspan"]
        twist=self.__cast_twist(n,twist_law)
        for i in range(2*n+1):
            eta[i,1]=self.__compute_eta(i,n,geom_parameters)
            
        for i in range(N):
            y=self.__compute_y(i,n,geom_parameters)
            Lloc=self.__compute_chord(i,n,geom_parameters)
            x=Lloc/2.#3quarters chord
            XYZ[i,:]=array([x,y,0])
            chords[i]=Lloc
            heights[i]=self.__compute_h(i,n,geom_parameters)
            wbs[i] = chords[i]/geom_parameters['box_chord_ratio']
     
        dw=Discrete_wing(XYZ,eta,chords,heights,wbs,twist,wingspan)
        
        return dw


    
    def __linear(self, x,y_min,y_max,x_min=0.,x_max=1.):
        return (y_max-y_min)*float(x-x_min)/(float(x_max-x_min))+y_min
    
    def __compute_eta(self,i,n,geom_parameters):
        ifl=float(i-n)
        wingspan=geom_parameters["wingspan"]
        return (ifl/float(n))*wingspan/2.
    
    def __compute_y(self,i,n,geom_parameters):
        r=float(i+0.5-n)/float(n)
        wingspan=geom_parameters["wingspan"]
        return r*wingspan/2.
        
    def __compute_chord(self,i,n,geom_parameters):
        root_chord=geom_parameters["root_chord"]
        
        if self.__wing_geometry_type=="Elliptic":
            r=float(i+0.5-n)/float(n)
            return root_chord*sqrt(1.-r**2)
        
        if self.__wing_geometry_type=="Rectangular":
            return root_chord
        
        if self.__wing_geometry_type=="Broken":
            b_root = geom_parameters['root_chord']
            b_break = geom_parameters['break_chord']
            b_tip = geom_parameters['tip_chord']
            p = geom_parameters['break_percent']/100.
            
            r=abs(float(i+0.5-n)/float(n))
            
            if r<= p:
                return self.__linear(r,b_root,b_break,0.,p)
            else :
                return self.__linear(r,b_break,b_tip,p,1.)
            
    def get_chord(self,i,n,geom_parameters):
        chord = self.__compute_chord(i, n, geom_parameters)
        return chord
            
    def __compute_h(self,i,n,geom_parameters):
        root_chord=geom_parameters["root_chord"]
        
        if self.__wing_geometry_type=="Broken":
            h_root = geom_parameters['root_height']
            h_break = geom_parameters['break_height']
            h_tip = geom_parameters['tip_height']
            p = geom_parameters['break_percent']/100.
            
            r=abs(float(i+0.5-n)/float(n))
            
            if r<= p:
                return self.__linear(r,h_root,h_break,0.,p)
            else :
                return self.__linear(r,h_break,h_tip,p,1.)
            
        else:
            h_root = geom_parameters['root_height']
            h_tip = geom_parameters['tip_height']
            
            r=abs((i+0.5-n)/float(n))
            
            return self.__linear(r,h_root,h_tip)
        
    def __compute_local_Sref(self,Lloc,n,wingspan):
        return Lloc*wingspan/(2.*float(n))
    
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
