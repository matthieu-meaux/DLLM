# Copyright: Airbus
# @version: 1.0
# @author: Francois Gallard
 
# - Local imports -
from numpy import array, zeros

class Discrete_wing():
    
    def __init__(self,XYZ,eta,chords,heights,twist,wingspan):
        '''
        Constructor for wings based on lifting line theory
        Store numpy arrays generated by wing_manufacture
        @param XYZ : the position of the airfoils
        @param eta : the coordinates of the slices of the wing along Y axis
        @param twist : the twist law of the airfoils
        '''

#TODO        if not (len(airfoils)==len(XYZ)==len(twist)):len(eta)-1
#            raise Exception,"The wing geometric discretisation is not coherent."
        
        if not (type(XYZ)==type(twist)==type(eta)==type(array([0.]))):
            raise Exception,"LiftingLineWing input parameters type must be : "+str(type(array([0.])))
        
        self._XYZ=XYZ
        self._eta=eta
        self._twist=twist
        self._wingspan=wingspan
        self._chords=chords
        self._N=len(twist)
        self._relative_thickness=zeros(self._N)
        self.set_heights(heights)
        
    def set_relative_thickness(self,relative_thickness):
        """
        Setter for the height of the sections
        """
        
        self._relative_thickness=relative_thickness
        self._heights=self._relative_thickness*self._chords
        
    def set_heights(self,heights):
        """
        Setter for the height of the sections
        """

#        self._relative_thickness=heights/self._chords
        self._heights=heights
        
    def get_relative_thickness(self):
        """
        Accessor for the relative thickness.
        """
        return self._relative_thickness
        
    def get_n_elements(self):
        """
        Accessor for the number of elements
        """
        return self._N

    def get_chords(self):
        """
        Accessor for the chords
        """
        return self._chords

    def get_heights(self):
        """
        Accessor for the airfoil heights 
        """
        return self._heights
    
    def get_wingspan(self):
        """
        Accessor for the wingspan
        """
        return self._wingspan

    def get_XYZ(self):
        """
        Accessor for the airfoil positions
        """
        return self._XYZ
    
    def get_eta(self):
        """
        Accessor for the sections positions
        """
        return self._eta
    
    def set_twist(self,twist):
        """
        Setter for the twist law.
        """
        self._twist=twist
            
    def get_twist(self):
        """
        Accessor for the twist law.
        """
        return self._twist

        
    

    