# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: EADS Innovation Works
# @version: 1.0
# @author: Regis LEBRUN
# @author: Andre BORTOLAZZI
try:
    from openturns import *
except:
    pass
import numpy as np

class SurrogateCoefs:
    """
    Import MetaModel
    Evaluate aerodynamic coefficients through the meta model
    Compute gradients of the meta model functions
    """
    def __init__(self, surrogate_model):
        
        fileName = surrogate_model
        study = Study()
        study.setStorageManager(XMLStorageManager(fileName))
        study.load()
        
        self.__meta_CL = NumericalMathFunction()
        study.fillObject("MetaModel_y4_deg_6", self.__meta_CL)
        
        self.__meta_CD = NumericalMathFunction()
        study.fillObject("MetaModel_y5_deg_6", self.__meta_CD)
        
        self.__meta_CM = NumericalMathFunction()
        study.fillObject("MetaModel_y6_deg_6", self.__meta_CM)
        
        distribution = Distribution()
        study.fillObject("distribution", distribution)
    
    def meta_Cl(self, thickness,  camber,  AoA, Mach):
        """
        Evaluate lift coefficient using the meta model
        """
        AoA = AoA*180./np.pi
        Cl = self.__meta_CL([thickness,  camber,  AoA, Mach])[0]
        return Cl
        
    def meta_Cd(self, thickness,  camber,  AoA, Mach):
        """
        Evaluate drag coefficient using the meta model
        """
        AoA = AoA*180./np.pi
        Cd = self.__meta_CD([thickness,  camber,  AoA, Mach])[0]
        return Cd
        
    def meta_Cm(self, thickness,  camber,  AoA, Mach):  
        """
        Evaluate moment coefficient using the meta model
        """
        AoA = AoA*180./np.pi
        Cm = self.__meta_CM([thickness,  camber,  AoA, Mach])[0]
        return Cm
    
    def compute_gradient(self, thickness,  camber,  AoA, Mach):
        """
        Compute gradients of all aerodynamic meta functions w.r.t. their entries
        """
        AoA = AoA*180./np.pi
        self.__gradient = np.array((0, 4))
        
        self.__gradient = np.vstack([self.__gradient, self.grad_Cl(thickness,  camber,  AoA, Mach)])
        self.__gradient = np.vstack([self.__gradient, self.grad_Cd(thickness,  camber,  AoA, Mach)])
        self.__gradient = np.vstack([self.__gradient, self.grad_Cm(thickness,  camber,  AoA, Mach)])
        for i in range(len(self.__gradient)):
            self.__gradient[i, 2] = self.__gradient[i, 2]*180./np.pi
        
    def grad_Cl(self, thickness,  camber,  AoA, Mach):
        """
        Evaluate gradient of Cl
        """
        AoA = AoA*180./np.pi
        grad_cl = np.transpose(np.array(self.__meta_CL.gradient([thickness,  camber,  AoA, Mach])))[0]
        grad_cl[2] = grad_cl[2]*180./np.pi
        
        return grad_cl
        
    def grad_Cd(self, thickness,  camber,  AoA, Mach):
        """
        Evaluate gradient of Cd
        """
        AoA = AoA*180./np.pi
        grad_cd = np.transpose(np.array(self.__meta_CD.gradient([thickness,  camber,  AoA, Mach])))[0]
        grad_cd[2] = grad_cd[2]*180./np.pi
        
        return grad_cd
        
        
    def grad_Cm(self, thickness,  camber,  AoA, Mach):
        """
        Evaluate gradient of Cm
        """
        AoA = AoA*180./np.pi
        grad_cm = np.transpose(np.array(self.__meta_CM.gradient([thickness,  camber,  AoA, Mach])))[0]
        grad_cm[2] = grad_cm[2]*180./np.pi
        
        return grad_cm
    
    
