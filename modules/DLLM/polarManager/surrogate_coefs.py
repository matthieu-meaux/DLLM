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
#  @author : Andre Mendes
#  @author : Regis Lebrun
#

try:
    from openturns import *
except:
    print 'OpenTURNS not in ENV'
import numpy as np

class SurrogateCoefs:
    """
    Import MetaModel
    Evaluate aerodynamic coefficients through the meta model
    Compute gradients of the meta model functions
    """
    def __init__(self, surrogate_model):
        
        fileName = surrogate_model
        try:
            study = Study()
        except:
            raise Exception,'OpenTURNS not available!'
        study.setStorageManager(XMLStorageManager(fileName))
        study.load()
        
        self.__meta_CL = NumericalMathFunction()
        study.fillObject("MetaModel_y4_deg_5", self.__meta_CL)
        
        self.__meta_CD = NumericalMathFunction()
        study.fillObject("MetaModel_y5_deg_5", self.__meta_CD)
        
        self.__meta_CM = NumericalMathFunction()
        study.fillObject("MetaModel_y6_deg_5", self.__meta_CM)
        
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
    
    
