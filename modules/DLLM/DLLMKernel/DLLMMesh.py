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
# @author : Francois Gallard
# @author : Matthieu MEAUX
#

import numpy
from numpy import array, transpose, outer, ones, zeros, copy, divide, diag, dot
from numpy.linalg import norm, solve

class DLLMMesh:
    """
    Class that deals with geometry for the lifting mine wing solver
    """
    def __init__(self, LLW,verbose = 0):
        self.__LLW = LLW
        self.__verbose = verbose
        self.__ndv = self.get_geom().get_ndv()
        self.__N       = None
        self.__K       = None
        self.__dK_dchi = None
        self.recompute()
    
    #-- Accessors
    def get_airfoils(self):
        return self.__LLW.get_airfoils()
    
    def get_tag(self):
        return self.__LLW.get_tag()
    
    def get_geom(self):
        return self.__LLW.get_geom()
    
    def get_OC(self):
        return self.__LLW.get_OC()
    
    def get_grad_active(self):
        return self.__LLW.get_grad_active()
    
    def get_K(self):
        return self.__K
    
    def get_dK_dchi(self):
        return self.__dK_dchi
    
    #-- Methods
    def recompute(self):
        self.__N   = self.get_geom().get_n_sect()
        
        # Set computational geometry
        self.__K       = None 
        self.__dK_dchi = None
        self.__setGeom()
         
    def __setGeom(self):
        '''
        Sets the geometry of the wing, builds the lifting line metric matrix
        '''
        #V   = self.get_OC().get_V()
        eta = self.get_geom().get_eta()[1,:]
        y   = self.get_geom().get_XYZ()[1,:]
        
        YminEta=transpose(outer(ones([self.__N+1]),y))-outer(ones([self.__N]),eta)
        Kmetric=divide(ones([self.__N,self.__N+1]),YminEta)
        Kmetric/=4.*numpy.pi
        
        DdGammaDy_DGamma = zeros([self.__N+1,self.__N])
        DdGammaDy_DGamma[0:self.__N,:]   = diag(ones([self.__N]))
        DdGammaDy_DGamma[self.__N,:]     = 0.0
        DdGammaDy_DGamma[1:self.__N+1,:]-= diag(ones([self.__N]))
        self.__K = - dot(Kmetric,DdGammaDy_DGamma)
        
        
        if self.get_grad_active():
            eta_grad = self.get_geom().get_eta_grad()[1,:,:]
            y_grad   = self.get_geom().get_XYZ_grad()[1,:,:]
            
            YminEta_grad=zeros((self.__N,self.__N+1,self.__ndv))
            for n in xrange(self.__ndv):
                YminEta_grad[:,:,n] = transpose(outer(ones([self.__N+1]),y_grad[:,n]))-outer(ones([self.__N]),eta_grad[:,n])
            
            dKmetric_dchi=zeros((self.__N,self.__N+1,self.__ndv))
    #         for n in xrange(self.__ndv):
    #             for j in xrange(self.__N+1):
    #                 for i in xrange(self.__N):
    #                     dKmetric_dchi[i,j,n]=-YminEta_grad[i,j,n]/YminEta[i,j]**2
    #                     
            for n in xrange(self.__ndv):
                dKmetric_dchi[:,:,n]=-YminEta_grad[:,:,n]/YminEta[:,:]**2
            dKmetric_dchi/=4.*numpy.pi
            
            self.__dK_dchi = zeros((self.__N,self.__N,self.__ndv))
            for n in xrange(self.__ndv):
                self.__dK_dchi[:,:,n]=-dot(dKmetric_dchi[:,:,n],DdGammaDy_DGamma)
                  

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
    