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

import string
from numpy import zeros, dot
from numpy import sin, cos
import matplotlib.pylab as plt

class DLLMPost:

    DEF_F_LIST_NAMES = [
        'Cl',
        'Cdi',
        'Cdvp',
        'Cdw',
        'Cdf',
        'Cdp',
        'Cd',
        'Lift',
        'Drag_Induced',
        'Drag_Viscous_pressure',
        'Drag_Wave',
        'Drag_Friction',
        'Drag_Pressure',
        'Drag',
        'LoD',
        'Sref']
    
    ERROR_MSG = 'ERROR in DLLMPost.'

    def __init__(self, LLW, verbose = 1):
        """
        Post-Processing module compatible with adjoint for DLLM
        """
        self.__LLW               = LLW
        self.__F_list_names_calc = self.DEF_F_LIST_NAMES

        self.__F_list_calc = None
        self.__F_list_names = self.DEF_F_LIST_NAMES
        self.__F_list_dim = 0
        self.__F_list = None
        self.__dpF_list_dpiAoA = None
        self.__dpF_list_dpchi = None
        self.__dpF_list_dpAoA = None
        self.__dpF_list_dpthetaY = None

        self.__Lift_distrib_res = None
        self.__dpLift_distrib_res_dpiAoA = None
        self.__dpLift_distrib_res_dpAoA = None
        
        #-- basic post-processing to recompute all information
        self.Cl              = None
        self.dpCl_dpiAoA     = None
        self.dpCl_dpAoA      = None
        self.dpCl_dpthethaY  = None
        self.dpCl_dpchi      = None
        
        self.Cdi              = None
        self.dpCdi_dpiAoA     = None
        self.dpCdi_dpAoA      = None
        self.dpCdi_dpthethaY  = None
        self.dpCdi_dpchi      = None
        
        self.Cdw              = None
        self.dpCdw_dpiAoA     = None
        self.dpCdw_dpAoA      = None
        self.dpCdw_dpthethaY  = None
        self.dpCdw_dpchi      = None
        
        self.Cdvp             = None
        self.dpCdvp_dpiAoA    = None
        self.dpCdvp_dpAoA     = None
        self.dpCdvp_dpthethaY = None
        self.dpCdvp_dpchi     = None
        
        self.Cdf              = None
        self.dpCdf_dpiAoA     = None
        self.dpCdf_dpAoA      = None
        self.dpCdf_dpthethaY  = None
        self.dpCdf_dpchi      = None
  
        self.Lift             = None
        self.dpLift_dpiAoA     = None
        self.dpLift_dpAoA      = None
        self.dpLift_dpthethaY  = None
        self.dpLift_dpchi      = None
        
        self.Drag             = None
        self.dpDrag_dpiAoA     = None
        self.dpDrag_dpAoA      = None
        self.dpDrag_dpthethaY  = None
        self.dpDrag_dpchi      = None
        
        self.LoD              = None
        self.dpLoD_dpiAoA     = None
        self.dpLoD_dpAoA      = None
        self.dpLoD_dpthethaY  = None
        self.dpLoD_dpchi      = None
        
        self.Cl_distrib       = None
        self.Cdi_distrib      = None
        self.Cdw_distrib      = None
        self.Cdvp_distrib     = None
        self.Cdf_distrib      = None

        self.Lift_distrib              = None
        self.dpLift_distrib_dpiAoA     = None
        self.dpLift_distrib_dpAoA      = None
        self.dpLift_distrib_dpthethaY  = None
        self.dpLift_distrib_dpchi      = None
        
        self.Drag_distrib              = None
        self.dpDrag_distrib_dpiAoA     = None
        self.dpDrag_distrib_dpAoA      = None
        self.dpDrag_distrib_dpthethaY  = None
        self.dpDrag_distrib_dpchi      = None

        self.__computed = False
        self.__target_loads = None
        
        self.__verbose = verbose
        
    #-- computed related methods
    def is_computed(self):
        return self.__computed

    def set_computed(self, comp=True):
        self.__computed = comp

    #-- Accessors
    def get_tag(self):
        return self.__LLW.get_tag()
    
    def get_grad_active(self):
        return self.__LLW.get_grad_active()

    def get_N(self):
        return self.get_geom().get_n_sect()

    def get_ndv(self):
        return self.get_geom().get_ndv()

    def get_F_list_names(self):
        return self.__F_list_names

    def get_F_list(self):
        return self.__F_list

    def get_dpF_list_dpiAoA(self):
        return self.__dpF_list_dpiAoA

    def get_dpF_list_dpchi(self):
        return self.__dpF_list_dpchi

    def get_dpF_list_dpAoA(self):
        return self.__dpF_list_dpAoA

    def get_dpF_list_dpthetaY(self):
        return self.__dpF_list_dpthetaY
    
    def get_geom(self):
        return self.__LLW.get_geom()

    def get_airfoils(self):
        return self.__LLW.get_airfoils()

    def get_OC(self):
        return self.__LLW.get_OC()

    def get_localAoA(self):
        return self.__LLW.get_localAoA()

    def get_dplocalAoA_dpiAoA(self):
        return self.__LLW.get_dplocalAoA_dpiAoA()

    def get_dplocalAoA_dpAoA(self):
        return self.__LLW.get_dplocalAoA_dpAoA()

    def get_dplocalAoA_dpchi(self):
        return self.__LLW.get_dplocalAoA_dpchi()

    def get_dplocalAoA_dpthetaY(self):
        return self.__LLW.get_dplocalAoA_dpthetaY()

    def get_iAoA(self):
        return self.__LLW.get_iAoA()

    def get_Lref(self):
        return self.__LLW.get_Lref()

    def get_Sref(self):
        return self.__LLW.get_Sref()

    def get_Sref_grad(self):
        return self.__LLW.get_Sref_grad()

    #-- Setters
    def set_F_list_names(self, F_list_names):
        self.__F_list_names = F_list_names
        
#  -- To update target loads capability
#     def set_target_loads_file(self, loads_file):
#         N = self.get_N()
# 
#         loads = zeros(N)
# 
#         fid = open(loads_file, 'r')
#         all_lines = fid.readlines()
#         fid.close()
#         i = 0
#         for line in all_lines:
#             words = string.split(line)
#             if len(words) == 2:
#                 loads[i] = eval(words[1])
#                 i += 1
# 
#         self.set_target_loads(loads)
#     def set_target_loads(self, loads):
#         self.__target_loads = loads
# 
#         if 'target_loads' not in self.__F_list_names_calc:
#             self.__F_list_names_calc.append('target_loads')
# 
#         if 'target_loads' not in self.__F_list_names:
#             self.__F_list_names.append('target_loads')
    
    #-- Run method
    def __init_run(self):
        N = self.get_N()
        ndv = self.get_ndv()
        self.__F_list_calc_dim = len(self.__F_list_names_calc)
        self.__F_list_calc     = zeros(self.__F_list_calc_dim)
        self.__F_list_dim      = len(self.__F_list_names)
        self.__F_list          = zeros(self.__F_list_dim)
        self.__dpF_list_dpAoA  = zeros(self.__F_list_dim)
        self.__dpF_list_dpiAoA = zeros((self.__F_list_dim, N))
        if self.get_grad_active():
            self.__dpF_list_dpchi  = zeros((self.__F_list_dim, ndv))
        else:
            self.__dpF_list_dpchi  = None
        self.__dpF_list_dpthetaY = zeros((self.__F_list_dim, N))

    def run(self, F_list_names=None):
        ERROR_MSG = self.ERROR_MSG + 'run: '
        N = self.get_N()
        grad_active = self.get_grad_active()
        self.__init_run()
        self.__basic_analysis()

        # Default analysis
        Pdyn     = self.get_OC().get_Pdyn()
        Sref     = self.get_Sref()
        if grad_active:
            Sref_grad = self.get_Sref_grad()
            
        for i, F_name in enumerate(self.__F_list_names_calc):
            if F_name == 'Cl':
                val = self.Cl
            elif F_name == 'Cdi':
                val = self.Cdi
            elif F_name == 'Cdvp':
                val = self.Cdvp
            elif F_name == 'Cdw':
                val = self.Cdw
            elif F_name == 'Cdf':
                val = self.Cdf
            elif F_name == 'Cdp':
                val = self.Cdi+self.Cdw+self.Cdvp
            elif F_name == 'Cd':
                val = self.Cdi+self.Cdw+self.Cdvp+self.Cdf
            elif F_name == 'Lift':
                val = Pdyn*Sref*self.Cl
            elif F_name == 'Drag_Induced':
                val = Pdyn*Sref*self.Cdi
            elif F_name == 'Drag_Viscous_pressure':
                val = Pdyn*Sref*self.Cdvp
            elif F_name == 'Drag_Wave':
                val = Pdyn*Sref*self.Cdw
            elif F_name == 'Drag_Friction':
                val = Pdyn*Sref*self.Cdf
            elif F_name == 'Drag_Pressure':
                val = Pdyn*Sref*(self.Cdi+self.Cdw+self.Cdvp)
            elif F_name == 'Drag':
                val = Pdyn*Sref*(self.Cdi+self.Cdw+self.Cdvp+self.Cdf)
            elif F_name == 'LoD':
                Cl = self.Cl
                Cd = (self.Cdi+self.Cdw+self.Cdvp+self.Cdf)
                val = Cl/Cd
                self.LoD = Cl/Cd
            elif F_name == 'Sref':
                val = Sref
            else:
                raise Exception(ERROR_MSG + ' unknown function ' + str(F_name))
            self.__F_list_calc[i] = val

        # Adjoint analysis
        if self.__verbose > 0:
            if grad_active :
                print "Post : partial derivatives for gradient assembly"
            else:
                print "Post : partial derivatives for other applications"
        for i, F_name in enumerate(self.__F_list_names):
            if F_name == 'Cl':
                val = self.Cl
                dpFdpiAoA   = self.dpCl_dpiAoA
                dpFdpAoA    = self.dpCl_dpAoA
                dpFdpthetaY = self.dpCl_dpthethaY	
                if grad_active:
                    dpFdpchi  =  self.dpCl_dpchi
            elif F_name == 'Cdi':
                val = self.Cdi
                dpFdpiAoA   = self.dpCdi_dpiAoA
                dpFdpAoA    = self.dpCdi_dpAoA
                dpFdpthetaY = self.dpCdi_dpthethaY    
                if grad_active:
                    dpFdpchi  =  self.dpCdi_dpchi
            elif F_name == 'Cdvp':
                val = self.Cdvp
                dpFdpiAoA   = self.dpCdvp_dpiAoA
                dpFdpAoA    = self.dpCdvp_dpAoA
                dpFdpthetaY = self.dpCdvp_dpthethaY    
                if grad_active:
                    dpFdpchi  =  self.dpCdvp_dpchi
            elif F_name == 'Cdw':
                val = self.Cdw
                dpFdpiAoA   = self.dpCdw_dpiAoA
                dpFdpAoA    = self.dpCdw_dpAoA
                dpFdpthetaY = self.dpCdw_dpthethaY    
                if grad_active:
                    dpFdpchi  =  self.dpCdw_dpchi
            elif F_name == 'Cdf':
                val = self.Cdf
                dpFdpiAoA   = self.dpCdf_dpiAoA
                dpFdpAoA    = self.dpCdf_dpAoA
                dpFdpthetaY = self.dpCdf_dpthethaY    
                if grad_active:
                    dpFdpchi  =  self.dpCdf_dpchi
            elif F_name == 'Cdp':
                val = (self.Cdi+self.Cdvp+self.Cdw)
                dpFdpiAoA   = (self.dpCdi_dpiAoA+self.dpCdvp_dpiAoA+self.dpCdw_dpiAoA)
                dpFdpAoA    = (self.dpCdi_dpAoA+self.dpCdvp_dpAoA+self.dpCdw_dpAoA)
                dpFdpthetaY = (self.dpCdi_dpthethaY+self.dpCdvp_dpthethaY+self.dpCdw_dpthethaY)    
                if grad_active:
                    dpFdpchi  =  (self.dpCdi_dpchi+self.dpCdvp_dpchi+self.dpCdw_dpchi)
            elif F_name == 'Cd':
                val = (self.Cdi+self.Cdvp+self.Cdw+self.Cdf)
                dpFdpiAoA   = (self.dpCdi_dpiAoA+self.dpCdvp_dpiAoA+self.dpCdw_dpiAoA+self.dpCdf_dpiAoA)
                dpFdpAoA    = (self.dpCdi_dpAoA+self.dpCdvp_dpAoA+self.dpCdw_dpAoA+self.dpCdf_dpAoA)
                dpFdpthetaY = (self.dpCdi_dpthethaY+self.dpCdvp_dpthethaY+self.dpCdw_dpthethaY+self.dpCdf_dpthethaY)    
                if grad_active:
                    dpFdpchi  =  (self.dpCdi_dpchi+self.dpCdvp_dpchi+self.dpCdw_dpchi+self.dpCdf_dpchi)
            elif F_name == 'Lift':
                self.Lift             = Pdyn*Sref*self.Cl
                self.dpLift_dpiAoA    = Pdyn*Sref*self.dpCl_dpiAoA
                self.dpLift_dpAoA     = Pdyn*Sref*self.dpCl_dpAoA
                self.dpLift_dpthethaY = Pdyn*Sref*self.dpCl_dpthethaY    
                if grad_active:
                    self.dpLift_dpchi  =  Pdyn*(Sref*self.dpCl_dpchi+Sref_grad*self.Cl)
                
                val         = self.Lift 
                dpFdpiAoA   = self.dpLift_dpiAoA  
                dpFdpAoA    = self.dpLift_dpAoA 
                dpFdpthetaY = self.dpLift_dpthethaY  
                if grad_active:
                    dpFdpchi  =  self.dpLift_dpchi
            elif F_name == 'Drag_Induced':
                val         = Pdyn*Sref*self.Cdi
                dpFdpiAoA   = Pdyn*Sref*self.dpCdi_dpiAoA
                dpFdpAoA    = Pdyn*Sref*self.dpCdi_dpAoA
                dpFdpthetaY = Pdyn*Sref*self.dpCdi_dpthethaY    
                if grad_active:
                    dpFdpchi  =  Pdyn*(Sref*self.dpCdi_dpchi+Sref_grad*self.Cdi)
            elif F_name == 'Drag_Viscous_pressure':
                val         = Pdyn*Sref*self.Cdvp
                dpFdpiAoA   = Pdyn*Sref*self.dpCdvp_dpiAoA
                dpFdpAoA    = Pdyn*Sref*self.dpCdvp_dpAoA
                dpFdpthetaY = Pdyn*Sref*self.dpCdvp_dpthethaY    
                if grad_active:
                    dpFdpchi  =  Pdyn*(Sref*self.dpCdvp_dpchi+Sref_grad*self.Cdvp)
            elif F_name == 'Drag_Wave':
                val         = Pdyn*Sref*self.Cdw
                dpFdpiAoA   = Pdyn*Sref*self.dpCdw_dpiAoA
                dpFdpAoA    = Pdyn*Sref*self.dpCdw_dpAoA
                dpFdpthetaY = Pdyn*Sref*self.dpCdw_dpthethaY    
                if grad_active:
                    dpFdpchi  =  Pdyn*(Sref*self.dpCdw_dpchi+Sref_grad*self.Cdw)
            elif F_name == 'Drag_Friction':
                val         = Pdyn*Sref*self.Cdf
                dpFdpiAoA   = Pdyn*Sref*self.dpCdf_dpiAoA
                dpFdpAoA    = Pdyn*Sref*self.dpCdf_dpAoA
                dpFdpthetaY = Pdyn*Sref*self.dpCdf_dpthethaY    
                if grad_active:
                    dpFdpchi  =  Pdyn*(Sref*self.dpCdf_dpchi+Sref_grad*self.Cdf)
            elif F_name == 'Drag_Pressure':
                val = Pdyn*Sref*(self.Cdi+self.Cdvp+self.Cdw)
                dpFdpiAoA   = Pdyn*Sref*(self.dpCdi_dpiAoA+self.dpCdvp_dpiAoA+self.dpCdw_dpiAoA)
                dpFdpAoA    = Pdyn*Sref*(self.dpCdi_dpAoA+self.dpCdvp_dpAoA+self.dpCdw_dpAoA)
                dpFdpthetaY = Pdyn*Sref*(self.dpCdi_dpthethaY+self.dpCdvp_dpthethaY+self.dpCdw_dpthethaY)    
                if grad_active:
                    dpFdpchi  =  Pdyn*(Sref*(self.dpCdi_dpchi+self.dpCdvp_dpchi+self.dpCdw_dpchi)+Sref_grad*(self.Cdi+self.Cdvp+self.Cdw))
            elif F_name == 'Drag':
                self.Drag             = Pdyn*Sref*(self.Cdi+self.Cdvp+self.Cdw+self.Cdf)
                self.dpDrag_dpiAoA    = Pdyn*Sref*(self.dpCdi_dpiAoA+self.dpCdvp_dpiAoA+self.dpCdw_dpiAoA+self.dpCdf_dpiAoA)
                self.dpDrag_dpAoA     = Pdyn*Sref*(self.dpCdi_dpAoA+self.dpCdvp_dpAoA+self.dpCdw_dpAoA+self.dpCdf_dpAoA)
                self.dpDrag_dpthethaY = Pdyn*Sref*(self.dpCdi_dpthethaY+self.dpCdvp_dpthethaY+self.dpCdw_dpthethaY+self.dpCdf_dpthethaY) 
                if grad_active:
                    self.dpDrag_dpchi = Pdyn*(Sref*(self.dpCdi_dpchi+self.dpCdvp_dpchi+self.dpCdw_dpchi+self.dpCdf_dpchi)+Sref_grad*(self.Cdi+self.Cdvp+self.Cdw+self.Cdf))
                
                val = self.Drag 
                dpFdpiAoA   = self.dpDrag_dpiAoA 
                dpFdpAoA    = self.dpDrag_dpAoA 
                dpFdpthetaY = self.dpDrag_dpthethaY   
                if grad_active:
                    dpFdpchi  = self.dpDrag_dpchi 
            elif F_name == 'LoD':
                Cd = (self.Cdi+self.Cdvp+self.Cdw+self.Cdf)
                dpCddpiAoA   = (self.dpCdi_dpiAoA+self.dpCdvp_dpiAoA+self.dpCdw_dpiAoA+self.dpCdf_dpiAoA)
                dpCddpAoA    = (self.dpCdi_dpAoA+self.dpCdvp_dpAoA+self.dpCdw_dpAoA+self.dpCdf_dpAoA)
                dpCddpthetaY = (self.dpCdi_dpthethaY+self.dpCdvp_dpthethaY+self.dpCdw_dpthethaY+self.dpCdf_dpthethaY)    
                if grad_active:
                    dpCddpchi  =  (self.dpCdi_dpchi+self.dpCdvp_dpchi+self.dpCdw_dpchi+self.dpCdf_dpchi)
                
                self.LoD    = self.Cl/Cd
                self.dpLoD_dpiAoA    = (self.dpCl_dpiAoA*Cd-self.Cl*dpCddpiAoA)/Cd**2
                self.dpLoD_dpAoA     = (self.dpCl_dpAoA*Cd-self.Cl*dpCddpAoA)/Cd**2
                self.dpLoD_dpthethaY = (self.dpCl_dpthethaY*Cd-self.Cl*dpCddpthetaY)/Cd**2
                if grad_active:
                    self.dpLoD_dpchi  = (self.dpCl_dpchi*Cd-self.Cl*dpCddpchi)/Cd**2
                
                val         = self.LoD
                dpFdpiAoA   = self.dpLoD_dpiAoA
                dpFdpAoA    = self.dpLoD_dpAoA  
                dpFdpthetaY = self.dpLoD_dpthethaY 
                if grad_active:
                    dpFdpchi  = self.dpLoD_dpchi 
                    
            elif F_name == 'Sref':
                val = self.get_Sref()
                dpFdpiAoA = zeros(N)
                dpFdpAoA  = 0.
                dpFdpthetaY = zeros(N)
                if grad_active:
                    dpFdpchi  = self.get_Sref_grad()
            else:
                raise Exception,ERROR_MSG+' unknown function '+str(F_name)
            
            self.__F_list[i]            = val
            self.__dpF_list_dpiAoA[i,:] = dpFdpiAoA[:]
            self.__dpF_list_dpAoA[i]    = dpFdpAoA
            self.__dpF_list_dpthetaY[i,:] = dpFdpthetaY[:]
            if grad_active:
                self.__dpF_list_dpchi[i,:]  = dpFdpchi[:]

        if self.__verbose > 0 :
            self.__display_info()
        self.set_computed(True)
        
    #-- basic analysis
    def __basic_analysis(self):
        grad_active = self.get_grad_active()
        N          = self.get_N()  
        iAoA       = self.get_iAoA()
        airfoils   = self.get_airfoils()
        dplocalAoA_dpiAoA   = self.get_dplocalAoA_dpiAoA()
        dplocalAoA_dpAoA    = self.get_dplocalAoA_dpAoA()
        dplocalAoA_dpthetaY = self.get_dplocalAoA_dpthetaY()
        Pdyn     = self.get_OC().get_Pdyn()
        Sref     = self.get_Sref()
        
        if grad_active:
            ndv = self.get_ndv()
            dlAoAdchi = self.get_dplocalAoA_dpchi()
            Sref_grad = self.get_Sref_grad()
        
        self.Cl_distrib   = zeros(N)
        self.Cdi_distrib  = zeros(N)
        self.Cdw_distrib  = zeros(N)
        self.Cdvp_distrib = zeros(N)
        self.Cdf_distrib  = zeros(N)
        self.Lift_distrib = zeros(N)
        self.Drag_distrib = zeros(N)
        self.Cl    = 0.0
        self.Cdi   = 0.0
        self.Cdw   = 0.0
        self.Cdvp  = 0.0
        self.Cdf   = 0.0
        
        self.dpCl_dpiAoA    = zeros(N)
        self.dpCdi_dpiAoA   = zeros(N)
        self.dpCdw_dpiAoA   = zeros(N)
        self.dpCdvp_dpiAoA  = zeros(N)
        self.dpCdf_dpiAoA   = zeros(N)
        self.dpLift_distrib_dpiAoA  = zeros((N,N))
        self.dpDrag_distrib_dpiAoA  = zeros((N,N))
                  
        self.dpCl_dpAoA      = 0.0
        self.dpCdi_dpAoA     = 0.0
        self.dpCdw_dpAoA     = 0.0
        self.dpCdvp_dpAoA    = 0.0
        self.dpCdf_dpAoA     = 0.0
        self.dpLift_distrib_dpAoA    = zeros(N)
        self.dpDrag_distrib_dpAoA    =zeros(N)
        
        self.dpCl_dpthethaY   = zeros(N)
        self.dpCdi_dpthethaY  = zeros(N)
        self.dpCdw_dpthethaY  = zeros(N)
        self.dpCdvp_dpthethaY = zeros(N)
        self.dpCdf_dpthethaY  = zeros(N)
        self.dpLift_distrib_dpthethaY  = zeros((N,N))
        self.dpDrag_distrib_dpthethaY  = zeros((N,N))   
        
        if grad_active:
            self.dpCl_dpchi   = zeros(ndv)
            self.dpCdi_dpchi  = zeros(ndv)
            self.dpCdw_dpchi  = zeros(ndv)
            self.dpCdvp_dpchi = zeros(ndv)
            self.dpCdf_dpchi  = zeros(ndv)
            self.dpLift_distrib_dpchi = zeros((N,ndv))
            self.dpDrag_distrib_dpchi = zeros((N,ndv))

        #Note: all airfoils are computed already during direct run
        for i in xrange(N):
            af = airfoils[i]
            surf_fact = af.get_Sref() / Sref
            if grad_active:
                dsurf_fact = (af.get_Sref_grad()*Sref-af.get_Sref()*Sref_grad)/Sref**2
                
            # Coefficients distributions
            self.Cl_distrib[i]   =  af.Cl * cos(iAoA[i])
            self.Cdi_distrib[i]  = -af.Cl * sin(iAoA[i])
            self.Cdw_distrib[i]  =  af.Cdw
            self.Cdvp_distrib[i] =  af.Cdvp
            self.Cdf_distrib[i]  =  af.Cdf
            
            self.Lift_distrib[i] = Pdyn*af.get_Sref()*af.Cl*cos(iAoA[i])
            locCd = -af.Cl*sin(iAoA[i])+af.Cdw+af.Cdvp+af.Cdf
            self.Drag_distrib[i] = Pdyn*af.get_Sref()*locCd
            
            # Coefficients
            self.Cl   += af.Cl * cos(iAoA[i]) * surf_fact
            self.Cdi  -= af.Cl * sin(iAoA[i]) * surf_fact
            self.Cdw  += af.Cdw  * surf_fact
            self.Cdvp += af.Cdvp * surf_fact
            self.Cdf  += af.Cdf  * surf_fact
            
            # Partial derivatives with respect to iAoA
            dlocCd_dAoA = -af.dCl_dAoA*sin(iAoA[i])+af.dCdw_dAoA+af.dCdvp_dAoA+af.dCdf_dAoA
            self.dpCl_dpiAoA[i]  -= af.Cl * sin(iAoA[i]) * surf_fact
            self.dpCl_dpiAoA     += af.dCl_dAoA * cos(iAoA[i]) * dplocalAoA_dpiAoA[i,:] * surf_fact
            
            self.dpCdi_dpiAoA[i] -= af.Cl * cos(iAoA[i]) * surf_fact
            self.dpCdi_dpiAoA    -= af.dCl_dAoA * sin(iAoA[i]) * dplocalAoA_dpiAoA[i,:] * surf_fact
            
            self.dpCdw_dpiAoA    += af.dCdw_dAoA  * dplocalAoA_dpiAoA[i,:] * surf_fact
            
            self.dpCdvp_dpiAoA   += af.dCdvp_dAoA * dplocalAoA_dpiAoA[i,:] * surf_fact
            
            self.dpCdf_dpiAoA    += af.dCdf_dAoA  * dplocalAoA_dpiAoA[i,:] * surf_fact
            
            self.dpLift_distrib_dpiAoA[i,i] -= Pdyn*af.get_Sref()*af.Cl*sin(iAoA[i])
            self.dpLift_distrib_dpiAoA[i,:] += Pdyn*af.get_Sref()*af.dCl_dAoA * cos(iAoA[i]) * dplocalAoA_dpiAoA[i,:]
            
            self.dpDrag_distrib_dpiAoA[i,i] -= Pdyn*af.get_Sref()*af.Cl*cos(iAoA[i])
            self.dpDrag_distrib_dpiAoA[i,:] += Pdyn*af.get_Sref()*dlocCd_dAoA*dplocalAoA_dpiAoA[i,:]
            
            # Partial derivative with respect to AoA
            self.dpCl_dpAoA   += af.dCl_dAoA * cos(iAoA[i]) * dplocalAoA_dpAoA[i] * surf_fact
            self.dpCdi_dpAoA  -= af.dCl_dAoA * sin(iAoA[i]) * dplocalAoA_dpAoA[i] * surf_fact
            self.dpCdw_dpAoA  += af.dCdw_dAoA * dplocalAoA_dpAoA[i] * surf_fact
            self.dpCdvp_dpAoA += af.dCdvp_dAoA * dplocalAoA_dpAoA[i] * surf_fact
            self.dpCdf_dpAoA  += af.dCdf_dAoA * dplocalAoA_dpAoA[i] * surf_fact
            self.dpLift_distrib_dpAoA[i] = Pdyn*af.get_Sref()*af.dCl_dAoA *cos(iAoA[i])*dplocalAoA_dpAoA[i]
            self.dpDrag_distrib_dpAoA[i] = Pdyn*af.get_Sref()*dlocCd_dAoA*dplocalAoA_dpAoA[i]
            
            # Partial derivatives with respect to thetaY
            self.dpCl_dpthethaY   += af.dCl_dAoA * cos(iAoA[i]) * dplocalAoA_dpthetaY[i,:] * surf_fact
            self.dpCdi_dpthethaY  -= af.dCl_dAoA * sin(iAoA[i]) * dplocalAoA_dpthetaY[i,:] * surf_fact
            self.dpCdw_dpthethaY  += af.dCdw_dAoA  * dplocalAoA_dpthetaY[i,:] * surf_fact
            self.dpCdvp_dpthethaY += af.dCdvp_dAoA * dplocalAoA_dpthetaY[i,:] * surf_fact
            self.dpCdf_dpthethaY  += af.dCdf_dAoA  * dplocalAoA_dpthetaY[i,:] * surf_fact
            self.dpLift_distrib_dpthethaY[i,:] = Pdyn*af.get_Sref()*af.dCl_dAoA*cos(iAoA[i])*dplocalAoA_dpthetaY[i,:]
            self.dpDrag_distrib_dpthethaY[i,:] = Pdyn*af.get_Sref()*dlocCd_dAoA*dplocalAoA_dpthetaY[i,:]
            
            # Partial derivatives with respect to chi
            if grad_active:
                self.dpCl_dpchi   += cos(iAoA[i])*(af.dCl_dAoA*dlAoAdchi[i, :]*surf_fact+af.dCl_dchi*surf_fact+af.Cl*dsurf_fact)
                self.dpCdi_dpchi  -= sin(iAoA[i])*(af.dCl_dAoA*dlAoAdchi[i, :]*surf_fact+af.dCl_dchi*surf_fact+af.Cl*dsurf_fact)
                self.dpCdw_dpchi  += (af.dCdw_dAoA*dlAoAdchi[i, :]*surf_fact+af.dCdw_dchi*surf_fact+af.Cdw*dsurf_fact)
                self.dpCdvp_dpchi += (af.dCdvp_dAoA*dlAoAdchi[i, :]*surf_fact+af.dCdvp_dchi*surf_fact+af.Cdvp*dsurf_fact)
                self.dpCdf_dpchi  += (af.dCdf_dAoA*dlAoAdchi[i, :]*surf_fact+af.dCdf_dchi*surf_fact+af.Cdf*dsurf_fact)
                self.dpLift_distrib_dpchi[i,:] = Pdyn*cos(iAoA[i])*(af.get_Sref()*af.dCl_dAoA *dlAoAdchi[i,:]+af.get_Sref()*af.dCl_dchi+af.get_Sref_grad()*af.Cl)

                dlocCd_dchi = -af.dCl_dchi*sin(iAoA[i])+af.dCdw_dchi+af.dCdvp_dchi+af.dCdf_dchi
                self.dpDrag_distrib_dpchi[i,:] = Pdyn*(af.get_Sref()*dlocCd_dAoA*dlAoAdchi[i, :]+af.get_Sref()*dlocCd_dchi+af.get_Sref_grad()*locCd)

    #-- Export methods
    def export_F_list(self, filename=None):    
        if filename is None:
            filename = self.get_tag()+'_F_list.dat'
        
        fid = open(filename,'w')
        for i,F_name in enumerate(self.__F_list_names_calc):
            unit=None
            if   F_name[0:4] in ['Lift','Drag']:
                unit = 'N'
            elif F_name[0:4]=='Sref':
                unit = 'm**2'
            else:
                unit = '-'
            line = F_name+' = '+str(self.__F_list_calc[i])+' ('+unit+') \n'
            fid.write(line)
        fid.close()
    
    #-- Display methods
    def plot(self):
        name = self.get_tag()
        Y_list = self.get_geom().get_XYZ()[1,:]
        
        #-- Plot Lift distrib
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('Local Lift')
        plt.plot(Y_list,self.Lift_distrib)
        plt.rc("font", size=14)
        plt.savefig(name+"_Lift_distrib.png",format='png')
        plt.close()
        
        #-- Plot Drag distrib
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('Local Drag')
        plt.plot(Y_list,self.Drag_distrib)
        plt.rc("font", size=14)
        plt.savefig(name+"_Drag_distrib.png",format='png')
        plt.close()
        
        #-- Plot Cl distrib
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('Local Cl')
        plt.plot(Y_list,self.Cl_distrib)
        plt.rc("font", size=14)
        plt.savefig(name+"_Cl_distrib.png",format='png')
        plt.close()
        
        #-- Plot Cdi distrib
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('Local Cdi')
        plt.plot(Y_list,self.Cdi_distrib)
        plt.rc("font", size=14)
        plt.savefig(name+"_Cdi_distrib.png",format='png')
        plt.close()
        
        #-- Plot Cdvp distrib
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('Local Cdvp')
        plt.plot(Y_list,self.Cdvp_distrib)
        plt.rc("font", size=14)
        plt.savefig(name+"_Cdvp_distrib.png",format='png')
        plt.close()
        
        #-- Plot Cdw distrib
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('Local Cdw')
        plt.plot(Y_list,self.Cdw_distrib)
        plt.rc("font", size=14)
        plt.savefig(name+"_Cdw_distrib.png",format='png')
        plt.close()
        
        #-- Plot Cdf distrib
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('Local Cdf')
        plt.plot(Y_list,self.Cdf_distrib)
        plt.rc("font", size=14)
        plt.savefig(name+"_Cdf_distrib.png",format='png')
        plt.close()
    
    def __display_info(self):
        print self.get_OC()
        print '\n*** aerodynamic functions and coefficients ***'
        print '  Sref  = ', self.get_Sref() ,'[m**2]'
        print '  Lref  = ', self.get_Lref() , '[m]'
        for i, func in enumerate(self.__F_list_names_calc):
            if func [0:4] in ['Lift', 'Drag']:
                unit = '[N]'
            else:
                unit = '[-]'
            if func == 'Sref':
                pass
            else:
                print '  ' + self.__F_list_names_calc[i] + '\t=\t' + str(self.__F_list_calc[i]) + ' ' + unit
