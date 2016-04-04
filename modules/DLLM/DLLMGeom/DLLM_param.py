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
# @author : Matthieu MEAUX
#

# - Local imports -
import string
from DLLM.DLLMGeom.DLLM_Geom import DLLM_Geom
from MDOTools.Controllers.Kernel.BCManager import BCManager
import numpy as np

class DLLM_param(DLLM_Geom):

    ERROR_MSG = 'ERROR in DLLM_param.'

    def __init__(self, tag, n_sect=20, grad_active=True):
        """
        Constructor: set main attributes
        """
        DLLM_Geom.__init__(self, tag, n_sect=n_sect, grad_active=grad_active)
        self.BC_manager = BCManager()
        
        self.__AoA_id = 'AoA'
        
    #-- Accessors
    def get_dv_array(self):
        return self.BC_manager.get_dv_array()

    def get_dv_id_list(self):
        return self.BC_manager.get_dv_id_list()

    def get_dv_info_list(self):
        return self.BC_manager.get_dv_info_list()

    def get_bounds_array(self):
        return self.BC_manager.get_bounds_array()
    
    #-- Setters
    def set_AoA_id(self, AoA_id):
        self.__AoA_id = AoA_id
    
    def set_value(self, Id, val):
        pt = self.__BC_manager.get_pt(Id)
        pt.set_value(val)
    
    #-- Methods
    def import_BC_from_file(self, filename):
        self.BC_manager.import_from_file(filename)
        self.BC_manager.update()
    
    def update_from_x_list(self, x):
        # print 'update wing_param with x=',x
        self.BC_manager.update_dv_from_x_list(x)
        self.update()
        
    def update(self):
        #-- Design variables update
        self.BC_manager.update()
        ndv = self.BC_manager.get_ndv()
        self.set_ndv(ndv)
        self.__update_AoA()
    
        #-- Update necessary information
        self.build_r_lists()
        
        #-- Build discretization
        self.build_discretization()
        
        DLLM_Geom.update(self)
        
    def build_discretization(self):
        print 'WARNING: build_discretization method has to be overloaded in child classes...'

    #-- Private methods
    def __update_AoA(self):
        Id = self.__AoA_id
        if Id in self.get_dv_id_list():
            pt = self.__BC_manager.get_pt(Id)
            deg_to_rad = np.pi / 180.
            AoA = pt.get_value() * deg_to_rad
            AoA_grad = pt.get_gradient() * deg_to_rad
        else:
            AoA = None
            AoA_grad = None
        self.set_AoA(AoA)
        self.set_AoA_grad(AoA_grad)
        
    def __repr__(self):
        DLLM_Geom.__repr__(self)
        info_string = '\n*** Wing param information ***'
        info_string += '\n  n_sect       : ' + str(self.get_n_sect())
        info_string += '\n  ndv          : ' + str(self.get_ndv())
        info_string += '\n  airfoil_type : ' + str(self.get_airfoil_type())
        info_string += '\n  --    parameters information section    --\n'
        for Id in self.BC_manager.get_list_id():
            pt = self.BC_manager.get_pt(Id)
            BC_Type = pt.get_BCType()
            if BC_Type == 'Variable':
                value = pt.get_value()
                info_string += "%30s" % Id + \
                    "%20s" % BC_Type + " %24.16e" % value + "\n"
            if BC_Type == 'DesignVariable':
                value = pt.get_value()
                bounds = pt.get_bounds()
                info_string += "%30s" % Id + "%20s" % BC_Type + \
                    " %24.16e" % value + " %24s" % str(bounds) + "\n"
            if BC_Type == 'Parameter':
                expr = pt.get_expr()
                info_string += "%30s" % Id + \
                    "%20s" % BC_Type + " %24s" % expr + "\n"
        info_string += '  -- end of parameters information section --\n'
        return info_string
    
#--   Methods to update with new way of working for parmaterization
#     def config_from_dict(self, OC, config_dict):
#         self.build_wing()
#         self.__config_airfoils(OC, config_dict)
#         self.__config_desc(config_dict)
#         self.update()
# 
#     def __config_airfoils(self, OC, config_dict):
#         in_keys_list = config_dict.keys()
# 
#         airfoil_type_key = self.get_tag() + '.airfoil.type'
#         if airfoil_type_key in in_keys_list:
#             airfoil_type = config_dict[airfoil_type_key]
#         else:
#             airfoil_type = 'simple'
#             
#         self.set_airfoil_type(airfoil_type)
#         
#         if airfoil_type == 'simple':
#             AoA0_key = self.get_tag() + '.airfoil.AoA0'
#             if AoA0_key in in_keys_list:
#                 AoA0 = config_dict[AoA0_key]
#             else:
#                 AoA0 = 0.
# 
#             self.build_linear_airfoil(OC, AoA0=AoA0, set_as_ref=True)
# 
#         elif airfoil_type == 'meta':
#             surrogate_model_key = self.__tag + '.airfoil.surrogate_model'
#             if surrogate_model_key in in_keys_list:
#                 surrogate_model = config_dict[surrogate_model_key]
#             else:
#                 surrogate_model = None
# 
#             self.build_meta_airfoil(OC, surrogate_model, set_as_ref=True)
# 
#         self.build_airfoils_from_ref()
# 
#     def __config_desc(self, config_dict):
#         in_keys_list=sorted(config_dict.keys())
#         existing_keys=deepcopy(self.__BC_manager.get_list_id())
#         
#         # Add user defined DesignVariable, Variable and Parameter
#         added_list=[]        
#         for in_key in in_keys_list:
#             words=in_key.split('.')
#             if len(words) >=4:
#                 test=string.join(words[:-2],'.')
#                 if test==self.get_tag()+'.desc':
#                     name=words[-2]
#                     Id=name
#                     Id_in=self.get_tag()+'.desc.'+name
#                     type=config_dict[Id_in+'.type']
#                     if Id not in existing_keys and Id not in added_list:
#                         if   type == 'DesignVariable':
#                             bounds = config_dict[Id_in+'.bounds']
#                             value  = config_dict[Id_in+'.value']
#                             self.__BC_manager.create_design_variable(Id,value,(bounds[0],bounds[1]))
#                         elif type == 'Variable':
#                             value  = config_dict[Id_in+'.value']
#                             self.__BC_manager.create_variable(Id,value)
#                         elif type == 'Parameter':
#                             fexpr = config_dict[Id_in+'.fexpr']
#                             self.__BC_manager.create_parameter(Id,fexpr)
#                         added_list.append(Id)
#         
#         # Convert pre-defined parameters
#         for existing_id in existing_keys:
#             ex_words=existing_id.split('.')
#             name=ex_words[-1]
#             Id_in=self.get_tag()+'.desc.'+name
#             Id=name
#             if Id_in+'.type' in in_keys_list:
#                 type=config_dict[Id_in+'.type']
#                 if   type == 'DesignVariable':
#                     bounds = config_dict[Id_in+'.bounds']
#                     value  = config_dict[Id_in+'.value']
#                     self.set_value(Id,value)
#                     self.convert_to_design_variable(Id,bounds)
#                 elif type == 'Variable':
#                     value  = config_dict[Id_in+'.value']
#                     self.set_value(Id,value)
#                     self.convert_to_variable(Id)
#                 elif type == 'Parameter':
#                     fexpr = config_dict[Id_in+'.fexpr']
#                     self.convert_to_parameter(Id,fexpr)
