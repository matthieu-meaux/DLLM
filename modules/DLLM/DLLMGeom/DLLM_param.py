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
        
        self.__AoA_id    = 'AoA'
        
        self.__common_OC = None
        
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
        
    def set_common_OC(self, OC):
        self.__common_OC = OC
    
    #-- Methods
    def import_BC_from_file(self, filename):
        self.BC_manager.import_from_file(filename)
        self.BC_manager.update()
        
    def config_from_dict(self, config_dict):
        config_dict_keys = config_dict.keys()
        #-- Import BCManager from filename defined in config_dict
        BC_filename_key=self.get_tag()+'.BCfilename'
        if BC_filename_key in config_dict_keys:
            BC_filename = config_dict[BC_filename_key]
            self.import_BC_from_file(BC_filename)
            
        #TBC: to be update with new way of working
        #-- Config airfoils
        airfoil_type_key = self.get_tag() + '.airfoil.type'
        if airfoil_type_key in config_dict_keys:
            airfoil_type = config_dict[airfoil_type_key]
        else:
            airfoil_type = 'simple'
             
        self.set_airfoil_type(airfoil_type)
         
        if airfoil_type == 'simple':
            AoA0_key = self.get_tag() + '.airfoil.AoA0'
            if AoA0_key in config_dict_keys:
                AoA0 = config_dict[AoA0_key]
            else:
                AoA0 = 0.
 
            self.build_linear_airfoil(self.__common_OC, AoA0=AoA0, set_as_ref=True)
 
        elif airfoil_type == 'meta':
            surrogate_model_key = self.__tag + '.airfoil.surrogate_model'
            if surrogate_model_key in config_dict_keys:
                surrogate_model = config_dict[surrogate_model_key]
            else:
                surrogate_model = None
 
            self.build_meta_airfoil(self.__common_OC, surrogate_model, set_as_ref=True)
 
        self.build_airfoils_from_ref()
        self.update()
    
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
            pt = self.BC_manager.get_pt(Id)
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
