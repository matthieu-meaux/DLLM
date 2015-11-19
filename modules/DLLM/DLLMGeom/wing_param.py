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
from numpy import zeros
from numpy import pi, sqrt, cos, sin
from copy import deepcopy

class Wing_param(DLLM_Geom):

    ERROR_MSG = 'ERROR in Wing_param.'
    POSSIBLE_GEOM_TYPES = ["Rectangular", "Elliptic", "Broken"]
    DISCRETE_ATTRIBUTES_LIST = [
        'span',
        'sweep',
        'break_percent',
        'root_chord',
        'break_chord',
        'tip_chord',
        'root_height',
        'break_height',
        'tip_height']

    def __init__(self, tag, geom_type='Broken', n_sect=20):
        """
        Constructor: set main attributes
        """
        DLLM_Geom.__init__(self, tag, n_sect=n_sect)
        
        self.__geom_type = None
        self.__BC_manager = BCManager()

        self.__distrib_type = 'cos_law'

        self.set_geom_type(geom_type)


        self.__AoA_id = 'AoA'

        self.__init_discrete_attributes()

    def __init_discrete_attributes(self):
        # -- discrete attributes
        self.__span = None
        self.__span_grad = None

        self.__sweep = None
        self.__sweep_grad = None

        self.__break_percent = None
        self.__break_percent_grad = None

        self.__root_chord = None
        self.__root_chord_grad = None

        self.__break_chord = None
        self.__break_chord_grad = None

        self.__tip_chord = None
        self.__tip_chord_grad = None

        self.__root_height = None
        self.__root_height_grad = None

        self.__break_height = None
        self.__break_height_grad = None

        self.__tip_height = None
        self.__tip_height_grad = None

        self.__heights = None
        self.__heights_grad = None

        self.__root_ToC = None
        self.__root_ToC_grad = None

        self.__break_ToC = None
        self.__break_ToC_grad = None

        self.__tip_ToC = None
        self.__tip_ToC_grad = None

        self.__AR = None
        self.__AR_grad = None

        self.__fuel = None
        self.__fuel_grad = None

    # -- Accessors
    def get_BC_manager(self):
        return self.__BC_manager

    # -- Accessors needed for pevious studies - can be removed
    def get_span(self):
        return self.__span
 
    def get_span_grad(self):
        return self.__span_grad
 
    def get_sweep(self):
        return self.__sweep
 
    def get_sweep_grad(self):
        return self.__sweep_grad
 
    def get_root_ToC(self):
        return self.__root_ToC
 
    def get_root_ToC_grad(self):
        return self.__root_ToC_grad
 
    def get_break_ToC(self):
        return self.__break_ToC
 
    def get_break_ToC_grad(self):
        return self.__break_ToC_grad
 
    def get_tip_ToC(self):
        return self.__tip_ToC
 
    def get_tip_ToC_grad(self):
        return self.__tip_ToC_grad

    # -- Setters
    def set_AoA_id(self, AoA_id):
        self.__AoA_id = AoA_id

    def set_geom_type(self, geom_type):
        """
        @param wing_geometry_type : Rectangular or Elliptic planforms
        @param type wing_geometry_type : String
        """
        if not geom_type in self.POSSIBLE_GEOM_TYPES:
            raise Exception("geom_type :" +
                            str(geom_type) +
                            " not in possible types : " +
                            str(self.POSSIBLE_GEOM_TYPES))
        self.__geom_type = geom_type



    def set_value(self, Id, val):
        pt = self.__BC_manager.get_pt(Id)
        pt.set_value(val)

    # -- controller management
    def convert_to_variable(self, Id):
        self.__BC_manager.convert_to_variable(Id)

    def convert_to_design_variable(self, Id, bounds):
        self.__BC_manager.convert_to_design_variable(Id, bounds)

    def convert_to_parameter(self, Id, fexpr):
        self.__BC_manager.convert_to_parameter(Id, fexpr)

    def add_AoA_design_variable(self, val, bounds):
        self.__BC_manager.create_design_variable(self.__AoA_id, val, bounds)

    def __update_AoA(self):
        Id = self.__AoA_id
        if Id in self.get_dv_id_list():
            pt = self.__BC_manager.get_pt(Id)
            deg_to_rad = pi / 180.
            AoA = pt.get_value() * deg_to_rad
            AoA_grad = pt.get_gradient() * deg_to_rad
        else:
            AoA = None
            AoA_grad = None
        self.set_AoA(AoA)
        self.set_AoA_grad(AoA_grad)

    def get_dv_array(self):
        return self.__BC_manager.get_dv_array()

    def get_dv_id_list(self):
        return self.__BC_manager.get_dv_id_list()

    def get_dv_info_list(self):
        return self.__BC_manager.get_dv_info_list()

    def get_bounds_array(self):
        return self.__BC_manager.get_bounds_array()

    def build_wing(self):
        self.__BC_manager.clean()
        self.__BC_manager.create_variable('span', 0.)
        if self.__geom_type not in ['Elliptic']:
            self.__BC_manager.create_variable('sweep', 0.)
        if self.__geom_type in ['Rectangular', 'Elliptic']:
            self.__BC_manager.create_variable('root_chord', 0.)
            self.__BC_manager.create_variable('root_height', 0.)
            self.__BC_manager.create_variable('tip_height', 0.)
        elif self.__geom_type == 'Broken':
            self.__BC_manager.create_variable('break_percent', 0.33)
            self.__BC_manager.create_variable('root_chord', 0.)
            self.__BC_manager.create_variable('break_chord', 0.)
            self.__BC_manager.create_variable('tip_chord', 0.)
            self.__BC_manager.create_variable('root_height', 0.)
            self.__BC_manager.create_variable('break_height', 0.)
            self.__BC_manager.create_variable('tip_height', 0.)

        N = self.get_n_sect()

        for i in xrange(N / 2):
            self.__BC_manager.create_design_variable('rtwist' + str(i), 0., (-10., 10.))
        for i in xrange(N / 2):
            # print
            # self.__tag+'.twist'+str(i),self.__tag+'.twist'+str(self.__n_sect/2+i),self.__tag+'.rtwist'+str(self.__n_sect/2-1-i),self.__tag+'.rtwist'+str(i)
            self.__BC_manager.create_parameter(
                'twist' + str(i), 'rtwist' + str(N / 2 - 1 - i))
            self.__BC_manager.create_parameter(
                'twist' + str(N  / 2 + i), 'rtwist' + str(i))
#         for i in xrange(self.__n_sect):
#             self.__BC_manager.create_design_variable(self.__tag+'.twist'+str(i),-25.,0.,25.)

    def config_from_dict(self, OC, config_dict):
        self.build_wing()
        self.__config_airfoils(OC, config_dict)
        self.__config_desc(config_dict)
        self.update()

    def __config_airfoils(self, OC, config_dict):
        in_keys_list = config_dict.keys()

        airfoil_type_key = self.__tag + '.airfoil.type'
        if airfoil_type_key in in_keys_list:
            airfoil_type = config_dict[airfoil_type_key]
        else:
            airfoil_type = 'simple'
            
        self.set_airfoil_type(airfoil_type)
        
        if airfoil_type == 'simple':
            AoA0_key = self.__tag + '.airfoil.AoA0'
            if AoA0_key in in_keys_list:
                AoA0 = config_dict[AoA0_key]
            else:
                AoA0 = 0.

            Cm0_key = self.__tag + '.airfoil.Cm0'
            if Cm0_key in in_keys_list:
                Cm0 = config_dict[Cm0_key]
            else:
                Cm0 = 0.

            self.build_linear_airfoil(OC, AoA0=AoA0, Cm0=Cm0, set_as_ref=True)

        elif airfoil_type == 'meta':
            surrogate_model_key = self.__tag + '.airfoil.surrogate_model'
            if surrogate_model_key in in_keys_list:
                surrogate_model = config_dict[surrogate_model_key]
            else:
                surrogate_model = None

            self.build_meta_airfoil(OC, surrogate_model, set_as_ref=True)

        self.build_airfoils_from_ref()

    def __config_desc(self, config_dict):
        in_keys_list=sorted(config_dict.keys())
        existing_keys=deepcopy(self.__BC_manager.get_list_id())
        
        # Add user defined DesignVariable, Variable and Parameter
        added_list=[]        
        for in_key in in_keys_list:
            words=in_key.split('.')
            if len(words) >=4:
                test=string.join(words[:-2],'.')
                if test==self.__tag+'.desc':
                    name=words[-2]
                    Id=name
                    Id_in=self.__tag+'.desc.'+name
                    type=config_dict[Id_in+'.type']
                    if Id not in existing_keys and Id not in added_list:
                        if   type == 'DesignVariable':
                            bounds = config_dict[Id_in+'.bounds']
                            value  = config_dict[Id_in+'.value']
                            self.__BC_manager.create_design_variable(Id,value,(bounds[0],bounds[1]))
                        elif type == 'Variable':
                            value  = config_dict[Id_in+'.value']
                            self.__BC_manager.create_variable(Id,value)
                        elif type == 'Parameter':
                            fexpr = config_dict[Id_in+'.fexpr']
                            self.__BC_manager.create_parameter(Id,fexpr)
                        added_list.append(Id)
        
        # Convert pre-defined parameters
        for existing_id in existing_keys:
            ex_words=existing_id.split('.')
            name=ex_words[-1]
            Id_in=self.__tag+'.desc.'+name
            Id=name
            if Id_in+'.type' in in_keys_list:
                type=config_dict[Id_in+'.type']
                if   type == 'DesignVariable':
                    bounds = config_dict[Id_in+'.bounds']
                    value  = config_dict[Id_in+'.value']
                    self.set_value(Id,value)
                    self.convert_to_design_variable(Id,bounds)
                elif type == 'Variable':
                    value  = config_dict[Id_in+'.value']
                    self.set_value(Id,value)
                    self.convert_to_variable(Id)
                elif type == 'Parameter':
                    fexpr = config_dict[Id_in+'.fexpr']
                    self.convert_to_parameter(Id,fexpr)

    def update(self):
        #-- Design variables update
        self.__BC_manager.update()
        ndv = self.__BC_manager.get_ndv()
        self.set_ndv(ndv)
        self.__update_AoA()
    
        #-- Update necessary information
        self.build_r_lists()
        
        #-- Build discretization
        self.__build_discretization()
        
        DLLM_Geom.update(self)

    def update_from_x_list(self, x):
        # print 'update wing_param with x=',x
        self.__BC_manager.update_dv_from_x_list(x)
        self.update()

    # -- discretization methods
    def __build_discretization(self):
        self.__get_and_set_parameters()
        self.__compute_ToCs()
        self.__build_chords_eta()
        self.__build_rel_thicks_eta()
        self.__build_sweep_eta()
        self.__build_eta()
        
    def __get_and_set_parameters(self):
        deg_to_rad = pi / 180.
        N   = self.get_n_sect()
        ndv = self.get_ndv()
        # -- span
        Id = 'span'
        pt = self.__BC_manager.get_pt(Id)
        val = pt.get_value()
        grad = pt.get_gradient()
        self.__span = val
        self.__span_grad = grad

        if self.__geom_type not in ['Elliptic']:
            # -- sweep
            Id = 'sweep'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__sweep = val * deg_to_rad
            self.__sweep_grad = grad * deg_to_rad
        else:
            self.__sweep = 0.
            self.__sweep_grad = zeros(self.get_ndv())

        if self.__geom_type in ['Rectangular', 'Elliptic']:
            Id = 'root_chord'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__root_chord = val
            self.__root_chord_grad = grad

            Id = 'root_height'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__root_height = val
            self.__root_height_grad = grad

            Id = 'tip_height'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_height = val
            self.__tip_height_grad = grad

        elif self.__geom_type == 'Broken':
            Id = 'break_percent'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__break_percent = val
            self.__break_percent_grad = grad

            Id = 'root_chord'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__root_chord = val
            self.__root_chord_grad = grad

            Id = 'break_chord'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__break_chord = val
            self.__break_chord_grad = grad

            Id = 'tip_chord'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_chord = val
            self.__tip_chord_grad = grad

            Id = 'root_height'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__root_height = val
            self.__root_height_grad = grad

            Id = 'break_height'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__break_height = val
            self.__break_height_grad = grad

            Id = 'tip_height'
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_height = val
            self.__tip_height_grad = grad

        twist = zeros(N)
        twist_grad = zeros((N, ndv))
        for i in xrange(N):
            Id = 'twist' + str(i)
            pt = self.__BC_manager.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            twist[i] = val * deg_to_rad
            twist_grad[i, :] = grad * deg_to_rad
        self.set_twist(twist)
        self.set_twist_grad(twist_grad)
        
    def __build_sweep_eta(self):
        N   = self.get_n_sect()
        ndv = self.get_ndv()
        
        sweep_eta = zeros(N+1)
        sweep_grad_eta = zeros((N+1,ndv))
        
        sweep_eta[:] = self.__sweep
        sweep_grad_eta[:] = self.__sweep_grad
        
        self.set_sweep_eta(sweep_eta)
        self.set_sweep_grad_eta(sweep_grad_eta)

    def __build_chords_eta(self):
        N   = self.get_n_sect()
        ndv = self.get_ndv()
        r_list_eta = self.get_r_list_eta()

        chords_eta = zeros((N + 1))
        chords_grad_eta = zeros((N+ 1, ndv))

        if self.__geom_type == 'Elliptic':
            for i, r in enumerate(r_list_eta):
                chords_eta[i] = self.__root_chord * \
                    sqrt(1. - (2. * r)**2)
                chords_grad_eta[i, :] = self.__root_chord_grad[:] \
                * sqrt(1. - (2. * r)**2)

        elif self.__geom_type == 'Rectangular':
            for i in xrange(N + 1):
                chords_eta[i] = self.__root_chord
                chords_grad_eta[i, :] = self.__root_chord_grad[:]

        elif self.__geom_type == 'Broken':
            p = self.__break_percent / 100.
            p_grad = self.__break_percent_grad / 100.

            for i, r in enumerate(r_list_eta):
                r = abs(2. * r)
                if r <= p:
                    coeff = r / p
                    dcoeff = -r * p_grad[:] / p**2
                    chords_eta[i] = (
                        self.__break_chord - self.__root_chord) * coeff + self.__root_chord
                    chords_grad_eta[i, :] = (self.__break_chord_grad[:] - self.__root_chord_grad[:]) * coeff  \
                        + (self.__break_chord - self.__root_chord) * dcoeff \
                        + self.__root_chord_grad[:]

                else:
                    coeff = (r - p) / (1. - p)
                    dcoeff = (r - 1) * p_grad[:] / (1. - p)**2
                    chords_eta[i] = (
                        self.__tip_chord - self.__break_chord) * coeff + self.__break_chord
                    chords_grad_eta[i, :] = (self.__tip_chord_grad[:] - self.__break_chord_grad[:]) * coeff  \
                        + (self.__tip_chord - self.__break_chord) * dcoeff \
                        + self.__break_chord_grad[:]
                        
        self.set_chords_eta(chords_eta)
        self.set_chords_grad_eta(chords_grad_eta)

    def __build_rel_thicks_eta(self):
        N   = self.get_n_sect()
        ndv = self.get_ndv()
        r_list_eta = self.get_r_list_eta()
        
        heights_eta = zeros(N+1)
        heights_grad_eta = zeros((N+1, ndv))

        rel_thicks_eta = zeros(N+1)
        rel_thicks_grad_eta = zeros((N+1, ndv))

        if self.__geom_type == 'Broken':
            p = self.__break_percent / 100.
            p_grad = self.__break_percent_grad / 100.

            for i, r in enumerate(r_list_eta):
                r = abs(2. * r)
                if r <= p:
                    coeff = r / p
                    dcoeff = -r * p_grad[:] / p**2
                    heights_eta[i] = (
                        (self.__break_height -
                         self.__root_height) *
                        coeff +
                        self.__root_height)
                    heights_grad_eta[i, :] = (self.__break_height_grad[:] - self.__root_height_grad[:]) * coeff \
                        + (self.__break_height - self.__root_height) * dcoeff \
                        + self.__root_height_grad[:]

                else:
                    coeff = (r - p) / (1. - p)
                    dcoeff = (r - 1) * p_grad[:] / (1. - p)**2
                    heights_eta[i] = (
                        (self.__tip_height -
                         self.__break_height) *
                        coeff +
                        self.__break_height)
                    heights_grad_eta[i, :] = (self.__tip_height_grad[:] - self.__break_height_grad[:]) * coeff \
                        + (self.__tip_height - self.__break_height) * dcoeff \
                        + self.__break_height_grad[:]
        else:
            for i, r in enumerate(r_list_eta):
                r = abs(2. * r)
                heights_eta[i] = (
                    (self.__tip_height -self.__root_height)*r + self.__root_height)
                heights_grad_eta[i,:] = (
                    (self.__tip_height_grad[:] -self.__root_height_grad[:])*r +self.__root_height_grad[:])

        #-- build rel_thiks
        chords_eta = self.get_chords_eta()
        chords_grad_eta = self.get_chords_grad_eta()
        for i in xrange(N+1):
            rel_thicks_eta[i] = heights_eta[i] / chords_eta[i]
            rel_thicks_grad_eta[i,:] = (heights_grad_eta[i,:] * chords_eta[i] - heights_eta[i] * chords_grad_eta[i,:]) / \
                                       (chords_eta[i])**2
                                       
        self.set_rel_thicks_eta(rel_thicks_eta)
        self.set_rel_thicks_grad_eta(rel_thicks_grad_eta)          

    def __build_eta(self):
        N          = self.get_n_sect()
        ndv        = self.get_ndv()
        r_list_eta = self.get_r_list_eta()

        eta      = zeros((3, N + 1))
        eta_grad = zeros((3, N + 1, ndv))
        
        chords_eta      = self.get_chords_eta()
        chords_grad_eta = self.get_chords_grad_eta()

        for i, r in enumerate(r_list_eta):
            abs_r = abs(r)
            eta[0, i] = abs_r * self.__span * sin(self.__sweep) + 0.25 * chords_eta[i]
            eta_grad[0, i, :] = abs_r * self.__span_grad[:] * sin(self.__sweep) + abs_r * self.__span * cos(self.__sweep) * self.__sweep_grad[:] + 0.25 * chords_grad_eta[i,:]

            eta[1, i] = r * self.__span
            eta_grad[1, i, :] = r * self.__span_grad[:]
            
        self.set_eta(eta)
        self.set_eta_grad(eta_grad)
        
    def __compute_ToCs(self):
        """
        Compute ToCs
        """
        self.__root_ToC = self.__root_height / self.__root_chord
        self.__root_ToC_grad = (self.__root_height_grad * self.__root_chord -
                                self.__root_height * self.__root_chord_grad) / self.__root_chord**2

        if self.__break_chord is not None:
            self.__break_ToC = self.__break_height / self.__break_chord
            self.__break_ToC_grad = (self.__break_height_grad * self.__break_chord -
                                     self.__break_height * self.__break_chord_grad) / self.__break_chord**2

        if self.__tip_chord is not None:
            self.__tip_ToC = self.__tip_height / self.__tip_chord
            self.__tip_ToC_grad = (self.__tip_height_grad * self.__tip_chord -
                                   self.__tip_height * self.__tip_chord_grad) / self.__tip_chord**2
                                   
    def __repr__(self):
        info_string = '\n*** Wing param information ***'
        info_string += '\n  geom_type    : ' + str(self.__geom_type)
        info_string += '\n  n_sect       : ' + str(self.get_n_sect())
        info_string += '\n  ndv          : ' + str(self.get_ndv())
        info_string += '\n  airfoil_type : ' + str(self.get_airfoil_type())
        info_string += '\n  --    parameters information section    --\n'
        for Id in self.__BC_manager.get_list_id():
            pt = self.__BC_manager.get_pt(Id)
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
