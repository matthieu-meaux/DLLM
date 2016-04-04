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
from DLLM.DLLMGeom.DLLM_param import DLLM_param
import numpy as np

class Wing_param(DLLM_param):

    ERROR_MSG = 'ERROR in Wing_param.'

    def __init__(self, tag, geom_type='Broken', n_sect=20, grad_active=True):
        """
        Constructor: set main attributes
        """
        DLLM_param.__init__(self, tag, n_sect=n_sect, grad_active=grad_active)
        

    # -- discretization methods
    def __build_discretization(self):
        self.__get_and_set_parameters()
        self.__build_chords_eta()
        self.__build_rel_thicks_eta()
        self.__build_sweep_eta()
        self.__build_eta()
        
    def __get_and_set_parameters(self):
        deg_to_rad = pi / 180.
        N   = self.get_n_sect()
        ndv = self.get_ndv()
        BCM = self.__BC_manager
        # -- span
        Id = 'span'
        pt = BCM.get_pt(Id)
        val = pt.get_value()
        grad = pt.get_gradient()
        self.__span = val
        self.__span_grad = grad

        if self.__geom_type not in ['Elliptic']:
            # -- sweep
            Id = 'sweep'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__sweep = val * deg_to_rad
            self.__sweep_grad = grad * deg_to_rad
        else:
            self.__sweep = 0.
            self.__sweep_grad = zeros(self.get_ndv())

        if self.__geom_type in ['Rectangular', 'Elliptic']:
            Id = 'root_chord'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__root_chord = val
            self.__root_chord_grad = grad

            Id = 'root_height'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__root_height = val
            self.__root_height_grad = grad

            Id = 'tip_height'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_height = val
            self.__tip_height_grad = grad

        elif self.__geom_type == 'Broken':
            Id = 'break_percent'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__break_percent = val
            self.__break_percent_grad = grad

            Id = 'root_chord'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__root_chord = val
            self.__root_chord_grad = grad

            Id = 'break_chord'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__break_chord = val
            self.__break_chord_grad = grad

            Id = 'tip_chord'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_chord = val
            self.__tip_chord_grad = grad

            Id = 'root_height'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__root_height = val
            self.__root_height_grad = grad

            Id = 'break_height'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__break_height = val
            self.__break_height_grad = grad

            Id = 'tip_height'
            pt = BCM.get_pt(Id)
            val = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_height = val
            self.__tip_height_grad = grad
        
        twist = zeros(N)
        twist_grad = zeros((N, ndv))
        for i in xrange(N / 2):
            pt1 = BCM.get_pt('twist'+str(N / 2 - 1 - i))
            twist[i] = pt1.get_value()*deg_to_rad
            twist_grad[i,:] = pt1.get_gradient()*deg_to_rad
            pt2 = BCM.get_pt('twist'+str(i))
            twist[N  / 2 + i] = pt2.get_value()*deg_to_rad
            twist_grad[N  / 2 + i,:] = pt2.get_gradient()*deg_to_rad 
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
            if chords_eta[i] ==0.:
                rel_thicks_eta[i] = 0.
                rel_thicks_grad_eta[i,:] = 0.*chords_eta[i]
            else:
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