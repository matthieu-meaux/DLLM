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
import sys
import string
import numpy as np
from DLLM.polarManager.analyticAirfoil import AnalyticAirfoil
from DLLM.polarManager.MetaAirfoil import MetaAirfoil
from numpy import zeros
from numpy import pi, sqrt, cos, sin
from copy import deepcopy
import matplotlib.pylab as plt

class DLLM_Geom(object):
    ERROR_MSG = 'ERROR in DLLM_Geom.'
    POS_DISTRIB = ['linear', 'cos_law']
    
    def __init__(self, tag, n_sect=20, grad_active=False):
        """
        Constructor
        """
        self.__tag         = tag

        self.__grad_active = grad_active
        
        #-- Necessary inputs for DLLM
        self.__n_sect              = None
        self.set_n_sect(n_sect)
        
        self.__distrib_type        = None
        self.set_distrib_type('cos_law')
#         self.__r_list_y        = None
        self.__r_list_eta          = None
        self.__ndv                 = None
        
        self.__perc_chord          = 0.25
        
        self.__Lref                = None
        self.__Lref_grad           = None
        
        self.__Sref                = None
        self.__Sref_grad           = None
        
        self.__AoA                 = None
        self.__AoA_grad            = None
        
        self.__twist               = None
        self.__twist_grad          = None
        
        self.__thetaY              = None
        
        self.__rel_thicks          = None
        self.__rel_thicks_grad     = None
        
        self.__rel_thicks_eta      = None
        self.__rel_thicks_grad_eta = None
        
        self.__chords              = None
        self.__chords_grad         = None

        self.__chords_eta          = None
        self.__chords_grad_eta     = None
        
        self.__sweep               = None
        self.__sweep_grad          = None
        
        self.__sweep_eta           = None
        self.__sweep_grad_eta      = None
        
        self.__XYZ                 = None
        self.__XYZ_grad            = None
        
        self.__eta                 = None
        self.__eta_grad            = None
        
        self.__airfoil_type        = None# for future use in case of meta airfoil with dynamic number of parameters
        self.__ref_airfoil         = None # reference airfoil, only used if the same airfoil is put for all sections
        self.__airfoils            = None # Airfoil list for each section
        self.__linked_airfoils     = None # Airfoil scaled to the the planform
        
    # -- Accessors
    def get_tag(self):
        return self.__tag
    
    def get_grad_active(self):
        return self.__grad_active
    
    def get_n_sect(self):
        return self.__n_sect
    
    def get_AoA(self):
        return self.__AoA
    
    def get_AoA_grad(self):
        return self.__AoA_grad
    
    def get_r_list_eta(self):
        return self.__r_list_eta
    
    def get_ndv(self):
        return self.__ndv
    
    def get_Lref(self):
        return self.__Lref
    
    def get_Lref_grad(self):
        return self.__Lref_grad
    
    def get_Sref(self):
        return self.__Sref
    
    def get_Sref_grad(self):
        return self.__Sref_grad
    
    def get_AR(self):
        return self.__AR

    def get_AR_grad(self):
        return self.__AR_grad

    #-- TBC: need to remove but not sure if it impacts some computations
    def get_fuel(self):
        return self.__fuel

    def get_fuel_grad(self):
        return self.__fuel_grad
    #-- end TBC
    
    def get_twist(self):
        return self.__twist
    
    def get_twist_grad(self):
        return self.__twist_grad
    
    def get_thetaY(self):
        return self.__thetaY
    
    def get_chords_eta(self):
        return self.__chords_eta
    
    def get_chords_grad_eta(self):
        return self.__chords_grad_eta
    
    def get_chords(self):
        return self.__chords
    
    def get_chords_grad(self):
        return self.__chords_grad
    
    def get_rel_thicks(self):
        return self.__rel_thicks
    
    def get_rel_thicks_grad(self):
        return self.__rel_thicks_grad
    
    def get_XYZ(self):
        return self.__XYZ
    
    def get_XYZ_grad(self):
        return self.__XYZ_grad
    
    def get_eta(self):
        return self.__eta
    
    def get_eta_grad(self):
        return self.__eta_grad
    
    def get_airfoil_type(self):
        return self.__airfoil_type
    
    def get_linked_airfoils(self):
        return self.__linked_airfoils
    
    #-- Setters
    def set_tag(self, tag):
        self.__tag = tag
        
    def set_grad_active(self, grad_active):
        self.__grad_active = grad_active
        
    def set_n_sect(self, n_sect):
        ERROR_MSG = self.ERROR_MSG + '.set_n_sect: ' + str(self.__tag) + ': '
        if n_sect % 2 != 0:
            raise Exception(
                ERROR_MSG +
                'The total number of elements in the wing must be even.')
        self.__n_sect = n_sect
        
    def set_distrib_type(self, distrib_type):
        if not distrib_type in self.POS_DISTRIB:
            raise Exception("distrib_type :" +
                            str(distrib_type) +
                            " not in possible types : " +
                            str(self.POS_DISTRIB))
        self.__distrib_type = distrib_type
        
    def set_ndv(self, ndv):
        self.__ndv = ndv
        
    def set_perc_chord(self, perc_chord):
        self.__perc_chord = perc_chord
        
    def set_AoA(self, AoA):
        self.__AoA = AoA
        
    def set_AoA_grad(self, AoA_grad):
        self.__AoA_grad = AoA_grad
    
    def set_twist(self, twist):
        self.__twist = twist
        
    def set_twist_grad(self, twist_grad):
        self.__twist_grad = twist_grad
        
    def set_thetaY(self, thetaY):
        self.__thetaY = thetaY
        
    def set_chords_eta(self, chords_eta):
        self.__chords_eta = chords_eta
        
    def set_chords_grad_eta(self, chords_grad_eta):
        self.__chords_grad_eta = chords_grad_eta
        
    def set_rel_thicks_eta(self, rel_thicks_eta):
        self.__rel_thicks_eta = rel_thicks_eta
        
    def set_rel_thicks_grad_eta(self, rel_thicks_grad_eta):
        self.__rel_thicks_grad_eta = rel_thicks_grad_eta
        
    def set_sweep_eta(self, sweep_eta):
        self.__sweep_eta = sweep_eta 
        
    def set_sweep_grad_eta(self, sweep_grad_eta):
        self.__sweep_grad_eta = sweep_grad_eta
        
    def set_eta(self, eta):
        self.__eta = eta 
        
    def set_eta_grad(self, eta_grad):
        self.__eta_grad = eta_grad
        
    def set_airfoil_type(self, airfoil_type):
        self.__airfoil_type = airfoil_type
        
    def set_ref_aifoil(self, ref_airfoil):
        self.__ref_airfoil = ref_airfoil

    def set_airfoils(self, airfoils):
        ERROR_MSG = self.ERROR_MSG + 'set_airfoils: '
        if len(airfoils) != self.__n_sect:
            print ERROR_MSG + 'number of airfoils is different than the number of sections. attribute not set'
            self.__airfoils = None
        else:
            self.__airfoils = airfoils

    #-- Methods     
    def build_r_lists(self, n_sect=None):
        if n_sect is not None:
            self.set_n_sect(n_sect)
            
        N = self.get_n_sect()

        if self.__distrib_type == 'linear':
            self.__r_list_eta = np.linspace(-0.5, 0.5, N + 1)
        else:
            theta_list = np.linspace(0., np.pi, N + 1)
            self.__r_list_eta = -0.5 + 0.5 * (1. - np.cos(theta_list))
            
    def build_airfoils_from_ref(self):
        ERROR_MSG = self.ERROR_MSG + 'build_airfoils_from_ref: '
        if self.__ref_airfoil is None:
            print ERROR_MSG + 'cannot build airfoils if reference airfoil is not defined.'
        else:
            airfoils = []
            for n in xrange(self.__n_sect):
                airfoils.append(self.__ref_airfoil.get_scaled_copy())

            self.set_airfoils(airfoils)
            
    def build_linear_airfoil(self, OC, AoA0=0., Sref=1.,
                             Lref=1., rel_thick=0., pcop=0.25, Ka=0.95, set_as_ref=True):
        self.__airfoil_type = 'simple'
        airfoil = AnalyticAirfoil(
            OC,
            AoA0=AoA0,
            Sref=Sref,
            Lref=Lref,
            rel_thick=rel_thick,
            pcop=pcop,
            Ka=Ka, grad_active=self.__grad_active)
        if set_as_ref:
            self.set_ref_aifoil(airfoil)
        return airfoil

    def build_meta_airfoil(self, OC, surrogate_model, relative_thickness=.12,
                           camber=0., Sref=1., Lref=1., sweep=.0, set_as_ref=True):
        self.__airfoil_type = 'meta'
        airfoil = MetaAirfoil(OC, surrogate_model, relative_thickness=relative_thickness,
                              camber=camber, Sref=Sref, Lref=Lref, sweep=sweep)
        if set_as_ref:
            self.set_ref_aifoil(airfoil)
        return airfoil
    
    def update(self):
        self.__check_thetaY()
        self.__build_data_from_eta()
        self.__check_airfoils_inputs()
        self.__link_airfoils_to_geom()
        self.__compute_Sref_Lref_AR_fuel() 
        
    def plot(self, prefix=None):
        N = self.get_n_sect()
        fact = self.__perc_chord
        if prefix is None:  
            name = self.get_tag()
        else:
            name = prefix+'_'+self.get_tag()
        
        #-- Planform mehses plot
        Y_list =   self.__eta[1,:]
        X_list = - self.__eta[0,:]
        
        Yp_list =   self.__XYZ[1,:]
        Xp_list = - self.__XYZ[0,:]
        
        lim = 1.1*max(max(Y_list),max(self.__chords_eta))
        
        # set limits to the graph
        plt.xlim(-lim, +lim)
        plt.ylim(-lim, +lim)
        
        # set graph labels (aero axes)
        plt.xlabel('y')
        plt.ylabel('-x')
        
        # plot eta points
        plt.plot(Y_list, X_list,'ro',markersize=3.0)
        
        # plot Y points
        plt.plot(Yp_list, Xp_list,'bo',markersize=3.0)

        # compute X for LE_points and TE_points
        X_LE_points = X_list[:] + fact*self.__chords_eta[:]
        X_TE_points = X_list[:] - (1.-fact)*self.__chords_eta[:]
        
        Xp_LE_points = Xp_list[:] + fact*self.__chords[:]
        Xp_TE_points = Xp_list[:] - (1-fact)*self.__chords[:]
        
        plt.plot(Y_list, X_LE_points,'-k')
        plt.plot(Y_list, X_TE_points,'-k')
        
        for i in xrange(N+1):
            plt.plot([Y_list[i],Y_list[i]],[X_LE_points[i],X_TE_points[i]],'-r')
            
        for i in xrange(N):
            plt.plot([Yp_list[i],Yp_list[i]],[Xp_LE_points[i],Xp_TE_points[i]],'--b')    
        

        plt.rc("font", size=14)
        plt.savefig(name+"_planform_meshes.png",format='png')
        plt.close()
        
        #-- planform plot
        # set limits to the graph
        plt.xlim(-lim, +lim)
        plt.ylim(-lim, +lim)
        
        # set graph labels (aero axes)
        plt.xlabel('y')
        plt.ylabel('-x')
        
        # compute X for LE_points and TE_points
        X_LE_points = X_list[:] + fact*self.__chords_eta[:]
        X_TE_points = X_list[:] - (1.-fact)*self.__chords_eta[:]
        
        plt.plot(Y_list, X_LE_points,'-k')
        plt.plot(Y_list, X_TE_points,'-k')
        plt.plot([Y_list[0],Y_list[0]],[X_LE_points[0],X_TE_points[0]],'-k')
        plt.plot([Y_list[-1],Y_list[-1]],[X_LE_points[-1],X_TE_points[-1]],'-k')     
 
        plt.rc("font", size=14)
        plt.savefig(name+"_planform.png",format='png')
        plt.close()
        
        
        #-- Chords plot
        plt.xlim(-lim, +lim)
        plt.xlabel('y')
        plt.ylabel('chords')
        plt.plot(Y_list,self.__chords_eta)
        plt.rc("font", size=14)
        plt.savefig(name+"_chords_distrib.png",format='png')
        plt.close()
        
        if self.__airfoil_type == 'Simple':
            #-- Toc plot
            plt.xlim(-lim, +lim)
            plt.xlabel('y')
            plt.ylabel('toc')
            plt.plot(Y_list,self.__rel_thicks_eta)
            plt.rc("font", size=14)
            plt.savefig(name+"_toc.png",format='png')
            plt.close()
            
            #-- Heights plot
            heights = self.__chords_eta*self.__rel_thicks_eta
            plt.xlim(-lim, +lim)
            plt.xlabel('y')
            plt.ylabel('heights')
            plt.plot(Y_list,heights)
            plt.rc("font", size=14)
            plt.savefig(name+"_heights_distrib.png",format='png')
            plt.close()
        
        #-- twist plot
        plt.xlim(-lim, +lim)
        plt.xlabel('y')
        plt.ylabel('twist (deg)')
        plt.plot(Yp_list,self.__twist*180./np.pi)
        plt.rc("font", size=14)
        plt.savefig(name+"_twist_distrib.png",format='png')
        plt.close()
        
    def __check_thetaY(self):
        # 6 dimensions for structural displacement:
        # dx,dy,dz,dthetax,dthetay,dthetaz at each section
        N = self.get_n_sect()
        if self.__thetaY is None:
            self.__thetaY = zeros(N)
            
    def __build_data_from_eta(self):
        N = self.get_n_sect()
        ndv = self.get_ndv()
        self.__XYZ      = np.zeros((3,N))
        self.__chords      = np.zeros(N)
        self.__rel_thicks      = np.zeros(N)
        self.__sweep      = np.zeros(N)
        
        self.__XYZ[:,:] = 0.5*(self.__eta[:,:-1]+self.__eta[:,1:])
        self.__chords[:] = 0.5*(self.__chords_eta[:-1]+self.__chords_eta[1:])
        self.__sweep[:]      = 0.5*(self.__sweep_eta[:-1]+self.__sweep_eta[1:])
        if self.__airfoil_type == 'Simple':
            self.__rel_thicks[:] = 0.5*(self.__rel_thicks_eta[:-1]+self.__rel_thicks_eta[1:])
        
        grad_active = self.get_grad_active()
        if grad_active:
            ndv = self.get_ndv()
            self.__XYZ_grad = np.zeros((3,N,ndv))
            self.__chords_grad = np.zeros((N,ndv))
            self.__rel_thicks_grad = np.zeros((N,ndv))
            self.__sweep_grad = np.zeros((N,ndv))
            self.__XYZ_grad[:,:,:] = 0.5*(self.__eta_grad[:,:-1,:]+self.__eta_grad[:,1:,:])
            self.__chords_grad[:,:] = 0.5*(self.__chords_grad_eta[:-1,:]+self.__chords_grad_eta[1:,:])
            self.__sweep_grad[:] = 0.5*(self.__sweep_grad_eta[:-1,:]+self.__sweep_grad_eta[1:,:])
            if self.__airfoil_type == 'Simple':
                self.__rel_thicks_grad[:,:] = 0.5*(self.__rel_thicks_grad_eta[:-1]+self.__rel_thicks_grad_eta[1:]) 

    def __link_airfoils_to_geom(self):
        ### ---- A modifier pour twist + modif DLLM Direct
        grad_active = self.get_grad_active()
        self.__linked_airfoils = []
        for i in range(self.__n_sect):
            LLoc, LLoc_grad, SLoc, SLoc_grad = self.__compute_local_info(i)
            linked_af = self.__airfoils[i].get_scaled_copy(Sref=SLoc, Lref=LLoc)
            linked_af.set_sweep(self.__sweep[i])
            linked_af.set_twist(self.__twist[i])
            if self.__airfoil_type == 'simple':
                linked_af.set_rel_thick(self.__rel_thicks[i])
            elif self.__airfoil_type == 'ref_CTA':
                linked_af.set_y_pos(self.__XYZ[1,i])
            if grad_active:
                linked_af.set_Sref_grad(SLoc_grad)
                linked_af.set_Lref_grad(LLoc_grad)
                linked_af.set_sweep_grad(self.__sweep_grad[i])
                linked_af.set_twist_grad(self.__twist_grad[i])
                if self.__airfoil_type == 'simple':
                    linked_af.set_rel_thick_grad(self.__rel_thicks_grad[i])

            self.__linked_airfoils.append(linked_af)
            
    #-- Private methods
    def __check_airfoils_inputs(self):
        ERROR_MSG = self.ERROR_MSG + \
            '__check_airfoils_inputs: ' + str(self.__tag) + ': '
        checked = True
        if self.__airfoils is None:
            checked = False
            print ERROR_MSG + 'airfoils attribute undefined. please set airfoils attribute.'
        elif len(self.__airfoils) != self.__n_sect:
            checked = False
            print ERROR_MSG + 'number of airfoils must equal to the number of geometrical sections.'

        if not checked:
            sys.exit(1)

    def __compute_local_info(self, i):
        LLoc = self.__chords[i]
        SLoc = LLoc * (self.__eta[1, i + 1] - self.__eta[1, i])
            
        grad_active = self.get_grad_active()
        if grad_active:
            LLoc_grad = self.__chords_grad[i]
            SLoc_grad = LLoc_grad * (self.__eta[1, i + 1] - self.__eta[1, i]) + \
                LLoc * (self.__eta_grad[1, i + 1, :] - self.__eta_grad[1, i, :])
        else:
            LLoc_grad = None
            SLoc_grad = None
            
        return LLoc, LLoc_grad, SLoc, SLoc_grad
    
    def __compute_Sref_Lref_AR_fuel(self):
        """
        Compute Lref and Sref from airfoils information, ToC and AR
        """
        N = self.get_n_sect()

        span      = self.__eta[1,-1]-self.__eta[1,0]

        self.__Sref = 0.
        self.__Lref = 0.
        self.__fuel = 0.

        for i, af in enumerate(self.__linked_airfoils):
            self.__Sref += af.get_Sref()
            self.__Lref += af.get_Lref()
            self.__fuel += self.__rel_thicks[i] * \
                self.__chords[i] * af.get_Sref() * 0.5

        self.__Lref /= N
        self.__AR = span**2 / self.__Sref
                 
        grad_active = self.get_grad_active()
        if grad_active:
            ndv=self.get_ndv()
            span_grad = self.__eta_grad[1,-1,:]-self.__eta_grad[1,0,:]
            self.__Sref_grad = zeros(ndv)
            self.__Lref_grad = zeros(ndv)
            self.__fuel_grad = zeros(ndv)
            for i, af in enumerate(self.__linked_airfoils):
                self.__Sref_grad += af.get_Sref_grad()
                self.__Lref_grad += af.get_Lref_grad()
                self.__fuel_grad += self.__rel_thicks_grad[i] * self.__chords[i] * af.get_Sref() * 0.5\
                    + self.__rel_thicks[i] * self.__chords_grad[i] * af.get_Sref() * 0.5\
                    + self.__rel_thicks[i] * self.__chords[i] * af.get_Sref_grad() * 0.5
                self.__Lref_grad /= N
                self.__AR_grad = (2. * span * span_grad *
                          self.__Sref - span**2 * self.__Sref_grad) / self.__Sref**2
            
    def __repr__(self):
        info_string = '\n*** Wing Geom Information ***'
        info_string += '\n  n_sect       : ' + str(self.get_n_sect())
        info_string += '\n  ndv          : ' + str(self.get_ndv())
        info_string += '\n  airfoil_type : ' + str(self.get_airfoil_type())
        return info_string
    
