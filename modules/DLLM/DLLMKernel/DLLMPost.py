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
#  http://github.com/TBD
#
import string
from numpy import zeros, dot
from numpy import sin, cos


class DLLMPost:

    DEF_F_LIST_NAMES = [
        'Lift',
        'Drag',
        'Drag_Pressure',
        'Drag_Induced',
        'Drag_Wave',
        'Drag_Friction',
        'Cl',
        'Cd',
        'Cdp',
        'Cdi',
        'Cdw',
        'Cdf',
        'LoD',
        'Sref']
    ERROR_MSG = 'ERROR in DLLMPost.'

    def __init__(self, LLW, verbose = 1):
        """
        Post-Processing module compatible with adjoint for DLLM
        """
        self.__LLW = LLW
        self.__F_list_names_calc = self.DEF_F_LIST_NAMES
        self.__F_list_calc = None
        self.__F_list_names = self.DEF_F_LIST_NAMES
        self.__F_list_dim = 0
        self.__F_list = None
        self.__dpF_list_dpiAoA = None
        self.__dpF_list_dpchi = None
        self.__dpF_list_dpAoA = None
        self.__computed = False

        self.__Lift_distrib_res = None
        self.__dpLift_distrib_res_dpiAoA = None
        self.__dpLift_distrib_res_dpAoA = None

        self.__target_loads = None
        
        self.__verbose = verbose
    #-- computed related methods
    def is_computed(self):
        return self.__computed

    def set_computed(self, bool=True):
        self.__computed = bool

    #-- Accessors
    def get_tag(self):
        return self.__LLW.get_tag()

    def get_N(self):
        return self.get_wing_param().get_n_sect()

    def get_ndv(self):
        return self.get_wing_param().get_ndv()

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

    def get_wing_param(self):
        return self.__LLW.get_wing_param()

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

    def set_target_loads_file(self, loads_file):
        N = self.get_N()

        loads = zeros(N)

        fid = open(loads_file, 'r')
        all_lines = fid.readlines()
        fid.close()
        i = 0
        for line in all_lines:
            words = string.split(line)
            if len(words) == 2:
                loads[i] = eval(words[1])
                i += 1

        self.set_target_loads(loads)

    def set_target_loads(self, loads):
        self.__target_loads = loads

        if 'target_loads' not in self.__F_list_names_calc:
            self.__F_list_names_calc.append('target_loads')

        if 'target_loads' not in self.__F_list_names:
            self.__F_list_names.append('target_loads')

    #-- Run method
    def __init_run(self):
        N = self.get_N()
        ndv = self.get_ndv()
        self.__F_list_calc_dim = len(self.__F_list_names_calc)
        self.__F_list_calc = zeros(self.__F_list_calc_dim)
        self.__F_list = zeros(self.__F_list_dim)
        self.__F_list_dim = len(self.__F_list_names)
        self.__F_list = zeros(self.__F_list_dim)
        self.__dpF_list_dpAoA = zeros(self.__F_list_dim)
        self.__dpF_list_dpiAoA = zeros((self.__F_list_dim, N))
        self.__dpF_list_dpchi = zeros((self.__F_list_dim, ndv))

    def run(self, F_list_names=None):
        ERROR_MSG = self.ERROR_MSG + 'run: '
        N = self.get_N()
        self.__init_run()

        # Default analysis
        for i, F_name in enumerate(self.__F_list_names_calc):
            if F_name == 'Cl':
                val = self.__Cl()
            elif F_name == 'Cd':
                val = self.__Cd()
            elif F_name == 'Cdp':
                val = self.__Cdp()
            elif F_name == 'Cdi':
                val = self.__Cdi()
            elif F_name == 'Cdw':
                val = self.__Cdw()
            elif F_name == 'Cdf':
                val = self.__Cdf()
#            Moments calculation are bugged for now
#             elif func == 'Cm':
#                 val = self.__Cm()
            elif F_name == 'Lift':
                val = self.__Lift()
            elif F_name == 'Drag':
                val = self.__Drag()
            elif F_name == 'Drag_Pressure':
                val = self.__Drag_Pressure()
            elif F_name == 'Drag_Induced':
                val = self.__Drag_Induced()
            elif F_name == 'Drag_Wave':
                val = self.__Drag_Wave()
            elif F_name == 'Drag_Friction':
                val = self.__Drag_Friction()
            elif F_name == 'LoD':
                val = self.__LoD()
            elif F_name == 'Sref':
                val = self.get_Sref()
            else:
                raise Exception(ERROR_MSG + ' unknown function ' + str(F_name))
            self.__F_list_calc[i] = val

        # Adjoint analysis
        if self.__verbose > 0 : print "Post : adjoint analysis"
        for i, F_name in enumerate(self.__F_list_names):
            if F_name == 'Cl':
                val = self.__Cl(distrib=True)
                dpFdpiAoA = self.__dpCl_dpiAoA()
                dpFdpchi = self.dpCl_dpchi()
                dpFdpAoA = self.dpCl_dpAoA()
            elif F_name == 'Cd':
                val = self.__Cd()
                dpFdpiAoA = self.__dpCd_dpiAoA()
                dpFdpchi = self.dpCd_dpchi()
                dpFdpAoA = self.dpCd_dpAoA()
            elif F_name == 'Cdi':
                val = self.__Cdi(distrib=True)
                dpFdpiAoA = self.__dpCdi_dpiAoA()
                dpFdpchi = self.dpCdi_dpchi()
                dpFdpAoA = self.dpCdi_dpAoA()
            elif F_name == 'Cdw':
                val = self.__Cdw()
                dpFdpiAoA = self.__dpCdw_dpiAoA()
                dpFdpchi = self.dpCdw_dpchi()
                dpFdpAoA = self.dpCdw_dpAoA()
            elif F_name == 'Cdp':
                val = self.__Cdp()
                dpFdpiAoA = self.__dpCdp_dpiAoA()
                dpFdpchi = self.dpCdp_dpchi()
                dpFdpAoA = self.dpCdp_dpAoA()
            elif F_name == 'Cdf':
                val = self.__Cdf()
                dpFdpiAoA = self.__dpCdf_dpiAoA()
                dpFdpchi = self.dpCdf_dpchi()
                dpFdpAoA = self.dpCdf_dpAoA()
#            Moments calculation are bugged for now
#             elif func == 'Cm':
#                 val = self.__Cm()
            elif F_name == 'Lift':
                val = self.__Lift()
                dpFdpiAoA = self.__dpLift_dpiAoA()
                dpFdpchi = self.dpLift_dpchi()
                dpFdpAoA = self.dpLift_dpAoA()
            elif F_name == 'Drag':
                val = self.__Drag()
                dpFdpiAoA = self.__dpDrag_dpiAoA()
                dpFdpchi = self.dpDrag_dpchi()
                dpFdpAoA = self.dpDrag_dpAoA()
            elif F_name == 'Drag_Pressure':
                val = self.__Drag_Pressure()
                dpFdpiAoA = self.__dpDrag_Pressure_dpiAoA()
                dpFdpchi = self.dpDrag_Pressure_dpchi()
                dpFdpAoA = self.dpDrag_Pressure_dpAoA()
            elif F_name == 'Drag_Induced':
                val = self.__Drag_Induced()
                dpFdpiAoA = self.__dpDrag_Induced_dpiAoA()
                dpFdpchi = self.dpDrag_Induced_dpchi()
                dpFdpAoA = self.dpDrag_Induced_dpAoA()
            elif F_name == 'Drag_Wave':
                val = self.__Drag_Wave()
                dpFdpiAoA = self.__dpDrag_Wave_dpiAoA()
                dpFdpchi = self.dpDrag_Wave_dpchi()
                dpFdpAoA = self.dpDrag_Wave_dpAoA()
            elif F_name == 'Drag_Friction':
                val = self.__Drag_Friction()
                dpFdpiAoA = self.__dpDrag_Friction_dpiAoA()
                dpFdpchi = self.dpDrag_Friction_dpchi()
                dpFdpAoA = self.dpDrag_Friction_dpAoA()
            elif F_name == 'LoD':
                val = self.__LoD()
                dpFdpiAoA = self.__dpLoD_dpiAoA()
                dpFdpchi = self.dpLoD_dpchi()
                dpFdpAoA = self.dpLoD_dpAoA()
            elif F_name == 'Sref':
                val = self.get_Sref()
                dpFdpiAoA = zeros(N)
                dpFdpchi = self.get_Sref_grad()
                dpFdpAoA = 0.
            else:
                raise Exception(ERROR_MSG + ' unknown function ' + str(F_name))
            self.__F_list[i] = val
            self.__dpF_list_dpiAoA[i, :] = dpFdpiAoA[:]
            self.__dpF_list_dpchi[i, :] = dpFdpchi[:]
            self.__dpF_list_dpAoA[i] = dpFdpAoA

        self.__Lift_distrib_res = self.__Lift_distrib(distrib=True)

        if self.__verbose > 0 :
            self.__display_info()
        self.set_computed(True)

    #-- Computation methods
    #-- Lift distrib related methods
    def comp_Lift_distrib(self, distrib=False):
        self.__Lift_distrib(distrib=distrib)
        return self.__Lift_distrib_res

    def comp_dpLift_distrib_dpiAoA(self):
        self.__dpLift_distrib_dpiAoA()
        return self.__dpLift_distrib_res_dpiAoA

    def comp_dpLift_distrib_dpAoA(self):
        self.__dpLift_distrib_dpAoA()
        return self.__dpLift_distrib_res_dpAoA

    def __Lift_distrib(self, distrib=False):
        N = self.get_N()
        self.__Lift_distrib_res = zeros(N)

        Pdyn = self.get_OC().get_Pdyn()
        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        airfoils = self.get_airfoils()

        Check_Lift = 0.

        for i in xrange(N):
            af = airfoils[i]
            self.__Lift_distrib_res[i] = af.Cl(
                localAoA[i], Mach) * cos(iAoA[i]) * af.get_Sref() * Pdyn
            Check_Lift += self.__Lift_distrib_res[i]

        # print 'CHECK LIFT = ',Check_Lift

        if distrib:
            y = self.get_wing_param().get_XYZ()[1, :]
            fid = open(self.get_tag() + '_Lift_distrib.dat', 'w')
            line = "#Slice\t%24s" % "y" + "\t%24s" % "Lift" + "\n"
            fid.write(line)
            for i in xrange(N):
                line = str(
                    i) + "\t%24.16e" % y[i] + "\t%24.16e" % self.__Lift_distrib_res[i] + "\n"
                fid.write(line)
            fid.close()

    def __dpLift_distrib_dpiAoA(self):
        N = self.get_N()
        self.__dpLift_distrib_res_dpiAoA = zeros((N, N))
        Pdyn = self.get_OC().get_Pdyn()

        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        dplocalAoA_dpiAoA = self.get_dplocalAoA_dpiAoA()
        airfoils = self.get_airfoils()
        for i in range(N):
            af = airfoils[i]
            self.__dpLift_distrib_res_dpiAoA[
                i, i] = -af.Cl(localAoA[i], Mach) * sin(iAoA[i]) * af.get_Sref() * Pdyn
            self.__dpLift_distrib_res_dpiAoA[i, :] = af.ClAlpha(
                localAoA[i], Mach=Mach) * cos(iAoA[i]) * dplocalAoA_dpiAoA[i] * af.get_Sref() * Pdyn

    def __dpLift_distrib_dpAoA(self):
        N = self.get_N()
        self.__dpLift_distrib_res_dpAoA = zeros(N)
        Pdyn = self.get_OC().get_Pdyn()

        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        dplocalAoA_dpAoA = self.get_dplocalAoA_dpAoA()
        airfoils = self.get_airfoils()
        for i in range(N):
            af = airfoils[i]
            self.__dpLift_distrib_res_dpAoA[i] = af.ClAlpha(
                localAoA[i], Mach=Mach) * cos(iAoA[i]) * af.get_Sref() * dplocalAoA_dpAoA[i] * Pdyn

    #-- Cl related methods
    def comp_Cl(self):
        Cl = self.__Cl()
        return Cl

    def comp_dpCl_dpiAoA(self):
        dpCl_dpiAoA = self.__dpCl_dpiAoA()
        return dpCl_dpiAoA

    def __Cl(self, distrib=False):
        N = self.get_N()

        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        airfoils = self.get_airfoils()
        Cl = 0.0
        Cl_loc = zeros(N)

        for i in xrange(N):
            af = airfoils[i]
            # change frame from local profile to wing profile and update
            # reference scale from profile to wing
            Cl_loc[i] = af.Cl(localAoA[i], Mach) * \
                cos(iAoA[i]) * af.get_Sref() / self.get_Sref()
            Cl += Cl_loc[i]

        if distrib:
            y = self.get_wing_param().get_XYZ()[1, :]
            fid = open(self.get_tag() + '_Cl_distrib.dat', 'w')
            line = "#Slice\t%24s" % "y" + "\t%24s" % "Cl" + "\n"
            fid.write(line)
            for i in xrange(N):
                line = str(
                    i) + "\t%24.16e" % y[i] + "\t%24.16e" % Cl_loc[i] + "\n"
                fid.write(line)
            fid.close()

        return Cl

    def __dpCl_dpiAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        dplocalAoA_dpiAoA = self.get_dplocalAoA_dpiAoA()
        airfoils = self.get_airfoils()
        dCl_diAoA = zeros(N)
        for i in range(N):
            af = airfoils[i]
            dCl_diAoA[i] -= af.Cl(localAoA[i], Mach) * \
                sin(iAoA[i]) * af.get_Sref()
            dCl_diAoA += af.ClAlpha(localAoA[i], Mach=Mach) * \
                cos(iAoA[i]) * dplocalAoA_dpiAoA[i] * af.get_Sref()

        dCl_diAoA /= self.get_Sref()

        return dCl_diAoA

    def dpCl_dpAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        dplocalAoA_dpAoA = self.get_dplocalAoA_dpAoA()
        airfoils = self.get_airfoils()
        dCl_dAoA = 0.
        for i in range(N):
            af = airfoils[i]
            dCl_dAoA += af.ClAlpha(localAoA[i], Mach=Mach) * \
                cos(iAoA[i]) * af.get_Sref() * dplocalAoA_dpAoA[i]
        dCl_dAoA /= self.get_Sref()

        return dCl_dAoA

    def dpCl_dpchi(self):
        N = self.get_N()
        ndv = self.get_ndv()
        Mach = self.get_OC().get_Mach()
        iAoA = self.get_iAoA()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        dlAoAdchi = self.get_dplocalAoA_dpchi()

        Cl = 0.0
        dCl = zeros(ndv)
        for i in xrange(N):
            af = airfoils[i]
            # change frame from local profile to wing profile and update
            # reference scale from profile to wing
            Cl += af.Cl(localAoA[i], Mach) * cos(iAoA[i]) * af.get_Sref()
            dCl += cos(iAoA[i]) * (af.ClAlpha(localAoA[i], Mach) * af.get_Sref() * dlAoAdchi[i, :] +
                                   # This term comes from local profile design
                                   # variables influence
                                   af.dCl_dchi(localAoA[i], Mach) * af.get_Sref() +
                                   af.Cl(localAoA[i], Mach) * af.get_Sref_grad())

        dCldchi = (
            dCl * self.get_Sref() - Cl * self.get_Sref_grad()) / (self.get_Sref()**2)
        return dCldchi

    #-- Cd related methods
    def __Cd(self):
        Cdp = self.__Cdp()
        Cdf = self.__Cdf()
        Cd = Cdp + Cdf
        return Cd

    def __dpCd_dpiAoA(self):
        dpCdp_dpiAoA = self.__dpCdp_dpiAoA()
        dpCdf_dpiAoA = self.__dpCdf_dpiAoA()
        dpCd_dpiAoA = dpCdp_dpiAoA + dpCdf_dpiAoA
        return dpCd_dpiAoA

    def dpCd_dpAoA(self):
        dpCdp_dpAoA = self.dpCdp_dpAoA()
        dpCdf_dpAoA = self.dpCdf_dpAoA()
        dpCd_dpAoA = dpCdp_dpAoA + dpCdf_dpAoA
        return dpCd_dpAoA

    def dpCd_dpchi(self):
        dpCdp_dpchi = self.dpCdp_dpchi()
        dpCdf_dpchi = self.dpCdf_dpchi()
        dpCd_dpchi = dpCdp_dpchi + dpCdf_dpchi
        return dpCd_dpchi

    #-- Cdp related methods
    def __Cdp(self):
        Cdi = self.__Cdi()
        Cdw = self.__Cdw()
        Cdp = Cdi + Cdw
        return Cdp

    def __dpCdp_dpiAoA(self):
        dpCdi_dpiAoA = self.__dpCdi_dpiAoA()
        dpCdw_dpiAoA = self.__dpCdw_dpiAoA()
        dpCdp_dpiAoA = dpCdi_dpiAoA + dpCdw_dpiAoA
        return dpCdp_dpiAoA

    def dpCdp_dpAoA(self):
        dpCdi_dpAoA = self.dpCdi_dpAoA()
        dpCdw_dpAoA = self.dpCdw_dpAoA()
        dpCdp_dpAoA = dpCdi_dpAoA + dpCdw_dpAoA
        return dpCdp_dpAoA

    def dpCdp_dpchi(self):
        dpCdi_dpchi = self.dpCdi_dpchi()
        dpCdw_dpchi = self.dpCdw_dpchi()
        dpCdp_dpchi = dpCdi_dpchi + dpCdw_dpchi
        return dpCdp_dpchi

    #-- Cdi related methods
    def __Cdi(self, distrib=False):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        Cdi_loc = zeros(N)
        Cdi = 0.0
        for i in xrange(N):
            af = airfoils[i]
            AoA = localAoA[i]
            Cdi_loc[i] = -af.Cl(AoA, Mach) * \
                sin(iAoA[i]) * af.get_Sref() / self.get_Sref()
            Cdi += Cdi_loc[i]

        if distrib:
            y = self.get_wing_param().get_XYZ()[1, :]
            fid = open(self.get_tag() + '_Cdi_distrib.dat', 'w')
            line = "#Slice\t%24s" % "y" + "\t%24s" % "Cdi" + "\n"
            fid.write(line)
            for i in xrange(N):
                line = str(
                    i) + "\t%24.16e" % y[i] + "\t%24.16e" % Cdi_loc[i] + "\n"
                fid.write(line)
            fid.close()
        return Cdi

    def __dpCdi_dpiAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        dplocalAoA_dpiAoA = self.get_dplocalAoA_dpiAoA()

        dCd_dAoA = zeros(N)
        for i in range(N):  # Dependance by induced angles
            af = airfoils[i]
            AoA = localAoA[i]
            dCd_dAoA[i] = -af.ClAlpha(AoA, Mach) * sin(iAoA[i]) * af.get_Sref()

        dCd_diAoA = dot(dCd_dAoA, dplocalAoA_dpiAoA)

        # Dependance by projection of Cl : the induced drag.
        for i in range(N):
            af = airfoils[i]
            dCd_diAoA[i] -= af.Cl(localAoA[i], Mach) * \
                af.get_Sref() * cos(iAoA[i])
        dCd_diAoA /= self.get_Sref()

        return dCd_diAoA

    def dpCdi_dpAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        dplocalAoA_dpAoA = self.get_dplocalAoA_dpAoA()

        dCd_dAoA = 0.
        for i in range(N):  # Dependance by induced angles
            af = airfoils[i]
            lAoA = localAoA[i]
            dCd_dAoA -= af.ClAlpha(lAoA, Mach) * \
                sin(iAoA[i]) * af.get_Sref() * dplocalAoA_dpAoA[i]
        dCd_dAoA /= self.get_Sref()

        return dCd_dAoA

    def dpCdi_dpchi(self):
        N = self.get_N()
        ndv = self.get_ndv()
        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        airfoils = self.get_airfoils()
        dlAoAdchi = self.get_dplocalAoA_dpchi()

        Cd = 0.0
        dCd = zeros(ndv)
        for i in range(N):
            af = airfoils[i]
            AoA = localAoA[i]
            Cdloc = -af.Cl(AoA, Mach) * sin(iAoA[i])
            dCdloc = -af.ClAlpha(AoA,
                                 Mach) * sin(iAoA[i]) * dlAoAdchi[i] - af.dCl_dchi(AoA,
                                                                                   Mach) * sin(iAoA[i])
            Cd += Cdloc * af.get_Sref()
            dCd += dCdloc * af.get_Sref() + Cdloc * af.get_Sref_grad()
        dCddchi = (
            dCd * self.get_Sref() - Cd * self.get_Sref_grad()) / (self.get_Sref()**2)

        return dCddchi

    #-- Cdw related methods
    def __Cdw(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        Cd = 0.0
        for i in xrange(N):
            af = airfoils[i]
            AoA = localAoA[i]
            Cdloc = af.Cdp(AoA, Mach)
            Cd += Cdloc * af.get_Sref()
        Cd /= self.get_Sref()

        return Cd

    def __dpCdw_dpiAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        dplocalAoA_dpiAoA = self.get_dplocalAoA_dpiAoA()

        dCd_dAoA = zeros(N)
        for i in range(N):  # Dependance by induced angles
            af = airfoils[i]
            AoA = localAoA[i]
            dCd_dAoA[i] = af.CdpAlpha(AoA, Mach) * af.get_Sref()

        dCd_diAoA = dot(dCd_dAoA, dplocalAoA_dpiAoA)
        dCd_diAoA /= self.get_Sref()

        return dCd_diAoA

    def dpCdw_dpAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        dplocalAoA_dpAoA = self.get_dplocalAoA_dpAoA()

        dCd_dAoA = 0.
        for i in range(N):  # Dependance by induced angles
            af = airfoils[i]
            lAoA = localAoA[i]
            dCd_dAoA += af.CdpAlpha(lAoA, Mach) * \
                af.get_Sref() * dplocalAoA_dpAoA[i]
        dCd_dAoA /= self.get_Sref()

        return dCd_dAoA

    def dpCdw_dpchi(self):
        N = self.get_N()
        ndv = self.get_ndv()
        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        airfoils = self.get_airfoils()
        dlAoAdchi = self.get_dplocalAoA_dpchi()

        Cd = 0.0
        dCd = zeros(ndv)
        for i in range(N):
            af = airfoils[i]
            AoA = localAoA[i]
            Cdloc = af.Cdp(AoA, Mach)
            dCdloc = af.CdpAlpha(
                AoA, Mach) * dlAoAdchi[i] + af.dCdp_dchi(AoA, Mach)
            Cd += Cdloc * af.get_Sref()
            dCd += dCdloc * af.get_Sref() + Cdloc * af.get_Sref_grad()
        dCddchi = (
            dCd * self.get_Sref() - Cd * self.get_Sref_grad()) / (self.get_Sref()**2)

        return dCddchi

    #-- Cdf related methods
    def __Cdf(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        Cd = 0.0
        for i in xrange(N):
            af = airfoils[i]
            AoA = localAoA[i]
            Cdloc = af.Cdf(AoA, Mach)
            Cd += Cdloc * af.get_Sref()
        Cd /= self.get_Sref()

        return Cd

    def __dpCdf_dpiAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        dplocalAoA_dpiAoA = self.get_dplocalAoA_dpiAoA()

        dCd_dAoA = zeros(N)
        for i in range(N):  # Dependance by induced angles
            af = airfoils[i]
            AoA = localAoA[i]
            dCd_dAoA[i] = af.CdfAlpha(AoA, Mach) * af.get_Sref()

        dCd_diAoA = dot(dCd_dAoA, dplocalAoA_dpiAoA)
        dCd_diAoA /= self.get_Sref()

        return dCd_diAoA

    def dpCdf_dpAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        dplocalAoA_dpAoA = self.get_dplocalAoA_dpAoA()

        dCd_dAoA = 0.
        for i in range(N):  # Dependance by induced angles
            af = airfoils[i]
            lAoA = localAoA[i]
            dCd_dAoA += af.CdfAlpha(lAoA, Mach) * \
                af.get_Sref() * dplocalAoA_dpAoA[i]
        dCd_dAoA /= self.get_Sref()

        return dCd_dAoA

    def dpCdf_dpchi(self):
        N = self.get_N()
        ndv = self.get_ndv()
        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA = self.get_iAoA()
        airfoils = self.get_airfoils()
        dlAoAdchi = self.get_dplocalAoA_dpchi()

        Cd = 0.0
        dCd = zeros(ndv)
        mySref = 0.
        for i in range(N):
            af = airfoils[i]
            lAoA = localAoA[i]
            Cdloc = af.Cdf(lAoA, Mach)
            dCdloc = af.CdfAlpha(lAoA, Mach) * dlAoAdchi[i] +\
                     af.dCdf_dchi(lAoA, Mach)
            Cd += Cdloc * af.get_Sref()
            dCd += dCdloc * af.get_Sref() + Cdloc * af.get_Sref_grad()
        dCddchi = (dCd * self.get_Sref() - Cd * self.get_Sref_grad()) / \
                  (self.get_Sref()**2)
        return dCddchi

    #-- Cm related methods (bugged for now)
    def __Cm(self):
        # Issue in moments calculation, no distances or sweep taken into
        # account...
        N = self.get_N()
        print 'WARNING: Cm calculation is not valid !!!...'
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        Cm = 0.0
        for i in xrange(N):
            af = airfoils[i]
            Cm += af.Cm(localAoA[i], Mach) * af.get_Sref()
        Cm /= self.get_Sref()
        return Cm

    def __dpCm_dpiAoA(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        dplocalAoA_dpiAoA = self.get_dplocalAoA_dpiAoA()
        dCm = 0.0
        for i in range(N):
            af = airfoils[i]
            dCm += af.CmAlpha(localAoA[i], Mach) * \
                af.getSref() * dplocalAoA_dpiAoA[i]
        dCm /= self.get_Sref()
        return dCm

    def dpCm_dpchi(self):
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        dlAoAdchi = self.get_dplocalAoA_dpchi()

        Cm = 0.0
        dCm = 0.0
        for i in range(N):
            af = airfoils[i]
            Cm += af.Cm(localAoA[i], Mach) * af.get_Sref()
            dCm += af.CmAlpha(localAoA[i],
                              Mach) * af.get_Sref() * dlAoAdchi[i] + af.dCm_dchi(localAoA[i],
                                                                                 Mach) * af.get_Sref() + af.Cm(localAoA[i],
                                                                                                               Mach) * af.get_Sref_grad()
        Cm /= self.get_Sref()
        dCmdchi = (
            dCm * self.get_Sref() - Cm * self.get_Sref_grad()) / (self.get_Sref()**2)

        return dCmdchi

    #-- Lift related methods
    def comp_Lift(self):
        Lift = self.__Lift()
        return Lift

    def get_Lift(self):
        return self.__Lift()

    def comp_dpLift_dpiAoA(self):
        dpLift_dpiAoA = self.__dpLift_dpiAoA()
        return dpLift_dpiAoA

    def __Lift(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Cl = self.__Cl()
        Lift = Pdyn * Sref * Cl
        return Lift

    def dpLift_dpAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCl_dAoA = self.dpCl_dpAoA()
        dLift_dAoA = Pdyn * Sref * dCl_dAoA
        return dLift_dAoA

    def __dpLift_dpiAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCl_diAoA = self.__dpCl_dpiAoA()
        dLift_diAoA = Pdyn * Sref * dCl_diAoA
        return dLift_diAoA

    def dpLift_dpchi(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Sref_grad = self.get_Sref_grad()
        Cl = self.__Cl()
        dpCl_dpchi = self.dpCl_dpchi()
        dpLift_dpchi = Pdyn * Sref * dpCl_dpchi + Pdyn * Sref_grad * Cl
        return dpLift_dpchi

    #-- Drag related methods
    def __Drag(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Cd = self.__Cd()
        Drag = Pdyn * Sref * Cd
        return Drag
    
    def get_Drag(self):
        return self.__Drag()

    def get_Drag_Friction(self):
        return self.__Drag_Friction()

    def get_Drag_Induced(self):
        return self.__Drag_Induced()

    def get_Drag_Wave(self):
        return self.__Drag_Wave()
    
    def get_Drag_Pressure(self):
        return self.__Drag_Pressure()
    
    def dpDrag_dpAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_dAoA = self.dpCd_dpAoA()
        dDrag_dAoA = Pdyn * Sref * dCd_dAoA
        return dDrag_dAoA

    def __dpDrag_dpiAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_diAoA = self.__dpCd_dpiAoA()
        dDrag_diAoA = Pdyn * Sref * dCd_diAoA
        return dDrag_diAoA

    def dpDrag_dpchi(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Sref_grad = self.get_Sref_grad()
        Cd = self.__Cd()
        dpCd_dpchi = self.dpCd_dpchi()
        dpDrag_dpchi = Pdyn * Sref * dpCd_dpchi + Pdyn * Sref_grad * Cd
        return dpDrag_dpchi

    #-- Pressure Drag related methods
    def __Drag_Pressure(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Cd = self.__Cdp()
        Drag = Pdyn * Sref * Cd
        return Drag

    def dpDrag_Pressure_dpAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_dAoA = self.dpCdp_dpAoA()
        dDrag_dAoA = Pdyn * Sref * dCd_dAoA
        return dDrag_dAoA

    def __dpDrag_Pressure_dpiAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_diAoA = self.__dpCdp_dpiAoA()
        dDrag_diAoA = Pdyn * Sref * dCd_diAoA
        return dDrag_diAoA

    def dpDrag_Pressure_dpchi(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Sref_grad = self.get_Sref_grad()
        Cd = self.__Cdp()
        dpCd_dpchi = self.dpCdp_dpchi()
        dpDrag_dpchi = Pdyn * Sref * dpCd_dpchi + Pdyn * Sref_grad * Cd
        return dpDrag_dpchi

    #-- Friction Drag related methods
    def __Drag_Friction(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Cd = self.__Cdf()
        Drag = Pdyn * Sref * Cd
        return Drag

    def dpDrag_Friction_dpAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_dAoA = self.dpCdf_dpAoA()
        dDrag_dAoA = Pdyn * Sref * dCd_dAoA
        return dDrag_dAoA

    def __dpDrag_Friction_dpiAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_diAoA = self.__dpCdf_dpiAoA()
        dDrag_diAoA = Pdyn * Sref * dCd_diAoA
        return dDrag_diAoA

    def dpDrag_Friction_dpchi(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Sref_grad = self.get_Sref_grad()
        Cd = self.__Cdf()
        dpCd_dpchi = self.dpCdf_dpchi()
        dpDrag_dpchi = Pdyn * Sref * dpCd_dpchi + Pdyn * Sref_grad * Cd
        return dpDrag_dpchi

    #-- Induced Drag related methods
    def __Drag_Induced(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Cd = self.__Cdi()
        Drag = Pdyn * Sref * Cd
        return Drag

    def dpDrag_Induced_dpAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_dAoA = self.dpCdi_dpAoA()
        dDrag_dAoA = Pdyn * Sref * dCd_dAoA
        return dDrag_dAoA

    def __dpDrag_Induced_dpiAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_diAoA = self.__dpCdi_dpiAoA()
        dDrag_diAoA = Pdyn * Sref * dCd_diAoA
        return dDrag_diAoA

    def dpDrag_Induced_dpchi(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Sref_grad = self.get_Sref_grad()
        Cd = self.__Cdi()
        dpCd_dpchi = self.dpCdi_dpchi()
        dpDrag_dpchi = Pdyn * Sref * dpCd_dpchi + Pdyn * Sref_grad * Cd
        return dpDrag_dpchi

    #-- Wave Drag related methods
    def __Drag_Wave(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Cd = self.__Cdw()
        Drag = Pdyn * Sref * Cd
        return Drag

    def dpDrag_Wave_dpAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_dAoA = self.dpCdw_dpAoA()
        dDrag_dAoA = Pdyn * Sref * dCd_dAoA
        return dDrag_dAoA

    def __dpDrag_Wave_dpiAoA(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        dCd_diAoA = self.__dpCdw_dpiAoA()
        dDrag_diAoA = Pdyn * Sref * dCd_diAoA
        return dDrag_diAoA

    def dpDrag_Wave_dpchi(self):
        Pdyn = self.get_OC().get_Pdyn()
        Sref = self.get_Sref()
        Sref_grad = self.get_Sref_grad()
        Cd = self.__Cdw()
        dpCd_dpchi = self.dpCdw_dpchi()
        dpDrag_dpchi = Pdyn * Sref * dpCd_dpchi + Pdyn * Sref_grad * Cd
        return dpDrag_dpchi

    #-- LoD related methods
    def __LoD(self):
        Cl = self.__Cl()
        Cd = self.__Cd()
        LoD = Cl / Cd
        return LoD

    def dpLoD_dpAoA(self):
        Cl = self.__Cl()
        Cd = self.__Cd()
        dpCl_dpAoA = self.dpCl_dpAoA()
        dpCd_dpAoA = self.dpCd_dpAoA()
        dpLoD_dpAoA = (dpCl_dpAoA * Cd - Cl * dpCd_dpAoA) / Cd**2
        return dpLoD_dpAoA

    def __dpLoD_dpiAoA(self):
        Cl = self.__Cl()
        Cd = self.__Cd()
        dpCl_dpiAoA = self.__dpCl_dpiAoA()
        dpCd_dpiAoA = self.__dpCd_dpiAoA()
        dpLoD_dpiAoA = (dpCl_dpiAoA * Cd - Cl * dpCd_dpiAoA) / Cd**2
        return dpLoD_dpiAoA

    def dpLoD_dpchi(self):
        Cl = self.__Cl()
        Cd = self.__Cd()
        dpCl_dpchi = self.dpCl_dpchi()
        dpCd_dpchi = self.dpCd_dpchi()
        dpLoD_dpchi = (dpCl_dpchi * Cd - Cl * dpCd_dpchi) / Cd**2
        return dpLoD_dpchi

    #-- Display methods
    def __display_info(self):
        print self.get_OC()
        print '\n*** aerodynamic functions and coefficients ***'
        print '  Sref  = ', self.get_Sref()
        print '  Lref  = ', self.get_Lref()
        for i, func in enumerate(self.__F_list_names_calc):
            if func in ['Lift', 'Drag']:
                unit = '(N)'
            else:
                unit = ''
            print '  ' + self.__F_list_names_calc[i] + '\t=\t' + str(self.__F_list_calc[i]) + ' ' + unit
