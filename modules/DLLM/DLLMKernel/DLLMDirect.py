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

import numpy as np
from numpy import array, transpose, outer, ones, zeros, copy, divide, diag, dot
from numpy.linalg import norm, solve
import matplotlib.pylab as plt

from MDOTools.Solvers.newton_raphson_problem import NewtonRaphsonProblem


class DLLMDirect:
    """
    Direct solver for the lifting line wing model
    """
    DEG_TO_RAD = np.pi / 180.
    RAD_TO_DEG = 180. / np.pi

    def __init__(self, LLW, verbose = 0):
        self.__LLW = LLW
        self.__verbose = verbose
        self.__computed = False
        self.__gamma_f_name = None

        # initialize local variables
        self.__init_local_variables()

        # initialize the Newton-Raphson problem
        self.__NRPb = None
        self.__init_Newton_Raphson()

    #-- Accessors
    def get_tag(self):
        return self.__LLW.get_tag()
    
    def get_grad_active(self):
        return self.__LLW.get_grad_active()

    def get_geom(self):
        return self.__LLW.get_geom()

    def get_airfoils(self):
        return self.__LLW.get_airfoils()

    def get_K(self):
        return self.__LLW.get_K()

    def get_dK_dchi(self):
        return self.__LLW.get_dK_dchi()

    def get_N(self):
        return self.get_geom().get_n_sect()

    def get_ndv(self):
        return self.get_geom().get_ndv()

    def get_OC(self):
        return self.__LLW.get_OC()

    def get_iAoA(self):
        return self.__iAoA

    def get_localAoA(self):
        return self.__localAoA

    def get_R(self):
        return self.__R

    def get_dpR_dpiAoA(self):
        return self.__dpR_dpiAoA

    def get_dpR_dpchi(self):
        return self.__dpR_dpchi

    def get_dplocalAoA_dpiAoA(self):
        return self.__dplocalAoA_dpiAoA

    def get_dplocalAoA_dpAoA(self):
        return self.__dplocalAoA_dpAoA

    def get_dplocalAoA_dpchi(self):
        return self.__dplocalAoA_dpchi

    def get_dplocalAoA_dpthetaY(self):
        return self.__dplocalAoA_dpthetaY
    
    def get_convergence_history(self):
        """
        Accessor to the last computation convergence history as a list of residuals normalized by the first iteration residual.
        """
        return self.__residuals_hist

    def is_computed(self):
        return self.__computed

    #-- Setters
    def set_computed(self, bool=True):
        self.__computed = bool

    def set_gamma_file_name(self, gamma_f_name):
        self.__gamma_f_name = gamma_f_name

    #-- Newton-Raphson related methods
    def __init_Newton_Raphson(self):
        N = self.get_N()
        iAoA0 = zeros(N)
        #self.comp_R(iAoA0)
        #self.comp_dpR_dpiAoA(iAoA0)
        self.__NRPb = NewtonRaphsonProblem(iAoA0,self.comp_R,self.comp_dpR_dpiAoA,verbose = self.__verbose)
        self.__NRPb.set_relax_factor(0.99)
        self.__NRPb.set_stop_residual(1.e-9)
        self.__NRPb.set_max_iterations(100)

    def set_relax_factor(self, relax_factor):
        self.__NRPb.set_relax_factor(relax_factor)

    def set_stop_residual(self, residual):
        self.__NRPb.set_stop_residual(residual)

    def set_max_iterations(self, max_it):
        self.__NRPb.set_max_iterations(max_it)

    def set_method(self, method):
        self.__NRPb.set_method(method)

    #-- Computation related methods
    def run(self, iAoA0=None):
        grad_active = self.get_grad_active()
        if iAoA0 is not None:
            self.__NRPb.set_W0(iAoA0)
        self.__NRPb.solve()
        self.set_computed(True)
        self.write_gamma_to_file()
        if grad_active:
            self.comp_dpR_dpchi()
            
    def __init_local_variables(self):
        # Initializing local variables for lifting line computations
        # Residual variables
        N = self.get_N()
        ndv = self.get_ndv()
        grad_active = self.get_grad_active()
        self.__R = zeros([N])
        self.__dpR_dpiAoA = None
        self.__dpR_dpchi = None
        self.__dpR_dpthetaY = None
        self.__dpR_dpAoA = None

        # Angle of attack variables
        self.__localAoA = zeros([N])
        self.__dplocalAoA_dpiAoA = diag(ones([N]))  # This is a constant matrix
        self.__dplocalAoA_dpthetaY = diag(ones([N]))  # This is a constant matrix
        self.__dplocalAoA_dpAoA = ones(N)
        if grad_active:
            self.__dplocalAoA_dpchi = zeros([N, ndv])
        else:
            self.__dplocalAoA_dpchi = None

        # Induced angle of attack variables
        self.__iAoA = None
        self.__iAoANew = None
        if grad_active:
            self.__dpiAoAnew_dpchi = zeros([N, ndv])
        else:
            self.__dpiAoAnew_dpchi = None
        self.__dpiAoAnew_dpiAoA = None

        # Circulation variables
        self.__gamma = zeros(N)
        self.__dpgamma_dpiAoA = None
        self.__dpgamma_dpthetaY = None
        self.__dpgamma_dpAoA = None
        self.__dpgamma_dplocalAoA = zeros([N, N])
        if grad_active:
            self.__dpgamma_dpchi = zeros([N, ndv])
        else:
            self.__dpgamma_dpchi = None
            
    #-- Residual related methods
    def comp_R(self, iAoA):
        self.__iAoA = iAoA
        print 'iAoA = ',iAoA
        self.__compute_localAoA()
        self.__compute_gamma()
        self.__compute_iAoAnew()

        self.__R = self.__iAoA - self.__iAoANew
        

        return self.__R

    def comp_dpR_dpiAoA(self, iAoA):
        N = self.get_N()
#        R = self.comp_R(iAoA)

        # dplocalAoA_dpiAoA is a constant matrix, no need to compute it
        # (self.__dplocalAoA_dpiAoA)
        self.__compute_dpgamma_dpiAoA()
        self.__compute_dpiAoAnew_dpiAoA()

        self.__dpR_dpiAoA = np.diag(ones([N])) - self.__dpiAoAnew_dpiAoA

        return self.__dpR_dpiAoA

    def comp_dpR_dpchi(self):
        self.__compute_dplocalAoA_dpchi()
        self.__compute_dpgamma_dpchi()
        self.__compute_dpiAoAnew_dpchi()
        self.__dpR_dpchi = -self.__dpiAoAnew_dpchi

        return self.__dpR_dpchi

    def comp_dpR_dpthetaY(self):
        K = self.get_K()
        self.__compute_dpgamma_dpthetaY()
        self.__dpR_dpthetaY = -dot(K, self.__dpgamma_dpthetaY)

        return self.__dpR_dpthetaY

    def comp_dpR_dpAoA(self):
        K = self.get_K()
        self.__compute_dpgamma_dpAoA()
        self.__dpR_dpAoA = -dot(K, self.__dpgamma_dpAoA)

        return self.__dpR_dpAoA

    def __compute_dpgamma_dpAoA(self):
        N = self.get_N()
        
        for i in xrange(N):
            self.__dpgamma_dplocalAoA[i,i] = self.get_airfoils()[i].dgamma_dAoA

        self.__dpgamma_dpAoA = dot(self.__dpgamma_dplocalAoA,self.__dplocalAoA_dpAoA)

    def __compute_dpgamma_dpthetaY(self):
        N = self.get_N()
        for i in xrange(N):
            self.__dpgamma_dplocalAoA[i,i] = self.get_airfoils()[i].dgamma_dAoA

        self.__dpgamma_dpthetaY = dot(self.__dpgamma_dplocalAoA,self.__dplocalAoA_dpthetaY)

    def __compute_dplocalAoA_dpchi(self):
        N = self.get_N()
        twist_grad = self.get_geom().get_twist_grad()
        AoA_grad = self.get_geom().get_AoA_grad()
        if AoA_grad is None:
            self.__dplocalAoA_dpchi =  twist_grad
        else:
            for i in xrange(N):
                self.__dplocalAoA_dpchi[i, :] = AoA_grad[:] + twist_grad[i, :]

    def __compute_dpgamma_dpchi(self):
        N = self.get_N()
        for i in xrange(N):
            self.__dpgamma_dpchi[i,:] = self.get_airfoils()[i].dgamma_dchi

        self.__dpgamma_dpchi = self.__dpgamma_dpchi + dot(self.__dpgamma_dplocalAoA, self.__dplocalAoA_dpchi)

    def __compute_dpiAoAnew_dpchi(self):
        K = self.get_K()
        dK_dchi = self.get_dK_dchi()
        ndv = self.get_ndv()
        self.__dpiAoAnew_dpchi = dot(K, self.__dpgamma_dpchi)
        for n in xrange(ndv):
            self.__dpiAoAnew_dpchi[:, n] += dot(dK_dchi[:, :, n], self.__gamma)

    def __compute_localAoA(self):
        N = self.get_N()
        Thetay = self.get_geom().get_thetaY()
        twist = self.get_geom().get_twist()
        AoA = self.get_geom().get_AoA()
        if AoA is None:
            AoA = self.get_OC().get_AoA_rad()
        else:
            self.get_OC().set_AoA(AoA * 180. / np.pi)

        # Why this formula ? twist increases the local airfoil angle of attack normally...
        # self.__localAoA=alpha-iaOa-self.get_wing_geom().get_twist()
        self.__localAoA = AoA + twist + self.__iAoA + Thetay
        print 'AoA = ',AoA
        print 'twist = ',twist
        print 'ThetaY = ',Thetay
        print 'iAoA = ',self.__iAoA
        print 'self.__localAoA = ',self.__localAoA

        for i in xrange(N):
            if self.__localAoA[i] > np.pi / 2. or self.__localAoA[i] < -np.pi / 2.:
                raise Exception("Local angle of attack out of bounds [-pi/2, pi/2]")

    def __compute_gamma(self):
        """
        Update the circulation
        """
        N = self.get_N()
        Mach = self.get_OC().get_Mach()
        for i in xrange(N):
            self.get_airfoils()[i].compute(self.__localAoA[i],Mach)
            self.__gamma[i] = self.get_airfoils()[i].gamma
            if np.isnan(self.get_airfoils()[i].gamma):
                af = self.get_airfoils()[i]
                print 'Lref = ',af.get_Lref()
                print 'Sref = ',af.get_Sref()
                print 'Cl = ',af.Cl

    def __compute_dpgamma_dpiAoA(self):
        N = self.get_N()
        for i in xrange(N):
            self.__dpgamma_dplocalAoA[i,i] = self.get_airfoils()[i].dgamma_dAoA

        self.__dpgamma_dpiAoA = dot(self.__dpgamma_dplocalAoA, self.__dplocalAoA_dpiAoA)

    def __compute_iAoAnew(self):
        '''
        Computes the induced angle on an airfoil for a given circulation on the wing.
        '''
        K = self.get_K()
        self.__iAoANew = dot(K, self.__gamma)

    def __compute_dpiAoAnew_dpiAoA(self):
        """
        Computes the derivative dpiAoAnew_dpiAoA
        """
        K = self.get_K()
        self.__dpiAoAnew_dpiAoA = dot(K, self.__dpgamma_dpiAoA)

    def write_gamma_to_file(self):
        '''
        Writes the circulation repartition in a file
        '''
        y = self.get_geom().get_XYZ()[1, :]
        if self.__gamma_f_name is None:
            gamma_f_name = 'gamma.dat'
        else:
            gamma_f_name = self.__gamma_f_name

        mod_gamma_f_name = self.get_tag() + '_' + gamma_f_name

        fid = open(mod_gamma_f_name, 'w')
        line = "#Slice\t%24s" % "y" + "\t%24s" % "Circulation" + "\n"
        fid.write(line)
        i = 0
        for i in range(len(self.__gamma)):
            line = str(i) + "\t%24.16e" % y[i] + \
                "\t%24.16e" % self.__gamma[i] + "\n"
            fid.write(line)
        fid.close()

        fid = open(self.get_tag() + '_iAoA.dat', 'w')
        line = "#Slice\t%24s" % "y" + "\t%24s" % "iAoA" + \
            "\t%24s" % "iAoA_deg" + "\n"
        fid.write(line)
        i = 0
        for i in range(len(self.__gamma)):
            line = str(i) + "\t%24.16e" % y[i] + "\t%24.16e" % self.__iAoA[
                i] + "\t%24.16e" % (self.__iAoA[i] * self.RAD_TO_DEG) + "\n"
            fid.write(line)
        fid.close()
    
    def plot(self):
        name = self.get_tag()
        Y_list = self.get_geom().get_XYZ()[1,:]
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('gamma')
        plt.plot(Y_list,self.__gamma)
        plt.rc("font", size=14)
        plt.savefig(name+"_gamma_distrib.png",format='png')
        plt.close()
        
        plt.xlim(1.1*Y_list[0], 1.1*Y_list[-1])
        plt.xlabel('y')
        plt.ylabel('iAoA')
        plt.plot(Y_list,self.__iAoA*self.RAD_TO_DEG)
        plt.rc("font", size=14)
        plt.savefig(name+"_iAoA_distrib.png",format='png')
        plt.close()
