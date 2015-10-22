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

from numpy import dot, zeros
from numpy.linalg import norm, solve

from MDOTools.Solvers.adjoint_problem import AdjointProblem

class DLLMAdjoint:
    def __init__(self, LLW, verbose = 1):
        """
        Adjoint module for the lifting line model
        """
        self.__LLW = LLW
        self.__verbose = verbose
        self.__computed = False

        self.__adj_list           = None
        self.__adj_conv_corr_list = None
        self.__dF_list_dchi       = None
        
    #-- accessors
    def get_tag(self):
        return self.__LLW.get_tag()
    
    def get_adjoint_list(self):
        return self.__adj_list
    
    def get_adjoint_convergence_correction_list(self):
        return self.__adj_conv_corr_list
    
    def get_dF_list_dchi(self):
        return self.__dF_list_dchi
        
    #-- computed related methods
    def is_computed(self):
        return self.__computed
    
    def set_computed(self, bool=True):
        self.__computed = bool
    
    #-- Residual accessors
    def get_R(self):
        return self.__LLW.get_R()
    
    def get_dpR_dpW(self):
        return self.__LLW.get_dpR_dpW()
    
    def get_dpR_dpchi(self):
        return self.__LLW.get_dpR_dpchi()
    
    #-- Post-processing accessors
    def get_F_list_names(self):
        return self.__LLW.get_F_list_names()
    
    def get_dpF_list_dpW(self):
        return self.__LLW.get_dpF_list_dpW()
    
    def get_dpF_list_dpchi(self):
        return self.__LLW.get_dpF_list_dpchi()
    
    #-- Run adjoint problem for each cost function
    def run(self):
        F_list_names = self.get_F_list_names()
        
        self.__adj_list           = []
        self.__adj_conv_corr_list = []
        self.__dF_list_dchi       = []
        
        R        = self.get_R()
        dpRdpW   = self.get_dpR_dpW()
        dpRdpchi = self.get_dpR_dpchi()
        
        dpFdpW_list   = self.get_dpF_list_dpW()
        dpFdpchi_list = self.get_dpF_list_dpchi()
        
        if self.__verbose == 0 :
            print "Running adjoint"
        for i, F_name in enumerate(F_list_names):
            if self.__verbose > 0 :
                print 'Run adjoint problem for func = '+str(F_name)
            dpFdpW   = dpFdpW_list[i]
            dpFdpchi = dpFdpchi_list[i]
            
            AdjPb = AdjointProblem(R, dpRdpW, dpRdpchi, dpFdpW, dpFdpchi)
            AdjPb.solve()
            
            self.__adj_list.append(AdjPb.get_adjoint_state())
            self.__adj_conv_corr_list.append(AdjPb.get_convergence_correction())
            self.__dF_list_dchi.append(AdjPb.get_dFdchi())
            
            if self.__verbose > 0 : print '  - Convergence adjoint correction for '+str(F_name)+' = '+str(self.__adj_conv_corr_list[i])
