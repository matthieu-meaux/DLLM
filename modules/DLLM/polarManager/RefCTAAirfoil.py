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
# @author : Matthieu Meaux

# - Local imports -
import string
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d

from airfoil import Airfoil

class RefCTAAirfoil(Airfoil):
    """
    An analytic airfoil based on linear theory
    """
    THICKNESS_CORRECTION = 0.7698
    
    def __init__(self, OC, Sref=1., Lref=1., y_pos=None, grad_active=False):
        '''
        Constructor for airfoils
        @param y_pos : span-wise position for interpolation
        @param Sref  : reference surface
        @param Lref  : reference length
        '''
        
        Airfoil.__init__(self, OC, Lref=Lref, Sref=Sref, grad_active=grad_active)
        
        self.y_pos = y_pos
        
        self.__y_def_list    = None
        self.__file_def_list = None
        self.__interp_list   = None
        
        self.__index_p       = None
        self.__index_m       = None
        self.__fact_p        = None
        self.__fact_m        = None

    #-- Setters
    def set_y_pos(self, y_pos):
        self.y_pos = y_pos
        
    def set_y_def_list(self, y_def_list):
        self.__y_def_list = y_def_list
        
    def set_file_def_list(self, file_def_list):
        self.__file_def_list = file_def_list
        
    def set_dict_info(self, dict_info):
        self.__dict_info = dict_info
    
    def init_interpolators(self):
#         print 'self.__y_def_list = ',self.__y_def_list 
#         print 'self.__file_def_list = ',self.__file_def_list
        self.__interp_list = []
        #-- Check length
        if len(self.__y_def_list) != len(self.__file_def_list):
            raise Exception,'RefCTAAirfoil: Inconsistent input data'
        for filename in self.__file_def_list:
            #-- Read file
            AoA_list ,CL_list, CDw_list, CDvp_list, CDf_list, Cmy_list = self.__read_file(filename)
            out_coeffs=np.array([CL_list,CDw_list,CDvp_list,CDf_list,Cmy_list])
            interp_coeffs = interp1d(AoA_list, out_coeffs, kind='cubic')
            self.__interp_list.append(interp_coeffs)
        
#         #-- Test plotting of polars
#         AoA_test=np.linspace(-6.0, 6.0, 50)
#         out_list = []
#         for i,y in enumerate(self.__y_def_list):
#             out_list.append(self.__interp_list[i](AoA_test))
#         
#         #-- CL plots
#         plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
#         for i,y in enumerate(self.__y_def_list):
#             out_data = out_list[i]
#             plt.plot(AoA_test ,out_data[0])
#         plt.xlabel('AoA (deg)')
#         plt.ylabel('CL')
#         plt.rc("font", size=14)
#         plt.savefig('CL_polars.png',format='png')
#         plt.close()
# 
#         #-- CDw plots
#         plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
#         for i,y in enumerate(self.__y_def_list):
#             out_data = out_list[i]
#             plt.plot(AoA_test ,out_data[1])
#         plt.xlabel('AoA (deg)')
#         plt.ylabel('CDw')
#         plt.rc("font", size=14)
#         plt.savefig('CDw_polars.png',format='png')
#         plt.close()
#         
#         #-- CDvp plots
#         plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
#         for i,y in enumerate(self.__y_def_list):
#             out_data = out_list[i]
#             plt.plot(AoA_test ,out_data[2])
#         plt.xlabel('AoA (deg)')
#         plt.ylabel('CDvp')
#         plt.rc("font", size=14)
#         plt.savefig('CDvp_polars.png',format='png')
#         plt.close()
#         
#         #-- CDf plots
#         plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
#         for i,y in enumerate(self.__y_def_list):
#             out_data = out_list[i]
#             plt.plot(AoA_test ,out_data[3])
#         plt.xlabel('AoA (deg)')
#         plt.ylabel('CDf')
#         plt.rc("font", size=14)
#         plt.savefig('CDf_polars.png',format='png')
#         plt.close()
#         
#         #-- Cmy_list plots
#         plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
#         for i,y in enumerate(self.__y_def_list):
#             out_data = out_list[i]
#             plt.plot(AoA_test ,out_data[4])
#         plt.xlabel('AoA (deg)')
#         plt.ylabel('Cmy')
#         plt.rc("font", size=14)
#         plt.savefig('Cmy_polars.png',format='png')
#         plt.close()
        
    def init_interp_factors(self):
        test_y = abs(self.y_pos)
        i=0
        test = 0
        while test_y > test:
            i = i+1
            test = self.__y_def_list[i]
        self.__index_m   = i-1
        self.__index_p   = i
#         print 'y_pos = ',self.y_pos, 'y_m = ',self.__y_def_list[self.__index_m],'y_p = ',self.__y_def_list[self.__index_p ]
        
        fact = (self.__y_def_list[self.__index_p ]-test_y)/(self.__y_def_list[self.__index_p]-self.__y_def_list[self.__index_m])
        self.__fact_m = fact
        self.__fact_p   = (1.-fact)

    #-- Methods to compute aero coefficients
    def comp_aero_coeffs(self, AoA, Mach):
        AoA_deg = AoA*180./np.pi
        #-- Mach is ignored for this airfoil
        out_data_m = self.__interp_list[self.__index_m](AoA_deg)
        out_data_p = self.__interp_list[self.__index_p](AoA_deg)
        out_data_m_d = self.__interp_list[self.__index_m](AoA_deg+0.05)
        out_data_p_d = self.__interp_list[self.__index_p](AoA_deg+0.05)
#         print 'out_data_m = ',out_data_m
#         print 'out_data_p = ',out_data_p
        self.Cl   = self.__fact_m*out_data_m[0] + self.__fact_p*out_data_p[0]
        Cl_d = self.__fact_m*out_data_m_d[0] + self.__fact_p*out_data_p_d[0]
        self.Cdw  = self.__fact_m*out_data_m[1] + self.__fact_p*out_data_p[1]
        self.Cdvp = self.__fact_m*out_data_m[2] + self.__fact_p*out_data_p[2]
        self.Cdf  = self.__fact_m*out_data_m[3] + self.__fact_p*out_data_p[3]
        
        self.dCl_dAoA = (Cl_d-self.Cl)/(0.05*np.pi/180.)
        
        Cmy = self.__fact_m*out_data_m[4] + self.__fact_p*out_data_p[4]
        
        self.pcop  = 0.25
        if self.Cl != 0:
            delta_d = Cmy/abs(self.Cl)
            
        self.pcop = self.pcop-delta_d
        
#         print 'delta_d = ',delta_d
#         
#         print 'Cl_d = ',Cl_d
#         print 'self.Cl = ',self.Cl
#         print 'self.Cdw = ',self.Cdw
#         print 'self.Cdvp = ',self.Cdvp
#         print 'self.Cdf = ',self.Cdf
#         
#         print 'self.dCl_dAoA = ',self.dCl_dAoA,2.*np.pi
#         
#         print 'Cmy = ',Cmy
#         print 'self.pcop = ',self.pcop
            
    def get_scaled_copy(self, OC=None, Sref=None, Lref=None, rel_thick=None, grad_active=True):
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        if rel_thick is None:
            rel_thick = self.get_rel_thick()
        if OC is None:
            OC = self.get_OC()
        return AnalyticAirfoil(OC, self.__AoA0_deg, Sref=Sref, Lref=Lref, rel_thick=rel_thick, Ka=self.__Ka, grad_active=self.is_grad_active())
    
    #-- Private methods
    def __read_file(self, filename):
        fid = open(filename,'r')
        all_lines = fid.readlines()
        n_lines = len(all_lines)
        #print 'n_lines = ',n_lines
        #-- get_tag_list
        tag_list = string.split(all_lines[0])
        #print 'tag_list = ',tag_list
        
        index_AoA = tag_list.index('AoA')
        #print 'index AoA = ',index_AoA
        
        index_CL = tag_list.index('CL')
        #print 'index CL = ',index_CL
        
        index_CDw = tag_list.index('CDw')
        #print 'index CDw = ',index_CDw
        
        index_CDvp = tag_list.index('CDvp')
        #print 'index CDvp = ',index_CDvp
        
        index_CDf = tag_list.index('CDf')
        #print 'index CDf = ',index_CDf
        
        index_Cmy = tag_list.index('Cmy')
        #print 'index Cmy = ',index_Cmy
        
        AoA_list  = []
        CL_list   = []
        CDw_list  = []
        CDvp_list = []
        CDf_list  = []
        Cmy_list  = []
        
        for i in xrange(1,n_lines):
            line = all_lines[i]
            values = string.split(line)
            AoA_list.append(float(values[index_AoA]))
            CL_list.append(float(values[index_CL]))
            CDw_list.append(float(values[index_CDw]))
            CDvp_list.append(float(values[index_CDvp]))
            CDf_list.append(float(values[index_CDf]))
            Cmy_list.append(float(values[index_Cmy]))
            
#         print 'AoA_list = ',AoA_list
#         print 'CL_list = ',CL_list
#         print 'CDw_list = ',CDw_list
#         print 'CDvp_list = ',CDvp_list
#         print 'CDf_list = ',CDf_list
#         print 'Cmy_list = ',Cmy_list
        
        return AoA_list, CL_list, CDw_list, CDvp_list, CDf_list, Cmy_list
            
            