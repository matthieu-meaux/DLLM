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

    #-- Setters
    def set_y_pos(self, y_pos):
        self.y_pos = y_pos
        
    def set_y_def_list(self, y_def_list):
        self.__y_def_list = y_def_list
        
    def set_file_def_list(self, file_def_list):
        self.__file_def_list = file_def_list
    
    def set_interp_list(self, interp_list):
        self.__interp_list = interp_list
    
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

    #-- Methods to compute aero coefficients
    def comp_aero_coeffs(self, AoA, Mach):
        AoA_deg = AoA*180./np.pi

        #-- Mach is ignored for this airfoil
        out_data_list = []
        for i,y in enumerate(self.__y_def_list):
            out_data = self.__interp_list[i](AoA_deg)
            out_data_list.append(out_data)
        out_data_list = np.array(out_data_list)
        out_data_list = out_data_list.T
        out_interp = interp1d(self.__y_def_list, out_data_list, kind='slinear')
        
        
#         #-- CL plots
#         plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
#         y_range = np.linspace(0.,19.0,100)
#         out_Cl = out_interp(y_range)[0]
#         plt.plot(y_range ,out_Cl)
#         plt.plot(self.__y_def_list, out_data_list[0])
#         plt.xlabel('y)')
#         plt.ylabel('CL')
#         plt.rc("font", size=14)
#         plt.savefig('CL_spanwise.png',format='png')
#         plt.close()
#         
        out_data_list_d = []
        for i,y in enumerate(self.__y_def_list):
            out_data_d = self.__interp_list[i](AoA_deg+0.05)
            out_data_list_d.append(out_data_d)
        out_data_list_d = np.array(out_data_list_d)
        out_data_list_d = out_data_list_d.T
        
        out_interp_d = interp1d(self.__y_def_list, out_data_list_d, kind='cubic')
        
        out_coeffs   = out_interp(abs(self.y_pos))
        out_coeffs_d = out_interp_d(abs(self.y_pos))
        
        Cmy = out_coeffs[4]   
        
        Cl_d = out_coeffs_d[0]
        self.dCl_dAoA = (out_coeffs_d[0]-out_coeffs[0])/(0.05*np.pi/180.)
        
        self.Cl   = out_coeffs[0]
        self.Cdw  = out_coeffs[1]*1e-4
        self.Cdvp = out_coeffs[2]*1e-4
        self.Cdf  = out_coeffs[3]*1e-4
        
        self.pcop  = 0.25
        if self.Cl != 0:
            delta_d = Cmy/self.Cl
            
        self.pcop = 0.25-delta_d
        
    def get_scaled_copy(self, OC=None, Sref=None, Lref=None, grad_active=False):
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        if OC is None:
            OC = self.get_OC()
            
        CopyAF = RefCTAAirfoil(OC, Sref=Sref, Lref=Lref, grad_active=self.is_grad_active())
        CopyAF.set_y_def_list(self.__y_def_list)
        CopyAF.set_interp_list(self.__interp_list)
        return CopyAF
    
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
            if len(values) >0:
                AoA_list.append(float(values[index_AoA]))
                CL_list.append(float(values[index_CL]))
                CDw_list.append(float(values[index_CDw]))
                CDvp_list.append(float(values[index_CDvp]))
                CDf_list.append(float(values[index_CDf]))
                Cmy_list.append(float(values[index_Cmy]))
        
        return AoA_list, CL_list, CDw_list, CDvp_list, CDf_list, Cmy_list
            
            