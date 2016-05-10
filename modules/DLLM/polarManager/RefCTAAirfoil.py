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
        
        self.__dict_info = {}
        
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
        print 'self.__y_def_list = ',self.__y_def_list 
        print 'self.__file_def_list = ',self.__file_def_list
        self.__interp_list = []
        #-- Check length
        if len(self.__y_def_list) != len(self.__file_def_list):
            raise Exception,'RefCTAAirfoil: Inconsistent input data'
        for filename in self.__file_def_list:
            #-- Read file
            AoA_list ,CL_list, CDw_list, CDvp_list, CDf_list, Cmy_list = self.__read_file(filename)
            out_coeffs=np.array([CL_list,CDw_list,CDvp_list,CDf_list,Cmy_list])
            print out_coeffs.shape
            interp_coeffs = interp1d(AoA_list, out_coeffs, kind='cubic')
            self.__interp_list.append(interp_coeffs)
            print 'CL_list = ',CL_list
        
        #-- Test plotting of polars
        AoA_test=np.linspace(-6.0, 6.0, 50)
        out_list = []
        for i,y in enumerate(self.__y_def_list):
            out_list.append(self.__interp_list[i](AoA_test))
        
        #-- CL plots
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
        for i,y in enumerate(self.__y_def_list):
            out_data = out_list[i]
            plt.plot(AoA_test ,out_data[0])
        plt.xlabel('AoA (deg)')
        plt.ylabel('CL')
        plt.rc("font", size=14)
        plt.savefig('CL_polars.png',format='png')
        plt.close()

        #-- CDw plots
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
        for i,y in enumerate(self.__y_def_list):
            out_data = out_list[i]
            plt.plot(AoA_test ,out_data[1])
        plt.xlabel('AoA (deg)')
        plt.ylabel('CDw')
        plt.rc("font", size=14)
        plt.savefig('CDw_polars.png',format='png')
        plt.close()
        
        #-- CDvp plots
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
        for i,y in enumerate(self.__y_def_list):
            out_data = out_list[i]
            plt.plot(AoA_test ,out_data[2])
        plt.xlabel('AoA (deg)')
        plt.ylabel('CDvp')
        plt.rc("font", size=14)
        plt.savefig('CDvp_polars.png',format='png')
        plt.close()
        
        #-- CDf plots
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
        for i,y in enumerate(self.__y_def_list):
            out_data = out_list[i]
            plt.plot(AoA_test ,out_data[3])
        plt.xlabel('AoA (deg)')
        plt.ylabel('CDf')
        plt.rc("font", size=14)
        plt.savefig('CDf_polars.png',format='png')
        plt.close()
        
        #-- Cmy_list plots
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3, mode="expand", borderaxespad=0.,prop={'size':12})
        for i,y in enumerate(self.__y_def_list):
            out_data = out_list[i]
            plt.plot(AoA_test ,out_data[4])
        plt.xlabel('AoA (deg)')
        plt.ylabel('Cmy')
        plt.rc("font", size=14)
        plt.savefig('Cmy_polars.png',format='png')
        plt.close()
        
    def init_interp_factors(self):
        test_y = abs(self.y_pos)
        i=0
        test = 0
        while test_y > test:
            i = i+1
            test = self.__y_def_list[i]
        print 'y_pos = ',self.y_pos, 'y_m = ',self.__y_def_list[i-1],'y_p = ',self.__y_def_list[i]
        
        fact = (self.__y_def_list[i]-self.y_pos)/(self.__y_def_list[i]-self.__y_def_list[i-1])
        self.__fact_p   = (1.-fact)
        self.__fact_m = fact
        print 'fact m = ',self.__fact_m
        print 'fact p = ',self.__fact_p

        
    #-- Methods to compute aero coefficients
    def comp_aero_coeffs(self, AoA, Mach):
        
        
        
        self.Cl   = 0.
        self.Cdw  = 0.
        self.Cdvp = 0.
        self.Cdf  = 0.
        #To fill: Cl, Cdw, Cdf, Cdvp
#         OC = self.get_OC()
#         c  = OC.get_c()
#         nu = OC.get_nu()
#         L     = self.get_Lref()
#         sweep = self.get_sweep()
#         toc   = self.get_rel_thick()
#         toc2  = self.get_rel_thick()/np.cos(sweep)
#         Mach_normal = Mach*np.cos(sweep)
#         
#         if self.is_grad_active():
#             dL  = self.get_Lref_grad()
#             dsweep = self.get_sweep_grad()
#             dtoc   = self.get_rel_thick_grad()
#             dtoc2   = (dtoc*np.cos(sweep)+toc*np.sin(sweep)*dsweep)/np.cos(sweep)**2
#             dMach_normal = -Mach*np.sin(sweep)*dsweep
# 
#         
#         #-- compute dCl_dAoA         
#         # Base slope with Prandtl correction
#         if Mach_normal<0.9:
#             base_slope=2.*np.pi/np.sqrt(1.-Mach_normal**2)
#             if self.is_grad_active():
#                 dbase_slope = 2.*np.pi*(Mach_normal*dMach_normal)/((1.-Mach_normal**2)*np.sqrt(1.-Mach_normal**2))
#         elif 0.9<=Mach_normal<=1.1:
#             s_sub = 2.*np.pi/np.sqrt(abs(1.-0.9**2))
#             s_sup = 4./np.sqrt(abs(1.1**2-1.))
#             fact = (1.1-Mach_normal)/0.2
#             base_slope = s_sub*fact+s_sup*(1.-fact)
#             if self.is_grad_active():
#                 dfact = -dMach_normal/0.2
#                 dbase_slope = s_sub*dfact+s_sup*(1.-dfact)
#         else:
#             base_slope = 4./np.sqrt(Mach_normal**2-1.)
#             if self.is_grad_active():
#                 dbase_slope = -4.*(Mach_normal*dMach_normal)/((1.-Mach_normal**2)*np.sqrt(1.-Mach_normal**2))
# 
#         # thickness correction
#         thick_corr = 1. + self.THICKNESS_CORRECTION*toc/np.cos(sweep)
#         if self.is_grad_active():
#             dthick_corr = self.THICKNESS_CORRECTION*(dtoc*np.cos(sweep)+toc*np.sin(sweep)*dsweep)/np.cos(sweep)**2
#         # sweep correction
#         sweep_corr = np.cos(sweep)**2
#         if self.is_grad_active():
#             dsweep_corr = -2.*np.sin(sweep)*np.cos(sweep)*dsweep
#             
#         self.dCl_dAoA = base_slope*thick_corr*sweep_corr
#         
#         if self.is_grad_active():
#             d2Cl_dAoAdchi = dbase_slope*thick_corr*sweep_corr\
#                           + base_slope*dthick_corr*sweep_corr\
#                           + base_slope*thick_corr*dsweep_corr
#                           
#         #-- Compute Cl
#         self.Cl = self.dCl_dAoA*(AoA-self.__AoA0)
#         
#         if np.isnan(self.Cl):
#             print 'dCl_dAoA = ',self.dCl_dAoA
#             print 'sweep_corr = ',sweep_corr
#             print 'thick_corr = ',thick_corr
#             print 'toc = ',toc
#         
#         #-- Compute dCl_dchi
#         if self.is_grad_active():
#             self.dCl_dchi = d2Cl_dAoAdchi*(AoA-self.__AoA0)
#     
#                       
#         #-- Compute Cdw
#         Mdd    = self.__Ka/np.cos(sweep) - toc/(np.cos(sweep))**2 - self.Cl/(10.*np.cos(sweep)**3)
#         Mcrit  = Mdd - (0.1/80)**(1./3.)
#         if Mach < Mcrit:
#             self.Cdw = 0.0
#         else:
#             self.Cdw = 20.*(Mach-Mcrit)**4
#             
#         
#         #-- Compute dCdw_dAoA
#         dMdd_dAoA   = -self.dCl_dAoA/(10.*np.cos(sweep)**3)
#         dMcrit_dAoA =  dMdd_dAoA
#         if Mach < Mcrit:
#             self.dCdw_dAoA = 0.0
#         else:
#             self.dCdw_dAoA = -80.*(Mach-Mcrit)**3*dMcrit_dAoA
#     
#         #-- compute dCdw_dchi
#         if self.is_grad_active():
#             dMdd   = self.__Ka*np.sin(sweep)*dsweep/(np.cos(sweep))**2 \
#                - dtoc/(np.cos(sweep))**2 - toc*(2.*np.sin(sweep)*dsweep)/np.cos(sweep)**3 \
#                - self.dCl_dchi/(10.*np.cos(sweep)**3) - self.Cl*(3.*np.sin(sweep)*dsweep)/(10.*np.cos(sweep)**4)
#             dMcrit = dMdd
#              
#             if Mach < Mcrit:
#                 self.dCdw_dchi = 0.0
#             else:
#                 self.dCdw_dchi = -80.*(Mach-Mcrit)**3*dMcrit
#                 
#         #-- compute Cdf
#         Re = Mach*np.cos(sweep)*c*L/nu
#         thick_coeff = 1.+2.1*toc2 # a rough guess since 1.21 for 10% relative thickness
#         if Re < 1.e-12:
#             # Prevent division by 0. Drag is null at zero Re number anyway
#             self.Cdf = 0.
#         elif Re < 1.e5:
#             # Laminar flow
#             self.Cdf=1.328/np.sqrt(Re)*thick_coeff
#         else:
#             # Turbulent flow
#             self.Cdf=0.074*Re**(-0.2)*thick_coeff
#             
#         #-- Compute dCdf_dAoA
#         self.dCdf_dAoA = 0.
#         
#         #-- Compute dCdf_dchi
#         if self.is_grad_active():
#             dRe = Mach*np.cos(sweep)*c*dL/nu-Mach*np.sin(sweep)*c*L/nu*dsweep
#             dthick_coeff = 2.1*dtoc2
#             if Re < 1.e-12:
#                 # Prevent division by 0. Drag is null at zero Re number anyway
#                 self.dCdf_dchi = 0.
#             elif Re < 1.e5:
#                 # Laminar flow
#                 self.dCdf_dchi=-0.664/(Re**1.5)*dRe*thick_coeff+1.328/np.sqrt(Re)*dthick_coeff
#             else:
#                 # Turbulent flow
#                 self.dCdf_dchi=-0.0148*Re**(-1.2)*dRe*thick_coeff+0.074*Re**(-0.2)*dthick_coeff
#                 
#         #-- compute Cdvp
#         Cdvp_min = 60.*toc2**4*self.Cdf
#         self.Cdvp = Cdvp_min + 2.*(Cdvp_min + self.Cdf)*self.Cl**2
#         #-- compute dCdvp_dAoA
#         # dCdvp_min_dAoA = 0. since dCdf_dAoA = 0. 
#         self.dCdvp_dAoA = 4.*(Cdvp_min + self.Cdf)*self.Cl*self.dCl_dAoA
#         
#         if self.is_grad_active():
#             dCdvp_min= 60.*(4.*dtoc2*toc2**3*self.Cdf+toc2**4*self.dCdf_dchi)
#             self.dCdvp_dchi = dCdvp_min + 2.*(dCdvp_min+self.dCdf_dchi)*self.Cl**2+4.*(Cdvp_min + self.Cdf)*self.Cl*self.dCl_dchi
            
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
        print 'n_lines = ',n_lines
        #-- get_tag_list
        tag_list = string.split(all_lines[0])
        print 'tag_list = ',tag_list
        
        index_AoA = tag_list.index('AoA')
        print 'index AoA = ',index_AoA
        
        index_CL = tag_list.index('CL')
        print 'index CL = ',index_CL
        
        index_CDw = tag_list.index('CDw')
        print 'index CDw = ',index_CDw
        
        index_CDvp = tag_list.index('CDvp')
        print 'index CDvp = ',index_CDvp
        
        index_CDf = tag_list.index('CDf')
        print 'index CDf = ',index_CDf
        
        index_Cmy = tag_list.index('Cmy')
        print 'index Cmy = ',index_Cmy
        
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
            
        print 'AoA_list = ',AoA_list
        print 'CL_list = ',CL_list
        print 'CDw_list = ',CDw_list
        print 'CDvp_list = ',CDvp_list
        print 'CDf_list = ',CDf_list
        print 'Cmy_list = ',Cmy_list
        
        return AoA_list, CL_list, CDw_list, CDvp_list, CDf_list, Cmy_list
            
            