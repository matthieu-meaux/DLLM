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
# - imports -
from MDOTools.OC.operating_condition import OperatingCondition
from DLLM.DLLMGeom.wing_broken import Wing_Broken
from DLLM.DLLMGeom.wing_elliptic import Wing_Elliptic
from DLLM.DLLMGeom.wing_straight import Wing_Straight
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from DLLM.DLLMKernel.DLLMTargetCl import DLLMTargetCl
from DLLM.DLLMKernel.DLLMTargetLift import DLLMTargetLift
import numpy
import cPickle

class DLLMWrapper():
    ERROR_MSG = 'ERROR in DLLMWrapper.'
    WARNING_MSG = 'WARNING in DLLMWrapper.'
    POS_SOLVER = ['Solver','TargetCl','TargetLift']
    POS_FMT = ['list','numpy']

    def __init__(self, tag, verbose=1):
        """
        Wrapper for the DLLM solver
        """
        self.__tag          = tag

        self.__verbose      = verbose
        self.__OC           = None # Operating condition
        self.__wing_param   = None # Wing_param class
        self.__DLLM_solver  = None # DLLM Solver class
        
        self.__config_dict  = None
        
        self.__out_format   = 'list'
        self.__grad_format  = 'list'
        
        self.__AoA_id       = 'AoA'
        
        self.__F_list       = None
        self.__F_list_grad  = None
        
    #-- Accessors
    def get_OC(self):
        return self.__OC
    
    def get_wing_param(self):
        return self.__wing_param
    
    def get_DLLM_solver(self):
        return self.__DLLM_solver
    
    def get_tags_x0_and_bounds(self):
        tags=self.__wing_param.get_dv_id_list()
        x0=self.__wing_param.get_dv_array()
        bounds=self.__wing_param.get_bounds_array()
        return tags,x0,bounds
    
    def get_x0_and_bounds(self):
        x0=self.__wing_param.get_dv_array()
        bounds=self.__wing_param.get_bounds_array()
        return x0,bounds
    
    def get_x0(self):
        return self.get_x()
    
    def get_x(self):
        x=self.__wing_param.get_dv_array()
        return x
    
    def get_F_list_names(self):
        return self.__DLLM_solver.get_F_list_names()
    
    def get_F_list(self):
        return self.__F_list
    
    def get_F_list_grad(self):
        return self.__F_list_grad
    
    def get_F_list_and_grad(self):
        return self.__F_list, self.__F_list_grad
    
    def get_F_value(self, F_id):
        index = self.get_F_list_names().index(F_id)
        return self.get_F_list()[index]
        
    #-- Setters
    def set_AoA_id(self, AoA_id):
        self.__AoA_id = AoA_id
        
    def set_out_format(self, format):
        WARNING_MSG=self.WARNING_MSG+'set_out_format: '
        if format not in self.POS_FMT:
            print WARNING_MSG+'format = '+str(format)+' not in '+str(self.POS_FMT)+'. Set to default out format = list'
            format='list'
        self.__out_format = format
    
    def set_grad_format(self, format):
        WARNING_MSG=self.WARNING_MSG+'set_grad_format: '
        if format not in self.POS_FMT:
            print WARNING_MSG+'format = '+str(format)+' not in '+str(self.POS_FMT)+'. Set to default grad format = list'
            format='list'
        self.__grad_format = format
    
    #-- Public methods
    def configure(self, config_dict):
        self.__config_dict = config_dict
        self.__config_OC()
        self.__config_param()
        self.__config_DLLM()
        
    def run(self, x):
        self.__wing_param.update_from_x_list(x)
        self.__DLLM_solver.set_geom(self.__wing_param)
        F_list=self.analysis()
        return F_list
    
    def run_grad(self, x):
        self.__wing_param.update_from_x_list(x)
        self.__DLLM_solver.set_geom(self.__wing_param)
        F_list_grad=self.analysis_grad()
        return F_list_grad
    
    def run_and_grad(self, x):
        self.__wing_param.update_from_x_list(x)
        self.__DLLM_solver.set_geom(self.__wing_param)
        F_list,F_list_grad = self.analysis_and_grad()
        return F_list,F_list_grad
    
    def analysis(self):
        if self.__verbose > 0:
            print self.__wing_param
        self.__DLLM_solver.run_direct()
        self.__DLLM_solver.run_post()
        F_list = self.__DLLM_solver.get_F_list()
        self.__DLLM_solver.export_F_list()
        if self.__out_format == 'list':
            F_list=F_list.tolist()
        self.__F_list = F_list
        return F_list
    
    def analysis_grad(self):
        if self.__verbose > 0:
            print self.__wing_param
        self.__DLLM_solver.run_direct()
        self.__DLLM_solver.run_post()
        self.__DLLM_solver.run_adjoint()
        F_list_grad=self.__DLLM_solver.get_dF_list_dchi()
        self.__DLLM_solver.export_dF_list_dchi()
        if self.__grad_format == 'numpy':
            F_list_grad=numpy.array(F_list_grad)
        self.__F_list_grad = F_list_grad
        return F_list_grad
    
    def analysis_and_grad(self):
        if self.__verbose > 0:
            print self.__wing_param
        self.__DLLM_solver.run_direct()
        self.__DLLM_solver.run_post()
        self.__DLLM_solver.run_adjoint()
        F_list = self.__DLLM_solver.get_F_list()
        F_list_grad=self.__DLLM_solver.get_dF_list_dchi()
        self.__DLLM_solver.export_F_list()
        self.__DLLM_solver.export_dF_list_dchi()
        if self.__out_format == 'list':
            F_list=F_list.tolist()
        if self.__grad_format == 'numpy':
            F_list_grad=numpy.array(F_list_grad)
        self.__F_list = F_list
        self.__F_list_grad = F_list_grad
        return F_list,F_list_grad
    
    def export_results(self):
        fid=open(self.__tag+'.res.dat','w')
        res=[self.__F_list,self.__F_list_grad]
        cPickle.dump(res, fid)
        fid.close()
        
    def import_results(self):
        fid=open(self.__tag+'.res.dat','r')
        res=cPickle.load(fid)
        fid.close()
        self.__F_list=res[0]
        self.__F_list_grad=res[1]
        
    #-- Private methods
    def __config_OC(self):
        """
        Set up the 
        """
        self.__OC = OperatingCondition(self.__tag+'.OC')
        self.__OC.config_from_dict(self.__config_dict)
        
    def __config_param(self):
        """
        Set up the __wing_param attribute.
        -The class associated to __wing_param depends on the chosen parameterization type
        -The method needs to be updated if a new parameterization is implemented
        """
        input_keys=self.__config_dict.keys()
        
        geom_type_key=self.__tag+'.param.geom_type'
        if geom_type_key in input_keys:
            geom_type = self.__config_dict[geom_type_key]
        else:
            geom_type= 'Broken'
            
        n_sect_key=self.__tag+'.param.n_sect'
        if n_sect_key in input_keys:
            n_sect = self.__config_dict[n_sect_key]
        else:
            n_sect = 20
        
        if   geom_type == 'Broken':
            self.__wing_param = Wing_Broken(self.__tag+'.param', n_sect=n_sect)
        elif geom_type == 'Elliptic':
            self.__wing_param = Wing_Elliptic(self.__tag+'.param', n_sect=n_sect)
        elif geom_type == 'Straight':
            self.__wing_param = Wing_Straight(self.__tag+'.param', n_sect=n_sect)
            
        self.__wing_param.set_AoA_id(self.__AoA_id)
        self.__wing_param.set_common_OC(self.__OC)
        self.__wing_param.config_from_dict(self.__config_dict)
        
    def __config_DLLM(self):
        ERROR_MSG=self.ERROR_MSG+'__config_DLLM: '
        input_keys=self.__config_dict.keys()
        type_key = self.__tag+'.DLLM.type'
        stype = self.__config_dict[type_key]
        
        if stype not in self.POS_SOLVER:
            raise Exception,ERROR_MSG+'solver_type = '+str(type)+' not in '+str(self.POS_SOLVER)+'. Set to default solver_type = '+str(type)

        if self.__verbose > 0:
            print 'Solver type = ',stype            
        if   stype == 'Solver':
            self.__DLLM_solver = DLLMSolver(self.__tag,self.__wing_param,self.__OC)          
        elif stype == 'TargetCl':
            self.__DLLM_solver = DLLMTargetCl(self.__tag,self.__wing_param,self.__OC)
            target_Cl_key = self.__tag+'.DLLM.target_Cl'
            target_Cl = self.__config_dict[target_Cl_key]
            self.__DLLM_solver.set_target_Cl(target_Cl)
        elif stype == 'TargetLift':
            self.__DLLM_solver = DLLMTargetLift(self.__tag,self.__wing_param,self.__OC)
            target_Lift_key = self.__tag+'.DLLM.target_Lift'
            target_Lift = self.__config_dict[target_Lift_key]
            self.__DLLM_solver.set_target_Lift(target_Lift)
        
        method_key = self.__tag+'.DLLM.method'
        if method_key in input_keys:
            method = self.__config_dict[method_key]
            self.__DLLM_solver.set_method(method)  
            
        relax_factor_key = self.__tag+'.DLLM.relax_factor'
        if relax_factor_key in input_keys:
            relax_factor = self.__config_dict[relax_factor_key]
            self.__DLLM_solver.set_relax_factor(relax_factor)
        
        stop_residual_key = self.__tag+'.DLLM.stop_residual'
        if stop_residual_key in input_keys:
            stop_residual = self.__config_dict[stop_residual_key]
            self.__DLLM_solver.set_stop_residual(stop_residual)
        
        max_iterations_key = self.__tag+'.DLLM.max_iterations'
        if max_iterations_key in input_keys:
            max_iterations = self.__config_dict[max_iterations_key]
            self.__DLLM_solver.set_max_iterations(max_iterations)
        
        gamma_file_name_key = self.__tag+'.DLLM.gamma_file_name'
        if gamma_file_name_key in input_keys:
            gamma_file_name = self.__config_dict[gamma_file_name_key]
            self.__DLLM_solver.set_gamma_file_name(gamma_file_name)
            
        F_list_names_key = self.__tag+'.DLLM.F_list_names'
        if F_list_names_key in input_keys:
            F_list_names = self.__config_dict[F_list_names_key]
            self.__DLLM_solver.set_F_list_names(F_list_names)
                