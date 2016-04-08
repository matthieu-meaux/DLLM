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
from DLLM.DLLMEval.DLLMWrapper import DLLMWrapper
from copy import deepcopy
import multiprocessing 
import numpy
import string
import sys
import os

class DLLMMP():
    ERROR_MSG = 'ERROR in DLLMMP.'
    WARNING_MSG = 'WARNING in DLLMMP.'
    POS_FMT = ['list','numpy']
    MP_KEYS = ['nb_conditions','condition_name']
    
    #-- Constructor
    def __init__(self, tag, verbose = 0):
        """
        Multi-point class for the DLLMSolver
        """
        self.__tag            = tag
        self.__verbose        = verbose
        
        self.__config_dict    = None
        self.__exclude_keys   = None
        self.__list_conf_dict = None
        self.__wrapper_list   = None
        self.__nb_cond        = None
        self.__cond_name      = None
        
        self.__out_format     = 'list'
        self.__grad_format    = 'list'
        
        self.__AoA_id_list    = None
        
    #-- Accessors
    def get_tags_x0_and_bounds(self):
        tags=[]
        Wrap0=self.__wrapper_list[0]
        wrap_tags,x0,bounds=Wrap0.get_tags_x0_and_bounds()
        for tag in wrap_tags:
            words=tag.split('.')
            if len(words)>1:
                next_index=self.__check_tag(words)
                new_tag=string.join(words[next_index+1:],'.')
            else:
                new_tag=words[0]
            tags.append(new_tag)
        return tags,x0,bounds
    
    def get_x0_and_bounds(self):
        tags,x0,bounds=self.get_tags_x0_and_bounds()
        return x0,bounds
    
    def get_x0(self):
        return self.get_x()
    
    def get_x(self):
        tags,x0,bounds=self.get_tags_x0_and_bounds()
        return x0
    
    def get_F_list_names(self):
        MP_F_list_names=[]
        for i in xrange(self.__nb_cond):
            cond_name = self.__cond_name+str(i+1)
            F_list_names=self.__wrapper_list[i].get_F_list_names()
            print 'MP F_list_names',F_list_names
            for name in F_list_names:
                MP_F_list_names.append(cond_name+'_'+name)
        return MP_F_list_names
                
    #-- Setters
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
        self.__init_configure()
        self.__set_list_config_dict()
        self.__configure_wrapper_list()
        
    def run(self, x):
        process_list=[]
        # Launch all sub-processes 
        for i in xrange(self.__nb_cond):
            process_list.append(multiprocessing.Process(target=self.__wrap_run, args=(i,x)))
            process_list[-1].start()
        # wait for all processes
        for i in xrange(self.__nb_cond):
            process_list[i].join()
        # import results
        self.__import_results()
        # Gather information
        MP_F_list=self.__build_MP_F_list()
        return MP_F_list
    
    def __wrap_run(self, i, x):
        std_sys=sys.stdout
        cond_name = self.__cond_name+str(i+1)
        sys.stdout = open(cond_name+'.log','a')
        print 'PID=',os.getpid()
        self.__wrapper_list[i].run(x)
        self.__wrapper_list[i].export_results()
        sys.stdout.close()
        sys.stdout=std_sys
    
    def run_grad(self, x):
        process_list=[]
        # Launch all sub-processes 
        for i in xrange(self.__nb_cond):
            process_list.append(multiprocessing.Process(target=self.__wrap_run_grad, args=(i,x)))
            process_list[-1].start()
        # wait for all processes
        for i in xrange(self.__nb_cond):
            process_list[i].join()
        # import results
        self.__import_results()
        # Gather information
        MP_F_list_grad=self.__build_MP_F_list_grad()
        return MP_F_list_grad
    
    def __wrap_run_grad(self, i, x):
        std_sys=sys.stdout
        cond_name = self.__cond_name+str(i+1)
        sys.stdout = open(cond_name+'.log','a')
        print 'PID=',os.getpid()
        self.__wrapper_list[i].run_grad(x)
        self.__wrapper_list[i].export_results()
        sys.stdout.close()
        sys.stdout=std_sys
    
    def run_and_grad(self, x):
        process_list=[]
        # Launch all sub-processes 
        for i in xrange(self.__nb_cond):
            process_list.append(multiprocessing.Process(target=self.__wrap_run_and_grad, args=(i,x)))
            process_list[-1].start()
        # wait for all processes
        for i in xrange(self.__nb_cond):
            process_list[i].join()
        # import results
        self.__import_results()
        # Gather information
        MP_F_list,MP_F_list_grad=self.__build_MP_F_list_and_grad()
        return MP_F_list, MP_F_list_grad
    
    def __wrap_run_and_grad(self, i, x):
        std_sys=sys.stdout
        cond_name = self.__cond_name+str(i+1)
        sys.stdout = open(cond_name+'.log','a')
        print 'PID=',os.getpid()
        self.__wrapper_list[i].run_and_grad(x)
        self.__wrapper_list[i].export_results()
        sys.stdout.close()
        sys.stdout=std_sys

    def analysis(self):
        process_list=[]
        # Launch all sub-processes 
        for i in xrange(self.__nb_cond):
            process_list.append(multiprocessing.Process(target=self.__wrap_analysis, args=(i,)))
            process_list[-1].start()
        # wait for all processes
        for i in xrange(self.__nb_cond):
            process_list[i].join()
        # import results
        self.__import_results()
        # Gather information
        MP_F_list=self.__build_MP_F_list()
        return MP_F_list
    
    def __wrap_analysis(self, i):
        std_sys=sys.stdout
        cond_name = self.__cond_name+str(i+1)
        sys.stdout = open(cond_name+'.log','a')
        print 'PID=',os.getpid()
        self.__wrapper_list[i].analysis()
        self.__wrapper_list[i].export_results()
        sys.stdout.close()
        sys.stdout=std_sys
    
    def analysis_grad(self):
        process_list=[]
        # Launch all sub-processes 
        for i in xrange(self.__nb_cond):
            process_list.append(multiprocessing.Process(target=self.__wrap_analysis_grad, args=(i,)))
            process_list[-1].start()
        # wait for all processes
        for i in xrange(self.__nb_cond):
            process_list[i].join()
        # import results
        self.__import_results()
        # Gather information
        MP_F_list_grad=self.__build_MP_F_list_grad()
        return MP_F_list_grad
    
    def __wrap_analysis_grad(self, i):
        std_sys=sys.stdout
        cond_name = self.__cond_name+str(i+1)
        sys.stdout = open(cond_name+'.log','a')
        print 'PID=',os.getpid()
        self.__wrapper_list[i].analysis_grad()
        self.__wrapper_list[i].export_results()
        sys.stdout.close()
        sys.stdout=std_sys
    
    def analysis_and_grad(self):
        process_list=[]
        # Launch all sub-processes 
        for i in xrange(self.__nb_cond):
            process_list.append(multiprocessing.Process(target=self.__wrap_analysis_and_grad, args=(i,)))
            process_list[-1].start()
        # wait for all processes
        for i in xrange(self.__nb_cond):
            process_list[i].join()
        # import results
        self.__import_results()
        # Gather information
        MP_F_list,MP_F_list_grad=self.__build_MP_F_list_and_grad()
        return MP_F_list, MP_F_list_grad
    
    def __wrap_analysis_and_grad(self, i):
        std_sys=sys.stdout
        cond_name = self.__cond_name+str(i+1)
        sys.stdout = open(cond_name+'.log','a')
        print 'PID=',os.getpid()
        self.__wrapper_list[i].analysis_and_grad()
        self.__wrapper_list[i].export_results()
        sys.stdout.close()
        sys.stdout=std_sys
        
    #-- Private methods
    def __import_results(self):
        for i in xrange(self.__nb_cond):
            self.__wrapper_list[i].import_results()
            
    def __build_MP_F_list(self):
        MP_F_list=[]
        for i in xrange(self.__nb_cond):
            F_list=self.__wrapper_list[i].get_F_list()
            MP_F_list+=F_list
        if self.__out_format == 'numpy':
            MP_F_list=numpy.array(MP_F_list)
        return MP_F_list
    
    def __build_MP_F_list_grad(self):
        MP_F_list_grad=[]
        for i in xrange(self.__nb_cond):
            F_list_grad=self.__wrapper_list[i].get_F_list_grad()
            MP_F_list_grad+=F_list_grad
        if self.__grad_format == 'numpy':
            MP_F_list_grad=numpy.array(MP_F_list_grad)
        return MP_F_list_grad
    
    def __build_MP_F_list_and_grad(self):
        MP_F_list=[]
        MP_F_list_grad=[]
        for i in xrange(self.__nb_cond):
            F_list,F_list_grad=self.__wrapper_list[i].get_F_list_and_grad()
            MP_F_list+=F_list
            MP_F_list_grad+=F_list_grad
        if self.__out_format == 'numpy':
            MP_F_list=numpy.array(MP_F_list)
        if self.__grad_format == 'numpy':
            MP_F_list_grad=numpy.array(MP_F_list_grad)
        return MP_F_list,MP_F_list_grad
    
    def __configure_wrapper_list(self):
        self.__wrapper_list=[]
        for i in xrange(self.__nb_cond):
            cond_name = self.__cond_name+str(i+1)
            self.__wrapper_list.append(DLLMWrapper(self.__tag+'.'+cond_name, verbose=self.__verbose))
            if self.__AoA_id_list is not None:
                AoA_id = self.__AoA_id_list[i]
                self.__wrapper_list[i].set_AoA_id(AoA_id)
            self.__wrapper_list[i].set_out_format('list')
            self.__wrapper_list[i].set_grad_format('list')
            self.__wrapper_list[i].configure(self.__list_config_dict[i])

            
    def __init_configure(self):
        config_keys=sorted(self.__config_dict.keys())
        self.__nb_cond     = self.__config_dict[self.__tag+'.nb_conditions']
        self.__cond_name   = self.__config_dict[self.__tag+'.condition_name']
        if self.__tag+'.AoA_id_list' in config_keys:
            self.__AoA_id_list = self.__config_dict[self.__tag+'.AoA_id_list']
        if self.__verbose > 0:
            print '*** Multi-conditions information ***'
            print '  Number of conditions = ',self.__nb_cond
            print '  Base condition name  = ',self.__cond_name
            print '  AoA_id list          = ',str(self.__AoA_id_list)
            print '***                              ***'
        
        self.__list_config_dict=[]
        self.__exclude_keys = deepcopy(self.MP_KEYS)
        
        for i in xrange(self.__nb_cond):
            cond_name = self.__cond_name+str(i+1)
            self.__exclude_keys.append(cond_name)
            self.__list_config_dict.append({})
    
    def __set_list_config_dict(self):
        config_keys=sorted(self.__config_dict.keys())
        
        for i in xrange(self.__nb_cond):
            cond_name = self.__cond_name+str(i+1)
            cond_key  = self.__tag+'.'+cond_name
            
            # Setup common keys
            for key in config_keys:
                words=key.split('.')
                next_index=self.__check_tag(words)
                next_word = words[next_index]
                # common keys
                if words[next_index] not in self.__exclude_keys:
                    common_key=string.join(words[next_index:],'.')
                    cond_spec_key = cond_key+'.'+common_key
                    self.__list_config_dict[i][cond_spec_key] = self.__config_dict[key]
                # condition specific key
                if words[next_index] == cond_name:
                    self.__list_config_dict[i][key]=self.__config_dict[key]
    
    def __check_tag(self, list, index=0):
        if string.join(list[0:index+1],'.') == self.__tag:
            next_index = index+1
            return next_index
        else:
            new_index=index+1
            return self.__check_tag(list, new_index)
        
        