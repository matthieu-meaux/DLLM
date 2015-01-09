# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus Group Innovations

import string
import traceback
import numpy as np

from DLLM.Controllers.Kernel.CManager import CManager
from DLLM.Controllers.Kernel.Variable import Variable
from DLLM.Controllers.Kernel.DesignVariable import DesignVariable
from DLLM.Controllers.Kernel.Parameter import Parameter

class BCManager(CManager):
    
    ERROR_MSG= 'ERROR BCManager.'
    
    def __init__(self, FatherObj=None, gradient_active = True):
        CManager.__init__(self, FatherObj, gradient_active=gradient_active)
        # Design variables handling
        self.__DV_id = []
        self.__DV_pt = []
        self.__ndv = 0
    
    #-- Update methods
    def update(self):
        # Get new list to update in case of dynamic item addition
        list_to_update=self.__get_list_to_update()
        if len(list_to_update) == 0:
            return
        else:
            for pt in list_to_update:
                if pt.is_to_update():
                    pt.update()
            self.update()
                
    def __get_list_to_update(self):
        list_to_update=[]
        for pt in self.get_list_pt():
            if not pt.is_updated():
                list_to_update.append(pt)
        return list_to_update

    #-- Methods 
    def get_print_header(self):
        info_string ='\n-o0 BCManager Information 0o-'
        return info_string
    
    def get_print_footer(self):
        info_string = '\n----------------------------------------------'
        info_string +='\n-o0 End of BCManager Information 0o-'
        return info_string
    
    def add(self,BCpt, check=True):
        """
        add a base controller
        """
        ERROR_MSG=self.ERROR_MSG+'.add: '
        CManager.add(self,BCpt, check=check)
        if BCpt.get_BCType() == 'DesignVariable':
            if BCpt.get_id() not in self.__DV_id:
                self.__DV_id.append(BCpt.get_id())
                self.__DV_pt.append(BCpt)
                self.__ndv+=1
                self.flag_update_all_grad_arrays()#Gradients dimension update
            else:
                raise Exception,ERROR_MSG+'cannot add design variable ID='+str(BCpt.get_id())+', design variable already exists!'
        
    def delete(self,Id):
        """
        delete a base controller
        """
        pointer=self.get_pt(Id)
        CManager.delete(self,Id)
        if pointer.get_BCType() == 'DesignVariable':
            index=self.__DV_id.index(Id)
            del self.__DV_id[index]
            del self.__DV_pt[index]
            self.__ndv-=1
            self.flag_update_all_grad_arrays()#Gradients dimension update
        
    def get_ndv(self):
        return self.__ndv
    
    def get_dv_pt(self,Id):
        """
        return pointer to element from base element id
        """
        ERROR_MSG=self.ERROR_MSG+'get_dv_pt: '
        if Id in self.__DV_id:
            index=self.__DV_id.index(Id)
        else:
            raise Exception,ERROR_MSG+'cannot find desing variable ID='+str(Id)
        return self.__DV_pt[index]
    
    def get_dv_pt_list(self):
        return self.__DV_pt
    
    def get_dv_id_list(self):
        return self.__DV_id
    
    def get_dv_index(self,Id):
        return self.__DV_id.index(Id)
    
    def get_dv_array(self,filter=None):
        dv_list=[]
        if filter is None:
            for pt in self.__DV_pt:
                dv_list.append(pt.get_value())
        else:
            for pt in self.__DV_pt:
                if pt.get_id() in filter:
                    dv_list.append(pt.get_value())
        return np.array(dv_list)
    
    def get_normalized_dv_array(self):
        dv_list=[]
        for pt in self.__DV_pt:
            dv_list.append(pt.get_normalized_value())
        return np.array(dv_list)
    
    def get_bounds_array(self,filter=None):
        bounds_list=[]
        if filter is None:
            for pt in self.__DV_pt:
                bounds_list.append(pt.get_bounds())
        else:
            for pt in self.__DV_pt:
                if pt.get_id() in filter:
                    bounds_list.append(pt.get_bounds())
        return bounds_list
    
    def get_variable_gradient(self):
        return self.__VARIABLE_GRADIENT
    
    # Object creation
    def create_variable(self,Id,value=0.0):
        try:
            return Variable(self,Id,value=value)
        except Exception,error:
            self.delete(Id)
            traceback.print_exc()
            raise Exception,error
        
        
    def create_design_variable(self,Id,value,bounds):
        try:
            return DesignVariable(self,Id,value,bounds)
        except Exception,error:
            self.delete(Id)
            traceback.print_exc()
            raise Exception,error
    
    def create_parameter(self,Id,fexpr,AliasDict=None):
        try:
            return Parameter(self,Id,fexpr,AliasDict=AliasDict)
        except Exception,error:
            self.delete(Id)
            traceback.print_exc()
            raise Exception,str(error)+' in '+Id
        
    # I/0 sections
    def display_design_variables(self):
        info_string ='\n-o0 Design variables information 0o-'
        info_string+='\n   Number of Design Variables : '+str(self.get_ndv())
        info_string+="\n   %20s|"%'VarId  '+"%24s|"%'Value  '+"%24s|"%'Bounds  '
        info_string+="\n      ---------------------------------------------------------------------------------------------"
        for i,varid in enumerate(self.__DV_id):
            pt=self.__DV_pt[i]
            value=pt.get_value()
            bounds=pt.get_bounds()
            info_string+="\n   %20s|"%varid+"%24.16e|"%value+"%24.16e|"%bounds
            
        print info_string
        
    def import_from_file(self,filename):
        """
        Create base controllers from a file
        """
        ERROR_MSG=self.ERROR_MSG+'import_from_file: '
        #- Open file and read all lines
        fid=open(filename,'r')
        all_lines=fid.readlines()
        fid.close()
        
        for line in all_lines:
            words  = string.split(line)
            nwords = len(words)
            if nwords > 0:
                if words[0][0] == '#':
                    pass
                else:
                    if   words[0] == 'Variable':
                        if nwords != 3:
                            raise Exception,ERROR_MSG+'invalid format for Variable!'
                        else:
                            Id    = words[1]
                            value = float(words[2])
                            self.create_variable(Id,value=value)
                    
                    elif words[0] == 'DesignVariable':
                        Id    = words[1]
                        value = float(words[2])
                        bounds  = eval(string.join(words[3:],' '))
                        self.create_design_variable(Id,value,bounds)
                            
                    elif words[0] == 'Parameter':
                        if nwords != 3:
                            raise Exception,ERROR_MSG+'invalid format for Parameter: '+str(words[1])+'!'
                        else:
                            Id    = words[1]
                            fexpr = words[2]
                            self.create_parameter(Id, fexpr, AliasDict=None)
        
        # De-allocate all_lines
        del all_lines
        
    def get_variables_id_list(self):
        """
        Returns the list of ids of variables in BC manager
        """
        var_list=[]
        for controller in self.get_list_pt():
            if controller.get_BCType() == 'Variable':
                var_list.append(controller.get_id())
        return var_list
        
    def export_to_file(self,filename):
        """
        Export BaseControllers to a file
        """
        # lists to store controllers by type
        variable_list = []
        design_variable_list = []
        parameter_list = []
        # loop over all controllers
        for controller in self.get_list_pt():
            if controller.get_BCType() == 'Variable':
                variable_list.append(controller)
            elif controller.get_BCType() == 'DesignVariable':
                design_variable_list.append(controller)
            elif controller.get_BCType() == 'Parameter':
                parameter_list.append(controller)
            else:
                print 'WARNING : BCManager export_to_file: ignoring controller export: '+str(controller.get_id())+' of BCType '+str(controller.get_BCType())
            
        #- Open file in write mode
        fid = open(filename,'w')
        if len(variable_list)>0:
            # Write variables part
            fid.write("# -----------------\n")
            fid.write("# Variables section\n")
            fid.write("# -----------------\n") 
            fid.write('#%23s'%'Type'+' %24s'%'Id'+' %24s'%'value'+'\n')
            for var_pt in variable_list:
                fid.write('%24s'%'Variable'+' %24s'%var_pt.get_id()+' %24.16e'%var_pt.get_value()+'\n')
        
        # Write design variable part
        if len(design_variable_list)>0:
            fid.write("# -----------------\n")
            fid.write("# Design Variables section\n")
            fid.write("# -----------------\n")
            fid.write('#%23s'%'Type'+' %24s'%'Id'+' %24s'%'value'+' %24s'%'bounds'+'\n')
            for design_var_pt in design_variable_list:
                bounds = design_var_pt.get_bounds()
                value = design_var_pt.get_value()
                fid.write('%24s'%'DesignVariable'+' %24s'%design_var_pt.get_id()+' %24.16e'%value+' %24s'%str(bounds)+'\n')
            
        # Write parameter part
        if len(parameter_list)>0:
            fid.write("# -----------------\n")
            fid.write("# Parameters section\n")
            fid.write("# -----------------\n")
            fid.write('#%23s'%'Type'+' %24s'%'Id'+' %24s'%'expression'+'\n')
            for parameter_pt in parameter_list:
                fid.write('%24s'%'Parameter'+' %24s'%parameter_pt.get_id()+' %24s'%parameter_pt.get_expr()+'\n')
        # close file    
        fid.close()
      
    def update_from_file(self,filename):
        """
        update base controller from file
        """
        ERROR_MSG=self.ERROR_MSG+'update_from_file: '
        #- Open file and store all lines in memory
        fid=open(filename,'r')
        all_lines=fid.readlines()
        fid.close()
        
        for line in all_lines:
            words  = string.split(line)
            nwords = len(words)
            if nwords > 0:
                if words[0][0] == '#':
                    pass
                else:
                    if   words[0] == 'Variable':
                        if nwords != 3:
                            raise Exception,ERROR_MSG+'invalid format for Variable!'
                        else:
                            Id    = words[1]
                            value = float(words[2])
                            BC=self.get_pt(Id)
                            BC.update_value(value)
                    
                    elif words[0] == 'DesignVariable':
                        Id    = words[1]
                        value = float(words[2])
                        bounds = eval(string.join(words[3,:]))
                        BC=self.get_pt(Id)
                        BC.set_bounds(bounds)
                        BC.update_value(value)
                            
                    elif words[0] == 'Parameter':
                        if nwords != 3:
                            raise Exception,ERROR_MSG+'invalid format for Parameter!'
                        else:
                            Id    = words[1]
                            fexpr = words[2]
                            BC=self.get_pt(Id)
                            BC.set_fexpr(fexpr)
                            
        self.get_CAD_model().update()
        #- De-allocate all_lines
        del all_lines
        
    def update_from_x(self, x):
        for i, pt in enumerate(self.get_dv_pt_list()):
            pt.set_value(x[i])
        
        self.update()
         
    def update_dv_from_x_list(self,x,filter=None):
        """
        Update design variable vector from x list
        @param filter : a list of Ids of variables to update. If the design variable is not in the filter, then it is not updated
        """
        ERROR_MSG=self.ERROR_MSG+'update_dv_from_x_list: '
        if len(x) != self.get_ndv() and filter is None:
            raise Exception,ERROR_MSG+'provided x vector does not have the required dimension'
        if filter is None:
            for Id,val in enumerate(x):
                self.__DV_pt[Id].update_value(val)
        else:
            for val,filter_id in zip(x,filter):
                self.get_dv_pt(filter_id).update_value(val)
        
        self.update()            
        
    def update_dv_from_normalized_x_list(self,x_norm):
        """
        Update design variable vector from normalized x list
        """
        ERROR_MSG=self.ERROR_MSG+'update_dv_from_normalized_x_list: '
        if len(x_norm) != self.get_ndv():
            raise Exception,ERROR_MSG+'provided x_norm vector does not have the required dimension'
        for Id,val in enumerate(x_norm):
            self.__DV_pt[Id].update_nomalized_value(val)
        
        self.update()
        
    def return_normalized_gradient(self,input_grad):
        output_grad=input_grad[:]
        for Id,pt in enumerate(self.__DV_pt):
            fact=pt.get_revert_fact()
            output_grad[Id]=input_grad[Id]*fact
        return output_grad
        
    def update_dv_from_file(self,filename):
        """
        Update design variables from a design variables file
        """
        ERROR_MSG=self.ERROR_MSG+'update_dv_from_file: '
        #- Open file and store all lines in memory
        fid=open(filename,'r')
        all_lines=fid.readlines()
        fid.close()
        
        for line in all_lines:
            words  = string.split(line)
            nwords = len(words)
            if nwords > 0:
                if words[0][0] == '#':
                    pass
                else:
                    if nwords == 4:
                            Id    = words[0]
                            lbnd  = float(words[1])
                            value = float(words[2])
                            ubnd  = float(words[3])
                            BC=self.get_pt(Id)
                            BC.set_bounds(lbnd,ubnd)
                            BC.update_value(value)
                    else:
                        raise Exception,ERROR_MSG+'unexpected file format :'+str(filename)+', please check line : '+str(line)
                            
        self.update()
        #- De-allocate all_lines
        del all_lines
            
    def export_dv_to_file(self,filename):
        """
        Export design variable to a file
        """
        #- File object creation
        all_lines=[]
        all_lines.append('# Design variables file generate from PADGE \n')
        for pt in self.__DV_pt:
            Id=pt.get_id()
            (lbnd,ubnd)=pt.get_bounds()
            value=pt.get_value()
            all_lines.append("%20s"%Id+" %24.16e"%lbnd+" %24.16e"%value+" %24.16e"%ubnd+"\n")
            
        fid=open(filename,'w')
        fid.writelines(all_lines)
        fid.close()
        
    def convert_to_variable(self,Id):
        ERROR_MSG=self.ERROR_MSG+'convert_to_variable: '
        try:
            pointer=self.get_pt(Id)
        except:
            raise Exception,ERROR_MSG+'Base controller id = '+str(Id)+' cannot be found, conversion failed!'
        
        BC_type = pointer.get_BCType()
        value   = pointer.get_value()
        if BC_type != 'Variable':
            influences_list=[]
            for pt in pointer.get_influences():
                influences_list.append(pt)
            for pt in influences_list:
                pt.del_dependance(pointer)
                    
            self.delete(Id)
            
            new_pointer=self.create_variable(Id, value)
            
            for pt in influences_list:
                pt.special_dependance_update(new_pointer)
            new_pointer.set_to_update()
            new_pointer.set_controllers_to_update()
            
            del influences_list
            
            
    def convert_to_design_variable(self,Id,bounds):
        ERROR_MSG=self.ERROR_MSG+'convert_to_design_variable: '
        try:
            pointer=self.get_pt(Id)
        except:
            raise Exception,ERROR_MSG+'Base controller id = '+str(Id)+' cannot be found, conversion failed!'
        
        BC_type = pointer.get_BCType()
        value   = pointer.get_value()
        if BC_type != 'DesignVariable':
            influences_list=[]
            for pt in pointer.get_influences():
                influences_list.append(pt)
            for pt in influences_list:
                pt.del_dependance(pointer)
                
            self.delete(Id)
            
            new_pointer=self.create_design_variable(Id,value,bounds)
            
            for pt in influences_list:
                pt.special_dependance_update(new_pointer)
            new_pointer.set_to_update()
            new_pointer.set_controllers_to_update()
            
            del influences_list
    
    def convert_to_parameter(self,Id,fexpr):
        ERROR_MSG=self.ERROR_MSG+'convert_to_parameter: '
        try:
            pointer=self.get_pt(Id)
        except:
            raise Exception,ERROR_MSG+'Base controller id = '+str(Id)+' cannot be found, conversion failed!'
        
        BC_type = pointer.get_BCType()
        if BC_type != 'Parameter':
            influences_list=[]
            for pt in pointer.get_influences():
                influences_list.append(pt)
            for pt in influences_list:
                pt.del_dependance(pointer)
                
            self.delete(Id)
            
            new_pointer=self.create_parameter(Id, fexpr)

            for pt in influences_list:
                pt.special_dependance_update(new_pointer)
            new_pointer.set_to_update()
            new_pointer.set_controllers_to_update()
            
            del influences_list
            
    def update_dv_from_info_list(self, info_list):
        for dv_info in info_list:
            id=dv_info[0]
            lbnd=dv_info[1]
            value=dv_info[2]
            ubnd=dv_info[3]
            local_dv=self.get_dv_pt(id)
            local_dv.set_bounds(lbnd,ubnd)
            local_dv.set_value(value)
            
    def get_dv_info_list(self):
        dv_info_list=[]
        for dv in self.get_dv_pt_list():
            id=dv.get_id()
            value=dv.get_value()
            bounds=dv.get_bounds()
            dv_info_list.append([id,bounds[0],value,bounds[1]])
        return dv_info_list
            
    def get_non_influent_DV_list(self):
        uninfl_dv=[]
        for id,pt in zip(self.__DV_id,self.__DV_pt):
            if not pt.is_influent():
                uninfl_dv.append(id)
                
        return uninfl_dv
    
    def get_tags_x0_bounds(self):
        tags=[]
        x0=[]
        bounds_list=[]
        for dv in self.get_dv_pt_list():
            tag = dv.get_id()
            value  = dv.get_value()
            bounds = dv.get_bounds()
            tags.append(tag)
            x0.append(value)
            bounds_list.append(bounds)
        return tags, x0, bounds_list
        
            
        
        