# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus

import string
from OpenDACE.evalcontrol.design_variable import DesignVariable

class DVHandler:
    """
    Class to handle design variable vector
    """
    ERROR_MSG       = 'ERROR DVHandler.'
    WARNING_MSG     = 'WARNING DVHandler.'
    
    def __init__(self):
        """
        Constructor
        """
        self.__DVId   = []
        self.__DVList = []
        
    def __repr__(self):
        info_string ='\n-o0 Design variables handler information 0o-'
        info_string+='\n   Number of design variables : '+str(self.get_nvar())
        info_string+="\n   %24s|"%'VarId  '+"%24s|"%'Lbnd  '+"%24s|"%'Value  '+"%24s|"%'Ubnd  '
        info_string+="\n      ---------------------------------------------------------------------------------------------"
        for DVpt in self.__DVList:
            name,lbnd,value,ubnd=DVpt.get_info()
            info_string+="\n   %24s|"%name+"%24.16e|"%lbnd+"%24.16e|"%value+"%24.16e|"%ubnd
        info_string+="\n      ---------------------------------------------------------------------------------------------"
        return info_string
    
    # accessors
    def get_DVId(self):
        return self.__DVId
        
    def add_DV(self,name,lbnd,value,ubnd):
        """
        Add a design variable to the handler
        """
        ERROR_MSG=self.ERROR_MSG+'add_DV: '
        
        if name not in self.__DVId:
            self.__DVId.append(name)
            self.__DVList.append(DesignVariable(name,lbnd,value,ubnd))
        else:
            raise Exception,ERROR_MSG+'variable named: '+name+' already exist in DVhandler!'
        
    def del_DV(self,name):
        WARNING_MSG = self.WARNING_MSG+'del_DV: '
        
        if name in self.__DVId:
            index=self.__DVId.index(name)
            del self.__DVId[index]
            del self.__DVList[index]
        else:
            warning_msg = WARNING_MSG+'Design variable named: '+name+' not deleted, design variable not found!'
            print warning_msg
            
    def replace_DV(self,name,new_name,lbound,new_value,ubound):
        WARNING_MSG = self.WARNING_MSG+'replace_DV: '
        if name in self.__DVId:
            indx=self.__DVId.index(name)
            self.__DVList[indx].set_name(new_name)
            self.__DVList[indx].set_all(lbound,new_value,ubound)
            self.__DVId[indx]=new_name
        else:
            warning_msg = WARNING_MSG+'Design variable named: '+name+' not replaced, design variable not found!'
            print warning_msg    
        
    def rename_DV(self,name,new_name):
        WARNING_MSG = self.WARNING_MSG+'rename_DV: '
        if name in self.__DVId:
            indx=self.__DVId.index(name)
            self.__DVList[indx].set_name(new_name)
            self.__DVId[indx]=new_name
        else:
            warning_msg = WARNING_MSG+'Design variable named: '+name+' not renamed, design variable not found!'
            print warning_msg
            
    def get_nvar(self):
        return len(self.__DVId)
            
    def get_x(self):
        x=[]
        for DVpt in self.__DVList:
            x.append(DVpt.get_value())
        return x

    def set_x(self, x):
        for DVpt,xi in zip(self.__DVList,x):
            DVpt.set_value(xi)
    
    def get_var(self,name):
        return self.get_variable(name).get_value()

    def get_variable(self,name):
        WARNING_MSG = self.WARNING_MSG+'get_var: '
        if name in self.__DVId:
            id=self.__DVId.index(name)
            return self.__DVList[id]
        else:
            warning_msg = WARNING_MSG+'Design variable named: '+name+' not found!'
            print warning_msg
    
    def set_var(self,name,value):
        WARNING_MSG = self.WARNING_MSG+'set_var: '
        if name in self.__DVId:
            id=self.__DVId.index(name)
            self.__DVList[id].set_value(value)
        else:
            warning_msg = WARNING_MSG+'Design variable named: '+name+' not found!'
            print warning_msg
        
    def is_var_in_handler(self,name):
        return self.__DVId.count(name)>0
    
    def get_x_names(self):
        x_names=[]
        for DVpt in self.__DVList:
            x_names.append(DVpt.get_name())
        return x_names
    
    def get_bounds(self):
        bounds=[]
        for DVpt in self.__DVList:
            bounds.append(DVpt.get_bounds())
        return bounds
    
    def get_var_bounds(self,name):
        return self.get_variable(name).get_nbounds()
    
    def get_x_bounds(self):
        x=self.get_x()
        bounds=self.get_bounds()
        return x,bounds
    
    def get_xid_x_bounds(self):
        xid=self.get_DVId()
        x=self.get_x()
        bounds=self.get_bounds()
        return xid,x,bounds
    
    def get_x_norm(self):
        xnorm=[]
        for DVpt in self.__DVList:
            xnorm.append(DVpt.get_nvalue())
        return xnorm
    
    def get_bounds_norm(self):
        boundsnorm=[]
        for DVpt in self.__DVList:
            boundsnorm.append(DVpt.get_nbounds())
        return boundsnorm

    def get_x_bounds_norm(self):
        xnorm=self.get_x_norm()
        boundsnorm=self.get_bounds_norm()
        return xnorm,boundsnorm


    # update methods
    def update_from_x(self,x):
        ERROR_MSG=self.ERROR_MSG+'update_from_x: '
        if len(x) != len(self.__DVId):
            raise Exception,ERROR_MSG+'provided x = '+str(x)+' does not have the correct dimension!'
        for index,value in enumerate(x):
            self.__DVList[index].set_value(value)

    def update_from_x_norm(self,xnorm):
        ERROR_MSG=self.ERROR_MSG+'update_from_x_norm: '
        if len(xnorm) != len(self.__DVId):
            raise Exception,ERROR_MSG+'provided xnorm = '+str(xnorm)+' does not have the correct dimension!'
        for index,value in enumerate(xnorm):
            self.__DVList[index].set_nvalue(value)
    
    #I/O methods
    def import_from_file(self,filename):
        """
        import design variables from a file
        """
        #- Open file and read all lines
        fid=open(filename,'r')
        all_lines=fid.readlines()
        fid.close()
        for line in all_lines:
            words=string.split(line)
            nwords=len(words)
            if nwords >0:
                if words[0][0] != '#':
                    if nwords !=4:
                        raise Exception,self.ERROR_MSG+'ImportFromFile: Invalid file!'
                    else:
                        name=words[0]
                        lbnd=eval(words[1])
                        value=eval(words[2])
                        ubnd=eval(words[3])
                        self.add_DV(name,lbnd,value,ubnd)
        fid.close()
        
    def export_to_file(self,filename):
        """
        export design variables to a file
        """
        #- Open file and read all lines
        fid=open(filename,'w')
        fid.write('# design variables file (generated for DVHandler)\n')
        fid.write("#%23s"%'Variable ID'+" %24s"%'Lower bound'+" %24s"%'Variable Value'+" %24s"%'Upper bound'+"\n")
        for DVpt in self.__DVList:
            name,lbnd,value,ubnd=DVpt.get_info()
            fid.write("%24s"%name+" %24.16e"%lbnd+" %24.16e"%value+" %24.16e"%ubnd+"\n")
        fid.close()




