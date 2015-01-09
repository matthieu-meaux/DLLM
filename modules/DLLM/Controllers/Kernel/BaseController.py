# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus Group Innovations
# Derived from PBaseController PADGE class
from DLLM.Controllers.Kernel.Controller import Controller

class BaseController(Controller):
    """
    BaseController Class
    """
    #--Class Variables
    ERROR_MSG='ERROR BaseController.'
    CLASSTYPE='BC'
    USE_GRADIENT_ARRAYS=True
    
    #--Constructor
    def __init__(self, BCManager, Id, value=0.0, BCType='Generic'):
        
        self.__BCType   = BCType
        self.__value    = None    #Initilizing attributes before updating the gradient in PController.handle_dv_changes
        self.__gradient = None
        
        Controller.__init__(self,BCManager,self.CLASSTYPE,Id)
        self.set_value(value)
        
    #--Private methods
    def __repr__(self):
        info_string  =Controller.__repr__(self)
        info_string += '\n   BC Type         : '+self.get_BCType()
        return info_string

    #--Accessors
    def get_BCType(self):
        return self.__BCType
    
    def get_value(self):
        return self.__value

    def get_gradient(self):
        return self.__gradient
    
    #--Setters        
    def set_value(self,value,flag_updates=True):
        self.__value = value
        if flag_updates:
            self.set_controllers_to_update()
        
    def update_value(self,value, flag_updates=True):
        if self.__value!=value:
            self.set_value(value,flag_updates=flag_updates)
        
    def set_gradient(self,gradient):
        self.__gradient = gradient
        self.set_controllers_to_update()
