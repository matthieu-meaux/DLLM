# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus Group Innovations

from DLLM.Controllers.Kernel.BaseController import BaseController

from numpy import zeros

class Variable(BaseController):
    """
    PVariable Class
    """
    #--Class variables
    CLASS_MSG = 'PVariable'
    ERROR_MSG = 'ERROR '+CLASS_MSG+'.'
    GRAD = None
    
    #--constructor
    def __init__(self,PBCManager,ID,value=0.0):
        BaseController.__init__(self, PBCManager, ID, value=value, BCType='Variable')
        
    #--Private methods
    def __repr__(self):
        """
        display some information about the variable
        """
        info_string=BaseController.__repr__(self)
        info_string += '\n   Value           :%24.16e'%self.get_value()
        info_string += '\n   Gradient        : '+str(self.get_gradient())
        return info_string

    #--Methods
    def handle_dv_changes(self):
        if self.__class__.GRAD is None or self.get_ndv() != len(self.__class__.GRAD):
            self.__class__.GRAD= zeros(self.get_ndv())
            
        self.set_gradient(self.__class__.GRAD)
        