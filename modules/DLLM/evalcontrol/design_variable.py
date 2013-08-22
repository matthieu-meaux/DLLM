# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus

class DesignVariable:
    """
    Class to describe a design variable
    """
    ERROR_MSG       = 'ERROR DVHandler.'
    WARNING_MSG     = 'WARNING DVHandler.'
    
    def __init__(self,name,lbnd,value,ubnd):
        """
        name: name of the design variable
        lbnd: lower boundary for the design variable
        value: value for the design variable
        ubnd: upper boundary for the design variable
        """
        self.__name   = name
        self.__lbnd   = float(lbnd)
        self.__value  = float(value)
        self.__ubnd   = float(ubnd)
        
        # normanized data
        self.__nlbnd  = 0.
        self.__nvalue = None
        self.__nubnd  = 1.
        
        self.__norm_value()
        
        self.__check()
        
    # accessors
    def get_name(self):
        return self.__name
    
    def get_lbnd(self):
        return self.__lbnd
    
    def get_value(self):
        return self.__value
    
    def get_ubnd(self):
        return self.__ubnd
    
    def get_all(self):
        return self.__lbnd,self.__value,self.__ubnd
    
    def get_all_norm(self):
        return self.__nlbnd,self.__nvalue,self.__nubnd
    
    def get_info(self):
        return self.__name,self.__lbnd,self.__value,self.__ubnd
    
    def get_bounds(self):
        return (self.__lbnd,self.__ubnd)
    
    def get_nvalue(self):
        return self.__nvalue
    
    def get_nbounds(self):
        return (self.__nlbnd,self.__nubnd)

    # set methods
    def set_name(self,name):
        self.__name  = name
    
    def set_lbnd(self,lbnd):
        self.__lbnd  = float(lbnd)
        self.__norm_value()
        
    def set_value(self,value):
        self.__value = float(value)
        self.__norm_value()
        
    def set_ubnd(self,ubnd):
        self.__ubnd  = float(ubnd)
        self.__norm_value()
        
    def set_all(self,lbnd,value,ubnd):
        self.__lbnd  = float(lbnd)
        self.__value = float(value)
        self.__ubnd  = float(ubnd)
        self.__norm_value()
        
    def set_nvalue(self,nvalue):
        self.__nvalue = float(nvalue)
        self.__revert_value()
        
    # private methods
    def __norm_value(self):
        self.__nvalue = (self.__value - self.__lbnd) / ( self.__ubnd - self.__lbnd)
        
    def __revert_value(self):
        self.__value = self.__lbnd + self.__nvalue * (self.__ubnd - self.__lbnd)
        
    def __check(self): 
        ERROR_MSG=self.ERROR_MSG+'.__check: '
        WARNING_MSG=self.WARNING_MSG+'.__check: '
        if self.__ubnd < self.__lbnd:
            raise Exception,ERROR_MSG+'lower bound must be lower than upper bound!'
        if self.__value < self.__lbnd:
            print WARNING_MSG+'current value is lower than lower bound!'
        if self.__value > self.__ubnd:
            print WARNING_MSG+'current value is upper than upper bound!'
        

