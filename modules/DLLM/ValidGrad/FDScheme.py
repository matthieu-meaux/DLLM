# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @author Francois Gallard

from numpy import shape

class FDScheme():
    order=None
    
    def __init__(self,fd_step):
        self.__fd_step=fd_step
        self.__samples=[]
        self.__n=0
        
    def get_order(self):
        return self.order
    
    def set_x(self,x):
        self.__x=x
        self.__n=shape(x)[0]
        
    def get_x(self):
        return self.__x
    
    def generate_samples(self):
        pass
    
    def get_fd_step(self):
        return self.__fd_step
    
    def get_samples(self):
        return self.__samples
    
    def set_samples(self,samples):
        self.__samples=samples
        
    def compute_grad(self,y_array):
        return None
    
    def get_grad_dim(self):
        return self.__n
    
