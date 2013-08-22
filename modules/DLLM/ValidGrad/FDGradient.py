# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @author Francois Gallard

from numpy import zeros, shape, array
from copy import deepcopy
class FDGradient():
    
    def __init__(self,scheme,f_pointer):
        self.__scheme=scheme
        self.__fpointer=f_pointer
        
    def get_scheme(self):
        return self.__scheme
    
    def grad_f(self,x):
        self.__scheme.set_x(x)
        self.__scheme.generate_samples()
        
        y=[]
        for x in self.__scheme.get_samples():
            y.append(deepcopy(self.__fpointer(x)))
            
        s=shape(y[0])
        if len(s)<2:
            y_array=array(y)
        elif len(s)==2:
            p=len(y)
            y_array=zeros((s[0],s[1],p))
            for i in range(p):
                y_array[:,:,i]=y[i]
        else:
            raise Exception, "Functional outputs of dimension >2 are not yet handled."
        return self.__scheme.compute_grad(y_array)
    
