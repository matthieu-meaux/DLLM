# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @author Francois Gallard

from FDScheme import FDScheme
from numpy import zeros, shape, array
from copy import deepcopy

class FDSecondOrderCentered(FDScheme):
    order=2
    
    def __init__(self,fd_step):
        FDScheme.__init__(self,fd_step)
        
    def generate_samples(self):
        x=self.get_x()
        e=self.get_fd_step()
        x_samples=[]
        
        for i in range(self.get_grad_dim()):
            x_c=deepcopy(x)
            x_c[i]+=e
            x_samples.append(x_c)
            
        for i in range(self.get_grad_dim()):
            x_c=deepcopy(x)
            x_c[i]-=e
            x_samples.append(x_c)            
            
        self.set_samples(x_samples)
        
    def compute_grad(self,y_array):
        if type(y_array) == type(zeros(1)):
            n=len(shape(y_array))
            p=self.get_grad_dim()
            if n < 3:
                grad=array(y_array[:p]-y_array[p:]).T/(2.*self.get_fd_step())
            elif n == 3:
                s=shape(y_array)
                grad=zeros((s[0],s[1],p))
                for i in range(p):
                    grad[:,:,i]=(y_array[:,:,i]-y_array[:,:,i+p])/(2.*self.get_fd_step())
            else:
                raise Exception, "Functional outputs of dimension > 2 are not yet handled"
                
            return grad
        else:
            raise Exception, "FDGradient for functional of type : "+str(type(y_array[0]))+" not implemented yet."
