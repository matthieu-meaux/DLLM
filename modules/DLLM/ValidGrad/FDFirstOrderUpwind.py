# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @author Francois Gallard

from FDScheme import FDScheme

from numpy import zeros, shape,dot ,eye
from copy import deepcopy

class FDFirstOrderUpwind(FDScheme):
    order=1
    
    def __init__(self,fd_step):
        FDScheme.__init__(self,fd_step)
        
    def generate_samples(self):
        x=self.get_x()
        x_samples=[deepcopy(x)]
        e=self.get_fd_step()
        
        for i in range(shape(x)[0]):
            x_c=deepcopy(x)
            x_c[i]+=e
            x_samples.append(x_c)
            
        self.set_samples(x_samples)
        
    def compute_grad(self,y_array):
        if type(y_array) == type(zeros(1)):
            n=len(shape(y_array))
            if n < 3:
                grad=(y_array[1:]-y_array[0]).T/self.get_fd_step()
            elif n == 3:
                s=shape(y_array)
                p=self.get_grad_dim()
                grad=zeros((s[0],s[1],p))
                for i in range(p):
                    grad[:,:,i]=(y_array[:,:,i+1]-y_array[:,:,0])/self.get_fd_step()
            else:
                raise Exception, "Functional outputs of dimension > 2 are not yet handled"
                
            return grad
        else:
            raise Exception, "FDGradient for functional of type : "+str(type(y_array))+" not implemented yet."