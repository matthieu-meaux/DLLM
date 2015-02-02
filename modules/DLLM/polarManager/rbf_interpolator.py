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
#  http://github.com/TBD
#
from scipy.interpolate import Rbf
from AeroElastAdj.ValidGrad.FDGradient import FDGradient
from AeroElastAdj.ValidGrad.FDSecondOrderCentered import FDSecondOrderCentered

class RBF_interpolator():
    EPSILON=1e-7
    
    def __init__(self,x,y,scale,epsilon=0.1):
        self.__x=x
        self.__y=y
        self.__scale=scale
        args=self.__get_args(x)
        args.append(y.T)
        fd2=FDSecondOrderCentered(self.EPSILON)
        self.__fd_grad=FDGradient(fd2,self.f)
        self.__rbfi=Rbf(*args,epsilon=epsilon,function="linear")

    def __get_args(self,x):
        args=[]
        for i in xrange(x.shape[1]):
            args.append(x[:,i]*self.__scale[i])
        return args
        
    def f(self,x):
        args=x*self.__scale
        return self.__rbfi(*args)
    
    def df(self,x,comp=0):
        return self.__fd_grad.grad_f(x)[comp]
    