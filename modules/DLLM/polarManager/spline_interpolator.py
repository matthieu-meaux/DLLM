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
#  https://github.com/matthieu-meaux/DLLM.git
#
# @author : Francois Gallard
#

from scipy.interpolate import bisplrep,bisplev

class Spline2D_interpolator():
    
    def __init__(self,x,y,scale):
        self.__x=x
        self.__y=y
        self.__scale=scale
        self.__tck = bisplrep(x=self.__x[:,0]*scale[0],y=self.__x[:,1]*scale[1],z=self.__y, kx=5, ky=5)
        
    def f(self,x):
        return bisplev(x[0]*self.__scale[0],x[1]*self.__scale[1],self.__tck)
    
    def df(self,x,comp=0):
        dx=1-comp
        dy=comp
        grad=bisplev(x[0]*self.__scale[0],x[1]*self.__scale[1],self.__tck,dx=dx,dy=dy)
        return grad*self.__scale[comp]

