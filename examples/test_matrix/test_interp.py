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
# @author : Matthieu MEAUX
#
import numpy
N = 8

gamma=numpy.array([1.,2.,3.,4.,4.,3.,2.,1.])

Mat = numpy.zeros([N+1,N])
Mat[0:N,:]   = numpy.diag(numpy.ones([N]))
Mat[N,:]     = 0.0
Mat[1:N+1,:]-= numpy.diag(numpy.ones([N]))

dgamma=numpy.dot(Mat,gamma)

print dgamma

file=open('interp.dat','w')
for i in xrange(N):
    file.write(str(i)+'\t'+str(dgamma[i])+'\n')
file.close()

Mat2 = numpy.zeros([N+1,N])
Mat2[0:N,:]   = numpy.diag(numpy.ones([N]))
Mat2[N,:]     = 0.0
Mat2[1:N+1,:]-= numpy.diag(numpy.ones([N])) 

Mat2[0,0] = -2
Mat2[0,1] = 3
Mat2[0,2] = -1 

Mat2[-1,-1] = 2
Mat2[-1,-2] = -3
Mat2[-1,-1] = 1 

dgamma=numpy.dot(Mat2,gamma)

file=open('interp2.dat','w')
for i in xrange(N):
    file.write(str(i)+'\t'+str(dgamma[i])+'\n')
file.close()