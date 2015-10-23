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
from numpy import array, transpose,outer, ones, zeros, copy, divide, diag, dot,log
from numpy.linalg import norm, solve

N = 4

eta=numpy.array([-2.,-1.,0.,1.,2.])
y=numpy.array([-1.5,-0.5,0.5,1.5])

# print 'y=',y
# print 'eta=',eta
#                 
# YminEta=transpose(outer(ones([N+1]),y))-outer(ones([N]),eta)
# print "outer(ones([N+1]),y)",outer(ones([N+1]),y)
# print 'YminEta=',YminEta
# 
# Kmetric=divide(ones([N,N+1]),YminEta)
# print "divide(ones([N,N+1]),YminEta)",divide(ones([N,N+1]),YminEta)
# Kmetric/=4.*numpy.pi
# 
# DdGammaDy_DGamma = zeros([N+1,N])
# DdGammaDy_DGamma[0:N,:]   = diag(ones([N]))
# DdGammaDy_DGamma[N,:]     = 0.0
# DdGammaDy_DGamma[1:N+1,:]-= diag(ones([N]))
# print "DdGammaDy_DGamma",DdGammaDy_DGamma
# 
# K = - dot(Kmetric,DdGammaDy_DGamma)
# 
# print "K",K

d_eta=eta[1:]-eta[:-1]

print 'd_eta=',d_eta

log_mat=zeros((N,N))
for i in xrange(N):
    for j in xrange(N):
        log_mat[i,j]=log(abs((eta[j+1]-y[i])/(eta[j]-y[i])))

DdGammaDy_DGamma = zeros([N,N])
for i in xrange(N):
        indexp=i+1
        indexm=i-1
        if indexp<=N-1:
            DdGammaDy_DGamma[i,indexp]=-1
        if indexm>=0:
            DdGammaDy_DGamma[i,indexm]=1.
        print indexp,indexm
        
    
print 'log_mat=',log_mat
print "DdGammaDy_DGamma",DdGammaDy_DGamma  
