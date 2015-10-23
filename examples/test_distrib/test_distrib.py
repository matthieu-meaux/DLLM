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

N=50

r_list_lin = numpy.linspace(0.,1.,N+1) 

theta_list = numpy.linspace(0.,numpy.pi,N+1)
r_list_cos = 0.5*(1.-numpy.cos(theta_list))

fid=open('distrib.dat','w')
for i in xrange(N+1):
    fid.write(str(0.)+' '+str(r_list_lin[i])+' '+str(r_list_cos[i])+'\n')
fid.close()
    

    