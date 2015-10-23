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
span=40.
gamma0=1.70524
N=50

theta_list = numpy.linspace(0.,numpy.pi,N+1)
r_list_eta = -0.5+0.5*(1.-numpy.cos(theta_list))

fid=open('ideal_gamma.dat','w')
for i,r in enumerate(r_list_eta):
    y=r*span
    gamma=numpy.sqrt(1-(2*r)**2)*gamma0
    fid.write(str(y)+' '+str(gamma)+'\n')
fid.close()
    
    