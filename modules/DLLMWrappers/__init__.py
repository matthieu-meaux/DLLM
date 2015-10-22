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

"""
B{DLLM}: B{D}ifferentiated B{L}ifting B{L}ine B{M}odel

  - Non linear lifting line model
  - Adjoint solver coded to provide sensitivities of post-processing functions with respect to design variables

B{DLLM} is using 1 programming language: B{Python}
  - B{Python} is used for the entire code an requires the numpy and scipy packages (U{http://www.python.org})
  - MDOTools library is required to be able to run the code
"""
__version__='0.1.0'