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
import string

list_file=['mesh_ONERAM6_inv_Mach_03.dat','mesh_ONERAM6_inv_Mach_06.dat','mesh_ONERAM6_inv_Mach_08.dat']
list_corr=[0.0010251106,0.0008860518,0.0012233352]
for i,file_name in enumerate(list_file):
    fid=open(file_name,'r')
    lines=fid.readlines()
    fid.close()
    
    fid=open('corr_'+file_name,'w')
    for line in lines:
        words=string.split(line)
        if len(words)> 0:
            if words[0] != '#':
                words[2]=str(eval(words[2])-list_corr[i])
        new_line = string.join(words,' ')
        fid.write(new_line+'\n')
    fid.close()
        