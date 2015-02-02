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
from DLLM.DLLMGeom.wing_param import Wing_param
from MDOTools.OC.operating_condition import OperatingCondition
import numpy
import string

OC=OperatingCondition('cond1')

OC.set_altitude(3000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)

wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
wing_param.build_wing()
wing_param.set_value('test_param.span',34.1)
wing_param.set_value('test_param.sweep',0.)
wing_param.set_value('test_param.break_percent',33.)
wing_param.set_value('test_param.root_chord',3.809036635014663)
wing_param.set_value('test_param.break_chord',3.809036635014663)
wing_param.set_value('test_param.tip_chord',3.809036635014663)
wing_param.set_value('test_param.root_height',0.8)
wing_param.set_value('test_param.break_height',0.8)
wing_param.set_value('test_param.tip_height',0.8)
wing_param.build_linear_airfoil(OC, AoA0=0., Cm0=-0.1, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

airfoils=wing_param.get_linked_airfoils()
Cla=airfoils[0].ClAlpha(0.0,0.3)

print "Cla=",Cla

#list_file = ['SimpleRect_3.00E-01.dat','SimpleRect_6.00E-01.dat','SimpleRect_8.00E-01.dat']
list_file = ['RectSweep0_3.00E-01.dat']
e=0.70
AR = wing_param.get_AR()
print 'AR=',AR

for i,file_name in enumerate(list_file):
    fid=open(file_name,'r')
    lines=fid.readlines()
    fid.close()
    
    fid=open('Cdp_'+file_name,'w')
    for j,line in enumerate(lines):
        words=string.split(line)
        new_words=words
        if j == 1:
            new_words = ['#','AoA','Cl','Cdp','Cl_theory','Cdp_theory','Cdp_elliptic']
        if len(words)> 0:
            if words[0] != '#':
                AoA=eval(words[0])
                Cl=eval(words[7])
                Cdp=eval(words[9])
                Cl_theory=AR*Cla/(AR+2.)*AoA*numpy.pi/180.
                Cdp_theory=Cl_theory**2/(numpy.pi*e*AR)
                Cdp_ell = Cl_theory**2/(numpy.pi*AR)
                new_words=[str(AoA),str(Cl),str(Cdp),str(Cl_theory),str(Cdp_theory),str(Cdp_ell)]
        
        line = string.join(new_words,' ')
        fid.write(line+'\n')
    fid.close()

fid=open('theory_iAoA.dat','w')
for i in xrange(20):
    fid.write(str(i)+'\t'+str(-Cl/(numpy.pi*AR))+'\t'+str(-Cl/(numpy.pi*e*AR))+'\n')
fid.close()