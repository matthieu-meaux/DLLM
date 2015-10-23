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
# @author : Andre Mendes
# @author : Matthieu MEAUX
#
from DLLM.DLLMGeom.wing_param import Wing_param
from MDOTools.OC.operating_condition import OperatingCondition
import numpy
import string

OC=OperatingCondition('cond1',atmospheric_model='simple')

OC.set_altitude(10000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)

wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
wing_param.build_wing()
wing_param.set_value('span',2.392)
#wing_param.set_value('sweep',34.)
wing_param.set_value('sweep',26.7)
wing_param.set_value('break_percent',33.)
wing_param.set_value('root_chord',0.806)
wing_param.set_value('break_chord',0.689)
wing_param.set_value('tip_chord',0.451)
#wing_param.set_value('root_height',1.28)
wing_param.set_value('root_height',0.0782)
#wing_param.set_value('break_height',0.97)
wing_param.set_value('break_height',0.0668)
#wing_param.set_value('tip_height',0.33)
wing_param.set_value('tip_height',0.0438)
#wing_param.convert_to_design_variable('span',10.,50.)
#wing_param.convert_to_design_variable('sweep',0.,40.)
#wing_param.convert_to_design_variable('break_percent',20.,40.)
#wing_param.convert_to_design_variable('root_chord',5.,7.)
#wing_param.convert_to_design_variable('break_chord',3.,5.)
#wing_param.convert_to_design_variable('tip_chord',1.,2.)
#wing_param.convert_to_design_variable('root_height',1.,1.5)
#wing_param.convert_to_design_variable('break_height',0.8,1.2)
#wing_param.convert_to_design_variable('tip_height',0.2,0.5)
#wing_param.convert_to_design_variable('root_height',0.7,1.)
#wing_param.convert_to_design_variable('break_height',0.45,0.8)
#wing_param.convert_to_design_variable('tip_height',0.10,0.26)
wing_param.build_linear_airfoil(OC, AoA0=0., Cm0=-0.1, set_as_ref=True)
#wing_param.build_meta_airfoil(OC, '../ONERA_D.xml', relative_thickness=.09, camber=0, Sref=1., Lref=1., sweep=.0, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()


list_file = ['MetaOneraM6_3.00E-01.dat','MetaOneraM6_6.00E-01.dat','MetaOneraM6_8.00E-01.dat']
e=1.
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
            new_words = ['#','AoA','Cl','Cdp','Cdp_theory']
        if len(words)> 0:
            if words[0] != '#':
                AoA=eval(words[0])
                Cl=eval(words[7])
                Cdp=eval(words[9])
                Cdp_theory=Cl**2/(numpy.pi*e*AR)
                new_words=[str(AoA),str(Cl),str(Cdp),str(Cdp_theory)]
        
        line = string.join(new_words,' ')
        fid.write(line+'\n')
    fid.close()