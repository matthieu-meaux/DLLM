# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
#  MDOTools (Multi-disciplinary Optimization Tools, open source software)
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
#  @author : Matthieu Meaux

#-- Imports
import numpy as np
from MDOTools.OC.operating_condition import OperatingCondition
from DLLM.DLLMGeom.DLLM_Geom import DLLM_Geom

#-- Build Operating conditions
OC=OperatingCondition('cond_tuto',atmospheric_model='ISA')
OC.set_Mach(0.8)
OC.set_AoA(3.5)
OC.set_altitude(10000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

#-- Build Geom class
span  = 40.
sweep = 34.

wing_geom = DLLM_Geom('tuto2_geom',n_sect=20, grad_active=False)
#-- initialize arrays
wing_geom.build_r_lists()
r_list = wing_geom.get_r_list_eta()
N      = wing_geom.get_n_sect()
twist          = np.zeros(N)
sweep_eta      = np.zeros(N+1)
chords_eta     = np.zeros(N+1)
rel_thicks_eta = np.zeros(N+1)
eta            = np.zeros((3,N+1))
#-- set data
twist[:]          = 0. # 0. twist for all sections
chords_eta[:]     = 5. # 5m chord for all sections
rel_thicks_eta[:] = 0.15 # same rel_thick for all sections
for i, r in enumerate(r_list):
    abs_r = abs(r)
    eta[0, i] = abs_r * span * np.sin(sweep) + 0.25 * chords_eta[i]
    eta[1, i] = r * span


wing_geom.set_twist(twist)
wing_geom.set_sweep_eta(sweep_eta)
wing_geom.set_chords_eta(chords_eta)
wing_geom.set_rel_thicks_eta(rel_thicks_eta)
wing_geom.set_eta(eta)
wing_geom.build_linear_airfoil(OC, AoA0=0.0, Cm0=-0.1, set_as_ref=True)
wing_geom.build_airfoils_from_ref()
wing_geom.update() 

print wing_geom

wing_geom.plot() 