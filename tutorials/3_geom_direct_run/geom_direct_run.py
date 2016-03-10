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
#  @author : Matthieu MEAUX
#

#-- Imports
import numpy as np
from MDOTools.OC.operating_condition import OperatingCondition
from DLLM.DLLMGeom.DLLM_Geom import DLLM_Geom
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver

#-- Build Operating conditions
OC=OperatingCondition('cond1',atmospheric_model='ISA')
OC.set_Mach(0.8)
OC.set_AoA(3.5)
OC.set_altitude(10000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

#-- Build Geom class
span  = 40.
sweep = 30.

sweep_rad = sweep*np.pi/180.

wing_geom = DLLM_Geom('tuto3_geom',n_sect=20, grad_active=False)
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
sweep_eta[:]      = sweep_rad
for i, r in enumerate(r_list):
    abs_r = abs(r)
    eta[0, i] = abs_r * span * np.sin(sweep_rad) + 0.25 * chords_eta[i]
    eta[1, i] = r * span

wing_geom.set_twist(twist)
wing_geom.set_sweep_eta(sweep_eta)
wing_geom.set_chords_eta(chords_eta)
wing_geom.set_rel_thicks_eta(rel_thicks_eta)
wing_geom.set_eta(eta)
wing_geom.build_linear_airfoil(OC, AoA0=0.0, set_as_ref=True)
wing_geom.build_airfoils_from_ref()
wing_geom.update()  
wing_geom.plot()
print wing_geom

DLLM = DLLMSolver('tuto3', wing_geom, OC, verbose=1, grad_active=False)
DLLM.run_direct()
DLLM.run_post()
DLLM.plot()
DLLM.export_F_list()

