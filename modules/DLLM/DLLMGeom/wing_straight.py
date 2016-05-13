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

# - Local imports -
from DLLM.DLLMGeom.DLLM_param import DLLM_param
import numpy as np

class Wing_Straight(DLLM_param):
    def __init__(self, tag, n_sect=20, grad_active=True):
        """
        Constructor: set main attributes
        """
        DLLM_param.__init__(self, tag, n_sect=n_sect, grad_active=grad_active)
        
        self.set_perc_chord(0.25)
        
    # -- discretization methods
    def build_discretization(self):
        deg_to_rad = np.pi / 180.
        N          = self.get_n_sect()
        ndv        = self.get_ndv()
        r_list_eta = self.get_r_list_eta()
        
        #-- recover some data
        span_pt          = self.BC_manager.get_pt('span')
        span             = span_pt.get_value()
        span_grad        = span_pt.get_gradient() 
        
        sweep_pt         = self.BC_manager.get_pt('sweep')
        sweep            = sweep_pt.get_value()*deg_to_rad
        sweep_grad       = sweep_pt.get_gradient()*deg_to_rad
        
        root_chord_pt    = self.BC_manager.get_pt('root_chord')
        root_chord       = root_chord_pt.get_value()
        root_chord_grad  = root_chord_pt.get_gradient()

        root_height_pt   = self.BC_manager.get_pt('root_height')
        root_height      = root_height_pt.get_value()
        root_height_grad = root_height_pt.get_gradient()
        
        tip_height_pt    = self.BC_manager.get_pt('tip_height')
        tip_height       = tip_height_pt.get_value()
        tip_height_grad  = tip_height_pt.get_gradient()
        
        #-- Build and set twist arrays
        twist = np.zeros(N)
        twist_grad = np.zeros((N, ndv))
        for i in xrange(N / 2):
            pt1 = self.BC_manager.get_pt('twist'+str(N / 2 - 1 - i))
            twist[i] = pt1.get_value()*deg_to_rad
            twist_grad[i,:] = pt1.get_gradient()*deg_to_rad
            pt2 = self.BC_manager.get_pt('twist'+str(i))
            twist[N  / 2 + i] = pt2.get_value()*deg_to_rad
            twist_grad[N  / 2 + i,:] = pt2.get_gradient()*deg_to_rad 
        self.set_twist(twist)
        self.set_twist_grad(twist_grad)
        
        #-- build and set sweep arrays
        sweep_eta = np.zeros(N+1)
        sweep_grad_eta = np.zeros((N+1,ndv))
        sweep_eta[:]      = sweep
        sweep_grad_eta[:] = sweep_grad
        self.set_sweep_eta(sweep_eta)
        self.set_sweep_grad_eta(sweep_grad_eta)
        
        #-- build and set chords arrays
        chords_eta         = np.zeros((N + 1))
        chords_grad_eta    = np.zeros((N+ 1, ndv))
        chords_eta[:]      = root_chord
        chords_grad_eta[:] = root_chord_grad
        self.set_chords_eta(chords_eta)
        self.set_chords_grad_eta(chords_grad_eta)

        #********* TBC TBC TBC : DEPENDS ON AIRFOILS TYPE UNKNOWN HOW TO DEAL WITH THAT*************
        #-- build and set rel_thicks arrays (TBC: Depends on Airfoils)
        rel_thicks_eta = np.zeros(N+1)
        rel_thicks_grad_eta = np.zeros((N+1, ndv))
        for i, r in enumerate(r_list_eta):
                r = abs(2. * r)
                heights_eta = (tip_height - root_height)*r + root_height
                heights_grad_eta= (tip_height_grad - root_height_grad)*r +root_height_grad
                rel_thicks_eta[i] = heights_eta / chords_eta[i]
                rel_thicks_grad_eta[i] = (heights_grad_eta * chords_eta[i] - heights_eta * chords_grad_eta[i]) / (chords_eta[i])**2      
        self.set_rel_thicks_eta(rel_thicks_eta)
        self.set_rel_thicks_grad_eta(rel_thicks_grad_eta)   

        #-- build and set eta arrays
        eta      = np.zeros((3, N + 1))
        eta_grad = np.zeros((3, N + 1, ndv))
        for i, r in enumerate(r_list_eta):
            abs_r = abs(r)
            eta[0, i] = abs_r * span * np.sin(sweep) + 0.25 * root_chord
            eta_grad[0, i, :] = abs_r * span_grad[:] * np.sin(sweep) + abs_r * span * np.cos(sweep) * sweep_grad[:] + 0.25 * root_chord_grad[:]

            eta[1, i] = r * span
            eta_grad[1, i, :] = r * span_grad[:]
        
        self.set_eta(eta)
        self.set_eta_grad(eta_grad)
