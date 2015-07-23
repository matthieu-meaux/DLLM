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
# @author : Damien Guenot
import unittest
from DLLMWrappers.OpenMDAOWrapper import AerodynamicOptimization


class TestDLLM_MDAO_opt(unittest.TestCase):

    def test_optim(self):
        N = 10
        Target_Lift = 606570.049598
        opt_problem = AerodynamicOptimization(N=N, Target_Lift=Target_Lift)
        Initial_drag = 15805.247796847883
        Initial_LoD = Target_Lift / Initial_drag
        print "Running optimization"
        opt_problem.run()
        self.assertAlmostEqual(Target_Lift,
                               opt_problem.llw_comp.Lift,
                               places=3)
        self.assertGreater(Initial_drag, opt_problem.llw_comp.Drag)
        self.assertGreater(opt_problem.llw_comp.LoD, Initial_LoD)
        return

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLM_MDAO_opt)
    unittest.TextTestRunner(verbosity=2).run(suite)
