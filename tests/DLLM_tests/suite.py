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

import unittest

from test_wing_param import TestWingParam
from test_DLLM_simple_base import TestDLLMSimpleBase
from test_DLLM_simple_R import TestDLLMSimpleR
from test_DLLM_simple_F import TestDLLMSimpleF
from test_DLLM_simple_Loads import TestDLLMSimpleLoads
from test_DLLM_simple_TCl_TLift import TestDLLMSimpleTClTLift
# from test_DLLM_meta import TestDLLMMeta

try:
    from openturns import *
    run_meta=True
except:
    run_meta=False

run_meta = False
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestWingParam))
    suite.addTest(unittest.makeSuite(TestDLLMSimpleBase))
    suite.addTest(unittest.makeSuite(TestDLLMSimpleR))
    suite.addTest(unittest.makeSuite(TestDLLMSimpleF))
    suite.addTest(unittest.makeSuite(TestDLLMSimpleLoads))
    suite.addTest(unittest.makeSuite(TestDLLMSimpleTClTLift))
    if run_meta:
        suite.addTest(unittest.makeSuite(TestDLLMMeta))
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=10).run(suite())