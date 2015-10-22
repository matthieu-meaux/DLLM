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

from DLLMWrappers.OpenMDAOWrapper import AerodynamicOptimization
        

import time
tt = time.time()
N = 10
opt_problem = AerodynamicOptimization(N)
opt_problem.configure()
opt_problem.DLLMOpenMDAO.execute()
print "Initial lift =", opt_problem.DLLMOpenMDAO.Lift
print "Initial drag =", opt_problem.DLLMOpenMDAO.Drag
print "Initial lift over drag =", opt_problem.DLLMOpenMDAO.LoD
print "Running optimization"
opt_problem.run()
print "*** Elapsed time: ", time.time() - tt, "seconds ***"

print "Number of function evaluations =", opt_problem.DLLMOpenMDAO.exec_count
print "Number of gradient evaluations =", opt_problem.DLLMOpenMDAO.derivative_exec_count
print "Final lift =", opt_problem.DLLMOpenMDAO.Lift
print "Final drag =", opt_problem.DLLMOpenMDAO.Drag
print "Final lift over drag =", opt_problem.DLLMOpenMDAO.LoD
for i_dv_info in opt_problem.DLLMOpenMDAO.get_dv_info_list():
    print 5 * " ", i_dv_info
