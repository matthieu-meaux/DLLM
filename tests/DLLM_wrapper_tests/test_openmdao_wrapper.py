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
# @author : Francois Gallard

from DLLMWrappers.OpenMDAOWrapper import DLLMOpenMDAOComponent
from MDOTools.ValidGrad.FDValidGrad import FDValidGrad
import unittest
import numpy as np
from DLLM.DLLMKernel.DLLMPost import DLLMPost
from MDOTools.OC.operating_condition import OperatingCondition


class TestDLLM_wrapper(unittest.TestCase):

    def __init_DLLM(self):
        N = 5
        Target_Lift = 606570.049598
        
        DLLMOpenMDAO = DLLMOpenMDAOComponent(N=5,
            Target_Lift=Target_Lift, verbose=0)
        return DLLMOpenMDAO, N

    def test_init_component(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        default_lift = 0.00
        self.assertEqual(DLLMOpenMDAO.Lift, default_lift)

    def test_change_OC(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        self.assertEqual(10000., DLLMOpenMDAO.OC.get_altitude())
        DLLMOpenMDAO.Mach = 0.1
        DLLMOpenMDAO.altitude = 0.
        DLLMOpenMDAO.T0 = OperatingCondition.T0 + 10
        DLLMOpenMDAO.P0 = OperatingCondition.P0 + 10
        DLLMOpenMDAO.execute()
        self.assertEqual(0., DLLMOpenMDAO.OC.get_altitude())
        self.assertEqual(0.1, DLLMOpenMDAO.OC.get_Mach())
        self.assertEqual(OperatingCondition.T0 + 10, DLLMOpenMDAO.OC.get_T0())
        self.assertEqual(OperatingCondition.P0 + 10, DLLMOpenMDAO.OC.get_P0())
        self.assertAlmostEqual(DLLMOpenMDAO.Target_Lift, DLLMOpenMDAO.Lift, places=2)
        self.assertAlmostEqual(335914.921978, DLLMOpenMDAO.Drag, places=2)


    def test_change_param(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        DLLMOpenMDAO.execute()
        span = 48.0
        twist = -1.
        sweep = 31.
        break_percent = 30.
        root_chord = 7.
        break_chord = 5.0
        tip_chord = 2.0
        root_height = 1.30
        break_height = 0.99
        tip_height = 0.38
        DLLMOpenMDAO.rtwist = twist * np.ones((5))
        DLLMOpenMDAO.span = span
        DLLMOpenMDAO.sweep = sweep
        DLLMOpenMDAO.break_percent = break_percent
        DLLMOpenMDAO.root_chord = root_chord
        DLLMOpenMDAO.break_chord = break_chord
        DLLMOpenMDAO.tip_chord = tip_chord
        DLLMOpenMDAO.root_height = root_height
        DLLMOpenMDAO.break_height = break_height
        DLLMOpenMDAO.tip_height = tip_height
        DLLMOpenMDAO.execute()
        self.assertEqual(span, DLLMOpenMDAO.get_dv_value('span'))
        self.assertEqual(sweep, DLLMOpenMDAO.get_dv_value('sweep'))
        self.assertEqual(break_percent, DLLMOpenMDAO.get_dv_value('break_percent'))
        self.assertEqual(root_chord, DLLMOpenMDAO.get_dv_value('root_chord'))
        self.assertEqual(break_chord, DLLMOpenMDAO.get_dv_value('break_chord'))
        self.assertEqual(tip_chord, DLLMOpenMDAO.get_dv_value('tip_chord'))
        self.assertEqual(root_height, DLLMOpenMDAO.get_dv_value('root_height'))
        self.assertEqual(break_height, DLLMOpenMDAO.get_dv_value('break_height'))
        self.assertEqual(tip_height, DLLMOpenMDAO.get_dv_value('tip_height'))
        for i in range(N) :
            self.assertEqual(twist, DLLMOpenMDAO.get_dv_value('rtwist%s' % i))
        self.assertAlmostEqual(DLLMOpenMDAO.Target_Lift, DLLMOpenMDAO.Lift, places=2)
        self.assertAlmostEqual(14013.884850035265, DLLMOpenMDAO.Drag, places=2)
        return

    def test_get_dv_id_list(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        for dv_id in DLLMOpenMDAO.get_dv_id_list() :
            if dv_id.startswith('rtwist') :
                dv_id = 'rtwist'
            if dv_id not in dir(DLLMOpenMDAO) :
                print "Can not find %s in %s" % (dv_id, dir(DLLMOpenMDAO))
                assert(False)

    def test_get_F_list_names(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        for f in DLLMOpenMDAO.get_F_list_names() :
            if f not in dir(DLLMOpenMDAO) :
                print "Can not find %s in %s" % (f, dir(DLLMOpenMDAO))
                assert(False)

    def test_get_dv_value(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        for dv_id in DLLMOpenMDAO.get_dv_id_list() :
            value = DLLMOpenMDAO.get_dv_value(dv_id)
            if dv_id.startswith('rtwist') :
                index_twist = int(dv_id.replace('rtwist', ''))
                rtwist = getattr(DLLMOpenMDAO, 'rtwist')
                self.assertEqual(value, rtwist[index_twist])
            else :
                self.assertEqual(value, getattr(DLLMOpenMDAO, dv_id))
#         DLLMOpenMDAO.execute()

    def test_get_F_value(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        DLLMOpenMDAO.execute()
        f_value = DLLMOpenMDAO.get_F_list()
        for i, f_name in enumerate(DLLMOpenMDAO.get_F_list_names()) :
            self.assertEqual(f_value[i], getattr(DLLMOpenMDAO, f_name))


    def update_dllm_from_x(self, DLLMOpenMDAO, x):
        DLLMOpenMDAO.rtwist = x[:5]
        DLLMOpenMDAO.span = x[5]
        DLLMOpenMDAO.sweep = x[6]
        DLLMOpenMDAO.break_percent = x[7]
        DLLMOpenMDAO.root_chord = x[8]
        DLLMOpenMDAO.break_chord = x[9]
        DLLMOpenMDAO.tip_chord = x[10]
        DLLMOpenMDAO.root_height = x[11]
        DLLMOpenMDAO.break_height = x[12]
        DLLMOpenMDAO.tip_height = x[13]
        return

    def test_exec(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        x0 = DLLMOpenMDAO.get_dv_array()
        self.update_dllm_from_x(DLLMOpenMDAO, x0)
        DLLMOpenMDAO.Mach = 0.4
        DLLMOpenMDAO.execute()
#         OC = DLLMOpenMDAO.OC
#         print "Final :"
#         print "    OC :"
#         print "        AoA :",OC.get_AoA()
#         print "        P0:",OC.get_P0()
#         print "        T0:",OC.get_T0()
#         print "        Alt.:",OC.get_altitude()
#         print "        Mach:",OC.get_Mach()
#         print "    Outputs :"
#         for i, varname in enumerate(DLLMPost.DEF_F_LIST_NAMES) :
#             print "        Varname : %s = %s" %(varname,DLLMOpenMDAO.get_F_list()[i])
#         print "    Inputs"
#         for i,varname in enumerate(DLLMOpenMDAO.wing_param.get_dv_id_list()):
# print "        Varname : %s=%s"
# %(varname,DLLMOpenMDAO.wing_param.get_dv_array()[i])
        self.assertAlmostEqual(
            DLLMOpenMDAO.Target_Lift,
            DLLMOpenMDAO.Lift,
            places=3)
        self.assertAlmostEqual(32951.826392216724, DLLMOpenMDAO.Drag, places=2)

    def test_list_deriv_vars(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        DLLMOpenMDAO.execute()
        self.assertAlmostEqual(
            DLLMOpenMDAO.Target_Lift,
            DLLMOpenMDAO.Lift,
            places=3)
        out_dvid_tuple, F_list_names_tuple = DLLMOpenMDAO.list_deriv_vars()
        self.assertTrue(isinstance(out_dvid_tuple, tuple))
        self.assertTrue(isinstance(F_list_names_tuple, tuple))
        F_list_names_ref = tuple(DLLMPost.DEF_F_LIST_NAMES)
        self.assertTupleEqual(F_list_names_ref, F_list_names_tuple)

    def test_jacobian(self):
        DLLMOpenMDAO, N = self.__init_DLLM()
        DLLMOpenMDAO.execute()
        out_dvid_tuple, F_list_names_tuple = DLLMOpenMDAO.list_deriv_vars()
        n_F = len(F_list_names_tuple)
        Jacobian = DLLMOpenMDAO.provideJ()
        n_dv = len(DLLMOpenMDAO.get_dv_id_list())
        self.assertEqual(n_dv, Jacobian.shape[0])
        self.assertEqual(n_F, Jacobian.shape[1])
        DLLMOpenMDAO.Mach = 0.666

        def d_obj(x):
            self.update_dllm_from_x(DLLMOpenMDAO, x)
            DLLMOpenMDAO.execute()
            return DLLMOpenMDAO.provideJ()

        def obj(x):
            self.update_dllm_from_x(DLLMOpenMDAO, x)
            DLLMOpenMDAO.execute()
            return DLLMOpenMDAO.get_F_list()

        treshold = 1e-2
        val_grad = FDValidGrad(2, obj, d_obj, fd_step=1.e-6)
        x0 = DLLMOpenMDAO.get_dv_array()
        ok, df_fd, df = val_grad.compare(
            x0, treshold=treshold, split_out=True, return_all=True, iprint=False)
        varnames, funcnames = DLLMOpenMDAO.list_deriv_vars()
        varnames = np.concatenate(
            (['tw0', 'tw1', 'tw2', 'tw3', 'tw4'], np.array(varnames[1:])))
        funcnames = np.array(funcnames)
        k = 0
        ok = True
        for fd, gr in zip(df_fd, df):
            err = np.where(((fd - gr) / fd > treshold) & (fd > 2e-5))[0]
            if err.size > 0:
                # Compute absolute error when gradient or fd is lower than
                # threshold
                inf_tol = err[np.where((np.abs(gr[err]) < treshold) | 
                                       (np.abs(fd[err]) < treshold))[0]]
                # Scale gradient by function value for the lift issue in
                # particular
                f_val = getattr(DLLMOpenMDAO, funcnames[k])
                actually_ok = np.where(
                    ((fd[inf_tol] - gr[inf_tol]) < treshold * f_val))[0]
                real_err = np.delete(err, actually_ok)
                if real_err.size > 0:
                    print "Errror in d", funcnames[k], " /d ", varnames[real_err]
                    print "df analytic = ", gr[real_err]
                    print "df finite d = ", fd[real_err]
                    print "Relative =", (fd[real_err] - gr[real_err]) / fd[real_err]
                    ok = False
            k = k + 1
        assert(ok)
        # OpenMDAO check grad method is not that reliable...
        # threshold is hard coded at 1e-5, only relative error is computed...
        # out = DLLMOpenMDAO.check_gradient(fd_form='central', fd_step=1e-6, fd_step_type='relative')
        # assert(len(out[-1])==0)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDLLM_wrapper)
    unittest.TextTestRunner(verbosity=2).run(suite)
