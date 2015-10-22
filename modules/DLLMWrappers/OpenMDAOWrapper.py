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
# @author : Francois Gallard
# @author : Damien Guenot

# pylint: disable-msg=E0611,F0401
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMTargetLift import DLLMTargetLift
from MDOTools.OC.operating_condition import OperatingCondition

from openmdao.main.api import Component
from openmdao.main.datatypes.api import Float, Array
from openmdao.lib.drivers.api import SLSQPdriver
from openmdao.main.api import Assembly

import numpy as np


class DLLMOpenMDAOComponent(Component):
    # set up interface to the framework
    # pylint: disable-msg=E1101
    # Outputs of lifting line problem
    """OpenMDAO component for DLLM implementation
    """
    
    Lift = Float(iotype='out', desc='Lift')
    Drag = Float(iotype='out', desc='Drag')
    Drag_Pressure = Float(iotype='out', desc='Drag_Pressure')
    Drag_Induced = Float(iotype='out', desc='Drag_Induced')
    Drag_Wave = Float(iotype='out', desc='Drag_Wave')
    Drag_Friction = Float(iotype='out', desc='Drag_Friction')
    Cd = Float(iotype='out', desc='Cd')
    Cdp = Float(iotype='out', desc='Cdp')
    Cdi = Float(iotype='out', desc='Cdi')
    Cdw = Float(iotype='out', desc='Cdw')
    Cdf = Float(iotype='out', desc='Cdf')
    Cl = Float(iotype='out', desc='Cl')
    LoD = Float(iotype='out', desc='LoD')
    Sref = Float(iotype='out', desc='Sref')

    # Design variables of lifting line problem
    rtwist = Array([], desc='rtwist', iotype="in")
    span = Float(desc='span', default_value=34., iotype="in")
    sweep = Float(desc='sweep', default_value=34., iotype="in")
    break_percent = Float(desc='break_percent', default_value=33., iotype="in")
    root_chord = Float(desc='root_chord', default_value=6.1, iotype="in")
    break_chord = Float(desc='break_chord', default_value=4.6, iotype="in")
    tip_chord = Float(desc='tip_chord', default_value=1.5, iotype="in")
    root_height = Float(desc='root_height', default_value=1.28, iotype="in")
    break_height = Float(desc='break_height', default_value=0.97, iotype="in")
    tip_height = Float(desc='tip_height', default_value=0.33, iotype="in")
    # Operating conditions variables
    Mach = Float(iotype='in', default_value=0.7, desc='Mach')
    altitude = Float(iotype='in', default_value=10000., desc='Altitude')
    T0 = Float(iotype='in',default_value=OperatingCondition.T0,
        desc='Ground ISA ref Temperature')
    P0 = Float(iotype='in',default_value=OperatingCondition.P0,
        desc='Ground ISA ref Pressure')

    def __init__(self,Target_Lift, N = 10, verbose=0):
        """Initialization of DLLM component.
        DLLM component use target lift capability of DLLM kernel
            @param Target_Lift : the targeted lift value (float)
            @param N : integer. Number of discrete section on 1/2 wing
            @param verbose : integer : verbosity level
        """
        try :
            float(Target_Lift)
        except:
            raise ValueError('You MUST define a float target lift value, get '+str(Target_Lift)+' instead.')

        self.Target_Lift = Target_Lift
        self.N = N
        self.OC = None
        self.rtwist = np.zeros(N)
        self.__display_wing_param = True
        self.__verbose = verbose
        self.__wing_param = None
        super(DLLMOpenMDAOComponent, self).__init__()        
        self.OC = OperatingCondition(tag='DLLMOC', atmospheric_model='ISA')
        self.__set_OC_values()
        self.__set_wing_param()
        self.__DLLM = DLLMTargetLift('test', self.__wing_param,
                                   self.OC, verbose=self.__verbose)
        self.__DLLM.set_target_Lift(Target_Lift)
        self.__DLLM.run_direct()
        self.__DLLM.run_post()

    def __set_wing_param(self, wing_param_name='test_param'):
        """Method for wing parameters setting : design variables initial values and bounds
        @param wing_param_name : wing parametrization names
        """
        self.__set_wing_param_values(wing_param_name=wing_param_name)

        self.__set_wing_param_bounds()
        self.__wing_param.build_linear_airfoil(self.OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        self.__wing_param.build_airfoils_from_ref()
        self.__wing_param.update()

        if self.__display_wing_param:
            self.__display_wing_param = False
            print self.__wing_param

    def __set_wing_param_values(self, wing_param_name='test_param'):
        """Method for wing parameters variables setting
        @param wing_param_name : wing parametrization names
        """
        self.__wing_param = Wing_param(wing_param_name,
                                     geom_type='Broken', n_sect=self.N * 2)
        self.__wing_param.build_wing()
        for i in xrange(self.N):
            self.__wing_param.set_value('rtwist%s' % i, 0.)
        for param in Wing_param.DISCRETE_ATTRIBUTES_LIST:
            self.__wing_param.set_value(param, getattr(self, param))

    def __set_wing_param_bounds(self):
        """Method for desing variables bounds settings
        Values are set to inf/-inf in DLLM component and their 'real' bounds
        are defined when optimization problem is set
        """
        for i in xrange(self.N):
            self.__wing_param.convert_to_design_variable(
                'rtwist%s' % i, (-float('inf'), float('inf')))
        for param in Wing_param.DISCRETE_ATTRIBUTES_LIST:
            self.__wing_param.convert_to_design_variable(param, (-float('inf'), float('inf')))
        return

    def __set_OC_values(self):
        """ Set default operating conditions"""
        self.OC.set_Mach(self.Mach)
        self.OC.set_altitude(self.altitude)
        self.OC.set_T0(self.T0)
        self.OC.set_P0(self.P0)
        self.OC.set_humidity(0.)
        self.OC.compute_atmosphere()
    
    def __update_wing_parame_values(self):
        for dv_id in self.__wing_param.get_dv_id_list():
            if dv_id.startswith('rtwist'):
                i_twist = int(dv_id.replace('rtwist', ''))
                self.__wing_param.set_value(dv_id, self.rtwist[i_twist])
            else:
                self.__wing_param.set_value(dv_id,getattr(self,dv_id))
        self.__wing_param.build_linear_airfoil(self.OC, AoA0=-2., Cm0=-0.1, set_as_ref=True)
        self.__set_wing_param_bounds()
        self.__wing_param.build_airfoils_from_ref()
        self.__wing_param.update()
        return

    def execute(self):
        """ Perform a DLLM computation with the """
        self.__update_wing_parame_values()        
        self.__set_OC_values()
        self.__DLLM.set_target_Lift(self.Target_Lift)
        self.__DLLM.run_direct()
        self.__DLLM.run_post()
        output = self.__DLLM.get_F_list()
        for f, f_name in zip(output, self.__DLLM.get_F_list_names()):
            setattr(self, f_name, f)

    def list_deriv_vars(self):
        """specify the inputs and outputs where derivatives are defined
        Specific treatment for twist : defined as rtwist0, rtwist1,... in DLLM
        but as an array rtwist in openmdao component"""

        out_dvid = []
        for dv_id in self.__wing_param.get_dv_id_list():
            if dv_id.startswith('rtwist0'):
                out_dvid.append('rtwist')
            elif not dv_id.startswith('rtwist'):
                out_dvid.append(dv_id)
        return tuple(out_dvid), tuple(self.__DLLM.get_F_list_names())

    def provideJ(self):
        """Calculate the Jacobian according inputs and outputs"""
        self.__DLLM.run_adjoint()
        return np.array(self.__DLLM.get_dF_list_dchi())

    def get_F_list(self):
        return self.__DLLM.get_F_list()

    def get_F_list_names(self):
        return self.__DLLM.get_F_list_names()
    
    def get_dv_array(self):
        return self.__wing_param.get_dv_array()
    
    def get_dv_id_list(self):
        return self.__wing_param.get_dv_id_list()
    
    def get_dv_info_list(self):
        return self.__wing_param.get_dv_info_list()
    
    def get_dv_value(self,dv_id):
        return self.get_dv_array()[self.get_dv_id_list().index(dv_id)]
    
    
class AerodynamicOptimization(Assembly):

    """Constrained optimization of the LLW Component."""

    def __init__(self, Target_Lift, N = 10, verbose = 0):
        """Initialization of OpenMDAO component.
        DLLM component use target lift capability of DLLM kernel
            @param Target_Lift : the targeted lift value (float)
            @param N : integer. Number of discrete section on 1/2 wing
            @param verbose : integer : verbosity level
        """
        try :
            float(Target_Lift)
        except:
            raise ValueError('You MUST define a float target lift value, get '+str(Target_Lift)+' instead.')

        self.Target_Lift = Target_Lift
        self.N = N
        self.verbose = verbose
        super(AerodynamicOptimization, self).__init__()

    def configure(self):
        print "Configuration openmdao optimizer"
        # Create Optimizer instance
        self.add('driver', SLSQPdriver())
        # Create Paraboloid component instances
        workflow_name = 'llw_comp'
        self.DLLMOpenMDAO = DLLMOpenMDAOComponent(N=self.N, Target_Lift = self.Target_Lift, verbose=self.verbose)
        self.DLLMOpenMDAO.Mach = 0.6
        #self.DLLMOpenMDAO.DLLM.set_target_Lift(self.Target_Lift)
        # DLLMOpenMDAO.check_gradient()
        self.add(workflow_name, self.DLLMOpenMDAO)

        # Iteration Hierarchy
        self.driver.workflow.add(workflow_name)

        # SLSQP Flags
        self.driver.iprint = 0
        self.driver.accuracy = 1e-6
        self.driver.maxiter = 20

        # Objective

        self.driver.add_objective('%s.Drag' % workflow_name)

        # Design Variables
        print "Adding desing variables to openMDAO problem "
        # print 5*" "+"Variable name"+20*' '+'Lower bound'+20*" "+'Upper bound'
        rtwist_lb_array = -12.*np.ones(self.N)
        rtwist_ub_array = 12.*np.ones(self.N)
        self.driver.add_parameter('%s.rtwist' %workflow_name,
            low=rtwist_lb_array,high=rtwist_ub_array)
        self.driver.add_parameter('%s.span'%workflow_name,10.,50.,34.)
        self.driver.add_parameter('%s.sweep'%workflow_name,0.,40.,34.)
        self.driver.add_parameter('%s.break_percent'%workflow_name,20.,40.,33.)
        self.driver.add_parameter('%s.root_chord'%workflow_name,5.,7.,6.1)
        self.driver.add_parameter('%s.break_chord'%workflow_name,3.,5.,4.6)
        self.driver.add_parameter('%s.tip_chord'%workflow_name,1.,2.,1.5)
        self.driver.add_parameter('%s.root_height'%workflow_name,1.,1.5,1.28)
        self.driver.add_parameter('%s.break_height'%workflow_name,0.8,1.2,0.97)
        self.driver.add_parameter('%s.tip_height'%workflow_name,0.2,0.5,0.33)        
