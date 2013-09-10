import sys
import numpy
from DLLM.polarManager.analyticAirfoil import AnalyticAirfoil
from DLLM.polarManager.airfoilPolar import AirfoilPolar
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver

class DLLMEvaluator:
    ERROR_MSG='ERROR in DLLMEvaluator.'
    
    #-- Constructor
    def __init__(self, wing_geom , OC):
        self.__wing_geom = wing_geom    # Wing geometry
        self.__OC        = OC           # Operating conditions
        
        self.__ref_airfoil     = None   # reference airfoil, only used if the same airfoil is put for all sections
        self.__airfoils        = None   # Airfoil list for each section
        self.__linked_airfoils = None   # Airfoil scaled to the the planform
        self.__DLLM            = None   # Differentiated lifting line model
        
    #-- Accessors
    def get_wing_geom(self):
        return self.__wing_geom
    
    def get_DLLM(self):
        return self.__DLLM
    
    #-- Setters
    def set_OC(self, OC):
        self.__OC = OC
        
    def set_ref_aifoil(self, ref_airfoil):
        self.__ref_airfoil = ref_airfoil
        
    def set_airfoils(self, airfoils):
        ERROR_MSG=self.ERROR_MSG+'set_airfoils: '
        if len(airfoils)!=self.__wing_geom.get_n_elements():
            print ERROR_MSG+'number of airfoils is different than the number of sections. attribute not set'
            self.__airfoils = None
        else:
            self.__airfoils = airfoils
            self.__link_airfoils_to_geom()
            
    #-- Public methods
    def setup(self):
        """
        Build lifting line model
        """
        print '\n--- DDLM Setup ---'
        sys.stdout.write('Checking inputs... ')
        self.__check_inputs()
        sys.stdout.write('Done.\n')
        sys.stdout.write('Build DLLM model... ')
        self.__build_DLLM_model()
        sys.stdout.write('Done.\n')
        
    def set_relax_factor(self, relax_factor):
        self.__DLLM.set_relax_factor(relax_factor)
        
    def set_stop_criteria(self, residual=None, n_it=None):
        self.__DLLM.set_stop_criteria(residual=residual, n_it=n_it)
        
    def run_direct(self):
        """
        Run the lifting line model
        """
        print '\n--- Run DLLM model - Direct mode ---'
        self.__DLLM.run_direct()
        
    def run_post(self, func_list=None):
        """
        Run the post processing for the lifting line model
        """
        print '\n--- Run DLLM model - Post mode ---'
        self.__DLLM.run_post(func_list=func_list)
        
    def run_adjoint(self):
        """
        Run the adjoint methods for the lifting line model
        """
        print '\n--- Run DLLM model - Adjoint mode ---'
        self.__DLLM.run_adjoint()
        
    def build_linear_airfoil(self, AoA0=0., Cd0=0., Cm0=0., relative_thickness=0., Sref=1., Lref=1., set_as_ref=True):
        degToRad = numpy.pi/180.
        airfoil  = AnalyticAirfoil(AoA0=degToRad*AoA0, Cd0=Cd0, Cm0=Cm0, relative_thickness= relative_thickness, Sref=Sref, Lref=Lref)
        if set_as_ref:
            self.set_ref_aifoil(airfoil)
        return airfoil
    
    def build_polar_airoil(self, database, Sref=1., Lref=1., interpolator='2DSpline', set_as_ref=True):
        # Why relative thickness usage ? The extraction from a polar should give us more freedom.
        airfoil = AirfoilPolar(database,relative_thickness=0.15, interpolator=interpolator, Sref=Sref, Lref=Lref)
        if set_as_ref:
            self.set_ref_aifoil(airfoil)
        return airfoil
    
    def build_airfoils_from_ref(self):
        ERROR_MSG=self.ERROR_MSG+'build_airfoils_from_ref: '
        if self.__ref_airfoil is None:
            print ERROR_MSG+'cannot build airfoils if reference airfoil is not defined.'
        else:
            airfoils=[]
            for n in xrange(self.__wing_geom.get_n_elements()):
                airfoils.append(self.__ref_airfoil.get_scaled_copy())
            
            self.set_airfoils(airfoils)

    #-- Private methods
    def __check_inputs(self):
        ERROR_MSG='\n'+self.ERROR_MSG+'__check_inputs: '
        checked=True
        if self.__airfoils is None:
            checked=False
            print ERROR_MSG+'airfoils attribute undefined. please set airfoils attribute.'
        elif len(self.__airfoils)!=self.__wing_geom.get_n_elements():
            checked=False
            print ERROR_MSG+'number of airfoils must equal to the number of geometrical sections.'
            
        if self.__linked_airfoils is None:
            checked=False
            print ERROR_MSG+'linked_airfoils attribute undefined.'
        
        if not checked:
            sys.exit(1)
            
    def __link_airfoils_to_geom(self):
        linked_airfoils=[]
        n=self.__wing_geom.get_n_elements()
        thicks=self.__wing_geom.get_relative_thickness()
        chords=self.__wing_geom.get_chords()
        for i in range(n):
            Lloc=chords[i]
            Sloc=self.__compute_local_Sref(Lloc)
            linked_airfoils.append(self.__airfoils[i].get_scaled_copy(thicks[i], Sloc, Lloc))
        self.__linked_airfoils = linked_airfoils
        
    def __compute_local_Sref(self,Lloc):
        wingspan = self.__wing_geom.get_wingspan()
        n = self.__wing_geom.get_n_elements()
        return Lloc*wingspan/float(n)
    
    def __build_DLLM_model(self):
        self.__DLLM=DLLMSolver(self.__wing_geom,self.__linked_airfoils, self.__OC)
    
        