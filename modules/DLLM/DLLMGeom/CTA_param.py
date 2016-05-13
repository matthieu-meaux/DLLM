# - Local imports -
import string
import numpy as np
from scipy.interpolate import interp1d
from DLLM.DLLMGeom.DLLM_param import DLLM_param
from DLLM.polarManager.RefCTAAirfoil import RefCTAAirfoil

class CTA_param(DLLM_param):

    def __init__(self, tag, n_sect=40, grad_active=False):
        """
        Constructor: set main attributes
        """
        DLLM_param.__init__(self, tag, n_sect=n_sect, grad_active=grad_active)
        
        self.set_perc_chord(0.25)
        
        self.__n_af = 9
        
    def set_n_af(self, n_af):
        self.__n_af = n_af
        
    # -- discretization methods
    def build_discretization(self):
        deg_to_rad = np.pi / 180.
        N          = self.get_n_sect()
        ndv        = self.get_ndv()
        r_list_eta = self.get_r_list_eta()
        
        #-- recover some data
        half_span_pt           = self.BC_manager.get_pt('half_span')
        half_span              = half_span_pt.get_value()
        
        y_fact      = np.zeros(self.__n_af)
        sweep_list  = np.zeros(self.__n_af)
        twist_list  = np.zeros(self.__n_af)
        chord_list  = np.zeros(self.__n_af)
        pos_25_list = np.zeros(self.__n_af)
        
        for i in xrange(self.__n_af):
            y_fact_pt = self.BC_manager.get_pt('section'+str(i)+'.y_fact')
            y_fact[i] = y_fact_pt.get_value()
            
            sweep_pt  = self.BC_manager.get_pt('section'+str(i)+'.sweep')
            sweep_list[i] = sweep_pt.get_value()
            
            twist_pt  = self.BC_manager.get_pt('section'+str(i)+'.twist')
            twist_list[i] = twist_pt.get_value()
            
            chord_pt  = self.BC_manager.get_pt('section'+str(i)+'.chord')
            chord_list[i] = chord_pt.get_value()
        
        y_list = y_fact*half_span
        
        pos_25_list[0] = 0.25*chord_list[0]
        int_dy = 0.
        for i in xrange(1, self.__n_af):
            dy = y_list[i]-y_list[i-1]
            int_dy +=dy
            pos_25_list[i] = pos_25_list[i-1]+np.sin(sweep_list[i-1]*np.pi/180.)*dy
            
        sweep_interp  = interp1d(y_list, sweep_list, kind='cubic')
        twist_interp  = interp1d(y_list, twist_list, kind='cubic')
        chord_interp  = interp1d(y_list, chord_list, kind='cubic')
        pos_25_interp = interp1d(y_list, pos_25_list, kind='cubic')
        
        
        r_list_y = 0.5*(r_list_eta[1:]+r_list_eta[:-1])
        Y_eval_list = 2.*r_list_y*half_span
        #-- Build and set twist arrays
        twist = twist_interp(abs(Y_eval_list))*np.pi/180.
        self.set_twist(twist)
         
        ETA_eval_list = 2.*r_list_eta*half_span
        
        #-- build and set sweep arrays
        sweep_eta = sweep_interp(abs(ETA_eval_list))
        self.set_sweep_eta(sweep_eta)
        
        chords_eta = chord_interp(abs(ETA_eval_list))
        self.set_chords_eta(chords_eta)
        
        pos_25_eta = pos_25_interp(abs(ETA_eval_list))
        
        #-- build and set eta arrays
        eta      = np.zeros((3, N + 1))
        eta[1, :] = ETA_eval_list[:]
        eta[0, :] = pos_25_eta[:]        
        self.set_eta(eta)
        
    def build_ref_CTA_airfoil(self, OC, Sref=1., Lref=1., y_def_list=None, file_def_list=None, set_as_ref=True):
        self.set_airfoil_type('ref_CTA')
        airfoil = RefCTAAirfoil(OC, Sref=Sref, Lref=Lref)
        airfoil.set_y_def_list(y_def_list)
        airfoil.set_file_def_list(file_def_list)
        airfoil.init_interpolators()

        if set_as_ref:
            self.set_ref_aifoil(airfoil)
        return airfoil
    
    def test_airfoils(self):
        airfoils = self.get_linked_airfoils()
        print 'INFO        Lref    Sref    Cl    Cdw    Cdvp    Cdf    pcop'
        for i,af in enumerate(airfoils):
            af.compute(0.0, 0.5)
            print 'AF '+str(i)+' data = ',af.get_Lref(), af.get_Sref(), af.Cl, af.Cdw, af.Cdvp, af.Cdf, af.pcop
