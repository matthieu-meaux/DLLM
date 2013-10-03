# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: EADS
# @version: 0.1
# @author: Matthieu Meaux

# - Local imports -
import sys
from padge.PCADEngine.Base.PCADModel import PCADModel
from DLLM.polarManager.analyticAirfoil import AnalyticAirfoil
from DLLM.polarManager.airfoilPolar import AirfoilPolar
from numpy import zeros, array
from numpy import pi, sqrt, cos,sin

class Wing_param():
    
    ERROR_MSG = 'ERROR in Wing_param.'
    POSSIBLE_GEOM_TYPES=["Rectangular","Elliptic","Broken"]
    
    def __init__(self, tag, geom_type='Broken', n_sect=20):
        """
        Constructor: set main attributes
        """
        self.__tag        = tag
        self.__geom_type  = None
        self.__n_sect     = None
        self.__PCADModel  = None
        self.__BC_manager = None
        
        self.set_geom_type(geom_type)
        self.set_n_sect(n_sect)
        
        self.__PCADModel  = PCADModel(self.__tag)
        self.__BC_manager = self.__PCADModel.get_BC_manager()
        self.__ndv        = 0
        
        self.__init_discrete_attributes()
        self.__init_airfoils_attributes()
        
        
    def __init_discrete_attributes(self):
        # -- discrete attributes
        self.__span               = None
        self.__span_grad          = None
        
        self.__sweep              = None
        self.__sweep_grad         = None
        
        self.__break_percent      = None
        self.__break_percent_grad = None
        
        self.__root_chord         = None
        self.__root_chord_grad    = None
        
        self.__break_chord        = None
        self.__break_chord_grad   = None
        
        self.__tip_chord          = None
        self.__tip_chord_grad     = None
        
        self.__root_height        = None
        self.__root_height_grad   = None
        
        self.__break_height       = None
        self.__break_height_grad  = None
        
        self.__tip_height         = None
        self.__tip_height_grad    = None
        
        self.__twist              = None
        self.__twist_grad         = None
        
        self.__chords             = None
        self.__chords_grad        = None
        
        self.__rel_thicks         = None
        self.__rel_thicks_grad    = None
        
        self.__XYZ                = None
        self.__XYZ_grad           = None
        
        self.__eta                = None
        self.__eta_grad           = None
        
    def __init_airfoils_attributes(self):
        self.__ref_airfoil     = None   # reference airfoil, only used if the same airfoil is put for all sections
        self.__airfoils        = None   # Airfoil list for each section
        self.__linked_airfoils = None   # Airfoil scaled to the the planform
        
    # -- Accessors
    def get_n_sect(self):
        return self.__n_sect
    
    def get_ndv(self):
        return self.__ndv
    
    def get_eta(self):
        return self.__eta
    
    def get_eta_grad(self):
        return self.__eta_grad
    
    def get_XYZ(self):
        return self.__XYZ
    
    def get_XYZ_grad(self):
        return self.__XYZ_grad
    
    def get_twist(self):
        return self.__twist
    
    def get_twist_grad(self):
        return self.__twist_grad
    
    def get_chords(self):
        return self.__chords
    
    def get_chords_grad(self):
        return self.__chords_grad
    
    def get_rel_thicks(self):
        return self.__rel_thicks
    
    def get_rel_thicks_grad(self):
        return self.__rel_thicks_grad
        
    # -- Setters
    def set_geom_type(self, geom_type):
        """
        @param wing_geometry_type : Rectangular or Elliptic planforms
        @param type wing_geometry_type : String
        """
        if not geom_type in self.POSSIBLE_GEOM_TYPES:
            raise Exception, "geom_type :"+str(geom_type)+" not in possible types : "+str(self.POSSIBLE_GEOM_TYPES)
        self.__geom_type = geom_type
        if self.__BC_manager is not None:
            self.reinitialize()
        
    def set_n_sect(self, n_sect):
        ERROR_MSG=self.ERROR_MSG+'.set_n_sect: '+str(self.__tag)+': '
        if n_sect%2!=0:
            raise Exception, ERROR_MSG+'The total number of elements in the wing must be even.'
        self.__n_sect = n_sect
        if self.__BC_manager is not None:
            self.reinitialize()
            
    def set_value(self, Id, val):
        pt=self.__BC_manager.get_pt(Id)
        pt.set_value(val)
        
    # -- controller management
    def convert_to_variable(self,Id):
        self.__BC_manager.convert_to_variable(Id)
        
    def convert_to_design_variable(self, Id, lbnd, ubnd):
        self.__BC_manager.convert_to_design_variable(Id, lbnd, ubnd)
        
    def convert_to_parameter(self, Id, fexpr):
        self.__BC_manager.convert_to_parameter(Id, fexpr)
        
    def get_dv_array(self):
        return self.__BC_manager.get_dv_array()
    
    def get_dv_id_list(self):
        return self.__BC_manager.get_dv_id_list()
        
    def build_wing(self):
        self.__BC_manager.clean()
        self.__BC_manager.create_variable(self.__tag+'.span',0.)
        if self.__geom_type not in ['Elliptic']:
            self.__BC_manager.create_variable(self.__tag+'.sweep',0.)
        if   self.__geom_type in ['Rectangular','Elliptic']:
            self.__BC_manager.create_variable(self.__tag+'.root_chord',0.)
            self.__BC_manager.create_variable(self.__tag+'.root_height',0.)
            self.__BC_manager.create_variable(self.__tag+'.tip_height',0.)
        elif self.__geom_type == 'Broken':
            self.__BC_manager.create_variable(self.__tag+'.break_percent',0.33)
            self.__BC_manager.create_variable(self.__tag+'.root_chord',0.)
            self.__BC_manager.create_variable(self.__tag+'.break_chord',0.)
            self.__BC_manager.create_variable(self.__tag+'.tip_chord',0.)
            self.__BC_manager.create_variable(self.__tag+'.root_height',0.)
            self.__BC_manager.create_variable(self.__tag+'.break_height',0.)
            self.__BC_manager.create_variable(self.__tag+'.tip_height',0.)
            
        for i in xrange(self.__n_sect):
            self.__BC_manager.create_design_variable(self.__tag+'.twist'+str(i),-25.,0.,25.)
    
    def update(self):
        self.__PCADModel.update()
        self.__ndv = self.__BC_manager.get_ndv()
        self.__build_discretization()
        self.__check_airfoils_inputs()
        self.__link_airfoils_to_geom()
        
    def update_from_x_list(self,x):
        #print 'update wing_param with x=',x
        self.__BC_manager.update_dv_from_x_list(x)
        self.update()
   
    def __repr__(self):
        info_string ='\n*** Wing param information ***'
        info_string+='\n  geom_type : '+str(self.__geom_type)
        info_string+='\n  n_sect    : '+str(self.__n_sect)
        info_string+='\n  ndv       : '+str(self.__ndv)
        info_string+='\n  --    parameters information section    --\n'
        for Id in self.__BC_manager.get_list_id():
            pt=self.__BC_manager.get_pt(Id)
            BC_Type = pt.get_BCType()
            if BC_Type == 'Variable':
                value   = pt.get_value()
                info_string+="%30s"%Id+"%20s"%BC_Type+" %24.16e"%value+"\n"
            if BC_Type == 'DesignVariable':
                value   = pt.get_value()
                bounds=pt.get_bounds()
                info_string+="%30s"%Id+"%20s"%BC_Type+" %24.16e"%value+" %24s"%str(bounds)+"\n"
        info_string+='  -- end of parameters information section --\n'
        return info_string
    
    # -- discretization methods
    def __build_discretization(self):
        self.__build_planform()
        self.__build_chords()
        self.__build_rel_thicks()
        self.__build_XYZ()
        self.__build_eta()
        

    def __build_planform(self):
        deg_to_rad = pi/180.
        
        # -- span
        Id   = self.__tag+'.span'
        pt   = self.__BC_manager.get_pt(Id)
        val  = pt.get_value()
        grad = pt.get_gradient()
        self.__span      = val
        self.__span_grad = grad
        
        if self.__geom_type not in ['Elliptic']:
            # -- sweep
            Id   = self.__tag+'.sweep'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__sweep      = val*deg_to_rad
            self.__sweep_grad = grad*deg_to_rad
        
        if   self.__geom_type in ['Rectangular','Elliptic']:
            Id   = self.__tag+'.root_chord'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__root_chord      = val
            self.__root_chord_grad = grad
            
            Id   = self.__tag+'.root_height'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__root_height      = val
            self.__root_height_grad = grad

            Id   = self.__tag+'.tip_height'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_height      = val
            self.__tip_height_grad = grad

        elif self.__geom_type == 'Broken':
            Id   = self.__tag+'.break_percent'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__break_percent      = val
            self.__break_percent_grad = grad
            
            Id   = self.__tag+'.root_chord'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__root_chord      = val
            self.__root_chord_grad = grad
            
            Id   = self.__tag+'.break_chord'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__break_chord      = val
            self.__break_chord_grad = grad
            
            Id   = self.__tag+'.tip_chord'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_chord      = val
            self.__tip_chord_grad = grad
            
            Id   = self.__tag+'.root_height'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__root_height      = val
            self.__root_height_grad = grad
            
            Id   = self.__tag+'.break_height'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__break_height      = val
            self.__break_height_grad = grad

            Id   = self.__tag+'.tip_height'
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__tip_height      = val
            self.__tip_height_grad = grad
            
        self.__twist      = zeros((self.__n_sect))
        self.__twist_grad = zeros((self.__n_sect,self.__ndv))
        for i in xrange(self.__n_sect):
            Id   = self.__tag+'.twist'+str(i)
            pt   = self.__BC_manager.get_pt(Id)
            val  = pt.get_value()
            grad = pt.get_gradient()
            self.__twist[i]        = val*deg_to_rad
            self.__twist_grad[i,:] = grad*deg_to_rad
            
    def __build_chords(self):

        self.__chords      = zeros((self.__n_sect))
        self.__chords_grad = zeros((self.__n_sect,self.__ndv))
        
        N = self.__n_sect
        n = N/2
        
        if   self.__geom_type == 'Elliptic':
            for i in xrange(N):
                r=float(i+0.5-n)/float(n)
                self.__chords[i]        = self.__root_chord*sqrt(1.-r**2)
                self.__chords_grad[i,:] = self.__root_chord_grad[:]*sqrt(1.-r**2)

        elif self.__geom_type == 'Rectangular':
            for i in xrange(N):
                self.__chords[i]        = self.__root_chord
                self.__chords_grad[i,:] = self.__root_chord_grad[:]

        elif self.__geom_type == 'Broken':
            p      = self.__break_percent/100.
            p_grad = self.__break_percent_grad/100.
            
            for i in xrange(N):
                r = abs(float(i+0.5-n)/float(n))
                
                if r <= p:
                    coeff = r/p
                    dcoeff = -r*p_grad[:]/p**2
                    self.__chords[i]        = (self.__break_chord - self.__root_chord)*coeff + self.__root_chord
                    self.__chords_grad[i,:] = (self.__break_chord_grad[:] - self.__root_chord_grad[:])*coeff  \
                                            + (self.__break_chord - self.__root_chord)*dcoeff \
                                            + self.__root_chord_grad[:]
                    
                else:
                    coeff = (r-p)/(1.-p)
                    dcoeff = (r-1)*p_grad[:]/(1.-p)**2
                    self.__chords[i]        = (self.__tip_chord - self.__break_chord)*coeff + self.__break_chord
                    self.__chords_grad[i,:] = (self.__tip_chord_grad[:] - self.__break_chord_grad[:])*coeff  \
                                            + (self.__tip_chord - self.__break_chord)*dcoeff \
                                            + self.__break_chord_grad[:]
                    
    def __build_rel_thicks(self):
        N = self.__n_sect
        n = N/2
        
        self.__rel_thicks      = zeros((self.__n_sect))
        self.__rel_thicks_grad = zeros((self.__n_sect,self.__ndv))
        
        if self.__geom_type == 'Broken':
            p      = self.__break_percent/100.
            p_grad = self.__break_percent_grad/100.
            for i in xrange(N):
                r=abs(float(i+0.5-n)/float(n))
                if r <= p:
                    coeff = r/p
                    term = ((self.__break_height - self.__root_height)*coeff + self.__root_height)
                    self.__rel_thicks[i] = term / self.__chords[i]

                    dcoeff = -r*p_grad[:]/p**2
                    term_grad = (self.__break_height_grad[:] - self.__root_height_grad[:])*coeff \
                              + (self.__break_height - self.__root_height)*dcoeff \
                              + self.__root_height_grad[:]
                    self.__rel_thicks_grad[i,:] = (term_grad[:]*self.__chords[i]-term*self.__chords_grad[i,:])/(self.__chords[i])**2
                
                else:
                    coeff = (r-p)/(1.-p)
                    term = ((self.__tip_height - self.__break_height)*coeff + self.__break_height)
                    self.__rel_thicks[i]  =  term / self.__chords[i]
                    
                    dcoeff = (r-1)*p_grad[:]/(1.-p)**2
                    term_grad = (self.__tip_height_grad[:] - self.__break_height_grad[:])*coeff \
                               + (self.__tip_height- self.__break_height)*dcoeff \
                               + self.__break_height_grad[:]
                    self.__rel_thicks_grad[i,:] = (term_grad[:]*self.__chords[i]-term*self.__chords_grad[i,:])/(self.__chords[i])**2         
                               
        else:
            for i in xrange(N):
                r=abs((i+0.5-n)/float(n)) 
                self.__rel_thicks[i]        = ((self.__tip_height - self.__root_height)*r + self.__root_height) / self.__chords[i]
                self.__rel_thicks_grad[i,:] = ( ((self.__tip_height_grad[:] - self.__root_height_grad[:])*r + self.__root_height_grad[:])*self.__chords[i] \
                                            - ((self.__tip_height - self.__root_height)*r + self.__root_height) * self.__chords_grad[i,:] ) / (self.__chords[i])**2

    def __build_XYZ(self):
        N = self.__n_sect
        n = N/2
        
        self.__XYZ      = zeros((3,N))
        self.__XYZ_grad = zeros((3,N,self.__ndv))
        
        for i in xrange(N):
            r=float(i+0.5-n)/float(n)
            abs_r=abs(r)
            self.__XYZ[0,i]        = 0.25*self.__chords[i] + abs_r*self.__span*cos(self.__sweep)/2.
            self.__XYZ_grad[0,i,:] = 0.25*self.__chords_grad[i] + abs_r*self.__span_grad[:]*cos(self.__sweep)/2.-abs_r*self.__span*sin(self.__sweep)*self.__sweep_grad[:]/2.
            
            self.__XYZ[1,i]        = r*self.__span/2.
            self.__XYZ_grad[1,i,:] = r*self.__span_grad[:]/2.
        
            self.__XYZ[2,i]        = 0.
            self.__XYZ_grad[2,i,:] = zeros(self.__ndv)
        
    def __build_eta(self):
        N = self.__n_sect
        n = N/2
        
        self.__eta      = zeros((3,N+1))
        self.__eta_grad = zeros((3,N+1,self.__ndv))
        
        for i in xrange(N+1):
            ifl=float(i-n)
            self.__eta[1,i]      = (ifl/float(n))*self.__span/2.
            self.__eta_grad[1,i] = (ifl/float(n))*self.__span_grad[:]/2.
            
    # -- Airfoils methods 
    def get_linked_airfoils(self):
        return self.__linked_airfoils
    
    def set_ref_aifoil(self, ref_airfoil):
        self.__ref_airfoil = ref_airfoil
        
    def set_airfoils(self, airfoils):
        ERROR_MSG=self.ERROR_MSG+'set_airfoils: '
        if len(airfoils)!=self.__n_sect:
            print ERROR_MSG+'number of airfoils is different than the number of sections. attribute not set'
            self.__airfoils = None
        else:
            self.__airfoils = airfoils
            
    def build_linear_airfoil(self, AoA0=0., Cd0=0., Cm0=0., Sref=1., Lref=1., rel_thick=0., sweep=0., set_as_ref=True):
        degToRad = pi/180.
        airfoil  = AnalyticAirfoil(AoA0=degToRad*AoA0, Cd0=Cd0, Cm0=Cm0, Sref=Sref, Lref=Lref, rel_thick=rel_thick, sweep=sweep)
        if set_as_ref:
            self.set_ref_aifoil(airfoil)
        return airfoil
    
    def build_polar_airoil(self, database, Sref=1., Lref=1., interpolator='2DSpline', set_as_ref=True):
        # Why relative thickness usage ? The extraction from a polar should give us more freedom.
        airfoil = AirfoilPolar(database,rel_thick=0.15, interpolator=interpolator, Sref=Sref, Lref=Lref)
        if set_as_ref:
            self.set_ref_aifoil(airfoil)
        return airfoil

    def build_airfoils_from_ref(self):
        ERROR_MSG=self.ERROR_MSG+'build_airfoils_from_ref: '
        if self.__ref_airfoil is None:
            print ERROR_MSG+'cannot build airfoils if reference airfoil is not defined.'
        else:
            airfoils=[]
            for n in xrange(self.__n_sect):
                airfoils.append(self.__ref_airfoil.get_scaled_copy())
            
            self.set_airfoils(airfoils)
            
    def __link_airfoils_to_geom(self):
        self.__linked_airfoils=[]
        for i in range(self.__n_sect):
            LLoc,LLoc_grad,SLoc,SLoc_grad=self.__compute_local_info(i)
            linked_af = self.__airfoils[i].get_scaled_copy(SLoc, LLoc)
            linked_af.set_Sref_grad(SLoc_grad)
            linked_af.set_Lref_grad(LLoc_grad)
            linked_af.set_rel_thick_grad(self.__rel_thicks_grad[i])
            linked_af.set_rel_thick(self.__rel_thicks[i])
            linked_af.set_sweep(self.__sweep)
            linked_af.set_sweep_grad(self.__sweep_grad)
            self.__linked_airfoils.append(linked_af)
        
    def __compute_local_info(self,i):
        LLoc      = self.__chords[i]
        LLoc_grad = self.__chords_grad[i]
        SLoc      = self.__span * LLoc / float(self.__n_sect)
        SLoc_grad = (self.__span_grad * LLoc + self.__span * LLoc_grad) / float(self.__n_sect)
        return LLoc, LLoc_grad, SLoc, SLoc_grad
    
    def __check_airfoils_inputs(self):
        ERROR_MSG=self.ERROR_MSG+'__check_airfoils_inputs: '+str(self.__tag)+': '
        checked=True
        if self.__airfoils is None:
            checked=False
            print ERROR_MSG+'airfoils attribute undefined. please set airfoils attribute.'
        elif len(self.__airfoils)!=self.__n_sect:
            checked=False
            print ERROR_MSG+'number of airfoils must equal to the number of geometrical sections.'
        
        if not checked:
            sys.exit(1)