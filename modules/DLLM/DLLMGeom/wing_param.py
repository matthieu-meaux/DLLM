# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus Group Innovations
# @version: 0.1
# @author: Matthieu Meaux

# - Local imports -
import sys
import string
import numpy
from DLLM.Controllers.Kernel.BCManager import BCManager
from DLLM.polarManager.analyticAirfoil import AnalyticAirfoil
from DLLM.polarManager.airfoilPolar import AirfoilPolar
from DLLM.polarManager.MetaAirfoil import MetaAirfoil
from numpy import zeros, array
from numpy import pi, sqrt, cos,sin
from copy import deepcopy

class Wing_param():
    
    ERROR_MSG = 'ERROR in Wing_param.'
    POSSIBLE_GEOM_TYPES=["Rectangular","Elliptic","Broken"]
    POS_DISTRIB=['linear','cos_law']
    
    def __init__(self, tag, geom_type='Broken', n_sect=20):
        """
        Constructor: set main attributes
        """
        self.__tag        = tag
        self.__geom_type  = None
        self.__n_sect     = None
        self.__PCADModel  = None
        self.__BC_manager = None
        
        self.__distrib_type = 'cos_law'
        
        self.set_geom_type(geom_type)
        self.set_n_sect(n_sect)
        
        self.__BC_manager = BCManager()

        self.__ndv        = 0
        
        self.__thetaY     = None
        
        self.__AoA_id     = 'AoA'
        
        self.__init_discrete_attributes()
        self.__init_airfoils_attributes()
        
        
    def __init_discrete_attributes(self):
        # -- OC design variable 
        self.__AoA                = None
        self.__AoA_grad           = None
        
        # -- discrete attributes
        self.__span               = None
        self.__span_grad          = None
        
        self.__r_list_y           = None
        self.__r_list_eta         = None
        
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
        
        self.__heights            = None
        self.__heights_grad       = None
        
        self.__rel_thicks         = None
        self.__rel_thicks_grad    = None
        
        self.__XYZ                = None
        self.__XYZ_grad           = None
        
        self.__eta                = None
        self.__eta_grad           = None
        
        self.__Lref               = None
        self.__Lref_grad          = None
        
        self.__Sref               = None
        self.__Sref_grad          = None
        
        self.__root_ToC           = None
        self.__root_ToC_grad      = None
        
        self.__break_ToC          = None
        self.__break_ToC_grad     = None
        
        self.__tip_ToC            = None
        self.__tip_ToC_grad       = None
        
        self.__AR                 = None
        self.__AR_grad            = None
        
        self.__fuel               = None
        self.__fuel_grad          = None
        
    def __init_airfoils_attributes(self):
        self.__airfoil_type    = None
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
    
    def get_thetaY(self):
        return self.__thetaY
    
    def get_AoA(self):
        return self.__AoA
    
    def get_AoA_grad(self):
        return self.__AoA_grad
    
    def get_span(self):
        return self.__span
    
    def get_span_grad(self):
        return self.__span_grad
    
    def get_sweep(self):
        return self.__sweep
    
    def get_sweep_grad(self):
        return self.__sweep_grad
    
    def get_Lref(self):
        return self.__Lref
    
    def get_Lref_grad(self):
        return self.__Lref_grad
    
    def get_Sref(self):
        return self.__Sref
    
    def get_Sref_grad(self):
        return self.__Sref_grad
    
    def get_root_ToC(self):
        return self.__root_ToC
    
    def get_root_ToC_grad(self):
        return self.__root_ToC_grad
    
    def get_break_ToC(self):
        return self.__break_ToC
    
    def get_break_ToC_grad(self):
        return self.__break_ToC_grad
    
    def get_tip_ToC(self):
        return self.__tip_ToC
    
    def get_tip_ToC_grad(self):
        return self.__tip_ToC_grad
    
    def get_AR(self):
        return self.__AR
    
    def get_AR_grad(self):
        return self.__AR_grad
    
    def get_fuel(self):
        return self.__fuel
    
    def get_fuel_grad(self):
        return self.__fuel_grad
        
    # -- Setters
    def set_AoA_id(self, AoA_id):
        self.__AoA_id = AoA_id
        
    def set_geom_type(self, geom_type):
        """
        @param wing_geometry_type : Rectangular or Elliptic planforms
        @param type wing_geometry_type : String
        """
        if not geom_type in self.POSSIBLE_GEOM_TYPES:
            raise Exception, "geom_type :"+str(geom_type)+" not in possible types : "+str(self.POSSIBLE_GEOM_TYPES)
        self.__geom_type = geom_type
        
    def set_n_sect(self, n_sect):
        ERROR_MSG=self.ERROR_MSG+'.set_n_sect: '+str(self.__tag)+': '
        if n_sect%2!=0:
            raise Exception, ERROR_MSG+'The total number of elements in the wing must be even.'
        self.__n_sect = n_sect
        
    def set_distrib_type(self, distrib_type):
        if not distrib_type in self.POS_DISTRIB:
            raise Exception, "distrib_type :"+str(distrib_type)+" not in possible types : "+str(self.POS_DISTRIB)
        self.__distrib_type = distrib_type
            
    def set_value(self, Id, val):
        pt=self.__BC_manager.get_pt(Id)
        pt.set_value(val)
        
    def set_thetaY(self, thetaY):
        self.__thetaY = thetaY
        
    # -- controller management
    def convert_to_variable(self,Id):
        self.__BC_manager.convert_to_variable(Id)
        
    def convert_to_design_variable(self, Id, bounds):
        self.__BC_manager.convert_to_design_variable(Id, bounds)
        
    def convert_to_parameter(self, Id, fexpr):
        self.__BC_manager.convert_to_parameter(Id, fexpr)
        
    def add_AoA_design_variable(self, val, bounds):
        self.__BC_manager.create_design_variable(self.__tag+'.'+self.__AoA_id, val, bounds)
        
    def __update_AoA(self):
        Id   = self.__tag+'.'+self.__AoA_id
        if Id in self.get_dv_id_list():
            pt   = self.__BC_manager.get_pt(Id)
            deg_to_rad = pi/180.
            self.__AoA      = pt.get_value()*deg_to_rad
            self.__AoA_grad = pt.get_gradient()*deg_to_rad
        
    def get_dv_array(self):
        return self.__BC_manager.get_dv_array()
    
    def get_dv_id_list(self):
        return self.__BC_manager.get_dv_id_list()
    
    def get_bounds_array(self):
        return self.__BC_manager.get_bounds_array()
        
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
            
        for i in xrange(self.__n_sect/2):
            self.__BC_manager.create_design_variable(self.__tag+'.rtwist'+str(i),0.,(-10.,10.))
        for i in xrange(self.__n_sect/2):
#             print self.__tag+'.twist'+str(i),self.__tag+'.twist'+str(self.__n_sect/2+i),self.__tag+'.rtwist'+str(self.__n_sect/2-1-i),self.__tag+'.rtwist'+str(i)
            self.__BC_manager.create_parameter(self.__tag+'.twist'+str(i),self.__tag+'.rtwist'+str(self.__n_sect/2-1-i))
            self.__BC_manager.create_parameter(self.__tag+'.twist'+str(self.__n_sect/2+i),self.__tag+'.rtwist'+str(i))
#         for i in xrange(self.__n_sect):
#             self.__BC_manager.create_design_variable(self.__tag+'.twist'+str(i),-25.,0.,25.)
    
    def config_from_dict(self, OC, config_dict):
        self.build_wing()
        self.__config_airfoils(OC, config_dict)
        self.__config_desc(config_dict)
        self.update()
        
    def __config_airfoils(self, OC, config_dict):
        in_keys_list=config_dict.keys()
        
        airfoil_type_key=self.__tag+'.airfoil.type'
        if airfoil_type_key in in_keys_list:
            airfoil_type=config_dict[airfoil_type_key]
        else:
            airfoil_type='simple'
        
        if airfoil_type == 'simple':
            AoA0_key=self.__tag+'.airfoil.AoA0'
            if AoA0_key in in_keys_list:
                AoA0=config_dict[AoA0_key]
            else:
                AoA0=0.
            
            Cm0_key=self.__tag+'.airfoil.Cm0'
            if Cm0_key in in_keys_list:
                Cm0=config_dict[Cm0_key]
            else:
                Cm0=0.
        
            self.build_linear_airfoil(OC, AoA0=AoA0, Cm0=Cm0, set_as_ref=True)
            
        elif airfoil_type == 'meta':
            surrogate_model_key=self.__tag+'.airfoil.surrogate_model'
            if surrogate_model_key in in_keys_list:
                surrogate_model = config_dict[surrogate_model_key]
            else:
                surrogate_model = None
                
            self.build_meta_airfoil(OC, surrogate_model, set_as_ref=True)
            
        self.build_airfoils_from_ref()
        
    def __config_desc(self, config_dict):
        in_keys_list=sorted(config_dict.keys())
        existing_keys=deepcopy(self.__BC_manager.get_list_id())
        
        # Convert pre-defined parameters
        for existing_id in existing_keys:
            ex_words=existing_id.split('.')
            name=ex_words[-1]
            Id_in=self.__tag+'.desc.'+name
            Id=self.__tag+'.'+name
            if Id_in+'.type' in in_keys_list:
                type=config_dict[Id_in+'.type']
                if   type == 'DesignVariable':
                    bounds = config_dict[Id_in+'.bounds']
                    value  = config_dict[Id_in+'.value']
                    self.set_value(Id,value)
                    self.convert_to_design_variable(Id,bounds)
                elif type == 'Variable':
                    value  = config_dict[Id_in+'.value']
                    self.set_value(Id,value)
                    self.convert_to_variable(Id)
                elif type == 'Parameter':
                    fexpr = config_dict[Id_in+'.fexpr']
                    self.convert_to_parameter(Id,fexpr)
        
        # Add user defined DesignVariable, Variable and Parameter
        added_list=[]        
        for in_key in in_keys_list:
            words=in_key.split('.')
            if len(words) >=4:
                test=string.join(words[:-2],'.')
                if test==self.__tag+'.desc':
                    name=words[-2]
                    Id=self.__tag+'.'+name
                    Id_in=self.__tag+'.desc.'+name
                    type=config_dict[Id_in+'.type']
                    if Id not in existing_keys and Id not in added_list:
                        if   type == 'DesignVariable':
                            bounds = config_dict[Id_in+'.bounds']
                            value  = config_dict[Id_in+'.value']
                            self.__BC_manager.create_design_variable(Id,value,(bounds[0],bounds[1]))
                        elif type == 'Variable':
                            value  = config_dict[Id_in+'.value']
                            self.__BC_manager.create_variable(Id,value)
                        elif type == 'Parameter':
                            fexpr = config_dict[Id_in+'.fexpr']
                            self.__BC_manager.create_parameter(Id,fexpr)
                        added_list.append(Id)
        
    def update(self):
        self.__BC_manager.update()
        self.__ndv = self.__BC_manager.get_ndv()
        self.__update_AoA()
        self.__check_thetaY()
        self.__build_discretization()
        self.__check_airfoils_inputs()
        self.__link_airfoils_to_geom()
        self.__compute_Sref_Lref_ToC_AR_fuel()
        
    def update_from_x_list(self,x):
        #print 'update wing_param with x=',x
        self.__BC_manager.update_dv_from_x_list(x)
        self.update()
   
    def __repr__(self):
        info_string ='\n*** Wing param information ***'
        info_string+='\n  geom_type    : '+str(self.__geom_type)
        info_string+='\n  n_sect       : '+str(self.__n_sect)
        info_string+='\n  ndv          : '+str(self.__ndv)
        info_string+='\n  airfoil_type : '+str(self.__airfoil_type)
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
            if BC_Type == 'Parameter':
                expr  = pt.get_expr()
                info_string+="%30s"%Id+"%20s"%BC_Type+" %24s"%expr+"\n"
        info_string+='  -- end of parameters information section --\n'
        return info_string
    
    # -- Structural displacement methods
    def __check_thetaY(self):
        # 6 dimensions for structural displacement: dx,dy,dz,dthetax,dthetay,dthetaz at each section
        if self.__thetaY is None:
            self.__thetaY = zeros(self.__n_sect) 
    
    # -- discretization methods
    def __build_discretization(self):
        self.__build_planform()
        self.__build_r_lists()
        self.__build_chords()
        self.__build_heights_rel_thicks()
        self.__build_XYZ_eta()
        
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
        else:
            self.__sweep      = 0.
            self.__sweep_grad = zeros(self.get_ndv())
        
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
    
    def __build_r_lists(self):
        N = self.__n_sect
        
        if self.__distrib_type == 'linear':
            self.__r_list_eta = numpy.linspace(-0.5,0.5,N+1)
        else:
            theta_list = numpy.linspace(0.,numpy.pi,N+1)
            self.__r_list_eta = -0.5+0.5*(1.-numpy.cos(theta_list))
        
        self.__r_list_y   = 0.5*(self.__r_list_eta[:-1]+self.__r_list_eta[1:])
        
    def __build_chords(self):
        self.__chords      = zeros((self.__n_sect))
        self.__chords_grad = zeros((self.__n_sect,self.__ndv))
        
        self.__chords_eta      = zeros((self.__n_sect+1))
        self.__chords_grad_eta = zeros((self.__n_sect+1,self.__ndv))
               
        if   self.__geom_type == 'Elliptic':
            for i,r in enumerate(self.__r_list_y):
                self.__chords[i]        = self.__root_chord*sqrt(1.-(2.*r)**2)
                self.__chords_grad[i,:] = self.__root_chord_grad[:]*sqrt(1.-(2.*r)**2)
            for i,r in enumerate(self.__r_list_eta):
                self.__chords_eta[i]        = self.__root_chord*sqrt(1.-(2.*r)**2)
                self.__chords_grad_eta[i,:] = self.__root_chord_grad[:]*sqrt(1.-(2.*r)**2)

        elif self.__geom_type == 'Rectangular':
            for i in xrange(self.__n_sect):
                self.__chords[i]        = self.__root_chord
                self.__chords_grad[i,:] = self.__root_chord_grad[:]
            for i in xrange(self.__n_sect+1):
                self.__chords_eta[i]        = self.__root_chord
                self.__chords_grad_eta[i,:] = self.__root_chord_grad[:]
                

        elif self.__geom_type == 'Broken':
            p      = self.__break_percent/100.
            p_grad = self.__break_percent_grad/100.
            
            for i,r in enumerate(self.__r_list_y):
                r=abs(2.*r)
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
                                            
            for i,r in enumerate(self.__r_list_eta):
                r=abs(2.*r)
                if r <= p:
                    coeff = r/p
                    dcoeff = -r*p_grad[:]/p**2
                    self.__chords_eta[i]        = (self.__break_chord - self.__root_chord)*coeff + self.__root_chord
                    self.__chords_grad_eta[i,:] = (self.__break_chord_grad[:] - self.__root_chord_grad[:])*coeff  \
                                                + (self.__break_chord - self.__root_chord)*dcoeff \
                                                + self.__root_chord_grad[:]
                    
                else:
                    coeff = (r-p)/(1.-p)
                    dcoeff = (r-1)*p_grad[:]/(1.-p)**2
                    self.__chords_eta[i]        = (self.__tip_chord - self.__break_chord)*coeff + self.__break_chord
                    self.__chords_grad_eta[i,:] = (self.__tip_chord_grad[:] - self.__break_chord_grad[:])*coeff  \
                                                + (self.__tip_chord - self.__break_chord)*dcoeff \
                                                + self.__break_chord_grad[:]
    
    def __build_heights_rel_thicks(self):        
        self.__heights         = zeros((self.__n_sect))
        self.__heights_grad    = zeros((self.__n_sect,self.__ndv))
        
        self.__rel_thicks      = zeros((self.__n_sect))
        self.__rel_thicks_grad = zeros((self.__n_sect,self.__ndv))
        
        if self.__geom_type == 'Broken':
            p      = self.__break_percent/100.
            p_grad = self.__break_percent_grad/100.
            
            for i,r in enumerate(self.__r_list_y):
                r=abs(2.*r)
                if r <= p:
                    coeff = r/p
                    dcoeff = -r*p_grad[:]/p**2
                    self.__heights[i]        = ((self.__break_height - self.__root_height)*coeff + self.__root_height)
                    self.__heights_grad[i,:] = (self.__break_height_grad[:] - self.__root_height_grad[:])*coeff \
                                            + (self.__break_height - self.__root_height)*dcoeff \
                                             + self.__root_height_grad[:]
                    
                else:
                    coeff = (r-p)/(1.-p)
                    dcoeff = (r-1)*p_grad[:]/(1.-p)**2
                    self.__heights[i]        = ((self.__tip_height - self.__break_height)*coeff + self.__break_height)
                    self.__heights_grad[i,:] = (self.__tip_height_grad[:] - self.__break_height_grad[:])*coeff \
                                             + (self.__tip_height- self.__break_height)*dcoeff \
                                             + self.__break_height_grad[:]
        else:
            for i,r in enumerate(self.__r_list_y):
                r=abs(2.*r)
                self.__heights[i]        = ((self.__tip_height - self.__root_height)*r + self.__root_height)
                self.__heights_grad[i,:] = ((self.__tip_height_grad[:] - self.__root_height_grad[:])*r + self.__root_height_grad[:])
                
        #-- build rel_thiks
        for i in xrange(self.__n_sect):
            self.__rel_thicks[i]         =  self.__heights[i] / self.__chords[i]
            self.__rel_thicks_grad[i,:]  = (self.__heights_grad[i,:]*self.__chords[i]-self.__heights[i]*self.__chords_grad[i,:])/(self.__chords[i])**2 
            
    def __build_XYZ_eta(self):       
        self.__XYZ      = zeros((3,self.__n_sect))
        self.__XYZ_grad = zeros((3,self.__n_sect,self.__ndv))
        
        self.__eta      = zeros((3,self.__n_sect+1))
        self.__eta_grad = zeros((3,self.__n_sect+1,self.__ndv))
        
        #print 'sweep=',self.__sweep,'sin sweep',sin(self.__sweep)
        
        for i,r in enumerate(self.__r_list_y):
            abs_r=abs(r)
            #self.__XYZ[0,i]        = 0.25*self.__chords[i] + abs_r*self.__span*sin(self.__sweep)
            self.__XYZ[0,i]        = abs_r*self.__span*sin(self.__sweep)
            #self.__XYZ_grad[0,i,:] = 0.25*self.__chords_grad[i] + abs_r*self.__span_grad[:]*sin(self.__sweep)+abs_r*self.__span*cos(self.__sweep)*self.__sweep_grad[:]
            self.__XYZ_grad[0,i,:] = abs_r*self.__span_grad[:]*sin(self.__sweep)+abs_r*self.__span*cos(self.__sweep)*self.__sweep_grad[:]
            
            self.__XYZ[1,i]        = r*self.__span
            self.__XYZ_grad[1,i,:] = r*self.__span_grad[:]
        
        for i,r in enumerate(self.__r_list_eta):
            abs_r=abs(r)
            #self.__eta[0,i]        = 0.25*self.__chords_eta[i] + abs_r*self.__span*sin(self.__sweep)
            self.__eta[0,i]        = abs_r*self.__span*sin(self.__sweep)
            #self.__eta_grad[0,i,:] = 0.25*self.__chords_grad_eta[i] + abs_r*self.__span_grad[:]*sin(self.__sweep)+abs_r*self.__span*cos(self.__sweep)*self.__sweep_grad[:]
            self.__eta_grad[0,i,:] = abs_r*self.__span_grad[:]*sin(self.__sweep)+abs_r*self.__span*cos(self.__sweep)*self.__sweep_grad[:]
            
            self.__eta[1,i]        = r*self.__span 
            self.__eta_grad[1,i,:] = r*self.__span_grad[:]

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
            
    def build_linear_airfoil(self, OC, AoA0=0., Cm0=0., Sref=1., Lref=1., rel_thick=0., sweep=0., Ka=0.95, set_as_ref=True):
        self.__airfoil_type = 'simple'
        degToRad = pi/180.
        airfoil  = AnalyticAirfoil(OC, AoA0=degToRad*AoA0, Cm0=Cm0, Sref=Sref, Lref=Lref, rel_thick=rel_thick, sweep=sweep, Ka=Ka)
        if set_as_ref:
            self.set_ref_aifoil(airfoil)
        return airfoil
    
#     def build_polar_airoil(self, OC, database, Sref=1., Lref=1., interpolator='2DSpline', set_as_ref=True):
#         # Why relative thickness usage ? The extraction from a polar should give us more freedom.
#         airfoil = AirfoilPolar(OC, database,rel_thick=0.15, interpolator=interpolator, Sref=Sref, Lref=Lref)
#         if set_as_ref:
#             self.set_ref_aifoil(airfoil)
#         return airfoil
        
    def build_meta_airfoil(self, OC, surrogate_model, relative_thickness=.12, camber=0., Sref=1., Lref=1., sweep=.0, set_as_ref=True):
        self.__airfoil_type = 'meta'
        airfoil = MetaAirfoil(OC, surrogate_model, relative_thickness=relative_thickness, camber=camber, Sref=Sref, Lref=Lref, sweep=sweep)
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
            linked_af = self.__airfoils[i].get_scaled_copy(Sref=SLoc, Lref=LLoc)
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
        SLoc      = LLoc*(self.__eta[1,i+1]-self.__eta[1,i])
        SLoc_grad = LLoc_grad*(self.__eta[1,i+1]-self.__eta[1,i])+LLoc*(self.__eta_grad[1,i+1,:]-self.__eta_grad[1,i,:])
#         SLoc      = self.__span * LLoc / float(self.__n_sect)
#         SLoc_grad = (self.__span_grad * LLoc + self.__span * LLoc_grad) / float(self.__n_sect)
        return LLoc, LLoc_grad, SLoc, SLoc_grad
    
    def __compute_Sref_Lref_ToC_AR_fuel(self):
        """
        Compute Lref and Sref from airfoils information, ToC and AR
        """
        N = self.__n_sect
        
        self.__Sref = 0.
        self.__Lref = 0.
        self.__fuel = 0.
        self.__Sref_grad = zeros(self.__ndv)
        self.__Lref_grad = zeros(self.__ndv)
        self.__fuel_grad = zeros(self.__ndv)

        for i,af in enumerate(self.__linked_airfoils):
            self.__Sref+=af.get_Sref()
            self.__Lref+=af.get_Lref()
            self.__fuel+=self.__rel_thicks[i]*self.__chords[i]*af.get_Sref()*0.5
            self.__Sref_grad+=af.get_Sref_grad()
            self.__Lref_grad+=af.get_Lref_grad()
            self.__fuel_grad+=self.__rel_thicks_grad[i]*self.__chords[i]*af.get_Sref()*0.5\
                            + self.__rel_thicks[i]*self.__chords_grad[i]*af.get_Sref()*0.5\
                            + self.__rel_thicks[i]*self.__chords[i]*af.get_Sref_grad()*0.5\
            
        self.__Lref/=N
        self.__Lref_grad/=N
        
        self.__root_ToC = self.__root_height/self.__root_chord
        self.__root_ToC_grad = (self.__root_height_grad*self.__root_chord - self.__root_height*self.__root_chord_grad)/self.__root_chord**2
        
        if self.__break_chord is not None:
            self.__break_ToC = self.__break_height/self.__break_chord
            self.__break_ToC_grad = (self.__break_height_grad*self.__break_chord-self.__break_height*self.__break_chord_grad)/self.__break_chord**2
            
        if self.__tip_chord is not None:
            self.__tip_ToC = self.__tip_height/self.__tip_chord
            self.__tip_ToC_grad = (self.__tip_height_grad*self.__tip_chord-self.__tip_height*self.__tip_chord_grad)/self.__tip_chord**2
            
        self.__AR = self.__span**2/self.__Sref
        self.__AR_grad = (2.*self.__span*self.__span_grad*self.__Sref - self.__span**2*self.__Sref_grad) / self.__Sref**2
    
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
            
#     def build_view(self):
#         from padgeGUI.ModelViewer.QTPadge import QTPadge
#         self.update()
#         LE_points=[]
#         TE_points=[]
#         n_eta=len(self.__r_list_eta)
#         for i in xrange(n_eta):
#             LE_points.append(self.__CE_manager.create_point(self.__tag+'.LE'+str(i)))
#             TE_points.append(self.__CE_manager.create_point(self.__tag+'.TE'+str(i)))
#         
#         for i in xrange(n_eta):
#             xyz_LE=[self.__eta[0,i],self.__eta[1,i],0.]
#             LE_points[i].set_xyz(xyz_LE)
#             xyz_grad_LE=zeros((3,self.__ndv))
#             xyz_grad_LE[1,:]=self.__eta_grad[0,i,:]
#             xyz_grad_LE[2,:]=self.__eta_grad[1,i,:]
#             LE_points[i].set_xyz_grad(xyz_grad_LE)
#             xyz_TE=[self.__eta[0,i]+self.__chords_eta[i],self.__eta[1,i],0.]
#             TE_points[i].set_xyz(xyz_TE)
#             xyz_grad_TE=zeros((3,self.__ndv))
#             xyz_grad_TE[1,:]=self.__eta_grad[0,i,:]+self.__chords_grad_eta[i,:]
#             xyz_grad_TE[2,:]=self.__eta_grad[1,i,:]
#             TE_points[i].set_xyz_grad(xyz_grad_TE)
#             
#         chords_curves=[]
#         for i in xrange(n_eta):
#             chords_curves.append(self.__CE_manager.create_line(self.__tag+'.chord'+str(i)))
#             chords_curves[i].set_controller('P1',LE_points[i].get_id())
#             chords_curves[i].set_controller('P2',TE_points[i].get_id())
#             
#         LE_curve1 = self.__CE_manager.create_curve_interp(self.__tag+'.LEcurve1', smoothing=True)
#         LE_curve1.set_tangent_method('bessel')
#         for i in xrange(0,(n_eta+1)/2):
#             pt=LE_points[i]
#             LE_curve1 .set_controller('P'+str(i+1),pt.get_id())
#             
#         LE_curve2 = self.__CE_manager.create_curve_interp(self.__tag+'.LEcurve2', smoothing=True)
#         LE_curve2.set_tangent_method('bessel')
#         for i in xrange((n_eta-1)/2,n_eta):
#             pt=LE_points[i]
#             LE_curve2.set_controller('P'+str(i+1-(n_eta-1)/2),pt.get_id())
#             
#         TE_curve1 = self.__CE_manager.create_curve_interp(self.__tag+'.TEcurve1', smoothing=True)
#         TE_curve1.set_tangent_method('bessel')
#         for i in xrange(0,(n_eta+1)/2):
#             pt=TE_points[i]
#             TE_curve1.set_controller('P'+str(i+1),pt.get_id())
#             
#         TE_curve2 = self.__CE_manager.create_curve_interp(self.__tag+'.TEcurve2', smoothing=True)
#         TE_curve2.set_tangent_method('bessel')
#         for i in xrange((n_eta-1)/2,n_eta):
#             pt=TE_points[i]
#             TE_curve2.set_controller('P'+str(i+1-(n_eta-1)/2),pt.get_id())
#             
#         self.__PCADModel.update()
#         
#         
#         GUI=QTPadge(self.__PCADModel, view_mode='2D')
#         GUI.set_plandef('xy')
#         GUI.launch()
