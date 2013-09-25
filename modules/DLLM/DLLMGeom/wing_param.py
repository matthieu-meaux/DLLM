# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: EADS
# @version: 0.1
# @author: Matthieu Meaux

# - Local imports -
from padge.PCADEngine.Base.PCADModel import PCADModel
from numpy import zeros, array
from numpy import pi, sqrt

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
        
        
    def __init_discrete_attributes(self):
        # -- discrete attributes
        self.__span               = None
        self.__span_grad          = None
        
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
        
    def build_wing(self):
        self.__BC_manager.clean()
        self.__BC_manager.create_variable(self.__tag+'.span',0.)
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
        # -- span
        Id   = self.__tag+'.span'
        pt   = self.__BC_manager.get_pt(Id)
        val  = pt.get_value()
        grad = pt.get_gradient()
        self.__span      = val
        self.__span_grad = grad
        
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
            
        deg_to_rad = pi/180.
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
                    self.__chords[i]        = (self.__break_chord - self.__root_chord)*r/p + self.__root_chord
                    self.__chords_grad[i,:] = (self.__break_chord_grad[:] - self.__root_chord_grad[:])*r/p + self.__root_chord_grad[:] \
                                            + (self.__break_chord - self.__root_chord)*r*p_grad[:]/p**2
                    
                else:
                    self.__chords[i]        = (self.__tip_chord - self.__break_chord)*(r-p)/(1.-p) + self.__break_chord
                    self.__chords_grad[i,:] = (self.__tip_chord_grad[:] - self.__break_chord_grad[:])*(r-p)/(1.-p) + self.__break_chord_grad[:] \
                                            + (self.__tip_chord - self.__break_chord)*(1-r)*p_grad[:]/(1.-p)**2
                    
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
                    term = ((self.__break_height - self.__root_height)*r/p + self.__root_height)
                    self.__rel_thicks[i] = term / self.__chords[i]

                    term_grad = (self.__break_height_grad[:] - self.__root_height_grad[:])*r/p + self.__root_height_grad[:] \
                              + (self.__break_height - self.__root_height)*r*p_grad[:]/p**2
                    self.__rel_thicks_grad[i,:] = (term_grad[:]*self.__chords[i]+term*self.__chords_grad[i,:])/(self.__chords[i])**2
                
                else:
                    term = ((self.__tip_height - self.__break_height)*(r-p)/(1.-p) + self.__break_height)
                    self.__rel_thicks[i]  =  term / self.__chords[i]
                    
                    term_grad = (self.__tip_height_grad[:] - self.__break_height_grad[:])*(r-p)/(1.-p) + self.__break_height_grad[:] \
                               + (self.__tip_height- self.__break_height)*(1-r)*p_grad[:]/(1.-p)**2
                    self.__rel_thicks_grad[i,:] = (term_grad[:]*self.__chords[i]+term*self.__chords_grad[i,:])/(self.__chords[i])**2         
                               
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
            self.__XYZ[0,i]        = 0.25*self.__chords[i]
            self.__XYZ_grad[0,i,:] = 0.25*self.__chords_grad[i]
            
            r=float(i+0.5-n)/float(n)
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
           
#         
#     def build_Discrete_wing(self,N,twist_law=None):
#         '''
#         Builds an elliptic liftingLineWing
#         @param reference_airfoil : The reference airfoil from which others will be built
#         @param N : the number of elements for the whole wing
#         @param twist_law : the twist law
#         '''
#         if N%2!=0:
#             raise Exception, "The total number of elements in the wing must be odd."
#         
#         n=int(N/2)
#         XYZ=zeros([N,3])
#         eta=zeros([N+1,3])
#         chords=zeros(N)
#         heights=zeros(N)
#         wbs = zeros(N)  # wingbox length
#         wingspan=self.__wingspan
#         twist=self.__cast_twist(n,twist_law)
#         for i in range(2*n+1):
#             eta[i,1]=self.__compute_eta(i,n)
#             
#         for i in range(N):
#             y=self.__compute_y(i,n)
#             Lloc=self.__compute_chord(i,n)
#             x=Lloc/2.#3quarters chord
#             XYZ[i,:]=array([x,y,0])
#             chords[i]=Lloc
#             heights[i]=self.__compute_h(i,n)
#      
#         dw=Discrete_wing(XYZ,eta,chords,heights,twist,wingspan)
#         
#         self.__discrete_wing = dw
#         
#         return dw
# 
#     def __linear(self, x,y_min,y_max,x_min=0.,x_max=1.):
#         return (y_max-y_min)*float(x-x_min)/(float(x_max-x_min))+y_min
#     

