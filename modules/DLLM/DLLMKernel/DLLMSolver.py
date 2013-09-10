# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard
# @author: Matthieu MEAUX (for refactoring)
 
# - Local imports -
from DLLM.DLLMKernel.DLLMMesh import DLLMMesh
from DLLM.DLLMKernel.DLLMDirect import DLLMDirect
from DLLM.DLLMKernel.DLLMPost import DLLMPost
from DLLM.DLLMKernel.DLLMAdjoint import DLLMAdjoint

class DLLMSolver:
    ERROR_MSG='ERROR in DLLMSolver.'
    def __init__(self, wing_geom, airfoils, OC):
        '''
        Constructor for wings based on lifting line theory
        @param wing_geom : the wing geometry
        @param Lref : Reference length for the moments computation
        @param relaxFactor : relaxation factor for induced angles computation
        @param stopCriteria : the stop criteria for the iterative method for computing the wing circulation.
        '''
        
        self.__wing_geom  = wing_geom
        self.__airfoils   = airfoils
        self.__OC         = OC
        
        self.__Lref       = 0.
        self.__Sref       = 0.
        
        self.__DLLMMesh   = DLLMMesh(self)
        self.__DLLMDirect = DLLMDirect(self)
        self.__DLLMPost   = DLLMPost(self)
        self.__DLLMAdjoint= DLLMAdjoint(self)
        
        self.set_OC(OC)
    
    #-- Accessors
    def get_wing_geom(self):
        return self.__wing_geom
    
    def get_airfoils(self):
        return self.__airfoils
    
    def get_OC(self):
        return self.__OC
    
    def get_Lref(self):
        return self.__Lref
    
    def get_Sref(self):
        return self.__Sref
    
    def get_K(self):
        return self.__DLLMMesh.get_K()
    
    def get_convergence_history(self):
        return self.__DLLMDirect.get_convergence_history()
    
    def get_localAoA(self):
        return self.__DLLMDirect.get_localAoA()
    
    def get_iAoA(self):
        return self.__DLLMDirect.get_iAoA()
    
    def get_R(self):
        return self.__DLLMDirect.get_R()
    
    def get_DR_DiAoA(self):
        return self.__DLLMDirect.get_DR_DiAoA()
    
    def get_DlocalAoA_DiAoA(self):
        return self.__DLLMDirect.get_DlocalAoA_DiAoA()
    
    def get_func_list(self):
        return self.__DLLMPost.get_func_list()
    
    def get_dfunc_diAoA(self):
        return self.__DLLMPost.get_dfunc_diAoA()
    
    def is_direct_computed(self):
        return self.__DLLMDirect.is_computed()
    
    def is_post_computed(self):
        return self.__DLLMPost.is_computed()
        
    #-- Setters
    def __reinit_modules(self):
        self.__DLLMDirect.set_computed(False)
        self.__DLLMPost.set_computed(False)
        self.__DLLMAdjoint.set_computed(False)
        
    def set_OC(self, OC):
        self.__OC = OC
        self.__reinit_modules()
        
    def set_wing_geom(self, wing_geom):
        self.__wing_geom  = wing_geom
        self.__DLLMMesh.recompute()
        self.__reinit_modules()
        
    def set_airfoil(self, airfoils):
        self.__airfoils = airfoils
        self.__DLLMMesh.recompute()
        self.__reinit_modules()
        
    def set_Lref(self, Lref):
        self.__Lref = Lref
        
    def set_Sref(self, Sref):
        self.__Sref = Sref
    
    def set_relax_factor(self, relax_factor):
        self.__DLLMDirect.set_relax_factor(relax_factor)
        
    def set_stop_criteria(self, residual=None, n_it=None):
        self.__DLLMDirect.set_stop_criteria(residual=residual, n_it=n_it)
        
    #-- Run methods
    def run_direct(self):
        self.__DLLMDirect.run()
        
    def run_post(self,func_list=None):
        ERROR_MSG=self.ERROR_MSG+'run_post: '
        if self.is_direct_computed():
            self.__DLLMPost.run(func_list=func_list)
        else:
            print ERROR_MSG+'Cannot run post-processing if solution is not computed'
            
    def run_adjoint(self):
        ERROR_MSG=self.ERROR_MSG+'run_adjoint: '
        if self.is_post_computed():
            self.__DLLMAdjoint.run()
        else:
            print ERROR_MSG+'Cannot run adjoint if post-processing is not computed'       
    
    def DJ_DTwist(self,dJ_dTwist,adjoint,alpha,Mach):
        '''
        Builds the gradient of J with the adjoint field
        @param dJ_dTwist : partial derivative of the objective function to the twist
        @param adjoint : the adjoint field
        @param alpha : the wing angle of attack
        '''
        self.__Iterate(alpha,Mach)
        dR=self.dRDTwist(self.__iAoA,alpha,Mach)
        
        return dJ_dTwist+dot(adjoint.T,dR)
    
    def DJ_DAoA(self,dJ_dAoA,adjoint,alpha,Mach):
        '''
        Builds the gradient of J with the adjoint field
        @param dJ_dTwist : partial derivative of the objective function to the twist
        @param adjoint : the adjoint field
        @param alpha : the wing angle of attack
        '''
        self.__Iterate(alpha,Mach)
        dR=self.dRDAoA(self.__iAoA,alpha,Mach)
        
        return dJ_dAoA+dot(adjoint.T,dR)

    def DJ_DThick(self,dJ_dThickness,adjoint,alpha,Mach):
        '''
        Builds the gradient of J with the adjoint field
        @param dJ_dThickness : partial derivative of the objective function to the twist
        @param adjoint : the adjoint field
        @param alpha : the wing angle of attack
        '''
        self.__Iterate(alpha,Mach)
        dR=self.dRDThickness(self.__iAoA,alpha,Mach)
        
        return dJ_dThickness+dot(adjoint.T,dR)
            
    def set_twist(self,twistLaw):
        '''
        Sets the twist law of the wing
        @param twistLaw : the twist law
        '''
        if type(twistLaw)==type([]):
            twist=array(twistLaw)
        elif type(twistLaw)==type(array([0.])):
            twist=twistLaw
        else:
            raise Exception, "Incorrect type for twistLaw : "+str(type(twistLaw))
        
        self.get_wing_geom().set_twist(twist)

    def set_relative_thickness(self,thickness):
        """
        Setter for the height of the airfoils
        """
        if type(thickness)==type([]):
            thick=array(thickness)
        elif type(thickness)==type(array([0.])):
            thick=thickness
        else:
            raise Exception, "Incorrect type for thickness : "+str(type(thickness))
        
        self.get_wing_geom().set_relative_thickness(thick)
        
        print "LLW set_relative_thickness thick = "+str(thick)
        for airfoil, thickness in zip(self.__airfoils,thick):
            airfoil.set_relative_thickness(thickness)
    
    def write_gamma_to_file(self,file_name,alpha,beta=0.0,Mach=0.0):
        '''
        Writes the circulation repartition in a file
        @param file_name : the file to write data
        @param type file_name : String
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        fid=open(file_name,'w')
        line="#Slice\t%24s"%"Circulation"+"\n"
        fid.write(line)
        i=0
        for i in range(len(self.__gamma)):
            line=str(i)+"\t%24.16e"%self.__gamma[i]+"\n"
            fid.write(line)
        fid.close()
