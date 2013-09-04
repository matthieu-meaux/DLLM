# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard
# @author: Matthieu MEAUX (for refactoring)
 
# - Local imports -
from DLLM.polarManager.differentiatedAeroShape import DifferentiatedAeroShape
from DLLM.DLLMKernel.DLLMGeom import DLLMGeom
from DLLM.DLLMKernel.DLLMDirect import DLLMDirect

class LiftingLineWing(DifferentiatedAeroShape):
    def __init__(self, wing_geom, airfoils, OC):
        '''
        Constructor for wings based on lifting line theory
        @param wing_geom : the wing geometry
        @param Lref : Reference length for the moments computation
        @param relaxFactor : relaxation factor for induced angles computation
        @param stopCriteria : the stop criteria for the iterative method for computing the wing circulation.
        '''
        DifferentiatedAeroShape.__init__(self,1.,1.) # Sref and Lref are computed by __set_Lref_Sref
        
        self.__wing_geom  = wing_geom
        self.__airfoils   = airfoils
        self.__OC         = OC
        
        self.__DLLMGeom   = DLLMGeom(self)
        self.__DLLMDirect = DLLMDirect(self)
    
    #-- Accessors
    def get_wing_geom(self):
        return self.__wing_geom
    
    def get_airfoils(self):
        return self.__airfoils
    
    def get_OC(self):
        return self.__OC
    
    def get_K(self):
        return self.__DLLMGeom.get_K()
    
    #-- Setters
    def set_relax_factor(self, relax_factor):
        self.__DLLMDirect.set_relax_factor(relax_factor)
        
    def set_stop_criteria(self, residual=None, n_it=None):
        self.__DLLMDirect.set_stop_criteria(residual=residual, n_it=n_it)
        
    #-- Run methods
    def run_direct(self):
        self.__DLLMDirect.run()
        
    #-- direct solver methods
    def __init_iterations(self):
        self.__iAoA=zeros([self.__N])
        R,DR_DiAoA=self.R(self.__iAoA,alpha,Mach)
        self.__R0       = norm(R)
        self.__residual = 1.0
        self.__residuals_hist = []
        
    def __sub_iteration(self):
        R,DR_DiAoA=self.R(self.__iAoA,alpha,Mach)
        #Newton zero search method
        self.__iAoA-=self.__relaxFactor*solve(DR_DiAoA,R)
        
        #Compute stop criteria
        self.__residual=norm(R)/R0
        self.__residuals_hist.append(self.__residual)
        print "||R||/||R0||= "+str(self.__residual)
        
    def __Iterate(self,alpha,Mach):
        '''
        For a given angle of attack, iterates to compute the circulation and the downwash angles.
        @param alpha : angle of attack
        '''
        self.__init_iterations()
        
        if self.__stop_criteria_type == self.STOP_CRITERIA_RESIDUAL:
            while(self.__residual>self.__stopCriteria):
                self.__sub_iteration()
        else:
            for i in xrange(self.__stopCriteria):
                self.__sub_iteration() 
         
    def get_convergence_history(self):
        """
        Accessor to the last computation convergence history as a list of residuals normalized by the first iteration residual.
        """
        return self.__residuals_hist
        
    def compute_localAoA(self,iaOa,alpha):
        '''
        Computes the local angle of attack = geometric AoA - induced downwash angle - twist
        @param alpha : the wing angle of attack
        '''
        self.__localAoA=alpha-iaOa-self.get_wing_geom().get_twist()
        
        for i in range(self.__N):
            if self.__localAoA[i] > math.pi/2. or self.__localAoA[i] < -math.pi/2.:
                raise Exception, "Local angle of attack out of bounds [-pi/2, pi/2]"

        return self.__localAoA,self.__DlocalAoA_DIaOa#last term is constant
        
    def dLocalAOa_DTwist(self):
        '''
        Computes the local angle of attack sensibility to twist
        '''
        return self.__DlocalAoA_DTwist #constant due to the expression of localAoA

    def dLocalAOa_DAoA(self):
        '''
        Computes the local angle of attack sensibility to twist
        '''
        return self.__DlocalAoA_DAoA #constant due to the expression of localAoA
            
    def compute_gamma(self,localAOa,Mach):
        '''
        Updates the circulation
        '''
        for i in range(self.__N):
            self.__gamma[i]=self.__airfoils[i].gamma(localAOa[i],Mach=Mach)
            self.__Dgamma_DlocalAoA[i,i]=self.__airfoils[i].dGammaDAoA(localAOa[i],Mach=Mach)
            self.__Dgamma_Dthickness[i,i]=self.__airfoils[i].dGammaDThickness(localAOa[i],Mach=Mach)
            
        return self.__gamma, self.__Dgamma_DlocalAoA,self.__Dgamma_Dthickness
    
    def compute_dGamma(self,gamma):
        '''
        Computes the circulation y derivation
        '''
        self.__dGamma[0:self.__N]=gamma
        self.__dGamma[self.__N]=0.0
        self.__dGamma[1:self.__N+1]-=gamma
        
        #Differenciation
        self.__DdGammaDy_DGamma[0:self.__N,:]=diag(ones([self.__N]))
        self.__DdGammaDy_DGamma[self.__N,:]=0.0
        self.__DdGammaDy_DGamma[1:self.__N+1,:]-=diag(ones([self.__N]))
        
        return self.__dGamma, self.__DdGammaDy_DGamma
    
    def compute_iAoA(self,dGamma):
        '''
        Computes the induced angle on an airfoil for a given circulation on the wing.
        '''
        self.__iAoANew=numpy.dot(self.__K,dGamma)
        self.__DiAoA_DdGamma=numpy.dot(self.__K,diag(ones(self.__N+1)))
        
        return self.__iAoANew, self.__DiAoA_DdGamma

    def R(self,iAoA,alpha,Mach):
        '''
        Computes the residual and its derivatives
        @param alpha : the wing angle of attack
        @param iAoA : the induced angles vector
        @return: Residual, Gradient of residual
        '''
        iAoAold=copy(iAoA)
        localAoA,DlocalAoA_DIaOa=self.compute_localAOa(iAoA,alpha)
        
        gamma, Dgamma_DlocalAoA,Dgamma_Dthickness=self.compute_gamma(localAoA,Mach)
        Dgamma_DiAoA=dot(Dgamma_DlocalAoA,DlocalAoA_DIaOa)
        
        dGamma, DdGammaDy_DGamma=self.compute_dGamma(gamma)
        DdGammaDy_DiAoA=dot(DdGammaDy_DGamma,Dgamma_DiAoA)
        
        iAoANew, DiAoA_DdGamma=self.compute_iAoA(dGamma)
        
        DiAoA_DiAoA=dot(DiAoA_DdGamma,DdGammaDy_DiAoA)
        
        self.__R=iAoAold-iAoANew
        self.__DR_DiAoA=numpy.diag(ones([self.__N]))-DiAoA_DiAoA
        
        return self.__R,self.__DR_DiAoA
    
    def dRDTwist(self,iAoA,alpha,Mach):
        '''
        Sensibility of residuals to twist.
        @param iAoA : the induced angles vector
        @param alpha : the wing angle of attack
        '''

        localAoA,DlocalAoA_DIaOa=self.compute_localAOa(iAoA,alpha)
        
        gamma, Dgamma_DlocalAoA,Dgamma_Dthickness=self.compute_gamma(localAoA,Mach)
        Dgamma_DTwist=dot(Dgamma_DlocalAoA,self.dLocalAOa_DTwist())
        
        dGamma, DdGammaDy_DGamma=self.compute_dGamma(gamma)
        DdGammaDy_DTwist=dot(DdGammaDy_DGamma,Dgamma_DTwist)
        
        iAoANew, DiAoA_DdGamma=self.compute_iAoA(dGamma)
        
        dRdTwist=-dot(DiAoA_DdGamma,DdGammaDy_DTwist)
        
        return dRdTwist

    def dRDAoA(self,iAoA,alpha,Mach):
        '''
        Sensibility of residuals to angle of attack
        @param iAoA : the induced angles vector
        @param alpha : the wing angle of attack
        '''

        localAoA,DlocalAoA_DIaOa=self.compute_localAOa(iAoA,alpha)
        
        gamma, Dgamma_DlocalAoA,Dgamma_Dthickness=self.compute_gamma(localAoA,Mach)
        Dgamma_DAoA=dot(Dgamma_DlocalAoA,self.dLocalAOa_DAoA())
        
        dGamma, DdGammaDy_DGamma=self.compute_dGamma(gamma)
        DdGammaDy_DAoA=dot(DdGammaDy_DGamma,Dgamma_DAoA)
        
        iAoANew, DiAoA_DdGamma=self.compute_iAoA(dGamma)
        
        dRAoA=-dot(DiAoA_DdGamma,DdGammaDy_DAoA)# (-) Because R=Gamma_input-Gamma_computed_from_input
        
        return dRAoA
    
    def dRDThickness(self,iAoA,alpha,Mach):
        '''
        Sensibility of residuals to relative thickness of airfoils.
        @param iAoA : the induced angles vector
        @param alpha : the wing angle of attack
        '''

        localAoA,DlocalAoA_DIaOa=self.compute_localAOa(iAoA,alpha)
        
        gamma, Dgamma_DlocalAoA,Dgamma_Dthickness=self.compute_gamma(localAoA,Mach)
        
        dGamma, DdGammaDy_DGamma=self.compute_dGamma(gamma)
        DdGammaDy_DThickness=dot(DdGammaDy_DGamma,Dgamma_Dthickness)
        
        iAoANew, DiAoA_DdGamma=self.compute_iAoA(dGamma)
        
        dRdThickness=-dot(DiAoA_DdGamma,DdGammaDy_DThickness)
        
        return dRdThickness
    
    def adjoint_correction(self,df_diAoA):
        '''
        Computes the adjoint correction
        @param df_diAoA : sensibility of the objective function to the state vector : induced aoa
        '''
        adj=self.adjoint(df_diAoA)
        
        return dot(adj.T,self.__R)
    
    def adj_correction_Cl(self,alpha, beta=0., Mach=0.):
        '''
        Computes the adjoint correction to Cl
        '''
        df_diAoA=self.dCl_dIaOa(alpha, beta, Mach)
        
        return self.adjoint_correction(df_diAoA)

    def adj_correction_Cd(self,alpha, beta=0., Mach=0.):
        '''
        Computes the adjoint correction to Cd
        '''
        df_diAoA=self.dCd_dIaOa(alpha, beta, Mach)
        
        return self.adjoint_correction(df_diAoA)
    
    def adj_correction_Cm(self,alpha, beta=0., Mach=0.):
        '''
        Computes the adjoint correction to Cm
        '''
        df_diAoA=self.dCm_dIaOa(alpha, beta, Mach)
        
        return self.adjoint_correction(df_diAoA)        
    
    def adjoint(self,dJ_diAoA):
        '''
        Computes the adjoint of the problem
        @param dJ_diAoA : sensibility of the objective function to the state vector : induced aoa
        '''

        return solve(self.__DR_DiAoA.T,-dJ_diAoA.T)
    
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
            
    def Cl(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        Cl=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            Cl+=af.Cl(self.__localAoA[i],Mach)*af.getSref()
        Cl/=self.getSref()
        
        return Cl
    
    def Cm(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float        
        '''
        self.__Iterate(alpha,Mach)
 
        Cm=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            Cm+=af.Cm(self.__localAoA[i],Mach)*af.getSref()
        Cm/=self.getSref()
        
        return Cm
 
    def Cd(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag function
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        #print "lifting line thickness = "+str(self._wing_geom.get_relative_thickness())
        Cd=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            aOa=self.__localAoA[i]
            Cdloc=af.Cl(aOa,Mach)*self.__iAoA[i]+af.Cd(aOa,Mach)
            Cd+=Cdloc*af.getSref()
 
        Cd/=self.getSref()
 
        return Cd
    
    def DCl_DTwist(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function Total derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCl_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dT=self.__d_part_Cl_d_part_Twist(alpha,Mach=Mach)
        
        return self.DJ_DTwist(dJ_dT,adj,alpha,Mach=Mach)
    
    def DCl_DAoA(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function Total derivation with respect to angle of attack
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCl_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dAoA=self.__d_part_Cl_d_part_AoA(alpha,Mach=Mach)
        
        return self.DJ_DAoA(dJ_dAoA,adj,alpha,Mach=Mach)
    
    def DCd_DTwist(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag coefficient function Total derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCd_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dT=self.__d_part_Cd_d_part_Twist(alpha,Mach=Mach)
        
        return self.DJ_DTwist(dJ_dT,adj,alpha,Mach=Mach)
    
    def DCd_DAoA(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag coefficient function Total derivation with respect to angle of attack
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCd_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dAoA=self.__d_part_Cd_d_part_AoA(alpha,Mach=Mach)
        
        return self.DJ_DAoA(dJ_dAoA,adj,alpha,Mach=Mach)
    
    def DCm_DTwist(self,alpha,beta=0.0,Mach=0.0):
        '''
        Pitch moment coefficient function Total derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCm_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dT=self.__d_part_Cm_d_part_Twist(alpha,Mach=Mach)
        
        return self.DJ_DTwist(dJ_dT,adj,alpha,Mach=Mach)
    
    def DCm_DAoA(self,alpha,beta=0.0,Mach=0.0):
        '''
        Moment coefficient function Total derivation with respect to angle of attack
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCm_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dAoA=self.__d_part_Cm_d_part_AoA(alpha,Mach=Mach)
        
        return self.DJ_DAoA(dJ_dAoA,adj,alpha,Mach=Mach)
    
    def DCl_DThickness(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function Total derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCl_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dT=self.__d_part_Cl_d_part_Thickness(alpha,Mach=Mach)
        
        return self.DJ_DThick(dJ_dT,adj,alpha,Mach=Mach)
    
    def DCd_DThickness(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag coefficient function Total derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCd_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dT=self.__d_part_Cd_d_part_Thickness(alpha,Mach=Mach)
        
        return self.DJ_DThick(dJ_dT,adj,alpha,Mach=Mach)
    
    def DCm_DThickness(self,alpha,beta=0.0,Mach=0.0):
        '''
        Pitch moment coefficient function Total derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dJ_diAoA=self.dCm_dIaOa(alpha,Mach=Mach)
        
        adj=self.adjoint(dJ_diAoA)
        
        dJ_dT=self.__d_part_Cm_d_part_Thickness(alpha,Mach=Mach)
        
        return self.DJ_DThick(dJ_dT,adj,alpha,Mach=Mach)
    
    def dCm_dIaOa(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function partial derivation with respect to iduced angles
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
 
        dCm=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            dCm+=af.CmAlpha(self.__localAoA[i],Mach)*af.getSref()*self.__DlocalAoA_DIaOa[i]
        dCm/=self.getSref()
        
        return dCm
     
    def dCl_dIaOa(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function partial derivation with respect to iduced angles
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        
        dCl=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            dCl+=af.ClAlpha(self.__localAoA[i],Mach=Mach)*af.getSref()*self.__DlocalAoA_DIaOa[i]
        dCl/=self.getSref()
        
        return dCl
    
    def dCd_dIaOa(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag function partial derivation with respect to iduced angles
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
 
        dCd_dAoA=zeros(self.__N)

        for i in range(self.__N):#Dependance by induced angles
            af=self.__airfoils[i]
            aOa=self.__localAoA[i]
            dCd_dAoA[i]=(af.ClAlpha(aOa,Mach)*self.__iAoA[i]+af.CdAlpha(aOa,Mach))*af.getSref()
 
        dCd=dot(dCd_dAoA,self.__DlocalAoA_DIaOa)
        
        for i in range(self.__N):#Dependance by projection of Cl : the induced drag.
            af=self.__airfoils[i]
            dCd[i]+=af.Cl(self.__localAoA[i],Mach)*af.getSref()
        
        dCd/=self.getSref()
        
        return dCd
    
    def __d_part_Cl_d_part_Twist(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function partial derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dlAoADT=self.dLocalAOa_DTwist()
        
        dCl=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            dCl+=af.ClAlpha(self.__localAoA[i],Mach)*af.getSref()*dlAoADT[i]
        dCl/=self.getSref()
        
        return dCl
    
    def __d_part_Cl_d_part_AoA(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function partial derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dlAoAdAoA=self.dLocalAOa_DAoA()
        
        dCl=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            dCl+=af.ClAlpha(self.__localAoA[i],Mach)*af.getSref()*dlAoAdAoA[i]
        dCl/=self.getSref()
        
        return dCl
    
    def __d_part_Cl_d_part_Thickness(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function partial derivation with respect to thickness
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dCl=zeros(self.__N)
        for i in range(self.__N):
            af=self.__airfoils[i]
            dCl[i]=af.dCl_dthickness(self.__localAoA[i],Mach)*af.getSref()
        dCl/=self.getSref()
        
        return dCl
    
    def __d_part_Cm_d_part_Twist(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function partial derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dlAoADT=self.dLocalAOa_DTwist()
        
        dCm=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            dCm+=af.CmAlpha(self.__localAoA[i],Mach)*af.getSref()*dlAoADT[i]
        dCm/=self.getSref()
        
        return dCm
    
    def __d_part_Cm_d_part_AoA(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function partial derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dlAoAdAoA=self.dLocalAOa_DAoA()
        
        dCm=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            dCm+=af.CmAlpha(self.__localAoA[i],Mach)*af.getSref()*dlAoAdAoA[i]
        dCm/=self.getSref()
        
        return dCm
    def __d_part_Cm_d_part_Thickness(self,alpha,beta=0.0,Mach=0.0):
        '''
        Lift coefficient function partial derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dCm=zeros(self.__N)
        for i in range(self.__N):
            af=self.__airfoils[i]
            dCm[i]=af.dCm_dthickness(self.__localAoA[i],Mach)*af.getSref()
        dCm/=self.getSref()
        
        return dCm
 
    def __d_part_Cd_d_part_Twist(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag function partial derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dlAoADT=self.dLocalAOa_DTwist()
 
        dCd=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            aOa=self.__localAoA[i]
            Cdloc=(af.ClAlpha(aOa,Mach)*self.__iAoA[i]+af.CdAlpha(aOa,Mach))*dlAoADT[i]
            dCd+=Cdloc*af.getSref()
 
        dCd/=self.getSref()
 
        return dCd
    
    def __d_part_Cd_d_part_AoA(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag function partial derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dlAoAdAoA=self.dLocalAOa_DAoA()
 
        dCd=0.0
        for i in range(self.__N):
            af=self.__airfoils[i]
            aOa=self.__localAoA[i]
            Cdloc=(af.ClAlpha(aOa,Mach)*self.__iAoA[i]+af.CdAlpha(aOa,Mach))*dlAoAdAoA[i]
            dCd+=Cdloc*af.getSref()
 
        dCd/=self.getSref()
 
        return dCd

    def __d_part_Cd_d_part_Thickness(self,alpha,beta=0.0,Mach=0.0):
        '''
        Drag function partial derivation with respect to twist law
        @param alpha: angle of Attack 
        @param type alpha : Float
        @param Mach : Mach number
        @param type Mach : Float
        @param beta : sideslip angle
        @param type beta : Float
        '''
        self.__Iterate(alpha,Mach)
        dCd=zeros(self.__N)
        for i in range(self.__N):
            af=self.__airfoils[i]
            aOa=self.__localAoA[i]
            dCd[i]=af.dCl_dthickness(aOa,Mach)*self.__iAoA[i]+af.dCd_dthickness(aOa,Mach)
            dCd[i]*=af.getSref()
        dCd/=self.getSref()
        
        return dCd  
       
    def ClAlpha(self,alpha, beta=0., Mach=0.):
        '''
        Sensibility of Cl to alpha.
        Done by adjoint method.
        @param alpha: angle of Attack 
        @param type alpha : Float
        '''
        return self.DCl_DAoA(alpha, beta, Mach)
    
    def CdAlpha(self,alpha, beta=0., Mach=0.):
        '''
        Sensibility of Cd to alpha.
        Done by adjoint method.
        @param alpha: angle of Attack 
        @param type alpha : Float
        '''
        return self.DCd_DAoA(alpha, beta, Mach)
    
    def CmAlpha(self,alpha, beta=0., Mach=0.):
        '''
        Sensibility of Cm to alpha.
        Done by adjoint method.
        @param alpha: angle of Attack 
        @param type alpha : Float
        '''
        return self.DCm_DAoA(alpha, beta, Mach)
    
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
        
