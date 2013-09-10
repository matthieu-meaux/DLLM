# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: Francois Gallard
# @author: Matthieu MEAUX (for refactoring)

from numpy import zeros, dot

class DLLMPost:
    def __init__(self, LLW):
        """
        Post-Processing module compatible with adjoint for DLLM
        """
        self.__LLW = LLW
        self.__func_list   = None
        self.__func_dim    = 0
        self.__func_values = None
        self.__dfunc_diAoA = None
        self.__N           = self.get_wing_geom().get_n_elements()
        self.__computed    = False
        
    #-- computed related methods
    def is_computed(self):
        return self.__computed
    
    def set_computed(self, bool=True):
        self.__computed = bool
        
    #-- Accessors
    def get_func_list(self):
        return self.__func_list
    
    def get_dfunc_diAoA(self):
        return self.__dfunc_diAoA
    
    def get_wing_geom(self):
        return self.__LLW.get_wing_geom()
    
    def get_airfoils(self):
        return self.__LLW.get_airfoils()
    
    def get_OC(self):
        return self.__LLW.get_OC()
    
    def get_localAoA(self):
        return self.__LLW.get_localAoA()
    
    def get_DlocalAoA_DiAoA(self):
        return self.__LLW.get_DlocalAoA_DiAoA()
    
    def get_iAoA(self):
        return self.__LLW.get_iAoA()
    
    def get_Sref(self):
        return self.__LLW.get_Sref()
    
    #-- Setters 
    def set_func_list(self, func_list):
        if func_list is None:
            #func_list = ['Lift', 'Drag', 'Moment', 'Cl', 'Cd', 'Cm']
            func_list = ['Lift', 'Drag', 'Cl', 'Cd']
        self.__func_list   = func_list
        self.__func_dim    = len(func_list)
        self.__func_values = [None]*self.__func_dim
        self.__dfunc_diAoA = [None]*self.__func_dim
        
    #-- Run method
    def run(self, func_list=None):
        self.set_func_list(func_list)
        for i,func in enumerate(self.__func_list):
            if   func == 'Cl':
                val     = self.__Cl()
                dfdiAoA = self.__dCl_diAoA()
            elif func == 'Cd':
                val = self.__Cd()
                dfdiAoA = self.__dCd_diAoA()
            elif func == 'Cm':
                val = self.__Cm()
            elif func == 'Lift':
                val = self.__Lift()
                dfdiAoA = self.__dLift_diAoA()
            elif func == 'Drag':
                val= self.__Drag()
                dfdiAoA = self.__dDrag_diAoA()
            else:
                val = None
            self.__func_values[i]=val
            self.__dfunc_diAoA[i]=dfdiAoA
        
        self.__display_info()
        self.set_computed(True)
                
    #-- Computation methods
    def __Cl(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        Cl=0.0
        for i in xrange(self.__N):
            af = airfoils[i]
            Cl+=af.Cl(localAoA[i],Mach)*af.get_Sref()
        Cl/=self.get_Sref()
        return Cl
    
    def __dCl_diAoA(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        DlocalAoA_DiAoA = self.get_DlocalAoA_DiAoA()
        airfoils = self.get_airfoils()
        dCl_diAoA=zeros(self.__N)
        for i in range(self.__N):
            af = airfoils[i]
            dCl_diAoA+=af.ClAlpha(localAoA[i],Mach=Mach)*af.get_Sref()*DlocalAoA_DiAoA[i]
        dCl_diAoA/=self.get_Sref()
        
        return dCl_diAoA
    
    def __Cd(self):
        Mach     = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        iAoA     = self.get_iAoA()
        Cd=0.0 
        for i in xrange(self.__N):
            af      = airfoils[i]
            AoA     = localAoA[i]
            iAoAloc = iAoA[i]
            Cdloc=af.Cl(AoA,Mach)*iAoAloc+af.Cd(AoA,Mach)
            Cd+=Cdloc*af.get_Sref()
 
        Cd/=self.get_Sref()
        return Cd
    
    def __dCd_diAoA(self):
        Mach     = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        iAoA     = self.get_iAoA()
        DlocalAoA_DiAoA = self.get_DlocalAoA_DiAoA()
        
        dCd_dAoA=zeros(self.__N)
        for i in range(self.__N):#Dependance by induced angles
            af=airfoils[i]
            AoA=localAoA[i]
            dCd_dAoA[i]=(af.ClAlpha(AoA,Mach)*iAoA[i]+af.CdAlpha(AoA,Mach))*af.get_Sref()
 
        dCd_diAoA=dot(dCd_dAoA,DlocalAoA_DiAoA)
        
        for i in range(self.__N):#Dependance by projection of Cl : the induced drag.
            af=airfoils[i]
            dCd_diAoA[i]+=af.Cl(localAoA[i],Mach)*af.get_Sref()
        
        dCd_diAoA/=self.get_Sref()
        
        return dCd_diAoA
    
    def __Cm(self):
        # Issue in moments calculation, no distances or sweep taken into account...
        print 'WARNING: Cm calculation to be checked...'
        Mach     = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        Cm=0.0
        for i in xrange(self.__N):
            af=airfoils[i]
            Cm+=af.Cm(localAoA[i],Mach)*af.get_Sref()
        Cm/=self.get_Sref()
        return Cm
    
    def __Lift(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        Cl = self.__Cl()
        Lift = Pdyn*Sref*Cl
        return Lift
    
    def __dLift_diAoA(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        dCl_diAoA=self.__dCl_diAoA()
        dLift_diAoA = Pdyn*Sref*dCl_diAoA
        return dLift_diAoA
        
    def __Drag(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        Cd = self.__Cd()
        Drag = Pdyn*Sref*Cd
        return Drag
    
    def __dDrag_diAoA(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        dCd_diAoA = self.__dCd_diAoA()
        dDrag_diAoA = Pdyn*Sref*dCd_diAoA
        return dDrag_diAoA
    
    #-- Display methods
    def __display_info(self):
        print self.get_OC()
        print '\n*** aerodynamic functions and coefficients ***'
        for i, func in enumerate(self.__func_list):
            if func in ['Lift', 'Drag']:
                unit = '(N)'
            else:
                unit = ''
            print self.__func_list[i]+'\t=\t'+str(self.__func_values[i])+' '+unit
    
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
       
#     def ClAlpha(self,alpha, beta=0., Mach=0.):
#         '''
#         Sensibility of Cl to alpha.
#         Done by adjoint method.
#         @param alpha: angle of Attack 
#         @param type alpha : Float
#         '''
#         return self.DCl_DAoA(alpha, beta, Mach)
#     
#     def CdAlpha(self,alpha, beta=0., Mach=0.):
#         '''
#         Sensibility of Cd to alpha.
#         Done by adjoint method.
#         @param alpha: angle of Attack 
#         @param type alpha : Float
#         '''
#         return self.DCd_DAoA(alpha, beta, Mach)
#     
#     def CmAlpha(self,alpha, beta=0., Mach=0.):
#         '''
#         Sensibility of Cm to alpha.
#         Done by adjoint method.
#         @param alpha: angle of Attack 
#         @param type alpha : Float
#         '''
#         return self.DCm_DAoA(alpha, beta, Mach)