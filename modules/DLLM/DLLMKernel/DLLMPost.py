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
        self.__dpJ_dpTwist = None
        self.__dpJ_dpAoA   = None
        self.__dpJ_dpThick = None
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
    
    def get_dpJ_dpTwist(self):
        return self.__dpJ_dpTwist
    
    def get_dpJ_dpAoA(self):
        return self.__dpJ_dpAoA
    
    def get_dpJ_dpThickness(self):
        return self.__dpJ_dpThick
    
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
    
    def get_DlocalAoA_DTwist(self):
        return self.__LLW.get_DlocalAoA_DTwist()
    
    def get_DlocalAoA_DAoA(self):
        return self.__LLW.get_DlocalAoA_DAoA()
    
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
        self.__dpJ_dpTwist = [None]*self.__func_dim
        self.__dpJ_dpAoA   = [None]*self.__func_dim
        self.__dpJ_dpThick = [None]*self.__func_dim
        
    #-- Run method
    def run(self, func_list=None):
        self.set_func_list(func_list)
        for i,func in enumerate(self.__func_list):
            if   func == 'Cl':
                val     = self.__Cl()
                dfdiAoA = self.__dCl_diAoA()
                dpJdpTwist = self.__dpCl_dpTwist()
                dpJdpAoA   = self.__dpCl_dpAoA()
                dpJdpThick = self.__dpCl_dpThickness()
            elif func == 'Cd':
                val = self.__Cd()
                dfdiAoA = self.__dCd_diAoA()
                dpJdpTwist = self.__dpCd_dpTwist()
                dpJdpAoA   = self.__dpCd_dpAoA()
                dpJdpThick = self.__dpCd_dpThickness()
#            Moments calculation are bugged for now
#             elif func == 'Cm':
#                 val = self.__Cm()
            elif func == 'Lift':
                val = self.__Lift()
                dfdiAoA = self.__dLift_diAoA()
                dpJdpTwist = self.__dpLift_dpTwist()
                dpJdpAoA   = self.__dpLift_dpAoA()
                dpJdpThick = self.__dpLift_dpThickness()
            elif func == 'Drag':
                val= self.__Drag()
                dfdiAoA = self.__dDrag_diAoA()
                dpJdpTwist = self.__dpDrag_dpTwist()
                dpJdpAoA   = self.__dpDrag_dpAoA()
                dpJdpThick = self.__dpDrag_dpThickness()
            else:
                val = None
            self.__func_values[i]=val
            self.__dfunc_diAoA[i]=dfdiAoA
            self.__dpJ_dpTwist[i]=dpJdpTwist
            self.__dpJ_dpAoA[i]  =dpJdpAoA
            self.__dpJ_dpThick[i]=dpJdpThick
        
        self.__display_info()
        self.set_computed(True)
                
    #-- Computation methods
    #-- Cl related methods
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
    
    def __dpCl_dpTwist(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        dlAoAdT=self.get_DlocalAoA_DTwist()
        
        dCl=0.0
        for i in range(self.__N):
            af=airfoils[i]
            dCl+=af.ClAlpha(localAoA[i],Mach)*af.get_Sref()*dlAoAdT[i]
        dCl/=self.get_Sref()
        
        return dCl
    
    def __dpCl_dpAoA(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        dlAoAdAoA= self.get_DlocalAoA_DAoA()
        
        dCl=0.0
        for i in range(self.__N):
            af=airfoils[i]
            dCl+=af.ClAlpha(localAoA[i],Mach)*af.get_Sref()*dlAoAdAoA[i]
        dCl/=self.get_Sref()
        
        return dCl
    
    def __dpCl_dpThickness(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        
        dCl=zeros(self.__N)
        for i in range(self.__N):
            af=airfoils[i]
            dCl[i]=af.dCl_dthickness(localAoA[i],Mach)*af.get_Sref()
        dCl/=self.get_Sref()
        
        return dCl
    
    #-- Cd related methods
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
    
    def __dpCd_dpTwist(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA     = self.get_iAoA()
        airfoils = self.get_airfoils()
        dlAoAdT=self.get_DlocalAoA_DTwist()
 
        dCd=0.0
        for i in range(self.__N):
            af=airfoils[i]
            AoA=localAoA[i]
            Cdloc=(af.ClAlpha(AoA,Mach)*iAoA[i]+af.CdAlpha(AoA,Mach))*dlAoAdT[i]
            dCd+=Cdloc*af.get_Sref()
        dCd/=self.get_Sref()
 
        return dCd
    
    def __dpCd_dpAoA(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA     = self.get_iAoA()
        airfoils = self.get_airfoils()
        dlAoAdAoA= self.get_DlocalAoA_DAoA()

        dCd=0.0
        for i in range(self.__N):
            af=airfoils[i]
            AoA=localAoA[i]
            Cdloc=(af.ClAlpha(AoA,Mach)*iAoA[i]+af.CdAlpha(AoA,Mach))*dlAoAdAoA[i]
            dCd+=Cdloc*af.get_Sref()
        dCd/=self.get_Sref()
 
        return dCd

    def __dpCd_dpThickness(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        iAoA     = self.get_iAoA()
        airfoils = self.get_airfoils()
        
        dCd=zeros(self.__N)
        for i in range(self.__N):
            af=airfoils[i]
            AoA=localAoA[i]
            dCd[i]=af.dCl_dthickness(AoA,Mach)*iAoA[i]+af.dCd_dthickness(AoA,Mach)
            dCd[i]*=af.get_Sref()
        dCd/=self.get_Sref()
        
        return dCd 
    
    #-- Cm related methods (bugged for now)
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
    
    def __dCm_diAoA(self):
        Mach     = self.get_OC().get_Mach()
        airfoils = self.get_airfoils()
        localAoA = self.get_localAoA()
        DlocalAoA_DiAoA = self.get_DlocalAoA_DiAoA()
        dCm=0.0
        for i in range(self.__N):
            af=airfoils[i]
            dCm+=af.CmAlpha(localAoA[i],Mach)*af.getSref()*DlocalAoA_DiAoA[i]
        dCm/=self.get_Sref()
        return dCm_diAoA
    
    def __dpCm_dpTwist(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        dlAoAdT=self.get_DlocalAoA_DTwist()
        
        dCm=0.0
        for i in range(self.__N):
            af=airfoils[i]
            dCm+=af.CmAlpha(localAoA[i],Mach)*af.get_Sref()*dlAoADT[i]
        dCm/=self.get_Sref()
        
        return dCm
    
    def __dpCm_dpAoA(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        dlAoAdAoA= self.get_DlocalAoA_DAoA()
        
        dCm=0.0
        for i in range(self.__N):
            af=airfoils[i]
            dCm+=af.CmAlpha(localAoA[i],Mach)*af.get_Sref()*dlAoAdAoA[i]
        dCm/=self.get_Sref()
        
        return dCm
    
    def __dpCm_dpThickness(self):
        Mach     = self.get_OC().get_Mach()
        localAoA = self.get_localAoA()
        airfoils = self.get_airfoils()
        dlAoAdAoA= self.get_DlocalAoA_DAoA()
        
        dCm=zeros(self.__N)
        for i in range(self.__N):
            af=airfoils[i]
            dCm[i]=af.dCm_dthickness(localAoA[i],Mach)*af.get_Sref()
        dCm/=self.get_Sref()
        
        return dCm
    
     #-- Lift related methods
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
    
    def __dpLift_dpTwist(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        dpCl_dpTwist=self.__dpCl_dpTwist()
        dpLift_dpTwist = Pdyn*Sref*dpCl_dpTwist
        return dpLift_dpTwist
    
    def __dpLift_dpAoA(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        dpCl_dpAoA=self.__dpCl_dpAoA()
        dpLift_dpAoA = Pdyn*Sref*dpCl_dpAoA
        return dpLift_dpAoA
    
    def __dpLift_dpThickness(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        dpCl_dpThickness=self.__dpCl_dpThickness()
        dpLift_dpThickness = Pdyn*Sref*dpCl_dpThickness
        return dpLift_dpThickness
        
    #-- Drag related methods
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
    
    def __dpDrag_dpTwist(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        dpCd_dpTwist = self.__dpCd_dpTwist()
        dpDrag_dpTwist = Pdyn*Sref*dpCd_dpTwist
        return dpDrag_dpTwist
    
    def __dpDrag_dpAoA(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        dpCd_dpAoA = self.__dpCd_dpAoA()
        dpDrag_dpAoA = Pdyn*Sref*dpCd_dpAoA
        return dpDrag_dpAoA
    
    def __dpDrag_dpThickness(self):
        Pdyn=self.get_OC().get_Pdyn()
        Sref=self.get_Sref()
        dpCd_dpThickness = self.__dpCd_dpThickness()
        dpDrag_dpThickness =  Pdyn*Sref*dpCd_dpThickness
        return dpDrag_dpThickness
    
    #-- Display methods
    def __display_info(self):
        print self.get_OC()
        print '\n*** aerodynamic functions and coefficients ***'
        for i, func in enumerate(self.__func_list):
            if func in ['Lift', 'Drag']:
                unit = '(N)'
            else:
                unit = ''
            print '  '+self.__func_list[i]+'\t=\t'+str(self.__func_values[i])+' '+unit 
