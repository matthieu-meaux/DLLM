# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: EADS Innovation Works
# @version: 1.0
# @author: Andre BORTOLAZZI

from airfoil import Airfoil
import numpy as np 
from surrogate_coefs import SurrogateCoefs

class MetaAirfoil(Airfoil):
    """
    A meta model based on 2D CFD data
    Provide lift, drag and moment coefficients
    Provide sensitivities of aerodynamic surrogate function w.r.t. input parameters
    """
    def __init__(self, OC, airfoil_model, relative_thickness=.12, camber=0., Sref=1., Lref=1., sweep=0., surrogate_fcs=None):

        Airfoil.__init__(self, OC, Sref, Lref, relative_thickness, sweep, camber)
        if not surrogate_fcs:
            self.__coefs = SurrogateCoefs(airfoil_model)
        else:
            self.__coefs = surrogate_fcs
        self.__airfoil_model = airfoil_model
        
    def Cl(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal = Mach*np.cos(sweep)
        Cl = self.__coefs.meta_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)
        sweep_corr = np.cos(self.get_sweep())**2
        
        return Cl*sweep_corr
    
    def ClAlpha(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCl_dAoA = self.__coefs.grad_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[2]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCl_dAoA*sweep_corr
        
    def Cd(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        Cd = self.__coefs.meta_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)  
        #sweep_corr = np.cos(self.get_sweep())**2
        
        #return Cd*sweep_corr
        return Cd
    
    def CdAlpha(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCd_dAoA = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[2]
        #sweep_corr = np.cos(self.get_sweep())**2
        
        return dCd_dAoA#*sweep_corr
        
    def Cm(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        Cm = self.__coefs.meta_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)  
        sweep_corr = np.cos(self.get_sweep())**2
        
        return Cm*sweep_corr
    
    def dCl_dchi(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal = Mach*np.cos(sweep)
        Cl = self.__coefs.meta_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)
        
        dsweep=self.get_sweep_grad()
        sweep_corr = np.cos(self.get_sweep())**2
        dsweep_corr = -2.*np.sin(sweep)*np.cos(sweep)*dsweep
        
        dCl_dthick = self.gradCl(alpha, Mach)[0]
        dthick_dchi = self.get_rel_thick_grad()
        
        dCl_dMach = self.gradCl(alpha, Mach)[3]
        dMach_dchi = -Mach*np.sin(sweep)*dsweep
        
    
        return (dCl_dthick*dthick_dchi + dCl_dMach*dMach_dchi)*sweep_corr + Cl*dsweep_corr
        
    def dCd_dchi(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal = Mach*np.cos(sweep)
        Cd = self.__coefs.meta_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)
        
        #sweep_corr = np.cos(self.get_sweep())**2
        dsweep=self.get_sweep_grad()
        #dsweep_corr = -2.*np.sin(sweep)*np.cos(sweep)*dsweep
        
        dCd_dthick = self.gradCd(alpha, Mach)[0]
        dthick_dchi = self.get_rel_thick_grad()
        
        dCd_dMach = self.gradCd(alpha, Mach)[3]
        dMach_dchi = -Mach*np.sin(sweep)*dsweep
        
        #return (dCd_dthick*dthick_dchi + dCd_dMach*dMach_dchi)*sweep_corr + Cd*dsweep_corr
        return (dCd_dthick*dthick_dchi + dCd_dMach*dMach_dchi)
        
    def gradCl(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        gradCl = self.__coefs.grad_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)
        gradCl[0] = gradCl[0]*100.
        return gradCl
        
    def gradCd(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        gradCd = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)
        gradCd[0] = gradCd[0]*100.
        return gradCd
        
    def gradCm(self, alpha, Mach):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        gradCm = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)
        gradCm[0] = gradCm[0]*100.
        return gradCm
    
    def dCl_dthickness(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCl_dthic = self.__coefs.grad_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[0]*100.
        
        return dCl_dthic*sweep_corr
        
    def dCl_dcamber(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCl_dcamb = self.__coefs.grad_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[1]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCl_dcamb*sweep_corr
        
    def dCl_dMach(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCl_dM = self.__coefs.grad_Cl(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[3]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCl_dM*sweep_corr
        
    def dCd_dthickness(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCd_dthic = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[0]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCd_dthic*sweep_corr
        
    def dCd_dcamber(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCd_dcamb = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[1]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCd_dcamb*sweep_corr
        
    def dCd_dalpha(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCd_dAoA = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[2]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCd_dAoA*sweep_corr
        
    def dCd_dMach(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCd_dM = self.__coefs.grad_Cd(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[3]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCd_dM*sweep_corr
        
    def dCm_dthickness(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCm_dthic = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[0]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCm_dthic*sweep_corr
        
    def dCm_dcamber(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCm_dcamb = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[1]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCm_dcamb*sweep_corr
        
    def dCm_dalpha(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCm_dAoA = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[2]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCm_dAoA*sweep_corr
        
    def dCm_dMach(self):
        sweep=self.get_sweep()
        Mach_normal= Mach*np.cos(sweep)
        dCm_dM = self.__coefs.grad_Cm(self.get_rel_thick()*100., self.get_camber(), alpha, Mach_normal)[3]
        sweep_corr = np.cos(self.get_sweep())**2
        
        return dCm_dM*sweep_corr
        
    def get_scaled_copy(self, OC=None, Sref=None,Lref=None, rel_thick=None, camber=None, sweep=None):
        if Sref is None:
            Sref=self.get_Sref()
        if Lref is None:
            Lref=self.get_Lref()
        if rel_thick is None:
            rel_thick = self.get_rel_thick()
        if camber is None:
            camber = self.get_camber()
        if sweep is None:
            sweep = self.get_sweep()
        if OC is None:
            OC = self.get_OC()
        return MetaAirfoil(OC, self.__airfoil_model, relative_thickness=rel_thick, camber=camber, Sref=Sref, Lref=Lref, sweep=sweep, surrogate_fcs=self.__coefs)
