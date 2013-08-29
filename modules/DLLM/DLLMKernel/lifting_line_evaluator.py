# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @author: François Gallard

from MDOTools.Components.Wing.wing_manufacture import Wing_manufacture
from MDOTools.Components.Wing.wing_geom_handler import Wing_geom_handler

from OpenDACE.database.Database import Database

from DLLM.polarManager.analyticAirfoil import AnalyticAirfoil
from DLLM.polarManager.airfoilPolar import AirfoilPolar

import math
from numpy import concatenate,zeros,ones,random,array,linalg
import numpy

class Lifting_line_evaluator():
    DATA_AERO_MODEL='Airfoils.Airfoil aerodynamic model'
    DATA_AOA0='Airfoils.AoA0'
    DATA_CD0='Airfoils.Cd0'
    DATA_CM='Airfoils.Cm'
    DATA_POLAR_DB_FILE='Airfoils.XML polar database'
    
    def __init__(self,data):
        self.__data=data
        self.__n=data[Wing_geom_handler.DATA_N]
        self.__thickness_active=False
        self.__N_iter_a_parameter=False
        self.__llw=self.__build_wing(data)
        self.__modes=None
        self.__Dtwist_Dmode=None
        self.__nb_modes=0
        
    def is_thickness_active(self):
        return self.__thickness_active
    
    def is_alpha_a_parameter(self):
        return self.__alpha_is_a_parameter
    
    def is_N_iter_a_parameter(self):
        return self.__N_iter_a_parameter
    
    def eval(self,data,DVH,alpha,Mach,noise_level=0.1):
        self.update_geom_from_DVH(DVH,data)
        wing=self.__llw
        
        if self.is_alpha_a_parameter():
            alpha=self.get_alpha_from_DVH(DVH)
            
        print "Evaluating aerodynamic coefficients in conditions :"
        print "Mach = "+str(Mach)
        print "Angle of Attack = "+str(alpha*180./math.pi)+" deg"
        
        data['y'] = {}
        data['dy'] = {}
        data['ddy'] = {}
        data['y']['Cl'] = wing.Cl(alpha,beta=0.,Mach=Mach)
        data['y']['Cd'] = wing.Cd(alpha,beta=0.,Mach=Mach)
        data['y']['Cm'] = wing.Cm(alpha,beta=0.,Mach=Mach)
        data['y']['LoD'] = data['y']['Cl']/data['y']['Cd']
        
        data['convergence_history']=wing.get_convergence_history()
        
        dClDTw= wing.DCl_DTwist(alpha,beta=0.,Mach=Mach)
        dCdDTw= wing.DCd_DTwist(alpha,beta=0.,Mach=Mach)
        dCmDTw= wing.DCm_DTwist(alpha,beta=0.,Mach=Mach)
        dClDTw=self.correct_mode_grad(dClDTw)
        dCdDTw=self.correct_mode_grad(dCdDTw)
        dCmDTw=self.correct_mode_grad(dCmDTw)
        
        dy={}
        if not self.is_thickness_active():
            dy['Cl'] = dClDTw.tolist()
            dy['Cd'] = dCdDTw.tolist()
            dy['Cm'] = dCmDTw.tolist()
            dy['LoD'] = (dClDTw/data['y']['Cd']-dCdDTw*data['y']['Cl']/(data['y']['Cd'])**2).tolist()
        else:
            dClDTh= wing.DCl_DThickness(alpha,beta=0.,Mach=Mach)
            dCdDTh= wing.DCd_DThickness(alpha,beta=0.,Mach=Mach)
            dCmDTh= wing.DCm_DThickness(alpha,beta=0.,Mach=Mach)
            dy['Cl'] = concatenate((dClDTw,dClDTh)).tolist()
            dy['Cd'] = concatenate((dCdDTw,dCdDTh)).tolist()
            dy['Cm'] = concatenate((dCmDTw,dCmDTh)).tolist()
            dy['LoD'] = (concatenate((dClDTw,dClDTh))/data['y']['Cd']\
                                 -concatenate((dCdDTw,dCdDTh))*(data['y']['Cl']/(data['y']['Cd'])**2)).tolist()
        
        if self.is_alpha_a_parameter():
            dClDAoA= wing.DCl_DAoA(alpha,beta=0.,Mach=Mach)
            dCdDAoA= wing.DCd_DAoA(alpha,beta=0.,Mach=Mach)
            dCmDAoA= wing.DCm_DAoA(alpha,beta=0.,Mach=Mach)
            dy['Cl'] .append(dClDAoA)
            dy['Cd'] .append(dCdDAoA)
            dy['Cm'] .append(dCmDAoA)
            dy['LoD'].append(dClDAoA/data['y']['Cd']-dCdDAoA*data['y']['Cl']/(data['y']['Cd'])**2)
        
        dy_sorted=self.__sort_derivatives(DVH,dy)
            
        for k,v in dy_sorted.items():
            data['dy'][k]=v
            
        if self.is_N_iter_a_parameter():
            data['y']['R'] = wing.get_convergence_history()
            data['dy']['R'] = zeros(DVH.get_nvar()).tolist()
            
        self.compute_adj_corr(data,alpha,beta=0.,Mach=Mach)
        
        if noise_level!=0.:
            self.noise_grads(data['dy'],noise_level)
        
        self.__write_output(data,alpha,Mach=Mach)
        
    def __sort_derivatives(self,DVH,dy):
        names=DVH.get_x_names()
        twist_names=[]
        thick_names=[]
        for i in range(self.__n):
            twist_names.append('twist'+str(i))
            thick_names.append('thick'+str(i))
        
        if self.is_use_modes():
            modes_names=[]
            for i in range(self.__nb_modes):
                modes_names.append('mode'+str(i))
        
        dy_sorted={}
        
        for func,dy_value in dy.items():
            dy_sorted[func]=[]
            for var in names:
                if var in twist_names:
                    dy_sorted[func].append(dy_value[twist_names.index(var)])
                elif var in thick_names:
                    dy_sorted[func].append(dy_value[len(twist_names)+thick_names.index(var)])
                elif self.is_use_modes() and var in modes_names:
                    dy_sorted[func].append(dy_value[modes_names.index(var)])
                elif var=='AoA':
                    if self.is_thickness_active():
                        if not self.is_use_modes():
                            dy_sorted[func].append(dy_value[2*self.__n])
                        else:
                            dy_sorted[func].append(dy_value[self.__nb_modes+self.__n])
                    else:
                        if not self.is_use_modes():
                            dy_sorted[func].append(dy_value[self.__n])
                        else:
                            dy_sorted[func].append(dy_value[self.__nb_modes])
                else:
                    dy_sorted[func].append(0.)
        return dy_sorted
            
    def is_use_modes(self):
        return self.__nb_modes>0
     
    def __write_output(self,data,alpha,Mach):
        file_msg="                    oOo Lifting line results oOo\n"
        file_msg+="\n ** Evaluation conditions ** \n\n"
        file_msg+="Angle of attack = "+str(alpha*180./math.pi)+" deg \n"
        file_msg+="Mach = "+str(Mach)+" \n"
        file_msg+="\n ** Aerodynamic coefficients ** \n\n"
        for funcname in data['y'].keys():
            msg=funcname+ ' = '+str(data['y'][funcname])
            file_msg+=msg+'\n'
        file_msg+="\n ** Aerodynamic gradients ** \n\n"
        for funcname in data['dy'].keys():
            msg='d'+funcname+'/dX = '+str(data['dy'][funcname])
            file_msg+=msg+'\n\n'
            
        file_msg+="\n ** Convergence error ** \n\n"  
        file_msg+="The adjoint correction to a function is a first order estimation of the convergence error of a function. See Arnaud Barhtet's PhD thesis.\n"
        file_msg+="Adjoint correction to Cl = "+str(data['Cl_adj_corr'])+"\n"
        file_msg+="Adjoint correction to Cd = "+str(data['Cd_adj_corr'])+"\n"
        file_msg+="Adjoint correction to Cm = "+str(data['Cm_adj_corr'])+"\n"
        
        print file_msg
        
        fid=open('eval.info','w')
        fid.write(file_msg)
        fid.close()
        data['evathicknesslinfo'] = 'eval.info'
        
        gamma_file_name='Circulation.dat'
        self.__llw.write_gamma_to_file(gamma_file_name,alpha=alpha,beta=0.,Mach=Mach)
        data['circulation'] = gamma_file_name
        
        
    def noise_grads(self,dy,level):
        for func,dy_value in dy.items():
            
            dya=array(dy_value)
            print "WARNING Noising artificially gradients "+str(func)+" value init="+str( dya)
            noised_dy=self.noise(dya, level).tolist()
            dy[func]=noised_dy
            print "noised value =                          "+str(array(noised_dy))
            print "err="+str(linalg.norm(dya-noised_dy)/linalg.norm(dya))
        return dy
            
    def noise(self,grad,level):
        
        ns=random.random_sample(grad.shape)
        ns*=level*linalg.norm(grad)/linalg.norm(ns)
        #print "noise ",ns
        return ns+grad
    
    def compute_adj_corr(self,data,alpha,beta,Mach):
        Cl_adj_corr=self.__llw.adj_correction_Cl(alpha,beta,Mach)
        Cd_adj_corr=self.__llw.adj_correction_Cd(alpha,beta,Mach)
        Cm_adj_corr=self.__llw.adj_correction_Cm(alpha,beta,Mach)
        
        data['Cl_adj_corr']=Cl_adj_corr
        data['Cd_adj_corr']=Cd_adj_corr
        data['Cm_adj_corr']=Cm_adj_corr
        
    def set_twist_from_modes(self,dvh,data):
        if "mode0" in dvh.get_DVId():
            twistLaw = zeros(self.__n)
            self.__modes=[]
            for p in range(len(dvh.get_DVId())):
                var="mode"+str(p)
                if dvh.is_var_in_handler(var):
                    self.__modes.append(dvh.get_var(var))
                    #print "Mode coeff "+str(p)+" = "+str(dvh.get_var(var))
            self.__Dtwist_Dmode=numpy.zeros((self.__n,len(self.__modes)))
            modes_num=numpy.loadtxt(data['Modes'])
            #print "Read modes = "+str(modes_num)
            self.__nb_modes=len(self.__modes)
            for p,mode in enumerate(self.__modes):
                twistLaw+=mode*modes_num[p,:]
                print "Twist = "+str(twistLaw)
                self.__Dtwist_Dmode[:,p]=modes_num[p,:]
            #print "D twist / D mode = "+str(self.__Dtwist_Dmode)
            self.__llw.set_twist(twistLaw)
            
    def correct_mode_grad(self,dy_twist):
        if self.__Dtwist_Dmode is not None:
            dy_corr= numpy.dot(numpy.atleast_2d(dy_twist),self.__Dtwist_Dmode)[0,:]
            return dy_corr
        return dy_twist
        
    def set_twit_from_dvh(self,dvh):
        twistLaw = zeros(self.__n)
        for i in range(self.__n):
            var="twist"+str(i)
            if dvh.is_var_in_handler(var):
                twistLaw[i]=dvh.get_var(var)
            else:
                print "Unparameterized twist no : "+str(i)
        self.__llw.set_twist(twistLaw)
        
    def update_geom_from_DVH(self,dvh,data):
        if "mode0" in dvh.get_DVId():
            self.set_twist_from_modes(dvh,data)
        else:
            self.set_twit_from_dvh(dvh)
        
        self.__thickness_active = dvh.is_var_in_handler("thick0")
        if self.__thickness_active:
            thickness = zeros(self.__n)
            for i in range(self.__n):
                thickness[i]=dvh.get_var("thick"+str(i))
            self.__llw.set_relative_thickness(thickness)
        else:
            thickness = 0.15*ones(self.__n)
            self.__llw.set_relative_thickness(thickness)
            print "WARNING thickness not active, using 0.15 by default"
            
        self.__alpha_is_a_parameter=dvh.is_var_in_handler("AoA")
        
        self.__N_iter_a_parameter=dvh.is_var_in_handler("Iterations_number")
        if self.__N_iter_a_parameter:
            nit=int(dvh.get_var("Iterations_number"))
            print "WARNING : OVERLOADING STOP CRITERIA WITH ITERATION NUMBER = "+str(nit)
            self.__llw.set_stop_criteria(n_it=nit)
            
    def get_alpha_from_DVH(self,dvh):
        return dvh.get_var("AoA")
        
    def __build_wing(self,data):
        wgh=Wing_geom_handler(data)
        wing_geometry=wgh.build_wing_geometry()
        lifting_line_wing=self.__build_Lifting_line_wing(data,wing_geometry)
        return lifting_line_wing
    
    def __build_Lifting_line_wing(self,data,wing_geometry):
        planform=data[Wing_geom_handler.DATA_PLANFORM]
        manufacture=Wing_manufacture(planform)
        
        Lref=data[Wing_geom_handler.DATA_ROOT_CHORD]
        Sref=data[Wing_geom_handler.DATA_WINGSPAN]*Lref
        
        airfoil_type=data[self.DATA_AERO_MODEL]
        if airfoil_type == 'Linear':
            AoA0=data[self.DATA_AOA0]
            Cd0=data[self.DATA_CD0]
            Cm=data[self.DATA_CM]
            reference_airfoil=self.__build_airfoil(airfoil_type,Sref=Sref,Lref=Lref,AoA0=AoA0,Cd0=Cd0,Cm=Cm,database=None)

        elif airfoil_type == "Polar":
            db=Database(data[self.DATA_POLAR_DB_FILE])
            reference_airfoil=self.__build_airfoil(airfoil_type,Sref=Sref,Lref=Lref,AoA0=None,Cd0=None,Cm=None,database=db)
        else:
            raise Exception, "Unknown airfoil type : "+str(airfoil_type)
        
        lifting_line_wing=manufacture.build_Lifting_line_wing(wing_geometry,reference_airfoil)
        
        if data['Numerical parameters.Stop criteria']=='Residual decrease':
            res=data['Numerical parameters.Decrease ratio']
            lifting_line_wing.set_stop_criteria(residual=res)
            
        elif data['Numerical parameters.Stop criteria']=='Iterations number':
            nit=data['Numerical parameters.Iterations number']
            lifting_line_wing.set_stop_criteria(n_it=nit)
                                                
        return lifting_line_wing
    
             
    def __build_airfoil(self,airfoil_type,Sref,Lref,AoA0=-2.,Cd0=0.01,Cm=-0.2,database=None):
        if airfoil_type=="Linear":
            degToRad=math.pi/180.
            return AnalyticAirfoil(Sref,Lref,AoA0*degToRad,Cd0,2.*math.pi,Cm)
        
        if airfoil_type=="Polar":
            return AirfoilPolar(Sref,Lref,database,relative_thickness=0.15)
        
        return None