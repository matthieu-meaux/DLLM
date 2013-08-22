# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @version: 1.0
# @author: François Gallard

# - Local imports -
from AeroElastAdj.wing.wing_manufacture import Wing_manufacture

class Wing_geom_handler():
        
    DATA_PLANFORM='Wing geometry.Planform'
    DATA_WINGSPAN='Wing geometry.Wingspan'
    DATA_ROOT_CHORD='Wing geometry.Root chord'
    DATA_BREAK_CHORD='Wing geometry.Break chord'
    DATA_TIP_CHORD='Wing geometry.Tip chord'
    DATA_ROOT_THICKNESS='Wing geometry.Root thickness'
    DATA_BREAK_THICKNESS='Wing geometry.Break thickness'
    DATA_TIP_THICKNESS='Wing geometry.Tip thickness'
    DATA_N='Wing geometry.Number of sections'
    DATA_BREAK_PERCENT='Wing geometry.Break %'
    
    def __init__(self,data):
        self.__data=data
        
    def build_wing_geometry(self):
        data=self.__data
        planform=data[self.DATA_PLANFORM]
        manufacture=Wing_manufacture(planform)
        
        N=data[self.DATA_N]
        geom_parameters=self.get_geom_parameters()
        
        wing_geometry=manufacture.build_Discrete_wing(N,geom_parameters,twist_law=None)
        
        return wing_geometry
    
    def get_geom_parameters(self):
        data=self.__data
        planform=data[self.DATA_PLANFORM]
        
        geom_parameters={}
        geom_parameters["wingspan"]=data[self.DATA_WINGSPAN]
        geom_parameters['root_chord']=data[self.DATA_ROOT_CHORD]
        geom_parameters['root_height']=data[self.DATA_ROOT_THICKNESS]
        geom_parameters['tip_height']=data[self.DATA_TIP_THICKNESS]
        geom_parameters['box_chord_ratio'] = 2.
        if planform =="Broken":
            geom_parameters['break_chord']=data[self.DATA_BREAK_CHORD]
            geom_parameters['tip_chord']=data[self.DATA_TIP_CHORD]
            geom_parameters['break_percent']=data[self.DATA_BREAK_PERCENT]
            geom_parameters['break_height']=data[self.DATA_BREAK_THICKNESS]
            
        return geom_parameters