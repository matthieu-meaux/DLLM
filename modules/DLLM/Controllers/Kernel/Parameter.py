# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus

from re import compile
from string import join
from DLLM.Controllers.Kernel.BaseController import BaseController
from DLLM.Controllers.Formula.Formula import Formula

class Parameter(BaseController):
    """
    PParameter Class
    """
    #--Class Variables
    CLASS_MSG = 'PParameter'
    ERROR_MSG = 'ERROR '+CLASS_MSG+'.'
    
    #--Constructor
    def __init__(self,PBCManager,Id,fexpr,AliasDict=None):

        BaseController.__init__(self, PBCManager, Id, value=0., BCType='Parameter')
        
        self.__fexpr       = None
        self.__AliasDict   = None
        self.__FormulaObj  = None
        self.__dep_id_list = []
        self.__dep_dict    = {}
        
        self.set_fexpr(fexpr,AliasDict=AliasDict)
        
    def __repr__(self):
        """
        Display some information about the variable
        """
        info_string =BaseController.__repr__(self)
        info_string+="\n   Formula         : "+self.__FormulaObj.get_formula()
        info_string+="\n   Grad Formula    : "+self.__FormulaObj.get_grad_formula()
        info_string+="\n   Value           :%24.16e"%self.get_value()
        info_string+="\n   Gradient        : "+str(self.get_gradient())
        return info_string
    
    def __replace_aliases(self):
        if self.__AliasDict is not None:
            alias_keys = self.__AliasDict.keys()
            sep = compile('([\^\*\+\[\]/\(\)-])')
            fexpr = sep.split(self.__fexpr)
            for indice in range(len(fexpr)):
                if fexpr[indice].strip() in alias_keys:
                    fexpr[indice] = self.__AliasDict[fexpr[indice].strip()]
            self.__fexpr=join(fexpr,'')
            
    def set_fexpr(self,fexpr,AliasDict=None):
        self.__fexpr       = fexpr
        self.__AliasDict   = AliasDict
        
        self.__replace_aliases()
        
        fgrad=self.is_gradient_active()
        
        self.__FormulaObj  = Formula(self.__fexpr,fgrad=fgrad) 
        
        dep_id_list = self.__FormulaObj.get_token_list()
        
        self.__update_dep_list(dep_id_list)
        
    def get_expr(self):
        return self.__fexpr
        
    def get_formula(self):
        return self.__FormulaObj.get_formula()
    
    def get_grad_formula(self):
        return self.__FormulaObj.get_grad_formula()
    
    def __update_dep_list(self,dep_id_list):
        BC_manager=self.get_manager()
        
        for Id in self.__dep_id_list:
            if Id not in dep_id_list:
                pt=BC_manager.get_pt(Id)
                self.del_dependance(pt)
                del self.__dep_dict[Id]
                
        for Id in dep_id_list:
            if Id not in self.__dep_id_list:
                pt=BC_manager.get_pt(Id)
                self.add_dependance(pt)
                self.__dep_dict[Id]=pt
            
        self.__dep_id_list=dep_id_list
        
    def special_dependance_update(self,pt):
        ERROR_MSG=self.ERROR_MSG+'special_dependance_update: '
        Id=pt.get_id()
        if Id not in self.__dep_id_list:
            raise Exception,ERROR_MSG+'method cannot be applied to non existing id ('+str(Id)+')'
        self.add_dependance(pt)
        self.set_to_update()
        self.__dep_dict[Id]=pt
    
    def update_specific(self):
        eval_dict={}
        for Id in self.__dep_id_list:
            eval_dict[Id]     = self.__dep_dict[Id].get_value()
        if self.is_gradient_active():
            self.__FormulaObj.set_grad(fgrad=True)
            for Id in self.__dep_id_list:
                eval_dict['d'+Id] = self.__dep_dict[Id].get_gradient()
        else:
            self.__FormulaObj.set_grad(fgrad=False)
        try:
            self.__FormulaObj.evaluate(eval_dict)
        except:
            ERROR_MSG=self.ERROR_MSG+' Failed to evaluate parameter '+self.get_id()+' of expression: '+self.__fexpr
            raise Exception,ERROR_MSG
        
        value    = self.__FormulaObj.get_value()
        self.set_value(value,flag_updates=False)
        
        if self.is_gradient_active():
            gradient = self.__FormulaObj.get_gradient()
            self.set_gradient(gradient)
