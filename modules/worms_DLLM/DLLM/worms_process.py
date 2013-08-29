# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @author: François Gallard

from worms.Workflow.worms_process_base import WormsProcessBase
from DLLM.DLLMKernel.lifting_line_evaluator import Lifting_line_evaluator
from OpenDACE.evalcontrol.dv_handler import DVHandler
from math import pi

class WormsProcess(WormsProcessBase):
    
    process_type='Evaluator'
    overloadable_inputs={'Mach':['Evaluation conditions.Mach','Float'],'AoA':['Evaluation conditions.AoA','Float']}
    
    DATA_AOA='Evaluation conditions.AoA'
    DATA_MACH='Evaluation conditions.Mach'
    DATA_DV='design_variables'
    
    def setup(self):
        self.__set_grad_hess()
        
    def __set_grad_hess(self):
        WF=self.get_workflow()
        WF.get_item( itemID = 'is_gradient_available', grammar_mode = 'input')._value=True
        WF.get_item( itemID = 'is_hessian_available' , grammar_mode = 'input')._value=False
        
    def run(self,data):
        print "oOo Lifting Line oOo"
        llw_evaluator = Lifting_line_evaluator(data)
        
        alpha=data[self.DATA_AOA]*pi/180.
        Mach=data[self.DATA_MACH]
        
        DVH=DVHandler()
        DVH.import_from_file(data[self.DATA_DV])
        
        llw_evaluator.eval(data,DVH,alpha,Mach,noise_level=0.)
        data['evalinfo']='eval.info'
        print "oOo END RUN Lifting Line oOo"
