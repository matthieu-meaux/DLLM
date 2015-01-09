# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
"""
The B{Formula} module contains the B{Formula} class

the B{Formula} class allows:
    - The analysis of a formula and its evaluation by defining the token values
    - The linearization of a formula and the gradient evaluattion by defining the token values

@author: Arnaud BARTHET
@author: Matthieu MEAUX
@author: Francois Gallard
"""
import math
from math import *
import numpy
from numpy import array
import types

import Fanalysis
import Fanalysislex

class Formula:
    """
    The B{Formula} class allows:
        - The analysis of a formula and its evaluation by defining the token values
        - The linearization of a formula and the gradient evaluation by defining the token values
        
    @ivar __input_expr: Formula expression (stores fexpr value)
    @ivar __fgrad: Option to compute gradient (stores fgrad value)
    @ivar __fexpr: Formula expression after analysis
    @ivar __fgradexpr: Linearized formula expression
    @ivar __value: Formula value after evaluation
    @ivar __gradient: Formula gradient value after evaluation 
    """
    
    REPLACE_DOT_BY='_DOT_'
    
    def __init__(self,fexpr,fgrad=True):
        """
        B{Formula} constructor
        
        @param fexpr: Formula expression
        @type fexpr: String 
        @param fgrad: option to compute gradient (B{True}: gradient active, B{False}: gradient inactive)
        @type fgrad: Boolean
        """
        # Formula Expressions
        self.__input_expr = fexpr
        self.__fgrad      = fgrad
        self.__fexpr      = None
        self.__fgradexpr  = None

        # Evaluation results
        self.__value      = None
        self.__gradient   = None
        
        # Initialise Expressions from input expression
        self.__init_expressions()
        
    # Private methods
    def __repr__(self):
        """
        Method that overload standard print command for the object.
        Displays information about the formula object by using the following command: 
            - C{>>>print myFormula}
        """
        tokenlist=self.get_token_list()
        if self.__fgrad:
            for token in self.get_token_list():
                tokenlist.append('d'+token)
        
        info_string  = '\n--o0 Formula Information 0o--'
        info_string += '\n  Active gradient     : %s'%self.__fgrad
        info_string += '\n  Formula             : %s'%self.__fexpr
        info_string += '\n  Gradient formula    : %s'%self.__fgradexpr
        info_string += '\n  Token list          : %s'%tokenlist
        info_string += '\n--o0 ------------------- 0o--'
        return info_string
    
    def __init_expressions(self):
        """
        Method that sets the following attributes:
            - __fexpr
            - __fgradexpr
        """
        # Value formula
        self.__fexpr = Fanalysis.exe(self.__input_expr)
        self.__init_formula(fgrad=False)
        # Gradient formula
        if self.__fgrad:
            gradient_expr    = '@('+self.__fexpr+')'
            self.__fgradexpr = Fanalysis.exe(gradient_expr)
            self.__init_formula(fgrad=True)
            
    def __clean_expr_tokens(self,expr,tokens):
        """
        Cleans the tokens to remove dots that cannot be contained in functions arguments 
        and replaces it by self.REPLACE_DOT_BY class attricute.
        Also cleans the formula in the same way
        @param expr : the expression
        @param type expr: string
        @param tokens: the tokens list
        @param type tokens: string list 
        """
        cln_tokens=[]
        cln_expr=expr
        for t in  tokens:
            cln_token=t.replace('.',self.REPLACE_DOT_BY)
            cln_tokens.append(cln_token)
            cln_expr=cln_expr.replace(t,cln_token)
        return   cln_expr,cln_tokens
    
    def __init_formula(self,fgrad=False):
        """
        Initializes the formula function self.f and self.df from the eval expressions
        Creates the tokens
        @param fgrad : if true df, the gradient of the function, is initialized, f otherwise
        @param type fgrad : Bool
        """
        tokens=self.get_token_list(fgrad=fgrad)
        if  fgrad:
            expr=self.__fgradexpr
            name='df'
        else:
            name='f'
            expr=self.__fexpr
                    
        cln_expr,cln_tokens=self.__clean_expr_tokens(expr, tokens)
        #Writes the function source code
        func= """
def %s(self,%s,**kwargs):
    return %s""" % (name,', '.join([t+'=None' for t in cln_tokens]), str(cln_expr))
        #Execute the code
        exec func
        if fgrad:#Attaches the interpreted method to the instance self
            self.df=types.MethodType(df,self)
        else:
            self.f=types.MethodType(f,self)
         
    # Some accessors
    def get_formula(self):
        """
        Accessor that returns __fexpr attribute
        @return: __fexpr
        """
        return self.__fexpr

    def get_grad_formula(self):
        """
        Accessor that returns __fgradexpr attribute (if fgrad is True)
        @return: __fgradexpr
        """
        if not self.__fgrad:
            print "Warning: gradient formula not generated!"
        return self.__fgradexpr
            
    def get_value(self):
        """
        Accessor that returns __value attribute
        @return: __value
        """
        return self.__value
    
    def get_gradient(self):
        """
        Accessor that returns __gradient attribute
        @return __gradient
        """
        return self.__gradient
    
    # Set methods
    def set_fexpr(self,fexpr,fgrad=True):
        """
        set __input_expr and __fgrad attribute
        @param fexpr: formula expression
        @type fexpr: String
        @param fgrad: option to compute gradient (B{True}: gradient active, B{False}: gradient inactive)
        @type fgrad: Boolean
        """
        self.__input_expr = fexpr
        self.__fgrad      = fgrad
        self.__init_expressions()
        
    def set_grad(self,fgrad=True):
        """
        set __fgrad attribute
        @param fgrad: option to compute gradient (B{True}: gradient active, B{False}: gradient inactive)
        @type fgrad: Boolean
        """
        if self.__fgrad != fgrad:
            self.__fgrad = fgrad
            self.__init_expressions()
    
    # Methods
    def get_token_list(self,fgrad=False,tokentype='VARIABLE'):
        """
        Returns the list of Token for a given token type
        @param fgrad: option to get tokens from the formula expression or from the linearized formula expression
        @type fgrad: Boolean
        @param tokentype: token type to be extracted
        @type tokentype: String
        """
        if  fgrad:
            expr=self.__fgradexpr
        else:
            expr=self.__fexpr

        lexer=Fanalysislex.create_lexer()
        lexer.input(expr)
        tokenlist=[]

        while 1:
            tok=lexer.token()
            # Stop if there is no token left
            if not tok: break
            
            if tok.type == tokentype :
            # Put the token found into the list
                if tok.value not in tokenlist:
                    tokenlist.append(tok.value)

        return tokenlist
    
    def evaluate(self,VarDict):
        """
        Set __var_dict attributes and evaluate expressions.
        @param VarDict: Dictionary of token values
        @type VarDict: Dictionary
        """
        #Replace the tokens by cleaned tokens
        eval_dict={}
        for k,v in VarDict.items():
            eval_dict[k.replace('.',self.REPLACE_DOT_BY)]=v
        #Call the f or df method implemented in __init_formula()
        self.__value = self.f(**eval_dict)
        if self.__fgrad:
            self.__gradient = self.df(**eval_dict)
        
    def display(self):
        """
        Display information about the B{Formula} object
        """
        print self
