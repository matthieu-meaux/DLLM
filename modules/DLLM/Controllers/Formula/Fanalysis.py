# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
"""
The B{Fanalysis} module contains B{exe} method that will analyze the constitency of a given formula.
if @(formula) is used it means the analysis will also linearize the formula

@author: Arnaud BARTHET
"""

import Fanalysisyacc

def exe(formula):
    """
    Execute the Formula analysis (@ for derivation)
    """
    testexpr='__init__'
    while str(testexpr) != str(formula) :
        testexpr=formula
        formula=Fanalysisyacc.parse(formula)
    return formula

##==============================================================================
## Main
##------------------------------------------------------------------------------
#if __name__ == '__main__':
#    formula = raw_input('expr > ')
#    analyzed=exe(formula)
#    print analyzed
