# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
"""
Formula analysis yacc

Describes the action to be done during formula analysis or during linearization of a formula.

@author: Arnaud BARTHET
"""
__version__=2.0

import yacc
import Fanalysislex
import math
tokens = Fanalysislex.tokens

# Parsing rules
# Priority
precedence = (
    ('left','PLUS','MINUS'),
    ('left','TIMES','DIVIDE'),
    ('right','POW'),
    ('right','UMINUS')
    )
# Rules
def p_statement_expr(t):
    '''statement : expr'''
    if t[1]=='ZERO':
        t[0]=str(0.)
    else :  
        t[0]=str(t[1])

# Definition of expression
def p_expr_opadd(t):
    '''expr       : expr PLUS term
                  | expr MINUS term
                  | ZERO PLUS expr
                  | ZERO MINUS expr
                  | expr PLUS ZERO
                  | expr MINUS ZERO
                  | ZERO PLUS ZERO
                  | ZERO MINUS ZERO
                  | term'''
                
    if (len(t) == 4):
        # PLUS
        if t[2] == '+'  : 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                t[0] = 'ZERO'
            elif t[1] == 'ZERO' :
                t[0] = str(t[3])
            elif t[3] == 'ZERO' :
                t[0] = str(t[1])
            else  :
                t[0] = str(t[1])+'+'+str(t[3])
        # MINUS
        elif t[2] == '-': 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                t[0] = 'ZERO'
            elif t[1] == 'ZERO' :
                t[0] = '-'+str(t[3])
            elif t[3] == 'ZERO' :
                t[0] = str(t[1])
            else  :
                t[0] = str(t[1])+'-'+str(t[3])
                
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of term
def p_term_opmult(t):
    '''term       : term TIMES fact
                  | term DIVIDE fact
                  | term TIMES ZERO
                  | term DIVIDE ZERO
                  | ZERO TIMES term
                  | ZERO DIVIDE term
                  | ZERO TIMES ZERO
                  | ZERO DIVIDE ZERO
                  | fact'''
                
    if (len(t) == 4):
        # TIMES
        if t[2] == '*'  : 
            if t[1] == 'ZERO' or t[3] == 'ZERO':
                t[0] = 'ZERO'
            else  :
                t[0] = str(t[1])+'*'+str(t[3])
        # DIVIDE
        elif t[2] == '/': 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                raise Exception,'Division by zero error'
            elif t[1] == 'ZERO' :
                t[0] = 'ZERO'
            elif t[3] == 'ZERO' :
                raise Exception,'Division by zero error'
            else  :
                t[0] = str(t[1])+'/'+str(t[3])
                
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of factor
def p_fact_opother(t):
    '''fact       : MINUS min
                  | PLUS min
                  | MINUS ZERO
                  | PLUS ZERO
                  | fact POW fact
                  | fact POW ZERO
                  | ZERO POW fact
                  | ZERO POW ZERO
                  | MAX LPAREN fact COMMA fact RPAREN
                  | MIN LPAREN fact COMMA fact RPAREN
                  | CMP LPAREN fact COMMA fact RPAREN
                  | EXP LPAREN fact RPAREN
                  | LOG LPAREN fact RPAREN
                  | COS LPAREN fact RPAREN
                  | SIN LPAREN fact RPAREN
                  | TAN LPAREN fact RPAREN
                  | ACOS LPAREN fact RPAREN
                  | ASIN LPAREN fact RPAREN
                  | ATAN LPAREN fact RPAREN
                  | INT LPAREN fact RPAREN
                  | ABS LPAREN fact RPAREN
                  | SQRT LPAREN fact RPAREN
                  | MAX LPAREN expr COMMA expr RPAREN
                  | MIN LPAREN expr COMMA expr RPAREN
                  | CMP LPAREN expr COMMA expr RPAREN
                  | EXP LPAREN expr RPAREN
                  | LOG LPAREN expr RPAREN
                  | COS LPAREN expr RPAREN
                  | SIN LPAREN expr RPAREN
                  | TAN LPAREN expr RPAREN
                  | ACOS LPAREN expr RPAREN
                  | ASIN LPAREN expr RPAREN
                  | ATAN LPAREN expr RPAREN
                  | INT LPAREN expr RPAREN
                  | ABS LPAREN expr RPAREN
                  | SQRT LPAREN expr RPAREN
                  | EXP LPAREN ZERO RPAREN
                  | LOG LPAREN ZERO RPAREN
                  | COS LPAREN ZERO RPAREN
                  | SIN LPAREN ZERO RPAREN
                  | TAN LPAREN ZERO RPAREN
                  | ACOS LPAREN ZERO RPAREN
                  | ASIN LPAREN ZERO RPAREN
                  | ATAN LPAREN ZERO RPAREN
                  | INT LPAREN ZERO RPAREN
                  | ABS LPAREN ZERO RPAREN
                  | SQRT LPAREN ZERO RPAREN
                  | min'''
    if (len(t) == 7):
        # MAX
        if t[1] == 'max' :
            t[0] = 'max('+str(t[3])+','+str(t[5])+')'
        # MIN
        if t[1] == 'min' :
            t[0] = 'min('+str(t[3])+','+str(t[5])+')'
        # CMP
        if t[1] == 'cmp' : 
            t[0] = 'cmp('+str(t[3])+','+str(t[5])+')'
            
    if (len(t) == 5):
        # EXP
        if t[1] == 'math.exp' :	
            if t[3] == 'ZERO' :
                t[0] = '1.0'
            else:
                t[0] = 'math.exp('+str(t[3])+')'
        # LOG
        if t[1] == 'math.log' :
            if t[3] == 'ZERO' :
                print 'log(0.0) undefined'
                raise Exception,'log(0.0) undefined'
            else:
                t[0] = 'math.log('+str(t[3])+')'
        # SIN
        if t[1] == 'math.sin' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.sin('+str(t[3])+')'
        # COS
        if t[1] == 'math.cos' :	
            if t[3] == 'ZERO' :
                t[0] = '1.0'
            else:
                t[0] = 'math.cos('+str(t[3])+')'
        # TAN
        if t[1] == 'math.tan' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.tan('+str(t[3])+')'
        # ASIN
        if t[1] == 'math.asin' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.asin('+str(t[3])+')'
        # ACOS
        if t[1] == 'math.acos' :	
            if t[3] == 'ZERO' :
                t[0] = str(math.pi/2)
            else:
                t[0] = 'math.acos('+str(t[3])+')'
        # ATAN
        if t[1] == 'math.atan' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.atan('+str(t[3])+')'
        # INT
        if t[1] == 'int' :
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'int('+str(t[3])+')'
        #ABS
        if t[1] == 'abs':
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'abs('+str(t[3])+')'
        #SQRT
        if t[1] == 'math.sqrt' : 
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.sqrt('+str(t[3])+')'
            
    if (len(t) == 4):
        # POWER
        if (t[2] == '**' or t[2] == '^') : 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                raise Exception,'zero pow zero error'
            elif t[1] == 'ZERO' :
                t[0] = 'ZERO'
            elif t[3] == 'ZERO' :
                t[0] = '1.'
            else  :	    
                t[0] = str(t[1])+'**'+str(t[3])
                
    if (len(t) == 3):
        # MINUS
        if t[1] == '-'  : 
            if t[2] == 'ZERO' :
                t[0] = 'ZERO'
            else :
                t[0] = '-'+str(t[2])
        # PLUS
        if t[1] == '+'  :
            if t[2] == 'ZERO' :
                t[0] = 'ZERO'
            else :
                t[0] = str(t[2])
                
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of init
def p_min_init(t):
    '''min        : NUMBER
                  | VARIABLE
                  | INTEGER
                  | DERIVATE ZERO
                  | LPAREN expr RPAREN
                  | LPAREN ZERO RPAREN
                  | DERIVATE LPAREN dexpr RPAREN
                  | DERIVATE LPAREN ZERO RPAREN'''
    if (len(t) == 5):
        if t[3] == 'ZERO' :
            t[0] = 'ZERO'
        else :	
            t[0] = '('+str(t[3])+')'
            
    if (len(t) == 4):	
        if t[2] == 'ZERO' :
            t[0] = 'ZERO'
        else :
            t[0] = '('+str(t[2])+')'
            
    if (len(t) == 3):	
        t[0] = 'ZERO'
        
    if (len(t) == 2):
        t[0] = str(t[1])

# DERIVATIVE COMPUTATION
def p_dexpr_opadd(t):
    '''dexpr     : expr2 PLUS term2
                 | expr2 MINUS term2
                 | ZERO PLUS expr2
                 | ZERO MINUS expr2
                 | expr2 PLUS ZERO
                 | expr2 MINUS ZERO
                 | ZERO PLUS ZERO
                 | ZERO MINUS ZERO
                 | dterm'''
    if (len(t) == 4):
        # PLUS
        if t[2] == '+'  : 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                t[0] = 'ZERO'
            elif t[1] == 'ZERO' :
                t[0] = '@('+str(t[3])+')'
            elif t[3] == 'ZERO' :
                t[0] = '@('+str(t[1])+')'
            else  :
                t[0] = '@('+str(t[1])+')+@('+str(t[3])+')'
        # MINUS
        if t[2] == '-': 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                t[0] = 'ZERO'
            elif t[1] == 'ZERO' :
                t[0] = '-@('+str(t[3])+')'
            elif t[3] == 'ZERO' :
                t[0] = '@('+str(t[1])+')'
            else  :
                t[0] = '@('+str(t[1])+')-@('+str(t[3])+')'
                
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of term
def p_dterm_opmult(t):
    '''dterm      : term2 TIMES fact2
                  | term2 DIVIDE fact2
                  | ZERO TIMES term2
                  | ZERO DIVIDE term2
                  | term2 TIMES ZERO
                  | term2 DIVIDE ZERO
                  | ZERO TIMES ZERO
                  | ZERO DIVIDE ZERO
                  | dfact'''

    if (len(t) == 4):
        # TIMES
        if t[2] == '*'  : 
            if t[1] == 'ZERO' or t[3] == 'ZERO':
                t[0] = 'ZERO'
            else  :
                t[0] = '@('+str(t[1])+')*'+str(t[3])+'+'+str(t[1])+'*@('+str(t[3])+')'
            
        # DIVIDE
        elif t[2] == '/': 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                raise Exception,'Division by zero error'
            elif t[1] == 'ZERO' :
                t[0] = 'ZERO'
            elif t[3] == 'ZERO' :
                raise Exception,'Division by zero error'
            else  :
                t[0] = '(@('+str(t[1])+')*'+str(t[3])+'-'+str(t[1])+'*@('+str(t[3])+'))'+'/('+str(t[3])+')^2'
                
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of dfactor
def p_dfact_opother(t):
    '''dfact      : MINUS min2
                  | PLUS min2
                  | MINUS ZERO
                  | PLUS ZERO
                  | fact2 POW fact2
                  | fact2 POW ZERO
                  | ZERO POW fact2
                  | ZERO POW ZERO
                  | MAX LPAREN fact2 COMMA fact2 RPAREN
                  | MIN LPAREN fact2 COMMA fact2 RPAREN
                  | CMP LPAREN fact2 COMMA fact2 RPAREN
                  | EXP LPAREN fact2 RPAREN
                  | LOG LPAREN fact2 RPAREN
                  | COS LPAREN fact2 RPAREN
                  | SIN LPAREN fact2 RPAREN
                  | TAN LPAREN fact2 RPAREN
                  | ACOS LPAREN fact2 RPAREN
                  | ASIN LPAREN fact2 RPAREN
                  | ATAN LPAREN fact2 RPAREN
                  | INT LPAREN fact2 RPAREN
                  | ABS LPAREN fact2 RPAREN
                  | SQRT LPAREN fact2 RPAREN
                  | MAX LPAREN expr2 COMMA expr2 RPAREN
                  | MIN LPAREN expr2 COMMA expr2 RPAREN
                  | CMP LPAREN expr2 COMMA expr2 RPAREN
                  | EXP LPAREN expr2 RPAREN
                  | LOG LPAREN expr2 RPAREN
                  | COS LPAREN expr2 RPAREN
                  | SIN LPAREN expr2 RPAREN
                  | TAN LPAREN expr2 RPAREN
                  | ACOS LPAREN expr2 RPAREN
                  | ASIN LPAREN expr2 RPAREN
                  | ATAN LPAREN expr2 RPAREN
                  | INT LPAREN expr2 RPAREN
                  | ABS LPAREN expr2 RPAREN
                  | SQRT LPAREN expr2 RPAREN
                  | EXP LPAREN ZERO  RPAREN
                  | LOG LPAREN ZERO  RPAREN
                  | COS LPAREN ZERO  RPAREN
                  | SIN LPAREN ZERO  RPAREN
                  | TAN LPAREN ZERO  RPAREN
                  | ACOS LPAREN ZERO  RPAREN
                  | ASIN LPAREN ZERO  RPAREN
                  | ATAN LPAREN ZERO  RPAREN
                  | INT LPAREN ZERO  RPAREN
                  | ABS LPAREN ZERO  RPAREN
                  | SQRT LPAREN ZERO RPAREN
                  | dmin'''
    if (len(t) == 7):
        # MAX
        if t[1] == 'max' :
#            t[0] = '@('+str(t[3])+')*int(min((max((1+'+str(t[3])+'-'+str(t[5])+'),0)),1))+@('+str(t[5])+')*int(min((max((1+'+str(t[5])+'-'+str(t[3])+'),0)),1))'
            t[0] = '@('+str(t[3])+')*(0.5+0.5*cmp('+str(t[3])+'-'+str(t[5])+',0.))+@('+str(t[5])+')*(0.5-0.5*cmp('+str(t[3])+'-'+str(t[5])+',0.))'
        # MIN
        if t[1] == 'min' :
#            t[0] = '@('+str(t[3])+')*int(min((max((1+'+str(t[5])+'-'+str(t[3])+'),0)),1))+@('+str(t[5])+')*int(min((max((1+'+str(t[3])+'-'+str(t[5])+'),0)),1))'
            t[0] = '@('+str(t[3])+')*(0.5-0.5*cmp('+str(t[3])+'-'+str(t[5])+',0.))+@('+str(t[5])+')*(0.5+0.5*cmp('+str(t[3])+'-'+str(t[5])+',0.))'
        # CMP
        if t[1] == 'cmp':
            t[0] = 'ZERO'
    if (len(t) == 5):
        # EXP
        if t[1] == 'math.exp' :	
            if t[3] == 'ZERO' :
                t[0] = '1.0'
            else:
                t[0] = '@('+str(t[3])+')*math.exp(('+str(t[3])+'))'
        # LOG
        if t[1] == 'math.log' :
            if t[3] == 'ZERO' :
                print 'log(0.0) undefined'
                raise Exception,'log(0.0) undefined'
            else:
                t[0] = '@('+str(t[3])+')/('+str(t[3])+')'
        # SIN
        if t[1] == 'math.sin' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = '@('+str(t[3])+')*math.cos('+str(t[3])+')'
        # TAN
        if t[1] == 'math.tan' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = '@('+str(t[3])+')*(1+(math.tan('+str(t[3])+'))^2)'
        # COS
        if t[1] == 'math.cos' :	
            if t[3] == 'ZERO' :
                t[0] = '1.0'
            else:
                t[0] = '@('+str(t[3])+')*(-1.0)*math.sin('+str(t[3])+')'

        # ASIN
        if t[1] == 'math.asin' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = '@('+str(t[3])+')/(math.sqrt(1-('+str(t[3])+')^2))'
        # ATAN
        if t[1] == 'math.atan' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = '@('+str(t[3])+')/(1+('+str(t[3])+')^2)'
        # ACOS
        if t[1] == 'math.acos' :	
            if t[3] == 'ZERO' :
                t[0] = str(math.pi/2)
            else:
                t[0] = '@(-1.0*'+str(t[3])+')/(math.sqrt(1-('+str(t[3])+')^2))'
        # INT
        if t[1] == 'int' :	
            t[0] = 'ZERO'
        # ABS 
        if t[1] == 'abs':
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = '@('+str(t[3])+')*cmp('+str(t[3])+',0)'
        #SQRT
        if t[1] == 'math.sqrt':
            if t[3] == 'ZERO' :
                raise Exception, 'Cannot compute derivative of SQRT at 0.'
            else:
                t[0] = '@('+str(t[3])+')/(2*math.sqrt('+str(t[3])+'))'
            
    if (len(t) == 4):
        # POWER
        if (t[2] == '**' or t[2] == '^')  : 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                raise Exception,'zero pow zero error'
            elif t[1] == 'ZERO' :
                t[0] = 'ZERO'
            elif t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else  :
                #t[0] = str(t[3])+'*@('+str(t[1])+')*'+str(t[1])+'**('+str(t[3])+'-1.)'
                t[0] = str(t[3])+'*@('+str(t[1])+')*'+str(t[1])+'**'+str(eval(t[3])-1)
            
    if (len(t) == 3):
        # MINUS
        if t[1] == '-'  : 
            if t[2] == 'ZERO' :
                t[0] = 'ZERO'
            else :
                t[0] = '-@('+str(t[2])+')'
        # PLUS
        if t[1] == '+'  : 
            if t[2] == 'ZERO' :
                t[0] = 'ZERO'
            else :
                t[0] = '@('+str(t[2])+')'
            
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of init
def p_dmin_init(t):
    '''dmin       : VARIABLE
                  | LPAREN expr2 RPAREN'''
    if (len(t) == 4):
        t[0] = '@('+str(t[2])+')'
    elif (len(t) == 2):
        t[0] = 'd'+str(t[1])

def p_dmin_init2(t):
    '''dmin       : NUMBER
                  | INTEGER'''
    t[0] = 'ZERO'


# Definition of expression2
def p_expr2_opadd(t):
    '''expr2      : expr2 PLUS term2
                  | expr2 MINUS term2
                  | ZERO PLUS expr2
                  | ZERO MINUS expr2
                  | expr2 PLUS ZERO
                  | expr2 MINUS ZERO
                  | ZERO PLUS ZERO
                  | ZERO MINUS ZERO
                  | term2'''
    if (len(t) == 4):
        # PLUS
        if t[2] == '+'  : 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                t[0] = 'ZERO'
            elif t[1] == 'ZERO' :
                t[0] = str(t[3])
            elif t[3] == 'ZERO' :
                t[0] = str(t[1])
            else  :
                t[0] = str(t[1])+'+'+str(t[3])
        # MINUS
        if t[2] == '-': 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                t[0] = 'ZERO'
            elif t[1] == 'ZERO' :
                t[0] = '-'+str(t[3])
            elif t[3] == 'ZERO' :
                t[0] = str(t[1])
            else  :
                t[0] = str(t[1])+'-'+str(t[3])
            
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of term
def p_term2_opmult(t):
    '''term2      : term2 TIMES fact2
                  | term2 DIVIDE fact2
                  | ZERO TIMES term2
                  | ZERO DIVIDE term2
                  | term2 TIMES ZERO
                  | term2 DIVIDE ZERO
                  | ZERO TIMES ZERO
                  | ZERO DIVIDE ZERO
                  | fact2'''
    if (len(t) == 4):
        # TIMES
        if t[2] == '*'  : 
            if t[1] == 'ZERO' or t[3] == 'ZERO':
                t[0] = 'ZERO'
            else  :
                t[0] = str(t[1])+'*'+str(t[3])
        # DIVIDE
        if t[2] == '/': 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                raise Exception,'Division by zero error'
            elif t[1] == 'ZERO' :
                t[0] = 'ZERO'
            elif t[3] == 'ZERO' :
                raise Exception,'Division by zero error'
            else  :
                t[0] = str(t[1])+'/'+str(t[3])
                
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of factor
def p_fact2_opother(t):
    '''fact2      : MINUS min2
                  | PLUS min2
                  | MINUS ZERO
                  | PLUS ZERO
                  | fact2 POW fact2
                  | fact2 POW ZERO
                  | ZERO POW fact2
                  | ZERO POW ZERO
                  | MAX LPAREN fact2 COMMA fact2 RPAREN
                  | MIN LPAREN fact2 COMMA fact2 RPAREN
                  | CMP LPAREN fact2 COMMA fact2 RPAREN
                  | EXP LPAREN fact2 RPAREN
                  | LOG LPAREN fact2 RPAREN
                  | COS LPAREN fact2 RPAREN
                  | SIN LPAREN fact2 RPAREN
                  | TAN LPAREN fact2 RPAREN
                  | ACOS LPAREN fact2 RPAREN
                  | ASIN LPAREN fact2 RPAREN
                  | ATAN LPAREN fact2 RPAREN
                  | INT LPAREN fact2 RPAREN
                  | ABS LPAREN fact2 RPAREN
                  | SQRT LPAREN fact2 RPAREN
                  | MAX LPAREN expr2 COMMA expr2 RPAREN
                  | MIN LPAREN expr2 COMMA expr2 RPAREN
                  | CMP LPAREN expr2 COMMA expr2 RPAREN
                  | EXP LPAREN expr2 RPAREN
                  | LOG LPAREN expr2 RPAREN
                  | COS LPAREN expr2 RPAREN
                  | SIN LPAREN expr2 RPAREN
                  | TAN LPAREN expr2 RPAREN
                  | ACOS LPAREN expr2 RPAREN
                  | ASIN LPAREN expr2 RPAREN
                  | ATAN LPAREN expr2 RPAREN
                  | INT LPAREN expr2 RPAREN
                  | ABS LPAREN expr2 RPAREN
                  | SQRT LPAREN expr2 RPAREN
                  | EXP LPAREN ZERO  RPAREN
                  | LOG LPAREN ZERO  RPAREN
                  | COS LPAREN ZERO  RPAREN
                  | SIN LPAREN ZERO  RPAREN
                  | TAN LPAREN ZERO  RPAREN
                  | ACOS LPAREN ZERO  RPAREN
                  | ASIN LPAREN ZERO  RPAREN
                  | ATAN LPAREN ZERO  RPAREN
                  | INT LPAREN ZERO  RPAREN
                  | ABS LPAREN ZERO  RPAREN
                  | SQRT LPAREN ZERO RPAREN
                  | min2'''
    if (len(t) == 7):
        # MAX
        if t[1] == 'max' :
            t[0] = 'max('+str(t[3])+','+str(t[5])+')'
        # MIN
        if t[1] == 'min' :
            t[0] = 'min('+str(t[3])+','+str(t[5])+')'
        # CMP
        if t[1] == 'cmp' : 
            t[0] = 'cmp('+str(t[3])+','+str(t[5])+')'
            
    if(len(t) == 5):
        # EXP
        if t[1] == 'math.exp' :	
            if t[3] == 'ZERO' :
                t[0] = '1.0'
            else:
                t[0] = 'math.exp('+str(t[3])+')'
        # LOG
        if t[1] == 'math.log' :
            if t[3] == 'ZERO' :
                print 'log(0.0) undefined'
                raise Exception,'log(0.0) undefined'
            else:
                t[0] = 'math.log('+str(t[3])+')'
        # SIN
        if t[1] == 'math.sin' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.sin('+str(t[3])+')'
        # TAN
        if t[1] == 'math.tan' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.tan('+str(t[3])+')'
        # COS
        if t[1] == 'math.cos' :	
            if t[3] == 'ZERO' :
                t[0] = '1.0'
            else:
                t[0] = 'math.cos('+str(t[3])+')'

        # ASIN
        if t[1] == 'math.asin' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.asin('+str(t[3])+')'
        # ATAN
        if t[1] == 'math.atan' :	
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'math.atan('+str(t[3])+')'
        # ACOS
        if t[1] == 'math.acos' :	
            if t[3] == 'ZERO' :
                t[0] = '1.5707963267948966'
            else:
                t[0] = 'math.acos('+str(t[3])+')'

        # INT
        if t[1] == 'int' :
            if t[3] == 'ZERO' :
                t[0] = 'ZERO'
            else:
                t[0] = 'int('+str(t[3])+')'
        # ABS
        if t[1] == 'abs' :
            if t[3] == 'ZERO' : 
                t[0] = 'ZERO'
            else:
                t[0] = 'abs('+str(t[3])+')'
        # SQRT
        if t[1] == 'math.sqrt' : 
            if t[3] == 'ZERO':
                t[0] = 'ZERO'
            else:
                t[0] = 'math.sqrt('+str(t[3])+')'
                
    if (len(t) == 4):
        # POWER
        if  (t[2] == '**' or t[2] == '^') : 
            if t[1] == 'ZERO' and t[3] == 'ZERO':
                raise Exception,'zero pow zero error'
            elif t[1] == 'ZERO' :
                t[0] = 'ZERO'
            elif t[3] == 'ZERO' :
                t[0] = '1.'
            else :
                t[0] = str(t[1])+'**'+str(t[3])
                
    if (len(t) == 3):
        # MINUS
        if t[1] == '-'  : 
            if t[2] == 'ZERO' :
                t[0] = 'ZERO'
            else :
                t[0] = '-'+str(t[2])
        elif t[1] == '+'  : 
            if t[2] == 'ZERO' :
                t[0] = 'ZERO'
            else :
                t[0] = str(t[2])
                
    if (len(t) == 2):
        t[0] = str(t[1])

# Definition of init
def p_min2_init(t):
    '''min2       : NUMBER
                  | VARIABLE
                  | INTEGER
                  | LPAREN expr2 RPAREN'''
    if (len(t) == 4):
        t[0] = '('+str(t[2])+')'
        
    elif (len(t) == 2):
        t[0] = str(t[1])

def p_error(t):
    print "Syntax error!"
    raise Exception,"Formula syntax error!"

# Create pointer lex and yacc
Fanalysisparser=yacc.yacc()
Fanalysislexer=Fanalysislex.create_lexer()
# Method to parse, used in the main program
def parse(data):
    # Useful for debug
    # return Fanalysisparser.parse(data,lexer=Fanalysislexer,debug=2)
    return Fanalysisparser.parse(data,lexer=Fanalysislexer)
