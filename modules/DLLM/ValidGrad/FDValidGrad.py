# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding: iso-8859-1 -*-
# Copyright: Airbus
# @author Francois Gallard

from numpy.linalg import norm
from numpy import savetxt
from FDGradient import FDGradient
from FDSecondOrderCentered import FDSecondOrderCentered
from FDFirstOrderUpwind import FDFirstOrderUpwind

class FDValidGrad():
    
    def __init__(self,scheme_order,f_pointer,df_pointer,fd_step=1e-6):
        self.__fpointer=f_pointer
        self.__df_pointer=df_pointer
        
        if scheme_order==1:
            self.__fd_scheme=FDFirstOrderUpwind(fd_step)
            
        elif scheme_order==2:
            self.__fd_scheme=FDSecondOrderCentered(fd_step)
        else:
            raise Exception, "Scheme of order" +str(scheme_order)+" not available now."
        
        self.__fd_grad=FDGradient(self.__fd_scheme,f_pointer)
        
    def compare(self,x,treshold=1e-4, args=None,force_print=False, grad_files_prefix="",iprint=True):
        df_fd=self.__compute_fd_grad(x)
        
        if args is None:
            df=self.__df_pointer(x)
        else:
            df=self.__df_pointer(x,args)
            
        error = self.__compute_error(df_fd, df,treshold)
        ok = error < treshold
        
        if not ok or force_print:
            if iprint:
                print "Gradient not in bounds, error = "+str(error)
                print "df =\n"+str(df)
                print 'df finite diff =\n'+str(df_fd)
            #for i in xrange(len(x)):
            #    savetxt(grad_files_prefix+'df_analytic_'+str(i)+'.txt',df[:,:,i].T)
            #    savetxt(grad_files_prefix+'df_FD_'+str(i)+'.txt',df_fd[:,:,i].T)
        return ok
        
    def __compute_fd_grad(self,x):
        return self.__fd_grad.grad_f(x)
    
    def __compute_error(self,df_fd,df,treshold):
        nfd=norm(df_fd)
        
        if norm(df_fd)<treshold:#In case df = 0
            return norm(df)
        
        return norm(df_fd-df)/nfd