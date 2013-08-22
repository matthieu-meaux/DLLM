from scipy.interpolate import Rbf
from AeroElastAdj.ValidGrad.FDGradient import FDGradient
from AeroElastAdj.ValidGrad.FDSecondOrderCentered import FDSecondOrderCentered

class RBF_interpolator():
    EPSILON=1e-7
    
    def __init__(self,x,y,scale,epsilon=0.1):
        self.__x=x
        self.__y=y
        self.__scale=scale
        args=self.__get_args(x)
        args.append(y.T)
        fd2=FDSecondOrderCentered(self.EPSILON)
        self.__fd_grad=FDGradient(fd2,self.f)
        self.__rbfi=Rbf(*args,epsilon=epsilon,function="linear")

    def __get_args(self,x):
        args=[]
        for i in xrange(x.shape[1]):
            args.append(x[:,i]*self.__scale[i])
        return args
        
    def f(self,x):
        args=x*self.__scale
        return self.__rbfi(*args)
    
    def df(self,x,comp=0):
        return self.__fd_grad.grad_f(x)[comp]
    