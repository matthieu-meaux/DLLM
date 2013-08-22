from scipy.interpolate import bisplrep,bisplev

class Spline2D_interpolator():
    
    def __init__(self,x,y,scale):
        self.__x=x
        self.__y=y
        self.__scale=scale
        self.__tck = bisplrep(x=self.__x[:,0]*scale[0],y=self.__x[:,1]*scale[1],z=self.__y, kx=5, ky=5)
        
    def f(self,x):
        return bisplev(x[0]*self.__scale[0],x[1]*self.__scale[1],self.__tck)
    
    def df(self,x,comp=0):
        dx=1-comp
        dy=comp
        grad=bisplev(x[0]*self.__scale[0],x[1]*self.__scale[1],self.__tck,dx=dx,dy=dy)
        return grad*self.__scale[comp]

