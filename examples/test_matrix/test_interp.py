import numpy
N = 8

gamma=numpy.array([1.,2.,3.,4.,4.,3.,2.,1.])

Mat = numpy.zeros([N+1,N])
Mat[0:N,:]   = numpy.diag(numpy.ones([N]))
Mat[N,:]     = 0.0
Mat[1:N+1,:]-= numpy.diag(numpy.ones([N]))

dgamma=numpy.dot(Mat,gamma)

print dgamma

file=open('interp.dat','w')
for i in xrange(N):
    file.write(str(i)+'\t'+str(dgamma[i])+'\n')
file.close()

Mat2 = numpy.zeros([N+1,N])
Mat2[0:N,:]   = numpy.diag(numpy.ones([N]))
Mat2[N,:]     = 0.0
Mat2[1:N+1,:]-= numpy.diag(numpy.ones([N])) 

Mat2[0,0] = -2
Mat2[0,1] = 3
Mat2[0,2] = -1 

Mat2[-1,-1] = 2
Mat2[-1,-2] = -3
Mat2[-1,-1] = 1 

dgamma=numpy.dot(Mat2,gamma)

file=open('interp2.dat','w')
for i in xrange(N):
    file.write(str(i)+'\t'+str(dgamma[i])+'\n')
file.close()