import numpy

N=50

r_list_lin = numpy.linspace(0.,1.,N+1) 

theta_list = numpy.linspace(0.,numpy.pi,N+1)
r_list_cos = 0.5*(1.-numpy.cos(theta_list))

fid=open('distrib.dat','w')
for i in xrange(N+1):
    fid.write(str(0.)+' '+str(r_list_lin[i])+' '+str(r_list_cos[i])+'\n')
fid.close()
    

    