import numpy
span=40.
gamma0=1.70524
N=50

theta_list = numpy.linspace(0.,numpy.pi,N+1)
r_list_eta = -0.5+0.5*(1.-numpy.cos(theta_list))

fid=open('ideal_gamma.dat','w')
for i,r in enumerate(r_list_eta):
    y=r*span
    gamma=numpy.sqrt(1-(2*r)**2)*gamma0
    fid.write(str(y)+' '+str(gamma)+'\n')
fid.close()
    
    