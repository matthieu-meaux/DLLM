from DLLM.DLLMGeom.wing_param import Wing_param
from MDOTools.OC.operating_condition import OperatingCondition
import numpy
import string

OC=OperatingCondition('cond1',atmospheric_model='simple')

OC.set_altitude(3000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)

wing_param=Wing_param('test_param',geom_type='Elliptic',n_sect=20)
wing_param.build_wing()
wing_param.set_value('test_param.span',40.)
wing_param.set_value('test_param.root_chord',4.)
wing_param.set_value('test_param.root_height',0.8)
wing_param.set_value('test_param.tip_height',0.1)
wing_param.build_linear_airfoil(OC, AoA0=0., Cm0=-0.1, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

airfoils=wing_param.get_linked_airfoils()
Cla=airfoils[0].ClAlpha(0.0,0.3)

print "Cla=",Cla

#list_file = ['SimpleRect_3.00E-01.dat','SimpleRect_6.00E-01.dat','SimpleRect_8.00E-01.dat']
list_file = ['Elliptic_3.00E-01.dat']
e=0.70
AR = wing_param.get_AR()
print 'AR=',AR

for i,file_name in enumerate(list_file):
    fid=open(file_name,'r')
    lines=fid.readlines()
    fid.close()
    
    fid=open('Cdp_'+file_name,'w')
    for j,line in enumerate(lines):
        words=string.split(line)
        new_words=words
        if j == 1:
            new_words = ['#','AoA','Cl','Cdp','Cl_theory','Cdp_theory','Cdp_elliptic']
        if len(words)> 0:
            if words[0] != '#':
                AoA=eval(words[0])
                Cl=eval(words[7])
                Cdp=eval(words[9])
                Cl_theory=AR*Cla/(AR+2.)*AoA*numpy.pi/180.
                Cdp_theory=Cl_theory**2/(numpy.pi*e*AR)
                Cdp_ell = Cl_theory**2/(numpy.pi*AR)
                new_words=[str(AoA),str(Cl),str(Cdp),str(Cl_theory),str(Cdp_theory),str(Cdp_ell)]
        
        line = string.join(new_words,' ')
        fid.write(line+'\n')
    fid.close()

fid=open('theory_iAoA.dat','w')
for i in xrange(20):
    fid.write(str(i)+'\t'+str(-Cl/(numpy.pi*AR))+'\t'+str(-Cl/(numpy.pi*e*AR))+'\n')
fid.close()