# Imports
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.OC.operating_condition import OperatingCondition
import numpy
import string

OC=OperatingCondition('cond1',atmospheric_model='simple')
#OC.set_Mach(0.8)
#OC.set_Mach(0.6)
#OC.set_AoA(3.5)
AoA_list = [float(xx) for xx in range(0, 10)]
AoA_list = [9.]
#Mach_list = [float(xx)/10. for xx in range(3, 9)]
Mach_list = [0.2]#,0.6,0.8]

OC.set_altitude(10000.)
OC.set_T0_deg(15.)
OC.set_P0(101325.)
OC.set_humidity(0.)

wing_param=Wing_param('test_param',geom_type='Rectangular',n_sect=51)
wing_param.build_wing()
wing_param.set_value('test_param.span',40.)
#wing_param.set_value('test_param.sweep',34.)
wing_param.set_value('test_param.sweep',0.)
wing_param.set_value('test_param.root_chord',1.0)
wing_param.set_value('test_param.root_height',0.1)
wing_param.set_value('test_param.tip_height',0.1)

wing_param.build_linear_airfoil(OC, AoA0=0., Cm0=-0.1, Ka=0.75, set_as_ref=True)
#wing_param.build_meta_airfoil(OC, '../ONERA_D.xml', relative_thickness=.09, camber=0, Sref=1., Lref=1., sweep=.0, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

print wing_param
    
def RunPolar(wing_param, OC, AoA_list, polar_name):
    
    polar_name =polar_name + '_%.2E' % OC.get_Mach()+'.dat'
    set_polar_file(OC,  polar_name)
    
    for i,alpha in enumerate(AoA_list):
        OC.set_AoA(alpha)
        output, F_list_names = runDLLM(wing_param, OC)
        tags=['AoA']+F_list_names
        add_to_polar(polar_name, numpy.hstack([[alpha],output]), tags, i)
        
def runDLLM(wing_param, OC):
    
        DLLM = DLLMSolver('SimpleRect',wing_param,OC)
        DLLM.run_direct()
        DLLM.run_post()
        
        output = DLLM.get_F_list()
        F_list_names = DLLM.get_F_list_names()
        return output, F_list_names      

def set_polar_file(OC, polar_name):
    
    with open(polar_name, 'w') as f:
	# Header is no up to date!!!
        header = '# polar for Mach ='+str(OC.get_Mach())
        f.write(header+'\n')
    f.closed
    
def add_to_polar(polar_name, values, tags, i):
    with open(polar_name, 'a') as f:
        if i == 0:
            line = "# "
            line+=string.join(tags,'\t')
            f.write(line+'\n')
        line  = ''
        for param in values:
            line = line+'%.10E    ' %  param 
        f.write(line+'\n')
    f.closed

for M in Mach_list:
    OC.set_Mach(M)
    OC.compute_atmosphere()
    RunPolar(wing_param,OC, AoA_list, 'SimpleRect')
