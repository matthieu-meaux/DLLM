from DLLM.DLLMKernel.lifting_line_evaluator import Lifting_line_evaluator
from OpenDACE.evalcontrol.dv_handler import DVHandler
import numpy

#Mach=
#Alpha=s
data={}
data['Evaluation conditions.Mach']=0.5
data['Evaluation conditions.AoA']=3.0
data['Airfoils.Airfoil aerodynamic model']='Linear'
data['Airfoils.AoA0']=-2.0
data['Airfoils.Cd0']=0.01
data['Airfoils.Cm']=-0.1
data['Numerical parameters.Stop criteria']='Residual decrease'
data['Numerical parameters.Relaxation factor']=0.95
data['Numerical parameters.Decrease ratio']=1e-06
data['Wing geometry.Planform']='Broken'
data['Wing geometry.Wingspan']=34.1
data['Wing geometry.Break %']=33.
data['Wing geometry.Root chord']=6.1
data['Wing geometry.Break chord']=4.6
data['Wing geometry.Tip chord']=1.5
data['Wing geometry.Root thickness']=1.28
data['Wing geometry.Break thickness']=0.97
data['Wing geometry.Tip thickness']=0.33
data['Wing geometry.Number of sections']=40
data['design_variables']='twist.dat'

llw_evaluator = Lifting_line_evaluator(data)

alpha=data['Evaluation conditions.AoA']*numpy.pi/180.
Mach=data['Evaluation conditions.Mach']

DVH=DVHandler()
DVH.import_from_file(data['design_variables'])

llw_evaluator.eval(data,DVH,alpha,Mach,noise_level=0.)
wing=llw_evaluator.get_wing()
print 'Cl=',wing.Cl(alpha,beta=0.,Mach=Mach)
print 'Cd=',wing.Cd(alpha,beta=0.,Mach=Mach)
print 'Cm=',wing.Cm(alpha,beta=0.,Mach=Mach)

