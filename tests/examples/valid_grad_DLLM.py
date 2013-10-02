from OpenDACE.ValidGrad.FDValidGrad import FDValidGrad
from DLLM.DLLMGeom.wing_param import Wing_param
from DLLM.DLLMKernel.DLLMSolver import DLLMSolver
from MDOTools.Base.operating_condition import OperatingCondition

wing_param=Wing_param('test_param',geom_type='Broken',n_sect=20)
wing_param.build_wing()
wing_param.set_value('test_param.span',34.1)
wing_param.set_value('test_param.break_percent',33.)
wing_param.set_value('test_param.root_chord',6.1)
wing_param.set_value('test_param.break_chord',4.6)
wing_param.set_value('test_param.tip_chord',1.5)
wing_param.set_value('test_param.root_height',1.28)
wing_param.set_value('test_param.break_height',0.97)
wing_param.set_value('test_param.tip_height',0.33)
wing_param.convert_to_design_variable('test_param.span',10.,50.)
wing_param.convert_to_design_variable('test_param.break_percent',20.,40.)
wing_param.convert_to_design_variable('test_param.root_chord',5.,7.)
wing_param.convert_to_design_variable('test_param.break_chord',3.,5.)
wing_param.convert_to_design_variable('test_param.tip_chord',1.,2.)
wing_param.convert_to_design_variable('test_param.root_height',1.,1.5)
wing_param.convert_to_design_variable('test_param.break_height',0.8,1.2)
wing_param.convert_to_design_variable('test_param.tip_height',0.2,0.5)
wing_param.build_linear_airfoil(AoA0=-2., Cd0=0.01, Cm0=-0.1, set_as_ref=True)
wing_param.build_airfoils_from_ref()
wing_param.update()

print wing_param

OC=OperatingCondition('cond1',atmospheric_model='simple')
OC.set_Mach(0.5)
OC.set_AoA(3.0)
OC.set_altitude(5000.)
OC.set_T0_deg(20.)
OC.set_P0(101325.)
OC.set_humidity(0.)
OC.compute_atmosphere()

#print OC
x0=wing_param.get_dv_array()
print 'dv array shape',x0.shape
print 'dv_array=',x0

def f(x):
    wing_param.update_from_x_list(x)
    DLLM = DLLMSolver(wing_param,OC)
    DLLM.run_direct()
    #DLLM.run_post(func_list=['Cl'])
    DLLM.run_post()
    func=DLLM.get_func_values()
    return func

def df(x):
    wing_param.update_from_x_list(x)
    DLLM = DLLMSolver(wing_param,OC)
    DLLM.run_direct()
    #DLLM.run_post(func_list=['Cl'])
    DLLM.run_post()
    DLLM.run_adjoint()
    func_grad=DLLM.get_dJ_dchi()
    return func_grad

val_grad=FDValidGrad(2,f,df,fd_step=1.e-3)
ok,df_fd,df=val_grad.compare(x0,treshold=1.e-2,return_all=True)

for j in xrange(len(df[:,0])):
    fid=open('gradient_file'+str(j)+'.dat','w')
    for i in xrange(len(x0)):
        fid.write(str(i)+' '+str(df_fd[j,i])+' '+str(df[j,i])+'\n')
    fid.close()

print '\n****************************************************'
if ok:
    print 'Gradients are valid.'
else:
    print 'Gradients are not valid!'
print '****************************************************'
