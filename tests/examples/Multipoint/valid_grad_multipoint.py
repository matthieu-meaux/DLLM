from DLLM.DLLMEval.DLLMMP import DLLMMP
from OpenDACE.ValidGrad.FDValidGrad import FDValidGrad
import os
from glob import glob

config_dict={}
config_dict['Case.nb_conditions']=3
config_dict['Case.condition_name']='cond'
# cond1 Operating condition information
config_dict['Case.cond1.OC.Mach']=0.8
config_dict['Case.cond1.OC.AoA']=3.5
config_dict['Case.cond1.OC.altitude']=10000.
# cond2 Operating condition information
config_dict['Case.cond2.OC.Mach']=0.6
config_dict['Case.cond2.OC.AoA']=4.5
config_dict['Case.cond2.OC.altitude']=5000.
# cond3 Operating condition information
config_dict['Case.cond3.OC.Mach']=0.4
config_dict['Case.cond3.OC.AoA']=6.
config_dict['Case.cond3.OC.altitude']=1000.

# Common parameterisation configuration
config_dict['Case.param.geom_type']='Broken'
config_dict['Case.param.n_sect']=20
config_dict['Case.param.airfoil.type']='simple'
config_dict['Case.param.airfoil.AoA0']=-2.
config_dict['Case.param.airfoil.Cm0']=-0.1
# config_dict['Case.param.airfoil.type']='meta'
# config_dict['Case.param.airfoil.surrogate_model']='../MetaModelCleaning.xml'
config_dict['Case.param.desc.span.type']='DesignVariable'
config_dict['Case.param.desc.span.value']=34.1
config_dict['Case.param.desc.span.bounds']=(10.,50.)
config_dict['Case.param.desc.sweep.type']='DesignVariable'
config_dict['Case.param.desc.sweep.value']=34.
config_dict['Case.param.desc.sweep.bounds']=(0.,40.)
config_dict['Case.param.desc.break_percent.type']='Variable'
config_dict['Case.param.desc.break_percent.value']=33.
#config_dict['cond1.param.desc.break_percent.bounds']=(20.,40.)
config_dict['Case.param.desc.root_chord.type']='DesignVariable'
config_dict['Case.param.desc.root_chord.value']=6.1
config_dict['Case.param.desc.root_chord.bounds']=(5.,7.)
config_dict['Case.param.desc.break_chord.type']='DesignVariable'
config_dict['Case.param.desc.break_chord.value']=4.6
config_dict['Case.param.desc.break_chord.bounds']=(3.,5.)
config_dict['Case.param.desc.tip_chord.type']='DesignVariable'
config_dict['Case.param.desc.tip_chord.value']=1.5
config_dict['Case.param.desc.tip_chord.bounds']=(1.,2.)
config_dict['Case.param.desc.root_height.type']='DesignVariable'
config_dict['Case.param.desc.root_height.value']=1.28
config_dict['Case.param.desc.root_height.bounds']=(1.,1.5)
config_dict['Case.param.desc.break_height.type']='DesignVariable'
config_dict['Case.param.desc.break_height.value']=0.97
config_dict['Case.param.desc.break_height.bounds']=(0.8,1.2)
config_dict['Case.param.desc.tip_height.type']='DesignVariable'
config_dict['Case.param.desc.tip_height.value']=0.33
config_dict['Case.param.desc.tip_height.bounds']=(0.2,0.5)

config_dict['Case.DLLM.type']='Solver'
config_dict['Case.DLLM.method']='inhouse'
config_dict['Case.DLLM.relax_factor']=0.99
config_dict['Case.DLLM.stop_residual']=1e-9
config_dict['Case.DLLM.max_iterations']=100

list_log=glob('*.log')
for log in list_log:
    os.remove(log)

MP=DLLMMP('Case')
MP.configure(config_dict)
MP.set_out_format('numpy')
MP.set_grad_format('numpy')

x0=MP.get_x0()

val_grad=FDValidGrad(2,MP.run,MP.run_grad,fd_step=1.e-8)
ok,df_fd,df=val_grad.compare(x0,treshold=1.e-6,split_out=True,return_all=True)

# for j in xrange(len(df[:,0])):
#     fid=open('gradient_file'+str(j)+'.dat','w')
#     for i in xrange(len(x0)):
#         fid.write(str(i)+' '+str(df_fd[j,i])+' '+str(df[j,i])+'\n')
#     fid.close()

print '\n****************************************************'
if ok:
    print 'DLLMMP gradients are valid.'
else:
    print 'DLLMMP gradients are not valid!'
print '****************************************************'
