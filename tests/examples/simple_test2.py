# Imports
from DLLM.DLLMGeom.wing_manufacture import Wing_manufacture
from OpenDACE.base.design_space import DesignSpace
from DLLM.DLLMEval.DLLM_evaluator import DLLMEvaluator
from MDOTools.Base.operating_condition import OperatingCondition

n_sect=20

# Build geometry
wmf=Wing_manufacture()
wmf.set_wing_geometry_type('Broken')
wmf.set_wingspan(34.1)
wmf.set_root_chord(6.1)
wmf.set_root_height(1.28)
wmf.set_break_percent(33.)
wmf.set_break_chord(4.6)
wmf.set_break_height(0.97)
wmf.set_tip_chord(1.5)
wmf.set_tip_height(0.33)

wing_geometry=wmf.build_Discrete_wing(n_sect,twist_law=None)

# Define design parameters
DS=DesignSpace()
for i in xrange(n_sect):
    tag='twist'+str(i)
    DS.add_DV(tag,0.,(-10.,10.))
print DS

# Update geometry with design parameters
# No geometry update since it's only the twist law... why? 

# Why not put all design parameters in the design space?
# Process should be (geometry+update)->(solver set up + run)->results

# Build operating condition object
OC=OperatingCondition('cond1',Mach=0.5,AoA=3.0)
print OC

# Set up DLLM evaluator
DLLMEval = DLLMEvaluator(wing_geometry,OC)
DLLMEval.build_linear_airfoil(AoA0=-2., Cd0=0.01, Cm0=-0.1, set_as_ref=True)
DLLMEval.build_airfoils_from_ref()
DLLMEval.setup()

# Set DLLM numerical parameters
DLLMEval.set_relax_factor(0.99)
DLLMEval.set_stop_criteria(residual=1.e-6)
# 
# # Run DLLM evaluator
DLLMEval.run_direct()

# Get outputs (Lift, Drag, loads,...)