# Imports
from DLLM.DLLMGeom.wing_manufacture import Wing_manufacture
from DLLM.DLLMGeom.wing_geom_handler import Wing_geom_handler

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

wing_geometry=wmf.build_Discrete_wing(40,twist_law=None)

# Define design parameters

# Update geometry with design parameters

# Build operating conditions object

# Run DLLM evaluator

# Get outputs (Lift, Drag, loads,...)