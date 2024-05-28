import sys
sys.path.append('/home/sh/WCY/T3/vali_Unsteady/validation')

from utils.config import VPM_ENV, VPM_Prop, VPM_Solver, VPM_Maneuver, VPM_Monitor, VPM_Simulation
from utils.process import initialize_VPMProcess

EXPName = 'Exp' 
EXP_stlName  = 'Exp_20inch.stl'
EXP_BEMName  = 'Exp_20inch.bem'

ModelPath_BEM = "0_model/BEM/"

## initialization
RPM                      = 5400
J                        = 0.01

ENV  = VPM_ENV()

BEM_Prop_Path = ModelPath_BEM + EXP_BEMName

Process = initialize_VPMProcess(BEM_Prop_Path, ENV, RPM, J)
## generation and calculation 
from utils.solver import *

SINGLE_VEHICLE_DEFINITION(Process)

SINGLE_MANEUVER_DEFINITION(Process)

