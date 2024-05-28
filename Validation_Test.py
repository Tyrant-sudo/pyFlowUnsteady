import sys
sys.path.append('/home/sh/WCY/T3/vali_Unsteady/validation')

from utils.config import  VPM_ENV, VPM_Prop, VPM_Solver, VPM_Maneuver, VPM_Monitor, VPM_Simulation
from utils.process import VPM_Process
from utils.process import initialize_VPMProcess

ModelPath_BEM            = '0_model/BEM/'
ModelName                = 'TestBEM_Prop.bem'

## initialization
RPM                      = 5400
J                        = 0.01

ENV  = VPM_ENV()

Process  = VPM_Process(ENV)

BEMProp_path = ModelPath_BEM + ModelName

Process = initialize_VPMProcess(BEMProp_path, ENV, RPM, J)
## generation and calculation
 
from utils.solver import *

SINGLE_VEHICLE_DEFINITION(Process)

SINGLE_MANEUVER_DEFINITION(Process)

SINGLE_SIMULATION_DEFINITION(Process)

SINGLE_MONITORS_DEFINITIONS(Process)

# SINGLE_RUN_SIMULATION(Process)