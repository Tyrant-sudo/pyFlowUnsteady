# transform stl to pointscloud
# from utils.config import VPM_process
from utils.visualization import plot_stl,easy_draw
from utils.geometry import build_PointsCloud_from_STL,build_BEM_from_PointsCloud
import os
import pandas as pd

ModelPath_STL            = '0_model/stl/'
ModelPath_PointCloud     = '0_model/PointsCloud/'
ModelPath_BEM            = '0_model/BEM/'
ModelPath_airfoil        = '0_model/airfoil/'


EXPName = 'Exp' 
EXP_stlName  = 'Exp_20inch.stl'
EXP_BEMName  = 'Exp_20inch.bem'
B        = 2

EXP_PointsCloud_Path = ModelPath_PointCloud + EXPName
PointsCloud_Direct    = [0,1,0]
Rtip_L, Rhub_L        = 1.0, 0.2
Num_Zsclices          = 20
Num_EachSlices        = 200

if not os.path.exists(EXP_PointsCloud_Path):
    PointsCloud_Result = build_PointsCloud_from_STL(ModelPath_STL + EXP_stlName,\
                        EXP_PointsCloud_Path ,\
                        PointsCloud_Direct    ,\
                        Rtip_L, Rhub_L, Num_Zsclices, Num_EachSlices)
    
# PointsCloud_Result = pd.read_csv(EXP_PointsCloud_Path)
# easy_draw(PointsCloud_Result.values,  ModelPath_PointCloud + EXPName + '.png',-45)

# 需要将 z 值划分的点云， 按顺序提取每一个截面的 airfoil[x,y], 
# Radius/R, Chord/R, Twist (deg), Rake/R, Skew/R, Sweep, t/c, CLi, Axial, Tangential 
EXP_BEM_Path  = ModelPath_BEM + EXP_BEMName
BEM_Result = build_PointsCloud_from_STL(EXP_PointsCloud_Path,\
                                             EXP_BEM_Path, 0.2,\
                                             True, 20, 2, 0.508 )




