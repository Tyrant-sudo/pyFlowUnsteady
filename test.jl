import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm

Rtip = 0.15000000; Rhub = 0.03000000; B = 3; chorddist = [[0.2        0.08      ]
 [0.27272727 0.11601803]
 [0.34545455 0.14698723]
 [0.41818182 0.17182569]
 [0.49090909 0.18945154]
 [0.56363636 0.19878287]
 [0.63636364 0.2       ]
 [0.70909091 0.2       ]
 [0.78181818 0.2       ]
 [0.85454545 0.2       ]
 [0.92727273 0.2       ]
 [1.         0.13      ]]; twistdist = [[ 0.2        46.75      ]
 [ 0.27272727 42.3174305 ]
 [ 0.34545455 38.15773854]
 [ 0.41818182 34.27092412]
 [ 0.49090909 30.65698723]
 [ 0.56363636 27.31592787]
 [ 0.63636364 24.24774606]
 [ 0.70909091 21.45244177]
 [ 0.78181818 18.93001503]
 [ 0.85454545 16.68046582]
 [ 0.92727273 14.70379414]
 [ 1.         13.        ]]; sweepdist = [[0.2        0.        ]
 [0.27272727 0.        ]
 [0.34545455 0.        ]
 [0.41818182 0.        ]
 [0.49090909 0.        ]
 [0.56363636 0.        ]
 [0.63636364 0.        ]
 [0.70909091 0.        ]
 [0.78181818 0.        ]
 [0.85454545 0.        ]
 [0.92727273 0.        ]
 [1.         0.        ]]; heightdist = [[0.2        0.        ]
 [0.27272727 0.        ]
 [0.34545455 0.        ]
 [0.41818182 0.        ]
 [0.49090909 0.        ]
 [0.56363636 0.        ]
 [0.63636364 0.        ]
 [0.70909091 0.        ]
 [0.78181818 0.        ]
 [0.85454545 0.        ]
 [0.92727273 0.        ]
 [1.         0.        ]]; airfoil_contours = [(0.0, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.09090908750000003, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.18181818750000006, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.27272727500000005, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.3636363625, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.45454544999999996, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.5454545500000001, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.6363636375, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.727272725, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.8181818125, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (0.9090909125000001, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" ), (1.0, "/home/sh/WCY/T3/vali_Unsteady/validation/0_model/airfoil/NACA 4412.csv" , "polor_file" )];                        pitch = 0.00000000; CW = false; n = 10; blade_r = 0.20000000;                         spline_k = 5; spline_s = 0.00100000; spline_bc = "extrapolate";                         rediscretize_airfoils = true;rfl_n_lower = 15 ; rfl_n_upper = 15; data_path = "/home/wangchenyu/T3/FLOWUnsteady/database"; read_polar = vlm.ap.read_polar; xfoil = true;                         alphas = [-20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]; ncrit = 9;                        ReD = 1722233.52348066; Matip = 0.24776670; altReD = (5400.00000000, 0.01000000, 0.00001478);                        verbose = true; verbose_xfoil = false; v_lvl = 1;                         plot_disc = true; save_polars = "1_data/testSingleRotor_BEM_TEST"; figsize_factor=0.66666667

 rotor = uns.generate_rotor(
Rtip, 
Rhub, 
B,
chorddist, 
twistdist, 
sweepdist, 
heightdist,
airfoil_contours;

# MORE ROTOR PARAMETERS
pitch = pitch,
CW = CW,

# DISCRETIZATION SETTINGS
n       = n, 
blade_r = blade_r,
rediscretize_airfoils = rediscretize_airfoils,
rfl_n_lower = rfl_n_lower, 
rfl_n_upper = rfl_n_upper,
spline_k    = spline_k, 
spline_s    = spline_s, 
splines_s   = nothing, 
spline_bc   = spline_bc,

# AIRFOIL PROCESSING
data_path  = data_path,
read_polar = read_polar,
xfoil      = xfoil,
alphas     = alphas, 
ncrit      = ncrit,
ReD        = ReD, 
altReD     = altReD, 
Matip      = Matip,

# OUTPUT OPTIONS
verbose       = verbose, 
verbose_xfoil = verbose_xfoil, 
v_lvl         = v_lvl,
plot_disc     = plot_disc, 
save_polars   = save_polars
)

system = vlm.WingSystem()
vlm.addwing(system, "Rotor", rotor)

rotors = [rotor]
rotor_systems = (rotors, )

wake_system = vlm.WingSystem()


vlm.addwing(wake_system, "Rotor", rotor)

VehicleType = uns.UVLMVehicle

vehicle = VehicleType(  system;
                        rotor_systems = rotor_systems,
                        wake_system   = wake_system
                     )