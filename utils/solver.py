from .process import VPM_Process

from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import PyCall

from julia import Main
# 导入 Julia 模块

Main.eval("import FLOWUnsteady as uns")    
Main.eval("import FLOWVLM as vlm")    
Main.eval("import FLOWVPM as vpm")


def  SINGLE_VEHICLE_DEFINITION(V:VPM_Process):
    print("Generating vehicle...")
    
    VPM_ENV = V.VPM_ENV
    Main.eval(V.VPM_defENV)

    VPM_Prop = V.VPM_Prop
    
    Main.eval(V.VPM_defProp)
    
    generate_rotor = """
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
    """

    Main.eval(generate_rotor)

    generate_vehicle = """
    system = vlm.WingSystem()
    vlm.addwing(system, "Rotor", rotor)

    rotors = [rotor]
    rotor_systems = (rotors, )

    wake_system = vlm.WingSystem()
    
    
    {}

    {}

    vehicle = VehicleType(  system;
                            rotor_systems = rotor_systems,
                            wake_system   = wake_system
                         )
    """.format(  "vlm.addwing(wake_system, \"Rotor\", rotor)" if V.VehicleType == "UVLMVehicle" else "", "VehicleType = uns.UVLMVehicle" if V.VehicleType == "UVLMVehicle" else "VehicleType = uns.QVLMVehicle")

    Main.eval(generate_vehicle)
    return 

def SINGLE_MANEUVER_DEFINITION(V:VPM_Process):
    print("Generating maneuver...")

    VPM_MAN  = V.MAN
    Main.eval(V.VPM_defMan)
    
    generate_maneuver = "maneuver = uns.KinematicManeuver(angles, RPMs, Vvehicle, anglevehicle)"
    Main.eval(generate_maneuver)

    return 

def SINGLE_SIMULATION_DEFINITION(V:VPM_Process):
    print("Generating simulation...")

    SIM      = V.SIM
    Solver   = V.Solver

    Main.eval(V.VPM_defSIM)
    
    generate_simulation = "simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, {:.8f}; Vinit = Vinit, Winit = Winit)".format(V.ttot)
    Main.eval(generate_simulation)

    return 

def SINGLE_MONITORS_DEFINITIONS(V:VPM_Process):
    print("Generating monitors...")
    
    Main.eval(V.VPM_defMON)

    generate_monitors = "monitors = uns.concatenate({})".format(V.monitors)
    Main.eval(generate_monitors)


    return 

def SINGLE_RUN_SIMULATION(V:VPM_Process):
    print("Running simulation...")
    
    Solver = V.Solver
    Main.eval(V.VPM_defSolver)

    run_simulation = """uns.run_simulation( 
            sim,                  
            nsteps;               

            Vinf            = V_inf,
            sound_spd       = sound_spd,
            rho             = rho,
            mu              = mu,
            tquit           = tquit,
            rand_RPM        = rand_RPM,
            extra_runtime_function = extra_runtime_function,
            max_particles   = max_particles,
            max_static_particles = max_static_particles
            p_per_step      = p_per_step,

            vpm_formulation = vpm_formulation,
            vpm_kernel      = vpm_kernel,
            vpm_UJ          = vpm_UJ,
            vpm_SFS         = vpm_SFS,
            vpm_integration = vpm_integration,
            vpm_transposed  = vpm_transposed,
            vpm_viscous     = vpm_viscous,
            vpm_fmm         = vpm_fmm,
            vpm_relaxation  = vpm_relaxation,
            vpm_surface     = vpm_surface,

            vlm_vortexsheet = vlm_vortexsheet,  
            vlm_vortexsheet_overlap     = vlm_vortexsheet_overlap,  
            vlm_vortexsheet_sigma_tbv   = vlm_vortexsheet_sigma_tbv,  
            vlm_rlx         = vlm_rlx, 
            vlm_init        = vlm_init,  
            hubtiploss_correction = hubtiploss_correction, 

            wake_coupled        = wake_coupled,
            shed_unsteady       = shed_unsteady,
            unsteady_shedcrit   = unsteady_shedcrit,  
            shed_starting       = shed_starting, 
            shed_boundarylayer  = shed_boundarylayer,
            boundarylayer_prescribedCd = boundarylayer_prescribedCd,   
            boundarylayer_d     = boundarylayer_d,
            omit_shedding       = omit_shedding, 

            sigma_vlm_solver    = sigma_vlm_solver,  
            sigma_vlm_surf      = sigma_vlm_surf,   
            sigma_rotor_surf    = sigma_rotor_surf,   
            sigmafactor_vpm     = sigmafactor_vpm,         
            sigmafactor_vpmonvlm = sigmafactor_vpmonvlm, 
            sigma_vpm_overwrite = sigma_vpm_overwrite,  

         
            restart_vpmfile = restart_vpmfile,   
            save_path       = save_path,          
            run_name        = run_name,
            create_savepath = create_savepath,             
            prompt          = prompt,             
            verbose         = verbose,             
            v_lvl           = v_lvl,                
            verbose_nsteps  = verbose_nsteps,               
            raisewarnings   = raisewarnings,             
            debug           = debug,            
            nsteps_save     = nsteps_save,                
            nsteps_restart  = nsteps_restart,               
            save_code       = save_code
            save_horseshoes = save_horseshoes,            
            save_static_particles = save_static_particles,       
            save_wopwopin   = save_wopwopin,            

        )
    """

    Main.eval(run_simulation)