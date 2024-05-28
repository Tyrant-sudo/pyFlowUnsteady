from .config import *

class VPM_Process:
    def __init__(self, VPM_ENV:VPM_ENV) -> None:

        self.state         = 'SingleRotor' #SingleRotor or CoaxialContra
        self.case_name     = 'TEST'
        self.file_type     = 'BEM'         #stl PointsCloud or BEM
        
        self.run_name      =  f"{self.state}_{self.file_type}_{self.case_name}"
        
        self.VehicleType   = "UVLMVehicle" # uns.UVLMVehicle UVLMVehicle for Unsteady solver, QLMVehicle for Quasi-steady solver
        
        self.VPM_ENV       = VPM_ENV
        self.rho           = VPM_ENV.rho
        self.mu            = VPM_ENV.mu
        self.speedofsound  = VPM_ENV.speedofsound
         
        self.data_path     = None     # the stored airfoils and propellers
        
        self.CW            = False    # Clock-wise rotation for the first 

        self.pitch         = 0.0      # (deg) collective pitch of blades
        
        self.AOA           = 0
        self.xfoil         = True     # Whether to run XFOIL  
        self.rand_RPM      = False    # (experimental) randomize RPM fluctuations


    def pre_single_rotor(self, prop:BEM_Prop, P:VPM_Prop, RPM , J ):

        R       = 1/2 * prop.Diameter 
        sec_num = prop.Num_Sections
        pi      = 3.14159

        P.Rtip = prop.Radius_R[-1] * R
        P.Rhub = prop.Radius_R[0] * R
        P.B    = prop.Num_Blade

        P.chorddist  = np.array([prop.Radius_R, prop.Chord_R]).T 
        P.twistdist  = np.array([prop.Radius_R, prop.Twist]).T
        P.sweepdist  = np.array([prop.Radius_R, prop.Sweep]).T
        P.heightdist = np.array([prop.Radius_R, prop.Rake_R]).T
        
        P.airfoil_contours = []

        for i in range(sec_num):
            # (r-Rhub)/(Rtip-Rhub)
            pos         = (prop.Radius_R[i] * R - P.Rhub)/(P.Rtip - P.Rhub)

            if prop.section_path[i]== "\"section_path\" ":
                contour     = prop.section_cod[i][::-1]
            else:
                contour     = prop.section_path[i]

            polorfile   = prop.polor_path[i]
            
            P.airfoil_contours.append((pos, contour, polorfile))
        
        string_contours   = ["(" + ", ".join(map(str, contour)) + ")" for contour in P.airfoil_contours]
        P.contours_string = "[" + ", ".join(string_contours) + "]"

        # P.n          = prop.Num_Sections

        if self.data_path != None:
            P.data_path = self.data_path
        
        P.pitch      = self.pitch
        P.CW         = self.CW
        P.xfoil      = self.xfoil
        omega        = RPM * 2 * pi / 60
        P.ReD        = (omega * R) * (2 * R)/ (self.mu/ self.rho)

        P.Matip      = 2*pi*RPM/60*R / self.speedofsound

        P.altReD     = (RPM, J, self.mu/ self.rho)
        
        P.save_polars = P.save_polars + self.run_name
        
        self.VPM_Prop       = P
        

        self.VPM_defProp = "Rtip = {:.8f}; Rhub = {:.8f}; B = {:d}; chorddist = {}; twistdist = {}; sweepdist = {}; heightdist = {}; airfoil_contours = {};\
                        pitch = {:.8f}; CW = {}; n = {:d}; blade_r = {:.8f}; \
                        spline_k = {:d}; spline_s = {:.8f}; spline_bc = \"{}\"; \
                        rediscretize_airfoils = {};rfl_n_lower = {:d} ; rfl_n_upper = {:d}; data_path = \"{}\"; read_polar = {}; xfoil = {}; \
                        alphas = {}; ncrit = {};\
                        ReD = {:.8f}; Matip = {:.8f}; altReD = ({:.8f}, {:.8f}, {:.8f});\
                        verbose = {}; verbose_xfoil = {}; v_lvl = {:d}; \
                        plot_disc = {}; save_polars = \"{}\"; figsize_factor={:.8f}"\
                        .format(P.Rtip, P.Rhub, P.B, P.chorddist, P.twistdist, P.sweepdist, P.heightdist, P.contours_string,\
                        P.pitch, julia_bool(P.CW), P.n, P.r,\
                        P.spline_k, P.spline_s, P.spline_bc, julia_bool(P.rediscretize_airfoils), P.rfl_n_lower, P.rfl_n_upper,\
                        P.data_path, P.read_polar, julia_bool(P.xfoil), \
                        P.alphas, P.ncrit,\
                        P.ReD, P.Matip, P.altReD[0], P.altReD[1], P.altReD[2],\
                        julia_bool(P.verbose), julia_bool(P.verbose_xfoil), P.v_lvl, julia_bool(P.plot_disc), P.save_polars, P.figsize_factor)
    
        ####### the simulation parameters
        self.R         = R
        self.RPM       = RPM
        self.J         = J
        self.magVinf   = J*RPM/60*(2*R)
        self.Vinf      = [self.magVinf * _ for _ in [np.cos(self.AOA), np.sin(self.AOA), 0]]
    
        self.VPM_defENV  = "rho = {:.8f}; mu = {:.8f}; speedofsound = {:.8f}; RPM = {:.8f}; J = {:.8f}".format(self.rho, self.mu, self.speedofsound, self.RPM, self.J)
    def pre_single_solver(self, S: VPM_Solver):
        
        S.nsteps          = S.nrevs*S.nsteps_per_rev # Number of time steps
        S.ttot            = S.nsteps/S.nsteps_per_rev / (self.RPM/60)       # (s) total simulation time
        S.max_particles   = ((2*self.VPM_Prop.n+1)*self.VPM_Prop.B)*S.nsteps*S.p_per_step + 1 # Maximum number of particles

        S.sigma_rotor_surf    = self.R / 40
        S.sigma_vpm_overwrite = S.lambda_vpm * 2*np.pi* self.R/(S.nsteps_per_rev* S.p_per_step)
        
        self.Solver       = S
        
        self.ttot = julia_none(S.ttot)
    
        if self.VPM_Prop.Matip != 0:
            sound_spd = "nothing"
        else:
            sound_spd = self.speedofsound
        

        if S.tquit > 1e8:
            tquit = "Inf"
        else:
            tquit = S.tquit

        S.run_name = "2_output/" + self.run_name
     
        SIMULATION_OPTIONS  = "Vinf = {}; sound_spd = {}; rho = {:.8f}; mu = {:.8f}; tquit = {} ; rand_RPM = {}; extra_runtime_function = monitors;"\
                      .format(self.Vinf,      sound_spd,      self.rho,     self.mu,     tquit , julia_bool(S.rand_RPM)  )

        # VPM Solver Options
        SOLVER_OPTIONS_vpm = "max_particles = {:d}; max_static_particles = {}; p_per_step = {:d}; vpm_formulation = {}; vpm_kernel = {}; vpm_UJ = {}; vpm_SFS = {}; vpm_integration = {}; vpm_transposed = {}; vpm_viscous = {}; vpm_fmm = {}; vpm_relaxation = {}; vpm_surface = {};"\
                    .format(S.max_particles       , julia_none(S.max_static_particles), S.p_per_step, S.vpm_formulation, S.vpm_kernel,   S.vpm_UJ,    S.vpm_SFS,    S.vpm_integration,     julia_bool(S.vpm_transposed), S.vpm_viscous, S.vpm_fmm, S.vpm_relaxation, julia_bool(S.vpm_surface) )
        
        # ASM/ALM Solver Options
  
        SOLVER_OPTIONS_asm = "vlm_vortexsheet = {}; vlm_vortexsheet_overlap = {:.8f}; vlm_vortexsheet_distribution = {}; vlm_vortexsheet_sigma_tbv = {}; vlm_rlx = {:.8f}; vlm_init = {}; hubtiploss_correction = {};"\
                    .format(julia_bool(S.vlm_vortexsheet), S.vlm_vortexsheet_overlap, S.vlm_vortexsheet_distribution,    julia_none(S.vlm_vortexsheet_sigma_tbv), S.vlm_rlx, julia_bool(S.vlm_init), S.hubtiploss_correction)

        # Wake Shedding Options
        SOLVER_OPTIONS_wake = "wake_coupled = {}; shed_unsteady = {}; unsteady_shedcrit = {:.8f}; shed_starting = {}; shed_boundarylayer = {}; boundarylayer_prescribedCd = {:.8f}; boundarylayer_d = {:.8f}; omit_shedding = {};"\
                       .format(julia_bool(S.wake_coupled), julia_bool(S.shed_unsteady), S.unsteady_shedcrit, julia_bool(S.shed_starting), julia_bool(S.shed_boundarylayer), S.boundarylayer_prescribedCd, S.boundarylayer_d, S.omit_shedding)
        
        # Regularization Options
        SOLVER_OPTIONS_regu = "sigma_vlm_solver = {:.8f}; sigma_vlm_surf = {:.8f}; sigma_rotor_surf = {:.8f}; sigmafactor_vpm = {:.8f}; sigmafactor_vpmonvlm = {:.8f}; sigma_vpm_overwrite = {};"\
                       .format(S.sigma_vlm_solver,        S.sigma_vlm_surf,        S.sigma_rotor_surf,        S.sigmafactor_vpm,        S.sigmafactor_vpmonvlm,        julia_none(S.sigma_vpm_overwrite))
        
        # Output Options
        SOLVER_OPTIONS_out = "restart_vpmfile = {};save_path = {}; run_name = {}; create_savepath = {}; prompt = {}; verbose = {}; v_lvl = {:.8f}; verbose_nsteps = {:d}; raisewarnings = {}; debug = {}; nsteps_save = {:d}; nsteps_restart = {:d}; save_code = {}; save_horseshoes = {}; save_static_particles = {}; save_wopwopin = {};"\
                       .format(julia_none(S.restart_vpmfile), julia_none(S.save_path), julia_none(S.run_name), julia_bool(S.create_savepath), julia_bool(S.prompt), julia_bool(S.verbose), S.v_lvl, S.verbose_nsteps, julia_bool(S.raisewarnings), julia_bool(S.debug), S.nsteps_save, S.nsteps_restart, S.save_code, julia_bool(S.save_horseshoes), julia_bool(S.save_static_particles), julia_bool(S.save_wopwopin) )
        
        self.VPM_defSolver     = SIMULATION_OPTIONS +  SOLVER_OPTIONS_vpm + SOLVER_OPTIONS_asm + SOLVER_OPTIONS_wake + SOLVER_OPTIONS_regu + SOLVER_OPTIONS_out
        

    def pre_single_maneuver(self, MAN: VPM_Maneuver):

        self.MAN        = MAN
        self.VPM_defMan = MAN.static_Maneuver() 
    
    def pre_single_simulation(self, SIM: VPM_Simulation):
        
        MAN        = self.MAN
        S          = self.Solver

        SIM.RPMref = self.RPM
        
        SIM.Vinit  = SIM.Vref * MAN.Vvehicle(0) 
        SIM.Winit  = np.pi/ 180 * ( MAN.anglevehicle(1e-6) - MAN.anglevehicle(0) )/ (1e-6 * S.ttot)

        self.SIM        = SIM
        self.VPM_defSIM = "Vref = {:f}; RPMref = {:f}; Vinit = {}; Winit = {}".format(\
        SIM.Vref, SIM.RPMref, SIM.Vinit, SIM.Winit\
        )

    def pre_single_monitor(self, MON: VPM_Monitor):
        
        MON.save_path = MON.save_path + self.run_name
        self.MON   = MON
        all_mon    = [ MON.set_monitor_rotor , MON.set_monitor_enstrophy , MON.set_monitor_Cd ]
        
        self.VPM_defMON_rotor = \
        """
        save_path = \"{}\";
        monitor_rotor = uns.generate_monitor_rotors(rotors, J, rho, RPM, nsteps;
                                            t_scale=RPM/60,       
                                            t_lbl="Revolutions",  
                                            save_path=save_path,
                                            run_name=run_name,
                                            figname="rotor monitor",
                                            )\n
        """.format(MON.save_path)
        
        self.VPM_defMON_enstrophy = \
        """
        save_path = \"{}\";
        monitor_enstrophy = uns.generate_monitor_enstrophy(;
                                                    save_path=save_path,
                                                    run_name=run_name,
                                                    figname="enstrophy monitor"
                                                    )\n
        """.format(MON.save_path)

        self.VPM_defMON_Cd = \
        """
        save_path = \"{}\";
        monitor_Cd = uns.generate_monitor_Cd(;
                                                    save_path=save_path,
                                                    run_name=run_name,
                                                    figname="Cd monitor"
                                                    )\n
        """.format(MON.save_path)

        result   = ""
        monitors = ""
        if MON.set_monitor_rotor:
            result   += self.VPM_defMON_rotor
            monitors += "monitor_rotor,"
        elif MON.set_monitor_enstrophy:
            result   += self.VPM_defMON_enstrophy
            monitors += "monitor_enstropy,"
        elif MON.set_monitor_Cd:
            result   += self.VPM_defMON_Cd
            monitors += "monitor_Cd"
        
        self.monitors = monitors

        self.VPM_defMON = result

    def pre_coaxial_rotor(prop1:VPM_Prop, prop2:VPM_Prop):
        1


class FWH_Process:

    software_path  = '/home/wangchenyu/T3/vali_Unsteady/validation/utils/WopWop3/wopwop3_linux'  



def initialize_VPMProcess(BEM_Prop_Path, ENV, RPM, J):

    Process  = VPM_Process(ENV)

    BEMProp  = build_BEMProp_from_BEMfile(BEM_Prop_Path)
    Prop     = VPM_Prop()
    Process.pre_single_rotor(BEMProp, Prop, RPM, J)

    Solver   = VPM_Solver()
    Process.pre_single_solver(Solver)

    Maneuver = VPM_Maneuver()
    Process.pre_single_maneuver(Maneuver)

    Simulation = VPM_Simulation()
    Process.pre_single_simulation(Simulation)

    Monitor    = VPM_Monitor()
    Process.pre_single_monitor(Monitor)

    print("Initilization finished")

    return Process

def build_BEMProp_from_BEMfile(filename, SectionPath_List = [], PolarPath_List = []):

    Class_BEM_Prop = BEM_Prop()

    with open(filename, 'r') as file:
        lines = file.readlines()

        # Extract general properties
        Class_BEM_Prop.Num_Sections = int(lines[1].split(': ')[1])
        Class_BEM_Prop.Num_Blade = int(lines[2].split(': ')[1])
        Class_BEM_Prop.Diameter = float(lines[3].split(': ')[1])/ 100 
        Class_BEM_Prop.Beta = float(lines[4].split(': ')[1])
        Class_BEM_Prop.Feather = float(lines[5].split(': ')[1])
        Class_BEM_Prop.Pre_Cone = float(lines[6].split(': ')[1])
        Class_BEM_Prop.Center = [float(x) for x in lines[7].split(': ')[1].split(',')]
        Class_BEM_Prop.Normal = [float(x) for x in lines[8].split(': ')[1].split(',')]
        
        # Extract properties for each section
        properties_start = 11
        properties_end = properties_start + Class_BEM_Prop.Num_Sections
        properties = [list(map(float, line.split(','))) for line in lines[properties_start:properties_end]]
        properties = np.array(properties).T

        properties_names = \
        ['Radius_R', 'Chord_R', 'Twist' , 'Rake_R', 'Skew_R', 'Sweep', 't_c', 'CLi', 'Axial', 'Tangential' ]
        
        for i, property in enumerate(properties):
            setattr(Class_BEM_Prop, properties_names[i], property)
        
        # Extract coordinates for each section
        Class_BEM_Prop.section_cod = []
        coord_start = properties_end + 2
        current_section = []
        i = coord_start
        while i < len(lines):
            line = lines[i]
            if line.strip() == '':  # 检查是否为空行或是否是文件的最后一行
                if current_section:
                    Class_BEM_Prop.section_cod.append(np.array(current_section))
                    current_section = []

                    if len(Class_BEM_Prop.section_cod) == Class_BEM_Prop.Num_Sections:
                        break

                if line.strip() == '':
                    i += 2  # 如果是空行，跳过接下来的两行
                else:
                    i += 1  # 如果是文件的最后一行，结束循环
            else:
                coords = list(map(float, line.split(',')))
                current_section.append(coords)
            
                if i == len(lines) - 1:
                    Class_BEM_Prop.section_cod.append(np.array(current_section))
                i += 1

        # Class_BEM_Prop.section_cod = np.array(Class_BEM_Prop.section_cod)
        if SectionPath_List:
            Class_BEM_Prop.section_path = SectionPath_List
        else:
            Class_BEM_Prop.section_path   = ["\"section_path\" " for _ in range(Class_BEM_Prop.Num_Sections)]  
        
        if PolarPath_List:
            Class_BEM_Prop.polar_path = PolarPath_List
        else:
            Class_BEM_Prop.polor_path   = ["\"polor_file\" " for _ in range(Class_BEM_Prop.Num_Sections)]  

    return Class_BEM_Prop

def build_BEMProp_from_STLfile(STLfile, output_file,  PointsCloud_Direct, Rtip_L, Rhub_L ,root_rR,  Num_Zsclices = 14, Num_EachSlices = 200, direct = True, Beta_34 = 20, Num_Blade = 2, Diameter = 0.5, save_PointsCloudfile = True ):
    """
    PointsCloud_Direct: the direct from hub to tip, which will be redirected to [0,0,1]
    Rtip_L, Rhub_L: the position of Rtip and Rhub at the lenth of stl model

    root_rR:  

    """
    from .geometry import build_PointsCloud_from_STL 
    
    

def create_geometry_name(g:VPM_Process):
    
    name = f'Geo_{g.state}_{g.case}_{g.Drotor}'
    return name
    
def julia_bool(B:bool):

    if B:
        return "true"
    else:
        return "false"

def julia_none(N):
    if N == None:
        return "nothing"
    else:
        return N