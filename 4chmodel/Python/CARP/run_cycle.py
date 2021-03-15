#!/usr/bin/env python
"""
Electromechanics benchmark for 4 chamber heart
"""

# --- import packages ---------------------------------------------------------
from __future__ import print_function
from os.path import join, isfile, dirname, realpath
from shutil import copy
from datetime import date

from carputils import settings
from carputils import tools
from carputils import model
from carputils import ep
from carputils.resources import petsc_options
from carputils.resources import petsc_block_options

# --- command line options ----------------------------------------------------
def parser_commands():
    """ Define parser commands """
    parser = tools.standard_parser()

    parser.add_argument('--t_peak',
                        type=float, default=100.0, # Ad hoc to achieve a higer EF. Marina if 100, but check (before it was higher) https://watermark.silverchair.com/cvq318.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAkwwggJIBgk qhkiG9w0BBwagggI5MIICNQIBADCCAi4GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM0GEFCgr3TmlYZomsAgEQgIIB_2zsbaVjOJgaYyP-Nr9PHVnX2fteJ9TRKx5UUYFj0oMvRzYqTY4CroW9Wlsx2afbwx7BHWc92b4hvKX4T3CxGix1BbQlS768x9K2WFhUX8_mAjlPskOaR-13T3AnAxDAZJNwp4-f0OkUuJpD3Lyw5pVFPgXBPMF0_neKwCVWJNULgohnqcyvIisII15XwxjkB4MXmhKX6MjvqiCeJ2nACtb7UKk6eqCorz6Qcowb7JGJrb4y8EZBguElgD8rBMjTuy8VoaGbmVFY7Y9aq-_OsUhIzcVQN6tEOceTnB2GJQ5rPAHd063PcwD6xPBRtAOt48CeEf-VVuYzJTzmxgumJQIwrc3uSwnQwD6U8kUwscuSruHfDmvqmbEO9YP82TKAQq4ti9GEr5_Wa4DzNvtVbmKxmNU5wL7YffLJHqdswjE87wFql_vdew4vZwKjKHlNupU3tHc7z0UIYXOzEqzpjYBzEcq4L5ZHuuFgodvztqAOJM8xg3VSgC2aGrCkc5cl9K7-z08-UxfcI6qMi6HG0M-U7pU0bWEOFeEEb6Itowlb97ByRoR_5FZZCY_0mrcd_VSl9DwcmiUxSxrqAQvIQ51JHAkEjeEU4SD3RrianCzF2G5-LHqKILr_aMG4AbaYHSiHh6Nrsuc_glyMNoKBjm6Q3MO0WmItDvoIxeOxTOs
                        help='peak value for active stress')

    parser.add_argument('--aortic_press',
			type=float,
			default=77,
            help='Diastolic Blood Pressure. Default extracted from https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0182173.t003 DBP age 50.')

    parser.add_argument('--case',
                        default='01HC',
                        help='Pick mesh')

    parser.add_argument('--duration',
                        type=float, default=-1,
                        help='duration of experiment in ms \
                              (default: -1: auto setting)')
    parser.add_argument('--experiment',
                        default='wk3',
                        choices=['wk3','unloading','free_contraction'])
    parser.add_argument('--loadStepping',
                        type=float, default=50.0)
    parser.add_argument('--mechDT',
                        type=float,
                        default=1.)
    parser.add_argument('--periOff',
                        action='store_true')
    parser.add_argument('--rvEDP_kPa',
                        type=float,
                        default=0.8)
    parser.add_argument('--lvEDP_kPA',
                        type=float,
                        default=1.6) #Relationships between beat-to-beat interval and the strength of contraction in the healthy and diseased human heart
    parser.add_argument('--scaling_Guccione',
                        type=float,
                        default=1.7)
    parser.add_argument('--sim_name',
			            default='')
    parser.add_argument('--stiffness_pericardium',
                        type=float,
                        default=0.001)
    parser.add_argument('--stiffness_LSPV',
                        type=float,
                        default=0.001)
    parser.add_argument('--stiffness_RSPV',
                        type=float,
                        default=0.001)
    parser.add_argument('--stiffness_SVC',
                        type=float,
                        default=0.001)
    parser.add_argument('--Tanh_time_contract',
                        type=float,
                        default=50.0,
                        help='Time constant of contraction for the Tanh model. If length-dependance is off, corresponds to the baseline time constant.')
    parser.add_argument('--Tanh_time_relax',
                        type=float,
                        default=50.0,
                        help='Time constant of relaxation for the Tanh model. Corresponds to tau_r in https://carpentry.medunigraz.at/carputils/examples/tutorials/tutorials.03_EM_single_cell.01_EM_coupling.run.html')
    return parser

def job_id(args):
    """
    Generate name of top level output directory.
    """
    if args.sim_name == '':
        jobname= '{}_{}_tpeak_{}'.format(args.case,args.experiment,int(args.t_peak))
    else:
        jobname = args.sim_name

    return jobname

def run(args, job):
    """ MAIN"""
    # Initial command
    cmd = tools.carp_cmd()
    # add mesh and simulation ID
    mesh_dir, mesh_name = mesh_file(args.case)
    cmd += ['-simID', job.ID,
            '-meshname', mesh_name]

    # --- configuration variables ---------------------------------------------
    cmd += setup_time_variables(args)

    # --- create eikonal activation -------------------------------------------
    #if args.experiment == 'wk3':
    #    cmd += set_propagation_opts(args, cmd, job)

    # --- boundary settings ---------------------------------------------------
    cmd += setup_bc(mesh_dir,args)

    # --- cardiovascular system settings --------------------------------------
    cmd += set_cv_sys(args,mesh_dir)

    # --- mechanical settings -------------------------------------------------
    # define materials for different regions
    cmd += setup_material(args)

    # set solver properties
    cmd += set_mechanic_options(args)

    # --- electrical settings ------------------------------------------------
    #if args.experiment == 'wk3':
        #cmd += setup_electrical_parameters()     # setup resolution dependent electrical parameters
    cmd += setup_stimuli(args.case)     # define stimuli
    cmd += setup_active(args)     # active stress setting
        #cmd += setup_em_coupling()     # electro mechanical coupling

    # visualization
    cmd += setup_visualization(args.visualize, args.postprocess)

    if args.experiment == 'unloading':
        cmd += setup_unloading(args)

    # Run simulation
    job.carp(cmd)

# =============================================================================
#    FUNCTIONS
# =============================================================================

def mesh_file(case):
    """ Get the full path of a file in the mesh directory
    """
    mesh_dir = join(settings.config.MESH_DIR, 'h4ckcl/meshes', case)
    mesh_name = join(mesh_dir, 'heart_case{}_800um'.format(case))
    return mesh_dir, mesh_name

# ----------------------------------------------------------------------------
def setup_stimuli(case):
    """ Set stimulus
    """
    act_seq_path = join(settings.config.MESH_DIR, 'h4ckcl/actSeq')
    act_seq_file = join(act_seq_path, '{}_vm_act_seq.dat'.format(case))
    # general stimulus options
    stm_opts = ['-num_stim', 1,
                '-stimulus[0].stimtype', 8,
                '-stimulus[0].data_file', act_seq_file,
		'-diffusionOn', 0]

    return stm_opts

# ----------------------------------------------------------------------------
"""
def set_propagation_opts(args,cmd,job):
    
   # pick propagation mode
    
    prop_opts = []

    # Eikonal propagation


    # if we don't get an explicit activation sequence, look for it in the base dir
    act_seq = args.act_seq
    if act_seq == '':
        # try to find activation sequence vector from previous eikonal run
        if isfile(act_seq_file):
            act_seq = act_seq_file

    if act_seq == '' and args.experiment != 'unloading':
        # no activation sequence given or found, run eikonal first to generate one
        print('No activation sequence given, run Eikonal first')

        ek_sim_id = join(job_id(args), 'eikonal')
        

        prop_opts += ['-experiment', 6,
                      '-localize_pts', 1,
                      '-pstrat', 1,
                      '-pstrat_i', 1,
                      '-mechDT', 0]
        ek_cmd = cmd + prop_opts + ['-simID', ek_sim_id]

        # run eikonal first
        job.carp(ek_cmd)

        # copy activation sequence vector to base dir including a mesh flag
        # pick the activation sequence without PS component as PS
        # is not supported yet in forced foot mode
        act_seq = './vm_act_seq_%s.dat' % (mesh_flag(args.mesh))
        copy('%s/vm_act_seq.dat'%(ek_sim_id), act_seq)

    # if we don't get an explicit activation sequence, look for it in the base dir
    stimulus = ep.fake_ep_opts(0, act_seq, pulse_file='nfoot')

    # in this case we turn off diffusion while propagation is ongoing
    diff = ['-diffusionOn', 1,
            '-num_diffusionOffSections', 1,
            '-diffusionOff[0].start', 0.0,
            '-diffusionOff[0].stop', 15.0]

    # forced foot allows larger time stepping
    d_t = 100.

    prop_opts = stimulus + diff + ['-dt', d_t] + ['-simID', job.ID]

    return prop_opts
"""
# ----------------------------------------------------------------------------
def setup_electrical_parameters():
    """ setup electrical parameters"""

    ep_opts = ['-bidomain', 0,
               '-diffusionOn', 0]

    return ep_opts

# ----------------------------------------------------------------------------
def set_cv_sys(args,mesh_dir):
    """ setup cardiovasculare system """

    mech_grid = model.mechanics.grid_id()
    vols = ['-numSurfVols', 4,
            '-surfVols[0].name', 'LV_endo',
            '-surfVols[0].surf_file', join(mesh_dir, 'LV_endo'),
            '-surfVols[0].grid', mech_grid,
            '-surfVols[1].name', 'RV_endo',
            '-surfVols[1].surf_file', join(mesh_dir, 'RV_endo'),
            '-surfVols[1].grid', mech_grid,
            '-surfVols[2].name', 'LA_endo',
            '-surfVols[2].surf_file', join(mesh_dir, 'LA_endo'),
            '-surfVols[2].grid', mech_grid,
            '-surfVols[3].name', 'RA_endo',
            '-surfVols[3].surf_file', join(mesh_dir, 'RA_endo'),
            '-surfVols[3].grid', mech_grid]

    if args.experiment == 'wk3':
        vols += ['-CV_coupling', 0,
            '-CV_FE_coupling', 1,
            '-CVS_mode', 0,
            '-num_cavities', 2,
            '-cavities[0].cav_type', 0,
            '-cavities[0].cavP', 6,
            '-cavities[0].tube', 2,
            '-cavities[0].cavVol', 0,
            '-cavities[0].p0_cav', 7.5*float(args.lvEDP_kPA), # mmHg, final pressure, 
            '-cavities[0].p0_in', 7.5*float(args.lvEDP_kPA), # mmHg, final pressure
            '-cavities[0].p0_out', float(args.aortic_press), 
            '-cavities[0].valve', 0,
            '-cavities[0].state', -1,
            '-lv_wk3.name', 'Aorta',
            '-lv_wk3.R1', 0.03,# https://link.springer.com/content/pdf/10.1007%2Fs10237-013-0523-y.pdf
            '-lv_wk3.R2', 0.63,# https://link.springer.com/content/pdf/10.1007%2Fs10237-013-0523-y.pdf
            '-lv_wk3.C', 5.16,# https://link.springer.com/content/pdf/10.1007%2Fs10237-013-0523-y.pdf
            '-mitral_valve.Rfwd', 0.05, # ad hoc
            '-mitral_valve.Rbwd', 1000.0, # ad hoc
            '-aortic_valve.Rfwd', 0.0, # ad hoc
            '-cavities[1].cav_type', 1,
            '-cavities[1].cavP', 7,
            '-cavities[1].tube', 2,
            '-cavities[1].cavVol', 1,
            '-cavities[1].p0_cav', 7.5*float(args.rvEDP_kPa), # mmHg https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4449250/pdf/PulmCirc-005-370.pdf
            '-cavities[1].p0_in', 7.5*float(args.rvEDP_kPa), #mmHg https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4449250/pdf/PulmCirc-005-370.pdf
            '-cavities[1].p0_out', 17.4, #Right ventricle dimensions and function in response to acute hypoxia in healthy human subjects
            '-cavities[1].valve', 0,
            '-cavities[1].state', -1,
            '-rv_wk3.name', 'PA',
            '-rv_wk3.R1', 0.5*0.03, # Scaled from Hyde al. Improvement of right ventricular... 2016
            '-rv_wk3.R2', 0.16*0.63,
            '-rv_wk3.C', 4.1*5.16,
            '-tricuspid_valve.Rfwd', 0.05,
            '-tricuspid_valve.Rbwd', 1000.0,
            '-pulmonic_valve.Rfwd', 0.0]

    return vols

# -----------------------------------------------------------------------------
def setup_solver_mech(args):
    """ adapt solver parameters """

    # for pt the switch is set automatically
    sopts = ['-pstrat', 1,
             '-pstrat_i', 1,
             '-krylov_tol_mech', 1e-10,
             '-krylov_norm_mech', 0,
             '-krylov_maxit_mech', 5000]

    sopts += ['-newton_atol_mech', 1e-6,
              '-newton_tol_mech', 1e-8,
              '-newton_adaptive_tol_mech', 2,
              '-newton_tol_cvsys', 1e-6,
              '-newton_line_search', 0,
              '-newton_maxit_mech', 20]

    solve_opts = ''
    # if args.flv == 'pt':
    #     solve_opts = join(dirname(realpath(__file__)),
    #                               'options/pt_mech_amg')

    if solve_opts != '':
        sopts += ['+F', solve_opts]

    return sopts

# -----------------------------------------------------------------------------
def set_mechanic_options(args):
    """ set mechanical options """

    # Track tissue volume
    mech_opts = ['-volumeTracking', 1,
                 '-numElemVols', 1,
                 '-elemVols[0].name', 'tissue',
                 '-elemVols[0].grid', model.mechanics.grid_id()]

    mech_opts += setup_solver_mech(args)

    return mech_opts

# -----------------------------------------------------------------------------
def setup_bc(mesh_dir,args):
    """ Set up Neumann boundary conditions """
    num_bs = 6 # Veins
    num_nbc = num_bs # Veins

    if args.experiment == 'unloading' and args.periOff:
        num_nbc += 2 # cavities
    elif args.experiment == 'free_contraction':
        num_bs += 1 # Pericardium
        num_nbc += 1
    elif args.experiment == 'wk3' or (args.experiment == 'unloading' and not args.periOff):
        num_nbc += 2 # Cavities
        num_nbc += 1 #Pericardium
        num_bs += 1

        
    nbc = ['-num_mechanic_nbc', num_nbc,
            '-num_mechanic_bs', num_bs]
    
    nbc += ['-mechanic_bs[0].value', args.stiffness_LSPV, # LSPV
            '-mechanic_bs[1].value', 0.0, # LIPV
            '-mechanic_bs[2].value', 0.0, # RIPV
            '-mechanic_bs[3].value', args.stiffness_RSPV, # RSPV
            '-mechanic_bs[4].value', args.stiffness_SVC, # SVC
            '-mechanic_bs[5].value', 0.0] # IVC

    nbc += ['-mechanic_nbc[0].name', 'LSPV',
            '-mechanic_nbc[0].surf_file', join(mesh_dir, 'LS_pulm_vein'),
            '-mechanic_nbc[0].spring_idx', 0,
            '-mechanic_nbc[1].name', 'LIPV',
            '-mechanic_nbc[1].surf_file', join(mesh_dir, 'LI_pulm_vein'),
            '-mechanic_nbc[1].spring_idx', 1,
            '-mechanic_nbc[2].name', 'RIPV',
            '-mechanic_nbc[2].surf_file', join(mesh_dir, 'RI_pulm_vein'),
            '-mechanic_nbc[2].spring_idx', 2,
            '-mechanic_nbc[3].name', 'RSPV',
            '-mechanic_nbc[3].surf_file', join(mesh_dir, 'RS_pulm_vein'),
            '-mechanic_nbc[3].spring_idx', 3,
            '-mechanic_nbc[4].name', 'SVC',
            '-mechanic_nbc[4].surf_file', join(mesh_dir, 'SVC'),
            '-mechanic_nbc[4].spring_idx', 4,
            '-mechanic_nbc[5].name', 'IVC',
            '-mechanic_nbc[5].surf_file', join(mesh_dir, 'IVC'),
            '-mechanic_nbc[5].spring_idx', 5]
    
    if args.experiment == 'unloading' or args.experiment == 'wk3':
        nbc += ['-mechanic_nbc[6].name', 'LV_endo',
                '-mechanic_nbc[6].surf_file', join(mesh_dir, 'LV_endo'),
                '-mechanic_nbc[6].pressure',args.lvEDP_kPA, #kPa
                '-mechanic_nbc[7].name', 'RV_endo',
                '-mechanic_nbc[7].surf_file', join(mesh_dir, 'RV_endo'),
                '-mechanic_nbc[7].pressure',args.rvEDP_kPa] #kPa

    if args.experiment == 'free_contraction' or args.experiment == 'wk3' or (args.experiment == 'unloading' and not args.periOff):
        ed_file = join(mesh_dir, 'PM_' + args.case)
        nbc += ['-num_mechanic_ed', 1,
                '-mechanic_ed[0].ncomp', 1,
                '-mechanic_ed[0].file', ed_file,
                '-mechanic_bs[6].edidx', 0,
                '-mechanic_bs[6].value', args.stiffness_pericardium,
                '-mechanic_nbc[' + str(num_nbc-1) + '].name', 'Pericardium',
                '-mechanic_nbc[' + str(num_nbc-1) + '].surf_file', join(mesh_dir, 'epicardium_ventricles'),
                '-mechanic_nbc[' + str(num_nbc-1) + '].spring_idx', 6,
                '-mechanic_nbc[' + str(num_nbc-1) + '].nspring_idx', 0,
                '-mechanic_nbc[' + str(num_nbc-1) + '].nspring_config',1]




    return nbc

# -----------------------------------------------------------------------------
def setup_active(args):
    """ setup active stress setting """

    if args.postprocess:
        return []
    if args.experiment == 'unloading': # All passive
        opts = ['-num_imp_regions', 1,
                '-num_stim', 0,
                '-imp_region[0].name','passive',
                '-imp_region[0].im','PASSIVE',
                '-imp_region[0].num_IDs',26,
                '-imp_region[0].ID[0]',1,
                '-imp_region[0].ID[1]',2,
                '-imp_region[0].ID[2]',3,
                '-imp_region[0].ID[3]',4,
                '-imp_region[0].ID[4]',5,
                '-imp_region[0].ID[5]',6,
                '-imp_region[0].ID[6]',7,
                '-imp_region[0].ID[7]',8,
                '-imp_region[0].ID[8]',9,
                '-imp_region[0].ID[9]',10,
                '-imp_region[0].ID[10]',11,
                '-imp_region[0].ID[11]',12,
                '-imp_region[0].ID[12]',13,
                '-imp_region[0].ID[13]',14,
                '-imp_region[0].ID[14]',15,
                '-imp_region[0].ID[15]',16,
                '-imp_region[0].ID[16]',17,
                '-imp_region[0].ID[17]',18,
                '-imp_region[0].ID[18]',19,
                '-imp_region[0].ID[19]',20,
                '-imp_region[0].ID[20]',21,
                '-imp_region[0].ID[21]',22,
                '-imp_region[0].ID[22]',23,
                '-imp_region[0].ID[23]',24,
                '-imp_region[0].ID[24]',25,
                '-imp_region[0].ID[25]',26]
    else:
        opts = ['-num_imp_regions', 3,
                '-imp_region[0].name', 'ventricles',
                '-imp_region[0].im', 'TT2',
                '-imp_region[0].num_IDs', 4,
                '-imp_region[0].ID[0]', 1,
                '-imp_region[0].ID[1]', 2,
                '-imp_region[0].ID[2]', 25,
                '-imp_region[0].ID[3]', 26,
                '-imp_region[1].name', 'atria',
                '-imp_region[1].im', 'COURTEMANCHE',
                '-imp_region[1].num_IDs', 2,
                '-imp_region[1].ID[0]', 3,
                '-imp_region[1].ID[1]', 4,
                '-imp_region[2].name', 'others',
                '-imp_region[2].im', 'PASSIVE',
                '-imp_region[2].num_IDs', 20,
                '-imp_region[2].ID[0]', 5,
                '-imp_region[2].ID[1]', 6,
                '-imp_region[2].ID[2]', 7,
                '-imp_region[2].ID[3]', 8,
                '-imp_region[2].ID[4]', 9,
                '-imp_region[2].ID[5]', 10,
                '-imp_region[2].ID[6]', 11,
                '-imp_region[2].ID[7]', 12,
                '-imp_region[2].ID[8]', 13,
                '-imp_region[2].ID[9]', 14,
                '-imp_region[2].ID[10]', 15,
                '-imp_region[2].ID[11]', 16,
                '-imp_region[2].ID[12]', 17,
                '-imp_region[2].ID[13]', 18,
                '-imp_region[2].ID[14]', 19,
                '-imp_region[2].ID[15]', 20,
                '-imp_region[2].ID[16]', 21,
                '-imp_region[2].ID[17]', 22,
                '-imp_region[2].ID[18]', 23,
                '-imp_region[2].ID[19]', 24]

        tanh_pars = set_tanh_stress_pars(args.t_peak, args.Tanh_time_relax, args.Tanh_time_contract)
        opts += ['-imp_region[0].plugins', 'TanhStress',
                 '-imp_region[0].plug_param', tanh_pars]
    return opts

# -----------------------------------------------------------------------------
def setup_em_coupling():
    """
    Setup electromechanical coupling
    """
    # setup weak coupling
    coupling = ['-mech_use_actStress', 1,
                '-mech_lambda_upd', 1,
                '-mech_deform_elec', 0]

    # add velocity dependence fudge factor
    veldep = 0  # fundge factor to attenuate force-velocity dependence in [0,1]
    coupling += ['-veldep', veldep]

    return coupling

# =============================================================================
#    Active stress models
# =============================================================================

# -----------------------------------------------------------------------------
def set_tanh_stress_pars(arg_t_peak, arg_tau_r, arg_tau_c0):
    """
    Tanh stress parameters
    Active stress model as used in Andrew's thesis and
    Niederer et al 2011 Cardiovascular Research 89
    """
    # current settings   #   default values (see TanhStress.c in LIMPET)
    # ------------------ # ----------------------------------------------------
    t_emd = 20           #Marina, but check Electro-mechanical delay, https://doi.org/10.1093/oxfordjournals.eurheartj.a060210
    t_peak = arg_t_peak  # 100.0 kPa peak isometric tension
    tau_c0 = arg_tau_c0      #  40.0 ms  time constant contraction (t_r)
    tau_r = arg_tau_r        # 110.0 ms  time constant relaxation (t_d)
    t_dur = 550.0        # 550.0 ms  duration of transient (t_max)
    lambda_0 = 0.7       #   0.7 -   sacomere length ratio (a_7)
    ld_deg = 6.0         #   5.0 -   degree of length dependence (a_6)
    ld_up = 500.0        # 500.0 ms  length dependence of upstroke time (a_4)
    ld_on = 0            #   0   -   turn on/off length dependence
    vm_thresh = -60.0    # -60.0 mV  threshold Vm for deriving LAT

    tpl = 't_emd={},Tpeak={},tau_c0={},tau_r={},t_dur={},lambda_0={},' + \
          'ld={},ld_up={},ldOn={},VmThresh={}'
    return tpl.format(t_emd, t_peak, tau_c0, tau_r, t_dur, lambda_0, ld_deg,
                      ld_up, ld_on, vm_thresh)

# -----------------------------------------------------------------------------
def setup_visualization(visualize, postprocess):
    """
    Visualization settings
    """
    # prevent output of large files
    vis_opts = ['-gridout_i', 0]

    if visualize:
        if postprocess:
            stress_val = 8  # principal stresses, elementwise
            strain_val = 4+8  # principal strains, elementwise
            vis_opts += ['-mech_output', 1+2,         # igb and vtk
                         '-vtk_output_mode', 3,       # vtu with zlib compr.
                         '-strain_value', strain_val,
                         '-stress_value', stress_val]
        else:
            vis_opts += ['-mech_output', 1+2,
                         '-vtk_output_mode', 3,       # vtu with zlib compr.
                         '-gzip_data', 0,
                         '-strain_value', 0,
                         '-stress_value', 0]
    else:
        vis_opts += ['-mech_output', 1,
                     '-strain_value', 0,
                     '-stress_value', 0]
    return vis_opts

# -----------------------------------------------------------------------------
def setup_time_variables(args):
    """
    Setup of time variables
    """
    time_opts = []
    if args.postprocess:
        output_step = 10
    else:
        output_step = 1

    #tstep = 1.
    tstep = args.mechDT
    mstep = tstep    # mech time step in ms
    estep = 100.     # EP time step in um
    tend = 1000.0    # end time in ms

    # duration not set by autosetting for the specific experiment
    if args.duration != -1:
        tend = args.duration
        if tend < output_step:
            output_step = tend

    time_opts += ['-timedt', tstep,
                  '-mechDT', mstep,
                  '-spacedt', output_step,
                  '-tend', tend,
                  #'-dt', estep,
                  ]
    if args.experiment == 'unloading' or args.experiment == 'wk3':
        time_opts += ['-loadStepping', args.loadStepping]

    return time_opts

# -----------------------------------------------------------------------------
def setup_material(args):
    """
    Material settings
    """
    mech_opts = []

    # set bulk modulus kappa depending on finite element
    kappa = 1000 if args.mech_element == 'P1-P0' else 1e100

    vtrcls = model.mechanics.GuccioneMaterial([1,2,25,26], 'Ventricles',
                                              kappa=kappa, a=args.scaling_Guccione, b_f=8.0,
                                              b_fs=4.0, b_t=3.0)
    atria = model.mechanics.NeoHookeanMaterial([3,4], 'Atria',
                                               kappa=kappa, c=7.45)
    vlvs = model.mechanics.NeoHookeanMaterial([7,8,9,10],
                                              'Valve planes', kappa=kappa, c=1000.0)
    inlets = model.mechanics.NeoHookeanMaterial([11,12,13,14,15,16,17],
                                              'Inlet planes', kappa=kappa, c=1000.0)
    aorta = model.mechanics.NeoHookeanMaterial([5], 'Aorta', kappa=kappa, c=26.66)
    pa = model.mechanics.NeoHookeanMaterial([6], 'Pulm_Artery', kappa=kappa, c=3.7)
    veins = model.mechanics.NeoHookeanMaterial([18,19,20,21,22,23,24], 'BC Veins',
                                               kappa=kappa, c=7.45)


    mech_opts += model.optionlist([vtrcls, atria, vlvs, inlets, aorta, pa, veins])

    mech_opts += ['-mech_vol_split_aniso', 1]

    return mech_opts

# -----------------------------------------------------------------------------
def setup_unloading(args):
    unload_opts = ['-experiment',   5,
                   '-loadStepping', args.loadStepping, 
                   '-unload_conv',  0, # 0  Volume based, 1 point based
                   '-unload_tol',   1e-3,
                   '-unload_err',   1,
                   '-unload_maxit', 10,
                   '-unload_stagtol',10.0]
    return unload_opts

# --- MAIN -------------------------------------------------------------------
if __name__ == '__main__':
    # teardown decorator function and assign it to constant function RUN
    RUN = tools.carpexample(parser_commands, job_id)(run)
    RUN()
