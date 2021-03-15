#!/usr/bin/env python2
# for heart in 21; do /home/crg17/Desktop/scripts/4chmodel/Python/CARP/run_4ch_healthy.py --electrode_name bottom_third --electrode_path /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$heart"/simulations" --experiment EP --mesh_name h_case$heart --mesh_path /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$heart"/meshing/1000um" --TestID $heart"_EP" --FEC_ratio=7 --np 20 --sim_folder /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$heart"/simulations" --fast; done;

import os
from datetime import date
from glob import glob

from carputils import settings
from carputils import tools
from carputils import ep
from carputils import model
#from carputils import testing
from carputils.resources import petsc_block_options

import itertools

MPI_EXEC       ='mpiexec'

#------------------------------------------------------------------
# For EP:
#
# for heart in 21; do /home/crg17/Desktop/scripts/4chmodel/Python/CARP/run_4ch_healthy.py --electrode_name bottom_third --electrode_path /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$heart"/simulations" --experiment EP --mesh_name h_case$heart --mesh_path /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$heart"/meshing/1000um" --TestID $heart"_EP" --FEC_ratio=7 --np 20 --sim_folder /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$heart"/simulations" --fast; done;
#
#------------------------------------------------------------------

#------------------------------------------------------------------
# For unloading:
#
# heart=01; /home/crg17/Desktop/scripts/4chmodel/Python/CARP/run_4ch_healthy.py --experiment unloading --mechBC_type Robin --mesh_name h_case$heart --mesh_path /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$heart/meshing/1000um --np 20 --TestID $heart"_unloading" --sim_folder /home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case$heart"/EP"
#
#------------------------------------------------------------------

def scripthelp():
  print('Script to run an EM simulation with Robin boundary conditions (pericardium).')
  print('Optional arguments: --experiment, --duration, --bulk_CV, --FEC_ratio, --FECRV_off, --C, --alpha, --ld_on, --atria_off, --fast, --coupling, --mechBC_type, --springVeins, --nspring, --stim_file, --TestID, --mv_resistance, --tv_resistance')
  print('Check the help of each argument separately for more information.')
def parser():
    parser = tools.standard_parser()
    parser.add_argument('--experiment',
                        choices=['EP','wk3','free','unloading','inflation'],
                        help='pick experiment type')
    parser.add_argument('--duration',
                        help='duration of experiment (default is 500 ms)')
    parser.add_argument('--bulk_CV',
                        help='CV of ventricular myocardium')
    parser.add_argument('--FEC_ratio',
                        help='Increase of CV in the FEC')
    parser.add_argument('--FECRV_off',
                        action='store_true',
                        help='switch off the FEC in the RV endocardium')
    parser.add_argument('--C',
                        help='Scaling parameter (C1 in Guccione model.) ventricular myocardium')
    parser.add_argument('--alpha',
                        help='Scaling factor for stiffness in f, s and n directions')
    parser.add_argument('--Tpeak',
                        help='Peak in active tension')
    parser.add_argument('--ld_on',
                        action='store_true',
                        help='switch on length-dependence of active force generation in ventricles')
    parser.add_argument('--atria_off',
                        action='store_true',
                        help='Turn off electrics in atria')
    parser.add_argument('--fast',
                        action='store_true',
                        help='Only eikonal solve')
    parser.add_argument('--coupling',
                        choices=['standard','circAdapt'],
                        help='pick coupling model')
    parser.add_argument('--mechBC_type',
                        choices=['Dirichlet','Robin'],
                        help='pick coupling model')
    parser.add_argument('--springVeins',
                        help='Penalty assigned to the displacement on the cropped veins')
    parser.add_argument('--nspring',
                        help='Penalty assigned to the displacement in the normal direction on the pericardium')
    parser.add_argument('--stim_file',
                        help='Give /path/to/file of the activation times file')
    parser.add_argument('--TestID',
                        help='name your test')
    parser.add_argument('--mv_resistance',
                        help='Forward resistance to flow of the mitral valve')
    parser.add_argument('--tv_resistance',
                        help='Forward resistance to flow of the tricuspid valve')
    parser.add_argument('--only_lv',
                        action='store_true',
                        help='Simulate only the left ventricle contraction. By default is a biv simulation.')
    parser.add_argument('--LV_FEC_tag',
                        help='Tag of the FEC layer in the LV.')
    parser.add_argument('--RV_FEC_tag',
                        help='Tag of the FEC layer in the RV.')
    parser.add_argument('--electrode_path')
    parser.add_argument('--electrode_name')
    parser.add_argument('--mesh_path')
    parser.add_argument('--mesh_name')
    parser.add_argument('--sim_folder')
    parser.add_argument('--act_threshold')
    parser.add_argument('--ar_myo')
    parser.add_argument('--restart')
    parser.add_argument('--CV_atria')
    parser.add_argument('--ar_atria')
    parser.add_argument('--LA_tag')
    parser.add_argument('--RA_tag')
    parser.add_argument('--App')
    parser.add_argument('--RIPV')
    parser.add_argument('--LIPV')
    parser.add_argument('--RSPV')
    parser.add_argument('--LSPV')
    parser.add_argument('--SVC')
    parser.add_argument('--IVC')
    parser.add_argument('--PM_name')
    parser.add_argument('--pericardium')
    parser.add_argument('--c_neohook_atria')
    parser.add_argument('--k_neohook_atria')
    parser.add_argument('--MV_tag')
    parser.add_argument('--PV_tag')
    parser.add_argument('--c_neohook_valves')
    parser.add_argument('--k_neohook_valves')
    parser.add_argument('--c_neohook_ao')
    parser.add_argument('--k_neohook_ao')
    parser.add_argument('--c_neohook_PA')
    parser.add_argument('--k_neohook_PA')
    parser.add_argument('--c_neohook_veins')
    parser.add_argument('--k_neohook_veins')
    parser.add_argument('--Ao_tag')
    parser.add_argument('--PA_tag')
    parser.add_argument('--LA_App_tag')
    parser.add_argument('--IVC_tag')
    parser.add_argument('--LA_App_BC_tag')
    parser.add_argument('--IVC_BC_tag')
    parser.add_argument('--newton_atol_mech')
    parser.add_argument('--loadStepping')

    return parser

def jobID(args):
    """
    Generate name of top level output directory.
    """
    args = parameters(args)
    return '{}/{}'.format(args.sim_folder,args.TestID)

@tools.carpexample(parser, jobID)
def run(args, job):

    meshdir = args.mesh_path
    meshname = '{}/{}'.format(args.mesh_path,args.mesh_name)

    cmd  = tools.carp_cmd()

    cmd += ['-simID', job.ID]

    cmd += ['-meshname', meshname]


    if args.experiment != 'EP':

      # set mechanics solver
      mechSolve = setMechanicOptions(args)
      cmd += mechSolve

      # set mechanics BCs
      actSeq = args.stim_file

      if args.experiment != 'unloading':
        if not os.path.isfile(actSeq):
          raise Exception('Activation sequence not provided')
        else:
          stimulus, nStim = ep.fake_ep_opts(0,actSeq)
          cmd += ['-num_stim', nStim]
          cmd += stimulus
          cmd += ['-diffusionOn', 0,
                  '-bidomain', 0 ]

      if args.mechBC_type == 'Dirichlet':
        mechDBC = setupMechDBC(args,meshdir)
        cmd += mechDBC


      mechNBC = setupMechNBC(args,meshdir)
      cmd += mechNBC

      # set mechanics parameters
      mregion = setupMechPar(args)
      cmd += mregion

      imp = impRegions(args)
      cmd += imp

      # set surfaces to keep track on LV and RV volume
      surfVols = setupVolSurf(args,meshdir)
      cmd += surfVols

      mechDT = 1.0

      if args.experiment == 'wk3':

        cav = setupCavities(args)
        cmd += cav

        cmd += ['-loadStepping',args.loadStepping]

      elif args.experiment == 'unloading':
        cmd += ['-experiment',   5,
                '-loadStepping', args.loadStepping, # Up to 100
                '-unload_conv',  0, # 0  Volume based, 1 point based
                '-unload_tol',   args.unload_tol,
                '-unload_err',   1,
                '-unload_maxit', args.unload_maxit,
                '-unload_stagtol',args.unload_stagtol]

    else:

      numStims, stims = setupStimuli(args, meshdir)
      propOpts, numStims, stims = setupPropagation('R-E+', numStims, stims)

      cmd += ['-bidomain',0,
              '-diffusionOn', 1*(1-args.fast),
              '-experiment', 0 + 6*(args.fast)]

      stimOpts = [ '-num_stim', numStims ] + stims
      cmd += stimOpts

      actTime = activationTime(args)
      cmd += actTime

      conduction = conductivity(args)
      cmd += conduction

      imp = impRegions(args)
      cmd += imp

      mechDT = 0.0

      # if args.experiment != 'unloading' and not(args.atria_off):

      #   # Atria/ventricle splitting
      #   cmd += ['-use_splitting', 1,
      #           '-gridout_i'    , 2 ]

      #   # needed for the splitting
      #   cmd += ['-localize_pts', 0]

    cmd += ['-mechDT', mechDT,
            '-dt',     args.dt,
            '-tend',   args.duration]

    # Set active contraction parameters
    act = setupActStress(args)
    cmd += act

    cmd += ['-timedt',       args.timedt,
            '-spacedt',      args.spacedt]

    cmd += ['-mech_output',  1,    # binary vtk
            '-stress_value', 0,    # principal stresses
            '-strain_value', 0]  # strain output

    if args.restart != '':
      cmd += ['-start_statef',args.restart]

    # Run main CARP simulation
    # ------------------------
    job.carp(cmd)


def parameters(args):
    auxiliar_params = {'Ao_tag'        :5,
                       'App'           :False,
                       'AV_tag'        :9,
                       'coupling'      :'standard',
                       'duration'      :1000.0,
                       'electrode_name':'',
                       'electrode_path':'',
                       'experiment'    :'unloading',
                       'IVC'           :False,
                       'IVC_BC_tag'    :24,
                       'IVC_tag'       :17,
                       'LA_App_BC_tag' :18,
                       'LA_App_tag'    :11,
                       'LA_tag'        :3,
                       'LIPV'          :False,
                       'LSPV'          :False,
                       'LV_FEC_tag'    :25,
                       'LV_tag'        :1,
                       'mechBC_type'   :'Dirichlet',
                       'mesh_name'     :'BiV',
                       'mesh_path'     :'/data/01/Mesh/800um/BiV',
                       'MV_tag'        :7,
                       'nspring'       :0.25,
                       'PA_tag'        :6,
                       'pericardium'   :True,
                       'PM_name'       :'UVC_elem',
                       'PV_tag'        :10,
                       'RA_tag'        :4,
                       'restart'       :'',
                       'RIPV'          :True,
                       'RSPV'          :True,
                       'RV_FEC_tag'    :26,
                       'RV_tag'        :2,
                       'sim_folder'    :'/data/simulations',
                       'springVeins'   :0.01, 
                       'stim_file'     :' ',
                       'SVC'           :True,
                       'TestID'        :' ',
                       'TV_tag'        :8
                       }

    unloading_params = {'loadStepping':50,
                        'lv_final_pressure': 1.453, # kPa, final pressure, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4694381/pdf/ijcem0008-18673.pdf
                        'rv_final_pressure' : 0.5, # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4449250/pdf/PulmCirc-005-370.pdf
                        'unload_maxit':10,
                        'unload_stagtol':10.0,
                        'unload_tol'  :1e-03
                        }

    eikonal_params = {'act_threshold':-60, # Doesn't seem to affect eikonal model
                      'ar_atria'     :1, # Don't care by now
                      'ar_FEC'       :1,  # doi: 10.1161/CIRCEP.111.965814,
                      'ar_myo'       :0.29,#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4608725
                      'bulk_CV'      :0.8, #https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4608725
                      'CV_atria'     :0, # Don't care by now
                      'Dstim'        :2.,
                      'Istim'        :150,
                      'FEC_ratio'    :7.0, 
                      'num_regions'  :4 # Ventricles, FECs, atria and others
                      }

    guccione_params = {'C'         :1.7, # https://link.springer.com/content/pdf/10.1007%2Fs10237-016-0865-3.pdf
                       'alpha'     :1,
                       'bf_scalar' :8,
                       'bfs_scalar':4,
                       'bt_scalar' :3,
                       'k_guccione':1000.
                       }

    neohookean_params = {'c_neohook_ao'    :26.66,#Marina's paper
                         'k_neohook_ao'    :1000.,
                         'c_neohook_atria' :7.45, 
                         'k_neohook_atria' :1000.,
                         'c_neohook_PA'    :3.7, 
                         'k_neohook_PA'    :1000.,
                         'c_neohook_valves':1000.,
                         'k_neohook_valves':1000.,
                         'c_neohook_veins' :7.45,
                         'k_neohook_veins' :1000.
                         }

    tanh_stress_params = {'Tpeak':125.0, #Marina, but check https://watermark.silverchair.com/cvq318.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAkwwggJIBgkqhkiG9w0BBwagggI5MIICNQIBADCCAi4GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM0GEFCgr3TmlYZomsAgEQgIIB_2zsbaVjOJgaYyP-Nr9PHVnX2fteJ9TRKx5UUYFj0oMvRzYqTY4CroW9Wlsx2afbwx7BHWc92b4hvKX4T3CxGix1BbQlS768x9K2WFhUX8_mAjlPskOaR-13T3AnAxDAZJNwp4-f0OkUuJpD3Lyw5pVFPgXBPMF0_neKwCVWJNULgohnqcyvIisII15XwxjkB4MXmhKX6MjvqiCeJ2nACtb7UKk6eqCorz6Qcowb7JGJrb4y8EZBguElgD8rBMjTuy8VoaGbmVFY7Y9aq-_OsUhIzcVQN6tEOceTnB2GJQ5rPAHd063PcwD6xPBRtAOt48CeEf-VVuYzJTzmxgumJQIwrc3uSwnQwD6U8kUwscuSruHfDmvqmbEO9YP82TKAQq4ti9GEr5_Wa4DzNvtVbmKxmNU5wL7YffLJHqdswjE87wFql_vdew4vZwKjKHlNupU3tHc7z0UIYXOzEqzpjYBzEcq4L5ZHuuFgodvztqAOJM8xg3VSgC2aGrCkc5cl9K7-z08-UxfcI6qMi6HG0M-U7pU0bWEOFeEEb6Itowlb97ByRoR_5FZZCY_0mrcd_VSl9DwcmiUxSxrqAQvIQ51JHAkEjeEU4SD3RrianCzF2G5-LHqKILr_aMG4AbaYHSiHh6Nrsuc_glyMNoKBjm6Q3MO0WmItDvoIxeOxTOs
                          't_emd':20, #Marina, but check Electro-mechanical delay, https://doi.org/10.1093/oxfordjournals.eurheartj.a060210
                          'tau_c0':130,#time constant of contraction
                          'tau_r':100, #time constant of relaxation, Inheritated
                          't_dur':550, # Inheritated, although it can be estimated by eye from the ref of the Tpeak
                          'lambda_0':0.7, # Tpeak
                          'ld':6, # Inheritated
                          'ld_up':500, # Inheritated
                          'VmThresh':-60
                          }

    windkessel_params = {'mv_resistance':0.05,
                         'tv_resistance':0.05,
                         'pv_resistance':0.0,
                         'Z_lv'         :0.030, # https://link.springer.com/content/pdf/10.1007%2Fs10237-013-0523-y.pdf
                         'R_lv'         :0.63, # As above
                         'C_lv'         :5.16, # As above
                         'p0_out_lv'    :67.5, 
                         'mv_bw_resistance':1000.,
                         'tv_bw_resistance':1000.,
                         'av_resistance':0.0,
                         'p0_out_rv':15.
                         }
    windkessel_params_rv = {'Z_rv':0.5*float(windkessel_params['Z_lv']), # Scaled from Hyde al. Improvement of right ventricular... 2016
                            'R_rv':0.16*float(windkessel_params['R_lv']), 
                            'C_rv':4.1*float(windkessel_params['C_lv'])} 

    numerical_method_params = {'dt'                         :100.0,
                               'parab_solve'                :1,
                               'pstrat'                     :1,
                               'pstrat_i'                   :1,
                               'pstrat_backpermute'         :1,
                               'mapping_mode'               :1,
                               'redist_extracted'           :0,
                               'mass_lumping'               :1,
                               'localize_pts'               :1,
                               'cg_tol_ellip'               :1e-8,
                               'cg_tol_parab'               :1e-8,
                               'cg_tol_mech'                :1e-8,
                               'cg_norm_ellip'              :0,
                               'cg_norm_parab'              :0,
                               'cg_norm_mech'               :0,
                               'cg_maxit_mech'              :1000,
                               'newton_stagnation_detection':0,
                               'newton_atol_mech'           :1e-06,
                               'newton_tol_mech'            :1e-08,
                               'newton_line_search'         :0,
                               'newton_maxit_mech'          :20,
                               'newton_adaptive_tol_mech'   :0,
                               'load_stiffening'            :1,
                               'mech_configuration'         :1,
                               'mech_tangent_comp'          :1,
                               'timedt'                     :1.,
                               'spacedt'                    :1.0
                               }

    params = dict(itertools.chain(auxiliar_params.items(),unloading_params.items(),eikonal_params.items(),guccione_params.items(),neohookean_params.items(),tanh_stress_params.items(),windkessel_params.items(),windkessel_params_rv.items(),numerical_method_params.items()))

    # Keep this parameters unless introduced by user
    for par in params:
      if not hasattr(args,par) or getattr(args,par) is None:
        setattr(args,par,params[par])

    return args



# ============================================================================
#    EP FUNCTIONS
# ============================================================================

def setupPropagation(propagation, numStims, stims):

    propOpts, numStims, stims = ep.setupPropagation(propagation, numStims, stims, '')

    return propOpts, numStims, stims

def setupSolverEP(args):

    sopts = [ '-parab_solve',        args.parab_solve,
              '-pstrat',             args.pstrat,
              '-pstrat_i',           args.pstrat_i,
              '-pstrat_backpermute', args.pstrat_backpermute,
              '-mapping_mode',       args.mapping_mode,
              '-redist_extracted',   args.redist_extracted,
              '-mass_lumping',       args.mass_lumping ]

    return sopts

# --- set conduction velocity and conductivities -----------------------------
def gregion(IDs,num_gregs,name,CV,ani_ratio):
  num_gregs = str(int(num_gregs))

  gregion = ['-gregion['+num_gregs+'].name',    name,
             '-gregion['+num_gregs+'].num_IDs', len(IDs)]

  for idx in range(0,len(IDs)):
    gregion += ['-gregion['+num_gregs+'].ID['+str(idx)+']', IDs[idx]]

  gregion += ['-ekregion['+num_gregs+'].ID',     num_gregs,
              '-ekregion['+num_gregs+'].vel_f',  CV,
              '-ekregion['+num_gregs+'].vel_s',  float(ani_ratio)*float(CV),
              '-ekregion['+num_gregs+'].vel_n',  float(ani_ratio)*float(CV)]

  num_gregs = int(num_gregs) + 1

  return num_gregs, gregion

def conductivity(args):

    num_regs = float(args.num_regions) - 1
    ar_myo = args.ar_myo
    ar_FEC = args.ar_FEC

    CV_a = args.CV_atria*(1-args.atria_off)+0*args.atria_off
    ar_a = args.ar_atria*(1-args.atria_off)+1.0*args.atria_off

    gregions = []

    # ventricles myocardium
    if args.only_lv:
      IDs_vmyo = [int(float(args.LV_tag))]
    elif args.FECRV_off:
      IDs_vmyo = [int(float(args.LV_tag)),int(float(args.RV_tag)),int(float(args.LV_FEC_tag))]
    else:
      IDs_vmyo = [int(float(args.LV_tag)),int(float(args.RV_tag)),int(float(args.LV_FEC_tag)),int(float(args.RV_FEC_tag))]
    num_regs, greg = gregion(IDs_vmyo,num_regs,'ventricles',args.bulk_CV,ar_myo)
    gregions += greg

    # ventricles FEC
    if args.FECRV_off or args.only_lv:
      IDs_FEC = [int(float(args.LV_FEC_tag))]
    else:
      IDs_FEC = [int(float(args.LV_FEC_tag)),int(float(args.RV_FEC_tag))]
    k_FEC = args.FEC_ratio
    num_regs, greg = gregion(IDs_FEC,num_regs,'FEC',float(k_FEC)*float(args.bulk_CV),ar_FEC)
    gregions += greg

    # atria
    num_regs, greg = gregion([int(float(args.LA_tag)),int(float(args.RA_tag))],num_regs,'atria',CV_a,ar_a)
    gregions += greg

    #others 
    num_regs, greg = gregion(range(5,25),num_regs,'others',0,1.0)
    gregions += greg 

    conduction = ['-num_gregions',  num_regs,
                  '-num_ekregions', num_regs] + gregions

    return conduction

# --- set cellular electrical dynamics ---------------------------------------
def impregion(IDs,num_imp_regs,name,reg_type):

  imp_region = ['-imp_region['+str(num_imp_regs)+'].name',    name,
                '-imp_region['+str(num_imp_regs)+'].im',      reg_type,
                '-imp_region['+str(num_imp_regs)+'].num_IDs', len(IDs)]

  for idx in range(0,len(IDs)):
    imp_region += ['-imp_region['+str(int(float(str(num_imp_regs))))+'].ID['+str(idx)+']', IDs[idx]]

  num_imp_regs += 1

  return num_imp_regs, imp_region

def impRegions(args):

    num_imp_regs = 0
    imp_regions = []

    if args.experiment == 'unloading' or args.experiment == 'inflation':
      num_imp_regs, imp_reg = impregion(range(1,27),num_imp_regs,'passive','PASSIVE')
      imp = ['-num_imp_regions',  num_imp_regs] + imp_reg + ['-num_stim', 0]

    else:

      # ventricles
      if args.only_lv:
        IDs_vmyo = [int(float(args.LV_tag))]
      elif args.FECRV_off:
        IDs_vmyo = [int(float(args.LV_tag)),int(float(args.RV_tag)),int(float(args.LV_FEC_tag))]
      else:
        IDs_vmyo = [int(float(args.LV_tag)),int(float(args.RV_tag)),int(float(args.LV_FEC_tag)),int(float(args.RV_FEC_tag))]
      if args.only_lv:
        IDs_valves = [args.MV_tag,args.AV_tag]
      else:
        IDs_valves = [args.MV_tag,args.AV_tag,args.TV_tag,args.PV_tag]

      num_imp_regs, imp_reg = impregion(IDs_vmyo,num_imp_regs,'ventricles','TT2')
      imp_regions += imp_reg

      # atria
      num_imp_regs, imp_reg = impregion([3,4],num_imp_regs,'atria','COURTEMANCHE')
      imp_regions += imp_reg

      # valve planes, arteries and veins
      num_imp_regs, imp_reg = impregion(range(5,25),num_imp_regs,'others','PASSIVE')
      imp_regions += imp_reg 

      imp = ['-num_imp_regions', num_imp_regs] + imp_regions

    return imp

# --- set activation time computation ----------------------------------------
def activationTime(args):

    activation = ['-num_LATs',            1,
                  '-lats[0].measurand',   0,      # Vm (0) or extracellular potential (1)
                  '-lats[0].all',         0,      # All activation times (1) or only the 1st (0)
                  '-lats[0].threshold', float(args.act_threshold),
                  '-lats[0].method',      1 ]     # Threshold crossing

    return activation

# --- set stimulus -----------------------------------------------------------
def stimulus(num_stims,stim_type,electrode_file,start,strength,duration):

    stim = ['-stimulus['+str(num_stims)+'].stimtype', stim_type,
            '-stimulus['+str(num_stims)+'].vtx_file', electrode_file,
            '-stimulus['+str(num_stims)+'].start',    start,
            '-stimulus['+str(num_stims)+'].strength', strength,
            '-stimulus['+str(num_stims)+'].duration', duration]

    num_stims += 1

    return num_stims, stim

def setupStimuli(args,meshdir):

    # Generate electrical trigger options
    Istim    = args.Istim  # strength
    Dstim    = args.Dstim    # duration
    num_stims = 0

    stims = []

    num_stims, stim = stimulus(num_stims,0,args.electrode_path + "/" + args.electrode_name,0.0,Istim*(not args.atria_off),Dstim)
    stims += stim

    return num_stims, stims

# ============================================================================
#    MECHANICS FUNCTIONS
# ============================================================================

# --- set up solver ----------------------------------------------------------
def setupSolverMech(args):

    # for pt the switch is set automatically
    sopts = [ '-localize_pts',args.localize_pts,
              '-cg_tol_ellip',args.cg_tol_ellip,
              '-cg_tol_parab',args.cg_tol_parab,
              '-cg_tol_mech',args.cg_tol_mech,
              '-cg_norm_ellip',args.cg_norm_ellip,
              '-cg_norm_parab',args.cg_norm_parab,
              '-cg_norm_mech',args.cg_norm_mech,
              '-cg_maxit_mech',args.cg_maxit_mech,
              '-newton_stagnation_detection',args.newton_stagnation_detection,
              '-newton_atol_mech',args.newton_atol_mech,
              '-newton_tol_mech',args.newton_tol_mech,
              '-newton_line_search',args.newton_line_search,
              '-newton_maxit_mech',args.newton_maxit_mech,
              '-newton_adaptive_tol_mech',args.newton_adaptive_tol_mech]

    return sopts

# --- set mechanical options -------------------------------------------------
def setMechanicOptions(args):

    mech_opts = ['-load_stiffening',     args.load_stiffening,
                 '-mech_configuration',  args.mech_configuration,
                 '-mech_tangent_comp',   args.mech_tangent_comp]

    mech_opts += setupSolverMech(args)

    return mech_opts

# --- set mechanics boundary conditions --------------------------------------
def mech_dbc(num_dbc,name,vtx_file):

    dbc = ['-mechanic_dbc['+str(num_dbc)+'].name',     name,
           '-mechanic_dbc['+str(num_dbc)+'].bctype',   0,
           '-mechanic_dbc['+str(num_dbc)+'].apply_ux', 1,
           '-mechanic_dbc['+str(num_dbc)+'].apply_uy', 1,
           '-mechanic_dbc['+str(num_dbc)+'].apply_uz', 1,
           '-mechanic_dbc['+str(num_dbc)+'].vtx_file', vtx_file]

    num_dbc += 1

    return num_dbc, dbc

def setupMechDBC(args,meshdir):
    num_mechanic_dbc = 0
    mechDBC = []

    # Pulmonary vein 1
    if args.App:
      num_mechanic_dbc, dbc = mech_dbc(num_mechanic_dbc,'LA Appendage',meshdir+'/BC/LA_App_BC.surf')
      mechDBC += dbc

    # Pulmonary vein 2
    if args.RIPV:
      num_mechanic_dbc, dbc = mech_dbc(num_mechanic_dbc,'RIPV',meshdir+'/BC/RIPV_BC.surf')
      mechDBC += dbc

    if args.RSPV:
      num_mechanic_dbc, dbc = mech_dbc(num_mechanic_dbc,'RSPV',meshdir+'/BC/RSPV_BC.surf')
      mechDBC += dbc

    # Pulmonary vein 3
    if args.LIPV:
      num_mechanic_dbc, dbc = mech_dbc(num_mechanic_dbc,'LIPV',meshdir+'/BC/LIPV_BC.surf')
      mechDBC += dbc

    # Pulmonary vein 4
    if args.LSPV:
      num_mechanic_dbc, dbc = mech_dbc(num_mechanic_dbc,'LSPV',meshdir+'/BC/LSPV_BC.surf')
      mechDBC += dbc

    # Superior vena cava
    if args.SVC:
      num_mechanic_dbc, dbc = mech_dbc(num_mechanic_dbc,'SVC',meshdir+'/BC/SVC_BC.surf')
      mechDBC += dbc

    # Inferior vena cava
    if args.IVC:
      num_mechanic_dbc, dbc = mech_dbc(num_mechanic_dbc,'IVC',meshdir+'/BC/IVC_BC.surf')
      mechDBC += dbc

    return ['-num_mechanic_dbc', num_mechanic_dbc] + mechDBC

def mech_nbc(num_nbc,name,pressure,surf_file,spring_idx=-1,nspring_idx=-1,load_spring_idx=-1,load_nspring_idx=-1,trace=""):

    num_nbc = int(num_nbc)


    nbc = ['-mechanic_nbc['+str(num_nbc)+'].name',      name,
           '-mechanic_nbc['+str(num_nbc)+'].surf_file', surf_file,
           '-mechanic_nbc['+str(num_nbc)+'].pressure',  pressure]

    if spring_idx+nspring_idx+load_spring_idx+load_nspring_idx > -4:
      nbc += ['-mechanic_nbc['+str(num_nbc)+'].spring_idx',        spring_idx,
              '-mechanic_nbc['+str(num_nbc)+'].nspring_idx',       nspring_idx,
              '-mechanic_nbc['+str(num_nbc)+'].dpspring_idx',      -1,
              '-mechanic_nbc['+str(num_nbc)+'].load_spring_idx',   load_spring_idx,
              '-mechanic_nbc['+str(num_nbc)+'].load_nspring_idx',  load_nspring_idx,
              '-mechanic_nbc['+str(num_nbc)+'].load_dpspring_idx', -1]

    if trace != "":
      nbc += ['-mechanic_nbc['+str(num_nbc)+'].trace', trace]

    num_nbc += 1

    return num_nbc, nbc

def elem_data(num_elem_data,n_components,elem_file):

    ed = ['-mechanic_ed['+str(num_elem_data)+'].ncomp', n_components,
          '-mechanic_ed['+str(num_elem_data)+'].file',  elem_file]

    num_elem_data += 1

    return num_elem_data, ed

def spring(num_s,value,elem_data_idx=-1,scale=1,active=1):

    s = ['-mechanic_bs['+str(num_s)+'].value',  value,
         '-mechanic_bs['+str(num_s)+'].edidx',  elem_data_idx,
         '-mechanic_bs['+str(num_s)+'].scale',  scale,
         '-mechanic_bs['+str(num_s)+'].active', active]

    num_s += 1

    return num_s, s

def setupSprings(args):

  num_springs = 0
  springs = []

  # Epicardium
  num_springs, s = spring(num_springs,args.nspring,elem_data_idx=0)
  springs += s

  # # Pulmonary vein 1
  # num_springs, s = spring(num_springs,float(args.springVeins)*float(args.App))
  # springs += s

  # Pulmonary vein 2
  num_springs, s = spring(num_springs,float(args.springVeins)*float(args.RIPV))
  springs += s

# Pulmonary vein 3
  num_springs, s = spring(num_springs,float(args.springVeins)*float(args.LIPV))
  springs += s

# Pulmonary vein 4
  num_springs, s = spring(num_springs,float(args.springVeins)*float(args.LSPV))
  springs += s

# Pulmonary vein 5
  num_springs, s = spring(num_springs,float(args.springVeins)*float(args.RSPV))
  springs += s

# Superior vena cava
  num_springs, s = spring(num_springs,float(args.springVeins)*float(args.SVC))
  springs += s

# Inferior vena cava
  num_springs, s = spring(num_springs,float(args.springVeins)*float(args.IVC))
  springs += s

  return ['-num_mechanic_bs', num_springs] + springs

def setupMechNBC(args,meshdir):

  num_mechanic_nbc = 2 - float(bool(args.only_lv))

  if args.only_lv:
    mechNBC = ['-mechanic_nbc[0].name',      'LV_endo',
              '-mechanic_nbc[0].pressure', float(args.lv_final_pressure) * float(bool(args.experiment == 'unloading')),
              '-mechanic_nbc[0].surf_file', meshdir+'/cavities/LV_endo_closed']
  else:
    mechNBC = ['-mechanic_nbc[0].name',      'LV_endo',
              '-mechanic_nbc[0].pressure', float(args.lv_final_pressure) * float(bool(args.experiment == 'unloading')),
              '-mechanic_nbc[0].surf_file', meshdir+'/cavities/LV_endo_closed',
              '-mechanic_nbc[1].name',      'RV_endo',
              '-mechanic_nbc[1].pressure', float(args.rv_final_pressure) * float(bool(args.experiment == 'unloading')),
              '-mechanic_nbc[1].surf_file', meshdir + '/cavities/RV_endo_closed']

  if args.mechBC_type == 'Robin':

    # define element data file
    num_elem_data, ed = elem_data(0,1,args.mesh_path + "/" + args.PM_name)
    el_data = ['-num_mechanic_ed', num_elem_data] + ed
    mechNBC += el_data

    # define the springs
    # if(args.experiment is not 'EP' and args.experiment is not 'unloading'):
    springs = setupSprings(args)
    mechNBC += springs

    # define nbc for the springs
    springidx = 0
    # Pericardium

    # Turn the pericardium on only if it is not unloading - for the unloading we turn it off
    if args.experiment != 'unloading' and args.pericardium:
      num_mechanic_nbc, nbc = mech_nbc(num_mechanic_nbc,'Pericardium',0.0,meshdir+'/BC/epicardium_ventricles',nspring_idx=springidx)
      mechNBC += nbc
      springidx += 1

  # Pulmonary vein 1
  if args.mechBC_type != 'Dirichlet':
      num_mechanic_nbc, nbc = mech_nbc(num_mechanic_nbc,'RIPV',0.0,meshdir+'/BC/RIPV_BC',spring_idx=springidx,load_spring_idx=springidx)
      mechNBC += nbc
      springidx += 1

    # Pulmonary vein 2
      num_mechanic_nbc, nbc = mech_nbc(num_mechanic_nbc,'LIPV',0.0,meshdir+'/BC/LIPV_BC',spring_idx=springidx,load_spring_idx=springidx)
      mechNBC += nbc
      springidx += 1

    # Pulmonary vein 3
      num_mechanic_nbc, nbc = mech_nbc(num_mechanic_nbc,'LSPV',0.0,meshdir+'/BC/LSPV_BC',spring_idx=springidx,load_spring_idx=springidx)
      mechNBC += nbc
      springidx += 1

    # Pulmonary vein 4
      num_mechanic_nbc, nbc = mech_nbc(num_mechanic_nbc,'RSPV',0.0,meshdir+'/BC/RSPV_BC',spring_idx=springidx,load_spring_idx=springidx)
      mechNBC += nbc
      springidx += 1

    # SVC
      num_mechanic_nbc, nbc = mech_nbc(num_mechanic_nbc,'SVC',0.0,meshdir+'/BC/SVC_BC',spring_idx=springidx,load_spring_idx=springidx)
      mechNBC += nbc
      springidx += 1

    # IVC
      num_mechanic_nbc, nbc = mech_nbc(num_mechanic_nbc,'IVC',0.0,meshdir+'/BC/IVC_BC',spring_idx=springidx,load_spring_idx=springidx)
      mechNBC += nbc
      springidx += 1

  return ['-num_mechanic_nbc',num_mechanic_nbc] + mechNBC

# --- set parameters of the constitutive law ----------------------------
def mregion(num_mregs,name,IDs,mreg_type,k,c,bf=None,bfs=None,bt=None):

    mreg = ['-mregion['+str(num_mregs)+'].name',    name,
            '-mregion['+str(num_mregs)+'].num_IDs', len(IDs)]

    for idx in range(0,len(IDs)):
      mreg += ['-mregion['+str(num_mregs)+'].ID['+str(idx)+']', int(float(IDs[idx]))]

    mreg += ['-mregion['+str(num_mregs)+'].type', mreg_type]

    if mreg_type == 7:
      mreg += ['-mregion['+str(num_mregs)+'].params', 'kappa='+str(k)+',C='+str(c)+',b_f='+str(bf)+',b_fs='+str(bfs)+',b_t='+str(bt)]
    elif mreg_type == 1:
      mreg += ['-mregion['+str(num_mregs)+'].params', 'kappa='+str(k)+',c='+str(c)]
    else:
      raise Exception('Unknown material type')

    num_mregs += 1

    return num_mregs, mreg

def setupMechPar(args):

    num_mregions = 0
    mregions = []

    alpha_float = float(args.alpha)

    bf = float(args.bf_scalar)*alpha_float
    bfs = float(args.bfs_scalar)*alpha_float
    bt = float(args.bt_scalar)*alpha_float

    # Ventricles
    if args.only_lv:
      ventr_mregions = [args.LV_tag,args.LV_FEC_tag]
    else:
      ventr_mregions = [args.LV_tag,args.LV_FEC_tag,args.RV_tag,args.RV_FEC_tag]
    num_mregions, mreg = mregion(num_mregions,'Ventricles',ventr_mregions,7,k=float(args.k_guccione),c=args.C,bf=bf,bfs=bfs,bt=bt)
    mregions += mreg

        # Atria
    num_mregions, mreg = mregion(num_mregions,'Atria',[args.LA_tag,args.RA_tag],1,k=args.k_neohook_atria,c=args.c_neohook_atria)
    mregions += mreg

    # Valve planes and veins inlets
    valve_planes = range(int(float(args.MV_tag)),int(float(args.IVC_tag)) + 1)
    num_mregions, mreg = mregion(num_mregions,'Valve planes',valve_planes,1,k=float(args.k_neohook_valves),c=float(args.c_neohook_valves))
    mregions += mreg

    # Aorta
    num_mregions, mreg = mregion(num_mregions,'Aorta',[float(args.Ao_tag)],1,k=float(args.k_neohook_ao),c=float(args.c_neohook_ao))
    mregions += mreg

    # Pulmonary artery
    num_mregions, mreg = mregion(num_mregions,'Pulm artery',[float(args.PA_tag)],1,k=float(args.k_neohook_PA),c=float(args.c_neohook_PA))
    mregions += mreg

    # SVC, IVC and pulmonary veins boundary conditions
    BC_vector = range(int(float(args.LA_App_BC_tag)),int(float(args.IVC_BC_tag)) + 1)
    num_mregions, mreg = mregion(num_mregions,'Boundary conditions',BC_vector,1,k=float(args.k_neohook_veins),c=float(args.c_neohook_veins))
    mregions += mreg

    return ['-num_mregions',num_mregions] + mregions

# --- set active stress model ------------------------------------------------
def setupActStress(args):

    act = ['-veldep', 0.]

    act += ['-imp_region[0].plugins',    'TanhStress',
            '-imp_region[0].plug_param', 't_emd='+str(args.t_emd)+',Tpeak='+str(args.Tpeak)+',tau_c0='+str(args.tau_c0)+',tau_r='+str(args.tau_r)+',t_dur='+str(args.t_dur)+',lambda_0='+str(args.lambda_0)+',ld='+str(args.ld)+',ld_up='+str(args.ld_up)+',VmThresh='+str(args.VmThresh)+',ldOn='+str(int(args.ld_on))]
            #'-imp_region[0].im_sv_init', './states/TT2_TanhStress_bcl_500_ms.sv']

    act += ['-mech_use_actStress', int(bool(args.experiment != 'unloading'))]

    return act

# --- set cavities -----------------------------------------------------------
def setupVolSurf(args,meshdir):

    numSurfVols = 4

    surfVols = ['-surfVols[0].name',      'LV_endo',
                '-surfVols[0].surf_file', meshdir+'/cavities/LV_endo_closed',
                '-surfVols[0].grid',      7,
                '-surfVols[1].name',      'RV_endo',
                '-surfVols[1].surf_file', meshdir+'/cavities/RV_endo_closed',
                '-surfVols[1].grid',      7,
                '-surfVols[2].name',      'LA_endo',
                '-surfVols[2].surf_file', meshdir+'/cavities/LA_endo_closed',
                '-surfVols[2].grid',      7,
                '-surfVols[3].name',      'RA_endo',
                '-surfVols[3].surf_file', meshdir+'/cavities/RA_endo_closed',
                '-surfVols[3].grid',      7]

    return ['-numSurfVols',numSurfVols] + surfVols

def setupCavities(args):

  numCavities = 2 - float(args.only_lv)



# Kerchoffs et al 2003
#  Z = 0.09
#  R = 0.9
#  C = 0.1867

#  Zr = 0.045
#  Rr = 0.1260
#  Cr = 0.7655

  # Ratio from Niederer et al 2010
#  Zr = 0.35 * Z
#  Rr = 0.125 * R
#  Cr = 4.5 * C

  state = -1

  cvs_mode    = float(args.coupling=='circAdapt')  # circAdapt based coupling
  valve_model = cvs_mode  # circAdapt valve


  cav = ['-CV_coupling', 0,
         '-CVS_mode',    cvs_mode]

  cav += ['-num_cavities',        numCavities,
         '-cavities[0].cav_type', 0,
         '-cavities[0].cavP',     0,
         '-cavities[0].tube',     2,
         '-cavities[0].cavVol',   0,
         '-cavities[0].p0_cav',   7.5*float(args.lv_final_pressure),
         '-cavities[0].p0_in',    7.5*float(args.lv_final_pressure),
         '-cavities[0].p0_out',   args.p0_out_lv,
         '-cavities[0].valve',    valve_model,
         '-cavities[0].state',    str(state)]

  cav += ['-lv_wk3.name', 'Aorta',
          '-lv_wk3.R1',   args.Z_lv,
          '-lv_wk3.R2',   args.R_lv,
          '-lv_wk3.C',    args.C_lv]

  # set valve resistance
  cav += ['-mitral_valve.Rfwd', args.mv_resistance,
 	        '-mitral_valve.Rbwd', args.mv_bw_resistance,
          '-aortic_valve.Rfwd', args.av_resistance]

  cav += ['-cavities[1].cav_type', 1,
          '-cavities[1].cavP',     1,
          '-cavities[1].tube',     2,
          '-cavities[1].cavVol',   1,
          '-cavities[1].p0_cav',   7.5*float(args.rv_final_pressure),
          '-cavities[1].p0_in',    7.5*float(args.rv_final_pressure),
          '-cavities[1].p0_out',   args.p0_out_rv,
          '-cavities[1].valve',    valve_model,
          '-cavities[1].state',    str(state)]

  cav += ['-rv_wk3.name', 'PA',
          '-rv_wk3.R1',   args.Z_rv,
          '-rv_wk3.R2',   args.R_rv,
          '-rv_wk3.C',    args.C_rv ]

  # set valve resistance
  cav += ['-tricuspid_valve.Rfwd', args.tv_resistance,
          '-tricuspid_valve.Rbwd', args.tv_bw_resistance,
          '-pulmonic_valve.Rfwd', args.pv_resistance]

  return cav

'''
test_basic = testing.Test('ring-basic', run, ['--duration', '100', '--np','4' ],
                             tags=[testing.tag.MEDIUM, testing.tag.PARALLEL])
test_basic.add_filecmp_check('x.dynpt', testing.max_error, 0.001)

__tests__ = [ test_basic ]
'''

if __name__ == '__main__':
    run()
