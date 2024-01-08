######################################################################
# ref_path and output have to be provided as a string
#
# Execute shepherd like this: ./shepherd.py ateles_params.py
#
#####################################################################
import os
import sys
import datetime
import shutil

from clone_and_build_function import *

# Set this switch to true to abort the recheck when the first job fails.
abort = False

templateFolder = './templates/'
machineFolder = './machines/'
apesFolder = os.path.expandvars('${HOME}/apes')

date = datetime.datetime.now().strftime("%Y-%m-%d__%X")
weekday = datetime.datetime.now().strftime("%A")

# Production directory, keep the past week as history.
prod_dir = 'ateles-runs_' + weekday

run_label = 'ATELES'

# Cleanup production directory before using it:
shutil.rmtree(prod_dir, ignore_errors=True)
loglevel = 'INFO'

# source for the mercurial functions
git_clone_source = 'https://github.com/apes-suite/'

# mail adress
from recheck import notify_atl, mail_server
mail_address = notify_atl
smtp_server = mail_server

# name of the shepherd log file
shepherd_out = 'shepherd.log'

# name of the log and rror file of the clone and build function
clone_build_out = 'clone_build.log'
clone_build_err = 'clone_build_error.log'

# Set this to true to have the current revision marked as working when all tests
# succeed.
create_tag_on = False
# Set this to true to have shepherd store the performance results in the loris
# repository.
grep_performance = True

loris_clone_url = apesFolder + 'loris/'

# path to the testsuite dir to shorten the string in the job_dict
atldir = os.path.join(apesFolder, 'ateles', 'atl', 'examples')


loris_clone_url = apesFolder + 'loris/'


shepherd_jobs = []

seeder_exe = clone_build( solver          = 'seeder',
                          revision        = 'main',
                          git_clone_source = git_clone_source+'seeder.git',
                          solver_dir      = 'seeder',
                          clone_build_out = clone_build_out,
                          clone_build_err = clone_build_err         )

ateles_exe = clone_build( solver          = 'ateles',
                          git_clone_source = git_clone_source+'ateles.git',
                          solver_dir      = 'ateles',
                          clone_build_out = clone_build_out,
                          clone_build_err = clone_build_err         )

## ATELES JOB 1
#Checked# (HK)
casedir = os.path.join( atldir, 'lineareuler', '3D',
                        'transientBackground' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        prefix = 'lineareuler_3D_transientBackground',
        label =  'lineareuler_3D_transientBackground_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_lineareuler_transientBackground_p00000.res' ),
        val_output_filename =
                'lineareuler_transientBackground_p00000.res'
    )
)

## ATELES JOB 2
casedir = os.path.join( atldir, 'euler', '3D',
                        'shear_layer_Q4' )

shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir, 'hyperfun.lua' ),
        prefix = 'shear_layer_Q4',
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        label = 'shearlayer_modg_Q4_valid_hyperfun',
    )
)
shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir,'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        depend = ['shearlayer_modg_Q4_valid_hyperfun'],
        label = 'shearlayer_modg_Q4_valid_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir,'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['shearlayer_modg_Q4_valid_seeder','shearlayer_modg_Q4_valid_hyperfun'],
        label = 'shearlayer_modg_Q4_valid_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_shear_layer_modg_probe_momentum_Q4_p00000.res' ),
        val_output_filename = 'shear_layer_modg_probe_momentum_Q4_p00000.res',
    )
)

## ATELES JOB 3
casedir = os.path.join( atldir, 'euler', '3D',
                        'shear_layer_Q8' )

shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir, 'hyperfun.lua' ),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        prefix = 'shearlayer_modg_Q8_valid',
        label = 'shearlayer_modg_Q8_valid_hyperfun',
    )
)
shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir,'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        depend = ['shearlayer_modg_Q8_valid_hyperfun'],
        label = 'shearlayer_modg_Q8_valid_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir,'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['shearlayer_modg_Q8_valid_seeder'],
        label = 'shearlayer_modg_Q8_valid_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_shear_layer_modg_probe_momentum_Q8_p00000.res'),
        val_output_filename = 'shear_layer_modg_probe_momentum_Q8_p00000.res',
    )
)
## ATELES JOB 4
casedir = os.path.join( atldir, 'euler', '3D',
                        'shear_layer_Q4_parallel' )

shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir, 'hyperfun.lua' ),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        prefix = 'shear_layer_Q4_parallel',
        label = 'parallel_shearlayer_modg_Q4_valid_hyperfun',
    )
)
shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir,'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        depend = ['parallel_shearlayer_modg_Q4_valid_hyperfun'],
        create_subdir = ['mesh'],
        label = 'parallel_shearlayer_modg_Q4_valid_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir,'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 2',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['parallel_shearlayer_modg_Q4_valid_seeder','parallel_shearlayer_modg_Q4_valid_hyperfun'],
        label = 'parallel_shearlayer_modg_Q4_valid_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_shear_layer_modg_probe_momentum_Q4_p00000.res' ),
        val_output_filename = 'shear_layer_modg_probe_momentum_Q4_p00000.res',
    )
)

## ATELES JOB 5
casedir = os.path.join( atldir, 'maxwell', 'periodic_oscillator_Q8_3D')

shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir,'posci_common.lua'),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        prefix = 'posci_modg_Q8_valid',
        label = 'posci_modg_Q8_valid_common',
    )
)
shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir,'valid_tracking.lua'),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        label = 'posci_modg_Q8_valid_tracking',
        depend = ['posci_modg_Q8_valid_common'],
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_Q8_valid.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        label = 'posci_modg_Q8_valid_ateles',
        depend = ['posci_modg_Q8_valid_tracking', 'posci_modg_Q8_valid_common'],
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_posci_modgQ8_probe_electricField_Q8_p00000.res' ),
        val_output_filename = 'posci_modgQ8_probe_electricField_Q8_p00000.res',
    )
)

## ATELES JOB 6
casedir = os.path.join( atldir, 'maxwell', 'periodic_oscillator_P8_3D')

shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir,'posci_common.lua'),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        prefix = 'posci_modg_P8_valid',
        label = 'posci_modg_P8_valid_common',
    )
)
shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir,'valid_tracking.lua'),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        label = 'posci_modg_P8_valid_tracking',
        depend = ['posci_modg_P8_valid_common'],
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_P8_valid.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        label = 'posci_modg_P8_valid_ateles',
        depend = ['posci_modg_P8_valid_common','posci_modg_P8_valid_tracking'],
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_posci_modgP8_probe_electricField_P8_p00000.res' ),
        val_output_filename = 'posci_modgP8_probe_electricField_P8_p00000.res',
    )
)

## ATELES JOB 7
casedir = os.path.join( atldir, 'maxwell', 'rectangular_wave_Q8_3D')

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'rectwave_modg_Q8',
        label = 'rectwave_modg_Q8_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_maxwell_modg_valid.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['rectwave_modg_Q8_seeder'],
        label = 'rectwave_modg_Q8_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_rectang_waveguide_maxwell_modg_probe_electricField_Q8_p00000.res' ),
        val_output_filename = 'rectang_waveguide_maxwell_modg_probe_electricField_Q8_p00000.res',
    )
)

## ATELES JOB 8
casedir = os.path.join( atldir, 'maxwell', 'pec_scatter_inhomogenous_matFun_Q8_3D')

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'pecScat_inhomMatFun_modg_Q8',
        label = 'pecScat_inhomMatFun_modg_Q8_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_maxwell_modg.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 4',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['pecScat_inhomMatFun_modg_Q8_seeder'],
        label = 'pecScat_inhomMatFun_modg_Q8_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_pec_scatter_maxwell_modg_probe_displacementField_Q8_p00000.res' ),
        val_output_filename = 'pec_scatter_maxwell_modg_probe_displacementField_Q8_p00000.res',
    )
)

## ATELES JOB 9
casedir = os.path.join( atldir, 'maxwell', 'pec_scatter_inhomogenous_matVar_Q8_3D')

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'pecScat_inhomMatVar_modg_Q8',
        label = 'pecScat_inhomMatVar_modg_Q8_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_maxwell_modg.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 4',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['pecScat_inhomMatVar_modg_Q8_seeder'],
        label = 'pecScat_inhomMatVar_modg_Q8_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_pec_scatter_maxwell_modg_probe_displacementField_Q8_p00000.res' ),
        val_output_filename = 'pec_scatter_maxwell_modg_probe_displacementField_Q8_p00000.res',
    )
)

## ATELES JOB 10
casedir = os.path.join( atldir, 'euler', '1D', 'toro1_x' )

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension='lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'toro1_x_modg_1d_Q4',
        label = 'toro1_x_modg_Q4_1d_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['toro1_x_modg_Q4_1d_seeder'],
        label = 'toro1_x_modg_Q4_1d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro1_x_euler_modg_1d_probe_density_Q4_toro_x_p00000.res' ),
        val_output_filename = 'toro1_x_euler_modg_1d_probe_density_Q4_toro_x_p00000.res',
    )
)

## ATELES JOB 11
casedir = os.path.join( atldir, 'euler', '2D', 'toro1_x' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'toro1_x_modg_Q4_2d',
        label = 'toro1_x_modg_Q4_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro1_x_euler_modg_2d_probe_density_Q4_toro_x_p00000.res' ),
        val_output_filename = 'toro1_x_euler_modg_2d_probe_density_Q4_toro_x_p00000.res',
    )
)
## ATELES JOB 12
casedir = os.path.join( atldir, 'euler', '2D', 'toro1_y' )

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'toro1_y_modg_Q4_2d',
        label = 'toro1_y_modg_Q4_2d_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['toro1_y_modg_Q4_2d_seeder'],
        label = 'toro1_y_modg_Q4_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro1_y_euler_modg_2d_probe_density_Q4_toro_y_p00000.res' ),
        val_output_filename = 'toro1_y_euler_modg_2d_probe_density_Q4_toro_y_p00000.res',
    )
)
## ATELES JOB 13
casedir = os.path.join( atldir, 'euler', '3D', 'toro1_x' )

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension='lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'toro1_x_modg_Q4',
        label = 'toro1_x_modg_Q4_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['toro1_x_modg_Q4_seeder'],
        label = 'toro1_x_modg_Q4_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro1_x_euler_modg_probe_density_Q4_toro_x_p00000.res' ),
        val_output_filename = 'toro1_x_euler_modg_probe_density_Q4_toro_x_p00000.res',
    )
)

## ATELES JOB 14
casedir = os.path.join( atldir, 'euler', '3D', 'toro1_y' )

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'toro1_y_modg_Q4',
        label = 'toro1_y_modg_Q4_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['toro1_y_modg_Q4_seeder'],
        label = 'toro1_y_modg_Q4_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro1_y_euler_modg_probe_density_Q4_toro_y_p00000.res' ),
        val_output_filename = 'toro1_y_euler_modg_probe_density_Q4_toro_y_p00000.res',
    )
)

## ATELES JOB 15
casedir = os.path.join( atldir, 'euler', '3D', 'toro1_z' )

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'toro1_z_modg_Q4',
        label = 'toro1_z_modg_Q4_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['toro1_z_modg_Q4_seeder'],
        label = 'toro1_z_modg_Q4_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro1_z_euler_modg_probe_density_Q4_toro_z_p00000.res' ),
        val_output_filename = 'toro1_z_euler_modg_probe_density_Q4_toro_z_p00000.res',
    )
)
## ATELES JOB 16
casedir = os.path.join( atldir, 'euler', '2D', 'toro2_x' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'toro2_modg_Q1_2d_SSPRK2',
        label = 'toro2_modg_Q1_2d_SSPRK2_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro2_euler_modg_2d_probe_density_Q1_toro_p00000.res' ),
        val_output_filename = 'toro2_euler_modg_2d_probe_density_Q1_toro_p00000.res',
    )
)
## ATELES JOB 17
casedir = os.path.join( atldir, 'euler', '2D', 'toro3_x' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'toro3_modg_Q1_2d',
        label = 'toro3_modg_Q1_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro3_euler_modg_2d_probe_density_toro_p00000.res' ),
        val_output_filename = 'toro3_euler_modg_2d_probe_density_toro_p00000.res',
    )
)
## ATELES JOB 18
casedir = os.path.join( atldir, 'euler', '2D', 'toro4_x' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'toro4_modg_Q1_2d',
        label = 'toro4_modg_Q1_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_toro4_euler_modg_2d_probe_density_toro_p00000.res' ),
        val_output_filename = 'toro4_euler_modg_2d_probe_density_toro_p00000.res',
    )
)


## ATELES JOB 19
casedir = os.path.join( atldir, 'euler', '3D', 'gauss_densitypulse_l2p' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'gaussian_pulse_l2p',
        label = 'gaussian_pulse_l2p_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_gPulseDens_euler_modg_track_momentum_l2p_p00000.res' ),
        val_output_filename = 'gPulseDens_euler_modg_track_momentum_l2p_p00000.res',
    )
)
## ATELES JOB 20
casedir = os.path.join( atldir, 'euler', '3D', 'gauss_densitypulse_fxt' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'gaussian_pulse_fxt',
        label = 'gaussian_pulse_fxt_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_gPulseDens_euler_fxt_track_momentum_fxt_p00000.res' ),
        val_output_filename = 'gPulseDens_euler_fxt_track_momentum_fxt_p00000.res',
    )
)

## ATELES JOB 21
casedir = os.path.join( atldir, 'euler', '3D', 'gauss_densitypulse_fpt' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'gaussian_pulse',
        label = 'gaussian_pulse_fpt_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_gPulseDens_euler_modg_track_momentum_p00000.res' ),
        val_output_filename = 'gPulseDens_euler_modg_track_momentum_p00000.res',
    )
)

## ATELES JOB 22
casedir = os.path.join( atldir, 'maxwelldivcorr', 'periodic_oscillator_P8_3D')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_divcor_p.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        label = 'maxwell_DivCor_ateles',
        prefix = 'maxDC_periodic_oscillator_P8_3D',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_maxwell_divcor_probe_electricField_P8_p00000.res' ),
        val_output_filename = 'maxwell_divcor_probe_electricField_P8_p00000.res',
    )
)

## ATELES JOB 23
casedir = os.path.join( atldir, 'maxwelldivcorr', 'periodic_oscillator_Q8_3D')

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder_source.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'maxwell_DivCor_source',
        label = 'maxwell_DivCor_source_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_source_Q8_recheck.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 6',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['maxwell_DivCor_source_seeder'],
        label = 'maxwell_DivCor_source_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 'ref_maxwell_source_divcor_source_probe_electricField_Q8_p00000.res' ),
        val_output_filename = 'maxwell_source_divcor_source_probe_electricField_Q8_p00000.res',
    )
)

## ATELES JOB 24
casedir = os.path.join( atldir, 'maxwellpml', 'cylinder_scattering_Q8_2D')

shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir, 'hyperfun.lua'),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        prefix = 'maxwell_circlematerial_pml_source',
        label = 'maxwell_circlematerial_pml_source_hyperfun',
    )
)
shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        depend = ['maxwell_circlematerial_pml_source_hyperfun'],
        label = 'maxwell_circlematerial_pml_source_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_recheck_Q8_maxwell_2d.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 5',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['maxwell_circlematerial_pml_source_seeder'],
        label = 'maxwell_circlematerial_pml_source_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_simulation_probe_electricField_Q8_pml_circMaterial_p00000.res' ),
        val_output_filename = 'simulation_probe_electricField_Q8_pml_circMaterial_p00000.res',
    )
)

## ATELES JOB 25
casedir = os.path.join( atldir, 'euler', '3D',
                        'gauss_densitypulse_deriv_quantity')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'euler_check_derVar',
        label = 'euler_check_derVar_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_gPulseDens_euler_modg_track_ke_p00000.res' ),
        val_output_filename = 'gPulseDens_euler_modg_track_ke_p00000.res',
    )
)
## ATELES JOB 26
casedir = os.path.join( atldir, 'euler', '3D',
                        'gauss_densitypulse_spongelayer')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'spongeLayer',
        label = 'spongeLayer_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_sponge_layer_modg_sponge_p00000.res' ),
        val_output_filename = 'sponge_layer_modg_sponge_p00000.res',
    )
)

## ATELES JOB 27
casedir = os.path.join( atldir, 'navierstokes', '2D',
                        'viscous_vortex')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'viscousVortex_Q4_2d',
        label = 'viscousVortex_Q4_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_simulation_track_momentum_p00000.res' ),
        val_output_filename = 'simulation_track_momentum_p00000.res',
    )
)

## ATELES JOB 28
casedir = os.path.join( atldir, 'heat', '3D',
                        'sinus_temperature')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension='lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'Heat_3D',
        label = 'Heat_3D_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_simulation_track_temp_p00000.res' ),
        val_output_filename = 'simulation_track_temp_p00000.res',
    )
)
## ATELES JOB 29
casedir = os.path.join( atldir, 'heat', '2D',
                        'sinus_temperature')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_heat_2d.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'Heat_2D',
        label = 'Heat_2D_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_simulation_track_temp_p00000.res' ),
        val_output_filename = 'simulation_track_temp_p00000.res',
    )
)

## ATELES JOB 30
casedir = os.path.join( atldir, 'heat', '1D',
                        'sinus_temperature')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_heat_1d.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'Heat_1D',
        label = 'Heat_1D_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_simulation_track_temp_p00000.res' ),
        val_output_filename = 'simulation_track_temp_p00000.res',
    )
)

## ATELES JOB 31
casedir = os.path.join( atldir, 'euler', '2D',
                        'gauss_densitypulse_spongelayer')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'spongeLayer_2d',
        label = 'spongeLayer_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_sponge_2d_modg_sponge_2d_p00000.res' ),
        val_output_filename = 'sponge_2d_modg_sponge_2d_p00000.res',
    )
)
## ATELES JOB 32
casedir = os.path.join( atldir, 'euler', '3D',
                        'vorticity')

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'vort',
        label = 'vort_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_vorticity_modg_vort_p00000.res' ),
        val_output_filename = 'vorticity_modg_vort_p00000.res',
    )
)

## ATELES JOB 33
casedir = os.path.join( atldir, 'euler', '2D',
                        'vorticity')
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'vort_2d',
        label = 'vort_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_vorticity_2d_modg_vort_p00000.res' ),
        val_output_filename = 'vorticity_2d_modg_vort_p00000.res',
    )
)

## ATELES JOB 34
casedir = os.path.join( atldir, 'acoustic', 'planar_wave_3D' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles_recheck.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'acoustic_3d',
        label = 'acoustic_3d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_plane_wave_XQ8_track_line_density_m7_p00000.res' ),
        val_output_filename = 'plane_wave_XQ8_track_line_density_m7_p00000.res',
    )
)
## ATELES JOB 35
casedir = os.path.join( atldir, 'acoustic', 'gauss_densitypulse_2D' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'acoustic_modg_2d_recheck.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'acoustic_2d',
        label = 'acoustic_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_simulation_track_density_p00000.res' ),
        val_output_filename = 'simulation_track_density_p00000.res',
    )
)
## ATELES JOB 36
casedir = os.path.join( atldir, 'euler', '3D',
                        'shock_stabilization_parallel' )

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension='lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'covolume_modg_Q4_3d_rktaylor',
        label = 'covolume_modg_Q4_3d_rktaylor_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 2',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['covolume_modg_Q4_3d_rktaylor_seeder'],
        label = 'covolume_modg_Q4_3d_rktaylor_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_euler_3d_probe_density_Q4_covolume_rktaylor_z_p00000.res' ),
        val_output_filename = 'euler_3d_probe_density_Q4_covolume_rktaylor_z_p00000.res',
    )
)
## ATELES JOB 37
casedir = os.path.join( atldir, 'euler', '3D',
                        'shock_stabilization_localrefinement' )
shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'covolume_modg_3d_localrefinement_serial',
        label = 'covolume_modg_3d_localrefinement_serial_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['covolume_modg_3d_localrefinement_serial_seeder'],
        label = 'covolume_modg_3d_localrefinement_serial_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_euler_3d_probe_density_covolume_multilevel_p00000.res' ),
        val_output_filename = 'euler_3d_probe_density_covolume_multilevel_p00000.res',
    )
)

## ATELES JOB 38
casedir = os.path.join( atldir, 'euler', '3D',
                        'shock_stabilization_periodic' )

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'covolume_modg_3d_localrefinement_parallel',
        label = 'covolume_modg_3d_localrefinement_parallel_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 3',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['covolume_modg_3d_localrefinement_parallel_seeder'],
        label = 'covolume_modg_3d_localrefinement_parallel_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_euler_3d_probe_density_Q4_periodic_covolume_z_p00000.res' ),
        val_output_filename = 'euler_3d_probe_density_Q4_periodic_covolume_z_p00000.res',
    )
)

## ATELES JOB 39
casedir = os.path.join( atldir, 'euler', '2D',
                        'shock_stabilization_parallel' )

shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        prefix = 'covolume_modg_Q4_2d',
        label = 'covolume_modg_Q4_2d_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 2',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['covolume_modg_Q4_2d_seeder'],
        label = 'covolume_modg_Q4_2d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_euler_2d_probe_density_Q4_covolume_rktaylor_y_p00000.res' ),
        val_output_filename = 'euler_2d_probe_density_Q4_covolume_rktaylor_y_p00000.res',
    )
)
## ATELES JOB 40
casedir = os.path.join( atldir, 'euler', '1D',
                        'shock_stabilization_parallel' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 2',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'covolume_modg_Q4_1d_rktaylor',
        label = 'covolume_modg_Q4_1d_rktaylor_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_euler_1d_probe_density_Q4_covolume_rktaylor_x_p00000.res' ),
        val_output_filename = 'euler_1d_probe_density_Q4_covolume_rktaylor_x_p00000.res',
    )
)

## ATELES JOB 41
casedir = os.path.join( atldir, 'euler', '3D',
                        'qCriterion_deriv_quantity' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'track_qCriterion_3d',
        label = 'track_qCriterion_3d_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_q_crit_3d_vort_p00000.res' ),
        val_output_filename = 'q_crit_3d_vort_p00000.res',
    )
)

## ATELES JOB 42
casedir = os.path.join( atldir, 'navierstokes', '2D',
                        'constant_state' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        prefix = 'track_NS_const_state',
        label =  'track_NS_const_state_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,'ref_simulation_track_const_state_l2p_p00000.res' ),
        val_output_filename = 'simulation_track_const_state_l2p_p00000.res',
    )
)

## ATELES JOB 43
casedir = os.path.join( atldir, 'euler', '3D',
                        'p_polynomials' )
shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir, 'hyperfun.lua' ),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        prefix = 'euler_3D_p_polynomials',
        label = 'euler_3D_p_polynomials_hyperfun',
    )
)
shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        create_subdir = ['mesh'],
        depend = ['euler_3D_p_polynomials_hyperfun'],
        label = 'euler_3D_p_polynomials_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['euler_3D_p_polynomials_seeder'],
        label = 'euler_3D_p_polynomials_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_shear_layer_modg_probe_momentum_P8_p00000.res', ),
        val_output_filename =
                'shear_layer_modg_probe_momentum_P8_p00000.res',
    )
)

## ATELES JOB 44
casedir = os.path.join( atldir, 'euler', '3D',
                        'multilevel' )
shepherd_jobs.append(
    dict(
        executable = seeder_exe,
        template = os.path.join( casedir, 'seeder.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        prefix = 'euler_3D_multilevel',
        label = 'euler_3D_multilevel_seeder',
        mail = False
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        depend = ['euler_3D_multilevel_seeder'],
        label = 'euler_3D_multilevel_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_matml_reflected_pulse_microphone_p00000.res', ),
        val_output_filename =
                'matml_reflected_pulse_microphone_p00000.res',
    )
)

## ATELES JOB 45
casedir = os.path.join( atldir, 'euler', '1D',
                        'modalEstimate' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'euler_1d_modalEstimate',
        label = 'euler_1d_modalEstimate_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_modalest_1d_point_series_p00000.res', ),
        val_output_filename =
                'modalest_1d_point_series_p00000.res'
    )
)

## ATELES JOB 46
casedir = os.path.join( atldir, 'euler', '2D',
                        'modalEstimate' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'euler_2d_modalEstimate',
        label = 'euler_2d_modalEstimate_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_modalest_2d_point_series_p00000.res', ),
        val_output_filename =
                'modalest_2d_point_series_p00000.res'
    )
)

## ATELES JOB 48
casedir = os.path.join( atldir, 'euler', '3D', 'overview' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'euler_3d_overview',
        label = 'euler_3d_overview_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_pp_PointProbe_p00000.res' ),
        val_output_filename =
                'pp_PointProbe_p00000.res'
    )
)

## ATELES JOB 49
casedir = os.path.join( atldir, 'navierstokes', '3D',
                        'shear_tube' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'navierstokes_3d_sheartube',
        label = 'navierstokes_3d_sheartube_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_shear_tube_track_momentum_p00000.res', ),
        val_output_filename =
                'shear_tube_track_momentum_p00000.res'
    )
)

## ATELES JOB 50
casedir = os.path.join( atldir, 'navierstokes', '3D',
                        'shear_hat' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'navierstokes_3d_shearhat',
        label = 'navierstokes_3d_shearhat_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_simulation_track_shearhat3D_p00000.res', ),
        val_output_filename =
                'simulation_track_shearhat3D_p00000.res'
    )
)

## ATELES JOB 51
casedir = os.path.join( atldir, 'navierstokes', '2D',
                        'shear_hat' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'navierstokes_2d_shearhat',
        label = 'navierstokes_2d_shearhat_ateles',
        attachment = True,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_simulation_track_shearhat2D_p00000.res', ),
        val_output_filename =
                'simulation_track_shearhat2D_p00000.res'
    )
)

## ATELES JOB 52
casedir = os.path.join( atldir, 'lineareuler', '3D',
                        'gradient' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        prefix = 'lineareuler_3D_gradient',
        label =  'lineareuler_3D_gradient_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_linearEuler_gradients_track_grads_p00000.res', ),
        val_output_filename =
                'linearEuler_gradients_track_grads_p00000.res',
    )
)

## ATELES JOB 53
casedir = os.path.join( atldir, 'lineareuler', '2D',
                        'transientBackground' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        prefix = 'lineareuler_2d_transientBackground',
        label =  'lineareuler_2d_transientBackground_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
            'ref_simulation_track_2d_density_temporalBackground_p00000.res', ),
        val_output_filename =
            'simulation_track_2d_density_temporalBackground_p00000.res',
    )
)

## ATELES JOB 54
casedir = os.path.join( atldir, 'lineareuler', '3D',
                        'gauss_pulse' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        prefix = 'lineareuler_3d_gausspulse',
        label =  'linearEuler_3d_gausspulse_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_simulation_track_density_p00000.res', ),
        val_output_filename = 'simulation_track_density_p00000.res',
    )
)

## ATELES JOB 55
casedir = os.path.join( atldir, 'lineareuler', '2D',
                        'gauss_pulse' )
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua' ),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        prefix = 'lineareuler_2d_gausspulse',
        label =  'linearEuler_2d_gausspulse_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir,
                'ref_simulation_track_2d_density_p00000.res', ),
        val_output_filename = 'simulation_track_2d_density_p00000.res',
    )
)
## ATELES JOB 56
casedir = os.path.join( atldir, 'euler', '1D', 'piston' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'piston_1D',
        label = 'piston_1D_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 
            'Ref_line_p00000_t5.000E-06.res' ),
        val_output_filename = 'ateles_line_p00000_t5.000E-06.res',
    )
)
## ATELES JOB 57
casedir = os.path.join( atldir, 'euler', '1D', 'piston_modereduction' )

shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        prefix = 'piston_1D_modereduction',
        label = 'piston_modereduction_1D_ateles',
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 
            'Ref_modereduction_line_p00000_t5.000E-06.res'),
        val_output_filename = 'ateles_line_p00000_t5.000E-06.res',
    )
)
## ATELES JOB 58
casedir = os.path.join( atldir, 'euler', '2D', 'wedge' )

shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir,'wedge.lua'),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        prefix = 'wedge_2D_configure',
        label = 'wedge_configure_ateles',
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 5',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        label = 'wedge_2D',
        depend = ['wedge_configure_ateles'],
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 
            'Ref_ateles_point_p00000.res' ),
        val_output_filename = 'ateles_point_p00000.res',
    )
)
## ATELES JOB 59 
casedir = os.path.join( atldir, 'euler', '2D', 'wedge_modereduction' )

shepherd_jobs.append(
    dict(
        executable = None,
        mail = False,
        template = os.path.join( casedir,'wedge.lua'),
        extension = 'lua',
        run_exec = False,
        abort_failure = abort,
        prefix = 'wedge_modereduction_2D_configure',
        label = 'wedge_modereduction_ateles',
    )
)
shepherd_jobs.append(
    dict(
        executable = ateles_exe,
        solver_name = 'ateles',
        template = os.path.join( casedir, 'ateles.lua'),
        extension = 'lua',
        run_exec = True,
        run_command = 'mpirun -np 5',
        abort_failure = abort,
        additional_params = dict(testsuite_path=atldir),
        label = 'wedge_modereduction_2D',
        depend = ['wedge_modereduction_ateles'],
        attachment = False,
        validation = True,
        val_method = 'difference',
        val_ref_path = os.path.join( casedir, 
            'Ref_ateles_point_p00000.res' ),
        val_output_filename = 'ateles_point_p00000.res',
    )
)

