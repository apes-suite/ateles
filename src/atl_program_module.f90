! Copyright (c) 2015-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016-2017, 2020 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
!
! Permission to use, copy, modify, and distribute this software for any
! purpose with or without fee is hereby granted, provided that the above
! copyright notice and this permission notice appear in all copies.
!
! THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
! WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
! ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
! WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
! ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
! OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
! **************************************************************************** !

!> author: Verena Krupp
!! This module has all the general function for ateles programm:
!! initialize, solve and finialze
!! gathering the process into 3 main routines for ateles to unified eg. for
!! apesmate
!!
module atl_program_module

  use mpi

  use env_module,                         only: rk, pathLen

  use aotus_module,                       only: close_config

  use treelmesh_module,                   only: treelmesh_type, free_treelmesh
  use tem_tools_module,                   only: tem_horizontalSpacer
  use tem_aux_module,                     only: tem_open_distconf_array, &
    &                                           tem_abort
  use tem_element_module,                 only: eT_fluid
  use tem_logging_module,                 only: logUnit
  use tem_simControl_module,              only: tem_simControl_clearStat,  &
    &                                           tem_simControl_syncUpdate, &
    &                                           tem_simControl_dump_now
  use tem_tracking_module,                only: tem_tracking_finalize
  use tem_restart_module,                 only: tem_restart_finalize
  use tem_timer_module,                   only: tem_startTimer, &
    &                                           tem_stopTimer
  use tem_trackmem_module,                only: tem_trackmem_load, &
    &                                           tem_trackmem
  use tem_status_module,                  only: tem_status_dump,        &
    &                                           tem_stat_interval,      &
    &                                           tem_status_run_end,     &
    &                                           tem_status_run_terminate
  use tem_timeControl_module,             only: tem_timeControl_dump
  use tem_convergence_module,             only: tem_convergence_check
  use tem_operation_var_module,           only: &
    & tem_opVar_reduction_transient_update

  ! Use precice interface:
  use tem_precice_module,                 only: precice_available,    &
    &                                           tem_precice_advance,  &
    &                                           tem_precice_finalize, &
    &                                           tem_precice_ongoing,  &
    &                                           precice_handle


  use atl_restart_module,                 only: atl_writeRestartIfNecessary
  use atl_initialize_module,              only: atl_initialize
  use atl_container_module,               only: atl_element_container_type
  use atl_equation_module,                only: atl_equations_type
  use atl_calc_time_module,               only: atl_get_timestep
  use atl_solver_param_module,            only: atl_solver_param_type,   &
    &                                           atl_load_solver_parameters
  use atl_source_module,                  only: atl_deallocate_sourceData
  use atl_scheme_module,                  only: atl_modg_scheme_prp, &
    &                                           atl_modg_2d_scheme_prp, &
    &                                           atl_modg_1d_scheme_prp
  use atl_kerneldata_module,              only: atl_kerneldata_update_estimates
  use atl_timer_module,                   only: atl_dumptimers, &
    &                                           atl_timerHandles, &
    &                                           atl_resetTimers
  use atl_weights_module,                 only: atl_dumpWeights
  use atl_global_time_integration_module, only: atl_initTimeStepInfo
  use atl_physCheck_module,               only: atl_check_val
  use atl_varSys_module,                  only: atl_varSys_solverData_type

  use ply_poly_project_module,            only: ply_poly_project_type
  use ply_sampled_tracking_module,        only: ply_sampled_track_output

  implicit none

  private


  public :: atl_load_config
  public :: atl_initialize_program
  public :: atl_solve_program
  public :: atl_finalize_program

contains

  ! ****************************************************************************
  ! routine which loads the config file and checks wether
  ! params%general%solver%configFile is given by apesmate or the filename is
  ! get from command line
  ! further loading can be done here
  subroutine atl_load_config(params, tree, configFile)
    ! --------------------------------------------------------------------------
    ! The structure that holds the solver parameter
    type(atl_solver_param_type), intent(inout) :: params

    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(out) :: tree

    !> Filename to use as default, when no argument provided.
    !!
    !! defaults to 'ateles.lua' if not provided.
    character(len=*), intent(in), optional :: configFile
    ! --------------------------------------------------------------------------
    integer :: nMaxThreads
    character(len=PathLen) :: filename
    ! --------------------------------------------------------------------------

    if (present(configFile)) then
      filename = trim(configFile)
    else
      filename = 'ateles.lua'
    end if

    ! Check whether params%general%solver%configFile is defined
    ! (coming from apesmate)
    ! If not get the configuration file as first parameter
    ! from the command line,
    ! if none is specified, assume the configuration to be located in the
    ! current working directory and named ateles.lua.
    if ( trim(params%general%solver%configFile) == '' ) then
      call get_command_argument(1,params%general%solver%configFile)
      if (trim(params%general%solver%configFile) == '') then
        ! Warn about usage of default configuration file name.
        ! Logging unit not initialized yet
        if (params%general%proc%rank == 0) then
          write(logUnit(1),*) 'No Lua configuration file specified, using'
          write(logUnit(1),*) trim(filename) // ' as default input file.'
        end if
        params%general%solver%configFile = trim(filename)
      end if
    end if


    if (params%general%proc%rank == 0) then
      write(logunit(1),*) "Reading configuration file " &
        & // trim(params%general%solver%configFile)
    end if

    nMaxThreads = 1
    allocate(params%general%solver%conf(nMaxThreads))


    ! Open the configuration file and store a handle to it in conf.
    call tem_open_distconf_array(                          &
      & L        = params%general%solver%conf,             &
      & filename = trim(params%general%solver%configFile), &
      & proc     = params%general%proc                     )

    ! Load the parameters of the solver
    ! This includes the information, if the simulation should restart
    ! from a previous run (readRestart).
    call atl_load_solver_parameters(params, tree)

  end subroutine atl_load_config
  ! ****************************************************************************

  ! ****************************************************************************
  subroutine atl_initialize_program( params, equation, tree, varSys_data,   &
    &                                nCellsNoBnd, element_container,        &
    &                                poly_proj_list, precice_dt             )
    ! --------------------------------------------------------------------------
    ! The structure that holds the solver parameter
    type(atl_solver_param_type), intent(inout) :: params
    ! Description of the equation system to solve
    type(atl_equations_type), intent(out) :: equation
    ! The treelmesh data structure
    type(treelmesh_type), intent(inout)           :: tree
    !> Data Infomation of the variable System
    type(atl_varSys_solverData_type), intent(out) :: varSys_data
    ! Number of cells on each levels
    integer, allocatable, intent(out)             :: nCellsNoBnd(:)
    ! Main data structure of Ateles describing the mesh elements
    type(atl_element_container_type), intent(out) :: element_container
    ! Desribe the projetion methods for the polynomials
    type(ply_poly_project_type), allocatable, intent(out) :: poly_proj_list(:)
    ! Timestep specified from  precice
    real(kind=rk), intent(out), optional :: precice_dt
    ! --------------------------------------------------------------------------
    ! Iterating over the various levels
    integer :: iList, iLevel
    integer :: firstLevel, lastlevel
    integer :: mindegree
    ! --------------------------------------------------------------------------
    ! start initialization timer
    call tem_startTimer( timerHandle = atl_timerHandles%init )


    call atl_initialize(                    &
      & params         = params,            &
      & container      = element_container, &
      & equation       = equation,          &
      & tree           = tree,              &
      & varSys_data    = varSys_data,       &
      & timestepLimit  = precice_dt,        &
      & poly_proj_list = poly_proj_list     )

    call tem_trackmem_load( me   = params%general%solver%trackmem_file, &
      &                     conf = params%general%solver%conf(1)        )

    call tem_trackmem(params%general%solver%trackmem_file, 0)

    nCellsNoBnd = element_container &
      &             %cubes          &
      &             %mesh_list(:)   &
      &             %descriptor     &
      &             %elem           &
      &             %nElems(eT_fluid)

    ! Write out the initial data if restart is active.
    call atl_writeRestartIfNecessary(              &
      &    restart    = params%general%restart,    &
      &    simControl = params%general%simControl, &
      &    equation   = equation,                  &
      &    tree       = tree,                      &
      &    mesh       = element_container%cubes    )

    ! Create the necessary info for the timestep control
    do iList = tree%global%minlevel, tree%global%maxLevel
      call atl_kerneldata_update_estimates(                              &
        &    statedata  = element_container%cubes%statedata_list(iList), &
        &    kerneldata = element_container%cubes%kerneldata_list(iList) )
      call atl_initTimeStepInfo(                                      &
        & equation   = equation,                                      &
        & mesh       = element_container%cubes%mesh_list(iList),      &
        & poly_proj  = poly_proj_list(element_container%cubes%        &
        &              poly_proj_pos(iList)),                         &
        & timestep   = element_container%time%elementSteps(iList),    &
        & control    = element_container%time%control,                &
        & statedata  = element_container%cubes%statedata_list(iList), &
        & kerneldata = element_container%cubes%kerneldata_list(iList) )
    end do

    if ( params%plySampleTrack%tracking%control%active ) then
      ! Write a tracking, if required.

      allocate(params%var_degree(equation%varsys%varname%nVals))
      allocate(params%lvl_degree(element_container%cubes%maxlevel))
      allocate(params%var_space(equation%varsys%varname%nVals))

      mindegree = -1
      firstLevel = lbound(element_container%cubes%scheme_list,1)
      lastLevel  = ubound(element_container%cubes%scheme_list,1)

      select case(element_container%cubes%scheme_list(firstLevel)%scheme)
      case (atl_modg_scheme_prp)
        params%var_degree = maxval( element_container%cubes%scheme_list(:) &
          &                                          %modg%maxPolyDegree   )
        do iLevel = firstLevel, lastLevel
          params%lvl_degree(iLevel) = element_container          &
            &                         %cubes%scheme_list(iLevel) &
            &                         %modg%maxPolyDegree
        end do
        mindegree = minval( element_container%cubes%scheme_list(:) &
          &                                  %modg%maxPolyDegree   )
        params%var_space = element_container%cubes%scheme_list(firstLevel) &
          &                                       %modg%basisType

      case (atl_modg_2d_scheme_prp)
        params%var_degree = maxval( element_container%cubes%scheme_list(:)  &
          &                                          %modg_2d%maxPolyDegree )
        do iLevel = firstLevel, lastLevel
          params%lvl_degree(iLevel) = element_container          &
            &                         %cubes%scheme_list(iLevel) &
            &                         %modg_2d%maxPolyDegree
        end do
        mindegree = minval( element_container%cubes%scheme_list(:)  &
          &                                  %modg_2d%maxPolyDegree )
        params%var_space = element_container%cubes%scheme_list(firstLevel) &
          &                                       %modg_2d%basisType

      case (atl_modg_1d_scheme_prp)
        params%var_degree = maxval( element_container%cubes%scheme_list(:)  &
          &                                          %modg_1d%maxPolyDegree )
        do iLevel = firstLevel, lastLevel
          params%lvl_degree(iLevel) = element_container          &
            &                         %cubes%scheme_list(iLevel) &
            &                         %modg_1d%maxPolyDegree
        end do
        mindegree = minval( element_container%cubes%scheme_list(:)  &
          &                                  %modg_1d%maxPolyDegree )
        params%var_space = element_container%cubes%scheme_list(firstLevel) &
          &                                       %modg_1d%basisType

      end select

      if (mindegree /= params%var_degree(1)) then
        write(logunit(1),*) 'ERROR: Tracking does not support different' &
          &                 // ' polynomial degrees!'
        write(logunit(1),*) 'See:'
        write(logunit(1),*) 'https://github.com/apes-suite/ateles-source/' &
          &                 // 'issues/6'
        write(logunit(1),*) 'Stopping...'
        call tem_abort()
      end if

      call ply_sampled_track_output( me         = params%plySampleTrack,    &
        &                            mesh       = tree,                     &
        &                            bc         = params%boundary,          &
        &                            solver     = params%general%solver,    &
        &                            proc       = params%general%proc,      &
        &                            varSys     = equation%varSys,          &
        &                            var_degree = params%var_degree,        &
        &                            lvl_degree = params%lvl_degree,        &
        &                            var_space  = params%var_space,         &
        &                            simControl = params%general%simControl )

    end if

    call tem_stopTimer( timerHandle = atl_timerHandles%init )

    ! Calculate the global timestep by a CFL condition including cfl
    ! calculation, mpi reduction to get the minimal timestep in
    ! parallel runs, checking for NAN timesteps and checking for
    ! the precice timestep.
    ! The timestep computation is done here in the initialization to
    ! have it available in apesmate if needed.
    if ( .not. precice_available ) then
      call atl_get_timestep(                                          &
        & tree              = tree,                                   &
        & mesh_list         = element_container%cubes%mesh_list,      &
        & scheme_list       = element_container%cubes%scheme_list,    &
        & material_list     = element_container%cubes%material_list,  &
        & equation          = equation,                               &
        & time              = element_container%time,                 &
        & statedata_list    = element_container%cubes%statedata_list, &
        & nCellsNoBnd       = nCellsNoBnd,                            &
        & general           = params%general,                         &
        & adaptive_timestep = equation%adaptive_timestep,             &
        & initial           = .true.                                  )
    else
      call atl_get_timestep(                                          &
        & tree              = tree,                                   &
        & mesh_list         = element_container%cubes%mesh_list,      &
        & scheme_list       = element_container%cubes%scheme_list,    &
        & material_list     = element_container%cubes%material_list,  &
        & equation          = equation,                               &
        & time              = element_container%time,                 &
        & statedata_list    = element_container%cubes%statedata_list, &
        & nCellsNoBnd       = nCellsNoBnd,                            &
        & general           = params%general,                         &
        & adaptive_timestep = equation%adaptive_timestep,             &
        & initial           = .true.,                                 &
        & precice_dt        = precice_dt                              )

    end if

    write(logUnit(1),*) 'Starting compute loop ...'
    write(logUnit(1),*) 'Step size (dt): ', &
      & element_container%cubes%scheme_list(tree%global%minLevel)%time%dt
    call tem_timeControl_dump(params%general%simControl%timeControl,logUnit(1))
    call tem_horizontalSpacer(funit = logUnit(1))
  end subroutine atl_initialize_program
  ! ****************************************************************************

  ! ****************************************************************************
  subroutine atl_solve_program( params, equation, tree, nCellsNoBnd,          &
    &                           element_container, poly_proj_list, precice_dt )
    ! --------------------------------------------------------------------------
    ! The structure that holds the solver parameter
    type(atl_solver_param_type), intent(inout) :: params
    ! Description of the equation system to solve
    type(atl_equations_type), intent(inout) :: equation
    ! The treelmesh data structure
    type(treelmesh_type), intent(in)           :: tree
    ! Number of cells on each levels
    integer, intent(inout)             :: nCellsNoBnd(:)
    ! Main data structure of Ateles describing the mesh elements
    type(atl_element_container_type), intent(inout) :: element_container
    ! Desribe the projetion methods for the polynomials
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    ! Timestep specified from  precice
    real(kind=rk), intent(inout), optional :: precice_dt
    ! --------------------------------------------------------------------------
    ! Iterating over the various levels
    integer :: iList

    !> integer which specifies is the coupling partner is still working
    integer :: is_ongoing=1
    ! --------------------------------------------------------------------------
    ! terminal output
    call tem_trackmem(params%general%solver%trackmem_file, 0)
    call tem_startTimer( timerHandle = atl_timerHandles%simLoop )

    ! !!!
    ! !!! now Compute !!!
    ! !!!  TIMELOOP   !!!
    ! !!!
    time_iter: do
      call tem_simControl_clearStat(params%general%simControl)

      ! Print iteration and timestep (useful to find too small timesteps/NaNs)
      write(logunit(10), *) 'Begin of iteration'
      call tem_simControl_dump_now(params%general%simControl, logUnit(10))


      ! Update the time stepping method, i.e. transfer the new timestep to it.
      do iList = tree%global%minLevel, tree%global%maxLevel
        call element_container%time%elementSteps(iList)%updateStep(    &
          & me           = element_container%time%elementSteps(iList), &
          & timestepInfo = element_container%cubes%scheme_list(iList)  &
          &                                       %time                )
      end do

      ! Do one timestep for the whole mesh. We call meshStep on the coarsest
      ! level and it will call itself recursively until it reaches the finest
      ! level.

      call element_container%time%meshStep(                        &
        &    minLevel       = tree%global%minLevel,                &
        &    maxLevel       = tree%global%maxLevel,                &
        &    currentLevel   = tree%global%minLevel,                &
        &    cubes          = element_container%cubes,             &
        &    tree           = tree,                                &
        &    timestep_list  = element_container%time%elementSteps, &
        &    nSteps         = element_container%time%nSteps,       &
        &    equation       = equation,                            &
        &    general        = params%general,                      &
        &    commStateTimer = atl_timerHandles%commState,          &
        &    poly_proj_list = poly_proj_list                       )


      ! Check if the result of the previous timestep is physical
      call tem_startTimer( timerHandle = atl_timerHandles%checkVal)
      call atl_check_val(                                          &
        & minlevel       = tree%global%minlevel,                   &
        & maxlevel       = tree%global%maxlevel,                   &
        & statedata_list = element_container%cubes%statedata_list, &
        & mesh_list      = element_container%cubes%mesh_list,      &
        & stat           = params%general%simControl%status,       &
        & equation       = equation,                               &
        & scheme_list    = element_container%cubes%scheme_list,    &
        & poly_proj_pos  = element_container%cubes%poly_proj_pos,  &
        & poly_proj_list = poly_proj_list,                         &
        & check          = element_container%physCheck,            &
        & iteration      = params%general%simcontrol%now%iter,     &
        & time           = element_container%time                  )
      call tem_stopTimer( timerHandle = atl_timerHandles%checkVal )
      do iList = tree%global%minLevel, tree%global%maxLevel
        call atl_kerneldata_update_estimates(                              &
          &    statedata  = element_container%cubes%statedata_list(iList), &
          &    kerneldata = element_container%cubes%kerneldata_list(iList) )
      end do

      ! check for convergence only if abortCriteria is set for steady_state
      if( params%general%simControl%abortCriteria%steady_state) then
        ! check for steady state convergence based on convergence criteria
        ! defined in convergence table
        call tem_startTimer( timerHandle = atl_timerHandles%convergeCheck)
        call tem_convergence_check(                                          &
          &                 me     = params%general%simControl%abortCriteria &
          &                                %convergence,                     &
          &                 time   = params%general%simControl%now,          &
          &                 status = params%general%simControl%status,       &
          &                 varSys = equation%varSys,                        &
          &                 tree   = tree                                    )
        call tem_stopTimer( timerHandle = atl_timerHandles%convergeCheck)
      end if

      ! Update reduction trasnient operation variables
      if (equation%redTransVarMap%varPos%nVals>0) then
        call tem_opVar_reduction_transient_update(                         &
          & redTransVarPos = equation%redTransVarMap%varPos                &
          &                  %val(1:equation%redTransVarMap%varPos%nVals), &
          & varSys         = equation%varSys,                              &
          & tree           = tree,                                         &
          & time           = params%general%simControl%now                 )
      end if

      ! for dt the timestep on the minLevel is choosen, since currently it is
      ! equal over alll levels
      !call tem_startTimer( timerHandle = atl_timerHandles%syncUpdate)
      call tem_simControl_syncUpdate(                        &
        & me      = params%general%simControl,               &
        & proc    = params%general%proc,                     &
        & dt      = element_container%cubes                  &
        &                 %scheme_list(tree%global%minlevel) &
        &                 %time                              &
        &                 %dt,                               &
        & outUnit = logUnit(1)                               )
     ! call tem_stopTimer( timerHandle = atl_timerHandles%syncUpdate)

      ! Sync times again.
      element_container%cubes%statedata_list(:)%local_time%sim &
        & = params%general%simControl%now%sim


      ! !!! Output related !!!
      call atl_writeRestartIfNecessary(           &
        & restart    = params%general%restart,    &
        & simControl = params%general%simControl, &
        & equation   = equation,                  &
        & tree       = tree,                      &
        & mesh       = element_container%cubes    )

      if (params%plySampleTrack%tracking%control%active) then
        call ply_sampled_track_output(                &
          &    me         = params%plySampleTrack,    &
          &    mesh       = tree,                     &
          &    bc         = params%boundary,          &
          &    solver     = params%general%solver,    &
          &    proc       = params%general%proc,      &
          &    varSys     = equation%varSys,          &
          &    var_degree = params%var_degree,        &
          &    lvl_degree = params%lvl_degree,        &
          &    var_space  = params%var_space,         &
          &    simControl = params%general%simControl )
      end if


      if ( tem_status_run_terminate(params%general%simControl &
        &                                         %status)    ) then
        exit
      end if
      if (tem_status_run_end(params%general%simControl%status)) exit

      if (params%general%simControl%status%bits(tem_stat_interval)) then
        call tem_trackmem( params%general%solver%trackmem_file, &
          &                params%general%SimControl%now%iter   )
        call tem_horizontalSpacer(funit = logUnit(1))
      end if

      ! create the necessary info for the timestep control
      call tem_startTimer( timerHandle = atl_timerHandles%TimeStepInfo)
      do iList = tree%global%minlevel, tree%global%maxLevel
        call atl_initTimeStepInfo(                                      &
          & equation   = equation,                                      &
          & mesh       = element_container%cubes%mesh_list(iList),      &
          & poly_proj  = poly_proj_list(element_container%cubes         &
          &                                    %poly_proj_pos(iList)),  &
          & timestep   = element_container%time%elementSteps(iList),    &
          & control    = element_container%time%control,                &
          & statedata  = element_container%cubes%statedata_list(iList), &
          & kerneldata = element_container%cubes%kerneldata_list(iList) )
      end do
      call tem_stopTimer( timerHandle = atl_timerHandles%TimeStepInfo)
      ! Calculate the global timestep for next timestep of a whole part of a
      ! cubic mesh by a CFL condition based on the state of the current timestep
      ! includingcfl calcualtion, mpi reduce to get minimal timestep,
      ! checking for NAN timestep and checking for precice timestep if
      ! it is avaiable
      if ( .not. precice_available ) then
        call tem_startTimer( timerHandle = atl_timerHandles%get_timestep)
        call atl_get_timestep(                                          &
          & tree              = tree,                                   &
          & mesh_list         = element_container%cubes%mesh_list,      &
          & scheme_list       = element_container%cubes%scheme_list,    &
          & material_list     = element_container%cubes%material_list,  &
          & equation          = equation,                               &
          & time              = element_container%time,                 &
          & statedata_list    = element_container%cubes%statedata_list, &
          & nCellsNoBnd       = nCellsNoBnd,                            &
          & general           = params%general,                         &
          & adaptive_timestep = equation%adaptive_timestep,             &
          & initial           = .false.                                        )
        call tem_stopTimer( timerHandle = atl_timerHandles%get_timestep)
      else

        ! Now the coupling tool precice is doing its work in the advance step
        ! input is the current timestep ( chosen the one on minLevel, since up to
        ! now all levels ), the input is overwritten in advance and output is the
        ! timestep steered by preCICE
        ! add the advance into if_precice statement tp avoid a extra if
        ! statement for the advance and for the timestep calculation
        ! Start preciceTimer
        if (precice_handle%use_RK2_inter) then
          call tem_startTimer( timerHandle = atl_timerHandles%preciceAdv )
          precice_dt = element_container%cubes                             &
            &                           %scheme_list(tree%global%minlevel) &
            &                                                   %time%dt
          call tem_precice_advance(precice_dt)
          ! Stop preciceTimer
          write(*,*) "PreCICE dt after call =", precice_dt
          call tem_stopTimer( timerHandle = atl_timerHandles%preciceAdv )
          call tem_startTimer( timerHandle = atl_timerHandles%get_timestep)
          call atl_get_timestep(                                          &
            & tree              = tree,                                   &
            & mesh_list         = element_container%cubes%mesh_list,      &
            & scheme_list       = element_container%cubes%scheme_list,    &
            & material_list     = element_container%cubes%material_list,  &
            & equation          = equation,                               &
            & time              = element_container%time,                 &
            & statedata_list    = element_container%cubes%statedata_list, &
            & nCellsNoBnd       = nCellsNoBnd,                            &
            & general           = params%general,                         &
            & adaptive_timestep = equation%adaptive_timestep,             &
            & initial           = .false.,                                &
            & precice_dt        = precice_dt                              )

          precice_dt = precice_dt *2
          call tem_precice_advance(precice_dt)
          ! Stop preciceTimer
          write(*,*) "PreCICE dt after call second timestep =", precice_dt
          call atl_get_timestep(                                          &
            & tree              = tree,                                   &
            & mesh_list         = element_container%cubes%mesh_list,      &
            & scheme_list       = element_container%cubes%scheme_list,    &
            & material_list     = element_container%cubes%material_list,  &
            & equation          = equation,                               &
            & time              = element_container%time,                 &
            & statedata_list    = element_container%cubes%statedata_list, &
            & nCellsNoBnd       = nCellsNoBnd,                            &
            & general           = params%general,                         &
            & adaptive_timestep = equation%adaptive_timestep,             &
            & initial           = .false.,                                &
            & precice_dt        = precice_dt                              )
          call tem_stopTimer( timerHandle = atl_timerHandles%get_timestep)
        else
          call tem_startTimer( timerHandle = atl_timerHandles%preciceAdv )
          precice_dt = element_container%cubes                             &
            &                           %scheme_list(tree%global%minlevel) &
            &                                                   %time%dt
          call tem_precice_advance(precice_dt)
          ! Stop preciceTimer
          call tem_stopTimer( timerHandle = atl_timerHandles%preciceAdv )
          call tem_startTimer( timerHandle = atl_timerHandles%get_timestep)
          call atl_get_timestep(                                          &
            & tree              = tree,                                   &
            & mesh_list         = element_container%cubes%mesh_list,      &
            & scheme_list       = element_container%cubes%scheme_list,    &
            & material_list     = element_container%cubes%material_list,  &
            & equation          = equation,                               &
            & time              = element_container%time,                 &
            & statedata_list    = element_container%cubes%statedata_list, &
            & nCellsNoBnd       = nCellsNoBnd,                            &
            & general           = params%general,                         &
            & adaptive_timestep = equation%adaptive_timestep,             &
            & initial           = .false.,                                &
            & precice_dt        = precice_dt                              )
          call tem_stopTimer( timerHandle = atl_timerHandles%get_timestep)
        end if
      end if
      ! Checking if precice is finializing the simLoop
      call tem_precice_ongoing(is_Ongoing)
      if (is_Ongoing == 0) then
        write(logUnit(2),*) 'Simulation is finalizing due to preCICE'
        exit
      end if
      !!NE ! Dump and reset the timer after each Iteration
      !!NE call atl_dumpTimers( general = params%general,                      &
      !!NE   &                  nElems  = tree%global%nElems,                  &
      !!NE   &                  nDofs   = element_container%cubes              &
      !!NE   &                              %scheme_list(tree%global%minLevel) &
      !!NE   &                              %nDofs,                            &
      !!NE   &                  nVars   = equation%varSys%nScalars             )
      !!NE call atl_resetTimers()

    end do time_iter
    write(logUnit(2),"(A)") 'Finished main loop.'

    call tem_stopTimer( timerHandle = atl_timerHandles%simLoop )

    call tem_horizontalSpacer(funit = logUnit(6))

  end subroutine atl_solve_program
  ! ****************************************************************************


  ! ****************************************************************************
  subroutine atl_finalize_program( params, equation, tree, nCellsNoBnd, &
    &                              element_container )
    ! --------------------------------------------------------------------------
    ! The structure that holds the solver parameter
    type(atl_solver_param_type), intent(inout) :: params
    ! Description of the equation system to solve
    type(atl_equations_type), intent(inout) :: equation
    ! The treelmesh data structure
    type(treelmesh_type), intent(inout)        :: tree
    ! Number of cells on each level
    integer, allocatable, intent(inout)        :: nCellsNoBnd(:)
    ! Main data structure of Ateles describing the mesh elements
    type(atl_element_container_type), intent(inout) :: element_container
    ! --------------------------------------------------------------------------
    integer :: iConf
    ! --------------------------------------------------------------------------
    write(logUnit(1),'(/A)') 'Done with simulation at:'
    call tem_simControl_dump_now(params%general%simControl, logUnit(1))
    if ( tem_status_run_terminate(params%general%simControl%status)    ) then
      write(logUnit(1),'(/A)') 'Abnormal termination with the following status:'
      call tem_status_dump(params%general%simControl%status, logUnit(1))
    else
      write(logUnit(1),'(/A)') 'SUCCESSFUL run!'
      write(logUnit(2),*) 'Status of the simulation now:'
      call tem_status_dump(params%general%simControl%status, logUnit(2))
    end if
    call tem_horizontalSpacer(funit = logUnit(1))

    deallocate(nCellsNoBnd)
    call atl_deallocate_sourceData(element_container%cubes%source)

    ! Output final simulation result after leaving the time iteration loop.
    call atl_writeRestartIfNecessary(           &
      & restart    = params%general%restart,    &
      & simControl = params%general%simControl, &
      & equation   = equation,                  &
      & tree       = tree,                      &
      & mesh       = element_container%cubes,   &
      & force      = .true.                     )

    ! Write out all timings for performance analysis.
    if ( .not. tem_status_run_terminate(params%general%simControl &
      &                                                %status) ) then
      call atl_dumptimers(                                  &
        & general = params%general,                         &
        & nElems  = tree%global%nElems,                     &
        & nDofs   = element_container%cubes                 &
        &                %scheme_list(tree%global%minLevel) &
        &                %nDofs,                            &
        & nVars   = equation%varSys%nScalars                )

      if ( .not. trim(tree%write_weights) == '') then
        call atl_dumpWeights( tree )
      end if

    end if

    do iConf=1,size(params%general%solver%conf)
      call close_config(params%general%solver%conf(iConf))
    end do
    call tem_restart_finalize(params%general%restart)
    call tem_tracking_finalize(params%plySampleTrack%tracking)

    ! finialize environment
    if ( precice_available ) call tem_precice_finalize()

    call free_treelmesh(tree)

  end subroutine atl_finalize_program
  ! ****************************************************************************


end module atl_program_module
