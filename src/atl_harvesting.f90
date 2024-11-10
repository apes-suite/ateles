! Copyright (c) 2015-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

!> ATELES Postprocessing
!!
!! This is the harvesting tool for Ateles, enabling the post-processing
!! simulation data from the solver.
!! (c) 2015 University of Siegen
!!
program atl_harvesting

  use mpi

  use aotus_module,             only: close_config

  use treelmesh_module,         only: treelmesh_type
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_general_module,       only: tem_start, tem_finalize
  use tem_logging_module,       only: logUnit
  use tem_tracking_module,      only: tem_tracking_finalize, &
    &                                 tem_tracking_print_last_VTK_files
  use tem_restart_module,       only: tem_restart_finalize
  use tem_timer_module,         only: tem_startTimer, tem_stopTimer

  use atl_aux_module,           only: atl_banner
  use atl_varSys_module,        only: atl_varSys_solverData_type
  use atl_restart_module,       only: atl_writeRestart
  use atl_initialize_module,    only: atl_initialize
  use atl_container_module,     only: atl_element_container_type
  use atl_equation_module,      only: atl_equations_type
  use atl_kerneldata_module,    only: atl_kerneldata_update_estimates
  use atl_solver_param_module,  only: atl_solver_param_type
  use atl_scheme_module,        only: atl_modg_scheme_prp, &
    &                                 atl_modg_2d_scheme_prp, &
    &                                 atl_modg_1d_scheme_prp
  use atl_stabilize_module,     only: atl_stabilize
  use atl_program_module,       only: atl_load_config
  use atl_cube_elem_module,     only: atl_get_numberOfElemsPerLevel
  use atl_timer_module,         only: atl_addTimers, atl_timerHandles
  use ply_poly_project_module,  only: ply_poly_project_type
  use ply_sampled_tracking_module, only: ply_sampled_track_output

  use aotus_module,       only: flu_State, aot_get_val

  use aot_table_module,   only: aot_table_open,   &
    &                           aot_table_close,  &
    &                           aot_table_set_val

  implicit none

  ! ------------------------------------------------------------------------ !
  ! The treelmesh data structure
  type(treelmesh_type) :: tree

  ! The structure that holds the solver parameter
  type(atl_solver_param_type) :: params

  ! Main data structure of Ateles describing the mesh elements
  type(atl_element_container_type) :: element_container

  ! Description of the equation system to solve
  type(atl_equations_type) :: equation

  ! data in the variable system
  type(atl_varSys_solverData_type) :: varSys_data

  ! Number of cells on each levels
  integer, allocatable :: nCellsNoBnd(:)

  integer :: iConf, iLevel
  integer :: iError
  integer :: firstlevel, lastLevel
  integer :: scheme_table, spatial_table
  logical :: use_post_filter

  type(flu_State) :: conf
  type(ply_poly_project_type), allocatable :: poly_proj_list(:)
  ! ------------------------------------------------------------------------ !

  ! Initialize the treelm environment
  call tem_start( codeName   = 'ATL_HARVESTING',         &
    &             version    = 'v0.1',                   &
    &             general    = params%general,           &
    &             simcontrol = params%general%simControl )

  call atl_banner(trim(params%general%solver%version))

  call atl_addTimers()

  call tem_startTimer( timerHandle = atl_timerHandles%init )

  ! !!! Initialize !!!
  !
  ! Get the configuration file as first parameter from the command line,
  ! if none is specified, assume the configuration to be located in the
  ! current working directory and named harvester.lua.
  call atl_load_config( params     = params,         &
    &                   tree       = tree,           &
    &                   configFile = 'harvester.lua' )

  ! Set the polynomial space to Q-space in the configuration file.
  ! For the case of a simulation with P-polynomials
  ! this will automatically map the P-polynomial data to Q-space.
  ! We should stick to this behaviour till we implemented
  ! a proper way of dealing with P-polynomials in
  ! adaptive subsampling.
  conf = params%general%solver%conf(1)
  call aot_table_open( L       = conf,         &
    &                  thandle = scheme_table, &
    &                  key     = 'scheme'      )
  call aot_table_open( L       = conf,          &
    &                  parent  = scheme_table,  &
    &                  thandle = spatial_table, &
    &                  key     = 'spatial'      )
  call aot_table_set_val( L       = conf,          &
    &                     thandle = spatial_table, &
    &                     key     = 'modg_space',  &
    &                     val     = 'Q'            )
  call aot_table_close( L       = conf,         &
    &                   thandle = spatial_table )
  call aot_table_close( L       = conf,        &
    &                   thandle = scheme_table )
  call aot_get_val( L       = conf,              &
    &               key     = 'use_post_filter', &
    &               val     = use_post_filter,   &
    &               default = .false.,           &
    &               ErrCode = iError             )

  call atl_initialize( params         = params,            &
    &                  container      = element_container, &
    &                  equation       = equation,          &
    &                  tree           = tree,              &
    &                  varSys_data    = varSys_data,       &
    &                  poly_proj_list = poly_proj_list     )

  call atl_get_numberOfElemsPerLevel(                                &
    &  descriptor = element_container%cubes%mesh_list(:)%descriptor, &
    &  nCells     = nCellsNoBnd,                                     &
    &  tree       = tree                                             )

  if (use_post_filter) then
    write(logUnit(3),*) 'Apply filters on loaded data'
    call atl_stabilize( minlevel            = tree%global%minlevel,       &
      &                 maxlevel            = tree%global%maxLevel,       &
      &                 statedata_list      = element_container           &
      &                                       %cubes                      &
      &                                       %statedata_list,            &
      &                 statedata_stab_list = element_container           &
      &                                       %cubes                      &
      &                                       %statedata_stab_list,       &
      &                 mesh_list           = element_container           &
      &                                       %cubes                      &
      &                                       %mesh_list,                 &
      &                 scheme_list         = element_container           &
      &                                       %cubes                      &
      &                                       %scheme_list,               &
      &                 poly_proj_pos       = element_container           &
      &                                       %cubes                      &
      &                                       %poly_proj_pos,             &
      &                 poly_proj_list      = poly_proj_list,             &
      &                 equation            = equation,                   &
      &                 tree                = tree,                       &
      &                 bc                  = element_container           &
      &                                       %cubes%bc,                  &
      &                 boundary            = element_container           &
      &                                       %cubes%boundary_list,       &
      &                 general             = params%general,             &
      &                 commStateTimer      = atl_timerHandles%commState, &
      &                 material_list       = element_container           &
      &                                       %cubes                      &
      &                                       %material_list              )
  end if

  write(logUnit(3),*) 'Obtain estimation for maximal deviation in each element'
  do iLevel = tree%global%minlevel, tree%global%maxLevel
    call atl_kerneldata_update_estimates(                               &
      &    statedata  = element_container%cubes%statedata_list(iLevel), &
      &    kerneldata = element_container%cubes%kerneldata_list(iLevel) )
  end do

  call tem_stopTimer( timerHandle = atl_timerHandles%init )

  write(logUnit(1),*) '... done with initialization!'
  call tem_horizontalSpacer(funit = logUnit(1))
  write(logUnit(1),*) 'Starting to process data ...'

  if ( params%plySampleTrack%tracking%control%active ) then

    allocate(params%var_degree(equation%varsys%varname%nVals))
    allocate(params%lvl_degree(element_container%cubes%maxlevel))
    allocate(params%var_space(equation%varsys%varname%nVals))
    firstLevel = lbound(element_container%cubes%scheme_list,1)
    lastLevel  = ubound(element_container%cubes%scheme_list,1)

    params%lvl_degree = 0

    select case(element_container%cubes%scheme_list(firstLevel)%scheme)
    case (atl_modg_scheme_prp)
      params%var_degree = maxval( element_container     &
        &                         %cubes%scheme_list(:) &
        &                         %modg%maxPolyDegree   )
      do iLevel = firstLevel, lastLevel
        params%lvl_degree(iLevel) = element_container          &
          &                         %cubes%scheme_list(iLevel) &
          &                         %modg%maxPolyDegree
      end do
      params%var_space = element_container              &
        &                %cubes%scheme_list(firstLevel) &
        &                %modg%basisType

    case (atl_modg_2d_scheme_prp)
      params%var_degree = maxval( element_container      &
        &                         %cubes%scheme_list(:)  &
        &                         %modg_2d%maxPolyDegree )
      do iLevel = firstLevel, lastLevel
        params%lvl_degree(iLevel) = element_container          &
          &                         %cubes%scheme_list(iLevel) &
          &                         %modg_2d%maxPolyDegree
      end do
      params%var_space = element_container              &
        &                %cubes%scheme_list(firstLevel) &
        &                %modg_2d%basisType

    case (atl_modg_1d_scheme_prp)
      params%var_degree = maxval( element_container      &
        &                         %cubes%scheme_list(:)  &
        &                         %modg_1d%maxPolyDegree )
      do iLevel = firstLevel, lastLevel
        params%lvl_degree(iLevel) = element_container          &
          &                         %cubes%scheme_list(iLevel) &
          &                         %modg_1d%maxPolyDegree
      end do
      params%var_space = element_container              &
        &                %cubes%scheme_list(firstLevel) &
        &                %modg_1d%basisType

    end select

    ! All data processing done via trackings now.
    ! Leaving out the simControl should ensure, that all trackings
    ! will be processed at this point.
    call ply_sampled_track_output(                 &
      & me         = params%plySampleTrack,        &
      & mesh       = tree,                         &
      & bc         = params%boundary,              &
      & solver     = params%general%solver,        &
      & proc       = params%general%proc,          &
      & varSys     = equation%varSys,              &
      & var_degree = params%var_degree,            &
      & lvl_degree = params%lvl_degree,            &
      & var_space  = params%var_space,             &
      & time       = params%general%simControl%now )

    call tem_tracking_print_last_VTK_files(params%plySampleTrack%tracking)

  end if

  ! Write out the initial data if restart is active.
  if (params%general%restart%controller%writeRestart) then

    ! This is only useful, if we did not read the restart data in the first
    ! place...
    if (.not. params%general%restart%controller%readRestart) then
      call atl_writeRestart(                                       &
        &    mesh           = element_container%cubes,             &
        &    restart        = params%general%restart,              &
        &    tree           = tree,                                &
        &    equation       = equation,                            &
        &    timing         = params%general%simControl%now,       &
        &    varSys         = equation%varSys,                     &
        &    timerHandle    = atl_timerHandles%wRestart            )
    end if

  end if

  call tem_horizontalSpacer(funit = logUnit(1))
  write(logUnit(1),*) ''
  write(logUnit(1),*) 'Done with Harvesting Ateles data!'
  write(logUnit(1),*) ''
  call tem_horizontalSpacer(funit = logUnit(1))

  ! !!! Finalization !!!
  do iConf=1,size(params%general%solver%conf)
    call close_config(params%general%solver%conf(iConf))
  end do
  call tem_restart_finalize(params%general%restart)
  call tem_tracking_finalize(params%plySampleTrack%tracking)

  ! finialize environment
  call tem_finalize(params%general)

end program atl_harvesting
