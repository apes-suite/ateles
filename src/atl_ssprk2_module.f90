! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017-2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> author: Jens Zudrop
!! Routines,functions, datatypes related to
!! second order, strong stability preserving explicit Runge Kutta method.
!!
!! This module implements all routines and datatypes that are
!! necessary for a second order, strong stability preserving
!! Runge-Kutta method.
!! The implementation is based on section 4.1 of:
!! Strong stability-preserving High-Order Time-Discretization Methods,
!! Gottlieb S. , Shu Chi-Wang,  Tadmor Etian ,
!! SIAM Review, Vol. 43, No. 1, pp. 89 - 112
module atl_ssprk2_module
  use env_module,                            only: rk
  use tem_aux_module,                        only: tem_abort
  use tem_element_module,                    only: eT_fluid
  use tem_general_module,                    only: tem_general_type
  use tem_logging_module,                    only: logUnit
  use treelmesh_module,                      only: treelmesh_type
  use tem_timer_module,                      only: tem_startTimer, tem_stopTimer
  use tem_precice_module,                    only: precice_available, &
   &                                               precice_handle

  use ply_poly_project_module,               only: ply_poly_project_type, &
    &                                              assignment(=)

  use atl_compute_module,                    only: atl_compute_rhs,    &
    &                                              atl_preprocess_rhs, &
    &                                              atl_postprocess_rhs
  use atl_elemental_time_integration_module, only: atl_timestep_type
  use atl_cube_elem_module,                  only: atl_cube_elem_type
  use atl_kerneldata_module,                 only: atl_kerneldata_type, &
    &                                              atl_statedata_type
  use atl_scheme_module,                     only: atl_scheme_type,         &
    &                                              atl_local_timestep_type
  use atl_source_types_module,               only: atl_source_type
  use atl_boundary_module,                   only: atl_level_boundary_type
  use atl_facedata_module,                   only: atl_facedata_type
  use atl_time_integration_module,           only: atl_global_timestep_type
  use atl_equation_module,                   only: atl_equations_type
  use atl_bc_header_module,                  only: atl_boundary_type
  use atl_stabilize_module,                  only: atl_stabilize
  use atl_materialPrp_module,                only: atl_material_type
  use atl_writePrecice_module,               only: atl_write_precice
  use atl_cube_container_module,             only: atl_cube_container_type
  use atl_penalization_module,               only: atl_penalizationData_type
  use atl_timer_module,                      only: atl_timerHandles

  implicit none
  private

  public :: atl_init_explicit_ssprk

contains


  ! ****************************************************************************
  !> Routine to initialize explicit runge kutta scheme for timestepping.
  subroutine atl_init_explicit_ssprk( me, minLevel, maxLevel, steps, &
    &                                 statedata_list                 )
    ! --------------------------------------------------------------------------
    !> The datatype to initialize.
    type(atl_global_timestep_type), intent(inout) :: me

    !> The minimum of level of the mesh.
    integer, intent(in) :: minLevel

    !> The maximum of level of the mesh.
    integer, intent(in) :: maxLevel

    !> The number of steps in the runge kutta procedure
    integer, intent(in) :: steps

    !> The state list used in your solver, for each level one entry.
    type(atl_statedata_type) , intent(in) :: statedata_list(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iStep
    ! --------------------------------------------------------------------------
    ! point to euler step function, we dont have to allocate any of the
    ! additional arrays for coefficients or buffering.
    allocate( me%elementSteps(minLevel:maxLevel) )
    select case(steps)
    case(2)
      do iLevel = minLevel, maxLevel
        if(allocated(statedata_list(iLevel)%state)) then
          ! We store one integer, indicating the elemental operation in which
          ! step we are.
          allocate(me%elementSteps(iLevel)%timestepInfoInteger(1))
          allocate(me%elementSteps(iLevel)%timestepData(2))
          do iStep = 1, 2
            allocate(me%elementSteps(iLevel)%timestepData(iStep) &
              &         %state(size(statedata_list(iLevel)%state, 1), &
              &                size(statedata_list(iLevel)%state, 2), &
              &                size(statedata_list(iLevel)%state, 3)) &
              &     )
          end do
          me%elementSteps(iLevel)%elemStep => elemental_timestep_ssprk2
          me%elementSteps(iLevel)%elemStep_vec => elemental_timestep_vec_ssprk2
          me%elementSteps(iLevel)%updateStep => update_timestep_ssprk2
        end if
      end do
      me%meshStep => mesh_timestep_ssprk2

    case default
      write(logUnit(1),*) 'init SSP Runge Kutta method for this number of steps not' &
        &            //' implemented, stopping...'
      call tem_abort()
    end select
  end subroutine atl_init_explicit_ssprk
  ! ****************************************************************************


  ! ****************************************************************************
  !> Levelwise updating of runge kutta of order 2
  subroutine update_timestep_ssprk2( me, timestepInfo )
    ! ------------------------------------------------------------------------
    !> The type of your timestepping.
    class(atl_timestep_type), intent(inout) :: me

    !> Local timestepping information for that part of the mesh
    type(atl_local_timestep_type), intent(in) :: timestepInfo
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! nothing to be done here.
  end subroutine update_timestep_ssprk2
  ! ****************************************************************************


  ! ****************************************************************************
  !> Elemental operation for timestepping of order 2.
  subroutine elemental_timestep_ssprk2( me, state, cell, dof, sideFlux)
    ! ------------------------------------------------------------------------
    !> Description of the timestep integration method.
    class(atl_timestep_type), intent(inout) :: me

    !> The state of all cells on this level. This field will be updated
    !! ad the cell position. See kerneldata type for more explanations.
    real(kind=rk), intent(inout) :: state(:,:,:)

    !> Position of the cell to update in the state vector.
    integer,intent(in) :: cell

    !> The degree of freedom to update
    integer, intent(in) :: dof

    !> The flux for one of the sides of this cell. The length of this array
    !! is the number of conservative variables of your equation.
    real(kind=rk), intent(in) :: sideFlux(:)
    ! ------------------------------------------------------------------------
    me%timestepData(me%timestepInfoInteger(1))%state(cell,dof,:) &
      & = me%timestepData(me%timestepInfoInteger(1))%state(cell,dof,:) &
      & + sideFlux(:)
  end subroutine elemental_timestep_ssprk2
  ! ****************************************************************************


  ! ****************************************************************************
  !> Elemental operation for timestepping of order 2.
  subroutine elemental_timestep_vec_ssprk2( me, state, kerneldata)
    ! ------------------------------------------------------------------------
    !> Description of the timestep integration method.
    class(atl_timestep_type), intent(inout) :: me

    !> The state of all cells on this level. This field will be updated
    !! ad the cell position. See kerneldata type for more explanations.
    real(kind=rk), intent(inout) :: state(:,:,:)

    !> Complete kerneldata to get the flux from with additional information.
    type(atl_kerneldata_type), intent(in) :: kerneldata
    ! ------------------------------------------------------------------------

    call compute_vec(kerneldata%nTotal, &
      &              kerneldata%nDofs, &
      &              kerneldata%nVars,  &
      &              me%timestepData(me%timestepInfoInteger(1))%state, &
      &              kerneldata%state_der)


  end subroutine elemental_timestep_vec_ssprk2
  ! ****************************************************************************


  ! ****************************************************************************
  ! put it in an own subroutine to allow vectorization with gcc
  subroutine compute_vec(nTotal, nDofs, nScalars, &
      &                  state, state_der         )
    ! ------------------------------------------------------------------------
    integer, intent(in) :: nTotal, nDofs, nScalars
    real(kind=rk), intent(inout) :: state(nTotal, nDofs, nScalars)
    real(kind=rk), intent(in) :: state_der(:,:,:)
    ! ------------------------------------------------------------------------
    state(:nTotal,:,:) = state(:nTotal,:,:) &
      & + state_der(:nTotal,:nDofs,:)

  end subroutine compute_vec
  ! ****************************************************************************


  ! ****************************************************************************
  !> Subroutine for timestepping with explicit runge kutta of order 2.
  subroutine mesh_timestep_ssprk2(minLevel, maxLevel, currentLevel,cubes, &
    &                tree, timestep_list, nSteps, equation, general,      &
    &                commStateTimer, poly_proj_list                       )
    ! --------------------------------------------------------------------------
    !> The minimum refinement level of the mesh.
    integer, intent(in)           :: minLevel

    !> The maximum refinement level of the mesh.
    integer, intent(in)           :: maxLevel

    !> The level the timestep has to be performed for.
    integer, intent(in)           :: currentLevel

    !> Container for the cubical elements.
    type(atl_cube_container_type), intent(inout) :: cubes

    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree

    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep_list( &
      &                       minLevel:)

    !> The number of steps of the time stepping scheme (assumed to be 2)
    integer, intent(in) :: nSteps

    !> The equation you are operating with.
    type(atl_equations_type),intent(inout) :: equation

    !> General treelm settings
    type(tem_general_type), intent(inout) :: general

    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer

    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iStep
    type(atl_statedata_type) :: statedata_list_temp(minLevel:maxLevel)
    ! --------------------------------------------------------------------------

    ! Init the arrays with the intermediate results
    do iLevel = minLevel, maxLevel
      do iStep = 1, 2
        timestep_list(iLevel)%timestepData(iStep)%state(:,:,:) = 0.0_rk
      end do
      allocate( statedata_list_temp(iLevel)        &
        &         %state(cubes%mesh_list(iLevel)   &
        &                     %descriptor          &
        &                     %elem                &
        &                     %nElems(eT_fluid),   &
        &                cubes%scheme_list(iLevel) &
        &                     %nDofs,              &
        &                equation%varSys           &
        &                        %nScalars)        )
      statedata_list_temp(iLevel)%state = 0.0_rk
    end do


    ! 1st runge kutta substep. Will call itself on the finer
    ! levels recursively.
    do iLevel = minLevel, maxLevel
      ! Store, that we are in the first step of the RK method.
      timestep_list(iLevel)%timestepInfoInteger(1) = 1
    end do
    call rk2_substep( minLevel              = minLevel,                    &
      &               maxLevel              = maxLevel,                    &
      &               currentLevel          = currentLevel,                &
      &               mesh_list             = cubes%mesh_list,             &
      &               tree                  = tree,                        &
      &               levelPointer          = cubes%levelPointer,          &
      &               kerneldata_list       = cubes%kerneldata_list,       &
      &               statedata_list        = cubes%statedata_list,        &
      &               facedata_list         = cubes%facedata_list,         &
      &               source                = cubes%source,                &
      &               penalizationdata_list = cubes%penalizationdata_list, &
      &               boundary_list         = cubes%boundary_list,         &
      &               bc                    = cubes%bc,                    &
      &               scheme_list           = cubes%scheme_list,           &
      &               poly_proj_pos         = cubes%poly_proj_pos,         &
      &               poly_proj_list        = poly_proj_list,              &
      &               timestep_list         = timestep_list,               &
      &               equation              = equation,                    &
      &               material_list         = cubes%material_list,         &
      &               general               = general,                     &
      &               commStateTimer        = commStateTimer               )
    ! Now, build the intermediate values for the RK2 method. Since the substep
    ! was carried out for all the levels in the previous step, we can continue
    ! in a simple loop over the levels.
    do iLevel = minLevel, maxLevel
      call compute_intermediate(                               &
        &    nTotal    = cubes%mesh_list(iLevel)               &
        &                     %descriptor                      &
        &                     %elem                            &
        &                     %nElems(eT_fluid),               &
        &    nDofs     = cubes%scheme_list(iLevel)             &
        &                     %nDofs,                          &
        &    nScalars  = equation%varSys                       &
        &                        %nScalars,                    &
        &    state_tmp = statedata_list_temp(iLevel)           &
        &                  %state,                             &
        &    state1    = cubes%statedata_list(iLevel)          &
        &                     %state,                          &
        &    state2    = timestep_list(iLevel)%timestepData(1) &
        &                                     %state,          &
        &    factor    = cubes%scheme_list(iLevel)             &
        &                     %time                            &
        &                     %dt                              )
    end do
    ! Stabilize the intermediate result of this stage of the RK2 scheme
    call atl_stabilize(                                     &
      &    minlevel            = minlevel,                  &
      &    maxlevel            = maxlevel,                  &
      &    statedata_list      = statedata_list_temp,       &
      &    statedata_stab_list = cubes%statedata_stab_list, &
      &    mesh_list           = cubes%mesh_list,           &
      &    scheme_list         = cubes%scheme_list,         &
      &    equation            = equation,                  &
      &    tree                = tree,                      &
      &    poly_proj_pos       = cubes%poly_proj_pos,       &
      &    poly_proj_list      = poly_proj_list,            &
      &    bc                  = cubes%bc,                  &
      &    boundary            = cubes%boundary_stab_list,  &
      &    general             = general,                   &
      &    material_list       = cubes%material_list,       &
      &    commStateTimer      = commStateTimer             )

    !if (precice_available) then
    if (precice_handle%use_RK2_inter) then
      ! After the runge kutta FULL step, write the results to precice
      ! Start precice timer with barrier
      call tem_startTimer( timerHandle = atl_timerHandles%preciceWrite )
        call atl_write_precice( stFunList = equation%stFunList,     &
          &                     varSys    = equation%varSys,        &
          &                     time      = general%simControl%now, &
          &                     tree      = tree                    )
      ! Stop precice timer with barrier
      call tem_stopTimer( timerHandle = atl_timerHandles%preciceWrite )
    end if
    !end if

    ! 2nd runge kutta substep. Will call itself on the finer
    ! levels recursively. Update again to the new timepoint of
    ! the next RK substep.
    do iLevel = minLevel, maxLevel
      cubes%statedata_list(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim &
        &  + cubes%scheme_list(iLevel)%time%dt
      ! update the intermediate time for next rk2 substep, this necessary to calculate
      ! source terms and time-dependent boundary values correctly.
      statedata_list_temp(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim
      ! indicates that we are in the second step of RK method
      timestep_list(iLevel)%timestepInfoInteger(1) = 2
    end do
    ! Just, the final substep, no need to store an intermediate
    ! state here.
    call rk2_substep( minLevel              = minLevel,                    &
      &               maxLevel              = maxLevel,                    &
      &               currentLevel          = currentLevel,                &
      &               mesh_list             = cubes%mesh_list,             &
      &               tree                  = tree,                        &
      &               levelPointer          = cubes%levelPointer,          &
      &               kerneldata_list       = cubes%kerneldata_list,       &
      &               statedata_list        = statedata_list_temp,         &
      &               facedata_list         = cubes%facedata_list,         &
      &               source                = cubes%source,                &
      &               penalizationdata_list = cubes%penalizationdata_list, &
      &               boundary_list         = cubes%boundary_list,         &
      &               bc                    = cubes%bc,                    &
      &               scheme_list           = cubes%scheme_list,           &
      &               poly_proj_pos         = cubes%poly_proj_pos,         &
      &               poly_proj_list        = poly_proj_list,              &
      &               timestep_list         = timestep_list,               &
      &               equation              = equation,                    &
      &               material_list         = cubes%material_list,         &
      &               general               = general,                     &
      &               commStateTimer        = commStateTimer               )


    ! No stabilizationwe stabilize after the final update step


    ! The final update step. No need to have a recusrive structure here,
    ! so, we can do it in parallel for the different levels.
    do iLevel = minLevel, maxLevel
      call rk2_update( cubes%statedata_list(iLevel), &
        &              cubes%scheme_list(iLevel), timestep_list(iLevel) )
    end do

    ! Stabilize the final result
    call atl_stabilize( minlevel = minlevel, maxlevel = maxlevel, &
                      & statedata_list = cubes%statedata_list,    &
                      & statedata_stab_list = cubes%statedata_stab_list,&
                      & mesh_list = cubes%mesh_list,              &
                      & scheme_list = cubes%scheme_list,          &
                      & equation = equation,                      &
                      & tree = tree,                              &
                      & poly_proj_pos= cubes%poly_proj_pos,       &
                      & poly_proj_list= poly_proj_list,           &
                      & bc = cubes%bc,                            &
                      & boundary = cubes%boundary_stab_list,      &
                      & general = general,                        &
                      & material_list= cubes%material_list,       &
                      & commStateTimer = commStateTimer           )


    if (precice_available) then
    ! After the runge kutta FULL step, write the results to precice
    ! Start precice timer with barrier
    call tem_startTimer( timerHandle = atl_timerHandles%preciceWrite )
      call atl_write_precice( stFunList = equation%stFunList,     &
        &                     varSys    = equation%varSys,        &
        &                     time      = general%simControl%now, &
        &                     tree      = tree                    )
    ! Stop precice timer with barrier
    call tem_stopTimer( timerHandle = atl_timerHandles%preciceWrite )
    end if

  end subroutine mesh_timestep_ssprk2
  ! ****************************************************************************



  ! ****************************************************************************
  ! own subroutine to allow vectorization with gcc
  subroutine compute_intermediate(nTotal, nDofs, nScalars, &
      &                           state_tmp, state1, state2, &
      &                           factor)
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nTotal, nDofs, nScalars
    real(kind=rk), intent(out) :: state_tmp(nTotal, nDofs, nScalars)
    real(kind=rk), intent(in) :: state1(nTotal, nDofs, nScalars)
    real(kind=rk), intent(in) :: state2(nTotal, nDofs, nScalars)
    real(kind=rk), intent(in) :: factor
    ! ---------------------------------------------------------------------------

    state_tmp(:,:,:) = state1(:,:,:) + factor*state2(:,:,:)

  end subroutine compute_intermediate
  ! ****************************************************************************


  ! ****************************************************************************
  !> Subroutine calculates the final update step of the Runge-Kutta method.
  !! It is performing levelwise.
  subroutine rk2_update( statedata_list, scheme_list, timestep_list )
    ! ---------------------------------------------------------------------------
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list
    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep_list
    ! ---------------------------------------------------------------------------
    integer :: nTotal, nDofs, nScalars
    ! ---------------------------------------------------------------------------

    nTotal = size(statedata_list%state, 1)
    nDofs = size(statedata_list%state, 2)
    nScalars = size(statedata_list%state, 3)

    call compute_up(nTotal, nDofs, nScalars, &
      &          statedata_list%state, &
      &          timestep_list%timestepData(1)%state, &
      &          timestep_list%timestepData(2)%state, &
      &          scheme_list%time%dt )



  end subroutine rk2_update
  ! ****************************************************************************


  ! ****************************************************************************
  ! own subroutine to allow vectorization with gcc
  subroutine compute_up(nTotal, nDofs, nScalars, &
      &              state, state1, state2, &
      &              dt)
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nTotal, nDofs, nScalars
    real(kind=rk), intent(inout) :: state(nTotal, nDofs, nScalars)
    real(kind=rk), intent(in) :: state1(nTotal, nDofs, nScalars)
    real(kind=rk), intent(in) :: state2(nTotal, nDofs, nScalars)
    real(kind=rk), intent(in) :: dt
    ! ---------------------------------------------------------------------------

    state(:,:,:) = state(:,:,:) + 0.5_rk*dt*state1(:,:,:) + 0.5_rk*dt*state2(:,:,:)

  end subroutine compute_up
  ! ****************************************************************************


  ! ****************************************************************************
  !> Subroutine calculates a substep of the Runge-Kutta timestepping scheme.
  !! Calls itself recursively for the finer levels until the finest level is reached.
  recursive subroutine rk2_substep(minLevel, maxLevel, currentLevel,       &
    & mesh_list, tree, levelPointer, kerneldata_list, statedata_list,      &
    & facedata_list, source, penalizationdata_list, boundary_list, bc,     &
    & scheme_list, poly_proj_pos, poly_proj_list, timestep_list, equation, &
    & material_list, general, commStateTimer                               )
    ! ---------------------------------------------------------------------------
    !> The minimum refinement level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum refinement level of the mesh.
    integer, intent(in) :: maxLevel
    !> The level the timestep has to be performed for.
    integer, intent(in) :: currentLevel
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type), intent(inout) :: mesh_list(minLevel:maxLevel)
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> Pointer for elements from global treeID list index to index in
    !! levelwise fluid lists
    integer, intent(in)  :: levelPointer(:)
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: &
      & kerneldata_list(minLevel:maxLevel)
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list(minLevel:maxLevel)
    !> List of faces states you want to calc the rhs for. For each level we have one.
    type(atl_facedata_type), intent(inout) :: facedata_list(minLevel:maxLevel)
    !> List of sources, for each level
    type(atl_source_type), intent(inout) :: source
    !> List of penalization data, for each level.
    type(atl_penalizationData_type), intent(inout) &
      & :: penalizationdata_list(minLevel:maxLevel)
    !> List of boundaries, for each level.
    type(atl_level_boundary_type), intent(inout) :: &
      & boundary_list(minLevel:maxLevel)
    !> Global description of the boundaries
    type(atl_boundary_type), intent(in) :: bc(:)
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list(minLevel:maxLevel)
    !> List of position pointer of projection, for each level.
    integer, intent(in) :: poly_proj_pos(minLevel:maxLevel)
    !> List of projection methods.
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep_list(minLevel:maxLevel)
    !> The equation you are operating with.
    type(atl_equations_type),intent(inout) :: equation
    !> Material parameter description.
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxlevel)
    !> General treelm settings
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    ! Before going to the next finer level, we do the following:
    !(0. Create modal representation of source terms)
    ! 1. Project modal representation to the faces.
    ! 2. Communicate the face representations (only the current level).
    ! 3. Interpolate face representations from the current level to
    !    the next finer level.
    ! All these steps are done in the preprocess step of the compute module.
    call atl_preprocess_rhs( minLevel        = minLevel,        &
      &                      maxLevel        = maxLevel,        &
      &                      currentLevel    = currentLevel,    &
      &                      mesh_list       = mesh_list,       &
      &                      tree            = tree,            &
      &                      statedata_list  = statedata_list,  &
      &                      facedata_list   = facedata_list,   &
      &                      boundary_list   = boundary_list,   &
      &                      bc              = bc,              &
      &                      scheme_list     = scheme_list,     &
      &                      poly_proj_list  = poly_proj_list,  &
      &                      equation        = equation,        &
      &                      material_list   = material_list,   &
      &                      general         = general          )

    ! Now, we call the timestep for the next level. I.e. the timestepping
    ! for the finer levels is called recursively unitl we reached the
    ! finest level (i.e. maxLevel) of the mesh.
    if(currentLevel .lt. maxLevel) then
      call rk2_substep( minLevel              = minLevel,              &
        &               maxLevel              = maxLevel,              &
        &               currentLevel          = currentLevel + 1,      &
        &               mesh_list             = mesh_list,             &
        &               tree                  = tree,                  &
        &               levelPointer          = levelPointer,          &
        &               kerneldata_list       = kerneldata_list,       &
        &               statedata_list        = statedata_list,        &
        &               facedata_list         = facedata_list,         &
        &               source                = source,                &
        &               penalizationdata_list = penalizationdata_list, &
        &               boundary_list         = boundary_list,         &
        &               bc                    = bc,                    &
        &               scheme_list           = scheme_list,           &
        &               poly_proj_pos         = poly_proj_pos,         &
        &               poly_proj_list        = poly_proj_list,        &
        &               timestep_list         = timestep_list,         &
        &               equation              = equation,              &
        &               material_list         = material_list,         &
        &               general               = general,               &
        &               commStateTimer        = commStateTimer         )
    end if

    ! After making the timestep on the finer level, we do the following:
    ! 1. Interpolate flux from finer level to the current level.
    ! 2. Calculate the fluxes for the current level.
    ! 3. Communicate the fluxes (only the current level).
    ! 4. Calculate the physical fluxes.
    ! 5. Project physical flux + numerical fluxes to test functions.
    ! 6. Multiply with inverse of (cell local) mass matrix.
    ! ... this call includes steps 1 to 5.
    call atl_compute_rhs( minLevel, maxLevel, currentLevel, mesh_list, tree, &
      & kerneldata_list, statedata_list, facedata_list, source,              &
      & penalizationdata_list, scheme_list,poly_proj_pos,poly_proj_list,     &
      & equation, material_list, general                                     )
    ! ... this call applies the inverse of the mass matrix.
    call atl_postprocess_rhs( mesh       = mesh_list(currentLevel),       &
      &                       kerneldata = kerneldata_list(currentLevel), &
      &                       statedata  = statedata_list(currentLevel),  &
      &                       scheme     = scheme_list(currentLevel),     &
      &                       timestep   = timestep_list(currentLevel),   &
      &                       equation   = equation                       )

  end subroutine rk2_substep
  ! ****************************************************************************

end module atl_ssprk2_module
