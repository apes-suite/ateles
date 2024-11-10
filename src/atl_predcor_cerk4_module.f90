! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> time integration approach for local timestepping with
!! local predictor and global corrector
!> local predictor: Continuous Extension Runge Kutta with 4 steps (order 3)
!! global corrector: gauss-quadrature with 2 points in time (order 3)
module atl_predcor_cerk4_module

  use env_module,                             only: rk
  use tem_aux_module,                         only: tem_abort
  use tem_general_module,                     only: tem_general_type
  use tem_logging_module,                     only: logUnit
  use treelmesh_module,                       only: treelmesh_type
  use tem_element_module,                     only: eT_fluid

  use ply_poly_project_module,                only: ply_poly_project_type

  use atl_elemental_time_integration_module,  only: atl_timestep_type
  use atl_cube_elem_module,                   only: atl_cube_elem_type
  use atl_kerneldata_module,                  only: atl_kerneldata_type, &
    &                                               atl_statedata_type
  use atl_scheme_module,                      only: atl_scheme_type,        &
    &                                               atl_local_timestep_type
  use atl_boundary_module,                    only: atl_level_boundary_type
  use atl_source_types_module,                only: atl_source_type
  use atl_facedata_module,                    only: atl_facedata_type
  use atl_time_integration_module,            only: atl_global_timestep_type
  use atl_equation_module,                    only: atl_equations_type
  use atl_bc_header_module,                   only: atl_boundary_type
  use atl_compute_local_module,               only: atl_preprocess_local_rhs
  use atl_materialPrp_module,                 only: atl_material_type
  use atl_cube_container_module,              only: atl_cube_container_type
  use atl_penalization_module,                only: atl_penalizationData_type

  implicit none
  private

  public :: atl_init_explicitLocalPredictorGlobalCorrector

contains


  ! ****************************************************************************
  !> Routine to initialize explicit runge kutta scheme for timestepping.
  subroutine atl_init_explicitLocalPredictorGlobalCorrector( &
    &          me, minLevel, maxLevel, steps, statedata_list)
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
    case(4) ! we have continuous extension runge kutta predictor with 4 steps, and third order
      do iLevel = minLevel, maxLevel
        if(allocated(statedata_list(iLevel)%state)) then
          ! We store one integer, indicating the elemental operation in which
          ! step we are.
          allocate(me%elementSteps(iLevel)%timestepInfoInteger(1))
          ! In a 4 step runge kutta method we have to store 5 intermediate
          ! results.
          allocate(me%elementSteps(iLevel)%timestepData(4))
          do iStep = 1, 4
            allocate( me%elementSteps(iLevel) &
              &         %timestepData(iStep) &
              &         %state( size(statedata_list(iLevel)%state, 1),  &
              &                 size(statedata_list(iLevel)%state, 2),  &
              &                 size(statedata_list(iLevel)%state, 3) ) )
          end do
          me%elementSteps(iLevel)%elemStep => elemental_timestep_predcor_cerk4
          me%elementSteps(iLevel)%elemStep_vec => &
            & elemental_timestep_vec_predcor_cerk4
          me%elementSteps(iLevel)%updateStep => update_timestep_predcor_cerk4
        end if
      end do
      me%meshStep => mesh_timestep_predcor_cerk4

    case default
      write(logUnit(1),*) 'init cerk-predictor--corrector method for this number of steps not' &
        &            //' implemented, stopping...'
      call tem_abort()
    end select
  end subroutine atl_init_explicitLocalPredictorGlobalCorrector
  ! ****************************************************************************


  ! ****************************************************************************
  !> Levelwise updating of runge kutta of order 4
  subroutine update_timestep_predcor_cerk4( me, timestepInfo )
    ! ------------------------------------------------------------------------
    !> The type of your timestepping.
    class(atl_timestep_type), intent(inout) :: me

    !> Local timestepping information for that part of the mesh
    type(atl_local_timestep_type), intent(in) :: timestepInfo
    ! ------------------------------------------------------------------------

    ! nothing to be done here.

  end subroutine update_timestep_predcor_cerk4
  ! ****************************************************************************


  ! ****************************************************************************
  !> Elemental operation for timestepping of order 4.
  subroutine elemental_timestep_predcor_cerk4( me, state, cell, dof, sideFlux)
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

    ! nothing to be done here.

  end subroutine elemental_timestep_predcor_cerk4
  ! ****************************************************************************


  ! ****************************************************************************
  !> Elemental operation for timestepping of order 4.
  subroutine elemental_timestep_vec_predcor_cerk4( me, state, kerneldata)
    ! ------------------------------------------------------------------------
    !> Description of the timestep integration method.
    class(atl_timestep_type), intent(inout) :: me

    !> The state of all cells on this level. This field will be updated
    !! ad the cell position. See kerneldata type for more explanations.
    real(kind=rk), intent(inout) :: state(:,:,:)

    !> Complete kerneldata to get the flux from with additional information.
    type(atl_kerneldata_type), intent(in) :: kerneldata
    ! ------------------------------------------------------------------------

    ! empty, all work is done by the mesh_timestep_predcor_cerk4 routine,
    ! as we cannot access the 'original' state from here
    ! (the argument 'state' above is a temporary array for the predicted state)

  end subroutine elemental_timestep_vec_predcor_cerk4
  ! ****************************************************************************


  ! ****************************************************************************
  !> Subroutine for timestepping with explicit runge kutta of order 4.
  subroutine mesh_timestep_predcor_cerk4(minLevel, maxLevel, currentLevel,    &
    &                         cubes, tree, timestep_list, nSteps, equation,   &
    &                         general, commStateTimer, poly_proj_list         )
    ! --------------------------------------------------------------------------
    !> The minimum refinement level of the mesh.
    integer, intent(in) :: minLevel

    !> The maximum refinement level of the mesh.
    integer, intent(in) :: maxLevel

    !> The level the timestep has to be performed for.
    integer, intent(in) :: currentLevel

    !> Container for the cubical elements.
    type(atl_cube_container_type), intent(inout) :: cubes

    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree

    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep_list( minLevel: )

    !> The number of steps of the time stepping scheme (assumed to be 4)
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
    integer :: iStep
    type(atl_statedata_type) :: statedata_list_temp(minLevel:maxLevel)
    integer :: nElems
    ! ---------------------------------------------------------------------------

    nElems = cubes%mesh_list(currentLevel)%descriptor%elem%nElems(eT_fluid)


    !============================ PREDICTOR ===================================

    ! Init the arrays with the intermediate results
    do iStep = 1, 4
      timestep_list(currentLevel)%timestepData(iStep)%state(:,:,:) = 0.0_rk
    end do
    allocate(statedata_list_temp(currentLevel)%state(       &
      & size(cubes%statedata_list(currentLevel)%state,1), &
      & size(cubes%statedata_list(currentLevel)%state,2), &
      & size(cubes%statedata_list(currentLevel)%state,3) ))

    statedata_list_temp(currentLevel)%state = &
      & cubes%statedata_list(currentLevel)%state


    ! 1st runge kutta substep
    timestep_list(currentLevel)%timestepInfoInteger(1) = 1
    call local_predictor_substep(                               &
      & mesh             = cubes%mesh_list(currentLevel),       &
      & stateData        = statedata_list_temp(currentLevel),   &
      & kernelData       = cubes%kerneldata_list(currentLevel), &
      & scheme           = cubes%scheme_list(currentLevel),     &
      & general          = general,                             &
      & poly_proj_list   = poly_proj_list,                      &
      & timestep         = timestep_list(currentLevel),         &
      & equation         = equation,                            &
      & boundary         = cubes%boundary_list(currentLevel),   &
      & material         = cubes%material_list(currentLevel)    )

    call compute_intermediate1(                                             &
      & nTotal    = size(cubes%mesh_list(currentLevel)%descriptor%total,1), &
      & nDofs     = cubes%scheme_list(currentLevel)%nDofs,                  &
      & nScalars  = equation%varSys%nScalars,                               &
      & state_tmp = statedata_list_temp(currentLevel)%state,                &
      & state     = cubes%statedata_list(currentLevel)%state,               &
      & state1    = timestep_list(currentLevel)%timestepData(1)%state,      &
      & dt        = cubes%scheme_list(currentLevel)%time%dt                 )

!! Stabilize the intermediate result of this stage of the predcor_cerk4 scheme
!call atl_stabilize( minlevel = minlevel, maxlevel = maxlevel, &
!                  & statedata_list = statedata_list_temp, &
!                  & mesh_list = mesh_list, &
!                  & scheme_list = scheme_list )




    ! 2nd runge kutta substep
    cubes%statedata_list(currentLevel)%local_time%sim &
      & = cubes%statedata_list(currentLevel)%local_time%sim &
      &   + (12.0_rk/23.0_rk) * cubes%scheme_list(currentLevel)%time%dt
    ! update the intermediate time for next predcor_cerk4 substep, this necessary to calculate
    ! source terms and time-dependent boundary values correctly.
    statedata_list_temp(currentLevel)%local_time%sim &
      & = cubes%statedata_list(currentLevel)%local_time%sim
    ! Store, that we are in the second step of the RK method.
    timestep_list(currentLevel)%timestepInfoInteger(1) = 2

    call local_predictor_substep(                               &
      & mesh             = cubes%mesh_list(currentLevel),       &
      & stateData        = statedata_list_temp(currentLevel),   &
      & kernelData       = cubes%kerneldata_list(currentLevel), &
      & scheme           = cubes%scheme_list(currentLevel),     &
      & general          = general,                             &
      & poly_proj_list   = poly_proj_list,                      &
      & timestep         = timestep_list(currentLevel),         &
      & equation         = equation,                            &
      & boundary         = cubes%boundary_list(currentLevel),   &
      & material         = cubes%material_list(currentLevel)    )

    call compute_intermediate2(                                             &
      & nTotal    = size(cubes%mesh_list(currentLevel)%descriptor%total,1), &
      & nDofs     = cubes%scheme_list(currentLevel)%nDofs,                  &
      & nScalars  = equation%varSys%nScalars,                               &
      & state_tmp = statedata_list_temp(currentLevel)%state,                &
      & state     = cubes%statedata_list(currentLevel)%state,               &
      & state1    = timestep_list(currentLevel)%timestepData(1)%state,      &
      & state2    = timestep_list(currentLevel)%timestepData(2)%state,      &
      & dt        = cubes%scheme_list(currentLevel)%time%dt                 )

!! Stabilize the intermediate result of this stage of the predcor_cerk4 scheme
!call atl_stabilize( minlevel = minlevel, maxlevel = maxlevel, &
!                  & statedata_list = statedata_list_temp, &
!                  & mesh_list = mesh_list, &
!                  & scheme_list = scheme_list )




    ! 3rd runge kutta substep
    cubes%statedata_list(currentLevel)%local_time%sim &
      & = cubes%statedata_list(currentLevel)%local_time%sim &
      &   + (4.0_rk/5.0_rk) * cubes%scheme_list(currentLevel)%time%dt
    ! update the intermediate time for next predcor_cerk4 substep, this necessary to calculate
    ! source terms and time-dependent boundary values correctly.
    statedata_list_temp(currentLevel)%local_time%sim &
      & = cubes%statedata_list(currentLevel)%local_time%sim
    ! Store, that we are in the third step of the RK method.
    timestep_list(currentLevel)%timestepInfoInteger(1) = 3

    call local_predictor_substep(                               &
      & mesh             = cubes%mesh_list(currentLevel),       &
      & stateData        = statedata_list_temp(currentLevel),   &
      & kernelData       = cubes%kerneldata_list(currentLevel), &
      & scheme           = cubes%scheme_list(currentLevel),     &
      & general          = general,                             &
      & poly_proj_list   = poly_proj_list,                      &
      & timestep         = timestep_list(currentLevel),         &
      & equation         = equation,                            &
      & boundary         = cubes%boundary_list(currentLevel),   &
      & material         = cubes%material_list(currentLevel)    )

    call compute_intermediate3(                                             &
      & nTotal    = size(cubes%mesh_list(currentLevel)%descriptor%total,1), &
      & nDofs     = cubes%scheme_list(currentLevel)%nDofs,                  &
      & nScalars  = equation%varSys%nScalars,                               &
      & state_tmp = statedata_list_temp(currentLevel)%state,                &
      & state     = cubes%statedata_list(currentLevel)%state,               &
      & state1    = timestep_list(currentLevel)%timestepData(1)%state,      &
      & state2    = timestep_list(currentLevel)%timestepData(2)%state,      &
      & state3    = timestep_list(currentLevel)%timestepData(3)%state,      &
      & dt        = cubes%scheme_list(currentLevel)%time%dt                 )

!! Stabilize the intermediate result of this stage of the predcor_cerk4 scheme
!call atl_stabilize( minlevel = minlevel, maxlevel = maxlevel, &
!                  & statedata_list = statedata_list_temp, &
!                  & mesh_list = mesh_list, &
!                  & scheme_list = scheme_list )




    ! 4th runge kutta substep
    cubes%statedata_list(currentLevel)%local_time%sim &
      & = cubes%statedata_list(currentLevel)%local_time%sim &
      &   + cubes%scheme_list(currentLevel)%time%dt
    ! update the intermediate time for next predcor_cerk4 substep, this necessary to calculate
    ! source terms and time-dependent boundary values correctly.
    statedata_list_temp(currentLevel)%local_time%sim &
      & = cubes%statedata_list(currentLevel)%local_time%sim
    ! indicates that we are in the fourth step of RK method
    timestep_list(currentLevel)%timestepInfoInteger(1) = 4

    ! Just, the final substep, no need to store an intermediate
    ! state here.
    call local_predictor_substep(                               &
      & mesh             = cubes%mesh_list(currentLevel),       &
      & stateData        = statedata_list_temp(currentLevel),   &
      & kernelData       = cubes%kerneldata_list(currentLevel), &
      & scheme           = cubes%scheme_list(currentLevel),     &
      & general          = general,                             &
      & poly_proj_list   = poly_proj_list,                      &
      & timestep         = timestep_list(currentLevel),         &
      & equation         = equation,                            &
      & boundary         = cubes%boundary_list(currentLevel),   &
      & material         = cubes%material_list(currentLevel)    )


    !============================ CORRECTOR ===================================


    ! first gauss-point in time
    call compute_prediction(                                                &
      & nTotal    = size(cubes%mesh_list(currentLevel)%descriptor%total,1), &
      & nDofs     = cubes%scheme_list(currentLevel)%nDofs,                  &
      & nScalars  = equation%varSys%nScalars,                               &
      & state_tmp = statedata_list_temp(currentLevel)%state,                &
      & state     = cubes%statedata_list(currentLevel)%state,               &
      & state1    = timestep_list(currentLevel)%timestepData(1)%state,      &
      & state2    = timestep_list(currentLevel)%timestepData(2)%state,      &
      & state3    = timestep_list(currentLevel)%timestepData(3)%state,      &
      & state4    = timestep_list(currentLevel)%timestepData(4)%state,      &
      & dt        = cubes%scheme_list(currentLevel)%time%dt,                &
      & theta     = 0.5_rk*(1-1/sqrt(3.0_rk))                               )


    ! apply corrector substep
    call global_corrector_substep(                           &
      & minLevel              = minLevel,                    &
      & maxLevel              = maxLevel,                    &
      & currentLevel          = currentLevel,                &
      & mesh_list             = cubes%mesh_list,             &
      & tree                  = tree,                        &
      & kerneldata_list       = cubes%kerneldata_list,       &
      & statedata_list        = statedata_list_temp,         &
      & facedata_list         = cubes%facedata_list,         &
      & source                = cubes%source,                &
      & penalizationData_list = cubes%penalizationdata_list, &
      & boundary_list         = cubes%boundary_list,         &
      & bc                    = cubes%bc,                    &
      & scheme_list           = cubes%scheme_list,           &
      & poly_proj_pos         = cubes%poly_proj_pos,         &
      & poly_proj_list        = poly_proj_list,              &
      & timestep_list         = timestep_list,               &
      & equation              = equation,                    &
      & material_list         = cubes%material_list,         &
      & general               = general                      )

    ! save old state
    statedata_list_temp(currentLevel)%state(:nElems,:,:) &
      & = cubes%statedata_list(currentLevel)%state(:nElems,:,:)
    ! add result to state
    call update_state(                                                      &
      & nTotal    = size(cubes%mesh_list(currentLevel)%descriptor%total,1), &
      & nDofs     = cubes%scheme_list(currentLevel)%nDofs,                  &
      & nScalars  = equation%varSys%nScalars,                               &
      & nElems    = nElems,                                                 &
      & state     = cubes%statedata_list(currentLevel)%state,               &
      & state_der = cubes%kerneldata_list(currentLevel)%state_der,          &
      & factor    = 0.5_rk*cubes%scheme_list(currentLevel)%time%dt          )

    ! second gauss-point in time
    call compute_prediction(                                                &
      & nTotal    = size(cubes%mesh_list(currentLevel)%descriptor%total,1), &
      & nDofs     = cubes%scheme_list(currentLevel)%nDofs,                  &
      & nScalars  = equation%varSys%nScalars,                               &
      & state_tmp = statedata_list_temp(currentLevel)%state,                &
      & state     = statedata_list_temp(currentLevel)%state,                &
      & state1    = timestep_list(currentLevel)%timestepData(1)%state,      &
      & state2    = timestep_list(currentLevel)%timestepData(2)%state,      &
      & state3    = timestep_list(currentLevel)%timestepData(3)%state,      &
      & state4    = timestep_list(currentLevel)%timestepData(4)%state,      &
      & dt        = cubes%scheme_list(currentLevel)%time%dt,                &
      & theta     = 0.5_rk*(1+1/sqrt(3.0_rk))                               )

    ! apply corrector substep
    call global_corrector_substep(                           &
      & minLevel              = minLevel,                    &
      & maxLevel              = maxLevel,                    &
      & currentLevel          = currentLevel,                &
      & mesh_list             = cubes%mesh_list,             &
      & tree                  = tree,                        &
      & kerneldata_list       = cubes%kerneldata_list,       &
      & statedata_list        = statedata_list_temp,         &
      & facedata_list         = cubes%facedata_list,         &
      & source                = cubes%source,                &
      & penalizationdata_list = cubes%penalizationdata_list, &
      & boundary_list         = cubes%boundary_list,         &
      & bc                    = cubes%bc,                    &
      & scheme_list           = cubes%scheme_list,           &
      & poly_proj_pos         = cubes%poly_proj_pos,         &
      & poly_proj_list        = poly_proj_list,              &
      & timestep_list         = timestep_list,               &
      & equation              = equation,                    &
      & material_list         = cubes%material_list,         &
      & general               = general                      )

    ! add result to state
    call update_state(                                                      &
      & nTotal    = size(cubes%mesh_list(currentLevel)%descriptor%total,1), &
      & nDofs     = cubes%scheme_list(currentLevel)%nDofs,                  &
      & nScalars  = equation%varSys%nScalars,                               &
      & nElems    = nElems,                                                 &
      & state     = cubes%statedata_list(currentLevel)%state,               &
      & state_der = cubes%kerneldata_list(currentLevel)%state_der,          &
      & factor    = 0.5_rk*cubes%scheme_list(currentLevel)%time%dt          )


!    ! Stabilize the final result
!    call atl_stabilize( minlevel = minlevel, maxlevel = maxlevel, &
!                      & statedata_list = statedata_list, &
!                      & mesh_list = mesh_list, &
!                      & scheme_list = scheme_list )


  contains


    ! own subroutine to allow vectorization with gcc
    subroutine compute_intermediate1( nTotal, nDofs, nScalars, state_tmp, &
      &                               state, state1, dt                   )
      ! ---------------------------------------------------------------------------
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(out) :: state_tmp(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state1(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: dt
      ! ---------------------------------------------------------------------------

      state_tmp = state + (12.0_rk/23.0_rk*dt)*state1

    end subroutine compute_intermediate1


    ! own subroutine to allow vectorization with gcc
    subroutine compute_intermediate2( nTotal, nDofs, nScalars, state_tmp, &
      &                               state, state1, state2, dt           )
      ! ---------------------------------------------------------------------------
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(out) :: state_tmp(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state1(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state2(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: dt
      ! ---------------------------------------------------------------------------

      state_tmp = state - (68.0_rk/375.0_rk*dt)*state1 + (368.0_rk/375.0_rk*dt)*state2

    end subroutine compute_intermediate2


    ! own subroutine to allow vectorization with gcc
    subroutine compute_intermediate3( nTotal, nDofs, nScalars, state_tmp, &
      &                               state, state1, state2, state3, dt   )
      ! ---------------------------------------------------------------------------
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(out) :: state_tmp(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state1(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state2(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state3(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: dt
      ! ---------------------------------------------------------------------------

      state_tmp = state                        &
        & + (31.0_rk/144.0_rk * dt) * state1   &
        & + (529.0_rk/1152.0_rk * dt) * state2 &
        & + (125.0_rk/384.0_rk * dt) * state3

    end subroutine compute_intermediate3


    ! own subroutine to allow vectorization with gcc
    subroutine compute_prediction( nTotal, nDofs, nScalars, state_tmp, state, &
      &                            state1, state2, state3, state4, dt, theta  )
      ! ---------------------------------------------------------------------------
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(out) :: state_tmp(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state1(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state2(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state3(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state4(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: dt, theta
      ! ---------------------------------------------------------------------------
      real(kind=rk) :: b(4)
      ! ---------------------------------------------------------------------------

      ! compute coefficients
      b(1) = dt*((41.0_rk/72.0_rk*theta - 65.0_rk/48.0_rk)*theta + 1)*theta
      b(2) = dt*(-529.0_rk/576.0_rk*theta + 529.0_rk/384.0_rk)*theta**2
      b(3) = dt*(-125.0_rk/192.0_rk*theta + 125.0_rk/128.0_rk)*theta**2
      b(4) = dt*(theta-1)*theta**2


      state_tmp = state   &
        & + b(1) * state1 &
        & + b(2) * state2 &
        & + b(3) * state3 &
        & + b(4) * state4

    end subroutine compute_prediction



    ! own subroutine to allow vectorization with gcc
    subroutine update_state( nTotal, nDofs, nScalars, nElems, state, &
      &                      state_der, factor                       )
      ! ---------------------------------------------------------------------------
      integer, intent(in) :: nTotal, nDofs, nScalars, nElems
      real(kind=rk), intent(inout) :: state(nTotal, nDofs, nScalars)
      !> The state derivatives, passed with assumed shape because of a potential
      !! padding applied for odd polynomial degrees.
      real(kind=rk), intent(in) :: state_der(:,:,:)
      real(kind=rk), intent(in) :: factor
      ! ---------------------------------------------------------------------------
      integer :: iVar
      ! ---------------------------------------------------------------------------


      do iVar = 1, nScalars
        !! Limit the second index by an upper bound of nDofs, as there
        !! potentially is a padding on this index.
        state(:nElems,:,iVar) = state(:nElems,:,iVar) &
          & + factor * state_der(:nElems,:nDofs,iVar)
      end do


    end subroutine update_state

  end subroutine mesh_timestep_predcor_cerk4
  ! ****************************************************************************



  ! ****************************************************************************
  !> Subroutine calculates a substep of the local predictor
  subroutine local_predictor_substep( mesh, kerneldata, statedata, scheme, &
    &                                 general, poly_proj_list, timestep,   &
    &                                 equation, material, boundary         )
    ! ---------------------------------------------------------------------------
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type), intent(inout) :: mesh
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme
    !> General treelm settings
    type(tem_general_type), intent(inout) :: general
    !> List of project, for each level.
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep
    !> The equation you are operating with.
    type(atl_equations_type),intent(inout) :: equation
    !> Material parameter description.
    type(atl_material_type), intent(inout) :: material
    type(atl_level_boundary_type), intent(in) :: boundary
    ! ---------------------------------------------------------------------------
    integer :: nElems
    ! ---------------------------------------------------------------------------

    !> TODO here need to be the source position pointer????
    call atl_preprocess_local_rhs(       &
      & mesh           = mesh,           &
      & statedata      = statedata,      &
      & boundary       = boundary,       &
      & scheme         = scheme,         &
      & general        = general,        &
      & poly_proj_list = poly_proj_list, &
      & equation       = equation,       &
      & material       = material        )

    nElems = mesh%descriptor%elem%nElems(eT_fluid)


    ! copy output to timestepData
    ! Limit the second index with an upper bound because of a potential padding
    timestep%timestepData(timestep%timestepInfoInteger(1))%state(:nElems,:,:) = &
     & kerneldata%state_der(:nElems,:kerneldata%nDofs,:)


  end subroutine local_predictor_substep
  ! ****************************************************************************



  ! ****************************************************************************
  !> Subroutine calculates a substep of corrector,
  !! this is the same as a usual substep of the RKDG
  recursive subroutine global_corrector_substep(minLevel, maxLevel,        &
    & currentLevel, mesh_list, tree, kerneldata_list, statedata_list,      &
    & facedata_list, source, penalizationdata_list, boundary_list, bc,     &
    & scheme_list, poly_proj_pos, poly_proj_list, timestep_list, equation, &
    & material_list, general                                               )
    ! ---------------------------------------------------------------------------
    ! Workaround ICE in Intel 16: If this use appears in the module-wide use
    !                             list above, we get an ICE in
    !                             atl_global_time_integration, by putting the use
    !                             here, we can avoid that error.
    use atl_compute_module, only: atl_compute_rhs,    &
      &                           atl_preprocess_rhs, &
      &                           atl_postprocess_rhs
    !> The minimum refinement level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum refinement level of the mesh.
    integer, intent(in) :: maxLevel
    !> The level the timestep has to be performed for.
    integer, intent(in) :: currentLevel
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type), intent(inout) :: mesh_list( minLevel:maxLevel )
    !> treelm mesh
    type(treelmesh_type) :: tree
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: &
      & kerneldata_list(minLevel:maxLevel)
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: &
      & statedata_list( minLevel:maxLevel )
    !> List of faces states you want to calc the rhs for. For each level we have
    !! one.
    type(atl_facedata_type), intent(inout) :: facedata_list( minLevel:maxLevel )
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
    type(atl_scheme_type), intent(inout) :: scheme_list( minLevel:maxLevel)
    !> List of projection pointer, for each level.
    integer, intent(inout) :: poly_proj_pos(minLevel:maxLevel)
    !> List of projection, for each level.
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep_list( minLevel:maxLevel )
    !> The equation you are operating with.
    type(atl_equations_type),intent(inout) :: equation
    !> Material parameter description.
    type(atl_material_type), intent(inout) :: material_list( minlevel:maxlevel )
    !> General treelm settings
    type(tem_general_type), intent(inout) :: general
    ! ---------------------------------------------------------------------------


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

    call atl_compute_rhs( minLevel, maxLevel, currentLevel, mesh_list, tree, &
      & kerneldata_list, statedata_list, facedata_list, source,              &
      & penalizationdata_list, scheme_list,poly_proj_pos, poly_proj_list,    &
      & equation, material_list, general                                     )

    call atl_postprocess_rhs( mesh       = mesh_list(currentLevel),       &
      &                       kerneldata = kerneldata_list(currentLevel), &
      &                       statedata  = statedata_list(currentLevel),  &
      &                       scheme     = scheme_list(currentLevel),     &
      &                       timestep   = timestep_list(currentLevel),   &
      &                       equation   = equation                       )


  end subroutine global_corrector_substep
  ! ****************************************************************************

end module atl_predcor_cerk4_module
