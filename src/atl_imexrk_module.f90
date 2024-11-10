! Copyright (c) 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Routines,functions, datatypes related to IMEX Runge Kutta timestepping methods.
!! \author{Jens Zudrop}
module atl_imexrk_module
  use env_module,               only: rk
  use tem_float_module,         only: operator(.fne.)
  use tem_aux_module,           only: tem_abort
  use tem_general_module,       only: tem_general_type
  use tem_element_module,       only: eT_fluid
  use tem_logging_module,       only: logUnit
  use treelmesh_module,         only: treelmesh_type
  use tem_precice_module,       only: precice_available

  use atl_elemental_time_integration_module, only: atl_timestep_type
  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_kerneldata_module,        only: atl_kerneldata_type, &
    &                                     atl_statedata_type
  use atl_scheme_module,            only: atl_scheme_type, &
    &                                     atl_local_timestep_type
  use atl_boundary_module,          only: atl_level_boundary_type
  use atl_source_types_module,      only: atl_source_type
  use atl_facedata_module,          only: atl_facedata_type
  use atl_time_integration_module,  only: atl_global_timestep_type
  use atl_equation_module,      only: atl_equations_type
  use atl_bc_header_module,     only: atl_boundary_type
  use atl_stabilize_module,     only: atl_stabilize
  use ply_poly_project_module,  only: ply_poly_project_type, &
    &                                 assignment(=),         &
    &                                 ply_poly_project_n2m,  &
    &                                 ply_poly_project_m2n
  use atl_materialPrp_module,   only: atl_material_type
  use atl_writePrecice_module,  only: atl_write_precice
  use atl_cube_container_module, only: atl_cube_container_type
  use atl_penalization_module,  only: atl_penalizationData_type
  use atl_reference_element_module, only: atl_refToPhysCoord

  implicit none
  private

  public :: atl_init_imexRungeKutta


contains


  ! ------------------------------------------------------------------------ !
  !> Routine to initialize IMEX Runge-Kutta scheme for timestepping.
  subroutine atl_init_imexRungeKutta( me, minLevel, maxLevel, &
    &                                 steps, statedata_list   )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    integer :: iLevel, iStep
    ! -------------------------------------------------------------------- !
    ! point to euler step function, we dont have to allocate any of the
    ! additional arrays for coefficients or buffering.
    allocate( me%elementSteps(minLevel:maxLevel) )
    select case(steps)
    case(4) ! we have runge kutta with 4 steps, and third order
      do iLevel = minLevel, maxLevel
        if(allocated(statedata_list(iLevel)%state)) then
          ! We store one integer, indicating the elemental operation in which
          ! step we are.
          allocate(me%elementSteps(iLevel)%timestepInfoInteger(1))
          ! In a 4 step runge kutta method we have to store 4 intermediate
          ! results.
          allocate(me%elementSteps(iLevel)%timestepData(7))
          do iStep = 1, 7
            allocate(me%elementSteps(iLevel)%timestepData(iStep) &
              &         %state(size(statedata_list(iLevel)%state, 1), &
              &                size(statedata_list(iLevel)%state, 2), &
              &                size(statedata_list(iLevel)%state, 3)) &
              &     )
          end do
          me%elementSteps(iLevel)%elemStep => elemental_timestep_imexrk
          me%elementSteps(iLevel)%elemStep_vec => elemental_timestep_vec_imexrk
          me%elementSteps(iLevel)%updateStep => update_timestep_imexrk
        end if
      end do
      me%meshStep => mesh_timestep_imexrk

    case default
      write(logUnit(1),*) 'IMEX Runge-Kutta method for ',  steps, &
        &                 ' steps not implemented!'
      write(logUnit(1),*) 'stopping...'
      call tem_abort()
    end select
  end subroutine atl_init_imexRungeKutta
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !



  ! ------------------------------------------------------------------------ !
  !> Levelwise update of IMEX Runge-Kutta
  subroutine update_timestep_imexrk( me, timestepInfo )
    !> The type of your timestepping.
    class(atl_timestep_type), intent(inout) :: me

    !> Local timestepping information for that part of the mesh
    type(atl_local_timestep_type), intent(in) :: timestepInfo
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    ! nothing to be done here.
  end subroutine update_timestep_imexrk
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !



  ! ------------------------------------------------------------------------ !
  !> Elemental operation for timestepping IMEX-RK.
  subroutine elemental_timestep_imexrk( me, state, cell, dof, sideFlux)
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    me%timestepData(me%timestepInfoInteger(1))%state(cell,dof,:)        &
      &  = me%timestepData(me%timestepInfoInteger(1))%state(cell,dof,:) &
      &  + sideFlux
  end subroutine elemental_timestep_imexrk
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Elemental operation for timestepping of order 4.
  subroutine elemental_timestep_vec_imexrk( me, state, kerneldata)
    ! -------------------------------------------------------------------- !
    !> Description of the timestep integration method.
    class(atl_timestep_type), intent(inout) :: me

    !> The state of all cells on this level. This field will be updated
    !! ad the cell position. See kerneldata type for more explanations.
    real(kind=rk), intent(inout) :: state(:,:,:)

    !> Complete kerneldata to get the flux from with additional information.
    type(atl_kerneldata_type), intent(in) :: kerneldata
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    call compute(kerneldata%nTotal,                                &
      &          kerneldata%nDofs,                                 &
      &          kerneldata%nVars,                                 &
      &          me%timestepData(me%timestepInfoInteger(1))%state, &
      &          kerneldata%state_der                              )

  contains

    ! put it in an own subroutine to allow vectorization with gcc
    subroutine compute(nTotal, nDofs, nScalars,  &
        &              state, state_der)
      ! ---------------------------------------------------------------- !
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(inout) :: state(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state_der(:,:,:)
      ! ---------------------------------------------------------------- !

      state(:nTotal,:,:) = state(:nTotal,:,:) + state_der(:nTotal,:nDofs,:)

    end subroutine compute

  end subroutine elemental_timestep_vec_imexrk
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine for timestepping with IMEX Runge-Kutta.
  subroutine mesh_timestep_imexrk( minLevel, maxLevel, currentLevel, cubes, &
    &                              tree, timestep_list, nSteps, equation,   &
    &                              general, commStateTimer, poly_proj_list  )
    ! -------------------------------------------------------------------- !
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
    type(atl_timestep_type), intent(inout) :: timestep_list(minLevel:)

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
    ! -------------------------------------------------------------------- !
    integer :: iLevel, iStep
    type(atl_statedata_type) :: statedata_list_temp(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !

    ! Init the arrays with the intermediate results
    do iLevel = minLevel, maxLevel
      do iStep = 1, 7
        timestep_list(iLevel)%timestepData(iStep)%state(:,:,:) = 0.0_rk
      end do
      allocate( statedata_list_temp(iLevel)         &
        &         %state( cubes%mesh_list(iLevel)   &
        &                      %descriptor          &
        &                      %elem                &
        &                      %nElems(eT_fluid),   &
        &                 cubes%scheme_list(iLevel) &
        &                      %nDofs,              &
        &                 equation%varSys           &
        &                         %nScalars )       )
      statedata_list_temp(iLevel)%state = 0.0_rk
    end do


    ! 1st runge kutta substep. Will call itself on the finer
    ! levels recursively.
    do iLevel = minLevel, maxLevel
      ! update the intermediate time for next imexrk substep, this necessary to
      ! calculate source terms and time-dependent boundary values correctly.
      statedata_list_temp(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim
      ! Store, that we are in the first step of the RK method.
      timestep_list(iLevel)%timestepInfoInteger(1) = 1
    end do
    call imexrk_substep(                                           &
      &    minLevel, maxLevel, currentLevel,                       &
      &    cubes%mesh_list, tree, cubes%levelpointer,              &
      &    cubes%kerneldata_list,                                  &
      &    cubes%statedata_list, cubes%facedata_list,              &
      &    cubes%source, cubes%penalizationdata_list,              &
      &    cubes%boundary_list, cubes%bc,                          &
      &    cubes%scheme_list, cubes%poly_proj_pos, poly_proj_list, &
      &    timestep_list, equation, cubes%material_list, general,  &
      &    commStateTimer                                          )

    ! Now, build the intermediate values for the RK4 method. Since the substep
    ! was carried out for all the levels in the previous step, we can continue
    ! in a simple loop over the levels.
    do iLevel = minLevel, maxLevel
      statedata_list_temp(iLevel)%state                      &
        & = cubes%statedata_list(iLevel)%state               &
        &   + cubes%scheme_list(iLevel)%time%dt              &
        &     * ( 0.4358665215_rk * timestep_list(iLevel)    &
        &                           %timestepData(1)%state   )
    end do

    ! make the first implicit solve
    do iLevel = minLevel, maxLevel
      timestep_list(iLevel)%timestepInfoInteger(1) = 5
    end do
    do iLevel = minLevel, maxLevel
      call implicit_update(                                                 &
        &    equation      = equation,                                      &
        &    statedata     = statedata_list_temp(iLevel),                   &
        &    timestep_list = timestep_list(iLevel),                         &
        &    poly_proj     = poly_proj_list(cubes%material_list(iLevel)     &
        &                                        %poly_proj_pos_state2Mat), &
        &    material      = cubes%material_list(iLevel),                   &
        &    dt            = cubes%scheme_list(iLevel)%time%dt,             &
        &    mesh          = cubes%mesh_list(iLevel),                       &
        &    weight        = 0.4358665215_rk                                )
    end do

    ! stabilize the result of the first implicit/explicit stage
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


    ! 2nd runge kutta substep. Will call itself on the finer
    ! levels recursively. Before, we enter the recursive call,
    ! update to the new intermediate timepoint of the RK4 method.
    do iLevel = minLevel, maxLevel
      cubes%statedata_list(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim &
        &  + 0.4358665215_rk * cubes%scheme_list(iLevel)%time%dt
      ! update the intermediate time for next imexrk substep, this necessary to
      ! calculate source terms and time-dependent boundary values correctly.
      statedata_list_temp(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim
      ! Store, that we are in the second step of the RK method.
      timestep_list(iLevel)%timestepInfoInteger(1) = 2
    end do

    call imexrk_substep(                                           &
      &    minLevel, maxLevel, currentLevel,                       &
      &    cubes%mesh_list, tree, cubes%levelpointer,              &
      &    cubes%kerneldata_list,                                  &
      &    statedata_list_temp, cubes%facedata_list,               &
      &    cubes%source, cubes%penalizationdata_list,              &
      &    cubes%boundary_list, cubes%bc,                          &
      &    cubes%scheme_list, cubes%poly_proj_pos, poly_proj_list, &
      &    timestep_list, equation, cubes%material_list, general,  &
      &    commStateTimer                                          )

    ! Now, build the intermediate state of the RK4 method.
    do iLevel = minLevel, maxLevel
      statedata_list_temp(iLevel)%state                      &
        & = cubes%statedata_list(iLevel)%state               &
        &   + cubes%scheme_list(iLevel)%time%dt              &
        &     * ( 0.3212788860_rk * timestep_list(iLevel)    &
        &                           %timestepData(1)%state   &
        &         + 0.3966543747_rk * timestep_list(iLevel)  &
        &                             %timestepData(2)%state &
        &         + 0.2820667392_rk * timestep_list(iLevel)  &
        &                             %timestepData(5)%state )
    end do

    ! make the implicit solve
    do iLevel = minLevel, maxLevel
      timestep_list(iLevel)%timestepInfoInteger(1) = 6
    end do
    do iLevel = minLevel, maxLevel
      call implicit_update(                                                 &
        &    equation      = equation,                                      &
        &    statedata     = statedata_list_temp(iLevel),                   &
        &    timestep_list = timestep_list(iLevel),                         &
        &    poly_proj     = poly_proj_list(cubes%material_list(iLevel)     &
        &                                        %poly_proj_pos_state2Mat), &
        &    material      = cubes%material_list(iLevel),                   &
        &    dt            = cubes%scheme_list(iLevel)%time%dt,             &
        &    mesh          = cubes%mesh_list(iLevel),                       &
        &    weight        = 0.4358665215_rk                                )
    end do

    ! Stabilize the intermediate result of this stage of the RK4 scheme
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


    ! 3rd runge kutta substep. Will call itself on the finer
    ! levels recursively. No update of the timepoint in the 3rd
    ! step.
    do iLevel = minLevel, maxLevel
      cubes%statedata_list(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim &
        &  + (0.7179332608_rk-0.4358665215_rk) * cubes%scheme_list(iLevel) &
        &                                             %time%dt
      ! update the intermediate time for next imexrk substep, this necessary to
      ! calculate source terms and time-dependent boundary values correctly.
      statedata_list_temp(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim
      ! Store, that we are in the third step of the RK method.
      timestep_list(iLevel)%timestepInfoInteger(1) = 3
    end do

    call imexrk_substep(                                           &
      &    minLevel, maxLevel, currentLevel,                       &
      &    cubes%mesh_list, tree, cubes%levelpointer,              &
      &    cubes%kerneldata_list,                                  &
      &    statedata_list_temp, cubes%facedata_list,               &
      &    cubes%source, cubes%penalizationdata_list,              &
      &    cubes%boundary_list, cubes%bc,                          &
      &    cubes%scheme_list, cubes%poly_proj_pos, poly_proj_list, &
      &    timestep_list, equation, cubes%material_list, general,  &
      &    commStateTimer                                          )

    ! Now, build the intermediate state of the RK4 method.
    do iLevel = minLevel, maxLevel
      statedata_list_temp(iLevel)%state                      &
        & = cubes%statedata_list(iLevel)%state               &
        &   + cubes%scheme_list(iLevel)%time%dt              &
        &     * ( - 0.1058582960_rk * timestep_list(iLevel)  &
        &                           %timestepData(1)%state   &
        &         + 0.5529291479_rk * timestep_list(iLevel)  &
        &                             %timestepData(2)%state &
        &         + 0.5529291479_rk * timestep_list(iLevel)  &
        &                             %timestepData(3)%state &
        &         + 1.2084966490_rk * timestep_list(iLevel)  &
        &                             %timestepData(5)%state &
        &         - 0.6443631710_rk * timestep_list(iLevel)  &
        &                             %timestepData(6)%state )
    end do

    ! make the implicit solve
    do iLevel = minLevel, maxLevel
      timestep_list(iLevel)%timestepInfoInteger(1) = 7
    end do
    do iLevel = minLevel, maxLevel
      call implicit_update(                                                 &
        &    equation      = equation,                                      &
        &    statedata     = statedata_list_temp(iLevel),                   &
        &    timestep_list = timestep_list(iLevel),                         &
        &    poly_proj     = poly_proj_list(cubes%material_list(iLevel)     &
        &                                        %poly_proj_pos_state2Mat), &
        &    material      = cubes%material_list(iLevel),                   &
        &    dt            = cubes%scheme_list(iLevel)%time%dt,             &
        &    mesh          = cubes%mesh_list(iLevel),                       &
        &    weight        = 0.4358665215_rk                                )
    end do

    ! Stabilize the intermediate result of this stage of the RK4 scheme
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


    ! 4th runge kutta substep. Will call itself on the finer
    ! levels recursively. Update again to the new timepoint of
    ! the next RK substep.
    do iLevel = minLevel, maxLevel
      cubes%statedata_list(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim &
        &  + (1.0_rk-0.7179332608_rk)*cubes%scheme_list(iLevel)%time%dt
      ! update the intermediate time for next imexrk substep, this necessary to
      ! calculate source terms and time-dependent boundary values correctly.
      statedata_list_temp(iLevel)%local_time%sim &
        &  = cubes%statedata_list(iLevel)%local_time%sim
      ! indicates that we are in the fourth step of RK method
      timestep_list(iLevel)%timestepInfoInteger(1) = 4
    end do
    ! Just, the final substep, no need to store an intermediate
    ! state here.

    call imexrk_substep(                                           &
      &    minLevel, maxLevel, currentLevel,                       &
      &    cubes%mesh_list, tree, cubes%levelpointer,              &
      &    cubes%kerneldata_list,                                  &
      &    statedata_list_temp, cubes%facedata_list,               &
      &    cubes%source, cubes%penalizationdata_list,              &
      &    cubes%boundary_list, cubes%bc,                          &
      &    cubes%scheme_list, cubes%poly_proj_pos, poly_proj_list, &
      &    timestep_list, equation, cubes%material_list, general,  &
      &    commStateTimer                                          )
    ! No stabilization here, as it is done after the imexrk update step.

    ! The final update step. No need to have a recusrive structure here,
    ! so, we can do it in parallel for the different levels.
    do iLevel = minLevel, maxLevel
      call imexrk_update( cubes%statedata_list(iLevel), &
        &                 cubes%scheme_list(iLevel),    &
        &                 timestep_list(iLevel)         )
    end do


    ! Stabilize the final result
    call atl_stabilize(                                     &
      &    minlevel            = minlevel,                  &
      &    maxlevel            = maxlevel,                  &
      &    statedata_list      = cubes%statedata_list,      &
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


    ! After the runge kutta FULL step, write the results to precice
    if (precice_available) then
      call atl_write_precice( stFunList = equation%stFunList,     &
        &                     varSys    = equation%varSys,        &
        &                     time      = general%simControl%now, &
        &                     tree      = tree                    )
    end if

  end subroutine mesh_timestep_imexrk
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine calculates the final update step of the Runge-Kutta method.
  !! It is performing levelwise.
  subroutine imexrk_update( statedata_list, scheme_list, timestep_list )
    ! -------------------------------------------------------------------- !
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list
    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep_list
    ! -------------------------------------------------------------------- !
    integer :: nTotal, nDofs, nScalars
    ! -------------------------------------------------------------------- !

    nTotal = size(statedata_list%state, 1)
    nDofs = size(statedata_list%state, 2)
    nScalars = size(statedata_list%state, 3)
    call compute(nTotal, nDofs, nScalars, &
      &          statedata_list%state, &
      &          timestep_list%timestepData(2)%state, &
      &          timestep_list%timestepData(3)%state, &
      &          timestep_list%timestepData(4)%state, &
      &          timestep_list%timestepData(5)%state, &
      &          timestep_list%timestepData(6)%state, &
      &          timestep_list%timestepData(7)%state, &
      &          scheme_list%time%dt )

  contains

    ! own subroutine to allow vectorization with gcc
    subroutine compute(nTotal, nDofs, nScalars, &
        &              state, state2, state3, state4, &
        &              state1_implicit, state2_implicit, state3_implicit, &
        &              dt)
      ! ---------------------------------------------------------------- !
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(inout) :: state(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state2(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state3(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state4(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state1_implicit(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state2_implicit(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state3_implicit(nTotal, nDofs, nScalars)
      real(kind=rk) :: dt
      ! ---------------------------------------------------------------- !

      state = state + (dt) * ( &
         ! & 0.0_rk * state1 &
         & + 1.208496649_rk*(state2+state1_implicit) &
         & + (-0.644363171_rk)*(state3+state2_implicit) &
         & + 0.4358665215_rk*(state4+state3_implicit) &
         & )

    end subroutine compute

  end subroutine imexrk_update
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine calculates a substep of the IMEX Runge-Kutta timestepping
  !!  scheme.
  !!
  !! Calls itself recursively for the finer levels until the finest level is
  !! reached.
  recursive subroutine imexrk_substep(minLevel, maxLevel, currentLevel,        &
    &                          mesh_list, tree, levelpointer, kerneldata_list, &
    &                          statedata_list, facedata_list,                  &
    &                          source, penalizationdata_list,                  &
    &                          boundary_list, bc,                              &
    &                          scheme_list, poly_proj_pos, poly_proj_list,     &
    &                          timestep_list,equation, material_list, general, &
    &                          commStateTimer                                  )
    ! -------------------------------------------------------------------- !
    ! Workaround ICE in Intel 16: If this use appears in the module-wide use
    !                             list above, we get an ICE in
    !                             atl_global_time_integration, by putting the use
    !                             here, we can avoid that error.
    use atl_compute_module,           only: atl_compute_rhs, &
      &                                     atl_preprocess_rhs, &
      &                                     atl_postprocess_rhs
    !> The minimum refinement level of the mesh.
    integer, intent(in)           :: minLevel
    !> The maximum refinement level of the mesh.
    integer, intent(in)           :: maxLevel
    !> The level the timestep has to be performed for.
    integer, intent(in)           :: currentLevel
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type), intent(inout) :: mesh_list( &
                            & minLevel:maxLevel)
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> Pointer for elements from global treeID list index to index in
    !! levelwise fluid lists
    integer, intent(in)  :: levelPointer(:)
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: kerneldata_list( &
                            & minLevel:maxLevel)
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list( &
      &                       minLevel:maxLevel)
    !> List of faces states you want to calc the rhs for.
    !! For each level we have one.
    type(atl_facedata_type), intent(inout) :: facedata_list( &
      &                       minLevel:maxLevel)
    !> List of sources, for each level
    type(atl_source_type), intent(inout) :: source
    !> List of penalization data, for each level.
    type(atl_penalizationData_type), intent(inout) :: penalizationdata_list( &
      &                       minLevel:maxLevel)
    !> List of boundaries, for each level.
    type(atl_level_boundary_type), intent(inout) :: &
      & boundary_list(minLevel:maxLevel)
    !> Global description of the boundaries
    type(atl_boundary_type), intent(in) :: bc(:)
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list( &
      &                       minLevel:maxLevel)
    !> List of position of projection method in unique projection list, for
    !! each level
    integer, intent(in) :: poly_proj_pos(minLevel:maxLevel)
    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep_list( &
      &                       minLevel:maxLevel)
    !> The equation you are operating with.
    type(atl_equations_type),intent(inout)              :: equation
    !> Material parameter description.
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxlevel)
    !> General treelm settings
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout)                        :: commStateTimer
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

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
    if (currentLevel < maxLevel) then
      call imexrk_substep( minLevel              = minLevel,              &
        &                  maxLevel              = maxLevel,              &
        &                  currentLevel          = currentLevel + 1,      &
        &                  mesh_list             = mesh_list,             &
        &                  tree                  = tree,                  &
        &                  levelPointer          = levelPointer,          &
        &                  kerneldata_list       = kerneldata_list,       &
        &                  statedata_list        = statedata_list,        &
        &                  facedata_list         = facedata_list,         &
        &                  source                = source,                &
        &                  penalizationdata_list = penalizationdata_list, &
        &                  boundary_list         = boundary_list,         &
        &                  bc                    = bc,                    &
        &                  scheme_list           = scheme_list,           &
        &                  poly_proj_pos         = poly_proj_pos,         &
        &                  poly_proj_list        = poly_proj_list,        &
        &                  timestep_list         = timestep_list,         &
        &                  equation              = equation,              &
        &                  material_list         = material_list,         &
        &                  general               = general,               &
        &                  commStateTimer        = commStateTimer         )
    end if

    ! After making the timestep on the finer level, we do the following:
    ! 1. Interpolate flux from finer level to the current level.
    ! 2. Calculate the fluxes for the current level.
    ! 3. Communicate the fluxes (only the current level).
    ! 4. Calculate the physical fluxes.
    ! 5. Project physical flux + numerical fluxes to test functions.
    ! 6. Multiply with inverse of (cell local) mass matrix.
    ! ... this call includes steps 1 to 5.
    call atl_compute_rhs( minlevel              = minLevel,              &
      &                   maxlevel              = maxLevel,              &
      &                   currentLevel          = currentLevel,          &
      &                   mesh_list             = mesh_list,             &
      &                   tree                  = tree,                  &
      &                   kerneldata_list       = kerneldata_list,       &
      &                   statedata_list        = statedata_list,        &
      &                   facedata_list         = facedata_list,         &
      &                   source                = source,                &
      &                   penalizationdata_list = penalizationdata_list, &
      &                   scheme_list           = scheme_list,           &
      &                   poly_proj_pos         = poly_proj_pos,         &
      &                   poly_proj_list        = poly_proj_list,        &
      &                   equation              = equation,              &
      &                   material_list         = material_list,         &
      &                   general               = general,               &
      &                   computePenalization   = .false.                )

    ! ... this call applies the inverse of the mass matrix.
    call atl_postprocess_rhs( mesh       = mesh_list(currentLevel),       &
      &                       kerneldata = kerneldata_list(currentLevel), &
      &                       statedata  = statedata_list(currentLevel),  &
      &                       scheme     = scheme_list(currentLevel),     &
      &                       timestep   = timestep_list(currentLevel),   &
      &                       equation   = equation                       )

  end subroutine imexrk_substep
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine implicit_update( equation, statedata,poly_proj, material, dt, &
    &                         timestep_list, weight, mesh                  )
    ! Workaround ICE in Intel 16: If this use appears in the module-wide use
    !                             list above, we get an ICE in
    !                             atl_global_time_integration, by putting the use
    !                             here, we can avoid that error.
    use ply_oversample_module, only: ply_convert2oversample,    &
      &                              ply_convertFromoversample, &
      &                              ply_convert2oversample,    &
      &                              ply_convertFromoversample
    use atl_eqn_maxwell_hlp_module, only: atl_eqn_maxwell_implicit_pen
    use atl_eqn_euler_hlp_module, only: atl_eqn_euler_implicit_pen
    type(atl_equations_type),intent(in) :: equation
    type(atl_statedata_type), intent(inout) :: statedata
    type(atl_material_type), intent(in) :: material
    type(ply_poly_project_type), intent(inout) :: poly_proj
    real(kind=rk), intent(in) :: dt
    type(atl_timestep_type), intent(inout) :: timestep_list
    real(kind=rk), intent(in) :: weight
    type(atl_cube_elem_type), intent(in) :: mesh
    ! -------------------------------------------------------------------- !
    integer :: iPoint, iElem, iMatElem, datIndex
    integer :: nElems
    real(kind=rk), allocatable :: modalCoeff(:,:), modalCoeff_cur(:,:)
    real(kind=rk), allocatable :: pointVal(:,:), cur(:,:)
    real(kind=rk), allocatable :: chebPhysCoord(:,:), scatterNodal(:,:)
    real(kind=rk) :: factor
    real(kind=rk) :: sigbyeps
    ! -------------------------------------------------------------------- !

    datIndex = timestep_list%timestepInfoInteger(1)

    select case(equation%eq_kind)
    case('pec_maxwell_scatter')
       !@todo Maxwell scattered equation system needs to be integrated again,
       !!     and the code from here should move into its hlp module.
       allocate(modalCoeff(poly_proj%body_3d%oversamp_dofs,3))
       allocate(pointVal(poly_proj%body_3d%nquadpoints,3))
       allocate(cur(poly_proj%body_3d%nquadpoints,3))
       allocate(modalCoeff_cur(poly_proj%body_3d%nquadpoints,3))
       allocate(chebPhysCoord(poly_proj%body_3d%nquadpoints,3))
       allocate(scatterNodal(poly_proj%body_3d%nquadpoints,3))

      nElems = material%material_desc%computeElems(1)%nElems

      !@todo Scatterfields need to be integrated again
      ! scatter_mag = sqrt(sum(equation%maxwell%scatter_k**2)) &
      !   &           * statedata%local_time%sim

      const3DElemsScatter: do iMatElem=1,nElems
        iElem = material%material_desc%computeElems(1)%totElemIndices(iMatElem)

        if ( material%material_dat%elemMaterialData(1) &
          & %materialDat(iMatElem,1,3) .fne. 0.0_rk ) then

          call atl_refToPhysCoord(                           &
            &    refpoints  = poly_proj%body_3d%nodes,       &
            &    nPoints    = poly_proj%body_3d%nquadpoints, &
            &    baryCoord  = mesh%bary_coord(iElem, :),     &
            &    elemLength = mesh%length,                   &
            &    physPoints = chebPhysCoord                  )

          !@todo Scatterfields need to be integrated again
          ! scatterNodal(:,1) = cos( scatter_mag                     &
          !   &                      - equation%maxwell%scatter_k(1) &
          !   &                        * chebPhysCoord(:,1)          &
          !   &                      - equation%maxwell%scatter_k(2) &
          !   &                        * chebPhysCoord(:,2)          &
          !   &                      - equation%maxwell%scatter_k(3) &
          !   &                        * chebPhysCoord(:,3)          )
          ! scatterNodal(:,3) = scatterNodal(:,1) &
          !   &                 * equation%maxwell%scatter_E_ampl(3)
          ! scatterNodal(:,2) = scatterNodal(:,1) &
          !   &                 * equation%maxwell%scatter_E_ampl(2)
          ! scatterNodal(:,1) = scatterNodal(:,1) &
          !   &                 * equation%maxwell%scatter_E_ampl(1)

          call ply_convert2oversample(state = statedata%state(iElem,:,:3), &
            &                         ndim = 3,                            &
            &                         poly_proj = poly_proj,               &
            &                         modalCoeffs = modalCoeff(:,:3)       )

          call ply_poly_project_m2n(me = poly_proj,         &
            &                       dim = 3,                &
            &                       nVars = 3,              &
            &                       nodal_data = pointVal,  &
            &                       modal_data = modalCoeff )

          ! Constant material, evaluate sigma by epsilon just once per element.
          sigbyeps = material%material_dat%elemMaterialData(1)       &
            &                             %materialDat(iMatElem,1,3) &
            &        / material%material_dat%elemMaterialData(1)     &
            &                               %materialDat(iMatElem,1,2)
          factor = 1 + weight*dt*sigbyeps
          factor = 1.0_rk / factor

          do iPoint = 1,poly_proj%body_3D%nquadpoints

            ! compute u_i
            pointVal(iPoint,1:3) = ( pointVal(iPoint,1:3)                &
              &                      - weight * dt * sigbyeps            &
              &                               * scatterNodal(iPoint,1:3) &
              &                    ) * factor

            ! compute g(u_i)
            cur(iPoint,1:3) = (-1.0_rk) * sigbyeps                     &
              &                         * ( pointVal(iPoint,1:3)       &
              &                             + scatterNodal(iPoint,1:3) )
          end do

          call ply_poly_project_n2m(me = poly_proj,         &
            &                       dim = 3,                &
            &                       nVars = 3,              &
            &                       nodal_data = pointVal,  &
            &                       modal_data = modalCoeff )

          ! write u_i
          call ply_convertFromOversample( modalCoeffs = modalCoeff(:,:3),      &
            &                             ndim = 3,                            &
            &                             poly_proj = poly_proj,               &
            &                             state = statedata%state(iElem,:,:3)  )

          ! write g(u_i)
          call ply_poly_project_n2m(me = poly_proj,             &
            &                       dim = 3 ,                   &
            &                       nVars = 3,                  &
            &                       nodal_data = cur,           &
            &                       modal_data = modalCoeff_cur )
          call ply_convertFromOversample(                     &
            &    modalCoeffs = modalCoeff_cur(:,:3),          &
            &    ndim = 3,                                    &
            &    poly_proj = poly_proj,                       &
            &    state = timestep_list%timestepData(datIndex) &
            &                         %state(iElem,:,:3)      )
          timestep_list%timestepData(datIndex)%state(iElem,:,4:) = 0.0_rk

        else

          ! No changes, as sigma=0
          ! u_i unchanged

          ! g(u_i) = 0
          timestep_list%timestepData(datIndex)%state(iElem,:,:) = 0.0_rk

        end if

      end do const3DElemsScatter

      nElems = material%material_desc%computeElems(2)%nElems
      var3DElemsScatter: do iMatElem=1,nElems
        iElem = material%material_desc%computeElems(2)%totElemIndices(iMatElem)

        call atl_refToPhysCoord(                           &
          &    refpoints  = poly_proj%body_3d%nodes,       &
          &    nPoints    = poly_proj%body_3d%nquadpoints, &
          &    baryCoord  = mesh%bary_coord(iElem, :),     &
          &    elemLength = mesh%length,                   &
          &    physPoints = chebPhysCoord                  )

        !@todo Scatterfields need to be integrated again
        !scatterNodal(:,1) = cos(scatter_mag                     &
        !  &                     - equation%maxwell%scatter_k(1) &
        !  &                       * chebPhysCoord(:,1)          &
        !  &                     - equation%maxwell%scatter_k(2) &
        !  &                       * chebPhysCoord(:,2)          &
        !  &                     - equation%maxwell%scatter_k(3) &
        !  &                       * chebPhysCoord(:,3)          )

        call ply_convert2oversample( state = statedata%state(iElem,:,:3), &
          &                          ndim = 3,                            &
          &                          poly_proj = poly_proj,               &
          &                          modalCoeffs = modalCoeff(:,:3)       )

        call ply_poly_project_m2n( me = poly_proj,         &
          &                        dim = 3 ,               &
          &                        nVars = 3,              &
          &                        nodal_data = pointVal,  &
          &                        modal_data = modalCoeff )

        do iPoint = 1,poly_proj%body_3D%nquadpoints
          sigbyeps = material%material_dat%elemMaterialData(2)            &
            &                             %materialDat(iMatElem,iPoint,3) &
            &        / material%material_dat%elemMaterialData(2)          &
            &                               %materialDat(iMatElem,iPoint,2)
          factor = 1 + weight*dt*sigbyeps
          factor = 1.0_rk / factor

          ! compute u_i
          pointVal(iPoint,1:3) = ( pointVal(iPoint,1:3)                &
            &                      - weight * dt * sigbyeps            &
            &                               * scatterNodal(iPoint,1:3) &
            &                    ) * factor

          ! compute g(u_i)
          cur(iPoint,1:3) = (-1.0_rk) * sigbyeps                     &
            &                         * ( pointVal(iPoint,1:3)       &
            &                             + scatterNodal(iPoint,1:3) )
        end do

        call ply_poly_project_n2m(me = poly_proj,         &
          &                       dim = 3 ,               &
          &                       nVars = 3,              &
          &                       nodal_data = pointVal,  &
          &                       modal_data = modalCoeff )

        ! write u_i
        call ply_convertFromOversample( modalCoeffs = modalCoeff(:,:3),     &
          &                             ndim = 3,                           &
          &                             poly_proj = poly_proj,              &
          &                             state = statedata%state(iElem,:,:3) )

        ! write g(u_i)
        call ply_poly_project_n2m(me = poly_proj,             &
          &                       dim = 3,                    &
          &                       nVars = 3,                  &
          &                       nodal_data = cur,           &
          &                       modal_data = modalCoeff_cur )
        call ply_convertFromOversample( &
          &    modalCoeffs = modalCoeff_cur(:,:3),                            &
          &    ndim = 3,                                                      &
          &    poly_proj = poly_proj,                                         &
          &    state = timestep_list%timestepData(datIndex)%state(iElem,:,:3) )
        timestep_list%timestepData(datIndex)%state(iElem,:,4:) = 0.0_rk

      end do var3DElemsScatter


    case('maxwell','pec_maxwell','pec_maxwell_pml')
      call atl_eqn_maxwell_implicit_pen(         &
        & material     = material,               &
        & weighted_dt  = weight*dt,              &
        & nDims        = 3,                      &
        & poly_proj    = poly_proj,              &
        & state        = statedata%state,        &
        & timestep_rhs = timestep_list           &
        &                %timestepData(datIndex) &
        &                %state                  )

    case('maxwell_2d','pec_maxwell_2d')
      call atl_eqn_maxwell_implicit_pen(         &
        & material     = material,               &
        & weighted_dt  = weight*dt,              &
        & nDims        = 2,                      &
        & poly_proj    = poly_proj,              &
        & state        = statedata%state,        &
        & timestep_rhs = timestep_list           &
        &                %timestepData(datIndex) &
        &                %state                  )

    case('euler', 'navier_stokes')
      call atl_eqn_euler_implicit_pen( material     = material,               &
        &                              eqn          = equation%euler,         &
        &                              weighted_dt  = weight*dt,              &
        &                              nDims        = 3,                      &
        &                              poly_proj    = poly_proj,              &
        &                              state        = statedata%state,        &
        &                              timestep_rhs = timestep_list           &
        &                                             %timestepData(datIndex) &
        &                                             %state                  )

    case('euler_2d', 'navier_stokes_2d')
      call atl_eqn_euler_implicit_pen( material     = material,               &
        &                              eqn          = equation%euler,         &
        &                              weighted_dt  = weight*dt,              &
        &                              nDims        = 2,                      &
        &                              poly_proj    = poly_proj,              &
        &                              state        = statedata%state,        &
        &                              timestep_rhs = timestep_list           &
        &                                             %timestepData(datIndex) &
        &                                             %state                  )

    case('euler_1d')
      call atl_eqn_euler_implicit_pen( material     = material,               &
        &                              eqn          = equation%euler,         &
        &                              weighted_dt  = weight*dt,              &
        &                              nDims        = 1,                      &
        &                              poly_proj    = poly_proj,              &
        &                              state        = statedata%state,        &
        &                              timestep_rhs = timestep_list           &
        &                                             %timestepData(datIndex) &
        &                                             %state                  )

    case default
      write(*,*) 'IMEX TIMESTEPPING IS NOT IMPLEMENTED FOR ', &
        &        trim(equation%eq_kind)
      write(*,*) ' stopping ...'
      call tem_abort()
    end select

  end subroutine implicit_update

end module atl_imexrk_module
