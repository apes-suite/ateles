! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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

!> author: Jens Zudrop
!! Routines,functions, datatypes related to Runge Kutta 4 (i.e. classical Runge
!! Kutta) timestepping method.
module atl_rk4_module
  use env_module,                            only: rk
  use tem_aux_module,                        only: tem_abort
  use tem_element_module,                    only: eT_fluid
  use tem_general_module,                    only: tem_general_type
  use tem_logging_module,                    only: logUnit
  use treelmesh_module,                      only: treelmesh_type
  use tem_precice_module,                    only: precice_available

  use atl_elemental_time_integration_module, only: atl_timestep_type
  use atl_cube_elem_module,                  only: atl_cube_elem_type
  use atl_kerneldata_module,                 only: atl_kerneldata_type, &
    &                                              atl_statedata_type
  use atl_scheme_module,                     only: atl_scheme_type,         &
    &                                              atl_local_timestep_type
  use atl_source_types_module,               only: atl_source_type
  use atl_compute_module,                    only: atl_compute_rhs,    &
    &                                              atl_preprocess_rhs, &
    &                                              atl_postprocess_rhs
  use atl_boundary_module,                   only: atl_level_boundary_type
  use atl_facedata_module,                   only: atl_facedata_type
  use atl_time_integration_module,           only: atl_global_timestep_type
  use atl_equation_module,                   only: atl_equations_type
  use atl_bc_header_module,                  only: atl_boundary_type
  use atl_stabilize_module,                  only: atl_stabilize
  use ply_poly_project_module,               only: ply_poly_project_type, &
    &                                              assignment(=)
  use atl_materialPrp_module,                only: atl_material_type
  use atl_writePrecice_module,               only: atl_write_precice
  use atl_cube_container_module,             only: atl_cube_container_type
  use atl_penalization_module,               only: atl_penalizationData_type
!!VK  use atl_timer_module,                      only: atl_timerHandles

  implicit none
  private

  public :: atl_init_explicitRungeKutta

contains

  !> Routine to initialize explicit runge kutta scheme for timestepping.
  subroutine atl_init_explicitRungeKutta( me, minLevel, maxLevel, steps, &
    &                                     statedata_list                 )
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
    type(atl_statedata_type), intent(in) :: statedata_list(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iStep
    ! --------------------------------------------------------------------------
    ! point to euler step function, we dont have to allocate any of the
    ! additional arrays for coefficients or buffering.
    allocate( me%elementSteps(minLevel:maxLevel) )
    select case(steps)
    case(4) ! we have runge kutta with 4 steps, and fourth order
      do iLevel = minLevel, maxLevel
        if(allocated(statedata_list(iLevel)%state)) then
          ! We store one integer, indicating the elemental operation in which
          ! step we are.
          allocate(me%elementSteps(iLevel)%timestepInfoInteger(1))
          ! In a 4 step runge kutta method we have to store 4 intermediate
          ! results.
          allocate(me%elementSteps(iLevel)%timestepData(4))
          do iStep = 1, 4
            allocate(me%elementSteps(iLevel)                          &
              &        %timestepData(iStep)                           &
              &         %state(size(statedata_list(iLevel)%state, 1), &
              &                size(statedata_list(iLevel)%state, 2), &
              &                size(statedata_list(iLevel)%state, 3)) )
          end do
          me%elementSteps(iLevel)%elemStep     => elemental_timestep_rk4
          me%elementSteps(iLevel)%elemStep_vec => elemental_timestep_vec_rk4
          me%elementSteps(iLevel)%updateStep   => update_timestep_rk4
        end if
      end do
      me%meshStep => mesh_timestep_rk4

    case default
      write(logUnit(1),*) 'init runge kutta method for this number of steps' &
        & //' not implemented, stopping...'
      call tem_abort()
    end select
  end subroutine atl_init_explicitRungeKutta



  !> Levelwise updating of runge kutta of order 4
  subroutine update_timestep_rk4( me, timestepInfo )
    ! --------------------------------------------------------------------------
    !> The type of your timestepping.
    class(atl_timestep_type), intent(inout) :: me

    !> Local timestepping information for that part of the mesh
    type(atl_local_timestep_type), intent(in) :: timestepInfo
    ! --------------------------------------------------------------------------
    ! nothing to be done here.
  end subroutine update_timestep_rk4



  !> Elemental operation for timestepping of order 4.
  subroutine elemental_timestep_rk4( me, state, cell, dof, sideFlux)
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    me%timestepData(me%timestepInfoInteger(1))%state(cell,dof,:)       &
      & = me%timestepData(me%timestepInfoInteger(1))%state(cell,dof,:) &
      &   + sideFlux(:)
  end subroutine elemental_timestep_rk4


  !> Elemental operation for timestepping of order 4.
  subroutine elemental_timestep_vec_rk4( me, state, kerneldata )
    ! --------------------------------------------------------------------------
    !> Description of the timestep integration method.
    class(atl_timestep_type), intent(inout) :: me

    !> The state of all cells on this level. This field will be updated
    !! ad the cell position. See kerneldata type for more explanations.
    real(kind=rk), intent(inout) :: state(:,:,:)

    !> Complete kerneldata to get the flux from with additional information.
    type(atl_kerneldata_type), intent(in) :: kerneldata
    ! --------------------------------------------------------------------------

    call compute(kerneldata%nTotal,                                &
      &          kerneldata%nDofs,                                 &
      &          kerneldata%nVars,                                 &
      &          me%timestepData(me%timestepInfoInteger(1))%state, &
      &          kerneldata%state_der                              )

  contains

    ! put it in an own subroutine to allow vectorization with gcc
    subroutine compute( nTotal, nDofs, nScalars, state, state_der )
      ! ------------------------------------------------------------------------
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(inout) :: state(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state_der(:,:,:)
      ! ------------------------------------------------------------------------

      state(:nTotal,:,:) = state(:nTotal,:,:) + state_der(:nTotal,:nDofs,:)

    end subroutine compute

  end subroutine elemental_timestep_vec_rk4


  !> Subroutine for timestepping with explicit runge kutta of order 4.
  !! It contains four Euler forwarding substeps:
  !! 0.1   temp = u0
  !! 1.1   u1   = rk4_substep(temp) = rhs(temp)
  !! 1.2   temp = u0 + dt/2 * u1
  !! 2.1   u2   = rk4_substep(temp) = rhs(temp)
  !! 2.2   temp = u0 + dt/2 * u2
  !! 3.1   u3   = rk4_substep(temp) = rhs(temp)
  !! 3.2   temp = u0 + dt   * u3
  !! 4.1   u4   = rk4_substep(temp) = rhs(temp)
  !! 4.2   u0   = u0 + dt/6 * ( u1 + 2*(u2+u3) + u4 )
  !! After each substep, the results are stabilized
  subroutine mesh_timestep_rk4( minLevel, maxLevel, currentLevel, cubes, tree, &
    &                           timestep_list, nSteps, equation, general,      &
    &                           commStateTimer, poly_proj_list                 )
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
    ! --------------------------------------------------------------------------
    integer :: iLevel, iStep
    type(atl_statedata_type) :: statedata_list_temp(minLevel:maxLevel)
    ! --------------------------------------------------------------------------

    ! Init the arrays with the intermediate results
    do iLevel = minLevel, maxLevel
      ! four arrays which store results of rk4_substep, i.e. u1,u2,u3,u4
      do iStep = 1, 4
        timestep_list(iLevel)%timestepData(iStep)%state(:,:,:) = 0.0_rk
      end do
      ! temp array which stores results of temporal forward steps
      allocate(statedata_list_temp(iLevel)%                               &
        & state(cubes%mesh_list(iLevel)%descriptor%elem%nElems(eT_Fluid), &
        & cubes%scheme_list(iLevel)%nDofs,                                &
        & equation%varSys%nScalars)                                       )
      statedata_list_temp(iLevel)%state = &
        & cubes%statedata_list(iLevel)%state
    end do

    ! Loop 4 RK substeps
    ! Each substep first calculates RHS, then Euler forward, finally stabilize
    do iStep = 1, 4
      call substep( iStep, currentLevel, minLevel, maxLevel, cubes, tree, &
        &           statedata_list_temp, poly_proj_list, timestep_list,   &
        &           equation, general, commStateTimer                     )
    end do

    ! copy temp data back to state
    do iLevel = minLevel, maxLevel
      cubes%statedata_list(iLevel)%state = statedata_list_temp(iLevel)%state
    end do

    ! After the runge kutta FULL step, write the results to precice
    !!VK      call tem_startTimer( timerHandle = atl_timerHandles%preciceWrite )
    if (precice_available) then
      call atl_write_precice( stFunList = equation%stFunList,     &
        &                     varSys    = equation%varSys,        &
        &                     time      = general%simControl%now, &
        &                     tree      = tree                    )
    end if
    !!VK      call tem_stopTimer( timerHandle = atl_timerHandles%preciceWrite )

  contains

    ! This subtouine performs the substep of RK4
    ! It first compute RHS, then compute temporal forward step, finally
    ! stabalize
    subroutine substep( iSubstep, currentLevel, minLevel, maxLevel, cubes, &
      &                 tree, statedata_list_temp, poly_proj_list,         &
      &                 timestep_list, equation, general, commStateTimer   )
      ! ------------------------------------------------------------------------
      !> the substep within the RK4
      integer, intent(in) :: iSubstep

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

      !> Temporary state data
      type(atl_statedata_type), intent(inout) &
        & :: statedata_list_temp(minLevel:maxLevel)

      !> List of levelwise timestepping algorihtms
      type(atl_timestep_type), intent(inout) :: timestep_list(minLevel:)

      !> The equation you are operating with.
      type(atl_equations_type),intent(inout) :: equation

      !> General treelm settings
      type(tem_general_type), intent(inout) :: general

      !> Timer for measuring the communication time inside this routine.
      integer, intent(inout) :: commStateTimer

      !> unique list for projection methods
      type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)

      ! ------------------------------------------------------------------------
      integer :: iLevel
      !> dt factor used in RK4 algorithm
      real(kind=rk) :: rk4_dt(4) = [ 0.5_rk, 0.5_rk, 1.0_rk, 1.0_rk ]
      !> dt factor used to forward local time
      real(kind=rk) :: local_dt(4) = [0.0_rk, 0.5_rk, 0.0_rk, 0.5_rk ]
      ! ------------------------------------------------------------------------

      do iLevel = minLevel, maxLevel
        cubes%statedata_list(iLevel)%local_time%sim                   &
          & = cubes%statedata_list(iLevel)%local_time%sim             &
          &   + local_dt(iSubstep) * cubes%scheme_list(iLevel)%time%dt
        ! update the intermediate time for next rk4 substep, this necessary to
        ! calculate source terms and time-dependent boundary values correctly.
        statedata_list_temp(iLevel)%local_time%sim &
          & = cubes%statedata_list(iLevel)%local_time%sim

        ! Store, in which step of the RK method we are
        timestep_list(iLevel)%timestepInfoInteger(1) = iSubstep
      end do

      !  Compute RHS
      !  timestep_list(iLevel)%timestepData(iSubstep)%state
      !   = rhs( statedata_list_temp )
      call rk4_substep( minLevel              = minLevel,                    &
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

      ! Now, build the intermediate values for the RK4 method. Since the substep
      ! was carried out for all the levels in the previous step, we can continue
      ! in a simple loop over the levels.
      if ( iSubstep == 4 ) then
        do iLevel = minLevel, maxLevel
          ! state = state + (dt/6) * ( state1 + 2*(state2+state3) + state4 )
          call compute_final(                                         &
            & nTotal   = size(cubes%statedata_list(iLevel)%state, 1), &
            & nDofs    = cubes%scheme_list(iLevel)%nDofs,             &
            & nScalars = equation%varSys%nScalars,                    &
            & state    = statedata_list_temp(iLevel)%state,           &
            & state0   = cubes%statedata_list(iLevel)%state,          &
            & state1   = timestep_list(iLevel)%timestepData(1)%state, &
            & state2   = timestep_list(iLevel)%timestepData(2)%state, &
            & state3   = timestep_list(iLevel)%timestepData(3)%state, &
            & state4   = timestep_list(iLevel)%timestepData(4)%state, &
            & dt       = cubes%scheme_list(iLevel)%time%dt            )
        end do
      else ! iSubstep = 1 or 2 or 3
        do iLevel = minLevel, maxLevel
          ! statedata_list_temp = cubes%statedata_list(iLevel)%state
          !   + factor * dt * timestep_list(iLevel)%timestepData(iSubstep)%state
          ! For 1st and 2nd substep, rk4_dt = 0.5
          ! For 3rd         substep, rk4_dt = 1.0
          call compute_intermediate(                                      &
            &    nTotal    = cubes%mesh_list(iLevel)                      &
            &                     %descriptor                             &
            &                     %elem                                   &
            &                     %nElems(eT_fluid),                      &
            &    nDofs     = cubes%scheme_list(iLevel)                    &
            &                     %nDofs,                                 &
            &    nScalars  = equation%varSys                              &
            &                        %nScalars,                           &
            &    state_tmp = statedata_list_temp(iLevel)                  &
            &                  %state,                                    &
            &    state1    = cubes%statedata_list(iLevel)                 &
            &                     %state,                                 &
            &    state2    = timestep_list(iLevel)%timestepData(iSubstep) &
            &                                     %state,                 &
            &    factor    = rk4_dt(iSubstep)                             &
            &                  * cubes%scheme_list(iLevel)                &
            &                         %time                               &
            &                         %dt                                 )
        end do
      end if ! iSubstep == 4

      ! Stabilize the intermediate result of this stage of the RK4 scheme
      call atl_stabilize( minlevel            = minlevel,                  &
        &                 maxlevel            = maxlevel,                  &
        &                 statedata_list      = statedata_list_temp,       &
        &                 statedata_stab_list = cubes%statedata_stab_list, &
        &                 mesh_list           = cubes%mesh_list,           &
        &                 scheme_list         = cubes%scheme_list,         &
        &                 equation            = equation,                  &
        &                 tree                = tree,                      &
        &                 poly_proj_pos       = cubes%poly_proj_pos,       &
        &                 poly_proj_list      = poly_proj_list,            &
        &                 bc                  = cubes%bc,                  &
        &                 boundary            = cubes%boundary_stab_list,  &
        &                 general             = general,                   &
        &                 material_list       = cubes%material_list,       &
        &                 commStateTimer      = commStateTimer             )

    end subroutine substep

    ! own subroutine to allow vectorization with gcc
    subroutine compute_intermediate( nTotal, nDofs, nScalars, state_tmp, &
      &                              state1, state2, factor              )
      ! ------------------------------------------------------------------------
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(out) :: state_tmp(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state1(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state2(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: factor
      ! ------------------------------------------------------------------------

      state_tmp = state1 + factor*state2

    end subroutine compute_intermediate

    ! own subroutine to allow vectorization with gcc
    subroutine compute_final( nTotal, nDofs, nScalars, state, state0, state1, &
      &                       state2, state3, state4, dt                      )
      ! ------------------------------------------------------------------------
      integer, intent(in) :: nTotal, nDofs, nScalars
      real(kind=rk), intent(out) :: state(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state0(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state1(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state2(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state3(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: state4(nTotal, nDofs, nScalars)
      real(kind=rk), intent(in) :: dt
      ! ------------------------------------------------------------------------

      state = state0 + (dt/6.0_rk) &
        & * ( state1 + 2.0_rk * (state2+state3) + state4 )

    end subroutine compute_final

  end subroutine mesh_timestep_rk4


  !> Subroutine calculates a substep of the Runge-Kutta timestepping scheme.
  !! Calls itself recursively for the finer levels until the finest level is
  !! reached.
  recursive subroutine rk4_substep( minLevel, maxLevel, currentLevel,       &
    & mesh_list, tree, levelPointer, kerneldata_list, statedata_list,       &
    & facedata_list, source, penalizationdata_list, boundary_list, bc,      &
    & scheme_list, poly_proj_pos, poly_proj_list, timestep_list, equation,  &
    & material_list, general, commStateTimer                                )
    ! --------------------------------------------------------------------------
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
    !> List of faces states you want to calc the rhs for. For each level we have
    !! one.
    type(atl_facedata_type), intent(inout) :: facedata_list(minLevel:maxLevel)
    !> List of sources, for each level
    type(atl_source_type), intent(inout) :: source
    !> List of penalization data, for each level.
    type(atl_penalizationData_type), intent(inout) :: &
      & penalizationdata_list(minLevel:maxLevel)
    !> List of boundaries, for each level.
    type(atl_level_boundary_type), intent(inout) :: &
      & boundary_list(minLevel:maxLevel)
    !> Global description of the boundaries
    type(atl_boundary_type), intent(in) :: bc(:)
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list(minLevel:maxLevel)
    !> List of position of projection method in unique projection list, for each
    !! level
    integer, intent(in) :: poly_proj_pos(minLevel:maxLevel)
    !> unique list for projection methods
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
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

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
      call rk4_substep( minLevel              = minLevel,              &
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
      &                   kerneldata_list, statedata_list, facedata_list,    &
      &                   source, penalizationdata_list, scheme_list,        &
      &                   poly_proj_pos, poly_proj_list, equation,           &
      &                   material_list, general                             )
    ! ... this call applies the inverse of the mass matrix.
    call atl_postprocess_rhs( mesh       = mesh_list(currentLevel),       &
      &                       kerneldata = kerneldata_list(currentLevel), &
      &                       statedata  = statedata_list(currentLevel),  &
      &                       scheme     = scheme_list(currentLevel),     &
      &                       timestep   = timestep_list(currentLevel),   &
      &                       equation   = equation                       )

  end subroutine rk4_substep

end module atl_rk4_module
