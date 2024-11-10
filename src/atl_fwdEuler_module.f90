! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

!> author: Jens Zudrop
!! Routines, functions and datatypes for Forward Euler (i.e. explicit Euler) timestepping.
module atl_fwdEuler_module
  use env_module,                             only: rk
  use tem_general_module,                     only: tem_general_type
  use treelmesh_module,                       only: treelmesh_type

  use ply_poly_project_module,                only: ply_poly_project_type, &
    &                                               assignment(=)

  use atl_compute_module,                     only: atl_compute_rhs,    &
    &                                               atl_preprocess_rhs, &
    &                                               atl_postprocess_rhs
  use atl_elemental_time_integration_module,  only: atl_timestep_type
  use atl_cube_container_module,              only: atl_cube_container_type
  use atl_scheme_module,                      only: atl_local_timestep_type
  use atl_kerneldata_module,                  only: atl_kerneldata_type
  use atl_time_integration_module,            only: atl_global_timestep_type
  use atl_equation_module,                    only: atl_equations_type
  use atl_stabilize_module,                   only: atl_stabilize

  implicit none
  private

  public :: atl_init_explicitEuler

contains

  !> Initialize explicit euler scheme for timestepping.
  subroutine atl_init_explicitEuler( me, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    !> The datatype to initialize.
    type(atl_global_timestep_type), intent(inout) :: me

    !> The minimal level in the mesh.
    integer, intent(in) :: minLevel

    !> The maximal level in the mesh.
    integer, intent(in) :: maxLevel
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------
    ! Point to euler step function, we do not have to allocate any of the
    ! additional arrays for coefficients or buffering.
    allocate( me%elementSteps(minLevel:maxLevel) )
    do iLevel = minLevel, maxLevel
      allocate(me%elementSteps(iLevel)%timestepInfoReal(1))
      me%elementSteps(iLevel)%elemStep => elemental_timestep_euler
      me%elementSteps(iLevel)%elemStep_vec => elemental_timestep_vec_euler
      me%elementSteps(iLevel)%updateStep => update_timestep_euler
    end do
    me%meshStep => mesh_timestep_euler
  end subroutine atl_init_explicitEuler



  !> Interface definition for levelwise updating of timestepping routine.
  subroutine update_timestep_euler( me, timestepInfo )
    ! ------------------------------------------------------------------------
    !> The type of your timestepping.
    class(atl_timestep_type), intent(inout) :: me

    !> Local timestepping information for that part of the mesh
    type(atl_local_timestep_type), intent(in) :: timestepInfo
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    me%timestepInfoReal(1) = timestepInfo%dt
  end subroutine update_timestep_euler



  !> Interface definition for elementwise timestepping routine.
  subroutine elemental_timestep_euler(me, state, cell, dof, sideFlux)
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

    ! we store the timestep in timestepInfoReal(1)
    state(cell,dof,:) = state(cell,dof,:) + me%timestepInfoReal(1)*sideFlux(:)

  end subroutine elemental_timestep_euler


  !> Interface definition for elementwise timestepping routine.
  subroutine elemental_timestep_vec_euler(me, state, kerneldata)
    ! ------------------------------------------------------------------------
    !> Description of the timestep integration method.
    class(atl_timestep_type), intent(inout) :: me

    !> The state of all cells on this level. This field will be updated
    !! ad the cell position. See kerneldata type for more explanations.
    real(kind=rk), intent(inout) :: state(:,:,:)

    !> Complete kerneldata to get the flux from with additional information.
    type(atl_kerneldata_type), intent(in) :: kerneldata
    ! ------------------------------------------------------------------------

    ! we store the timestep in timestepInfoReal(1)
    state(:kerneldata%nTotal, :, :) =                                &
      & state(:kerneldata%nTotal, :, :)                              &
      & + me%timestepInfoReal(1)                                     &
      & * kerneldata%state_der(:kerneldata%nTotal,:kerneldata%nDofs,:)

  end subroutine elemental_timestep_vec_euler



  !> summary: subroutine for timestepping with explicit euler.
  recursive subroutine mesh_timestep_euler(minLevel, maxLevel, currentLevel,  &
    &                  cubes, tree, timestep_list, nSteps, equation, general, &
    &                  commStateTimer, poly_proj_list                         )
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

    !> The number of steps of the time stepping scheme (assumed to be 1)
    integer, intent(in) :: nSteps

    !> The equation you are operating with.
    type(atl_equations_type),intent(inout) :: equation

    !> General treelm settings
    type(tem_general_type), intent(inout) :: general

    !> Timer for measuring the communication time inside this routine.
    integer, intent(inout) :: commStateTimer

    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    ! --------------------------------------------------------------------------

    ! Before going to the next finer level, we do the following:
    !(0. Create modal representation of source terms)
    ! 1. Project modal representation to the faces.
    ! 2. Communicate the face representations (only the current level).
    ! 3. Interpolate face representations from the current level to
    !    the next finer level.
    ! All these steps are done in the preprocess step of the compute module.
    call atl_preprocess_rhs(                     &
      & minLevel        = minLevel,              &
      & maxLevel        = maxLevel,              &
      & currentLevel    = currentLevel,          &
      & mesh_list       = cubes%mesh_list,       &
      & tree            = tree,                  &
      & statedata_list  = cubes%statedata_list,  &
      & facedata_list   = cubes%facedata_list,   &
      & boundary_list   = cubes%boundary_list,   &
      & bc              = cubes%bc,              &
      & scheme_list     = cubes%scheme_list,     &
      & poly_proj_list  = poly_proj_list,        &
      & equation        = equation,              &
      & material_list   = cubes%material_list,   &
      & general         = general                )
    ! Now, we call the timestep for the next level. I.e. the timestepping
    ! for the finer levels is called recursively unitl we reached the
    ! finest level (i.e. maxLevel) of the mesh.
    if(currentLevel .lt. maxLevel) then
      call mesh_timestep_euler(minLevel, maxLevel, currentLevel+1, cubes, &
        &                      tree,timestep_list, nSteps,                &
        &                      equation, general, commStateTimer,         &
        &                      poly_proj_list                             )
    end if

    ! After making the timestep on the finer level, we do the following:
    ! 1. Interpolate flux from finer level to the current level.
    ! 2. Calculate the fluxes for the current level.
    ! 3. Communicate the fluxes (only the current level).
    ! 4. Calculate the physical fluxes.
    ! 5. Project physical flux + numerical fluxes to test functions.
    ! 6. Multiply with inverse of (cell local) mass matrix.
    ! ... this call includes steps 1 to 5.
    call atl_compute_rhs( minLevel              = minLevel,                    &
      &                   maxLevel              = maxLevel,                    &
      &                   currentLevel          = currentLevel,                &
      &                   mesh_list             = cubes%mesh_list,             &
      &                   tree                  = tree,                        &
      &                   kerneldata_list       = cubes%kerneldata_list,       &
      &                   statedata_list        = cubes%statedata_list,        &
      &                   facedata_list         = cubes%facedata_list,         &
      &                   source                = cubes%source,                &
      &                   penalizationdata_list = cubes%penalizationdata_list, &
      &                   scheme_list           = cubes%scheme_list,           &
      &                   poly_proj_pos         = cubes%poly_proj_pos,         &
      &                   poly_proj_list        = poly_proj_list,              &
      &                   equation              = equation,                    &
      &                   material_list         = cubes%material_list,         &
      &                   general               = general                      )
    ! ... this call applies the inverse of the mass matrix.
    call atl_postprocess_rhs(                             &
      & mesh       = cubes%mesh_list(currentLevel),       &
      & kerneldata = cubes%kerneldata_list(currentLevel), &
      & statedata  = cubes%statedata_list(currentLevel),  &
      & scheme     = cubes%scheme_list(currentLevel),     &
      & timestep   = timestep_list(currentLevel),         &
      & equation   = equation                             )

    ! Update to the new time
    cubes%statedata_list(currentLevel)%local_time%sim &
      &  = cubes%statedata_list(currentLevel)%local_time%sim &
      &    + cubes%scheme_list(currentLevel)%time%dt

    ! Stabilize the final result
    call atl_stabilize( minlevel = minlevel, maxlevel = maxlevel, &
                      & statedata_list = cubes%statedata_list,    &
                      & statedata_stab_list = cubes%statedata_stab_list, &
                      & mesh_list = cubes%mesh_list,              &
                      & scheme_list = cubes%scheme_list,          &
                      & equation = equation,                      &
                      & tree = tree,                              &
                      & poly_proj_pos = cubes%poly_proj_pos,      &
                      & poly_proj_list= poly_proj_list,           &
                      & bc = cubes%bc,                            &
                      & boundary = cubes%boundary_stab_list,      &
                      & general = general,                        &
                      & material_list= cubes%material_list,       &
                      & commStateTimer = commStateTimer           )

  end subroutine mesh_timestep_euler

end module atl_fwdEuler_module
