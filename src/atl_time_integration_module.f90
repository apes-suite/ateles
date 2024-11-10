! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2013-2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2015, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!! This module keeps the generic routines and datatypes for the timestepping
!! methods.
!!
!! It provides the methods to perform both, single step time integration and
!! multistep in a single framework.
module atl_time_integration_module
  use env_module,                             only: rk

  use tem_general_module,                     only: tem_general_type
  use treelmesh_module,                       only: treelmesh_type

  use atl_cube_container_module,              only: atl_cube_container_type
  use atl_equation_module,                    only: atl_equations_type
  use atl_elemental_time_integration_module,  only: atl_timestep_type
  use ply_poly_project_module,                only: ply_poly_project_type

  implicit none

  private

  public :: atl_global_timestep_type, atl_timestep_control_type,          &
    &       explicitEuler, explicitRungeKutta,                            &
    &       explicitSSPRungeKutta, explicitLocalPredictorGlobalCorrector, &
    &       explicitRungeKuttaTaylor, imexRungeKutta

  !> Explicit euler in time.
  integer, parameter :: explicitEuler = 1

  !> Explicit Runge Kutta in time.
  integer, parameter :: explicitRungeKutta = 2

  !> Explicit Strong-Stability-Preserving Runge Kutta in time.
  integer, parameter :: explicitSSPRungeKutta = 3

  !> Explicit Local-Predictor Global-Corrector Approach in time
  integer, parameter :: explicitLocalPredictorGlobalCorrector = 4

  !> Explicit Runge-Kutta-Taylor in time.
  integer, parameter :: explicitRungeKuttaTaylor = 5

  !> IMEX Runge-Kutta in time.
  integer, parameter :: imexRungeKutta = 6

  !> Interface definition for meshwise timestepping routine.
  !!
  !! This step is used to update the complete mesh.
  abstract interface
    subroutine mesh_timestep( minLevel, maxLevel, currentLevel, cubes, tree, &
      &                       timestep_list, nSteps, equation, general,      &
      &                       commStateTimer, poly_proj_list                 )
      ! ------------------------------------------------------------------------
      import :: rk, atl_equations_type, atl_timestep_type, tem_general_type, &
        &       atl_cube_container_type, ply_poly_project_type, treelmesh_type
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
      type(atl_timestep_type), intent(inout) :: timestep_list(:)

      !> The number of steps of the time stepping scheme
      integer, intent(in) :: nSteps

      !> The equation you are operating with.
      type(atl_equations_type),intent(inout) :: equation

      !> General treelm settings
      type(tem_general_type), intent(inout) :: general

      !> Timer for measuring the communication time inside this routine.
      integer,intent(inout) :: commStateTimer

      !> unique list for projection methods
      type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
      ! ------------------------------------------------------------------------
    end subroutine mesh_timestep
  end interface

  !> Datatype to represent the control of timesteps.
  type atl_timestep_control_type
    !> The cfl number that has to be satisfied by the
    !! timestep (for the convective part of the equation).
    real(kind=rk) :: cfl

    !> The cfl number that has to be satisfied by the
    !! timestep (for the viscous part of the equation).
    real(kind=rk) :: cfl_visc

    !> Flag to indicate whether to use modal min/max estimations in
    !! computation of timestep limitation.
    logical :: use_modal_estimate = .false.

    !> Fixed timestep length
    !!
    !! A negative setting indicates that adaptive timesteps via CFL should
    !! be used.
    !! @note This feature is currently implemented in a minimal change kind of
    !!       approach and still allocates temporary arrays for the timestep
    !!       computation. Those might still be necessary to compute the
    !!       dynamic CFL for the fixed dt, but this is currently also not
    !!       implemented.
    real(kind=rk) :: fixed_dt
  end type


  type atl_global_timestep_type
    !> The type of the timestepping. See the parameters above:
    !!
    !! - explicitEuler
    !! - explicitRungeKutta
    !! - explicitRungeKuttaTaylor
    !! - explicitSSPRungeKutta
    integer :: timestepType = 0

    !> The number of steps in the time-stepping scheme
    integer :: nSteps = 0

    !> The data controlling the timestep sizes.
    type(atl_timestep_control_type) :: control

    !> The final timestep stage.
    procedure(mesh_timestep), pointer, nopass :: meshStep => null()

    !> Additional storage for multisptep/multistage timestepping schemes where
    !! additional coefficients are required.
    !! The usage and the dimensions of this array is left to the solver and is
    !! therefore varying for different timestepping schemes.
    real(kind=rk), allocatable :: timestepCoeff(:)

    !> Datatypes for elementwise timestepping. One for each level.
    type(atl_timestep_type), allocatable :: elementSteps(:)
  end type atl_global_timestep_type

end module atl_time_integration_module
