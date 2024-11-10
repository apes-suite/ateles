! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
!! Module keeping interfaces
module atl_elemental_time_integration_module
  use env_module,            only: rk
  use atl_kerneldata_module, only: atl_kerneldata_type, atl_statedata_type
  use atl_scheme_module,     only: atl_local_timestep_type

  implicit none

  private

  public :: atl_timestep_type
  public :: atl_elemental_timestep
  public :: atl_elemental_timestep_vec
  public :: atl_update_timestep

  !> Timestep information for Euler equation.
  type atl_eulerTimestep_type
    !> Modulus of the maximum velocity (for all elements on the current level).
    real(kind=rk), allocatable :: maxVel(:)
    !> Maximum speed of sound (for all elements on the current level).
    real(kind=rk), allocatable :: speedOfSound(:)
  end type atl_eulerTimestep_type

  !> Timestep information for local linear Euler equation
  type atl_LoclinEulerTimestep_type
    !>mean velocity for all elements
    real(kind=rk), allocatable :: meanVel(:)
    !> mean speed of sound for all elements
    real(kind=rk),allocatable  :: speedOfSound(:)
  end type atl_LoclinEulerTimestep_type

  type atl_timestep_type
    !> Additional storage for timestepping algorithm.
    !! The usage is optional and some solvers do not use is, while others have
    !! to use it (e.g.  multistep/multistage methods).
    !! This field is always passed as an argument to elementwise
    !! timestepping rotine (specified above as a pointer).
    !! The dimension of this array
    !! and its usage is left to the stepping algorithm and is therefore
    !! different for different timestepping schemes.
    type(atl_statedata_type), allocatable :: timestepData(:)

    !> Information about the propagation speeds of the 3D Euler equation.
    type(atl_eulerTimestep_type) :: euler

    !> Information about the propagation speeds of the 2D Euler equation.
    type(atl_eulerTimestep_type) :: euler_2d

    !> Information about the propagation speeds of the 1D Euler equation.
    type(atl_eulerTimestep_type) :: euler_1d

    !> Information about the propagation speed of the local linear Euler
    !! equation
    type(atl_LoclinEulerTimestep_type) :: LoclinEuler

    !> Pointer to a subroutine implementing the elemental procedure of a
    !! timestepping algorithm.
    procedure(atl_elemental_timestep), pointer, nopass :: elemStep => null()
    procedure(atl_elemental_timestep_vec), pointer, nopass &
      & :: elemStep_vec => null()

    !> Storage for additional information of the timestepping method.
    !!
    !! Could be something like the iteration in a multistep scheme.
    !! The exact definition and usage is left to the timestepping scheme.
    integer, allocatable :: timestepInfoInteger(:)
    real(kind=rk), allocatable :: timestepInfoReal(:)

    !> Pointer to bring information from timestep information to the
    !! timestepping scheme itself.
    procedure(atl_update_timestep), pointer, nopass :: updateStep => null()
  end type

  !> Interface definition for elementwise timestepping routine.
  abstract interface
    subroutine atl_elemental_timestep( me, state, cell, dof, sideFlux)
      ! ------------------------------------------------------------------------
      import :: rk, atl_timestep_type
      !> The type of your timestepping.
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
    end subroutine atl_elemental_timestep
  end interface

  !> Interface definition for vectorized elementwise timestepping routine.
  abstract interface
    subroutine atl_elemental_timestep_vec( me, state, kerneldata)
      ! ------------------------------------------------------------------------
      import :: rk, atl_timestep_type, atl_kerneldata_type
      !> The type of your timestepping.
      class(atl_timestep_type), intent(inout) :: me

      !> The state of all cells on this level. This field will be updated
      !! ad the cell position. See kerneldata type for more explanations.
      real(kind=rk), intent(inout) :: state(:,:,:)

      !> Complete kerneldata to get the flux from with additional information.
      type(atl_kerneldata_type), intent(in) :: kerneldata
      ! ------------------------------------------------------------------------
    end subroutine atl_elemental_timestep_vec
  end interface

  !> Interface definition for updating the timestep width.
  abstract interface
    subroutine atl_update_timestep( me, timestepInfo )
      ! ------------------------------------------------------------------------
      import :: rk, atl_timestep_type, atl_local_timestep_type
      !> The type of your timestepping.
      class(atl_timestep_type), intent(inout) :: me

      !> Local timestepping information for that part of the mesh
      type(atl_local_timestep_type), intent(in) :: timestepInfo
      ! ------------------------------------------------------------------------
    end subroutine atl_update_timestep
  end interface

end module atl_elemental_time_integration_module
