! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> module that holds all routines to calculate the flux for
!! hyperbolic linearzied gas dynamic equations.

module atl_acoustic_physflux_module
  ! Treelm modules
  use env_module,         only: rk

  use atl_eqn_acoustic_module, only: atl_acoustic_type

  implicit none

  private

  public :: atl_acoustic_physFlux


contains


! ******************************************************************************
  !> Function for physical flux of the acoustic equation F, 1D?
  !! Since it is 1d, there need to be passed the correct background velocity (u0
  !! for F - flux in x direction, v0 for G - flux in y direction, w0 for H -
  !! flux in z direction)
  function atl_acoustic_physFlux(state, acoustic, idir) result(flux)
    ! ---------------------------------------------------------------------------
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    !> State to compute the fluxes (rho, u, v, w)
    ! the size of array differ for 2d and 3d, hence it is always dimension+1
    real(kind=rk), intent(in) :: state(4)
    !> Direction of flux, used fot  background velocity
    integer, intent(in)  :: iDir
    !> The resulting flux in x direction
    real(kind=rk) :: flux(4)
    ! ---------------------------------------------------------------------------

    ! 1....disturbance in density
    flux(1) = acoustic%velocity_0(idir)*state(1)+acoustic%density_0*state(2)
    ! 2....disturbance velocity in x direction
    flux(2) = acoustic%velocity_0(idir)*state(2)+(acoustic%speedOfSound**2/ &
      &       acoustic%density_0)*state(1)
    ! 2....disturbance velocity in y direction
    flux(3) = acoustic%velocity_0(idir)*state(3)
    ! 2....disturbance velocity in z direction
    flux(4) = acoustic%velocity_0(idir)*state(4)

  end function atl_acoustic_physFlux
! ******************************************************************************


end module atl_acoustic_physflux_module
