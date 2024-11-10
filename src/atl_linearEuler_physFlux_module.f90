! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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

module atl_LinearEuler_physflux_module
  ! Treelm modules
  use env_module,         only: rk
  ! use tem_aux_module,     only: tem_abort

  use atl_eqn_LinearEuler_module, only: atl_LinearEuler_type

  implicit none

  private

  public :: atl_LinearEuler_physFlux


contains


! ******************************************************************************
  !> Function for physical flux of the LinearEuler equation F, 1D?
  !! Since it is 1d, there need to be passed the correct background velocity (u0
  !! for F - flux in x direction, v0 for G - flux in y direction, w0 for H -
  !! flux in z direction)
  function atl_LinearEuler_physFlux(state, LinearEuler, idir) result(flux)
    ! ---------------------------------------------------------------------------
    !> Datatype for LinearEuler equation include all background data
    type(atl_LinearEuler_type), intent(in) :: LinearEuler
    !> State to compute the fluxes (rho, u, v, w)
    ! the size of array differ for 2d and 3d, hence it is always dimension+1
    real(kind=rk), intent(in) :: state(5)
    !> Direction of flux, used fot  background velocity
    integer, intent(in)  :: iDir
    !> The resulting flux in x direction
    real(kind=rk) :: flux(5)
    ! ---------------------------------------------------------------------------

    ! 1....disturbance in density
    flux(1) = LinearEuler%velocity_0(idir)*state(1)+LinearEuler%density_0*state(2)
    ! 2....disturbance velocity in x direction
    flux(2) = LinearEuler%velocity_0(idir)*state(2)+state(5)/LinearEuler%density_0
    ! 2....disturbance velocity in y direction
    flux(3) = LinearEuler%velocity_0(idir)*state(3)
    ! 2....disturbance velocity in z direction
    flux(4) = LinearEuler%velocity_0(idir)*state(4)
    !5....distrubance in pressure
    flux(5) = LinearEuler%velocity_0(idir)*state(5)+ &
      &       LinearEuler%isen_coef*LinearEuler%pressure_0*state(2)

  end function atl_LinearEuler_physFlux
! ******************************************************************************

end module atl_LinearEuler_physflux_module
