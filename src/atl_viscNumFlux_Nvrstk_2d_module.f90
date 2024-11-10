! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> author: Jens Zudrop
!! Module provides viscous flux for the compressible Navier-Stokes equation.
module atl_viscNumFlux_Nvrstk_2d_module
  ! Treelm modules
  use env_module,                   only: rk
  ! Ateles modules
  use atl_physFluxNvrstk_2d_module, only: atl_viscPhysFluxNavierStokes_2d

  implicit none
  private

  public :: atl_viscNavierStokes_2d

contains

  subroutine atl_viscNavierStokes_2d(left, left_grad, right, right_grad, &
                                    & mu, lambda, thermCond, heatCap,    &
                                    & penaltyIP, flux)
    ! ---------------------------------------------------------------------------
    !> The state on the face from its left limit (in conservative variables).
    real(kind=rk), intent(in) :: left(4)
    !> The gradient state on the face from its left limit (in conservative variables).
    real(kind=rk), intent(in) :: left_grad(4,2)
    !> The state on the face from its right limit (in conservative variables).
    real(kind=rk), intent(in) :: right(4)
    !> The gradient state on the face from its right limit (in conservative variables).
    real(kind=rk), intent(in) :: right_grad(4,2)
    !> Resulting flux for the left element (in conservative variables).
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The penalty parameter
    real(kind=rk), intent(in) :: penaltyIP
    real(kind=rk), intent(out) :: flux(4)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: physical_left(4), physical_right(4)
    ! ---------------------------------------------------------------------------

    ! Physical flux on the left face
    physical_left = atl_viscPhysFluxNavierStokes_2d( state = left,               &
                                                  & state_gradient = left_grad, &
                                                  & mu = mu,                    &
                                                  & lambda = lambda,            &
                                                  & thermCond = thermCond,      &
                                                  & heatCap =heatCap            )

    ! Physical flux on the right face
    physical_right = atl_viscPhysFluxNavierStokes_2d( state = right,             &
                                                  & state_gradient = right_grad,&
                                                  & mu = mu,                    &
                                                  & lambda = lambda,            &
                                                  & thermCond = thermCond,      &
                                                  & heatCap =heatCap            )

    ! Flux for density
    flux(:) = ((physical_left(:) + physical_right(:))/2.0_rk) - penaltyIP*mu*( left(:) - right(:) )

  end subroutine atl_viscNavierStokes_2d


end module atl_viscNumFlux_Nvrstk_2d_module
