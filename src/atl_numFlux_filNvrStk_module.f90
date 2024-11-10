! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

!> Module provides numerical fluxes for the RANS equation.
module atl_numFlux_filNvrStk_module

  ! Treelm modules
  use env_module,                   only: rk
  ! Ateles modules
  use atl_physFluxFilNvrstk_module, only: atl_viscPhysFluxRans_2d
  use atl_eqn_nvrstk_module, only: atl_Navier_stokes_rans_type

  implicit none
  private

  public :: atl_viscRans_2d

contains

  subroutine atl_viscRans_2d(left, left_grad, right, right_grad, isen_coeff, &
    & mu, lambda, thermCond, heatCap, penaltyIP, rans_params,  flux          )
    ! ---------------------------------------------------------------------------
    !> The state on the face from its left limit (in conservative variables).
    real(kind=rk), intent(in) :: left(6)
    !> The gradient state on the face from its left limit (in conservative variables).
    real(kind=rk), intent(in) :: left_grad(6,2)
    !> The state on the face from its right limit (in conservative variables).
    real(kind=rk), intent(in) :: right(6)
    !> The gradient state on the face from its right limit (in conservative variables).
    real(kind=rk), intent(in) :: right_grad(6,2)
    !> The adaiabtic expansion coefficient
    real(kind=rk), intent(in) :: isen_coeff
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
    !> The constants for the Rans eqn
    type(atl_Navier_stokes_rans_type), intent(in) :: rans_params
    real(kind=rk), intent(out) :: flux(6)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: physical_left(6), physical_right(6)
    ! ---------------------------------------------------------------------------

    ! Physical flux on the left face
    physical_left = atl_viscPhysFluxRans_2d( state = left,                      &
                                                  & state_gradient = left_grad, &
                                                  & isenCoeff = isen_coeff,     &
                                                  & mu = mu,                    &
                                                  & lambda = lambda,            &
                                                  & thermCond = thermCond,      &
                                                  & rans_params = rans_params,  &
                                                  & heatCap =heatCap            )

    ! Physical flux on the right face
    physical_right = atl_viscPhysFluxRans_2d( state = right,                    &
                                                  & state_gradient = right_grad,&
                                                  & isenCoeff = isen_coeff,     &
                                                  & mu = mu,                    &
                                                  & lambda = lambda,            &
                                                  & thermCond = thermCond,      &
                                                  & rans_params = rans_params,  &
                                                  & heatCap =heatCap            )

    ! Flux for density
    flux(:) = ((physical_left(:) + physical_right(:))/2.0_rk) - penaltyIP*mu*( left(:) - right(:) )

  end subroutine atl_viscRans_2d


end module atl_numFlux_filNvrStk_module
