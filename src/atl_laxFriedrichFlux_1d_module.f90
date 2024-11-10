! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!! Module collects all Lax-Friedrich flux for different types of equations.
module atl_laxFriedrichFlux_1d_module
  ! Treelm modules
  use env_module,                   only: rk
  ! Ateles modules
  use atl_physFluxEuler_1d_module,  only: atl_physFluxEuler_1d
  use atl_eqn_euler_module,         only: atl_euler_type

  implicit none
  private

  public :: atl_laxFriedEuler_1d

contains

  !> Lax-Friedrich flux (in fully conservative variables) for the 1D Euler
  !! equation
  !!
  !! This interface has to match the abstract definition compute_numflux in the
  !! atl_equation_module.
  subroutine atl_laxFriedEuler_1d( euler, state_left, state_right,         &
    &                              material_left, material_right, nPoints, &
    &                              flux                                    )
    ! ---------------------------------------------------------------------------
    class(atl_euler_type), intent(in) :: euler
    !> The state on the face from its left limit (in conservative variables).
    real(kind=rk), intent(in) :: state_left(:,:)
    !> The state on the face from its right limit (in conservative variables).
    real(kind=rk), intent(in) :: state_right(:,:)
    !> The left value of the characteristic function (stemming from
    !! penalization)
    real(kind=rk), intent(in) :: material_left(:,:)
    !> The right value of the characteristic function (stemming from
    !! penalization)
    real(kind=rk), intent(in) :: material_right(:,:)
    !> Number of points to evaluate the flux at.
    integer, intent(in) :: nPoints
    !> Resulting flux for the left element (in conservative variables).
    real(kind=rk), intent(out) :: flux(:,:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: maxSpeed, soundSpeed, pressureLeft, pressureRight
    real(kind=rk) :: isen_coeff, icm1
    real(kind=rk) :: volLeft, volRight
    real(kind=rk) :: pfluxl(3), pfluxr(3)
    real(kind=rk) :: penalty_left, penalty_right
    real(kind=rk) :: U_o_left, U_o_right
    integer :: iPoint, matpoint, mm
    ! ---------------------------------------------------------------------------

    isen_coeff = euler%isen_coef
    icm1 = isen_coeff-1.0_rk
    mm = ubound(material_left,1)

    do iPoint = 1, nPoints
      matpoint = min(iPoint, mm)
      penalty_left = material_left(matpoint,1) &
        &              * (1.0_rk/euler%porosity - 1.0_rk)
      penalty_right = material_right(matpoint,1) &
        &               * (1.0_rk/euler%porosity - 1.0_rk)
      U_o_left = material_left(matpoint,2)
      U_o_right = material_right(matpoint,2)
      volLeft = 1.0_rk/state_left(iPoint,1)
      volRight = 1.0_rk/state_right(iPoint,1)
      pressureLeft = icm1*( state_left(iPoint,3)               &
        &                   - 0.5_rk*(state_left(iPoint,2)**2) &
        &                           * volLeft )
      pressureRight = icm1*( state_right(iPoint,3)               &
        &                    - 0.5_rk*(state_right(iPoint,2)**2) &
        &                            * volRight )

      soundSpeed = max( sqrt( isen_coeff * pressureLeft * volLeft ),  &
        &               sqrt( isen_coeff * pressureRight * volRight ) )

      ! the maximum propagation speed of the waves
      maxSpeed = abs(soundSpeed)                          &
        &      + max( abs(state_left(iPoint,2))*volLeft,  &
        &             abs(state_right(iPoint,2))*volRight )

      pfluxl = atl_physFluxEuler_1d( state           = state_left(iPoint,:), &
        &                            isenCoeff       = isen_coeff,           &
        &                            U_o             = U_o_left,             &
        &                            penalty_scaling = penalty_left          )
      pfluxr = atl_physFluxEuler_1d( state           = state_right(iPoint,:), &
        &                            isenCoeff       = isen_coeff,            &
        &                            U_o             = U_o_right,             &
        &                            penalty_scaling = penalty_right          )

      ! now, calculate the flux
      flux(iPoint,:) = 0.5_rk*( maxSpeed*( state_left(iPoint,:)     &
        &                                 - state_right(iPoint,:) ) &
        &                       + (pfluxl + pfluxr) )
    end do

  end subroutine atl_laxFriedEuler_1d

end module atl_laxFriedrichFlux_1d_module
