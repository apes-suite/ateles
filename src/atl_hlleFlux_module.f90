! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011 Daniel Harlacher <daniel.harlacher@uni-siegen.de>
! Copyright (c) 2011-2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2013, 2015, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Module to collecting all data types and subroutines which are related
!! to flux calculations.
module atl_hlleFlux_module
  use env_module,           only: rk
  use atl_eqn_euler_module, only: atl_euler_type

  implicit none

  private

  public :: atl_HLLEuler
  public :: atl_HLLEuler2D
  public :: atl_HLLEuler1D


contains


  !> Calculate the HLL flux given the left an right state.
  !!
  !! @todo Properly treat porosity and penalty terms. Currently they are
  !!       completely ignored!
  subroutine atl_HLLEuler(euler, state_left, state_right, &
    &                     material_left, material_right, nPoints, flux)
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
    real(kind=rk) :: fcl(5),fcr(5), priml(5), primr(5)
    real(kind=rk) :: hl, hr, sqrt_rho_r, sqrt_rho_l !, rho_mean
    real(kind=rk) :: u_mean, v_mean, w_mean, h_mean, c_mean, sr, sl
    real(kind=rk) :: cr, cl
    real(kind=rk) :: isen_coeff
    real(kind=rk) :: icm1
    integer :: iPoint
    !integer :: matpoint
    integer :: mm
    integer :: mr

    isen_coeff = euler%isen_coef
    icm1       = isen_coeff - 1.0_rk
    mm = ubound(material_left,1)
    mr = ubound(material_right,1)

    do iPoint = 1, nPoints
      !matpoint = min(iPoint, mm)
      ! calculate the left primitive states
      priml(1)   = state_left(iPoint,1)
      priml(2) = state_left(iPoint,2) / priml(1)
      priml(3) = state_left(iPoint,3) / priml(1)
      priml(4) = state_left(iPoint,4) / priml(1)
      priml(5)   = icm1*( state_left(iPoint,5)                  &
        &                 - 0.5_rk*(  state_left(iPoint,2)**2   &
        &                           + state_left(iPoint,3)**2   &
        &                           + state_left(iPoint,4)**2 ) &
        &                         / priml(1)                    )

      ! calculate the right primitive states
      primr(1)   = state_right(iPoint,1)
      primr(2) = state_right(iPoint,2) / primr(1)
      primr(3) = state_right(iPoint,3) / primr(1)
      primr(4) = state_right(iPoint,4) / primr(1)
      primr(5)   = icm1*( state_right(iPoint,5)                  &
        &                 - 0.5_rk*(  state_right(iPoint,2)**2   &
        &                           + state_right(iPoint,3)**2   &
        &                           + state_right(iPoint,4)**2 ) &
        &                         / primr(1)                     )

      ! Physical convective fluxes left and right
      fcl(1) = priml(1)*priml(2)
      fcl(2) = priml(1)*priml(2)*priml(2) + priml(5)
      fcl(3) = priml(1)*priml(2)*priml(3)
      fcl(4) = priml(1)*priml(2)*priml(4)
      fcl(5) = priml(2)*(state_left(iPoint,5)+priml(5))
      !
      fcr(1) = primr(1)*primr(2)
      fcr(2) = primr(1)*primr(2)*primr(2) + primr(5)
      fcr(3) = primr(1)*primr(2)*primr(3)
      fcr(4) = primr(1)*primr(2)*primr(4)
      fcr(5) = primr(2)*(state_right(iPoint,5)+primr(5))

      ! speed of sound
      cl   = sqrt(isen_coeff*priml(5)/priml(1))
      cr   = sqrt(isen_coeff*primr(5)/primr(1))

      ! enthalpy
      hl = (state_left(iPoint,5) + priml(5)) / priml(1)
      hr = (state_right(iPoint,5) + primr(5)) / primr(1)

      ! aux. vars
      sqrt_rho_r = sqrt(primr(1))
      sqrt_rho_l = sqrt(priml(1))

      ! mean values
      !rho_mean = sqrt_rho_r*sqrt_rho_l
      u_mean   = (sqrt_rho_r*primr(2) + sqrt_rho_l*priml(2)) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      v_mean   = (sqrt_rho_r*primr(3) + sqrt_rho_l*priml(3)) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      w_mean   = (sqrt_rho_r*primr(4) + sqrt_rho_l*priml(4)) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      h_mean   = (sqrt_rho_r*hr + sqrt_rho_l*hl) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      c_mean   = sqrt((isen_coeff-1.0_rk) &
        &        * (h_mean - 0.5_rk * (u_mean**2 + v_mean**2 + w_mean**2) ))

      sr = max(primr(2)+cr, u_mean+c_mean, 0.0_rk)
      sl = min(priml(2)-cl, u_mean-c_mean, 0.0_rk)

      ! HLLE-Flux
      flux(iPoint,1) = (sr*fcl(1) - sl*fcr(1)) / (sr-sl) &
        &  +  sr*sl / (sr - sl)*(state_right(iPoint,1) - state_left(iPoint,1))
      flux(iPoint,2) = (sr*fcl(2) - sl*fcr(2)) / (sr-sl) &
        &  +  sr*sl / (sr - sl)*(state_right(iPoint,2) - state_left(iPoint,2))
      flux(iPoint,3) = (sr*fcl(3) - sl*fcr(3)) / (sr-sl) &
        &  +  sr*sl / (sr - sl)*(state_right(iPoint,3) - state_left(iPoint,3))
      flux(iPoint,4) = (sr*fcl(4) - sl*fcr(4)) / (sr-sl) &
        &  +  sr*sl / (sr - sl)*(state_right(iPoint,4) - state_left(iPoint,4))
      flux(iPoint,5) = (sr*fcl(5) - sl*fcr(5)) / (sr-sl) &
        &  +  sr*sl / (sr - sl)*(state_right(iPoint,5) - state_left(iPoint,5))
    end do

  end subroutine atl_HLLEuler


  !> Calculate the 2D HLL flux given the left an right state.
  !!
  !! @todo Properly treat porosity and penalty terms. Currently they are
  !!       completely ignored!
  subroutine atl_HLLEuler2D(euler, state_left, state_right, &
    &                       material_left, material_right, nPoints, flux)
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
    real(kind=rk) :: fcl(4),fcr(4), priml(4), primr(4)
    real(kind=rk) :: hl, hr, sqrt_rho_r, sqrt_rho_l !, rho_mean
    real(kind=rk) :: u_mean, v_mean, h_mean, c_mean, sr, sl
    real(kind=rk) :: cr, cl
    real(kind=rk) :: isen_coeff
    real(kind=rk) :: icm1
    integer :: iPoint
    !integer :: matpoint
    integer :: mm
    integer :: mr

    isen_coeff = euler%isen_coef
    icm1       = isen_coeff - 1.0_rk
    mm = ubound(material_left,1)
    mr = ubound(material_right,1)

    do iPoint = 1, nPoints
      !matpoint = min(iPoint, mm)
      ! calculate the left primitive states
      priml(1)   = state_left(iPoint,1)
      priml(2:3) = state_left(iPoint,2:3) / priml(1)
      priml(4)   = icm1*( state_left(iPoint,4)                    &
        &                 - 0.5_rk*sum(state_left(iPoint,2:3)**2) &
        &                         / priml(1)                      )

      ! calculate the right primitive states
      primr(1)   = state_right(iPoint,1)
      primr(2:3) = state_right(iPoint,2:3) / primr(1)
      primr(4)   = icm1*( state_right(iPoint,4)                    &
        &                 - 0.5_rk*sum(state_right(iPoint,2:3)**2) &
        &                         / primr(1)                       )

      ! Physical convective fluxes left and right
      fcl(1) = priml(1)*priml(2)
      fcl(2) = priml(1)*priml(2)*priml(2) + priml(4)
      fcl(3) = priml(1)*priml(2)*priml(3)
      fcl(4) = priml(2)*(state_left(iPoint,4)+priml(4))
      !
      fcr(1) = primr(1)*primr(2)
      fcr(2) = primr(1)*primr(2)*primr(2) + primr(4)
      fcr(3) = primr(1)*primr(2)*primr(3)
      fcr(4) = primr(2)*(state_right(iPoint,4)+primr(4))

      ! speed of sound
      cl   = sqrt(isen_coeff*priml(4)/priml(1))
      cr   = sqrt(isen_coeff*primr(4)/primr(1))

      ! enthalpy
      hl = (state_left(iPoint,4) + priml(4)) / priml(1)
      hr = (state_right(iPoint,4) + primr(4)) / primr(1)

      ! aux. vars
      sqrt_rho_r = sqrt(primr(1))
      sqrt_rho_l = sqrt(priml(1))

      ! mean values
      !rho_mean = sqrt_rho_r*sqrt_rho_l
      u_mean   = (sqrt_rho_r*primr(2) + sqrt_rho_l*priml(2)) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      v_mean   = (sqrt_rho_r*primr(3) + sqrt_rho_l*priml(3)) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      h_mean   = (sqrt_rho_r*hr + sqrt_rho_l*hl) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      c_mean   = sqrt((isen_coeff-1.0_rk) &
        &        * (h_mean - 0.5_rk * (u_mean**2 + v_mean**2) ))

      sr = max(primr(2)+cr, u_mean+c_mean, 0.0_rk)
      sl = min(priml(2)-cl, u_mean-c_mean, 0.0_rk)

      ! HLLE-Flux
      flux(iPoint,:) = (sr*fcl(:) - sl*fcr(:)) / (sr-sl) &
        &  +  sr*sl / (sr - sl)*(state_right(iPoint,:) - state_left(iPoint,:))
    end do

  end subroutine atl_HLLEuler2D


  !> Calculate the 1D HLL flux given the left an right state.
  !!
  !! @todo Properly treat porosity and penalty terms. Currently they are
  !!       completely ignored!
  subroutine atl_HLLEuler1D(euler, state_left, state_right, &
    &                       material_left, material_right, nPoints, flux)
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
    real(kind=rk) :: fcl(3),fcr(3), priml(3), primr(3)
    real(kind=rk) :: hl, hr, sqrt_rho_r, sqrt_rho_l !, rho_mean
    real(kind=rk) :: u_mean, h_mean, c_mean, sr, sl
    real(kind=rk) :: cr, cl
    real(kind=rk) :: isen_coeff
    real(kind=rk) :: icm1
    integer :: iPoint
    !integer :: matpoint
    integer :: mm
    integer :: mr

    isen_coeff = euler%isen_coef
    icm1       = isen_coeff - 1.0_rk
    mm = ubound(material_left,1)
    mr = ubound(material_right,1)

    do iPoint = 1, nPoints
      !matpoint = min(iPoint, mm)
      ! calculate the left primitive states
      priml(1) = state_left(iPoint,1)
      priml(2) = state_left(iPoint,2) / priml(1)
      priml(3) = icm1*( state_left(iPoint,3)               &
        &               - 0.5_rk*(state_left(iPoint,2)**2) &
        &                       / priml(1)                 )

      ! calculate the right primitive states
      primr(1) = state_right(iPoint,1)
      primr(2) = state_right(iPoint,2) / primr(1)
      primr(3) = icm1*( state_right(iPoint,3)               &
        &               - 0.5_rk*(state_right(iPoint,2)**2) &
        &                       / primr(1)                  )

      ! Physical convective fluxes left and right
      fcl(1) = priml(1)*priml(2)
      fcl(2) = priml(1)*priml(2)*priml(2) + priml(3)
      fcl(3) = priml(2)*(state_left(iPoint,3)+priml(3))
      !
      fcr(1) = primr(1)*primr(2)
      fcr(2) = primr(1)*primr(2)*primr(2) + primr(3)
      fcr(3) = primr(2)*(state_right(iPoint,3)+primr(3))

      ! speed of sound
      cl   = sqrt(isen_coeff*priml(3)/priml(1))
      cr   = sqrt(isen_coeff*primr(3)/primr(1))

      ! enthalpy
      hl = (state_left(iPoint,3) + priml(3)) / priml(1)
      hr = (state_right(iPoint,3) + primr(3)) / primr(1)

      ! aux. vars
      sqrt_rho_r = sqrt(primr(1))
      sqrt_rho_l = sqrt(priml(1))

      ! mean values
      !rho_mean = sqrt_rho_r*sqrt_rho_l
      u_mean   = (sqrt_rho_r*primr(2) + sqrt_rho_l*priml(2)) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      h_mean   = (sqrt_rho_r*hr + sqrt_rho_l*hl) &
        &        / (sqrt_rho_r + sqrt_rho_l)
      c_mean   = sqrt((isen_coeff-1.0_rk) &
        &        * (h_mean - 0.5_rk * (u_mean**2) ))

      sr = max(primr(2)+cr, u_mean+c_mean, 0.0_rk)
      sl = min(priml(2)-cl, u_mean-c_mean, 0.0_rk)

      ! HLLE-Flux
      flux(iPoint,:) = (sr*fcl(:) - sl*fcr(:)) / (sr-sl) &
        &  +  sr*sl / (sr - sl)*(state_right(iPoint,:) - state_left(iPoint,:))
    end do

  end subroutine atl_HLLEuler1D

end module atl_hlleFlux_module
