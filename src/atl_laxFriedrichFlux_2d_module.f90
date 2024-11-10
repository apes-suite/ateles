! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!! Module collects all Lax-Friedrich flux for different types of equations.
module atl_laxFriedrichFlux_2d_module
  ! Treelm modules
  use env_module,                         only: rk
  ! Ateles modules
  use atl_eqn_euler_module,               only: atl_euler_type
  use atl_eqn_LinearEuler_module,         only: atl_LinearEuler_type
  use atl_physFluxEuler_2d_module,        only: atl_physFluxEuler_2d
  use atl_LinearEuler_2d_physFlux_module, only: atl_LinearEuler_2d_physFlux

  implicit none
  private

  public :: atl_laxFriedEuler_2d
  public :: atl_laxFriedLinearEuler_2d


contains


  !> Lax-Friedrich flux (in fully conservative variables) for the 2D Euler
  !! equation
  !!
  !! This interface has to match the abstract definition compute_numflux in the
  !! atl_equation_module.
  subroutine atl_laxFriedEuler_2d(euler, state_left, state_right, &
    &                             material_left, material_right, nPoints, flux)
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
    real(kind=rk) :: pfluxl(4), pfluxr(4)
    integer :: iPoint, matpoint, mm
    ! ---------------------------------------------------------------------------

    isen_coeff = euler%isen_coef
    icm1 = isen_coeff-1.0_rk
    mm = ubound(material_left,1)

    do iPoint = 1, nPoints
      matpoint = min(iPoint, mm)
      volLeft = 1.0_rk/state_left(iPoint,1)
      volRight = 1.0_rk/state_right(iPoint,1)

      pressureLeft = icm1*( state_left(iPoint,4) &
        &                   - 0.5_rk*sum(state_left(iPoint,2:3)**2) &
        &                           * volLeft )
      pressureRight = icm1*( state_right(iPoint,4) &
        &                    - 0.5_rk*sum(state_right(iPoint,2:3)**2) &
        &                            * volRight )

      soundSpeed = max( sqrt( isen_coeff * pressureLeft * volLeft ),  &
        &               sqrt( isen_coeff * pressureRight * volRight ) )

      ! the maximum propagation speed of the waves
      maxSpeed = abs(soundSpeed)                          &
        &      + max( abs(state_left(iPoint,2))*volLeft,  &
        &             abs(state_right(iPoint,2))*volRight )

      pfluxl = atl_physFluxEuler_2d(                          &
        &          state = state_left(iPoint,:),              &
        &          isenCoeff = isen_coeff,                    &
        &          penalty_char = material_left(matpoint,1),  &
        &          U_o = material_left(matpoint,2),           &
        &          porosity = euler%porosity                  )
      pfluxr = atl_physFluxEuler_2d(                          &
        &          state = state_right(iPoint,:),             &
        &          isenCoeff = isen_coeff,                    &
        &          penalty_char = material_right(matpoint,1), &
        &          U_o = material_right(matpoint,2),          &
        &          porosity = euler%porosity                  )

      ! now, calculate the flux
      flux(iPoint,:) = 0.5_rk*( maxSpeed*( state_left(iPoint,:)     &
        &                                 - state_right(iPoint,:) ) &
        &                       + (pfluxl + pfluxr) )
    end do

  end subroutine atl_laxFriedEuler_2d
  ! *****************************************************************************


  ! ****************************************************************************
  !> Lax-Friedrich flux (in fully conservative variables) for the
  !  linear euler equation
  subroutine atl_laxFriedLinearEuler_2d( nSides, nFaceDofs, faceRep, faceFlux, &
    &                                    leftPos, rightPos, var, LinearEuler,  &
    &                                    iDir                                  )
    ! --------------------------------------------------------------------------
    !> Datatype for acoustic equation include all background data
    type(atl_LinearEuler_type), intent(in) :: LinearEuler
    integer, intent(in) :: nFaceDofs, nSides
    real(kind=rk), intent(in) :: faceRep(:, :, :, :)
    real(kind=rk), intent(inout) :: faceFlux(:, :, :, :)
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    integer, intent(in) :: var(:)
    !> Direction of the flow, used for background velocity
    integer, intent(in) :: idir
    ! --------------------------------------------------------------------------
    integer :: iSide, left, right, iDof, iter
    real(kind=rk) :: flux(4)
    real(kind=rk) :: leftstate(4), rightstate(4)
    real(kind=rk) :: max_vel
    ! --------------------------------------------------------------------------

    ! loop over all dofs
    do iter=1,nFaceDofs*nSides
      iDof = (iter-1)/(nSides)+1
      iSide = mod(iter-1,nSides)+1

      ! The position of the left and right element in the state
      ! vector.
      left = leftPos(iSide)
      right = rightPos(iSide)

      leftstate = faceRep(left,iDof,var,2)
      rightstate = faceRep(right,iDof,var,1)

      ! the maximum propagation speed of the waves
      max_vel = sqrt(sum(LinearEuler%velocity_0**2)) &
        &       + abs(LinearEuler%speedofSound)

      ! now, calculate the flux
      flux(:) = (max_vel/2.0_rk)*( leftstate(:) - rightstate(:) )          &
        & + ( atl_LinearEuler_2d_physFlux(leftstate(:), LinearEuler, idir) &
        & + atl_LinearEuler_2d_physFlux(rightstate(:), LinearEuler,idir) ) &
        & / 2.0_rk

      !Assign the same flux for both adjacent elements
      faceFlux(left,iDof,var,2) = flux
      faceFlux(right,iDof,var,1) = flux
    end do

  end subroutine atl_laxFriedLinearEuler_2d
  ! ****************************************************************************


end module atl_laxFriedrichFlux_2d_module
