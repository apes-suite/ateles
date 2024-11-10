! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013, 2015, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015, 2017 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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
module atl_laxFriedrichFlux_module
  ! Treelm modules
  use env_module,                       only: rk
  ! Ateles modules
  use atl_eqn_euler_module,             only: atl_euler_type
  ! use atl_physFluxEuler_module,     only: atl_physFluxEuler
  use atl_eqn_acoustic_module,          only: atl_acoustic_type
  use atl_acoustic_physflux_module,     only: atl_acoustic_physFlux
  use atl_eqn_LinearEuler_module,       only: atl_LinearEuler_type
  use atl_LinearEuler_physflux_module,  only: atl_LinearEuler_physFlux
  use atl_physFluxFilNvrStk_module,     only: atl_PhysFluxRans, &
    &                                         atl_PhysFluxRans_2d
  implicit none
  private

  public :: atl_laxFriedEuler
  public :: atl_laxFriedAcoustic
  public :: atl_laxFriedLinearEuler
  public :: atl_laxFriedRans
  public :: atl_laxFriedRans_2D


contains


  !> Lax-Friedrich flux (in fully conservative variables) for the Euler equation
  !!
  !! This interface has to match the abstract definition compute_numflux in the
  !! atl_equation_module.
  subroutine atl_laxFriedEuler(euler, state_left, state_right, &
    &                          material_left, material_right, nPoints, flux)
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
    real(kind=rk) :: pfluxl(5), pfluxr(5)
    real(kind=rk) :: sl(5), sr(5)
    real(kind=rk) :: porosity, penalty_charL, penalty_charR, U_oL, U_oR
    integer :: iPoint, matpoint, mm
    ! ---------------------------------------------------------------------------

    isen_coeff = euler%isen_coef
    icm1 = isen_coeff-1.0_rk
    mm = ubound(material_left,1)
    porosity = euler%porosity

    !$NEC ivdep
    do iPoint = 1, nPoints

      matpoint = min(iPoint, mm)
      penalty_charL = material_left(matpoint,1)
      penalty_charR = material_right(matpoint,1)
      U_oR = material_right(matpoint,2)
      U_oL = material_left(matpoint,2)

      !$NEC unroll(5)
      sl = state_left(iPoint,1:5)
      !$NEC unroll(5)
      sr = state_right(iPoint,1:5)

      volLeft  = 1.0_rk/sl(1)
      volRight = 1.0_rk/sr(1)

      pressureLeft = icm1*( sl(5) &
        &                   - 0.5_rk*(sl(2)*sl(2)+sl(3)*sl(3)+sl(4)*sl(4))*volLeft )
      pressureRight = icm1*( sr(5) &
        &                   - 0.5_rk*(sr(2)*sr(2)+sr(3)*sr(3)+sr(4)*sr(4))*volRight)

      soundSpeed = max( sqrt( isen_coeff * pressureLeft * volLeft ),  &
        &               sqrt( isen_coeff * pressureRight * volRight ) )

      ! the maximum propagation speed of the waves
      maxSpeed = abs(soundSpeed)                          &
        &      + max( abs(sl(2))*volLeft,  abs(sr(2))*volRight )

      pfluxl(1) = sl(2) + ((1.0_rk/porosity)-1.0_rk)*penalty_charL &
        &                 * (sl(2)*volLeft*U_oL)*sl(1)
      pfluxl(2) = pressureLeft + sl(2)*sl(2)*volLeft
      pfluxl(3) = sl(2)*sl(3)*volLeft
      pfluxl(4) = sl(2)*sl(4)*volLeft
      pfluxl(5) = sl(2)*volLeft * ( sl(5) + pressureLeft )

      pfluxr(1) = sr(2) + ((1.0_rk/porosity)-1.0_rk)*penalty_charR &
        &                 * (sr(2)*volRight*U_oR)*sr(1)
      pfluxr(2) = pressureRight + sr(2)*sr(2)*volRight
      pfluxr(3) = sr(2)*sr(3)*volRight
      pfluxr(4) = sr(2)*sr(4)*volRight
      pfluxr(5) = sr(2)*volRight * ( sr(5) + pressureRight )

      flux(iPoint,1:5) = 0.5_rk*( maxSpeed*( sl(1:5) - sr(1:5) ) &
        &                       + (pfluxl(1:5) + pfluxr(1:5) ) )
    end do

  end subroutine atl_laxFriedEuler
  ! *****************************************************************************


  ! ****************************************************************************
  !> Lax-Friedrich flux (in fully conservative variables) for the
  !  linear euler equation
  subroutine atl_laxFriedLinearEuler( nSides, nFaceDofs, faceRep, faceFlux, &
    &                                 leftPos, rightPos, var, LinearEuler,  &
    &                                 iDir                                  )
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
    real(kind=rk) :: flux(5)
    real(kind=rk) :: leftstate(5), rightstate(5)
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
      flux(:) = (max_vel/2.0_rk)*( leftstate(:) - rightstate(:) )       &
        & + ( atl_LinearEuler_physFlux(leftstate(:), LinearEuler, idir) &
        & + atl_LinearEuler_physFlux(rightstate(:), LinearEuler,idir) ) &
        & / 2.0_rk

      !Assign the same flux for both adjacent elements
      faceFlux(left,iDof,var,2) = flux
      faceFlux(right,iDof,var,1) = flux
    end do

  end subroutine atl_laxFriedLinearEuler
  ! ****************************************************************************


  ! *****************************************************************************
  !> Lax-Friedrich flux (in fully conservative variables) for the Acoustic equation
  subroutine atl_laxFriedAcoustic(left, right, acoustic, flux, iDir)
    ! ---------------------------------------------------------------------------
    !> The state on the face from its left limit (in conservative variables).
    real(kind=rk), intent(in) :: left(4)
    !> The state on the face from its right limit (in conservative variables).
    real(kind=rk), intent(in) :: right(4)
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    !> Resulting flux for the left element (in conservative variables).
    real(kind=rk), intent(out) :: flux(4)
    !> Direction of the flow, used for background velocity
    integer, intent(in) :: idir
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: max_vel
    ! ---------------------------------------------------------------------------

    ! the maximum propagation speed of the waves
    max_vel = sqrt(sum(acoustic%velocity_0**2)) + abs(acoustic%speedofSound)

    ! now, calculate the flux
    flux(:) = (max_vel/2.0_rk)*( left(:) - right(:) )&
      & + ( atl_acoustic_physFlux(left(:), acoustic, idir) &
         & + atl_acoustic_physFlux(right(:), acoustic,idir) )/2.0_rk

  end subroutine atl_laxFriedAcoustic
  ! *****************************************************************************

  !> Lax-Friedrich flux (in fully conservative variables) for the Euler equation
  !!
  !! This interface has to match the abstract definition compute_numflux in the
  !! atl_equation_module.
  subroutine atl_laxFriedRans(euler, state_left, state_right, &
    &                          material_left, material_right, nPoints, flux)
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
    real(kind=rk) :: pfluxl(7), pfluxr(7)
    integer :: iPoint, matpoint, mm
    ! ---------------------------------------------------------------------------

    isen_coeff = euler%isen_coef
    icm1 = isen_coeff-1.0_rk
    mm = ubound(material_left,1)

    do iPoint = 1, nPoints
      matpoint = min(iPoint, mm)
      volLeft = 1.0_rk/state_left(iPoint,1)
      volRight = 1.0_rk/state_right(iPoint,1)

      pressureLeft = icm1*( state_left(iPoint,5) &
        &                   - 0.5_rk*sum(state_left(iPoint,2:4)**2)*volLeft    &
        &                   - state_left(iPoint,6) )
      pressureRight = icm1*( state_right(iPoint,5) &
        &                    - 0.5_rk*sum(state_right(iPoint,2:4)**2)*volRight &
        &                    - state_right(iPoint,6) )

      soundSpeed = max( sqrt( isen_coeff * pressureLeft * volLeft ),  &
        &               sqrt( isen_coeff * pressureRight * volRight ) )

      ! the maximum propagation speed of the waves
      maxSpeed = abs(soundSpeed)                          &
        &      + max( abs(state_left(iPoint,2))*volLeft,  &
        &             abs(state_right(iPoint,2))*volRight )

      pfluxl = atl_physFluxRans( state = state_left(iPoint,:),               &
        &                         isenCoeff = isen_coeff,                    &
        &                         penalty_char = material_left(matpoint,1),  &
        &                         porosity = euler%porosity                  )
      pfluxr = atl_physFluxRans( state = state_right(iPoint,:),              &
        &                         isenCoeff = isen_coeff,                    &
        &                         penalty_char = material_right(matpoint,1), &
        &                         porosity = euler%porosity                  )

      ! now, calculate the flux
      flux(iPoint,:) = 0.5_rk*( maxSpeed*( state_left(iPoint,:)     &
        &                                 - state_right(iPoint,:) ) &
        &                       + (pfluxl + pfluxr) )
    end do

  end subroutine atl_laxFriedRans
  ! *****************************************************************************

  subroutine atl_laxFriedRans_2D(euler, state_left, state_right,        &
    &                          material_left, material_right, nPoints, flux)
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
    real(kind=rk) :: pfluxl(6), pfluxr(6)
    integer :: iPoint, matpoint, mm
    ! ---------------------------------------------------------------------------

    isen_coeff = euler%isen_coef
    icm1 = isen_coeff-1.0_rk
    mm = ubound(material_left,1)

    !$NEC ivdep
    do iPoint = 1, nPoints
      matpoint = min(iPoint, mm)
      volLeft = 1.0_rk/state_left(iPoint,1)
      volRight = 1.0_rk/state_right(iPoint,1)

      pressureLeft = icm1*( state_left(iPoint,4) &
        &                   - 0.5_rk*sum(state_left(iPoint,2:3)**2)*volLeft    &
        &                   - state_left(iPoint,5) )
      pressureRight = icm1*( state_right(iPoint,4) &
        &                    - 0.5_rk*sum(state_right(iPoint,2:3)**2)*volRight &
        &                    - state_right(iPoint,5) )

      soundSpeed = max( sqrt( isen_coeff * pressureLeft * volLeft ),  &
        &               sqrt( isen_coeff * pressureRight * volRight ) )

      ! the maximum propagation speed of the waves
      maxSpeed = abs(soundSpeed)                          &
        &      + max( abs(state_left(iPoint,2))*volLeft,  &
        &             abs(state_right(iPoint,2))*volRight )

      pfluxl = atl_physFluxRans_2d(                              &
        &              state        = state_left(iPoint,:),      &
        &              isenCoeff    = isen_coeff,                &
        &              penalty_char = material_left(matpoint,1), &
        &              porosity     = euler%porosity             )

      pfluxr = atl_physFluxRans_2d(                              &
        &              state        = state_right(iPoint,:),     &
        &              isenCoeff    = isen_coeff,                &
        &              penalty_char = material_right(matpoint,1),&
        &              porosity     = euler%porosity             )

      ! now, calculate the flux
      flux(iPoint,:) = 0.5_rk*( maxSpeed*( state_left(iPoint,:)     &
        &                                 - state_right(iPoint,:) ) &
        &                       + (pfluxl + pfluxr) )
    end do

  end subroutine atl_laxFriedRans_2D

end module atl_laxFriedrichFlux_module
