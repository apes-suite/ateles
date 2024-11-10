! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2018 Harald Klimach <harald.klimach@uni-siegen.de>
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

module atl_LinearEuler_numflux_module
  ! Treelm modules
  use env_module,         only: rk

  use atl_eqn_LinearEuler_module,      only: atl_LinearEuler_type

  implicit none

  private

  public :: atl_LinearEuler_numflux_subleft
  public :: atl_LinearEuler_numflux_subright
  public :: atl_LinearEuler_numflux_superleft
  public :: atl_LinearEuler_numflux_superright

contains


! ******************************************************************************
subroutine atl_LinearEuler_numflux_subright( nSides, nFaceDofs,  faceRep,      &
    &                                        faceFlux, leftPos, rightPos, var, &
    &                                        LinearEuler, iDir                 )
    ! --------------------------------------------------------------------------
    !> Datatype for LinearEuler equation include all background data
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
    real(kind=rk) :: leftstate(5), rightstate(5)
    !> intermediated state on the interface ( in the * region)
    real(kind=rk) :: dens, velX, velY, velZ, p
    real(kind=rk) :: c, c_inv, c_cube_inv, rho0, rho0_inv, u0
    ! --------------------------------------------------------------------------

    c = LinearEuler%SpeedOfSound
    c_inv = 1 / c
    c_cube_inv = 1 / (c ** 2)
    rho0 = LinearEuler%density_0
    rho0_inv = 1 / rho0
    u0 = LinearEuler%velocity_0(iDir)

    ! loop over all dofs
    do iter=1,nFaceDofs*nSides
      iDof = (iter-1)/(nSides)+1
      iSide = mod(iter-1,nSides)+1

      ! The position of the left and right element in the state
      ! vector.
      left = leftPos(iSide)
      right = rightPos(iSide)

      leftstate(1) = faceRep(left, iDof, var(1), 2)
      leftstate(2) = faceRep(left, iDof, var(2), 2)
      leftstate(3) = faceRep(left, iDof, var(3), 2)
      leftstate(4) = faceRep(left, iDof, var(4), 2)
      leftstate(5) = faceRep(left, iDof, var(5), 2)
      rightstate(1) = faceRep(right, iDof, var(1), 1)
      rightstate(2) = faceRep(right, iDof, var(2), 1)
      rightstate(3) = faceRep(right, iDof, var(3), 1)
      rightstate(4) = faceRep(right, iDof, var(4), 1)
      rightstate(5) = faceRep(right, iDof, var(5), 1)

      ! density
      dens = 0.5 * c_cube_inv * ( rightstate(5) + leftstate(5)   &
        &    + c * rho0 * (leftstate(2) - rightstate(2))) &
        &    + rightstate(1) - rightstate(5) * c_cube_inv
      ! velocity
      velX = 0.5 * c_inv * rho0_inv * (leftstate(5) - rightstate(5) &
        &    + c * rho0 * (rightstate(2) + leftstate(2)))
      ! pressure
      p = 0.5 * (rightstate(5) + leftstate(5) &
        & + c * rho0 * (leftstate(2) - rightstate(2)))

      velY = rightstate(3)
      velZ = rightstate(4)

      ! 1....disturbance in density
      faceFlux(left, iDof, var(1), 2) = u0 * dens + rho0 * velX
      ! 2....disturbance velocity in x direction
      faceFlux(left, iDof, var(2), 2) = u0 * velX + p * rho0_inv
      ! 2....disturbance velocity in y direction
      faceFlux(left, iDof, var(3), 2) = u0 * velY
      ! 2....disturbance velocity in z direction
      faceFlux(left, iDof, var(4), 2) = u0 * velZ
      ! 5....distrubance in pressure
      faceFlux(left, iDof, var(5), 2) = u0 * p + LinearEuler%isen_coef &
        &                               * LinearEuler%pressure_0 * velX

      ! Assign the same flux for both adjacent elements
      faceFlux(right, iDof, var(1), 1) = faceFlux(left, iDof, var(1), 2)
      faceFlux(right, iDof, var(2), 1) = faceFlux(left, iDof, var(2), 2)
      faceFlux(right, iDof, var(3), 1) = faceFlux(left, iDof, var(3), 2)
      faceFlux(right, iDof, var(4), 1) = faceFlux(left, iDof, var(4), 2)
      faceFlux(right, iDof, var(5), 1) = faceFlux(left, iDof, var(5), 2)

    end do

end subroutine atl_LinearEuler_numFlux_subright
! ******************************************************************************


! ******************************************************************************
subroutine atl_LinearEuler_numflux_subleft( nSides, nFaceDofs, faceRep,        &
    &                                       faceFlux, leftPos,  rightPos, var, &
    &                                       LinearEuler, iDir                  )
    ! --------------------------------------------------------------------------
    !> Datatype for LinearEuler equation include all background data
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
    real(kind=rk) :: leftstate(5), rightstate(5)
    !> intermediated state on the interface ( in the * region)
    real(kind=rk) :: dens, velX, velY, velZ, p
    real(kind=rk) :: c, c_inv, c_cube_inv, rho0, rho0_inv, u0
    ! --------------------------------------------------------------------------

    c = LinearEuler%SpeedOfSound
    c_inv = 1 / c
    c_cube_inv = 1 / (c ** 2)
    rho0 = LinearEuler%density_0
    rho0_inv = 1 / rho0
    u0 = LinearEuler%velocity_0(iDir)

    ! loop over all dofs
    do iter=1,nFaceDofs*nSides
      iDof = (iter-1)/(nSides)+1
      iSide = mod(iter-1,nSides)+1

      ! The position of the left and right element in the state
      ! vector.
      left = leftPos(iSide)
      right = rightPos(iSide)

      leftstate(1) = faceRep(left, iDof, var(1), 2)
      leftstate(2) = faceRep(left, iDof, var(2), 2)
      leftstate(3) = faceRep(left, iDof, var(3), 2)
      leftstate(4) = faceRep(left, iDof, var(4), 2)
      leftstate(5) = faceRep(left, iDof, var(5), 2)
      rightstate(1) = faceRep(right, iDof, var(1), 1)
      rightstate(2) = faceRep(right, iDof, var(2), 1)
      rightstate(3) = faceRep(right, iDof, var(3), 1)
      rightstate(4) = faceRep(right, iDof, var(4), 1)
      rightstate(5) = faceRep(right, iDof, var(5), 1)

      ! density
      dens = 0.5 * c_cube_inv * ( rightstate(5) + leftstate(5) &
        &    + c * rho0 * ( leftstate(2) - rightstate(2) ) )   &
        &    + leftstate(1) - leftstate(5) * c_cube_inv
      ! velocity
      velX = 0.5 * c_inv * rho0_inv * ( leftstate(5) - rightstate(5) &
        &    + c * rho0 * ( rightstate(2) + leftstate(2) ) )
      ! pressure
      p = 0.5 * ( rightstate(5) + leftstate(5) &
        & + c * rho0 * ( leftstate(2) - rightstate(2) ) )

      velY = leftstate(3)
      velZ = leftstate(4)

      ! 1....disturbance in density
      faceFlux(left, iDof, var(1), 2) = u0 * dens + rho0 * velX
      ! 2....disturbance velocity in x direction
      faceFlux(left, iDof, var(2), 2) = u0 * velX + p * rho0_inv
      ! 2....disturbance velocity in y direction
      faceFlux(left, iDof, var(3), 2) = u0 * velY
      ! 2....disturbance velocity in z direction
      faceFlux(left, iDof, var(4), 2) = u0 * velZ
      ! 5....distrubance in pressure
      faceFlux(left, iDof, var(5), 2) = u0 * p + LinearEuler%isen_coef &
        &                               * LinearEuler%pressure_0 * velX

      ! Assign the same flux for both adjacent elements
      faceFlux(right, iDof, var(1), 1) = faceFlux(left, iDof, var(1), 2)
      faceFlux(right, iDof, var(2), 1) = faceFlux(left, iDof, var(2), 2)
      faceFlux(right, iDof, var(3), 1) = faceFlux(left, iDof, var(3), 2)
      faceFlux(right, iDof, var(4), 1) = faceFlux(left, iDof, var(4), 2)
      faceFlux(right, iDof, var(5), 1) = faceFlux(left, iDof, var(5), 2)

    end do

end subroutine atl_LinearEuler_numFlux_subleft
! ******************************************************************************


! ******************************************************************************
subroutine atl_LinearEuler_numflux_superright( nSides, nFaceDofs, faceRep,  &
    &                                          faceFlux, leftPos, rightPos, &
    &                                          var, LinearEuler, iDir       )
    ! --------------------------------------------------------------------------
    !> Datatype for LinearEuler equation include all background data
    type(atl_LinearEuler_type), intent(in) :: LinearEuler
    integer, intent(in) :: nFaceDofs, nSides
    real(kind=rk), intent(in) :: faceRep(:, :, :, :)
    real(kind=rk), intent(inout) :: faceFlux(:, :, :, :)
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    integer, intent(in) :: var(:)
    !> Direction of the flow, used for background velocity
    integer, intent(in) :: idir
    ! --------------------------------------------------------------------------
    integer :: iSide, right, left, iDof, iter
    !> intermediated state on the interface ( in the * region)
    real(kind=rk) :: dens, velX, velY, velZ, p
    real(kind=rk) :: rho0, rho0_inv, u0
    ! --------------------------------------------------------------------------

    rho0 = LinearEuler%density_0
    rho0_inv = 1 / rho0
    u0 = LinearEuler%velocity_0(iDir)

    ! loop over all dofs
    do iter=1,nFaceDofs*nSides
      iDof = (iter-1)/(nSides)+1
      iSide = mod(iter-1,nSides)+1

      ! The position of the left and right element in the state
      ! vector.
      right = rightPos(iSide)
      left = leftPos(iSide)

      ! this is the case when all characteristics are going to the left
      ! hence the state is always transport from right to left
      ! which implies the state becomes the state from the right
      dens = faceRep(right, iDof, var(1), 1)
      velX = faceRep(right, iDof, var(2), 1)
      velY = faceRep(right, iDof, var(3), 1)
      velZ = faceRep(right, iDof, var(4), 1)
      p = faceRep(right, iDof, var(5), 1)

      ! 1....disturbance in density
      faceFlux(left, iDof, var(1), 2) = u0 * dens + rho0 * velX
      ! 2....disturbance velocity in x direction
      faceFlux(left, iDof, var(2), 2) = u0 * velX + p * rho0_inv
      ! 2....disturbance velocity in y direction
      faceFlux(left, iDof, var(3), 2) = u0 * velY
      ! 2....disturbance velocity in z direction
      faceFlux(left, iDof, var(4), 2) = u0 * velZ
      ! 5....distrubance in pressure
      faceFlux(left, iDof, var(5), 2) = u0 * p + LinearEuler%isen_coef &
        &                               * LinearEuler%pressure_0 * velX

      ! Assign the same flux for both adjacent elements
      faceFlux(right, iDof, var(1), 1) = faceFlux(left, iDof, var(1), 2)
      faceFlux(right, iDof, var(2), 1) = faceFlux(left, iDof, var(2), 2)
      faceFlux(right, iDof, var(3), 1) = faceFlux(left, iDof, var(3), 2)
      faceFlux(right, iDof, var(4), 1) = faceFlux(left, iDof, var(4), 2)
      faceFlux(right, iDof, var(5), 1) = faceFlux(left, iDof, var(5), 2)

    end do

end subroutine atl_LinearEuler_numFlux_superright
! ******************************************************************************


! ******************************************************************************
subroutine atl_LinearEuler_numflux_superleft( nSides, nFaceDofs, faceRep,  &
    &                                         faceFlux, leftPos, rightPos, &
    &                                         var, LinearEuler, iDir       )
    ! --------------------------------------------------------------------------
    !> Datatype for LinearEuler equation include all background data
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
    !> intermediated state on the interface ( in the * region)
    real(kind=rk) :: dens, velX, velY, velZ, p
    real(kind=rk) :: rho0, rho0_inv, u0
    ! --------------------------------------------------------------------------

    rho0 = LinearEuler%density_0
    rho0_inv = 1 / rho0
    u0 = LinearEuler%velocity_0(iDir)

    ! loop over all dofs
    do iter=1,nFaceDofs*nSides
      iDof = (iter-1)/(nSides)+1
      iSide = mod(iter-1,nSides)+1

      ! The position of the left and right element in the state
      ! vector.
      right = rightPos(iSide)
      left = leftPos(iSide)

      ! this is the case when all characteristics are going to the right
      ! hence the state is always transport from left to right
      ! which implies the state becomes the state from the left
      dens = faceRep(left, iDof, var(1), 2)
      velX = faceRep(left, iDof, var(2), 2)
      velY = faceRep(left, iDof, var(3), 2)
      velZ = faceRep(left, iDof, var(4), 2)
      p = faceRep(left, iDof, var(5), 2)

      ! 1....disturbance in density
      faceFlux(left, iDof, var(1), 2) = u0 * dens + rho0 * velX
      ! 2....disturbance velocity in x direction
      faceFlux(left, iDof, var(2), 2) = u0 * velX + p * rho0_inv
      ! 2....disturbance velocity in y direction
      faceFlux(left, iDof, var(3), 2) = u0 * velY
      ! 2....disturbance velocity in z direction
      faceFlux(left, iDof, var(4), 2) = u0 * velZ
      ! 5....distrubance in pressure
      faceFlux(left, iDof, var(5), 2) = u0 * p + LinearEuler%isen_coef &
        &                               * LinearEuler%pressure_0 * velX

      ! Assign the same flux for both adjacent elements
      faceFlux(right, iDof, var(1), 1) = faceFlux(left, iDof, var(1), 2)
      faceFlux(right, iDof, var(2), 1) = faceFlux(left, iDof, var(2), 2)
      faceFlux(right, iDof, var(3), 1) = faceFlux(left, iDof, var(3), 2)
      faceFlux(right, iDof, var(4), 1) = faceFlux(left, iDof, var(4), 2)
      faceFlux(right, iDof, var(5), 1) = faceFlux(left, iDof, var(5), 2)

    end do

end subroutine atl_LinearEuler_numFlux_superleft
! ******************************************************************************

end module atl_LinearEuler_numflux_module
