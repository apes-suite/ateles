! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
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

module atl_acoustic_2d_flux_module
  ! Treelm modules
  use env_module,         only: rk

  use atl_eqn_acoustic_module, only: atl_acoustic_type
!!VK  use atl_laxFriedrichFlux_module, only: atl_laxFriedAcoustic

  implicit none

  private

  !> Interface for fluxes of acoustic equations.
  interface atl_acoustic_2d_numflux
    module procedure atl_acoustic_2d_numflux_cube_vec
  end interface atl_acoustic_2d_numflux

  public :: atl_acoustic_2d_numflux
!  public :: atl_acoustic_numflux_oneDir
  public :: atl_acoustic_2d_physFlux


contains


! ******************************************************************************
  !> Function for physical flux of the acoustic equation F, 1D?
  !! Since it is 1d, there need to be passed the correct background velocity (u0
  !! for F - flux in x direction, v0 for G - flux in y direction, w0 for H -
  !! flux in z direction)
  function atl_acoustic_2d_physFlux(state, acoustic, idir) result(flux)
    ! ---------------------------------------------------------------------------
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    !> State to compute the fluxes (rho, u, v, w)
    ! the size of array differ for 2d and 3d, hence it is always dimension+1
    real(kind=rk), intent(in) :: state(acoustic%ndims+1)
    !> Direction of flux, used fot  background velocity
    integer, intent(in)  :: iDir
    !> The resulting flux in x direction
    real(kind=rk) :: flux(acoustic%ndims+1)
    ! ---------------------------------------------------------------------------

    ! 1....disturbance in density
    flux(1) = acoustic%velocity_0(idir)*state(1)+acoustic%density_0*state(2)
    ! 2....disturbance velocity in x direction
    flux(2) = acoustic%velocity_0(idir)*state(2)+(acoustic%speedOfSound**2/ &
      &       acoustic%density_0)*state(1)
    ! 2....disturbance velocity in y direction
    flux(3) = acoustic%velocity_0(idir)*state(3)

  end function atl_acoustic_2d_physFlux
! ******************************************************************************


! ******************************************************************************
  !> summary: calculate flux of pure acoustic  equation directly on the face-vector
  !!
  subroutine atl_acoustic_2d_numflux_cube_vec(nTotalFaces, nSides, nFaceDofs, &
      &                             faceRep, faceFlux, &
      &                             leftPos, rightPos, var, acoustic, iDir )
    ! --------------------------------------------------------------------------
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    integer, intent(in) :: nTotalFaces,  nFaceDofs, nSides
    real(kind=rk), intent(in) :: faceRep(nTotalFaces,nFaceDofs,3,2)
    real(kind=rk), intent(inout) :: faceFlux(nTotalFaces,nFaceDofs,3,2)
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    integer, intent(in) :: var(3)
    !> Direction of the flow, used for background velocity
    integer, intent(in) :: idir
    ! --------------------------------------------------------------------------
    integer :: iSide, left, right, iDof
    real(kind=rk) :: leftstate(3), rightstate(3)
    ! --------------------------------------------------------------------------

    ! AT THE MOMENT ONLY FOR U0 < SPEEDOFSOUND...THUS WE DONT HAVE PURE UPWIND


    ! loop over all dofs
    do iDof = 1, nFaceDofs

      do iSide = 1, nSides

        ! The position of the left and right element in the state
        ! vector.
        left = leftPos(iSide)
        right = rightPos(iSide)

        leftstate = faceRep(left,iDof,var,2)
        rightstate = faceRep(right,iDof,var,1)

!!VK        ! call lax friedrich flux
!!VK        call atl_laxFriedAcoustic(left = leftstate,                 &
!!VK          &                       right = rightstate,               &
!!VK          &                       acoustic = acoustic,              &
!!VK          &                       flux = flux,                      &
!!VK          &                       iDir = idir                       )
!!VK
!!VK        faceFlux(left,iDof,var,2) = flux

        faceFlux(left, iDof,var,2)= &
          &  atl_acoustic_2d_numFlux_oneDir( left= faceRep(left,iDof,var,2),   &
          &                               right= faceRep(right,iDof,var,1), &
          &                               acoustic = acoustic,              &
          &                               idir = idir                       )

        ! Assign the same flux for both adjacent elements
        faceFlux(right,iDof,:,1) = faceFlux(left,iDof,:,2)

      end do

    end do

  end subroutine atl_acoustic_2d_numflux_cube_vec
! ******************************************************************************


! ******************************************************************************
function atl_acoustic_2d_numFlux_oneDir(left,right, acoustic,idir) result(flux)
    ! ---------------------------------------------------------------------------
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    !> The resulting flux in x direction
    real(kind=rk) :: flux(3)
    !> State to compute the fluxes (rho, u, v, w)
    real(kind=rk), intent(in) :: left(3), right(3)
    !> Direction of flux, used fot  background velocity
    integer, intent(in)  :: iDir
    ! ---------------------------------------------------------------------------


    ! the flux for rho
     flux(1) = acoustic%velocity_0(iDir) * 0.5_rk *                            &
     &        ( right(1) + left(1) - acoustic%density_0/acoustic%speedOfSound* &
     &        ( right(2) - left(2) ) )                                         &
     &        + acoustic%density_0 * 0.5_rk *                                  &
     &        ( right(2) + left(2) - acoustic%speedOfSound /acoustic%density_0*&
     &        ( right(1) - left(1) ) )


    ! the flux for velX
     flux(2) = acoustic%velocity_0(iDir) * 0.5_rk *                            &
     &        ( right(2) + left(2) - acoustic%density_0/acoustic%speedOfSound* &
     &        ( right(1) - left(1) ) )                                         &
     &        + acoustic%speedOfSound**2/acoustic%density_0 * 0.5_rk *                                  &
     &        ( right(1) + left(1) - acoustic%density_0/acoustic%SpeedOfSound* &
     &        ( right(2) - left(2) ) )

    ! flux for velY, velZ
    flux(3) = 0.0_rk

end function atl_acoustic_2d_numFlux_oneDir
! ******************************************************************************


end module atl_acoustic_2d_flux_module
