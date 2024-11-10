! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014, 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

module atl_acoustic_numflux_module
  ! Treelm modules
  use env_module,         only: rk

  use atl_eqn_acoustic_module, only: atl_acoustic_type
!!VK  use atl_laxFriedrichFlux_module, only: atl_laxFriedAcoustic

  implicit none

  private

  !> Interface for fluxes of acoustic equations.
  interface atl_acoustic_numflux
    module procedure atl_acoustic_numflux_cube_vec
  end interface atl_acoustic_numflux

  public :: atl_acoustic_numflux
  public :: atl_acoustic_numflux_oneDir


contains


! ******************************************************************************
  !> summary: calculate flux of pure acoustic  equation directly on the face-vector
  !!
  subroutine atl_acoustic_numflux_cube_vec(nTotalFaces, nSides, nFaceDofs, &
      &                             faceRep, faceFlux, &
      &                             leftPos, rightPos, var, acoustic, iDir )
    ! --------------------------------------------------------------------------
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    integer, intent(in) :: nTotalFaces,  nFaceDofs, nSides
    real(kind=rk), intent(in) :: faceRep(nTotalFaces,nFaceDofs,4,2)
    real(kind=rk), intent(inout) :: faceFlux(nTotalFaces,nFaceDofs,4,2)
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    integer, intent(in) :: var(4)
    !> Direction of the flow, used for background velocity
    integer, intent(in) :: idir
    ! --------------------------------------------------------------------------
    integer :: iSide, left, right, iDof, iter
!!VK    real(kind=rk) :: leftstate(4), rightstate(4)
!!VK    real(kind=rk) :: flux(4)
    ! --------------------------------------------------------------------------

    ! AT THE MOMENT ONLY FOR background velocity is 0, hence U0< SPEEDOFSOUND


    ! loop over all dofs
    do iter=1,nFaceDofs*nSides
      iDof = (iter-1)/(nSides)+1
      iSide = mod(iter-1,nSides)+1

      ! The position of the left and right element in the state
      ! vector.
      left = leftPos(iSide)
      right = rightPos(iSide)

!!VK      leftstate = faceRep(left,iDof,var,2)
!!VK      rightstate = faceRep(right,iDof,var,1)
!!VK
!!VK        ! call lax friedrich flux
!!VK        call atl_laxFriedAcoustic(left = leftstate,                 &
!!VK          &                       right = rightstate,               &
!!VK          &                       acoustic = acoustic,              &
!!VK          &                       flux = flux,                      &
!!VK          &                       iDir = idir                       )
!!VK
!!VK        faceFlux(left,iDof,var,2) = flux

      faceFlux(left, iDof,var,2)= &
        &  atl_acoustic_numFlux_oneDir( left= faceRep(left,iDof,var,2),   &
        &                               right= faceRep(right,iDof,var,1), &
        &                               acoustic = acoustic,              &
        &                               idir = idir                       )

      ! Assign the same flux for both adjacent elements
      faceFlux(right,iDof,:,1) = faceFlux(left,iDof,:,2)

    end do
!!VK !upwind case (at the moment not needed, since U0=0)
!!VK      faceFlux(left, iDof,var,2)=  faceRep(left,iDof,var,2)
!!VK      faceFlux(right,iDof,var,1) = faceFlux(left,iDof,var,2)
!!VK

  end subroutine atl_acoustic_numflux_cube_vec
! ******************************************************************************


! ******************************************************************************
function atl_acoustic_numFlux_oneDir(left,right, acoustic,idir) result(flux)
    ! ---------------------------------------------------------------------------
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    !> The resulting flux in x direction
    real(kind=rk) :: flux(4)
    !> State to compute the fluxes (rho, u, v, w)
    real(kind=rk), intent(in) :: left(4), right(4)
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
    if (acoustic%velocity_0(iDir) >= 0 ) then
      flux(3) = left(3)
      flux(4) = left(4)
    else
      flux(3) = right(3)
      flux(4) = right(4)
    end if


end function atl_acoustic_numFlux_oneDir
! ******************************************************************************

end module atl_acoustic_numflux_module
