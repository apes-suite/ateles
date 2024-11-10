! Copyright (c) 2020 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach for University of Siegen.
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
!
!> Collection of node sets to use in the nodal representation of the solution.
!!
!! Node distributions on the reference interval [-1,1] are to be described by
!! functions that satisfy the [[ply_nodeset_coords]] interface and return an
!! array of 1D coordinates for the given number of points.
!!
!! The following sets of nodes are available:
!!
!! * [[ply_nodeset_legendre]] the nodes of the Gauss-Legendre integration
!! * [[ply_nodeset_chebyshev]] the nodes of the Chebyshev integration
!! * [[ply_nodeset_chebyloba]] the nodes of the Chebyshev-Lobatto integration
!!   (includes the interval boundaries)
module ply_nodeset_module
  use env_module, only: rk
  use tem_param_module, only: PI

  implicit none

  private

  interface
    function ply_nodeset_coords( nPoints ) result(x)
      use env_module, only: rk
      !> The number of points to create.
      integer, intent(in) :: nPoints
      !> The coordinates of the nodeset in the interval [-1,1].
      !! The array has to have the length nPoints.
      real(kind=rk) :: x(nPoints)
    end function ply_nodeset_coords
  end interface

  public :: ply_nodeset_coords
  public :: ply_nodeset_legendre
  public :: ply_nodeset_chebyshev
  public :: ply_nodeset_chebyloba


contains


  ! ------------------------------------------------------------------------ !
  !> Compute Gauss-Legendre integration points on the interval [-1,1].
  function ply_nodeset_legendre( nPoints ) result(x)
    ! -------------------------------------------------------------------- !
    !> The number of integration points.
    integer, intent(in) :: nPoints
    !> The coordinates of the Legendre points on the interval [-1,1].
    !! The array has to have the length nPoints.
    real(kind=rk) :: x(nPoints)
    ! -------------------------------------------------------------------- !
    ! some working variables
    real(kind=rk) :: z1, z, pp, p3, p2, p1
    ! Precision to find roots and stop Newton iterations.
    real(kind=rk) :: EPS
    integer :: m, i, j
    ! -------------------------------------------------------------------- !

    EPS = epsilon(z)*16
    m = nPoints/2

    ! Set all symmetric points (center is done seperately for odd nPoints)
    PointLoop: do i=1,m

      z = cos(PI*((i-1)+0.75_rk)/(nPoints+0.5_rk))

      NewtonLoop: do
        p1 = 1.0_rk
        p2 = 0.0_rk
        PolyEval: do j=0,nPoints-1
          p3 = p2
          p2 = p1
          p1 = (real(2*j + 1, kind=rk)*z*p2 - j*p3) / real(j+1, kind=rk)
        end do PolyEval
        pp = nPoints * (z*p1 - p2) / (z*z - 1.0_rk)
        z1 = z
        z = z1 - p1/pp
        if ( abs(z-z1) < EPS ) EXIT NewtonLoop
      end do NewtonLoop

      x(i) = -z
      x(nPoints-i+1) = z

    end do PointLoop

    ! For odd number of points the center is an integration point.
    if (mod(nPoints,2) > 0) then
      i = (nPoints+1)/2
      x(i) = 0.0_rk
    end if

  end function ply_nodeset_legendre
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Generates a given number of Chebyshev points on the unit interval [-1;+1].
  function ply_nodeset_chebyshev( nPoints ) result(x)
    ! -------------------------------------------------------------------- !
    !> The number of points to generate
    integer, intent(in) :: nPoints
    !> The coordinates of the Chebyshev points on the interval [-1,1].
    !! The array has to have the length nPoints.
    real(kind=rk) :: x(nPoints)
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    integer :: iBack
    integer :: nHalf
    ! -------------------------------------------------------------------- !

    nHalf = nPoints/2

    ! Set the middle point for odd number of points.
    ! Will be overwritten for even number of points.
    x(nHalf+1) = 0.0_rk

    do iPoint=1,nHalf
      iBack = nPoints - iPoint + 1
      x(iBack) = cos( PI / nPoints * ( (iPoint - 1) + 0.5_rk ) )
      x(iPoint) = -x(iBack)
    end do

  end function ply_nodeset_chebyshev
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Generates a given number of Chebyshev-Lobatto points on the unit interval
  !! [-1;+1].
  function ply_nodeset_chebyloba( nPoints ) result(x)
    ! -------------------------------------------------------------------- !
    !> The number of points to generate
    integer, intent(in) :: nPoints
    !> The coordinates of the Chebyshev-Lobatto points on the interval [-1,1].
    !! The array has to have the length nPoints.
    real(kind=rk) :: x(nPoints)
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    integer :: iBack
    integer :: nHalf
    ! -------------------------------------------------------------------- !

    nHalf = nPoints/2

    ! Set the middle point for odd number of points.
    ! Will be overwritten for even number of points.
    x(nHalf+1) = 0.0_rk

    ! Both interval limits can be set directly and don't have to be
    ! computed.
    x(1) = 1.0_rk
    x(nPoints) = -1.0_rk

    do iPoint=2,nHalf
      iBack = nPoints - iPoint + 1
      x(iPoint) = cos( (iPoint - 1) * PI / real(nPoints - 1, kind=rk) )
      x(iBack) = -x(iPoint)
    end do

  end function ply_nodeset_chebyloba
  ! ------------------------------------------------------------------------ !

end module ply_nodeset_module
