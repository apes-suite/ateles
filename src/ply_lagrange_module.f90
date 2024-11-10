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
!> Lagrange polynomial representation.
!!
!! Lagrange polynomial series by the values at given nodes.
!! The nodes to be used are to be given in form of [[ply_nodeset_coords]].
module ply_lagrange_module
  use env_module, only: rk

  use ply_nodeset_module, only: ply_nodeset_coords

  implicit none

  private

  type ply_lagrange_type
    !> Number of points to represent the Lagrange polynomials
    integer :: nPoints

    !> Coordinates of the points where the nodes are to be found.
    real(kind=rk), allocatable :: coords(:)

    !> Values of the function at all coords.
    real(kind=rk), allocatable :: values(:)
  end type ply_lagrange_type

  public :: ply_lagrange_type
  public :: ply_lagrange_define
  public :: ply_lagrange_eval
  public :: ply_lagrange_mode_at
  public :: ply_lagrange_1D


contains


  ! ------------------------------------------------------------------------ !
  !> Define a new polynomial in the Lagrange basis.
  function ply_lagrange_define(nPoints, nodeset, values) result(me)
    ! -------------------------------------------------------------------- !
    !> Number of points to define the polynomial.
    integer, intent(in) :: nPoints

    !> The set of nodes where the function assumes the given values.
    procedure(ply_nodeset_coords) :: nodeset

    !> Function values at all nPoints of the nodeset.
    real(kind=rk), intent(in) :: values(nPoints)

    !> The newly created Lagrange series describing the polynomial function.
    type(ply_lagrange_type) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    me%nPoints = nPoints
    allocate(me%coords(nPoints))
    me%coords = nodeset(nPoints)

    allocate(me%values(nPoints))
    me%values = values

  end function ply_lagrange_define
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Evaluate a polynomial in the Lagrange basis at some point x.
  function ply_lagrange_eval(me, x) result(f)
    ! -------------------------------------------------------------------- !
    !> The polynomial in Lagrange basis to evaluate at point x.
    type(ply_lagrange_type), intent(in) :: me

    !> Coordinate at which the function is to be evaluated.
    real(kind=rk), intent(in) :: x

    !> Value of the polynomial at coordinate x.
    real(kind=rk) :: f
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    ! -------------------------------------------------------------------- !

    f = me%values(1) * ply_lagrange_mode_at(me = me, mode = 1, x = x)
    do iPoint=2,me%nPoints
      f = f + me%values(iPoint) &
        &     * ply_lagrange_mode_at(me = me, mode = iPoint, x = x)
    end do

  end function ply_lagrange_eval
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Evaluate the given Lagrangian mode (which is 1 at coord(mode) and 0 in
  !! all other points) at a given point x.
  function ply_lagrange_mode_at(me, mode, x) result(f)
    ! -------------------------------------------------------------------- !
    !> The polynomial in Lagrange basis.
    type(ply_lagrange_type), intent(in) :: me

    !> Mode to evaluate at x.
    !!
    !! Here mode identifies the polynomial that is 1 in me%coord(mode) and
    !! 0 in all other nodes.
    integer, intent(in) :: mode

    !> Coordinate at which the mode is to be evaluated.
    real(kind=rk), intent(in) :: x

    !> Value of the polynomial at coordinate x.
    real(kind=rk) :: f
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    real(kind=rk) :: x_m
    ! -------------------------------------------------------------------- !

    f = 1.0_rk
    x_m = me%coords(mode)

    do iPoint=1,mode-1
      f = f * (x - me%coords(iPoint)) &
        &   / (x_m - me%coords(iPoint))
    end do
    do iPoint=mode+1,me%nPoints
      f = f * (x - me%coords(iPoint)) &
        &   / (x_m - me%coords(iPoint))
    end do

  end function ply_lagrange_mode_at
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! Compute each Lagrange mode for every given point.
  function ply_lagrange_1D(me, points) result(pointval)
    ! -------------------------------------------------------------------- !
    !> Definition of the Lagrange polynomial basis to evaluate at points.
    type(ply_lagrange_type), intent(in) :: me

    !> List of points at which the polynomials are to be evaluated.
    real(kind=rk), intent(in) :: points(:)

    !> Resulting Lagrange values at all points.
    !!
    !! First dimension holds the Lagrange modes, second dimension the
    !! points.
    real(kind=rk) :: pointval(me%nPoints, size(points))
    ! -------------------------------------------------------------------- !
    integer :: nPoints
    integer :: iPoint
    integer :: iMode
    ! -------------------------------------------------------------------- !
    nPoints = size(points)
    do iPoint=1,nPoints
      do iMode=1,me%nPoints
        pointval(iMode, iPoint) = ply_lagrange_mode_at(   &
          &                         me   = me,            &
          &                         mode = iMode,         &
          &                         x    = points(iPoint) )
      end do
    end do

  end function ply_lagrange_1D
  ! ------------------------------------------------------------------------ !

end module ply_lagrange_module
