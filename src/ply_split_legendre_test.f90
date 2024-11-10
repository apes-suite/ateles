! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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

program ply_split_legendre_test_prog
  use env_module, only: rk
  use ply_split_legendre_module, only: ply_split_legendre_test, &
    &                                  ply_split_legendre_matrix
  use ply_modg_basis_module, only: ply_legendre_1d

  implicit none

  integer, parameter :: nModes = 30
  integer, parameter :: nPoints = nModes
  integer, parameter :: nOps = nModes**2
  real(kind=rk), parameter :: tolerance = nOps*epsilon(1.0_rk)

  ! ------------------------------------------------------------------------ !
  logical :: success
  logical :: matching_right
  logical :: matching_left

  integer :: iMode
  integer :: iPoint

  real(kind=rk) :: parentval
  real(kind=rk) :: childval

  ! Polynomial modes
  real(kind=rk) :: parentmode(nModes)
  real(kind=rk) :: childmode(nModes)

  ! Mode values at a set of points
  real(kind=rk) :: legparent(nModes,nPoints)
  real(kind=rk) :: legchild(nModes,nPoints)
  real(kind=rk) :: xi(nPoints)
  real(kind=rk) :: x(nPoints)
  real(kind=rk) :: split_matrix(nModes, nModes)
  ! ------------------------------------------------------------------------ !

  ! Check the basic functionality with the test provided by the split_legendre
  ! module.
  call ply_split_legendre_test(success)


  ! Check whether the polynomials in the subintervals return the same
  ! values as the original polynomial for arbitrary modes and point locations.

  ! First compute the split matrix.
  split_matrix = ply_split_legendre_matrix(nModes)

  ! Create some random parent polynomial.
  call random_number(parentmode)

  ! Create some random sample points in the subinterval.
  call random_number(xi)

  ! Compute all legendre modes at those points.
  legchild = ply_legendre_1D(xi, nModes-1)

  ! Right half:
  ! Do the coordinate transform to obtain the point values in the parent, and
  ! compute all legendre modes for these points.
  x = 0.5_rk*xi + 0.5_rk
  legparent = ply_legendre_1D(x, nModes-1)

  ! Transformation is achieved by multiplication with the upper triangular
  ! matrix. For the right half this is stored in the format (row, column).
  do iMode=1,nModes
    childmode(iMode) = sum( split_matrix(iMode, iMode:) &
      &                     * parentmode(iMode:)        )
  end do

  ! Do the polynomials have the same value at coinciding points?
  matching_right = .true.
  do iPoint=1,nPoints
    parentval = sum(parentmode*legparent(:,iPoint))
    childval = sum(childmode*legchild(:,iPoint))
    matching_right = matching_right &
      &              .and. ( abs(parentval - childval) < tolerance )
  end do

  if (.not. matching_right) write(*,*) &
    &                       'Parent and child do not match in right half!'

  ! Left half:
  ! Do the coordinate transform to obtain the point values in the parent, and
  ! compute all legendre modes for these points.
  x = 0.5_rk*xi - 0.5_rk
  legparent = ply_legendre_1D(x, nModes-1)

  ! Transformation is achieved by multiplication with the lower triangular
  ! matrix. For the left half this is stored in the format (column, row).
  do iMode=1,nModes
    childmode(iMode) = sum( split_matrix(iMode:, iMode) &
      &                     * parentmode(iMode:)        )
  end do

  ! Do the polynomials have the same value at coinciding points?
  matching_left = .true.
  do iPoint=1,nPoints
    parentval = sum(parentmode*legparent(:,iPoint))
    childval = sum(childmode*legchild(:,iPoint))
    matching_left = matching_left &
      &             .and. ( abs(parentval - childval) < tolerance )
  end do

  if (.not. matching_left) write(*,*) &
    &                      'Parent and child do not match in left half!'

  success = success .and. matching_right .and. matching_left

  if (success) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program ply_split_legendre_test_prog
