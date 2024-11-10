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
program ply_lagrange_test
  use env_module, only: rk
  use ply_nodeset_module, only: ply_nodeset_legendre,  &
    &                           ply_nodeset_chebyshev, &
    &                           ply_nodeset_chebyloba
  use ply_lagrange_module, only: ply_lagrange_type,   &
    &                            ply_lagrange_define, &
    &                            ply_lagrange_eval,   &
    &                            ply_lagrange_mode_at

  implicit none

  integer :: iSample
  real(kind=rk), allocatable :: linear(:)
  real(kind=rk), allocatable :: quad(:)
  real(kind=rk) :: sample_x
  real(kind=rk) :: f
  real(kind=rk) :: diff
  real(kind=rk) :: maxdiff
  real(kind=rk) :: lobval(3)
  type(ply_lagrange_type) :: poly
  type(ply_lagrange_type) :: lobpoly

  write(*,*) 'Testing the Lagrange Polynomial module...'

  ! Define a linear polynomial with f(x) = x
  allocate(linear(5))
  linear = ply_nodeset_legendre(5)

  poly = ply_lagrange_define( nPoints = 5,                    &
    &                         nodeset = ply_nodeset_legendre, &
    &                         values  = linear                )

  maxdiff = 0.0_rk
  do iSample=1,20
    call random_number(sample_x)
    f = ply_lagrange_eval(me = poly, x = sample_x)
    diff = abs(f - sample_x)
    maxdiff = max(diff, maxdiff)
  end do

  ! Define a quadric polynomial with f(x) = x*x
  allocate(quad(10))
  quad = ply_nodeset_chebyshev(10)
  quad = quad*quad

  poly = ply_lagrange_define( nPoints = 10,                    &
    &                         nodeset = ply_nodeset_chebyshev, &
    &                         values  = quad                   )

  do iSample=1,40
    call random_number(sample_x)
    f = ply_lagrange_eval(me = poly, x = sample_x)
    diff = abs(f - sample_x*sample_x)
    maxdiff = max(diff, maxdiff)
  end do

  ! Define a reference element with quadratic elements on the
  ! Chebyshev-Lobatto nodes in the interval [-1,1].
  call random_number(lobval) ! The values at the nodes do not matter.
  lobpoly = ply_lagrange_define( nPoints = 3,                     &
    &                            nodeset = ply_nodeset_chebyloba, &
    &                            values  = lobval                 )

  ! The three basis polynomials associated with the respective modes
  ! should be: 0.5*(x*x + x); 1-x*x; 0.5*(x*x - x)
  do iSample=1,60
    call random_number(sample_x)

    f = ply_lagrange_mode_at(me = lobpoly, mode = 1, x = sample_x)
    diff = abs(f - 0.5_rk*(sample_x*sample_x + sample_x))
    maxdiff = max(diff, maxdiff)

    f = ply_lagrange_mode_at(me = lobpoly, mode = 2, x = sample_x)
    diff = abs(f - (1.0_rk - sample_x*sample_x))
    maxdiff = max(diff, maxdiff)

    f = ply_lagrange_mode_at(me = lobpoly, mode = 3, x = sample_x)
    diff = abs(f - 0.5_rk*(sample_x*sample_x - sample_x))
    maxdiff = max(diff, maxdiff)
  end do

  if (maxdiff < 16*epsilon(maxdiff)) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'Found a maximal deviation of ', maxdiff
    write(*,*) 'from the expected value.'
    write(*,*) 'FAILED'
  end if

end program ply_lagrange_test
