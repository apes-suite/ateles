! Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> Testing Piessens algorithm implemented in ply_legser_module with the
!! data provided in his description in [1] Communications of the ACM,
!! January 1974, Volume 17, Number 1, page 25.
!! Used is the expansion of the function f(x) = 1 / (2-x).
program ply_legser_test
  use env_module, only: rk
  use ply_legser_module, only: ply_legser

  implicit none

  real(kind=rk) :: chebymodes(0:20)
  real(kind=rk) :: legmodes(0:20)

  real(kind=rk) :: refleg(0:20)
  real(kind=rk),parameter :: rt3fourth = sqrt(0.75_rk)
  real(kind=rk) :: cterm, ci
  real(kind=rk) :: maxdiff

  integer :: i

  ! Exact Legendre coefficients for expansion of f(x)=1/(2-x)
  ! as provided in [1].
  refleg( 0) = 0.549294e00_rk
  refleg( 1) = 0.295830e00_rk
  refleg( 2) = 0.105917e00_rk
  refleg( 3) = 0.340972e-1_rk
  refleg( 4) = 0.104495e-1_rk
  refleg( 5) = 0.311269e-2_rk
  refleg(10) = 0.601250e-5_rk
  refleg(15) = 0.101339e-7_rk
  refleg(20) = 0.161332e-10_rk

  ! The Chebyshev expansion:
  cterm = 2.0_rk*(1.0_rk - rt3fourth)
  ci = 1.0_rk
  do i=0,20
    chebymodes(i) = ci / rt3fourth
    ci = cterm*ci
  end do

  call ply_legser( A = chebymodes, &
    &              B = legmodes,   &
    &              n = 21          )

  maxdiff = 0.0_rk
  do i=0,5
    maxdiff = max(abs(legmodes(i) - refleg(i)), maxdiff)
    write(*,*) 'Absolute error for mode ', i, ':', &
      &        abs(legmodes(i) - refleg(i))
  end do
  do i=10,20,5
    maxdiff = max(abs(legmodes(i) - refleg(i)), maxdiff)
    write(*,*) 'Absolute error for mode ', i, ':', &
      &        abs(legmodes(i) - refleg(i))
  end do
  if (maxdiff < 0.122e-4_rk) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program ply_legser_test
