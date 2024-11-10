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

!> Testing the integration of Legendre Polynomials.
program integrateLeg_test
  use env_module, only: rk

  use ply_modg_basis_module, only: ply_integrateLeg, ply_legendre_1d

  implicit none

  real(kind=rk), allocatable :: integrand(:)
  real(kind=rk), allocatable :: integral(:)
  real(kind=rk), allocatable :: evalatends(:,:)
  real(kind=rk), allocatable :: asserted_half(:)
  real(kind=rk) :: halfint
  real(kind=rk) :: signfact

  integer, parameter :: maxdegree = 10
  integer :: iMode
  integer :: nEven
  logical :: check_ok

  allocate(integrand(maxdegree+1))
  allocate(asserted_half(maxdegree+1))
  allocate(integral(maxdegree+2))
  allocate(evalatends(maxdegree+2,2))

  ! The expected integrals on the half interval for the first Legendre
  ! polynomials.
  ! See
  ! http://math.stackexchange.com/questions/35804/integrating-legendre-polynomials-over-half-range
  ! for a nice summary on this.
  asserted_half = 0.0_rk
  asserted_half(1) = 1.0_rk ! int(l_0) = int(1) = x -> int_0^1 l_0 = 1
  asserted_half(2) = 0.5_rk ! int(l_1) = int(x) = x^2/2 -> int_0^1 l_1 = 1/2

  ! All even modes except the first one yield 0.
  ! For the odd modes the following holds:
  ! int_0^1 l_(2i+1) = -1^i / 2^(2i+1)*(i+1) * (2i over i)
  signfact = 1.0_rk
  nEven = (maxdegree)/2 + mod(maxdegree,2) -1
  do iMode=1,nEven
    signfact = - signfact
    asserted_half(2*iMode+2) = signfact*bico(2*iMode, iMode) &
      &                       / real(2**(2*iMode+1) * (iMode+1), kind=rk)
  end do

  check_ok = .true.

  ! Check integration of Legendre polynomials up to maximal polynomial degree
  ! of maxdegree.
  do iMode=1,maxdegree+1
    integrand = 0.0_rk
    integral = 0.0_rk
    integrand(iMode) = 1.0_rk
    integral(:iMode+1) = ply_integrateLeg( integrand = integrand(:iMode), &
      &                                maxdegree = iMode              )
    evalatends(:iMode+1,:) = ply_legendre_1d( points = [0.0_rk, 1.0_rk], &
      &                                   degree = iMode              )
    halfint = sum( integral(:iMode+1)*(evalatends(:iMode+1,2) &
      &            - evalatends(:iMode+1,1)) )
    write(*,*) halfint, asserted_half(iMode)
    check_ok = check_ok &
      &        .and. (abs(halfint - asserted_half(iMode)) < epsilon(halfint))
  end do

  if (check_ok) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if


contains


  !> Compute the binomial n over k
  elemental function bico(n, k) result(b)
    integer, intent(in) :: n, k
    real(kind=rk) :: b

    integer :: iFact

    b = 1.0_rk
    if ((k >= n-k) .and. (k < n)) then
      do iFact=k+1,n-k
        b = b * real(iFact, kind=rk)
      end do
      do iFact=n-k+1,n
        b = b * (real(iFact, kind=rk) / real(iFact-n+k, kind=rk))
      end do
    else if ((k > 0) .and. (k < n)) then
      do iFact=n-k+1,k
        b = b * real(iFact, kind=rk)
      end do
      do iFact=k+1,n
        b = b * (real(iFact, kind=rk) / real(iFact-k, kind=rk))
      end do
    end if

  end function bico

end program integrateLeg_test
