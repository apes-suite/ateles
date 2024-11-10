!> This module provides Piessens Algorithm 473 from the Communications of the
!! ACM, January 1974, Volume 17, Number 1.
!!
!! Its unrestricted use in a computer is permitted.
!!
!! This Algorithm transforms a known Chebyshev expansion into a Legendre
!! expansion in N^2 time complexity.
module ply_legser_module
  use env_module, only: rk

  implicit none

  private

  public :: ply_legser

contains

  ! ************************************************************************ !
  !> Subroutine to convert Chebyshev (A) to Legendre (B) coefficients.
  !!
  !! This is Piessens Algorithm, which can be found as Algorithm 473 in the
  !! Communications of the ACM, January 1974, Volume 17, Number 1, page 25.
  !! It is slightly modified to account for Fortran 90 practices.
  !!
  !! Algorithm makes use of the recurrence relation, by which the integral
  !! of the product of the nth Legendre polynomial with the kth Chebyshev
  !! polynomial I_{n,k} is given by:
  !!
  !! I_{n,k+2} = {[(k-1)*k - n(n+1)]*(k+2) / ([(k+3)*(k+2) - n*(n+1)]*k)}
  !!           * I_{n,k}
  !! with I_{0,0} and I_{n,k} = 0 if k < n. For I_{n,n} we have:
  !! I_{n,n} = 2^{2n}*(n!)^2/(2n+1)! if n > 0
  subroutine ply_legser(A, B, n)
    ! -------------------------------------------------------------------- !
    !> Number of coefficients.
    integer, intent(in) :: n

    !> Known coefficients of the Chebyshev approximation.
    real(kind=rk), intent(in) :: A(n)

    !> Computed corresponding coefficients of the Legendre approximation.
    real(kind=rk), intent(out) :: B(n)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ak, al, bb, c, d
    integer :: k, l, ll
    ! -------------------------------------------------------------------- !

    ak = 0.0_rk

    ! Calculation of the first Legendre coefficient
    b(1) = 0.5_rk * a(1)
    do k=3,n,2
      ak = ak + 2.0_rk
      b(1) = b(1) - a(k)/(ak*ak - 1.0_rk)
    end do
    c = 2.0_rk / 3.0_rk
    al = 0.0_rk

    ! Start main loop (remaining Legendre coefficients)
    do l=2,n
      ! Calculation of the Lth coefficient
      ll = l+2
      al = al + 1.0_rk
      bb = c*a(l)
      d = c
      ak = al
      do k=ll,n,2
        d = ( (ak-1.0_rk)*ak - al*(al+1.0_rk) ) &
          &   * (ak+2.0_rk) &
          &   * d &
          &   / ( ((ak+3.0_rk)*(ak+2.0_rk) - al*(al+1.0_rk)) * ak )
        bb = bb + a(k)*d
        ak = ak + 2.0_rk
      end do
      c = 4.0_rk * c * (al+1.0_rk)*(al+1.0_rk) &
        &   / ( (al+al+3.0_rk)*(al+al+2.0_rk) )
      b(l) = (al+0.5_rk)*bb
    end do

  end subroutine ply_legser
  ! ************************************************************************ !

end module ply_legser_module
