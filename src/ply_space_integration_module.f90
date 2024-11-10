! Copyright (c) 2011-2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp <v.krupp@grs-sim.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Vyacheslav Korchagin, Melven
! Zoellner and Harald Klimach for German Research School for Simulation
! Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Verena Krupp, Peter Vitt,
! Daniel Petró and Nikhil Anand for University of Siegen.
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

!> Spatial integration with the Gauss-Legendre numerical integration.
module ply_space_integration_module
  use env_module,         only: rk
  use tem_param_module,   only: PI

  implicit none

  private

  public :: ply_gaussLegPoints


contains


  ! ------------------------------------------------------------------------ !
  !> Create Gauss-Legendre integration points and weights for one-dimensional
  !! integration on the interval [x1,x2].
  subroutine ply_gaussLegPoints( x1, x2, x, w, nIntP )
    ! -------------------------------------------------------------------- !
    !> lower limit of integration interval
    real(kind=rk), intent(in) :: x1
    !> upper limit of integration interval
    real(kind=rk), intent(in) :: x2
    !> The coordinates of the gauss points on the interval [x1,x2].
    !! The array has the length nIntP.
    real(kind=rk), intent(out) :: x(:)
    !> The quadrature weights. The array has the length nIntP.
    real(kind=rk), intent(out) :: w(:)
    !> The number of integration points.
    integer, intent(in) :: nIntP
    ! -------------------------------------------------------------------- !
    !> some working variables
    real(kind=rk) :: z1,z,xm,xl,pp,p3,p2,p1;
    !> the relative precision of the points
    real(kind=rk) :: EPS
    integer :: m, i, j
    ! -------------------------------------------------------------------- !

    EPS = 1.0 / (10.0**(PRECISION(1.0_rk)-2) )
    m = (nIntP+1)/2
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)

    do i = 1, m

      z = cos(PI*((i-1)+0.75_rk)/(nIntP+0.5_rk))

      loopToExit: do
        p1 = 1.0_rk
        p2 = 0.0_rk
        do j=0 , nIntP-1
          p3 = p2
          p2 = p1
          p1 = ((2.0_rk*j+1.0_rk)*z*p2-j*p3) / (j+1.0_rk)
        end do
        pp = nIntP*(z*p1-p2)/(z*z-1.0_rk)
        z1 = z
        z = z1-p1/pp
        if ( abs(z-z1) < EPS ) EXIT loopToExit
      end do loopToExit

      x(i) = xm-xl*z
      x(nIntP-i+1) = xm+xl*z
      w(i) = 2.0_rk*xl/((1.0_rk-z*z)*pp*pp)
      w(nIntp-i+1) = w(i)

    end do

  end subroutine ply_gaussLegPoints
  ! ------------------------------------------------------------------------ !

end module ply_space_integration_module
