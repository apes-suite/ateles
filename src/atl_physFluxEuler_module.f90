! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> author: Jens Zudrop
!! Collects all functions related to the physical fluxes of the Euler equations.
module atl_physFluxEuler_module

  use env_module,            only: rk

  implicit none
  private

  public :: atl_physFluxEuler
  public :: atl_physFluxEuler_vec

contains


  !> Physical flux calculation along x direction for Euler equation.
  function atl_physFluxEuler(state, isenCoeff, penalty_char, porosity, U_o) &
    &                       result(physFlux)
    ! ---------------------------------------------------------------------------
    !> The state in nodal space. Dimension is the number of vars, i.e. 5 for Euler
    real(kind=rk), intent(in) :: state(:)
    !> Adiabatice index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff
    !> The value of the characteristic function (stemming from penalization)
    real(kind=rk), intent(in) :: penalty_char
    !> The porosity at the current point
    real(kind=rk), intent(in) :: porosity
    !> Velocity of the obstacle
    real(kind=rk), intent(in) :: U_o
    !> The physical flux along the x axis for all variables
    real(kind=rk) :: physFlux(5)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: pressure, velocity(1:3)
    ! ---------------------------------------------------------------------------

    ! calculate pressure
    pressure = (isenCoeff-1.0_rk) * ( &
                 & state(5) - 0.5_rk*(sum(state(2:4)**2))/state(1) &
                                    & )

    !> @todo JZ: here, we divide by a polynomial, we should be careful! We are leaving
    !! the polynomial space here!
    velocity(1:3) = state(2:4)/state(1)

    ! calculate the nonlinear term for different varibales now.
    ! ... density
    physFlux(1) = state(2) + ((1.0_rk/porosity)-1.0_rk)*penalty_char * (velocity(1)-U_o) * state(1)
    ! ... x-velocity
    physFlux(2) = pressure + state(2)*velocity(1)
    ! ... y-velocity
    physFlux(3) = state(2)*velocity(2)
    ! ... z-velocity
    physFlux(4) = state(2)*velocity(3)
    ! ... total energy
    physFlux(5) = velocity(1) * ( state(5) + pressure )

  end function atl_physFluxEuler

  !> Physical flux calculation along x direction for Euler equation.
  subroutine atl_physFluxEuler_vec(state, isenCoeff, penalty_char, porosity, nPoints, &
    &                              rot, physFlux, U_o)
    ! ---------------------------------------------------------------------------
    !> number of points
    integer, intent(in) :: nPoints
    !> The state in nodal space. Dimension is the number of vars, i.e. 5 for Euler
    real(kind=rk), intent(in) :: state(:,:)
    !> Adiabatice index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff
    !> The value of the characteristic function (stemming from penalization)
    real(kind=rk), intent(in) :: penalty_char(nPoints)
    !> The velocity of the obstacle
    real(kind=rk), intent(in) :: U_o(nPoints)
    !> The porosity at the current point
    real(kind=rk), intent(in) :: porosity
    !> rotation
    integer, intent(in) :: rot(5)
    !> The physical flux along the x axis for all variables
    real(kind=rk), intent(inout) :: physFlux(:,:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: inv_s1, s2, s3, s4, s5, s1
    real(kind=rk) :: pressure, vel1, cPor
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    cPor = (1.0_rk/porosity) - 1.0_rk
    do iPoint = 1, nPoints

      inv_s1 = 1.0_rk / state( iPoint, rot(1) )
      s1 = state( iPoint, rot(1) )
      s2 = state( iPoint, rot(2) )
      s3 = state( iPoint, rot(3) )
      s4 = state( iPoint, rot(4) )
      s5 = state( iPoint, rot(5) )

      ! calculate pressure
      pressure = (isenCoeff-1.0_rk) * ( s5 - 0.5_rk*(s2*s2+s3*s3+s4*s4) * inv_s1 )

      !> @todo JZ: here, we divide by a polynomial, we should be careful! We are leaving
      !! the polynomial space here!
      vel1 = s2 * inv_s1
      ! velocity(2) = s3 * inv_s1
      ! velocity(3) = s4 * inv_s1

      ! calculate the nonlinear term for different varibales now.
      ! ... density
      ! physFlux(1) = (1.0_rk + ((1.0_rk/porosity)-1.0_rk)*penalty_char) * s2
      ! physFlux(3) = s2*velocity(2)
      ! physFlux(4) = s2*velocity(3)
      physFlux( iPoint, rot(1) ) = s2 + cPor * penalty_char(iPoint)*(vel1 - U_o(iPoint))*s1
      physFlux( iPoint, rot(2) ) = pressure + s2*vel1
      physFlux( iPoint, rot(3) ) = s2 * s3 * inv_s1
      physFlux( iPoint, rot(4) ) = s2 * s4 * inv_s1
      physFlux( iPoint, rot(5) ) = vel1 * ( s5 + pressure )
    end do

  end subroutine atl_physFluxEuler_vec

end module atl_physFluxEuler_module
