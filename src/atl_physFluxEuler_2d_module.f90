! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
module atl_physFluxEuler_2d_module
  ! Treelm modules
  use env_module,            only: rk

  implicit none
  private


  public :: atl_physFluxEuler_2d

contains

  !> Physical flux calculation along x direction for Euler equation.
  function atl_physFluxEuler_2d(state, isenCoeff, penalty_char, porosity, U_o) result(physFlux)
    ! ------------------------------------------------------------------------------------
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
    real(kind=rk), allocatable :: physFlux(:)
    ! ------------------------------------------------------------------------------------
    real(kind=rk) :: pressure, velocity(1:2)
    ! ------------------------------------------------------------------------------------

    ! allocate the memory for the output
    allocate(physFlux(4))


    ! calculate pressure
    pressure = (isenCoeff-1.0_rk) * ( &
                 & state(4) - 0.5_rk*(sum(state(2:3)**2))/state(1) &
                                    & )

    !> @todo JZ: here, we devide by a polynomial, we should be careful! We are leaving
    !! the polynomial space here!
    velocity(1:2) = state(2:3)/state(1)

    ! calculate the nonlinear term for different varibales now, for the first component consider
    ! also the velocity of the obstacle.
    ! ... density
    physFlux(1) = state(2) + ((1.0_rk/porosity)-1.0_rk)*penalty_char * (velocity(1)-U_o) * state(1)
    ! ... x-velocity
    physFlux(2) = pressure + state(2)*velocity(1)
    ! ... y-velocity
    physFlux(3) = state(2)*velocity(2)
    ! ... total energy
    physFlux(4) = velocity(1) * ( state(4) + pressure )


  end function atl_physFluxEuler_2d

end module atl_physFluxEuler_2d_module
