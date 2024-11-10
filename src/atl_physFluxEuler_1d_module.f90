! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2018 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Physical flux implementation for 1D Euler equations.
!!
!! author: Jens Zudrop
!! Collects all functions related to the physical fluxes of the Euler equations.
module atl_physFluxEuler_1d_module
  ! Treelm modules
  use env_module,            only: rk

  implicit none

  private

  public :: atl_physFluxEuler_1d


contains


  ! ------------------------------------------------------------------------ !
  !> Physical flux calculation along x direction for Euler equation.
  function atl_physFluxEuler_1d(state, isenCoeff, penalty_scaling, U_o) &
    &        result(physFlux)
    ! -------------------------------------------------------------------- !
    !> The state in nodal space. Dimension is the number of vars, i.e. 3 for
    !! Euler 1d
    real(kind=rk), intent(in) :: state(:)

    !> Adiabatic index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff

    !> The scaling of the mass flux due to the penality in the porous
    !! material.
    !!
    !! It is given by $1 + (\frac{1}{\phi} - 1) \Chi$
    !! Thus, if there is no porous media, this factor has to be 1.
    real(kind=rk), intent(in) :: penalty_scaling
    !> Obstacle velocity
    real(kind=rk), intent(in) :: U_o
    !> The physical flux along the x axis for all variables
    real(kind=rk) :: physFlux(3)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: pressure, velocity
    ! -------------------------------------------------------------------- !

    ! calculate pressure
    pressure = (isenCoeff-1.0_rk) &
      &        * ( state(3) - 0.5_rk*(state(2)**2)/state(1) )

    ! the polynomial space here!
    velocity = state(2)/state(1)

    ! calculate the nonlinear term for different varibales now.
    ! ... density
    physFlux(1) = state(2) + penalty_scaling * (velocity-U_o) * state(1)
    ! ... x-velocity
    physFlux(2) = pressure + state(2)*velocity
    ! ... total energy
    physFlux(3) = velocity * ( state(3) + pressure )


  end function atl_physFluxEuler_1d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module atl_physFluxEuler_1d_module
