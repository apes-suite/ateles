! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2017-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Module for routines and datatypes of Modal Discontinuous Galerkin (MODG)
!! scheme for the LinearEuler equation. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
module atl_modg_1d_LoclinEuler_kernel_module
  use env_module,                       only: rk

  use ply_poly_project_module,          only: ply_poly_project_type
  use atl_equation_module,              only: atl_equations_type

  implicit none

  private

  public ::  atl_modg_1d_LoclinEuler_physFlux


contains


  ! ************************************************************************ !
  !> Calculate the physical flux for the MODG scheme and
  !! Linearized euler equation.
  subroutine atl_modg_1d_LoclinEuler_physFlux(equation, res, state, poly_proj )
    ! -------------------------------------------------------------------- !
    !> The equation system we are working with
    type(atl_equations_type), intent(in) :: equation
    !> The result in the modal form
    real(kind=rk), intent(inout)     :: res(:,:)
    !> The state in the modal form
    real(kind=rk), intent(in), optional :: state(:,:)
    !> The projection datatype for the projection information
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The material information
    ! -------------------------------------------------------------------- !
    ! Loop var for all the dof in an element
    integer :: iDof, nDofs
    ! -------------------------------------------------------------------- !

    real(kind=rk) :: state_0(3)
    real(kind=rk) :: vsq, R_inv, R_inv_sq, isen, isencoef
    real(kind=rk) :: A_21, A_22, A_23
    real(kind=rk) :: A_31, A_32, A_33

    ! get the rotation for the physical flux calculation
    nDofs    = poly_proj%body_1d%ndofs

    state_0  = state(1,:)
    !defining the coefficients for matrix vector multiplication
    isencoef = equation%euler%isen_coef
    vsq      = state_0(2)**2
    R_inv    = 1._rk / state_0(1)
    R_inv_sq = R_inv**2
    isen     = isencoef - 1

    A_21 = R_inv_sq * ((isencoef - 3) / 2) * vsq
    A_22 = (3 - isencoef) * state_0(2) * R_inv
    A_23 = isen

    A_31 = (isen * vsq - isencoef * state_0(3) * state_0(1)) &
      &      * state_0(2)                                    &
      &      * (R_inv**(3))
    A_32 = R_inv_sq                                                        &
      &      * (isencoef * state_0(3) * state_0(1) - (isen / 2) * (3 * vsq))
    A_33 = isencoef * state_0(2) * R_inv

    ! non linear flux for the first mode

    res(1,1) = state_0(2)
    res(1,2) = R_inv * state_0(2)**2                        &
      &          + isen * (state_0(3) - 0.5_rk * R_inv * vsq)
    res(1,3) = R_inv * state_0(2) * state_0(3) * isencoef    &
      &          - 0.5_rk * isen * R_inv_sq * state_0(2) * vsq

    ! linear flux calculation for the higher modes
    dofLoop: do iDof = 2, ndofs

      res(iDof,1) = state(iDof,2)

      res(iDof,2) = A_21 * state(iDof,1)     &
        &             + A_22 * state(iDof,2) &
        &             + A_23 * state(iDof,3)


      res(iDof,3) = A_31 * state(iDof,1)     &
        &             + A_32 * state(iDof,2) &
        &             + A_33 * state(iDof,3)

    end do dofLoop

  end subroutine atl_modg_1d_LoclinEuler_physFlux
  ! ************************************************************************ !


end module atl_modg_1d_LoclinEuler_kernel_module
