! Copyright (c) 2016-2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <Michael.Gaida@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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
module atl_modg_LoclinEuler_kernel_module
  use env_module,                       only: rk

  use ply_poly_project_module,          only: ply_poly_project_type
  use atl_equation_module,              only: atl_equations_type
  use atl_scheme_module,                only: atl_scheme_type
  use atl_penalization_module,          only: atl_penalizationData_type
  use atl_materialPrp_module,           only: atl_material_type

  implicit none
  private

  public :: atl_modg_LoclinEuler_physFlux


contains


  ! ****************************************************************************
  !> Calculate the physical flux for the MODG scheme and
  !! Linearized euler equation.
  subroutine atl_modg_LoclinEuler_physFlux( equation, res, state, iElem, iDir, &
    &                                       penalizationData, poly_proj,       &
    &                                       material, nodal_data,              &
    &                                       nodal_gradData, nodal_res,         &
    &                                       elemLength, scheme_min,            &
    &                                       scheme_current                     )
    ! --------------------------------------------------------------------------
    !> The equation system we are working with
    type(atl_equations_type), intent(in) :: equation
    !> The result in the modal form
    real(kind=rk), intent(inout)     :: res(:,:)
    !> The state in the modal form
    real(kind=rk), intent(in), optional :: state(:,:)
    !> The current element index
    integer, intent(in) :: iElem
    !> The current direction
    integer, intent(in) :: iDir
    !> The Penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> The projection datatype for the projection information
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The material information
    type(atl_material_type), intent(inout) :: material
    !> The state data in the nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    real(kind=rk), intent(in), optional :: nodal_GradData(:,:,:)
    !> The result in the nodal form
    real(kind=rk), intent(inout)     :: nodal_res(:,:)
    !> The length of the current element
    real(kind=rk), intent(in) :: ElemLength
    !> The scheme information of the min level (This is needed for the temp
    ! buffer array for evaluating the physical fluxes )
    type(atl_scheme_type), intent(inout) :: scheme_min
    !> Information about the current level
    type(atl_scheme_type), intent(inout) :: scheme_current
    ! --------------------------------------------------------------------------!
    ! Loop var for all the dof in an element
    integer :: iDof, nDofs
    ! Rotation indices for physical flux calculation
    integer :: rot(5)
    ! --------------------------------------------------------------------------!

    real(kind=rk)             :: state_0(5)

    real(kind=rk)             :: vsq,R_inv,R_inv_sq,isen,isencoef
    real(kind=rk)             :: A_21,A_22,A_23,A_24,A_25
    real(kind=rk)             :: A_31,A_32,A_33
    real(kind=rk)             :: A_41,A_42,A_44
    real(kind=rk)             :: A_51,A_52,A_53,A_54,A_55


    ! get the rotation for the physical flux calculation
    rot      = equation%varRotation(iDir)%varTransformIndices(1:5)
    nDofs    = poly_proj%body_3d%ndofs

    state_0  = state(1,rot)

    !defining the coefficients for matrix vector multiplication
    isencoef = equation%euler%isen_coef
    vsq      = (state_0(2)**2)+(state_0(3)**2)+(state_0(4)**2)
    R_inv    = 1._rk/state_0(1)
    R_inv_sq = R_inv**2
    isen     = (isencoef-1)

    A_21     = R_inv_sq*(0.5*isen*vsq - state_0(2)**2)
    A_22     = (3-isencoef)*state_0(2)*R_inv
    A_23     = -isen*state_0(3)*R_inv
    A_24     = -isen*state_0(4)*R_inv
    A_25     = isen

    A_31     = -state_0(2)*state_0(3)*R_inv_sq
    A_32     = R_inv*state_0(3)
    A_33     = R_inv*state_0(2)

    A_41     = -state_0(2)*state_0(4)*R_inv_sq
    A_42     = state_0(4)*R_inv
    A_44     = A_33

    A_51     = (isen*vsq-isencoef*state_0(5)*state_0(1))*state_0(2)*(R_inv**(3))
    A_52     =  R_inv_sq*(isencoef*state_0(5)*state_0(1)-(isen/2)*(vsq+2*(state_0(2)**2)))
    A_53     = -isen*state_0(2)*state_0(3)*R_inv_sq
    A_54     = -isen*state_0(2)*state_0(4)*R_inv_sq
    A_55     = isencoef*A_33

    ! non linear flux for the first mode

    res(1,rot(1))      = state_0(2)
    res(1,rot(2))      = R_inv*state_0(2)**2 + isen *(state_0(5)   &
      &                - 0.5_rk*R_inv*vsq )
    res(1,rot(3))      = R_inv*state_0(2)*state_0(3)
    res(1,rot(4))      = R_inv*state_0(2)*state_0(4)
    res(1,rot(5))      = R_inv*state_0(2)*state_0(5)*isencoef      &
      &                - 0.5_rk*isen*R_inv_sq*state_0(2)*vsq

    ! linear flux calculation for the higher modes
    dofLoop: do iDof   = 2, ndofs

      res(iDof,rot(1)) = state(iDof,rot(2))

      res(iDof,rot(2)) = A_21*state(iDof,1)      &
        &              + A_22*state(iDof,rot(2)) &
        &              + A_23*state(iDof,rot(3)) &
        &              + A_24*state(iDof,rot(4)) &
        &              + A_25*state(iDof,5)

      res(iDof,rot(3)) = A_31*state(iDof,1)      &
        &              + A_32*state(iDof,rot(2)) &
        &              + A_33*state(iDof,rot(3))

      res(iDof,rot(4)) = A_41*state(iDof,1)      &
        &              + A_42*state(iDof,rot(2)) &
        &              + A_44*state(iDof,rot(4))

      res(iDof,rot(5)) = A_51*state(iDof,1)      &
        &              + A_52*state(iDof,rot(2)) &
        &              + A_53*state(iDof,rot(3)) &
        &              + A_54*state(iDof,rot(4)) &
        &              + A_55*state(iDof,5)


    end do dofLoop


  end subroutine atl_modg_LoclinEuler_physFlux
  ! ****************************************************************************


end module atl_modg_LoclinEuler_kernel_module
