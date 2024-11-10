! Copyright (c) 2018-2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2018-2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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

module atl_eqn_sponge_module
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer, c_ptr
  use env_module,                  only: rk, labelLen

  use tem_time_module,             only: tem_time_type
  use tem_varSys_module,           only: tem_varSys_type
  use ply_poly_project_module,     only: ply_poly_project_type, &
    &                                    assignment(=)
  use atl_equation_source_module,  only: atl_equation_evaluate_source_nodal, &
    &                                    atl_compute_source_interface
  use atl_source_types_module,     only: atl_source_op_type
  use atl_cube_elem_module,        only: atl_cube_elem_type

  implicit none

  private

  public :: atl_eval_sponge
  public :: atl_eval_source_spongeLayer
  public :: atl_eval_source_spongeLayer_2d


contains


! ******************************************************************************
  subroutine atl_eval_sponge(rhs, source, state, constants)
    ! --------------------------------------------------------------------------
    !> The Right Hand side to be updated
    real(kind=rk), intent(inout) :: rhs(:,:)
    !> The source data to be used
    real(kind=rk), intent(in) :: source(:,:)
    !> The state in the modal form
    real(kind=rk), intent(in) :: state(:,:)
    !> the constants required for the evaluation of source
    real(kind=rk), intent(in) :: constants(:)
    ! --------------------------------------------------------------------------
    integer :: iComp, nComps
    integer :: iPoint, nQuadpoints
    ! --------------------------------------------------------------------------

    nQuadpoints = size(state,1)
    nComps =size(state,2)
    rhs = 0.0_rk
    ! Compute RHS using the nodal values of source and state
    do iComp = 1, nComps
      do iPoint = 1, nQuadpoints
        rhs(iPoint, iComp) = rhs(iPoint, iComp) &
          & - source(iPoint, 1) &
          & * ( state(iPoint, iComp) &
          &   - source(iPoint, iComp + 1) )
      end do
    end do

  end subroutine atl_eval_sponge
! ******************************************************************************


! ******************************************************************************
  subroutine atl_eval_source_spongeLayer( fun, varSys, time, mesh, poly_proj, &
    &                                     currentLevel, state, material,      &
    &                                     sourcedata                          )
    !---------------------------------------------------------------------------

    !> Description of method to update source
    class(atl_source_op_type), intent(in) :: fun

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> Current level mesh information
    type(atl_cube_elem_type), intent(in) :: mesh

    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> current level
    integer, intent(in) :: currentLevel

    !> The state in modal space.
    !! This is needed for several source terms that have to be applied to the
    !! current state
    real(kind=rk), intent(in) :: state(:,:,:)

    !> Material description for the complete domain. Used for evaluation of some
    !! source terms.
    real(kind=rk), intent(in) :: material(:)

    !> The source data to update. When all source terms are added to this
    !! buffer, it is applied to the state.
    real(kind=rk), intent(inout) :: sourcedata(:,:,:)
    ! --------------------------------------------------------------------------
    procedure(atl_compute_source_interface) , pointer:: evaluate_source
    ! --------------------------------------------------------------------------

    !> @todo PV: Create a unit test for this routine and compare it to the
    !! version before the new varSys

    ! Set the function pointer for the evaluation of spongeLayer
    evaluate_source => atl_eval_sponge

    ! Call the common function for updating the sourceData
    call atl_equation_evaluate_source_nodal(    &
      & fun          = fun,                     &
      & varSys       = varSys,                  &
      & currentLevel = currentLevel,            &
      & nDim         = 3,                       &
      & time         = time,                    &
      & eval_rhs     = evaluate_source,         &
      & state        = state,                   &
      & poly_proj    = poly_proj,               &
      & polyProjBody = poly_proj%body_3d,       &
      & sourceData   = sourceData               )


  end subroutine atl_eval_source_spongeLayer
! ******************************************************************************

! *******************************************************************************
  subroutine atl_eval_source_spongeLayer_2d( fun, varSys, time, mesh, poly_proj, &
    &                                    currentLevel, state, material,      &
    &                                    sourcedata                          )
    !---------------------------------------------------------------------------

    !> Description of method to update source
    class(atl_source_op_type), intent(in) :: fun

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> Current level mesh information
    type(atl_cube_elem_type), intent(in) :: mesh

    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> current level
    integer, intent(in) :: currentLevel

    !> The state in modal space.
    !! This is needed for several source terms that have to be applied to the
    !! current state
    real(kind=rk), intent(in) :: state(:,:,:)

    !> Material description for the complete domain. Used for evaluation of some
    !! source terms.
    real(kind=rk), intent(in) :: material(:)

    !> The source data to update. When all source terms are added to this
    !! buffer, it is applied to the state.
    real(kind=rk), intent(inout) :: sourcedata(:,:,:)
    ! ---------------------------------------------------------------------------
    procedure(atl_compute_source_interface), pointer:: evaluate_source
    ! ---------------------------------------------------------------------------

    !> @todo PV: Create a unit test for this routine and compare it to the
    !! version before the new varSys

    ! Set the function pointer for the evaluation of spongeLayer_2d
    evaluate_source => atl_eval_sponge

    ! Call the common function for updating the sourceData
    call atl_equation_evaluate_source_nodal(    &
      & fun          = fun,                     &
      & varSys       = varSys,                  &
      & currentLevel = currentLevel,            &
      & nDim         = 2,                       &
      & time         = time,                    &
      & eval_rhs     = evaluate_source,         &
      & state        = state,                   &
      & poly_proj    = poly_proj,               &
      & polyProjBody = poly_proj%body_2d,       &
      & sourceData   = sourceData               )

  end subroutine atl_eval_source_spongeLayer_2d
! *******************************************************************************
end module atl_eqn_sponge_module

