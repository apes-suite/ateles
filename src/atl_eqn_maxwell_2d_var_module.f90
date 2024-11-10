! Copyright (c) 2015-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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

!> summary: module to configure information about the variables of the maxwell
!! equations
module atl_eqn_maxwell_2d_var_module

  use, intrinsic :: iso_c_binding,      only: c_loc, c_ptr
  use env_module,                       only: rk

  use tem_time_module,                  only: tem_time_type
  use tem_varSys_module,                only: tem_varSys_type,              &
    &                                         tem_varSys_init,              &
    &                                         tem_varSys_append_stateVar,   &
    &                                         tem_varSys_proc_point,        &
    &                                         tem_varSys_proc_element,      &
    &                                         tem_varSys_proc_setparams,    &
    &                                         tem_varSys_proc_getparams,    &
    &                                         tem_varSys_proc_setupIndices, &
    &                                         tem_varSys_proc_getValOfIndex
  use tem_varMap_module,                only: tem_possible_variable_type, &
    &                                         init, append

  use ply_poly_project_module,          only: ply_poly_project_type, &
    &                                         assignment(=)

  use atl_equation_module,              only: atl_equations_type
  use atl_varSys_module,                only: atl_varSys_solverData_type,    &
    &                                         atl_varSys_getStateForElement, &
    &                                         atl_varSys_getStateForPoint,   &
    &                                         atl_get_new_varSys_data_ptr
  use atl_cube_elem_module,             only: atl_cube_elem_type
  use atl_source_types_module,          only: atl_eqn_sourceMap_type, &
    &                                         atl_source_op_type
  use atl_equation_source_module,               &
    & only: atl_equation_evaluate_source_nodal, &
    &       atl_equation_evaluate_source_modal, &
    &       atl_compute_source_interface

  implicit none

  private

  public :: atl_init_maxwell_2d_vars
  public :: atl_append_maxwell_2d_vars
  public :: atl_init_maxwell_2d_sourceTerms
  public :: atl_eval_source_currentDensity_2d


contains


! ******************************************************************************
  !> summary: init the variables for maxwell equation (2D, TE-mode formulation).
  subroutine atl_init_maxwell_2d_vars( equation, methodData )
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type) :: methodData
    ! --------------------------------------------------------------------------
    ! initialize variable system
    call tem_varSys_init( me = equation%varSys, systemName = 'maxwell_2d' )

    allocate(equation%stateVar(4))

    ! Append conservative Variables to variable system
    call atl_append_maxwell_2d_vars(  equation, methodData  )

    equation%hasPrimitiveVariables = .false.

    ! Set values fro allocating temp flux arrays
    equation%temp%overSamp = 1
    equation%temp%modal = 0
    equation%temp%nodal = 2
    equation%temp%nScal = equation%varsys%nScalars
  end subroutine atl_init_maxwell_2d_vars
! ******************************************************************************


! ******************************************************************************
  !> summary: append the variables for electrodynamic simulations
  subroutine atl_append_maxwell_2d_vars(equation, methodData)
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), target :: methodData
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer &
      & :: get_valOfIndex => NULL()
    ! --------------------------------------------------------------------------

    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    set_params => null()
    get_params => null()
    setup_indices => null()
    get_valOfIndex => null()

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'displacement_field',                    &
      & nComponents    = 2,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(1)                     )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'magnetic_field',                        &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(2)                     )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'pml_P',                                 &
      & nComponents    = 2,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(3)                     )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'pml_Q',                                 &
      & nComponents    = 2,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(4)                     )

  end subroutine atl_append_maxwell_2d_vars
! ******************************************************************************


! ******************************************************************************
  subroutine eval_currentDensity_2d(rhs, source, state, constants)
    ! --------------------------------------------------------------------------
    !> The Right Hand side to be updated
    real(kind=rk), intent(inout) :: rhs(:,:)
    !> The source data to be used
    real(kind=rk), intent(in) :: source(:,:)
    !> The state in the modal form
    real(kind=rk), intent(in) :: state(:,:)
    !> the constants required for the evaluation of source
    real(kind = rk ), intent(in) :: constants(:)
    ! --------------------------------------------------------------------------
    rhs = 0.0_rk

    ! Compute RHS using the nodal values of source and state
    ! current density is relevant in Ampere's law:
    rhs(:,1:2) = - source(:,1:2)

  end subroutine eval_currentDensity_2d
! ******************************************************************************


! ******************************************************************************
  !> summary: evaluate "currentDensity" source
  subroutine atl_eval_source_currentDensity_2d( fun, varSys, time, mesh,    &
    &                                           poly_proj, currentLevel,    &
    &                                           state, material, sourcedata )
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
    ! Set the function pointer for the evaluation of spongeLayer_2d
    evaluate_source => eval_currentDensity_2d

    ! Call the common function for updating the sourceData
    call atl_equation_evaluate_source_modal(    &
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


  end subroutine atl_eval_source_currentDensity_2d
! ******************************************************************************


! ******************************************************************************
  subroutine eval_pml_2d(rhs, source, state, constants)
    ! --------------------------------------------------------------------------
    !> The Right Hand side to be updated
    real(kind=rk), intent(inout) :: rhs(:,:)
    !> The source data to be used
    real(kind=rk), intent(in) :: source(:,:)
    !> The state in the modal form
    real(kind=rk), intent(in) :: state(:,:)
    !> the constants required for the evaluation of source
    real(kind = rk ), intent(in) :: constants(:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: permittivity, permeability
    ! --------------------------------------------------------------------------
    rhs = 0.0_rk
    permeability = constants(1)
    permittivity = constants(2)

    ! Compute RHS using the nodal values of source and state

    ! ... pml for d_t D_x = ...
    !              - 2 kappa_2 (1/sqrt(eps mu)) D_1 - kappa_2 sqrt(eps/mu) P_2
    RHS(:,1) = RHS(:,1)                                              &
      & - 2.0_rk * source(:,2) * state(:,1)                          &
      &   / sqrt(permittivity * permeability)                        &
      & - source(:,2) * sqrt(permittivity / permeability) * state(:,5)

    ! ... pml for d_t D_y = ...
    !               - 2 kappa_1 (1/sqrt(eps mu)) D_2 - kappa_1 sqrt(eps/mu) P_1
    RHS(:,2) = RHS(:,2)                                              &
      & - 2.0_rk * source(:,1) * state(:,2)                          &
      &   / sqrt(permittivity * permeability)                        &
      & - source(:,1) * sqrt(permittivity / permeability) * state(:,4)

    ! ... pml for d_t B_z = ... + (d_x kappa_1) Q_1 - (d_y kappa_2) Q_2
    RHS(:,3) = RHS(:,3) + source(:,3)*state(:,6) - source(:,4)*state(:,7)

    ! ... pml for d_t P_1 = ... + kappa_1 * (1/sqrt(eps mu)) * (D_2/eps)
    RHS(:,4) = RHS(:,4)                                         &
      & + source(:,1) * (1 / sqrt(permittivity * permeability)) &
      &   * (state(:,2)/permittivity)

    ! ... pml for d_t P_2 = ... + kappa_2 * (1/sqrt(eps mu)) * (D_1/eps)
    RHS(:,5) = RHS(:,5)                                         &
      & + source(:,2) * (1 / sqrt(permittivity * permeability)) &
      &   * (state(:,1) / permittivity)

    ! ... pml for d_t Q_1 = ...
    !         - kappa_1 * (1/sqrt(eps mu)) * Q_1 - (1/sqrt(eps mu)) * (D_2/eps)
    RHS(:,6) = RHS(:,6)                                                       &
      & - source(:,1) * (1 / sqrt(permittivity * permeability)) * state(:,6)  &
      & - (1 / sqrt(permittivity * permeability)) * (state(:,2) / permittivity)

    ! ... pml for d_t Q_2 = ...
    !         - kappa_2 * (1/sqrt(eps mu)) * Q_2 - (1/sqrt(eps mu)) * (D_1/eps)
    RHS(:,7) = RHS(:,7)                                                       &
      & - source(:,2) * (1 / sqrt(permittivity * permeability)) * state(:,7)  &
      & - (1 / sqrt(permittivity * permeability)) * (state(:,1) / permittivity)

  end subroutine eval_pml_2d
! ******************************************************************************


! ******************************************************************************
  !> summary: evaluate "pml" source (uniaxial PML)
  subroutine eval_source_pml_2d( fun, varSys, time, mesh, poly_proj,       &
    &                            currentLevel, state, material, sourcedata )
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
    ! Set the function pointer for the evaluation of spongeLayer_2d
    evaluate_source => eval_pml_2d

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
      & consts       = material,                &
      & sourceData   = sourceData               )

  end subroutine eval_source_pml_2d
! ******************************************************************************


! ******************************************************************************
  !> summary: init source terms for electrodynamic simulations.
  subroutine atl_init_maxwell_2d_sourceTerms(possVars, eval_source)
    ! --------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(out) :: possVars
    type(atl_eqn_sourceMap_type), allocatable, intent(out) :: eval_source(:)
    ! --------------------------------------------------------------------------
    integer :: pos
    ! --------------------------------------------------------------------------

    allocate(eval_source(2))
    call init(me = possVars, length = 2)

    ! ... current density
    call append( me          = possVars,          &
      &          varName     = 'current_density', &
      &          nComponents = 2,                 &
      &          pos         = pos                )
    eval_source(pos)%do => atl_eval_source_currentDensity_2d

    ! ... pml layer (returns 2-vector with normal-projection onto
    ! the starting plane of the PML and the derivative of the first entry
    ! with regards to x and the derivative of the second entry with regards
    ! to y.
    call append( me          = possVars, &
      &          varName     = 'pml',    &
      &          nComponents = 4,        &
      &          pos         = pos       )
    eval_source(pos)%do => eval_source_pml_2d

  end subroutine atl_init_maxwell_2d_sourceTerms
! ******************************************************************************

end module atl_eqn_maxwell_2d_var_module
