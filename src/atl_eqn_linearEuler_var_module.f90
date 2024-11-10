! Copyright (c) 2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> summary: module to configure information about the variables of the acoustic
!! equations
module atl_eqn_LinearEuler_var_module

  use, intrinsic :: iso_c_binding,       only: c_loc
  use env_module,                        only: labelLen

  use tem_aux_module,                    only: tem_abort
  use tem_varSys_module,                 only: tem_varSys_init,              &
    &                                          tem_varSys_type,              &
    &                                          tem_varSys_append_stateVar,   &
    &                                          tem_varSys_append_derVar,     &
    &                                          tem_varSys_proc_point,        &
    &                                          tem_varSys_proc_element,      &
    &                                          tem_varSys_proc_setparams,    &
    &                                          tem_varSys_proc_getparams,    &
    &                                          tem_varSys_proc_setupIndices, &
    &                                          tem_varSys_proc_getValOfIndex
  use tem_logging_module,                only: logUnit
  use atl_equation_module,               only: atl_equations_type
  use atl_eqn_linearEuler_derive_module, only: atl_linEuler_completState_getElement,&
    &                                          atl_linEuler_completState_getPoint
  use atl_varSys_module,                 only: atl_varSys_solverData_type,    &
    &                                          atl_varSys_getStateForElement, &
    &                                          atl_varSys_getStateForPoint,   &
    &                                          atl_get_new_varSys_data_ptr,   &
    &                                          atl_varSys_setupStateIndices,  &
    &                                          atl_varSys_getStateValofIndex
  use tem_varMap_module,                 only: tem_possible_variable_type, &
    &                                          init, append
  use atl_source_types_module,           only: atl_eqn_sourceMap_type
  use atl_eqn_sponge_module,             only: atl_eval_source_spongeLayer

  implicit none

  private

  public :: atl_init_LinearEuler_vars
  public :: atl_append_LinearEuler_vars
  public :: atl_init_lineuler_sourceTerms

contains


! *******************************************************************************
  !> summary: init the variables for LinearEuler equation
  subroutine atl_init_LinearEuler_vars( equation, methodData )
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), intent(in) :: methodData
    ! --------------------------------------------------------------------------
    ! Initialize variable system
    call tem_varSys_init( me = equation%varSys, systemName = 'LinearEuler' )

    ! Append conservative Variables to variable system
    call atl_append_LinearEuler_vars(equation, methodData)

    ! Append derived Variables to variable system
    call atl_append_LinearEuler_derivedVars( equation%varSys, methodData )
  end subroutine atl_init_LinearEuler_vars
! *******************************************************************************


! *******************************************************************************
  !> summary: append the variables for LinearEuler simulations
  !! the acoustic equation only has 'primitive' variables or different speaking,
  !! the equation describe the pertubation in primitive variables
  subroutine atl_append_linearEuler_vars(equation, methodData)
    ! ---------------------------------------------------------------------------
    !> The equation type.
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), target :: methodData
    ! ---------------------------------------------------------------------------
    procedure(tem_varSys_proc_point), pointer :: get_point => null()
    procedure(tem_varSys_proc_element), pointer :: get_element => null()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    ! ---------------------------------------------------------------------------

    allocate(equation%stateVar(3))

    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    setup_indices => atl_varSys_setupStateIndices
    get_valOfIndex => atl_varSys_getStateValofIndex

    ! Append variables to varSys
    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'density',                               &
      & nComponents    = 1,                                       &
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
      & varName        = 'velocity',                              &
      & nComponents    = 3,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                         &
      & pos            = equation%stateVar(2)                    )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'pressure',                              &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(3)                     )

  end subroutine atl_append_LinearEuler_vars
! *******************************************************************************


! ******************************************************************************
  !> summary: append / set data of derived quantities
  subroutine atl_append_lineareuler_derivedVars( varSys, methodData )
    ! --------------------------------------------------------------------------
    !> The Euler variable system to modify. It has to contain the conservative
    !! and primitive variables already.
    type(tem_varSys_type), intent(inout) :: varSys
    !> the pointer to the data required for the varsys to fulfill all operations
    !! and derivations on the variables
    type(atl_varSys_solverData_type), target :: methodData
    ! --------------------------------------------------------------------------
    integer :: nDerivedVars, iVar, nComponents
    character(len=20), allocatable :: derVarName(:)
    character(len=labelLen), allocatable :: invar_name(:)
    character(len=labelLen) :: varname
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer &
      & :: get_valOfIndex => NULL()
    ! --------------------------------------------------------------------------
    nDerivedVars = 1
    allocate(derVarName(nDerivedVars))
    derVarName    = [ 'completeState' ]

    do iVar = 1, nDerivedVars
      varname = trim(adjustl(derVarName(iVar)))

      select case(varname)

      ! Define the complet State meaning adding perturbation + background state
      case ('completeState')
        get_point => atl_linEuler_completState_getPoint
        get_element => atl_linEuler_completState_getElement
        set_params => null()
        get_params => null()
        setup_indices => null()
        get_valOfIndex => null()
        nComponents = 5
        allocate(invar_name(3))
        invar_name(1) = 'density'
        invar_name(2) = 'velocity'
        invar_name(3) = 'pressure'

      case default
        write(logUnit(1),*) 'WARNING: Unknown variable: '//trim(varname)
        cycle !go to next variable

      end select

      ! append variable to varSys
      call tem_varSys_append_derVar(                                &
        & me             = varSys,                                  &
        & varName        = varname,                                 &
        & operType       = 'st_fun',                                &
        & nComponents    = nComponents,                             &
        & input_varname  = invar_name,                              &
        & method_data    = atl_get_new_varSys_data_ptr(methodData), &
        & get_point      = get_point,                               &
        & get_element    = get_element,                             &
        & set_params     = set_params,                              &
        & get_params     = get_params,                              &
        & setup_indices  = setup_indices,                           &
        & get_valOfIndex = get_valOfIndex,                          &
        & wasAdded       = wasAdded                                 )

      if (wasAdded) then
        write(logUnit(10),*) 'Appended variable: ' // trim(varname)
      else
        call tem_abort( 'Error: variable ' // trim(varname) &
          & // ' is not added to variable system'           )
      end if

      deallocate(invar_name)

    end do

    deallocate(derVarName)

  end subroutine atl_append_lineareuler_derivedVars
! ******************************************************************************

! ******************************************************************************
  !> Init source terms for flow simulations.
  !> This routine initializes possible source variables and returns the filled
  !! up list of the poss_srcVars
  subroutine atl_init_lineuler_sourceTerms(possVars, eval_source)
    ! --------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(inout)  :: possVars
    type(atl_eqn_sourceMap_type), allocatable, intent(out) :: eval_source(:)
    ! --------------------------------------------------------------------------
    integer :: pos
    ! --------------------------------------------------------------------------

    allocate(eval_source(3))
    call init( me = possVars, length = 3 )

    ! Add the spongeLayer Source term
    call append( me          = possVars,      &
      &          varName     = 'spongelayer', &
      &          nComponents = 6,             &
      &          pos         = pos            )
    eval_source(pos)%do => atl_eval_source_spongeLayer

  end subroutine atl_init_lineuler_sourceTerms
! ******************************************************************************


end module atl_eqn_LinearEuler_var_module
