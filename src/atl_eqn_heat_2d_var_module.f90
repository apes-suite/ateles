! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> Module to configure the variables of the Heat equations.
module atl_eqn_heat_2d_var_module
  use, intrinsic :: iso_c_binding,  only: c_loc

  use tem_varSys_module,            only: tem_varSys_init,              &
    &                                     tem_varSys_append_stateVar,   &
    &                                     tem_varSys_proc_element,      &
    &                                     tem_varSys_proc_point,        &
    &                                     tem_varSys_proc_setparams,    &
    &                                     tem_varSys_proc_getparams,    &
    &                                     tem_varSys_proc_setupIndices, &
    &                                     tem_varSys_proc_getValOfIndex

  use atl_equation_module,          only: atl_equations_type
  use atl_varSys_module,            only: atl_varSys_solverData_type,    &
    &                                     atl_varSys_getStateForElement, &
    &                                     atl_varSys_getStateForPoint,   &
    &                                     atl_get_new_varSys_data_ptr,   &
    &                                     atl_varSys_setupStateIndices,  &
    &                                     atl_varSys_getStateValofIndex


  implicit none

  private

  public :: atl_init_heat_2d_vars


contains


  !> Init the variable system for simulations of Heat equation.
  !!
  !! The variable system describes, which variables are to be used and how
  !! they are organized in the memory.
  !! The first few variables up to the sys_mark are those, describing the
  !! state, and thus describe the output for regular restart files.
  !! Here these are the conservative variables density, momentum and energy.
  !! After the mark, there are additional values described that can be deduced
  !! from the state variables.
  subroutine atl_init_heat_2d_vars(  equation, methodData  )
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type) :: methodData
    ! --------------------------------------------------------------------------

    ! Initialize variable system
    call tem_varSys_init( me = equation%varSys, systemName = 'heat_2d' )

    ! Append the conservative variables
    call append_heat_2d_consVars(equation, methodData)

    ! Append primitive variables
    equation%hasPrimitiveVariables = .false.

    ! Set values fro allocating temp flux arrays
    equation%temp%overSamp = 0
    equation%temp%modal = 2
    equation%temp%nodal = 0
    equation%temp%nScal = 1
  end subroutine atl_init_heat_2d_vars


  !> Append conservative variables for Heat equations.
  subroutine append_heat_2d_consVars(  equation, methodData  )
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
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    ! --------------------------------------------------------------------------

    allocate(equation%stateVar(1))

    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    setup_Indices => atl_varSys_setupStateIndices
    get_valOfindex => atl_varSys_getStateValOfIndex

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'temperature',                           &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(1)                     )

  end subroutine append_heat_2d_consVars

end module atl_eqn_heat_2d_var_module
