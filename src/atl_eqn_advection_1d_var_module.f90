! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

!> Module to configure the variables of the Euler equations.
module atl_eqn_advection_1d_var_module
  use, intrinsic :: iso_c_binding,    only: c_loc

  use tem_varSys_module,              only: tem_varSys_type,              &
    &                                       tem_varSys_init,              &
    &                                       tem_varSys_append_stateVar,   &
    &                                       tem_varSys_proc_element,      &
    &                                       tem_varSys_proc_point,        &
    &                                       tem_varSys_proc_setparams,    &
    &                                       tem_varSys_proc_getparams,    &
    &                                       tem_varSys_proc_setupIndices, &
    &                                       tem_varSys_proc_getValOfIndex

  use atl_varSys_module,              only: atl_varSys_solverData_type, &
    &                                       atl_get_new_varSys_data_ptr
  use atl_varSys_module,              only: atl_varSys_getStateForElement, &
    &                                       atl_varSys_getStateForPoint,   &
    &                                       atl_varSys_setupStateIndices,  &
    &                                       atl_varSys_getStateValofIndex


  implicit none

  private

  public :: atl_init_advection_1d_vars


contains


  !> Init the variable system for simulations of advection equation.
  !!
  !! The variable system describes, which variables are to be used and how
  !! they are organized in the memory.
  !! The first few variables up to the sys_mark are those, describing the
  !! state, and thus describe the output for regular restart files.
  !! Here these are the conservative variables density, momentum and energy.
  !! After the mark, there are additional values described that can be deduced
  !! from the state variables.
  subroutine atl_init_advection_1d_vars( varSys, hasPrimitiveVariables, &
      &                                  methodData                     )
    ! --------------------------------------------------------------------------
    !> The resulting variable system used in the Euler equations
    type(tem_varSys_type), intent(out) :: varSys

    !> A flag indicating, that this system has primitive variables.
    logical, intent(out) :: hasPrimitiveVariables

    !> Data required for the varsys to perform different operations of it's
    !! variables
    type(atl_varSys_solverData_type) :: methodData
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Initialize variable system
    call tem_varSys_init( me = varSys, systemName = 'advection_1d')

    ! Append the conservative variables
    call append_advection_1d_consVars(varSys, methodData)

    ! Append primitive variables
    hasPrimitiveVariables = .false.

  end subroutine atl_init_advection_1d_vars


  !> Append conservative variables for Euler equations.
  !!
  !! These are density, momentum and energy here.
  subroutine append_advection_1d_consVars(varSys, methodData)
    ! --------------------------------------------------------------------------
    !> The variable system to append the variables to.
    type(tem_varSys_type), intent(inout)  :: varSys
    !> Data required for the varsys to perform different operations of it's
    !! variables
    type(atl_varSys_solverData_type), target    :: methodData
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    ! --------------------------------------------------------------------------

    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    setup_indices => atl_varSys_setupStateIndices
    get_valOfIndex => atl_varSys_getStateValOfIndex

    call tem_varSys_append_stateVar(                              &
      & me             = varSys,                                  &
      & varName        = 'density',                               &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex                           )

    get_element => null()
    get_point => null()
    set_params => null()
    get_params => null()
    setup_indices => null()
    get_valOfIndex => null()

  end subroutine append_advection_1d_consVars

end module atl_eqn_advection_1d_var_module
