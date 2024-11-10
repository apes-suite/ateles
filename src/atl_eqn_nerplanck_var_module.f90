! Copyright (c) 2013-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

!> summary: module to configure information about the variables of the
!! Nernst-Planck equations
module atl_eqn_nerplanck_var_module
  use, intrinsic :: iso_c_binding,  only: c_loc

  use tem_varSys_module,            only: tem_varSys_type,              &
    &                                     tem_varSys_init,              &
    &                                     tem_varSys_append_stateVar,   &
    &                                     tem_varSys_proc_element,      &
    &                                     tem_varSys_proc_point,        &
    &                                     tem_varSys_proc_element,      &
    &                                     tem_varSys_proc_setparams,    &
    &                                     tem_varSys_proc_getparams,    &
    &                                     tem_varSys_proc_setupIndices, &
    &                                     tem_varSys_proc_getValOfIndex

  use atl_varSys_module,            only: atl_varSys_solverData_type,    &
    &                                     atl_varSys_getStateForElement, &
    &                                     atl_varSys_getStateForPoint,   &
    &                                     atl_get_new_varSys_data_ptr,   &
    &                                     atl_varSys_setupStateIndices,  &
    &                                     atl_varSys_getStateValofIndex

  implicit none

  private

  public :: atl_init_nerplanck_vars
  public :: atl_append_nerplanck_vars


contains


  !> summary: init the variables for nernst-planck equation
  subroutine atl_init_nerplanck_vars(varSys, hasPrimitiveVariables, methodData)
    ! --------------------------------------------------------------------------
    type(tem_varSys_type), intent(out) :: varSys
    logical, intent(out)               :: hasPrimitiveVariables
    type(atl_varSys_solverData_type)   :: methodData
    ! --------------------------------------------------------------------------
      ! initialize variable system
      call tem_varSys_init( me = varSys, systemName = 'nerplanck' )

      ! Append conservative Variables to variable system
      call atl_append_nerplanck_vars(varSys, methodData)

      hasPrimitiveVariables = .false.

  end subroutine atl_init_nerplanck_vars


  !> summary: append the variables for diffusion simulations
  subroutine atl_append_nerplanck_vars(varSys, methodData)
    ! --------------------------------------------------------------------------
    type(tem_varSys_type), intent(inout) :: varSys
    type(atl_varSys_solverData_type), target   :: methodData
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
    setup_indices => atl_varSys_setupStateIndices
    get_valOfIndex => atl_varSys_getStateValofIndex

    call tem_varSys_append_stateVar(                              &
      & me             = varSys,                                  &
      & varName        = 'concentration',                         &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex                           )

    call tem_varSys_append_stateVar(                              &
      & me             = varSys,                                  &
      & varName        = 'diffusiveFlux',                         &
      & nComponents    = 3,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex                           )

  end subroutine atl_append_nerplanck_vars


end module atl_eqn_nerplanck_var_module
