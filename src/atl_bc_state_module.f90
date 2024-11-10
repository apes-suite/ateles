! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> A module to extend tem_bc_state with Ateles specific information.
module atl_bc_state_module
  use aotus_module,                   only: flu_State

  use tem_bc_module,                  only: tem_bc_state_type, &
    &                                       tem_load_bc_state
  use tem_stringKeyValuePair_module,  only: grw_stringKeyValuePairArray_type
  use tem_varSys_module,              only: tem_varSys_type, &
    &                                       tem_varSys_solverData_evalElem_type

  implicit none

  private

  public :: atl_load_bc_state, atl_bc_state_set_frompoint

  type(tem_varSys_solverData_evalElem_type) :: solverData_evalElem


contains


  ! ************************************************************************ !
  !> Define the method to set the solverData_evalElem routine for stfuns.
  !!
  !! This needs to be done before any atl_load_bc_state call is made.
  subroutine atl_bc_state_set_fromPoint( stfun_solverelem )
    ! -------------------------------------------------------------------- !
    !> Datatype describing the setter callback function and the varsys data
    !! we need to do the projection.
    type(tem_varSys_solverData_evalElem_type), intent(in) :: stfun_solverelem
    ! -------------------------------------------------------------------- !

    ! Just set the module variable here to allow its proper usage in the
    ! atl_load_bc_state calls.
    solverData_evalElem = stfun_solverelem

  end subroutine atl_bc_state_set_fromPoint
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load the boundary condition for state variables.
  !!
  !! This is only a thin wrapper around tem_load_bc_state, to always provide
  !! the solverData_evalElem routine.
  subroutine atl_load_bc_state( bc, state_name, nComp, style, conf, &
    &                           bc_handle, varDict, varSys, ErrCode )
    ! -------------------------------------------------------------------- !
    !> The boundary to fill
    type(tem_bc_state_type), intent(inout) :: bc
    !> The state variable to set with this boundary condition
    character(len=*), intent(in) :: state_name
    !> Number of Components in this boundary variable.
    integer, intent(in), optional :: nComp
    !> Style of this boundary condition
    !! dirichlet = set value itself
    !! neumann = set derivative of value
    character(len=*), optional, intent(in) :: style
    type(flu_State),intent(in) :: conf !< Lua state
    !> Handle to the table describing the boundary
    integer, intent(in) :: bc_handle
    !> The dictionary that contains the mapping between expected variables
    !! and the actual variables defined by the user.
    type(grw_stringKeyValuePairArray_type), intent(inout) :: varDict
    type(tem_varSys_type), intent(inout) :: varSys
    !> Error code
    integer, optional, intent(out) :: ErrCode
    ! -------------------------------------------------------------------- !

    ! Insert solverData_evalElem from the module variable here
    call tem_load_bc_state(                           &
      &    bc                  = bc,                  &
      &    state_name          = state_name,          &
      &    nComp               = nComp,               &
      &    style               = style,               &
      &    conf                = conf,                &
      &    bc_handle           = bc_handle,           &
      &    varDict             = varDict,             &
      &    varSys              = varSys,              &
      &    solverData_evalElem = solverData_evalElem, &
      &    ErrCode             = ErrCode              )

  end subroutine atl_load_bc_state
  ! ************************************************************************ !


end module atl_bc_state_module
