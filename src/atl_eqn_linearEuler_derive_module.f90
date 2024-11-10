! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> Routines to derive quantities from the state in the linearized Euler equation
! system.
module atl_eqn_lineareuler_derive_module
  use, intrinsic :: iso_c_binding,  only: c_f_pointer
  use env_module,                   only: rk

  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use tem_varSys_module,            only: tem_varSys_type,   &
    &                                     tem_varSys_op_type

  use atl_eqn_linearEuler_module,   only: atl_LinearEuler_type
  use atl_varSys_module,            only: atl_varSys_data_type

  implicit none

  private

  public :: atl_linEuler_completState_getElement
  public :: atl_linEuler_completState_getPoint

contains


! ******************************************************************************
  subroutine atl_linEuler_completState_getPoint( fun, varsys, point, time,tree, &
      &                                          nPnts, res                     )
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    type(atl_linearEuler_type), pointer :: linEulerPtr
    real(kind=rk) :: density(nPnts), velocity(3*nPnts), pressure(nPnts)
    ! ---------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )
    linEulerPtr => fPtr%solverData%equationPtr%LinearEuler

    ! get the state values
    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = density                                  )
    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = velocity                                 )
    call varSys%method%val(fun%input_varPos(3))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = pressure                                 )

    ! add state + background
    density  = density + linEulerPtr%density_0
    velocity(1:3*nPnts:3) = velocity(1:3*nPnts:3) + linEulerPtr%velocity_0(1)
    velocity(2:3*nPnts:3) = velocity(2:3*nPnts:3) + linEulerPtr%velocity_0(2)
    velocity(3:3*nPnts:3) = velocity(3:3*nPnts:3) + linEulerPtr%velocity_0(3)
    pressure = pressure + linEulerPtr%pressure_0

    ! transform it into res
    res(1:fun%nComponents*nPnts:fun%nComponents) = density
    res(2:fun%nComponents*nPnts:fun%nComponents) = velocity(1:3*nPnts:3)
    res(3:fun%nComponents*nPnts:fun%nComponents) = velocity(2:3*nPnts:3)
    res(4:fun%nComponents*nPnts:fun%nComponents) = velocity(3:3*nPnts:3)
    res(5:fun%nComponents*nPnts:fun%nComponents) = pressure

  end subroutine atl_linEuler_completState_getPoint

  subroutine atl_linEuler_completState_getElement(fun, varsys, elempos, time, &
      &                                           tree, nElems, nDofs, res    )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    type(atl_linearEuler_type), pointer :: linEulerPtr
    real(kind=rk) :: density(nElems), velocity(3*nElems), pressure(nElems)
    ! ---------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )
    linEulerPtr => fPtr%solverData%equationPtr%LinearEuler

    ! get the state values
    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = density                                    )
    call varSys%method%val(fun%input_varPos(2))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = velocity                                   )
    call varSys%method%val(fun%input_varPos(3))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = pressure                                   )

    ! add state + background
    density  = density + linEulerPtr%density_0
    velocity(1::3) = velocity(1::3) + linEulerPtr%velocity_0(1)
    velocity(2::3) = velocity(2::3) + linEulerPtr%velocity_0(2)
    velocity(3::3) = velocity(3::3) + linEulerPtr%velocity_0(3)
    pressure = pressure + linEulerPtr%pressure_0

    ! transform it into res
    res(1::fun%nComponents) = density
    res(2::fun%nComponents) = velocity(1::3)
    res(3::fun%nComponents) = velocity(2::3)
    res(4::fun%nComponents) = velocity(3::3)
    res(5::fun%nComponents) = pressure

    ! transform it into res
  end subroutine atl_linEuler_completState_getElement

end module atl_eqn_linearEuler_derive_module
