! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
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

!> Routines to derive quantities from the state in the Euler equation system.
module atl_eqn_euler_1d_derive_module
  
  use, intrinsic :: iso_c_binding,    only: c_f_pointer
  use env_module,                     only: rk

  use atl_equation_module,            only: atl_equations_type
  use tem_time_module,                only: tem_time_type
  use tem_varSys_module,              only: tem_varSys_type,   &
    &                                       tem_varSys_op_type
  use treelmesh_module,               only: treelmesh_type
  use tem_logging_module,             only: logUnit

  use atl_equation_module,            only: atl_equations_type
  use atl_varSys_module,              only: atl_varSys_data_type
  use atl_derive_module,              only: atl_derive_inputVar_type,        &
    &                                       atl_derive_fromModalData,        &
    &                                       atl_generic_fromModal_getElement
  implicit none

  private

  public :: atl_eqn_euler_1d_cons2prim
  public :: atl_eqn_euler_1d_prim2cons
  public :: atl_eqn_euler_1d_cons2prim_elems
  public :: atl_eqn_euler_1d_prim2cons_elems
  public :: atl_pressure_1d_getPoint
  public :: atl_pressure_1d_getIndex
  public :: atl_derivePressure_1d
  public :: atl_pressure_1d_getElement

contains

! ******************************************************************************!
  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_1d_prim2cons(equation, instate, outstate, material)
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type),   intent(in)    :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional, intent(in)    :: material(:,:)
    ! --------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,1) * instate(:,2)
      outstate(:,3) = instate(:,3) / (equation%euler%isen_coef - 1.0_rk)         &
        & + 0.5_rk * instate(:,1) * ( instate(:,2)**2 )
    else
      instate(:,3) = instate(:,3) / (equation%euler%isen_coef - 1.0_rk)         &
        & + 0.5_rk * instate(:,1) * ( instate(:,2)**2 )
      instate(:,2) = instate(:,1)*instate(:,2)
    end if
  end subroutine atl_eqn_euler_1d_prim2cons
! ******************************************************************************!


! ******************************************************************************!
  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_1d_cons2prim(equation, instate, outstate, material)
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type),   intent(in)    :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional, intent(in)    :: material(:,:)

    ! --------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,2)/instate(:,1)
      outstate(:,3) = (equation%euler%isen_coef - 1.0_rk) * ( instate(:,3) &
        & - 0.5_rk * instate(:,1) * (outstate(:,2)**2 ) )
    else
      instate(:,2) = instate(:,2)/instate(:,1)
      instate(:,3) = (equation%euler%isen_coef - 1.0_rk) * ( instate(:,3) &
        & - 0.5_rk * instate(:,1) * ( instate(:,2)**2 ) )
    end if
  end subroutine atl_eqn_euler_1d_cons2prim
! ******************************************************************************!

! ******************************************************************************!
  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_1d_prim2cons_elems( equation, instate, outstate, &
    &                                          nElems                       )
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type),   intent(in)    :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:,:)

    !> Number of elements to act on (first index in the state arrays).
    integer,                 intent(in)    :: nElems
    ! --------------------------------------------------------------------------!

    if (present(outstate)) then
      outstate(1:nElems,:,1) = instate(1:nElems,:,1)
      outstate(1:nElems,:,2) = instate(1:nElems,:,1)*instate(1:nElems,:,2)
      outstate(1:nElems,:,3) =                                            &
        & instate(1:nElems,:,3) / (equation%euler%isen_coef-1.0_rk)       &
        &   + 0.5_rk * instate(1:nElems,:,1) * ( instate(1:nElems,:,2)**2 )
    else
      instate(1:nElems,:,3) =                                             &
        & instate(1:nElems,:,3) / (equation%euler%isen_coef-1.0_rk)       &
        &   + 0.5_rk * instate(1:nElems,:,1) * ( instate(1:nElems,:,2)**2 )
      instate(1:nElems,:,2) = instate(1:nElems,:,1)*instate(1:nElems,:,2)
    end if
  end subroutine atl_eqn_euler_1d_prim2cons_elems
! ******************************************************************************!

! ******************************************************************************!
  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_1d_cons2prim_elems( equation, instate, outstate, &
    &                                          nElems                       )
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type),   intent(in)    :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:,:)

    !> Number of elements to act on (first index in the state arrays).
    integer,                 intent(in)    :: nElems
    ! --------------------------------------------------------------------------!

    if (present(outstate)) then
      outstate(1:nElems,:,1) = instate(1:nElems,:,1)
      outstate(1:nElems,:,2) = instate(1:nElems,:,2)/instate(1:nElems,:,1)
      outstate(1:nElems,:,3) = (equation%euler%isen_coef - 1.0_rk) &
        & * ( instate(1:nElems,:,3) - 0.5_rk * instate(1:nElems,:,1) &
        &   * (outstate(1:nElems,:,2)**2) )
    else
      instate(1:nElems,:,2) = instate(1:nElems,:,2)/instate(1:nElems,:,1)
      instate(1:nElems,:,3) = (equation%euler%isen_coef - 1.0_rk) &
        & * ( instate(1:nElems,:,3) - 0.5_rk * instate(1:nElems,:,1) &
        &   * (instate(1:nElems,:,2)**2) )
    end if
  end subroutine atl_eqn_euler_1d_cons2prim_elems
! ******************************************************************************!

! ******************************************************************************!
  subroutine atl_pressure_1d_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! ---------------------------------------------------------------------------
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
    real(kind=rk) :: density(nPnts), momentum(nPnts), energy(nPnts)
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )

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
      & res     = momentum                                 )
    call varSys%method%val(fun%input_varPos(3))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = energy                                   )

    res = (fPtr%solverData%equationPtr%Euler%isen_coef - 1.0_rk) &
      &   * ( energy                                             &
      &        - ( 0.5_rk / density * momentum**2 )              )

  end subroutine atl_pressure_1d_getPoint

! ******************************************************************************!
  subroutine atl_pressure_1d_getIndex( fun, varSys, time, iLevel, &
    &                               idx, idxLen, nVals,  res   )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables,
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: n
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: density(nVals), momentum(nVals), energy(nVals)
    ! ---------------------------------------------------------------------------
    write(logUnit(4),*) 'Get the values of indices for derived variable ',  &
      &                  trim(varSys%varname%val(fun%myPos))

    ! Initialize the momentum to zero
    density = 0.0_rk
    momentum = 0.0_rk
    energy = 0.0_rk

    call C_F_POINTER( fun%method_Data, fPtr )

    call varSys%method%val(fun%input_varPos(1))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = fPtr%opData%input_pntIndex(1)                 &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = density                                       )
    call varSys%method%val(fun%input_varPos(2))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = fPtr%opData%input_pntIndex(2)                 &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = momentum                                      )
    call varSys%method%val(fun%input_varPos(3))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = fPtr%opData%input_pntIndex(3)                 &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = energy                                        )

    res = (fPtr%solverData%equationPtr%Euler%isen_coef - 1.0_rk)  &
      &   * ( energy                                              &
      &      - ( 0.5_rk / density * momentum**2)                  )

  end subroutine atl_pressure_1d_getIndex

! ******************************************************************************!
  subroutine atl_derivePressure_1d(fun, varsys, tree, iElem, elemPos,    &
    &                                        nodalInput, nodalRes     )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> The Current element index
    integer, intent(in) :: iElem

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> The input data. nodalInput contains one entry for each input variable.
    !! This entry itself contains the nodal data for the dofs and components of
    !! the input variable. These nodal data has to be gained by oversampling
    !! and projecting the modal state into nodal space.
    type(atl_derive_inputVar_type) :: nodalInput(:)
    !> The result in nodal space
    real(kind=rk), allocatable :: nodalRes(:,:)
    ! ---------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer, parameter :: density = 1, momentum = 2, energy = 3
    ! ---------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    nodalRes(:,1) = (fPtr%solverData%equationPtr%Euler%isen_coef - 1.0_rk) &
      &   * (nodalInput(energy)%data(:,1)                                  &
      &   - 0.5_rk                                                         &
      &   / nodalInput(density)%data(:,1)                                  &
      &   * nodalInput(momentum)%data(:,1)**2                              )

  end subroutine atl_derivePressure_1d

! ******************************************************************************!
  subroutine atl_pressure_1d_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                    nDofs, res                               )
    ! -----------------------------------------------------------------------------
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
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! ---------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    fnCalcPtr => atl_derivePressure_1d

    call atl_generic_fromModal_getELement( &
      & fun        = fun,                  &
      & varsys     = varsys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr%solverData,      &
      & res        = res                   )

  end subroutine atl_pressure_1d_getElement

end module atl_eqn_euler_1d_derive_module
