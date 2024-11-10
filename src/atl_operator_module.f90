! Copyright (c) 2013-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014, 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2016-2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Jana Gericke <jana.gericke@student.uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
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

! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015, 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!--------------------------------------------
!    A O S - Array of structures layout new
!-------------------------------------------
! Access to get_point value output
! Access to get_element value output
! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt and Harald Klimach for
! University of Siegen.
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
!
! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the number of degrees of freedom for Q polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * A variable to store the number of degrees of freedom for a P tensor product
!   polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.

! The x ansatz degree is turned into the degree of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x ansatz degree is turned into the degree of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! **************************************************************************** !
!> This module provides the routine for applying operators. Currently it is
!! only implemented for 3D and needs to be extended to 2d
module atl_operator_module
  use, intrinsic :: iso_c_binding,  only: c_f_pointer, c_ptr
  use env_module,                   only: rk, long_k

  use tem_logging_module,           only: logUnit, lldebug
  use tem_varSys_module,            only: tem_varSys_type,                     &
    &                                     tem_varSys_op_type,                  &
    &                                     tem_varSys_solverData_evalElem_type
  use tem_aux_module,               only: tem_abort
  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use tem_topology_module,          only: tem_coordOfId, &
    &                                     tem_levelOf, tem_IdOfCoord
  use tem_geometry_module,          only: tem_CoordOfReal, tem_ElemSize, &
    &                                     tem_PosofId, tem_BaryOfId
  use tem_operation_var_module,     only: tem_opVar_fill_inputIndex, &
    &                                     tem_varSys_op_Data_type
  use tem_timer_module,             only: tem_startTimer, &
    &                                     tem_stopTimer
  use tem_grow_array_module,        only: grw_intArray_type, init, &
    &                                     append, truncate
  use tem_dyn_array_module,         only: dyn_intArray_type, init, append, &
    &                                     truncate

  use ply_leg_diff_module,          only: ply_calcdiff_leg,       &
    &                                     ply_calcdiff_leg_2d,    &
    &                                     ply_calcdiff_leg_1d,    &
    &                                     ply_calcDiff_leg_x_vec, &
    &                                     ply_calcDiff_leg_y_vec, &
    &                                     ply_calcDiff_leg_z_vec
  use ply_modg_basis_module,        only: ply_evalLegendreTensPoly

  use atl_varSys_module,            only: atl_varSys_data_type,          &
    &                                     atl_varSys_solverData_type,    &
    &                                     atl_get_new_varSys_data_ptr
  use atl_derive_module,            only: atl_generic_fromModal_getElement, &
    &                                     atl_derive_inputVar_type,         &
    &                                     atl_derive_fromModalData
  use atl_modg_1d_basis_module,     only: atl_evalLegendreTensPoly1d
  use atl_modg_2d_basis_module,     only: atl_evalLegendreTensPoly2d
  use atl_timer_module,             only: atl_cpl_elemTimers,    &
    &                                     atl_timerHandles

  implicit none

  private

  public :: atl_op_meansquare_forElement
  public :: atl_op_local_L2_mean_forElement
  public :: atl_op_divideVecByScal_forElement
  public :: atl_op_gradient_forElement, atl_op_gradient_forPoint
  public :: atl_op_GradientX_forElement
  public :: atl_op_GradientY_forElement
  public :: atl_op_GradientZ_forElement
  public :: atl_op_gradient_fromIndex
  public :: atl_set_opVar_getElement
  public :: atl_opVar_setupIndices


contains


  ! ************************************************************************ !
  subroutine atl_set_opVar_getElement( solData_evalElem, fun )
    ! ------------------------------------------------------------------ !
    !> Description on how to set the element retrieval function for stfuns.
    class(tem_varSys_solverData_evalElem_type), intent(in) :: solData_evalElem

    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    type(tem_varSys_op_type), intent(inout) :: fun
    ! ------------------------------------------------------------------ !
    type(tem_varSys_op_Data_type), pointer :: fptr
    type(atl_varSys_solverData_type), pointer :: fSDptr
    ! ------------------------------------------------------------------ !

    write(logunit(lldebug),*) "Setting different solver_bundle and ", &
      & fun%myPos
    call C_F_Pointer(fun%method_data, fptr)
    call c_f_pointer(solData_evalElem%solver_bundle, fSDptr)
    fptr%solver_bundle = atl_get_new_varSys_data_ptr( fSDptr )

    select case(trim(fun%operType))
    case ('division', 'div')
      fun%get_element => atl_op_division_forElement
    case ('divide_vector_by_scalar')
      fun%get_element => atl_op_divideVecByScal_forElement
    case ('gradient', 'grad')
      fun%get_element    => atl_op_Gradient_forElement
      fun%get_point      => atl_op_Gradient_forPoint
      fun%get_valOfIndex => atl_op_gradient_fromIndex
    case ('gradientX', 'gradX')
      fun%get_element    => atl_op_GradientX_forElement
      fun%get_point      => atl_op_Gradient_forPoint
      fun%get_valOfIndex => atl_op_gradient_fromIndex
    case ('gradientY', 'gradY')
      fun%get_element    => atl_op_GradientY_forElement
      fun%get_point      => atl_op_Gradient_forPoint
      fun%get_valOfIndex => atl_op_gradient_fromIndex
    case ('gradientZ', 'gradZ')
      fun%get_element    => atl_op_GradientZ_forElement
      fun%get_point      => atl_op_Gradient_forPoint
      fun%get_valOfIndex => atl_op_gradient_fromIndex
    case ('meansquare')
      fun%get_element => atl_op_meansquare_forElement
    case ('deviation')
      fun%get_element => atl_op_deviation_forElement
    case ('locall2mean')
      fun%get_element => atl_op_local_L2_mean_forElement
    case default
      write(logUnit(4),*) 'operType: '   &
        & // trim(fun%operType)          &
        & // ' not supported by ateles.'
    end select

  end subroutine atl_set_opVar_getElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_division( fun, varsys, tree, iElem, elemPos, nodalInput, &
    &                      nodalRes )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Current Element Index
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
    ! -------------------------------------------------------------------- !
    integer, parameter :: dividend = 1, divisor = 2
    integer :: iComp
    ! -------------------------------------------------------------------- !

    do iComp = 1, fun%nComponents
      nodalRes(:,iComp) =                     &
        & nodalInput(dividend)%data(:,iComp)  &
        & / nodalInput(divisor)%data(:,iComp)
    end do

  end subroutine atl_division
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_op_division_forElement( fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr_atl
    type(tem_varSys_op_Data_type), pointer :: fptr
    ! -------------------------------------------------------------------- !
    call C_F_POINTER(fun%method_data, fPtr)
    call C_F_POINTER(fPtr%solver_bundle, fPtr_atl)

    fnCalcPtr => atl_division

    call atl_generic_fromModal_getElement(  &
      & fun        = fun,                   &
      & varSys     = varSys,                &
      & elempos    = elempos,               &
      & time       = time,                  &
      & tree       = tree,                  &
      & nElems     = nElems,                &
      & nDofs      = nDofs,                 &
      & fnCalcPtr  = fnCalcPtr,             &
      & solverData = fPtr_atl%solverData,   &
      & res        = res                    )

  end subroutine atl_op_division_forElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_divideVecByScal( fun, varsys, tree, iElem, elemPos, &
    &                             nodalInput, nodalRes               )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Current Element Index
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
    ! -------------------------------------------------------------------- !
    integer, parameter :: dividend = 1, divisor = 2
    integer :: iComp
    ! -------------------------------------------------------------------- !

    do iComp = 1, fun%nComponents
      nodalRes(:,iComp) =                                                   &
        & nodalInput(dividend)%data(:,iComp) / nodalInput(divisor)%data(:, 1)
    end do

  end subroutine atl_divideVecByScal
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_op_divideVecByScal_forElement( fun, varsys, elempos, time, &
    &                                           tree, nElems, nDofs, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

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
    ! -------------------------------------------------------------------- !
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr_atl
    type(tem_varSys_op_Data_type), pointer :: fptr
    ! -------------------------------------------------------------------- !
    call C_F_POINTER(fun%method_data, fPtr)
    call C_F_POINTER(fPtr%solver_bundle, fPtr_atl)

    fnCalcPtr => atl_divideVecByScal

    call atl_generic_fromModal_getElement( &
      & fun        = fun,                  &
      & varSys     = varSys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr_atl%solverData,  &
      & res        = res                   )

  end subroutine atl_op_divideVecByScal_forElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_op_Gradient_forPoint( fun, varsys, point, time, tree, nPnts, &
    &                                  res                                    )
    ! -------------------------------------------------------------------- !
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
    type(tem_time_type), intent(in) :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(atl_varSys_data_type), pointer :: fPtr
    type(tem_varSys_op_data_type), pointer :: fPtr_temOp
    integer :: coord(4)
    integer(kind=long_k):: treeID
    integer :: iPoint, elemPos(1), level, iDof, iComp
    integer ::  nDofs,  maxPolyDegree, nComp
    real(kind=rk), allocatable :: input(:), inputModal(:,:)
    real(kind=rk), allocatable :: polyVal(:,:), val(:)
    real(kind=rk) :: bary(3), extent, halfwidth, MappedCoord(1,3)
    ! -------------------------------------------------------------------- !
    write(logUnit(10),*) 'Get the values of indices for gradient in' &
      & // 'variable '                                               &
      & // trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr_temOp )
    call C_F_POINTER( fPtr_temOp%solver_bundle, fPtr )

    nComp = fun%nComponents
    allocate(val(fun%nComponents))

    ! Loop over points
    do iPoint = 1, nPnts
      val = 0.0_rk

      ! Find out which element this point belongs to
      coord = tem_CoordOfReal(tree, point(iPoint,:), tree%global%maxLevel )
      treeId = tem_IdOfCoord(coord)
      elemPos(1) = tem_PosofId(treeId, tree%treeID)
      if (elemPos(1) == 0 ) then
        write(*,*) 'ERROR in atl_op_Gradient_forPoint: no element for this ', &
          &  'position ', point(iPoint,:), 'found, aborting...'
        call tem_abort()
      end if

      ! start element wise timer for LB weights
      call tem_startTimer( me          = atl_cpl_elemTimers, &
        &                  timerHandle = elemPos(1)          )
      call tem_startTimer( timerHandle = atl_timerHandles%gradient)

      level = tem_levelOf(treeId)
      bary = tem_BaryOfId(tree,TreeID)
      extent = tem_ElemSize( tree, TreeID)
      halfwidth = extent/2

      select case(fPtr%solverData%equationPtr%nDimensions)
        case(1)
          maxPolyDegree = fPtr%solverData%scheme_listPtr(level)%modg_1d     &
            &                                                  %maxpolydegree
        case(2)
          maxPolyDegree = fPtr%solverData%scheme_listPtr(level)%modg_2d     &
            &                                                  %maxpolydegree
        case(3)
          maxPolyDegree = fPtr%solverData%scheme_listPtr(level)%modg        &
            &                                                  %maxpolydegree
      end select

      nDofs = (maxPolyDegree + 1 )**fPtr%solverData%equationPtr%nDimensions

      allocate( polyval( nDofs, 1))
      allocate(input(nDofs*nComp))
      allocate(inputModal(ndofs, nComp))
      polyVal = 0_rk

      ! Get the gradient of element where the point is laying in,
      ! call get_elemt routin for that variable since it is pointing to
      ! atl_op_Gradient_forElement, hence in input are the modal gradient
      ! for that elemet
      call varSys%method%val(fun%myPos)%get_element( varSys  = varSys,  &
        &                                            elemPos = elemPos, &
        &                                            time    = time,    &
        &                                            tree    = tree,    &
        &                                            nelems  = 1,       &
        &                                            nDofs   = nDofs,   &
        &                                            res     = input    )
      ! Arrange the data of modal gradient
      do iComp = 1, nComp
        do iDof = 1, nDofs
            inputModal(iDof, iComp) =                    &
              & input((( 1-1)* ncomp* ndofs+( idof-1)* ncomp+icomp))
        end do
      end do

      ! So, get the Mapped Coordinates in reference element from the
      ! physical coordinates
      MappedCoord(1,1) = (point(iPoint,1)- bary(1))/halfwidth
      MappedCoord(1,2) = (point(iPoint,2)- bary(2))/halfwidth
      MappedCoord(1,3) = (point(iPoint,3)- bary(3))/halfwidth

      ! Now we Evaluate the polynomial at the given point using the modal
      ! values
      select case (fPtr%solverData%equationPtr%nDimensions)
      case(1)
        call atl_evalLegendreTensPoly1d(                                  &
          & coords        = MappedCoord,                                  &
          & ncoords       = 1,                                            &
          & maxPolyDegree = maxpolyDegree,                                &
          & basisType     = fPtr%solverData                               &
          &                         %scheme_listPtr(tree%global%minlevel) &
          &                         %modg_1d                              &
          &                         %basisType,                           &
          & polyVal       = polyVal                                       )
      case(2)
        call atl_evalLegendreTensPoly2d(                                  &
          & coords        = MappedCoord,                                  &
          & ncoords       = 1,                                            &
          & maxPolyDegree = maxpolyDegree,                                &
          & basisType     = fPtr%solverData                               &
          &                         %scheme_listPtr(tree%global%minlevel) &
          &                         %modg_2d                              &
          &                         %basisType,                           &
          & polyVal       = polyVal                                       )
      case(3)
        call ply_evalLegendreTensPoly(                                    &
          & coords        = MappedCoord,                                  &
          & ncoords       = 1,                                            &
          & maxPolyDegree = maxpolyDegree,                                &
          & basisType     = fPtr%solverData                               &
          &                         %scheme_listPtr(tree%global%minlevel) &
          &                         %modg                                 &
          &                         %basisType,                           &
          & polyVal       = polyVal                                       )
      end select

      ! Multiply with modal gradient values to get the exact  point value
      do iComp = 1, nComp
        do iDof = 1, nDofs
          val(icomp) = val(icomp)                      &
            & + polyVal(iDof,1) * inputModal(iDof, iComp)
        end do
      end do

      ! Store it in res
      res((iPoint-1)*nComp+1:(iPoint-1)*nComp+nComp) = val(:)

      deallocate( polyval)
      deallocate(input)
      deallocate(inputModal)

      call tem_stopTimer( me          = atl_cpl_elemTimers, &
        &                 timerHandle = elemPos(1)          )
      call tem_stopTimer( timerHandle = atl_timerHandles%gradient )
    end do

  end subroutine atl_op_Gradient_forPoint
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine takes in a variable and differentiates it in a modal way
  ! and returns the tensor
  subroutine atl_op_Gradient_forElement( fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    type(atl_varSys_data_type), pointer :: fPtr
    type(tem_varSys_op_data_type), pointer :: fPtr_temOp
    real(kind=rk), allocatable :: input(:), inputModalState(:,:,:)
    real(kind=rk), allocatable :: DiffInpTemp(:,:,:), DiffInp(:,:)
    integer :: nCompInp, iComp, iDof, iElem, i, j, coord(4), level
    real(kind=rk) :: elemLength
    integer :: maxPolyDegree
    ! This variable would have the actual number of degrees of freedom
    ! of the simulation
    integer :: numDofs
    integer :: nDimensions
    ! -------------------------------------------------------------------- !
    write(logUnit(4),*) 'Getting gradient by element for variable ', &
      &                  trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, fPtr_temOp )
    call C_F_POINTER( fPtr_temOp%solver_bundle, fPtr )

    nDimensions = fPtr%solverData%equationPtr%nDimensions
    select case(nDimensions)
    case(1)
      numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_1d%nDofs)
    case(2)
      numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_2d%nDofs)
    case(3)
      numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_3d%nDofs)
    end select

    nCompInp = varSys%method%val(fun%input_varPos(1))%nComponents
    allocate(input(numDofs*nElems*nCompInp))
    allocate(inputModalState(nElems, numdofs, nCompInp))
    allocate(DiffInpTemp(numdofs,nCompInp,nDimensions ))
    allocate(DiffInp(numdofs, fun%nComponents))

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = numDofs,                                   &
      & res     = input                                      )

    ! Arrange the data
    do iComp = 1, nCompInp
      do iDof = 1, numDofs
        do iElem = 1, nElems
          inputModalState(iElem,iDof, iComp) =                    &
            & input((( ielem-1)* ncompinp* numdofs+( idof-1)* ncompinp+icomp))
        end do
      end do
    end do

    !Differentiate the modal values
    do iElem = 1, nElems

      call tem_startTimer( timerHandle = atl_timerHandles%gradient )
      level = tem_levelOf(tree%treeID(elemPos(iElem)))

      ! The physical length of the element
      coord = tem_coordOfId( tree%treeid(iElem) )
      elemLength = tree%global%BoundingCubeLength / ( 2**(coord(4)) )

      select case (nDimensions)
      case(1)
        maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
          &                            %modg_1d               &
          &                            %maxpolydegree
        call ply_calcDiff_leg_1d( legCoeffs     = inPutModalState(iElem,:,:),  &
          &                       legCoeffsDiff = diffInpTemp(:,:,1),          &
          &                       maxPolyDegree = maxPolyDegree,               &
          &                       elemLength    = elemLength                   )

      case(2)
        maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
          &                            %modg_2d               &
          &                            %maxpolydegree
        call ply_calcDiff_leg_2d( legCoeffs     = inPutModalState(iElem,:,:),  &
          &                       legCoeffsDiff = diffInpTemp,                 &
          &                       maxPolyDegree = maxPolyDegree,               &
          &                       nVars         = nCompInp,                    &
          &                       elemLength    = elemLength                   )

      case(3)
        maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
          &                            %modg                  &
          &                            %maxpolydegree
        call ply_calcDiff_leg( legCoeffs     = inPutModalState(iElem,:,:),  &
          &                    legCoeffsDiff = diffInpTemp,                 &
          &                    maxPolyDegree = maxPolyDegree,               &
          &                    nVars         = nCompInp,                    &
          &                    elemLength    = elemLength                   )
      end select

      ! Put it in the output array
      do i = 1, nCompInp
        do j = 1, nDimensions
          DiffInp(:,(i-1)*nCompInp +j) = diffInpTemp( :,i,j)
        end do
      end do

      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res(( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp) = &
            & DiffInp(iDof, iComp)
        end do
      end do
      call tem_stopTimer( timerHandle = atl_timerHandles%gradient )

    end do

  end subroutine atl_op_Gradient_forElement
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> This routine takes in a variable and differentiates it in a modal way
  ! and returns the tensor
  subroutine atl_op_GradientX_forElement( fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    type(atl_varSys_data_type), pointer :: fPtr
    type(tem_varSys_op_data_type), pointer :: fPtr_temOp
    real(kind=rk), allocatable :: input(:), inputModalState(:,:,:)
    real(kind=rk), allocatable :: DiffInp(:,:)
    integer :: nCompInp, iComp, iDof, iElem, coord(4), level
    real(kind=rk) :: elemLength
    integer :: maxPolyDegree
    ! This variable would have the actual number of degrees of freedom
    ! of the simulation
    integer :: numDofs
    integer :: nDimensions
    ! -------------------------------------------------------------------- !
    write(logUnit(4),*) 'Getting gradientX normal by element for variable ', &
      &                  trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, fPtr_temOp )
    call C_F_POINTER( fPtr_temOp%solver_bundle, fPtr )

    nDimensions = fPtr%solverData%equationPtr%nDimensions
    select case(nDimensions)
    case(1)
      numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_1d%nDofs)
    case(2)
      numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_2d%nDofs)
    case(3)
      numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_3d%nDofs)
    end select
    !> The number of input components is equal to the number of output components
    !! This means nCompInp * nDimension = fun%nComponents
    !! In this routine both are the same, as we consider just one dimension
    nCompInp = varSys%method%val(fun%input_varPos(1))%nComponents
    allocate(input(numDofs*nElems*nCompInp))
    allocate(inputModalState(nElems, numdofs, nCompInp))
    allocate(DiffInp(numdofs, fun%nComponents))

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      &    varSys  = varSys,                                 &
      &    elemPos = elemPos,                                &
      &    time    = time,                                   &
      &    tree    = tree,                                   &
      &    nElems  = nElems,                                 &
      &    nDofs   = numDofs,                                &
      &    res     = input                                   )

    ! Arrange the data
    do iComp = 1, nCompInp
      do iDof = 1, numDofs
        do iElem = 1, nElems
          inputModalState(iElem,iDof, iComp) =                    &
            & input((( ielem-1)* ncompinp* numdofs+( idof-1)* ncompinp+icomp))
        end do
      end do
    end do

    !Differentiate the modal values
    Elemloop: do iElem = 1, nElems
      call tem_startTimer( timerHandle = atl_timerHandles%gradient )

      level = tem_levelOf(tree%treeID(elemPos(iElem)))

      ! The physical length of the element
      coord = tem_coordOfId( tree%treeid(iElem) )
      elemLength = tree%global%BoundingCubeLength / ( 2**(coord(4)) )

      select case (nDimensions)
      case(1)
        maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
          &                            %modg_1d               &
          &                            %maxpolydegree
        call ply_calcDiff_leg_1d( legCoeffs     = inPutModalState(iElem,:,:), &
          &                       legCoeffsDiff = diffInp,                    &
          &                       maxPolyDegree = maxPolyDegree,              &
          &                       elemLength    = elemLength                  )

      !!case(2)
      !!  maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
      !!    &                            %modg_2d               &
      !!    &                            %maxpolydegree
      !!  call ply_calcDiff_leg_2d(legCoeffs     = inPutModalState(iElem,:,:), &
      !!    &                      legCoeffsDiff = diffInp,                    &
      !!    &                      maxPolyDegree = maxPolyDegree,              &
      !!    &                      nVars         = nCompInp,                   &
      !!    &                      elemLength    = elemLength                  )

      case(3)
        maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
          &                            %modg                  &
          &                            %maxpolydegree
        call ply_calcDiff_leg_x_vec(                       &
          &    legCoeffs     = inPutModalState(iElem,:,:), &
          &    legCoeffsDiff = diffInp,                    &
          &    mPd           = maxPolyDegree,              &
          &    nVars         = nCompInp,                   &
          &    elemLength    = elemLength                  )
      end select

      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res(( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp) = &
            & DiffInp(iDof, iComp)
        end do
      end do
      call tem_stopTimer( timerHandle = atl_timerHandles%gradient )

    end do Elemloop

  end subroutine atl_op_GradientX_forElement
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> This routine takes in a variable and differentiates it in a modal way
  ! and returns the tensor
  subroutine atl_op_GradientY_forElement( fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    type(atl_varSys_data_type), pointer :: fPtr
    type(tem_varSys_op_data_type), pointer :: fPtr_temOp
    real(kind=rk), allocatable :: input(:), inputModalState(:,:,:)
    real(kind=rk), allocatable ::  DiffInp(:,:)
    integer :: nCompInp, iComp, iDof, iElem, coord(4), level
    real(kind=rk) :: elemLength
    integer :: maxPolyDegree
    ! This variable would have the actual number of degrees of freedom
    ! of the simulation
    integer :: numDofs
    integer :: nDimensions
    ! -------------------------------------------------------------------- !
    write(logUnit(4),*) 'Getting gradientY normal by element for variable ', &
      &                  trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, fPtr_temOp )
    call C_F_POINTER( fPtr_temOp%solver_bundle, fPtr )

    nDimensions = fPtr%solverData%equationPtr%nDimensions
    select case(nDimensions)
    case(2)
      numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_2d%nDofs)
    case(3)
      numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_3d%nDofs)
    end select

    nCompInp = varSys%method%val(fun%input_varPos(1))%nComponents
    allocate(input(numDofs*nElems*nCompInp))
    allocate(inputModalState(nElems, numdofs, nCompInp))
    allocate(DiffInp(numdofs, fun%nComponents))

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      &    varSys  = varSys,                                 &
      &    elemPos = elemPos,                                &
      &    time    = time,                                   &
      &    tree    = tree,                                   &
      &    nElems  = nElems,                                 &
      &    nDofs   = numDofs,                                &
      &    res     = input                                   )

    ! Arrange the data
    do iComp = 1, nCompInp
      do iDof = 1, numDofs
        do iElem = 1, nElems
          inputModalState(iElem,iDof, iComp) =                    &
            & input((( ielem-1)* ncompinp* numdofs+( idof-1)* ncompinp+icomp))
        end do
      end do
    end do

    !Differentiate the modal values
    do iElem = 1, nElems

      call tem_startTimer( timerHandle = atl_timerHandles%gradient )
      level = tem_levelOf(tree%treeID(elemPos(iElem)))

      ! The physical length of the element
      coord = tem_coordOfId( tree%treeid(iElem) )
      elemLength = tree%global%BoundingCubeLength / ( 2**(coord(4)) )

      ! We do not need to compute the gardients for a 1D test case
      ! In y direction. Therefore there is no 1D test case here
      ! TODO change the routine for the 2D gradient computation
      select case (nDimensions)
      !!case(2)
      !!  maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
      !!    &                            %modg_2d               &
      !!    &                            %maxpolydegree
      !!  call ply_calcDiff_leg_2d( legCoeffs     = inPutModalState(iElem,:,:), &
      !!    &                       legCoeffsDiff = diffInp,                    &
      !!    &                       maxPolyDegree = maxPolyDegree,              &
      !!    &                       nVars         = nCompInp,                   &
      !!    &                       elemLength    = elemLength                  )

      case(3)
        maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
          &                            %modg                  &
          &                            %maxpolydegree
        call ply_calcDiff_leg_y_vec(                       &
          &    legCoeffs     = inPutModalState(iElem,:,:), &
          &    legCoeffsDiff = diffInp,                    &
          &    mPd           = maxPolyDegree,              &
          &    nVars         = nCompInp,                   &
          &    elemLength    = elemLength                  )
      case default
        call tem_abort( 'GradientY is just provided in 2D and 3D. ' &
          & // ' Please your configuration, stopping ...'           )

      end select

      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res(( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp) = &
            & DiffInp(iDof, iComp)
        end do
      end do
      call tem_stopTimer( timerHandle = atl_timerHandles%gradient )

    end do

  end subroutine atl_op_GradientY_forElement
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> This routine takes in a variable and differentiates it in a modal way
  ! and returns the tensor
  subroutine atl_op_GradientZ_forElement( fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
    ! -------------------------------------------------------------------- !
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
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(atl_varSys_data_type), pointer :: fPtr
    type(tem_varSys_op_data_type), pointer :: fPtr_temOp
    real(kind=rk), allocatable :: input(:), inputModalState(:,:,:)
    real(kind=rk), allocatable :: DiffInp(:,:)
    integer :: nCompInp, iComp, iDof, iElem, coord(4), level
    real(kind=rk) :: elemLength
    integer :: maxPolyDegree
    ! This variable would have the actual number of degrees of freedom
    ! of the simulation
    integer :: numDofs
    integer :: nDimensions
    ! -------------------------------------------------------------------- !
    write(logUnit(6),*) 'Getting gradientZ normal by element for variable ', &
      &                  trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, fPtr_temOp )
    call C_F_POINTER( fPtr_temOp%solver_bundle, fPtr )

    nDimensions = fPtr%solverData%equationPtr%nDimensions
    numDofs = maxval(fPtr%solverData%polyProjectPtr(:)%body_3d%nDofs)

    nCompInp = varSys%method%val(fun%input_varPos(1))%nComponents
    allocate(input(numDofs*nElems*nCompInp))
    allocate(inputModalState(nElems, numdofs, nCompInp))
    allocate(DiffInp(numdofs, fun%nComponents))

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      &    varSys  = varSys,                                 &
      &    elemPos = elemPos,                                &
      &    time    = time,                                   &
      &    tree    = tree,                                   &
      &    nElems  = nElems,                                 &
      &    nDofs   = numDofs,                                &
      &    res     = input                                   )

    ! Arrange the data
    do iComp = 1, nCompInp
      do iDof = 1, numDofs
        do iElem = 1, nElems
          inputModalState(iElem,iDof, iComp) =                    &
            & input((( ielem-1)* ncompinp* numdofs+( idof-1)* ncompinp+icomp))
        end do
      end do
    end do

    !Differentiate the modal values
    Elemloop: do iElem = 1, nElems
      call tem_startTimer( timerHandle = atl_timerHandles%gradient )

      level = tem_levelOf(tree%treeID(elemPos(iElem)))

      ! The physical length of the element
      coord = tem_coordOfId( tree%treeid(iElem) )
      elemLength = tree%global%BoundingCubeLength / ( 2**(coord(4)) )

      if (nDimensions == 3) then
        maxPolyDegree = fPtr%solverData%scheme_listPtr(level) &
           &                            %modg                 &
           &                            %maxpolydegree
        call ply_calcDiff_leg_z_vec(                       &
          &    legCoeffs     = inPutModalState(iElem,:,:), &
          &    legCoeffsDiff = diffInp,                    &
          &    mPd           = maxPolyDegree,              &
          &    nVars         = nCompInp,                   &
          &    elemLength    = elemLength                  )

        do iDof = 1, nDofs
          do iComp = 1, fun%nComponents
            res(( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp) = &
              & DiffInp(iDof, iComp)
          end do
        end do

      else
        call tem_abort( 'GradientZ is just supported in 3D. '   &
          & // ' Please check your configuration, stopping ...' )
      end if
    end do Elemloop
    call tem_stopTimer( timerHandle = atl_timerHandles%gradient )

  end subroutine atl_op_GradientZ_forElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_op_gradient_fromIndex( fun, varsys, time, iLevel, idx, &
    &                                   idxLen, nVals, res              )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: most times nVals, if contiguous arrays are used it depends
    !! on the number of first indices
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(atl_varSys_data_type), pointer :: fPtr
    type(atl_varSys_data_type), pointer :: fPtr_in
    type(tem_varSys_op_data_type), pointer :: fPtr_temOp
    integer :: iPoint, iDof, iComp, iELem
    integer ::  nDofs,  maxPolyDegree, nComp, elemPos(1)
    real(kind=rk), allocatable :: input(:), inputModal(:,:)
    real(kind=rk), allocatable :: polyVal(:,:), val(:)
    real(kind=rk) :: coord(1,3)
    logical :: found
    type(dyn_intArray_type) :: unique_elemPos
    type(grw_intArray_type) :: grw_nCoordPerElem
    logical :: wasAdded
    integer :: loc_level
    integer :: pos, global_count, buf_start, buf_end
    integer, allocatable :: pnt_pos(:)
    ! -------------------------------------------------------------------- !
    write(logUnit(10),*) 'Get the values of indices for gradient in' &
      & // 'variable '                                               &
      & // trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr_temOp )
    call C_F_POINTER( fPtr_temOp%solver_bundle, fPtr )

    ! now we search recursivly for the input variable where the main info
    ! like elemPos and points are stored
    ! this variable should not have further input variables and it either
    ! a state variable of spacetime function
    call get_statePtr( varPos = fun%myPos, &
      &                varSys = varSys,    &
      &                found = found,      &
      &                fPtr = fPtr_in      )

    nComp = fun%nComponents
    allocate(val(fun%nComponents))

    ! sort the point array elementwise
    ! loop over points and store them per element
    call init( unique_elemPos )
    call init( grw_nCoordPerElem )
    do iPoint = 1, nVals
      elemPos = fPtr_in%pointData%pntLvl(iLevel)%elemPos%val(idx(iPoint))
      ! get an unqiue array of elempos
      call append( me       = unique_elemPos, &
        &          val      = elemPos(1),     &
        &          pos      = pos,            &
        &          wasAdded = wasAdded        )
      ! count points per elempos
      if (wasAdded) then
        ! add new ElemPos
        call append(grw_nCoordPerElem, 1)
      else
        ! count up the elempos which is already there
        grw_nCoordPerElem%val(Pos) = grw_nCoordPerElem%val(Pos) + 1
      end if
    end do
    call truncate( unique_elemPos )
    call truncate( grw_nCoordPerElem )

    ! now, store the position in the coord array accoding to the sequence in the
    ! unquie element array
    allocate(pnt_pos(nVals))
    pnt_pos = -1
    global_count = 1
    do iElem = 1, unique_elemPos%nVals
      do iPoint = 1, nVals
        if (unique_elemPos%val(iElem) == fPtr_in%pointData%pntLvl(iLevel)     &
          &                              %elemPos%val(idx(iPoint))       ) then
          pnt_pos(global_count) = iPoint
          global_count = global_count +1
        end if
      end do
    end do

    ! now we go the sequence of the element
    Elemloop: do iElem = 1, unique_elemPos%nVals
      call tem_startTimer( timerHandle = atl_timerHandles%gradient )

      elemPos(1) = unique_elemPos%val(iElem)
      ! start element wise timer for LB weights
      call tem_startTimer( me          = atl_cpl_elemTimers, &
        &                  timerHandle = elemPos(1)          )
      ! upper and lower bound of array access according to elemPos offset/
      ! counter grw_nCoordPerElem
      buf_start = sum(grw_nCoordPerElem%val(1:iElem-1))+1
      buf_end = sum(grw_nCoordPerElem%val(1:iElem))

      ! get teh correct local level for the points
      loc_level = fPtr_in%pointData%pntLvl(iLevel)%pntLevel%val(idx(pnt_pos(buf_start)))

      ! get the correct degree to allocate the correct arrray size
      select case(fPtr%solverData%equationPtr%nDimensions)
        case(1)
          maxPolyDegree = fPtr%solverData%scheme_listPtr(loc_level)%modg_1d &
            &                       %maxpolydegree
        case(2)
          maxPolyDegree = fPtr%solverData%scheme_listPtr(loc_level)%modg_2d &
            &                       %maxpolydegree
        case(3)
          maxPolyDegree = fPtr%solverData%scheme_listPtr(loc_level)%modg &
            &                       %maxpolydegree
      end select
      nDofs = (maxPolyDegree + 1 )**fPtr%solverData%equationPtr%nDimensions

      allocate(polyval( nDofs, 1))
      allocate(input(nDofs*nComp))
      allocate(inputModal(ndofs, nComp))
      input = 0_rk
      inputModal = 0_rk


      ! Get the gradient of element where the point is laying in,
      ! call get_elemt routin for that variable since it is pointing to
      ! atl_op_Gradient_forElement, hence in input are the modal gradient
      ! for that elemet
      call varSys%method%val(fun%mypos)%get_element( &
        &  varSys  = varSys,                         &
        &  elemPos = elemPos,                        &
        &  time    = time,                           &
        &  tree    = fPtr%solverData%tree,           &
        &  nelems  = 1,                              &
        &  nDofs   = nDofs,                          &
        &  res     = input                           )
      ! Arrange the data of modal gradient
      do iComp = 1, nComp
        do iDof = 1, nDofs
            inputModal(iDof, iComp) =                    &
              & input((( 1-1)* ncomp* ndofs+( idof-1)* ncomp+icomp))
        end do
      end do

      Pointloop: do iPoint = buf_start, buf_end
        polyVal = 0_rk
        val = 0.0_rk

        ! So, get the all physical coordinates
        ! distinguish if we have an array of index or we have contingous memory
        ! acces where index are always first entries!
        if ( .not. present(idxLen) ) then
          coord(1,1)=fPtr_in%pointData%pntLvl(iLevel)%grwPnt%coordX &
            &                    %val(idx(pnt_Pos(iPoint)))
          coord(1,2)=fPtr_in%pointData%pntLvl(iLevel)%grwPnt%coordY &
            &                    %val(idx(pnt_Pos(iPoint)))
          coord(1,3)=fPtr_in%pointData%pntLvl(iLevel)%grwPnt%coordZ &
            &                    %val(idx(pnt_Pos(iPoint)))
        else ! idxLen is present
          ! integer which stores the position in the coord array
          !> to do
          call tem_abort( 'idxLen in atl_op_gradient_fromIndex not yet' &
            &             // 'implemented'                              )
        end if
        ! Now we Evaluate the polynomial at the given point using the modal

        ! values
        select case (fPtr%solverData%equationPtr%nDimensions)
        case(1)
          call atl_evalLegendreTensPoly1d(                              &
            & coords        = coord,                                    &
            & ncoords       = 1,                                        &
            & maxPolyDegree = maxpolyDegree,                            &
            & basisType     = fPtr%solverData%scheme_listPtr(loc_level) &
            &                     %modg_1d%basisType,                   &
            & polyVal       = polyVal                                   )
        case(2)
          call atl_evalLegendreTensPoly2d(                              &
            & coords        = coord,                                    &
            & ncoords       = 1,                                        &
            & maxPolyDegree = maxpolyDegree,                            &
            & basisType     = fPtr%solverData%scheme_listPtr(loc_level) &
            &                         %modg_2d%basisType,               &
            & polyVal       = polyVal                                   )
        case(3)
          call ply_evalLegendreTensPoly(                                &
            & coords        = coord,                                    &
            & ncoords       = 1,                                        &
            & maxPolyDegree = maxpolyDegree,                            &
            & basisType     = fPtr%solverData%scheme_listPtr(loc_level) &
            &                         %modg%basisType,                  &
            & polyVal       = polyVal                                   )
        end select

        ! Multiply with modal gradient values to get the exact  point value
        do iComp = 1, nComp
          do iDof = 1, nDofs
            val(icomp) = val(icomp)                      &
              & + polyVal(iDof,1) * inputModal(iDof, iComp)
          end do
        end do

        ! Store it in res
        res( (pnt_Pos(iPoint)-1)*nComp +1 :               &
          &     (pnt_Pos(iPoint)-1)*nComp +nComp ) = val(:)

      end do Pointloop

      call tem_stopTimer( timerHandle = atl_timerHandles%gradient )
      call tem_stopTimer( me          = atl_cpl_elemTimers, &
        &                 timerHandle = elemPos(1)          )
      deallocate(polyval)
      deallocate(inputModal)
      deallocate(input)
    end do Elemloop

    deallocate(val)

  end subroutine atl_op_gradient_fromIndex
  ! ************************************************************************ !


  ! ************************************************************************ !
  ! This routine looks recursivley for the input variable which is either a
  ! state or a spacetimefunction, the criteria is that this variable depends on
  ! no more input variables and there for we get the pointer to find infomation
  ! stored there.
  recursive subroutine get_statePtr(varPos, varSys, found, fPtr)
    ! -------------------------------------------------------------------- !
  integer, intent(in) :: varPos
  !> The variable system to obtain the variable from.
  type(tem_varSys_type), intent(in) :: varSys
  logical, intent(out) :: found
  type(atl_varSys_data_type), pointer, intent(out) :: fPtr
    ! -------------------------------------------------------------------- !
  type(tem_varSys_op_type) :: fun
  integer :: nInput
    ! -------------------------------------------------------------------- !
  found = .false.

  fun = varSys%method%val(varPos)

  if (fun%nInputs > 0) then
    do nInput=1, fun%nInputs
      call get_statePtr(fun%input_varPos(ninput), varSys, found, fPtr)
      if (found) exit
    end do
  else
    write(logUnit(10),*) 'state variable is ',  trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( varSys%method%val(fun%myPos)%method_data, fPtr )
    found = .true.
  end if
  end subroutine
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine takes points coordinates, pass them tp the input variables
  !! the opertaion depends and return indices
  subroutine atl_opVar_setupIndices( fun, varSys, point, offset_bit, &
    &                                iLevel, tree, nPnts, idx        )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> List of space coordinate points to store as growing array in
    !! method_data
    real(kind=rk), intent(in) :: point(:,:)

    !> Offset bit encoded as character for every point.
    !!
    !! Offset integer coord(3) is converted into a character with
    !! offset_bit = achar( (coord(1)+1) + (coord(2)+1)*4 + (coord(3)+1)*16 )
    !! Backward transformation form character to 3 integer:
    !! coord(1) = mod(ichar(offset_bit),4) - 1
    !! coord(2) = mod(ichar(offset_bit),16)/4 - 1
    !! coord(3) = ichar(offset_bit)/16 - 1
    !!
    !! If not present default is to center i.e offset_bit = achar(1+4+16)
    character, optional, intent(in) :: offset_bit(:)

    !> Level to which input points belong to
    integer, intent(in) :: iLevel

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of points to add in method_data of this variable
    integer, intent(in) :: nPnts

    !> Index of points in the growing array and variable val array.
    !! Size: nPoints
    !!
    !! This must be stored in boundary or source depends on who
    !! calls this routine.
    !! This index is required to return a value using getValOfIndex.
    integer, intent(out) :: idx(:)
    ! -------------------------------------------------------------------- !
    type(atl_varSys_data_type), pointer :: fPtr
    integer :: iPnt, iDep
    type(grw_intArray_type), allocatable :: inputIndex_loc(:)
    integer, allocatable :: idxPerPnt(:)
    ! -------------------------------------------------------------------- !
    write(logUnit(10),*) 'setup indices for the points of operation' &
      & // ' variable '                                              &
      & // trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr )

    ! allcoate the index array for all inpits
    if (.not. allocated(fPtr%opData%input_pntIndex)) then
      allocate( fPtr%OpData%input_pntIndex(fun%nInputs) )
    end if

    ! allocate temporary inputIndex with size of nInputs and initialize
    ! growing array with length nPnts
    allocate(inputIndex_loc(fun%nInputs))

    ! Now fill in the index arrays for the inputs
    call tem_opVar_fill_inputIndex( fun        = fun,           &
      &                             varSys     = varSys,        &
      &                             point      = point,         &
      &                             offset_bit = offset_bit,    &
      &                             iLevel     = iLevel,        &
      &                             tree       = tree,          &
      &                             nPnts      = nPnts,         &
      &                             inputIndex = inputIndex_loc )

    ! fill the index array of the derived variable, it starts with the first
    ! entry in this call = nVals_prev and is continguous until nVals_prev+nVals
    allocate(idxPerPnt(fun%nInputs))
    idx = 0
    do iPnt = 1, nPnts
      do iDep = 1, fun%nInputs
        idxPerPnt(iDep) = inputIndex_loc(iDep)%val(iPnt)
      end do
      ! set index only when any of dependent variable has valid index
      if (any(idxPerPnt > 0)) then
        do iDep = 1, fun%nInputs
          call append(me  = fPtr%opData%input_pntIndex(iDep)%indexLvl(iLevel), &
            &         val = inputIndex_loc(iDep)%val(iPnt)                     )
        end do!iDep
        ! set index to last position in input_pntIndex of dep var 1 of
        ! indexLvl of iLevel
        idx(iPnt) = fPtr%opData%input_pntIndex(1)%indexLvl(iLevel)%nVals
      end if
    end do!iPnt

    do iDep = 1, fun%nInputs
      call truncate (fPtr%opData%input_pntIndex(iDep)%indexLvl(iLevel) )
    end do
  end subroutine atl_opVar_setupIndices
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_op_meansquare_forElement (fun, varsys, elempos, time, tree, &
    &                                       nElems, nDofs, res )
    ! ---------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    ! Integer for the loop
    integer       :: iDof,iElem,iComp
    integer :: nComps
    integer :: firstdof

    !(ijk) indices of the coeficients
    integer       :: ansFuncX,ansFuncY,ansFuncZ

    ! maximal degree of the polynomials
    integer       :: maxdegree,PolyOrd
    ! -------------------------------------------------------------------- !

    nComps = varSys%method%val(fun%input_varPos(1))%nComponents

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = res                                        )
    PolyOrd   = nint(nDofs**(1.0_rk/3.0_rk))
    maxdegree = PolyOrd - 1

    do iElem    = 1 , nElems
      do iComp = 1,nComps
        firstdof = ((ielem-1)*ncomps*ndofs+(1-1)*ncomps+icomp)

        ansFuncX    = 1
        ansFuncY    = 1
        ansFuncZ    = 1

        ! Square the first degree of freedom for each element and
        ! component.
        res(firstdof) = res(firstdof)**2

        ! Get to the second degree of freedom
  if (ansfuncx .ne. (maxdegree + 1)) then
    ! next x index
    ansfuncx = ansfuncx + 1
  elseif (ansfuncy .ne. (maxdegree + 1)) then
    ! next y index
    ansfuncx = 1
    ansfuncy = ansfuncy + 1
  else
    ! next z index
    ansfuncx = 1
    ansfuncy = 1
    ansfuncz = ansfuncz + 1
  end if

        do iDof = 2, nDofs
          res(firstdof) = res(firstdof) &
            &           + ( res(((ielem-1)*ncomps*ndofs+(idof-1)*ncomps+icomp))**2 ) &
            &            / ( (2*ansFuncX - 1)                                  &
            &               *(2*ansFuncY - 1)                                  &
            &               *(2*ansFuncZ - 1) )
          ! Set higher mode to 0 after adding it to the (reduced) first
          ! mode:
          res(((ielem-1)*ncomps*ndofs+(idof-1)*ncomps+icomp)) = 0.0_rk
  if (ansfuncx .ne. (maxdegree + 1)) then
    ! next x index
    ansfuncx = ansfuncx + 1
  elseif (ansfuncy .ne. (maxdegree + 1)) then
    ! next y index
    ansfuncx = 1
    ansfuncy = ansfuncy + 1
  else
    ! next z index
    ansfuncx = 1
    ansfuncy = 1
    ansfuncz = ansfuncz + 1
  end if
        end do
      end do
    end do
  end subroutine atl_op_meansquare_forElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_op_deviation_forElement (fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res )
    ! ---------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: mean
    ! Integer for the loop
    integer       :: iDof,iElem,iComp
    integer :: nComps
    integer :: firstdof

    !(ijk) indices of the coeficients
    integer       :: ansFuncX,ansFuncY,ansFuncZ

    ! maximal degree of the polynomials
    integer       :: maxdegree,PolyOrd
    ! -------------------------------------------------------------------- !

    nComps = varSys%method%val(fun%input_varPos(1))%nComponents

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = res                                        )
    PolyOrd   = nint(nDofs**(1.0_rk/3.0_rk))
    maxdegree = PolyOrd - 1

    do iElem    = 1 , nElems
      do iComp = 1,nComps
        firstdof = ((ielem-1)*ncomps*ndofs+(1-1)*ncomps+icomp)

        ansFuncX    = 1
        ansFuncY    = 1
        ansFuncZ    = 1

        mean = max(res(firstdof), epsilon(mean))

        ! Get to the second degree of freedom
  if (ansfuncx .ne. (maxdegree + 1)) then
    ! next x index
    ansfuncx = ansfuncx + 1
  elseif (ansfuncy .ne. (maxdegree + 1)) then
    ! next y index
    ansfuncx = 1
    ansfuncy = ansfuncy + 1
  else
    ! next z index
    ansfuncx = 1
    ansfuncy = 1
    ansfuncz = ansfuncz + 1
  end if

        res(firstdof) = 0.0_rk
        do iDof = 2, nDofs
          res(firstdof) = res(firstdof) &
            &           + abs(res(((ielem-1)*ncomps*ndofs+(idof-1)*ncomps+icomp)))
          ! Set higher mode to 0 after adding it to the (reduced) first
          ! mode:
          res(((ielem-1)*ncomps*ndofs+(idof-1)*ncomps+icomp)) = 0.0_rk
  if (ansfuncx .ne. (maxdegree + 1)) then
    ! next x index
    ansfuncx = ansfuncx + 1
  elseif (ansfuncy .ne. (maxdegree + 1)) then
    ! next y index
    ansfuncx = 1
    ansfuncy = ansfuncy + 1
  else
    ! next z index
    ansfuncx = 1
    ansfuncy = 1
    ansfuncz = ansfuncz + 1
  end if
        end do
        res(firstdof) = res(firstdof)/mean
      end do
    end do
  end subroutine atl_op_deviation_forElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_op_local_L2_mean_forElement( fun, varsys, elempos, time,  &
    &                                         tree, nElems, nDofs, res     )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !

    call atl_op_meansquare_forElement (fun     = fun,       &
      &                                varsys  = varsys,    &
      &                                elempos = elempos,   &
      &                                time    = time,      &
      &                                tree    = tree,      &
      &                                nElems  = nElems,    &
      &                                nDofs   = nDofs,     &
      &                                res     = res       )

    res = sqrt(res)

  end subroutine atl_op_local_L2_mean_forElement
  ! ************************************************************************ !

end module atl_operator_module
