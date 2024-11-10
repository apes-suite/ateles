! Copyright (c) 2014-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016, 2018-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
!> This module contains types and routines to work with the ateles variable
!! system. This variable system is build upon the treelm routines, but contains
!! solver specific enhancements.
module atl_varSys_module


  use, intrinsic :: iso_c_binding,  only: c_ptr, c_f_pointer, c_loc, c_null_ptr
  use env_module,                   only: rk, long_k, labelLen

  use aotus_module,                 only: flu_State, aot_get_val

  use tem_logging_module,           only: logUnit, lldebug
  use tem_varSys_module,            only: tem_varSys_type,               &
    &                                     tem_varSys_op_type,            &
    &                                     tem_varSys_solverData_evalElem_type
  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use tem_geometry_module,          only: tem_CoordOfReal, &
    &                                     tem_PosofId,     &
    &                                     tem_BaryOfId,    &
    &                                     tem_ElemSize
  use tem_topology_module,          only: tem_IdOfCoord, &
    &                                     tem_levelOf
  use tem_shape_module,             only: tem_global_shape
  use tem_variable_module,          only: tem_variable_type, &
    &                                     tem_append_solverVar_method
  use tem_aux_module,               only: tem_abort
  use tem_grow_array_module,        only: grw_intArray_type,  &
    &                                     append, init, truncate, destroy
  use tem_dyn_array_module,         only: dyn_intArray_type, init, append, &
    &                                     truncate, destroy
  use tem_spacetime_fun_module,     only: tem_st_fun_listelem_type
  use tem_pointData_module,         only: tem_pointData_list_type, &
    &                                     init, append, truncate
  use tem_operation_var_module,     only: tem_varSys_op_data_type
  use tem_timer_module,             only: tem_startTimer, &
    &                                     tem_stopTimer

  use ply_poly_project_module,      only: ply_poly_project_type
  use ply_oversample_module,        only: ply_convertFromOversample
  use ply_poly_project_module,      only: ply_poly_project_n2m
  use ply_modg_basis_module,        only: ply_evalLegendreTensPoly

  use atl_equation_module,          only: atl_equations_type
  use atl_modg_1d_basis_module,     only: atl_evalLegendreTensPoly1d
  use atl_modg_2d_basis_module,     only: atl_evalLegendreTensPoly2d
  use atl_scheme_module,            only: atl_scheme_type
  use atl_kerneldata_module,        only: atl_statedata_type, &
    &                                     atl_kerneldata_type
  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_reference_element_module, only: atl_ref_in_elempos
  use atl_timer_module,             only: atl_cpl_elemTimers
  use atl_legpolyvar_module,        only: atl_legpolyvar_type, &
    &                                     atl_legpolyvar_load, &
    &                                     atl_legpolyvar_append
  use ply_dof_module,               only: Q_space, &
    &                                     P_space

  implicit none

  private

  public :: atl_varSys_data_type
  public :: atl_varSys_solverData_type
  public :: atl_varSys_getStateForElement
  public :: atl_varSys_getStateForPoint
  public :: atl_varSys_setupStateIndices
  public :: atl_varSys_getStateValofIndex
  public :: atl_varSys_load_user
  public :: atl_init_varSys_solverData
  public :: atl_get_new_varSys_data_ptr
  public :: atl_create_fortranVar
  public :: atl_set_stFun_getElement

  !> Solver-specific structure for solver and source term variables.
  type atl_varSys_data_type
    type(atl_varSys_solverData_type), pointer :: solverData
    !>the point_datas need to be stored levelwise
    type(tem_pointData_list_type) :: pointData
    !> data array for operation or derived varibales
    !! consists the index arrys for points stored in the
    !! poingtData of input variable
    !! size is number of input variables
    type(tem_varSys_op_data_type) :: opData
  end type

  !> Solver-specific structure for solver and source term variables.
  !!
  !! Method data structures are used to provide access to solver structures via
  !! function pointers routed through the variable system. This structure is
  !! also used to provide access to structures needed to evaluate source
  !! terms.
  type atl_varSys_solverData_type
    type(atl_scheme_type), pointer :: scheme_listPtr(:) => null()
    type(atl_statedata_type), pointer :: statedata_listPtr(:) => null()
    type(atl_kerneldata_type), pointer :: kerneldata_listPtr(:) => null()
    type(atl_equations_type), pointer :: equationPtr => null()
    type(ply_poly_project_type), pointer :: polyProjectPtr(:) => null()
    !> Levelwise mesh information used for creating chebychev coordinates.
    type(atl_cube_elem_type), pointer :: mesh_listPtr(:) => null()
    integer, pointer :: poly_proj_posPtr(:) => null()
    integer, pointer :: levelPointer(:) => null()
    type(treelmesh_type), pointer :: tree => null()
  end type atl_varSys_solverData_type

  type atl_varSys_solverVar_type
    character(len=labelLen) :: vartype
    type(atl_legpolyvar_type) :: legpolyvar
  end type atl_varSys_solverVar_type


contains


  ! ************************************************************************ !
  !> Routine to get a pointer to a new instance of atl_varSys_data_type to be
  !! used as method data for a variable in the variable system.
  !!
  !! A new instance is allocated and a c_ptr to this type is returned. We
  !! currently don't need a reference to this instance besides the one that
  !! will be stored in the variable syste, thus we don't keep one here.
  function atl_get_new_varSys_data_ptr( solverData ) result(resPtr)
    ! -------------------------------------------------------------------- !
    !> The prototype is used to initialize the new instance.
    type(atl_varSys_solverData_type), intent(in), target :: solverData
    !> Pointer to the newly created instance.
    type(c_ptr) :: resPtr
    ! -------------------------------------------------------------------- !
    !> Local variable to allocate a new instance.
    type(atl_varSys_Data_type), pointer :: res
    ! -------------------------------------------------------------------- !

    allocate(res)
    res%solverData => solverData

    !!! the point_data_type need to be levelwise
    !!init( res%pointData%pointX )
    !!init( res%pointData%pointY )
    !!init( res%pointData%pointZ )
    !!init( res%pointData%elemPos )
    resPtr = c_loc(res)

  end function atl_get_new_varSys_data_ptr
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_init_varSys_solverData( this, tree, equation, schemeList,   &
    &                                    statedataList, kerneldataList,      &
    &                                    levelPointer,                       &
    &                                    meshList, polyProjPos, polyProjList )
    ! -------------------------------------------------------------------- !
    type(atl_varSys_solverData_type), intent(inout) :: this
    type(treelmesh_type),target, intent(in) :: tree
    type(atl_equations_type), target, intent(in) :: equation
    type(atl_scheme_type), target, intent(in) &
      & :: schemeList(tree%global%minLevel:)
    type(atl_statedata_type), target, intent(in) :: statedataList( &
      & tree%global%minLevel:)
    type(atl_kerneldata_type), target, intent(in) :: kerneldataList( &
      & tree%global%minLevel:)
    !> Levelwise mesh information used for creating chebychev coordinates.
    type(atl_cube_elem_type), target, intent(in) :: meshList(tree%global%minLevel:)
    integer, target, intent(in) :: polyProjPos(tree%global%minLevel:)
    integer, target, intent(in) :: levelPointer(:)
    type(ply_poly_project_type), target, intent(in) :: polyProjList(:)
    ! -------------------------------------------------------------------- !
    this%equationPtr => equation
    this%scheme_listPtr => schemeList
    this%statedata_listPtr => statedataList
    this%kerneldata_listPtr => kerneldataList
    this%mesh_listPtr => meshList
    this%polyProjectPtr => polyProjList
    this%poly_proj_posPtr => polyProjPos
    this%levelPointer => levelPointer
    this%tree => tree
  end subroutine atl_init_varSys_solverData
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_set_stFun_getElement(solData_evalElem, fun)
    ! -------------------------------------------------------------------- !
    !> Description on how to set the element retrieval function for stfuns.
    class(tem_varSys_solverData_evalElem_type), intent(in) :: solData_evalElem

    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    type(tem_varSys_op_type), intent(inout) :: fun
    ! -------------------------------------------------------------------- !
    type(tem_st_fun_listElem_type), pointer :: fptr
    type(atl_varSys_solverData_type), pointer :: fSDptr
    ! -------------------------------------------------------------------- !

    write(logunit(lldebug),*) "Setting different solver_bundle and" &
      & // " get_element routine for variable at position ",        &
      & fun%myPos
    call C_F_Pointer(fun%method_data, fptr)
    call c_f_pointer(solData_evalElem%solver_bundle, fSDptr)
    fptr%solver_bundle = atl_get_new_varSys_data_ptr( fSDptr )

    if (all(fptr%val(:)%fun_kind == 'const')) then
      fun%get_element => atl_generic_fromConst_getElement
    else
      fun%get_element => atl_generic_fromNodal_getElement
    end if

  end subroutine atl_set_stFun_getElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> To obtain values of a given variable, it is necessary to state the
  !! treeID and time at which the variable should be evaluated.
  !! The interface is nDofs values to cover the all degrees of freedoms
  !! in the element.
  !! Of course the variable system itself also needs to be passed in, to
  !! allow the computation of other derived quantities as needed.
  !! The method description itself is passed in automatically, and has not
  !! to be provided explicitly.
  subroutine atl_varSys_getStateForElement(fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
    ! -------------------------------------------------------------------- !

    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in)     :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in)                   :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)       :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in)      :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in)                   :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in)                   :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out)            :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDof, iComp, nScalars, level
    integer(kind=long_k) :: treeID
    type(atl_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, fPtr )

    ! number of scalars in state array
    nScalars = varSys%nScalars

    res = 0.0_rk
    do iElem = 1, nElems
      treeID = fptr%solverData%levelPointer(elemPos(iElem))
      level = tem_levelOf(tree%treeID(elemPos(iElem)))
      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res( ( iElem - 1 ) * fun%nComponents * nDofs &
            &   + ( iDof - 1 ) * fun%nComponents       &
            &   + iComp )                              &
            & = fPtr%solverData%statedata_listPtr(level)%state(   &
            &   treeID, iDof, fun%state_varPos(iComp)  )
        end do !iComp
      end do !iDof
    end do !iElem

  end subroutine atl_varSys_getStateForElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Interface description for a variable access method (single point).
  !!
  !! To obtain values of a given variable, it is necessary to state the
  !! space and time at which the variable should be evaluated.
  !! The interface is vectorized and provides n values at n different space
  !! locations at the same time.
  !! Of course the variable system itself also needs to be passed in, to
  !! allow the computation of other derived quantities as needed.
  !! The method description itself is passed in automatically, and has not
  !! to be provided explicitly.
  subroutine atl_varSys_getStateForPoint( fun, varsys, point, time, tree, &
      &                                   nPnts, res)
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
    ! -------------------------------------------------------------------- !
    integer :: iPoint, nComp, iDof, nDofs, level, iComp, maxPolyDegree, elemPos
    integer :: coord(4), levelPos
    integer(kind=long_k):: treeID
    integer :: basisType
    integer :: iDir
    real(kind=rk) :: bary(3)
    real(kind=rk) :: extent
    real(kind=rk) :: halfwidth
    real(kind=rk) :: MappedCoord(1,3)
    real(kind=rk), allocatable :: polyval(:,:), val(:)
    type(atl_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, fPtr )

    allocate( val( fun%nComponents ))

    nComp = fun%nComponents

    ! step 1 : Loop over total number of Points
    do iPoint = 1, nPnts
      val = 0.0_rk

      ! Step 2 : Get the treeID of the element and then get the element number
      coord =  tem_CoordOfReal(tree, point(iPoint,:), tree%global%maxLevel )
      treeId = tem_IdOfCoord(coord)
      ! get the position of treeid or position of the parent treeid
      elemPos = abs(tem_PosofId(treeId, tree%treeID))

      if( elemPos == 0 ) then
        call tem_abort( 'Element in multilevel mesh not on the finest level.' &
          & // ' Fix the position determination here.'                        )
      endif

      ! start element wise timer for LB weights
      call tem_startTimer( me          = atl_cpl_elemTimers, &
        &                  timerHandle = elemPos             )

      levelPos = fPtr%solverData%levelPointer(elemPos)
      level = tem_levelOf( tree%treeID( elemPos ) )

      bary = tem_BaryOfId(tree,TreeID)
      extent = tem_ElemSize( tree, TreeID)
      halfwidth = extent/2


      select case(fPtr%solverData%equationPtr%nDimensions)
      case(1)
        maxPolyDegree = fptr%solverData%scheme_listPtr(level)%modg_1d     &
          &                                                  %maxpolydegree
        basisType     = fPtr%solverdata            &
          &                 %scheme_listPtr(level) &
          &                 %modg_1d               &
          &                 %basisType
      case(2)
        maxPolyDegree = fptr%solverData%scheme_listPtr(level)%modg_2d     &
          &                                                  %maxpolydegree
        basisType     = fPtr%solverdata            &
          &                 %scheme_listPtr(level) &
          &                 %modg_2d               &
          &                 %basisType
      case(3)
        maxPolyDegree = fptr%solverData%scheme_listPtr(level)%modg        &
          &                                                  %maxpolydegree
        basisType     = fPtr%solverdata            &
          &                 %scheme_listPtr(level) &
          &                 %modg                  &
          &                 %basisType
      end select

      select case(basisType)
      case(Q_space)
        nDofs = (maxPolyDegree + 1)**fPtr%solverData%equationPtr%nDimensions
      case(P_space)
        nDofs = (maxPolyDegree + 1)
        do iDir = 2,fPtr%solverData%equationPtr%nDimensions
          nDofs = (nDofs * (maxPolyDegree + iDir)) / iDir
        end do
      end select

      allocate( polyval( nDofs, 1))
      polyVal = 0.0_rk

      ! Step 3 : Get the Mapped Coordinate from the physical coordinates
      MappedCoord(1,1) = (point(iPoint,1)- bary(1))/halfwidth
      MappedCoord(1,2) = (point(iPoint,2)- bary(2))/halfwidth
      MappedCoord(1,3) = (point(iPoint,3)- bary(3))/halfwidth

      ! Step 4 : Get the nodal coefficient
      select case (fPtr%solverData%equationPtr%nDimensions)
      case(1)
        call atl_evalLegendreTensPoly1d(      &
          &    coords        = MappedCoord,   &
          &    ncoords       = 1,             &
          &    maxPolyDegree = maxpolyDegree, &
          &    basisType     = basisType,     &
          &    polyVal       = polyVal        )
      case(2)
        call atl_evalLegendreTensPoly2d(      &
          &    coords        = MappedCoord,   &
          &    ncoords       = 1,             &
          &    maxPolyDegree = maxpolyDegree, &
          &    basisType     = basisType,     &
          &    polyVal       = polyVal        )
      case(3)
        call ply_evalLegendreTensPoly(        &
          &    coords        = MappedCoord,   &
          &    ncoords       = 1,             &
          &    maxPolyDegree = maxpolyDegree, &
          &    basisType     = basisType,     &
          &    polyVal       = polyVal        )
      end select

      ! Step 5 : Multiply with state to get the exact  point value
      do iComp = 1, nComp
        do iDof = 1, nDofs
          val(icomp) = val(icomp)                     &
            & + polyVal(iDof,1)                       &
            &   * fPtr%solverdata                     &
            &         %statedata_listPtr(level)       &
            &         %state( levelPos,               &
            &                 iDof,                   &
            &                 fun%state_varPos(iComp) )
        end do
      end do

      ! Step 6 : Store it in res
      res((iPoint-1)*nComp+1 : iPoint*nComp) = val(:)

      deallocate(polyval)

      call tem_stopTimer( me          = atl_cpl_elemTimers, &
        &                 timerHandle = elemPos             )
    end do ! loop over points

  end subroutine atl_varSys_getStateForPoint
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine takes points coordinates, stores them in the method_data and
  !! return indices where points are located in the growing array of points or
  !! values ( sometimes we do not need to store the points )
  !! It is need to setup points for every variable. Points will be provided by
  !! boundaries or sources depends on what uses the variable. This points do not
  !! change with time . This indices will be stored in corresponding boundary
  !! or source to evaluate a value on that point later using
  !! tem_varSys_proc_getValOfIndex.
  subroutine atl_varSys_setupStateIndices( fun, varSys, point, offset_bit, &
      &                                    iLevel, tree, nPnts, idx        )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, for this routine
    !! we need the location where to store the points.
    class(tem_varSys_op_type), intent(in) :: fun
    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys
    !>  arrays of points for which the indices are returned
    real(kind=rk), intent(in) :: point(:,:)
    !> Offset bit encoded as character for every point.
    !! If not present default is to center i.e offset_bit = achar(1+4+16)
    character, optional, intent(in) :: offset_bit(:)
    !> the point data need to be loaded levelwise, we need the current level
    integer, intent(in) :: iLevel
    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree
    !> number of points
    integer, intent(in) :: nPnts
    integer, intent(out) :: idx(:)
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    integer :: coord(4)
    integer(kind=long_k) :: treeID
    integer :: elemPos
    integer :: loc_level
    real(kind=rk) :: extent
    real(kind=rk) :: halfwidth
    real(kind=rk) :: MappedCoord(3)
    real(kind=rk) :: bary(3)
    type(atl_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !

    write(logUnit(10),*) 'setup indices for the points of state variable ', &
      &                 trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr )

    ! transform the point coordinates into coordintes on the reference element
    do iPoint = 1, nPnts

      ! Get the treeID of the element and then get the element number, use the
      ! maxLevel of the tree to get the coordinates of the elements.
      coord =  tem_CoordOfReal( tree, point(iPoint,:), tree%global%maxLevel )

      treeId = tem_IdOfCoord(coord)
      elemPos = tem_PosofId(treeId, tree%treeID)
      loc_level = tem_levelOf( tree%treeID( elemPos ) )
      !based on the local level get the correct treeid required for barr, extend
      coord =  tem_CoordOfReal( tree, point(iPoint,:), loc_level)
      treeId = tem_IdOfCoord(coord)
      ! append grw array for this point with element position
      call append(fPtr%pointData%pntLvl(iLevel)%elemPos, elemPos)

      ! Get the Mapped Coordinate from the physical coordinates
      bary = tem_BaryOfId(tree,TreeID)
      extent = tem_ElemSize(tree, TreeID)
      halfwidth = extent/2

      MappedCoord(1) = ( point(iPoint,1) - bary(1) ) / halfwidth
      MappedCoord(2) = ( point(iPoint,2) - bary(2) ) / halfwidth
      MappedCoord(3) = ( point(iPoint,3) - bary(3) ) / halfwidth

      ! Append the growing array for points
      call append(fPtr%pointData%pntLvl(iLevel)%grwPnt, MappedCoord)
      ! Append the growing array for local level of each point
      call append(fPtr%pointData%pntLvl(iLevel)%pntLevel, loc_level)

      ! store the index,
      idx(iPoint) = fPtr%pointData%pntLvl(iLevel)%elemPos%nVals
    end do

    call truncate(fPtr%pointData%pntLvl(iLevel)%grwPnt)
    call truncate(fPtr%pointData%pntLvl(iLevel)%elemPos)
    call truncate(fPtr%pointData%pntLvl(iLevel)%pntLevel)

  end subroutine atl_varSys_setupStateIndices
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine for gettint the actual value for a given array of indices.
  !! The indices belong to the grwarray of points storing levelwise in
  !! Pointdata%pntLvl(iLevel).
  !! Hence this routines takes the indeices as input, can refer to the pointData
  !! and evaluate the variable and returns the values
  subroutine atl_varSys_getStateValOfIndex( fun, varSys, time, iLevel, &
    &                                       idx, idxLen, nVals,  res   )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: coord(:,:)
    real(kind=rk), allocatable :: polyval(:,:), val(:)
    integer :: maxpolyDegree, nDofs, iDof, nComp, iComp, iPoint, ii, iElem
    integer :: elemPos, coordPos
    type(atl_varSys_data_type), pointer :: fPtr
    type(dyn_intArray_type) :: unique_elemPos
    type(grw_intArray_type) :: grw_nCoordPerElem
    logical :: wasAdded
    integer :: pos, global_count, buf_start, buf_end, loc_count
    integer :: loc_level, levelPos
    integer, allocatable :: pnt_pos(:)
    ! -------------------------------------------------------------------- !
    write(logUnit(4),*) 'Get the values of indices for state variable ',  &
      &                  trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, fPtr )


    allocate(coord(nVals,3))
    allocate(val(fun%nComponents))

    ! distinguish if we have an array of index or we have contingous memory
    ! acces where index are always first entries!
    if ( .not. present(idxLen) ) then
      coord(:,1)=fPtr%pointData%pntLvl(iLevel)%grwPnt%coordX%val(idx)
      coord(:,2)=fPtr%pointData%pntLvl(iLevel)%grwPnt%coordY%val(idx)
      coord(:,3)=fPtr%pointData%pntLvl(iLevel)%grwPnt%coordZ%val(idx)
    else ! idxLen is present
      ! integer which stores the position in the coord array
      coordPos = 1
      do ii = 1, size(idxLen)
        coord(coordPos:coordPos+idxLen(ii)-1,1) =  &
          & fPtr%pointData                         &
          &     %pntLvl(iLevel)                    &
          &     %grwPnt                            &
          &     %coordX                            &
          &     %val( idx(ii):idx(ii)+idxLen(ii)-1 )
        coord(coordPos:coordPos+idxLen(ii)-1,2) =  &
          & fPtr%pointData                         &
          &     %pntLvl(iLevel)                    &
          &     %grwPnt                            &
          &     %coordY                            &
          &     %val( idx(ii):idx(ii)+idxLen(ii)-1 )
        coord(coordPos:coordPos+idxLen(ii)-1,3) =  &
          & fPtr%pointData                         &
          &     %pntLvl(iLevel)                    &
          &     %grwPnt                            &
          &     %coordZ                            &
          &     %val( idx(ii):idx(ii)+idxLen(ii)-1 )
        coordPos = coordPos+idxLen(ii)
      end do
    end if

    ! sort the point array elementwise
    ! loop over points and store them per element
    call init( unique_elemPos )
    call init( grw_nCoordPerElem )
    do iPoint = 1, nVals
      elemPos = fPtr%pointData%pntLvl(iLevel)%elemPos%val(idx(iPoint))
      ! get an unqiue array of elempos
      call append( me       = unique_elemPos, &
        &          val      = elemPos,        &
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
        if (unique_elemPos%val(iElem) == fPtr%pointData%pntLvl(iLevel)%elemPos%val(idx(iPoint)) ) then
          pnt_pos(global_count) = iPoint
          global_count = global_count +1
        end if
      end do
    end do

    if ( any(pnt_pos == -1) ) then
        write(*,*) 'problem in getvalOfIndex: -1 in pnt_pos occured!!!'
        call tem_abort()
    end if

    ! now we go the sequence of the element
    do iElem = 1, unique_elemPos%nVals
      ! start element wise timer for LB weights
      call tem_startTimer( me          = atl_cpl_elemTimers,        &
        &                  timerHandle = unique_elemPos%val(iElem)  )

      ! upper and lower bound of array access according to elemPos offset/
      ! counter grw_nCoordPerElem
      buf_start = sum(grw_nCoordPerElem%val(1:iElem-1))+1
      buf_end = sum(grw_nCoordPerElem%val(1:iElem))

      ! get the correct levelPos
      levelPos = fPtr%solverData%levelPointer(unique_elemPos%val(iElem))

      ! get teh correct local level for the points
      ! when coupling with apesmate, it could happen that the requested
      ! level is different from local level
      ! since elementwise access it is sufficent to ask the first point
      ! for local level
      loc_level = fPtr%pointData%pntLvl(iLevel)%pntLevel%val(idx(pnt_pos(buf_start)))

      ! Coords are the coordinates on the reference element, hence we can get
      ! directly the nodal coefficients
      select case (fPtr%solverData%equationPtr%nDimensions)
      case(1)
        maxPolyDegree = fptr%solverData%scheme_listPtr(loc_level)%modg_1d     &
          &                                                      %maxpolydegree
        nDofs = (maxPolyDegree + 1)
        allocate( polyval( nDofs, grw_nCoordPerElem%val(iElem)))
        call atl_evalLegendreTensPoly1d(                          &
          & coords        = Coord(pnt_pos(buf_start:buf_end),:),  &
          & ncoords       = grw_nCoordPerElem%val(iElem),         &
          & maxPolyDegree = maxpolyDegree,                        &
          & basisType     = fPtr%solverdata                       &
          &                     %scheme_listPtr(loc_level)        &
          &                     %modg_1d                          &
          &                     %basisType,                       &
          & polyVal       =  polyVal                              )
      case(2)
        maxPolyDegree = fptr%solverData%scheme_listPtr(loc_level)%modg_2d     &
          &                                                     %maxpolydegree
        nDofs = (maxPolyDegree + 1)**2
        allocate( polyval( nDofs, grw_nCoordPerElem%val(iElem)))
        call atl_evalLegendreTensPoly2d(                          &
          & coords        = Coord(pnt_pos(buf_start:buf_end),:),  &
          & ncoords       = grw_nCoordPerElem%val(iElem),         &
          & maxPolyDegree = maxpolyDegree,                        &
          & basisType     = fPtr%solverdata                       &
          &                     %scheme_listPtr(loc_level)        &
          &                     %modg_2d                          &
          &                     %basisType,                       &
          & polyVal       =  polyVal                              )
      case(3)
        maxPolyDegree = fptr%solverData%scheme_listPtr(loc_level)%modg%maxpolydegree
        nDofs = (maxPolyDegree + 1)**3
        allocate( polyval( nDofs, grw_nCoordPerElem%val(iElem)))
        call ply_evalLegendreTensPoly(                            &
          & coords        = Coord(pnt_pos(buf_start:buf_end),:),  &
          & ncoords       = grw_nCoordPerElem%val(iElem),         &
          & maxPolyDegree = maxpolyDegree,                        &
          & basisType     = fPtr%solverData                       &
          &                     %scheme_listPtr(loc_level)        &
          &                     %modg                             &
          &                     %basisType,                       &
          & polyVal       =  polyVal                              )
      end select

      ! Multiply with state to get the exact  point value
      ! and copi it to the output array
      nComp = fun%nComponents
      if (present(idxlen)) then
        !! @todo Implement blocks of contiguous indices.
        call tem_abort(                                                   &
          & 'atl_varSys_getStateValOfIndex is not implemented for idxlen' )
      else
        loc_count = 0
        do iPoint = buf_start, buf_end
          loc_count = loc_count + 1
          val = 0.0_rk
          ! now we go in sequence of the element
          elemPos = fPtr%pointData%pntLvl(iLevel)%elemPos%val(idx(pnt_Pos(iPoint)))
          do iComp = 1, nComp
            do iDof = 1, nDofs
              val(icomp) = val(icomp)                                           &
                &        + polyVal(iDof,loc_count)                              &
                &          * fPtr%solverdata                                    &
                &                %statedata_listPtr(loc_level)                  &
                &                %state( levelPos, iDof, fun%state_varPos(iComp) )
            end do
          end do
          res((pnt_pos(iPoint)-1)*nComp+1: (pnt_pos(iPoint)-1)*nComp+nComp)  = val(:)
        end do!iPoint
      end if
      deallocate(polyVal)
      call tem_stopTimer( me          = atl_cpl_elemTimers,        &
        &                 timerHandle =unique_elemPos%val(iElem)   )
    end do!ielem

    deallocate(coord)
    deallocate(val)
    call destroy( unique_elemPos )
    call destroy( grw_nCoordPerElem )

  end subroutine atl_varSys_getStateValofIndex
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine creates hard coded variables which are required by ateles.
  !!
  !! For example: Lot of boundary condition variables require constant zero as
  !! spacetime-function, thus we create a variable "zero_const" here to be
  !! able to refer to it in boundary condition definitions.
  subroutine atl_create_fortranVar( me )
    ! -------------------------------------------------------------------- !
    !> Contains list of hard coded variables
    type(tem_variable_type), allocatable, intent(out) :: me(:)
    ! -------------------------------------------------------------------- !
    allocate( me(1) )
    ! create variable which return constant zero
    me(1)%label = 'zero_const'
    me(1)%nComponents = 1
    me(1)%varType = 'st_fun'
    me(1)%evalType = 'add'
    allocate( me(1)%st_fun(1) )
    allocate( me(1)%st_fun(1)%const(1) )
    me(1)%st_fun(1)%fun_kind = 'const'
    me(1)%st_fun(1)%const = 0.0_rk
    allocate( me(1)%st_fun(1)%geom(1) )
    me(1)%st_fun(1)%geom(1)%kind = 'all'
    me(1)%st_fun(1)%geom(1)%shapeID = tem_global_shape

  end subroutine atl_create_fortranVar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to obtain a modal representation for a variable, which is only
  !! available in nodal form, like space-time functions.
  subroutine atl_generic_fromNodal_getElement( fun, varsys, elempos, time, &
    &                                          tree, nElems, nDofs, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of element in tree%treeID to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of elements to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (nComponents of resulting variable) x (nDegrees of freedom) x (nElems)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(tem_st_fun_listElem_type), pointer :: md_ptr
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk), allocatable :: modalRes(:,:)
    real(kind=rk), allocatable :: nodalRes(:,:), res2d(:,:)
    real(kind=rk), allocatable :: physPoints(:,:)
    integer :: level
    integer :: pos
    integer :: iElem, iDof, iComp, iPoint
    integer :: nPoints
    ! -------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, md_Ptr )
    call C_F_POINTER( md_Ptr%solver_bundle, fPtr )

    allocate(res2d(nDofs,fun%nComponents))

    do iElem = 1, nElems
      level = tem_levelOf(tree%treeID(elemPos(iElem)))
      pos = fptr%solverData%poly_proj_posPtr(level)

      ! Obtain the point coordinates to get values for
      select case(fPtr%solverData%equationPtr%nDimensions)
      case(1)
        nPoints = fPtr%solverData%polyProjectPtr(pos)%body_1d%nQuadPoints
        ! Reallocate the physical points if they are different in size from the
        ! previous ones.
        if (.not. allocated(physPoints)) then
          allocate(physPoints(nPoints,3))
        else
          if (size(physPoints,1) /= nPoints) then
            deallocate(physPoints)
            allocate(physPoints(nPoints, 3))
          end if
        end if

        ! Get the physical points to use in the get_point call
        call atl_ref_in_elempos( refPoints  = fPtr%solverData          &
          &                                       %polyProjectPtr(pos) &
          &                                       %body_1d             &
          &                                       %nodes,              &
          &                      tree       = tree,                    &
          &                      elemPos    = elemPos(iElem),          &
          &                      physPoints = physPoints               )

      case(2)
        nPoints = fPtr%solverData%polyProjectPtr(pos)%body_2d%nQuadPoints
        ! Reallocate the physical points if they are different in size from the
        ! previous ones.
        if (.not. allocated(physPoints)) then
          allocate(physPoints(nPoints,3))
        else
          if (size(physPoints,1) /= nPoints) then
            deallocate(physPoints)
            allocate(physPoints(nPoints, 3))
          end if
        end if

        ! Get the physical points to use in the get_point call
        call atl_ref_in_elempos( refPoints  = fPtr%solverData          &
          &                                       %polyProjectPtr(pos) &
          &                                       %body_2d             &
          &                                       %nodes,              &
          &                      tree       = tree,                    &
          &                      elemPos    = elemPos(iElem),          &
          &                      physPoints = physPoints               )

      case(3)
        nPoints = fPtr%solverData%polyProjectPtr(pos)%body_3d%nQuadPoints
        if (.not. allocated(physPoints)) then
          allocate(physPoints(nPoints,3))
        else
          if (size(physPoints,1) /= nPoints) then
            deallocate(physPoints)
            allocate(physPoints(nPoints, 3))
          end if
        end if

        ! Get the physical points to use in the get_point call
        call atl_ref_in_elempos( refPoints  = fPtr%solverData          &
          &                                       %polyProjectPtr(pos) &
          &                                       %body_3d             &
          &                                       %nodes,              &
          &                      tree       = tree,                    &
          &                      elemPos    = elemPos(iElem),          &
          &                      physPoints = physPoints               )
      end select

      !!@todo Compute nodalRes by using getpoint.
      if (.not. allocated(nodalRes)) then
        allocate(nodalRes(nPoints,fun%nComponents))
      else
        if (size(nodalRes,1) /= nPoints) then
          deallocate(nodalRes)
          allocate(nodalRes(nPoints, fun%nComponents))
        end if
      end if
      if (.not. allocated(modalRes)) then
        allocate(modalRes(nPoints,fun%nComponents))
      else
        if (size(modalRes,1) /= nPoints) then
          deallocate(modalRes)
          allocate(modalRes(nPoints, fun%nComponents))
        end if
      end if

      do iPoint=1,nPoints
        call fun%get_point( varsys = varsys,                      &
          &                 point  = physPoints(iPoint:iPoint,:), &
          &                 time   = time,                        &
          &                 tree   = tree,                        &
          &                 nPnts  = 1,                           &
          &                 res    = nodalRes(iPoint,:)           )
      end do

      call ply_poly_project_n2m( &
        & me         = fPtr%solverData%polyProjectPtr(pos),     &
        & dim        = fPtr%solverData%equationPtr%nDimensions, &
        & nVars      = fun%nComponents,                         &
        & nodal_data = nodalRes,                                &
        & modal_data = modalRes                                 )

      call ply_convertFromOversample(                            &
        & poly_proj   = fPtr%solverData%polyProjectPtr(pos),     &
        & modalCoeffs = modalRes,                                &
        & nDim        = fPtr%solverData%equationPtr%nDimensions, &
        & state       = res2d                                    )

      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res(( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp) &
            & = res2d(iDof, iComp)
        end do
      end do

    end do

  end subroutine atl_generic_fromNodal_getElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to obtain a modal representation for a variable, which is only
  !! available in nodal form, like space-time functions.
  subroutine atl_generic_fromConst_getElement( fun, varsys, elempos, time, &
    &                                          tree, nElems, nDofs, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of element in tree%treeID to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of elements to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (nComponents of resulting variable) x (nDegrees of freedom) x (nElems)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(tem_st_fun_listElem_type), pointer :: md_ptr
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: nodalRes(fun%nComponents)
    real(kind=rk) :: physPoints(1,3)
    integer :: iElem, iComp
    ! -------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, md_Ptr )
    call C_F_POINTER( md_Ptr%solver_bundle, fPtr )

    ! Set all (higher) modes to 0.
    res = 0.0_rk

    do iElem = 1, nElems
      physPoints(1,:) = tem_baryofID(tree, tree%treeID(elempos(iElem)))

      ! Get the value at the barycenter to account for possibly overlapping
      ! stfun definitions and use the correct overlapping approach.
      call fun%get_point(         &
        &    varsys = varsys,     &
        &    point  = physPoints, &
        &    time   = time,       &
        &    tree   = tree,       &
        &    nPnts  = 1,          &
        &    res    = nodalRes    )

      ! Put values from barycenters into first mode.
      do iComp = 1, fun%nComponents
        res(( ielem-1)* fun%ncomponents* ndofs+( 1-1)* fun%ncomponents+icomp) &
          &  = nodalres(iComp)
      end do

    end do

  end subroutine atl_generic_fromConst_getElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Method to load user defined variables from the configuration.
  !!
  !! This adds the possibility to load ateles specific variable types.
  subroutine atl_varSys_load_user(L, parent, specifics, appender, iError)
    !> Lua script to load the variable data from.
    type(flu_State) :: L

    !> Parent table in the Lua script to read the variable from.
    integer, intent(in) :: parent

    !> Data to read for the solver specific variable
    type(c_ptr), intent(out) :: specifics

    !> Function pointer to use for appending the solver variable.
    procedure(tem_append_solverVar_method), pointer :: appender

    !> Indication whether the attempted reading was successful.
    integer, intent(out) :: iError
    ! -------------------------------------------------------------------- !
    type(atl_varSys_solverVar_type), pointer :: user_var
    ! -------------------------------------------------------------------- !

    appender => NULL()
    specifics = c_null_ptr
    allocate(user_var)

    call aot_get_val( L       = L,                &
      &               thandle = parent,           &
      &               val     = user_var%varType, &
      &               ErrCode = iError,           &
      &               default = 'none',           &
      &               key     = 'vartype'         )

    if (iError==0) then
      select case (trim(user_var%varType))
      case('legpoly')
        call atl_legpolyvar_load( L          = L,                  &
          &                       parent     = parent,             &
          &                       legpolyvar = user_var%legpolyvar )
        appender => atl_legpolyvar_append
        specifics = c_loc(user_var)
      case default
        deallocate(user_var)
        iError = 1
      end select
    end if

  end subroutine atl_varSys_load_user

end module atl_varSys_module
