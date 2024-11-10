! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> Some variables provided by Navier-Stokes system.
module atl_eqn_nvrstk_var_module
  use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr
  use env_module, only: rk, long_k, labelLen
  use treelmesh_module, only: treelmesh_type
  use tem_aux_module, only: tem_abort
  use tem_time_module, only: tem_time_type
  use tem_logging_module, only: logUnit
  use tem_varSys_module, only: tem_varSys_type,              &
    &                          tem_varSys_op_type,           &
    &                          tem_varSys_append_derVar,     &
    &                          tem_varSys_proc_point,        &
    &                          tem_varSys_proc_element,      &
    &                          tem_varSys_proc_setparams,    &
    &                          tem_varSys_proc_getparams,    &
    &                          tem_varSys_proc_setupIndices, &
    &                          tem_varSys_proc_getValOfIndex
  use tem_topology_module, only: tem_coordOfId, tem_IDofCoord, &
    &                            tem_levelOf
  use tem_geometry_module, only: tem_CoordOfReal, &
    &                            tem_PosofId
  use atl_varSys_module, only: atl_varSys_data_type, &
    &                          atl_varSys_solverData_type, &
    &                          atl_get_new_varSys_data_ptr
  use atl_operator_module, only: atl_opVar_setupIndices
  use atl_eqn_nvrstk_module, only: atl_navierstokes_type

  implicit none

  private

  public :: atl_append_nvrstk_derivedVars


contains

  !> Append / set methods and data to compute derived quantities to the
  !! variable system.
  !!
  !! Available quantities are:
  !!
  !! * inviscindicator: Indication whether the viscous terms are neglected in
  !!                    the element. It will be one where the viscous fluxes
  !!                    are computed and 0, where they'll be neglected.
  !! * shearestimate:   Basis for the inviscindicator, provides the estimate
  !!                    for the maximal shear in the element. This is based
  !!                    on the derivative estimation, which is only computed
  !!                    if adaptivity for viscous terms is activated!
  !! * all quantities from the euler system, see
  !!                   [[atl_eqn_euler_var_module:atl_append_euler_derivedVars]]
  subroutine atl_append_nvrstk_derivedVars(varSys, solverData)
    !> The variable system to modify. It has to contain the conservative
    !! and primitive variables already.
    type(tem_varSys_type), intent(inout) :: varSys
    !> the pointer to the data required for the varsys to fulfill all operations
    !! and derivations on the variables
    type(atl_varSys_solverData_type), target :: solverData
    ! -------------------------------------------------------------------- !
    logical :: wasAdded
    integer :: nComponents
    character(len=labelLen), allocatable :: invar_name(:)
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer &
      & :: get_valOfIndex => NULL()
    type(c_ptr) :: method_data
    ! -------------------------------------------------------------------- !

    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    ! Appending the inviscid indicator to the variable system.
    get_point => atl_viscindicator_getPoint
    get_element => atl_viscindicator_getElement
    setup_indices => atl_opVar_setupIndices
    method_data = atl_get_new_varSys_data_ptr(solverData)
    nComponents = 1

    allocate(invar_name(0))
    call tem_varSys_append_derVar( me             = varSys,            &
      &                            varName        = 'inviscindicator', &
      &                            nComponents    = nComponents,       &
      &                            input_varname  = invar_name,        &
      &                            method_data    = method_data,       &
      &                            get_point      = get_point,         &
      &                            get_element    = get_element,       &
      &                            set_params     = set_params,        &
      &                            get_params     = get_params,        &
      &                            setup_indices  = setup_indices,     &
      &                            get_valOfIndex = get_valOfIndex,    &
      &                            wasAdded       = wasAdded           )

    if (wasAdded) then
      write(logUnit(10),*) 'Appended variable: inviscindicator'
    else
      call tem_abort( 'Error: variable inviscindicator' &
        & // ' is not added to variable system'         )
    end if

    ! Appending the shear estimation as variable.
    get_point => atl_shearestimate_getPoint
    get_element => atl_shearestimate_getElement
    setup_indices => atl_opVar_setupIndices
    method_data = atl_get_new_varSys_data_ptr(solverData)
    nComponents = 1

    call tem_varSys_append_derVar( me             = varSys,          &
      &                            varName        = 'shearestimate', &
      &                            nComponents    = nComponents,     &
      &                            input_varname  = invar_name,      &
      &                            method_data    = method_data,     &
      &                            get_point      = get_point,       &
      &                            get_element    = get_element,     &
      &                            set_params     = set_params,      &
      &                            get_params     = get_params,      &
      &                            setup_indices  = setup_indices,   &
      &                            get_valOfIndex = get_valOfIndex,  &
      &                            wasAdded       = wasAdded         )
    if (wasAdded) then
      write(logUnit(10),*) 'Appended variable: shearestimate'
    else
      call tem_abort( 'Error: variable shearestimate' &
        & // ' is not added to variable system'       )
    end if
    deallocate(invar_name)

  end subroutine atl_append_nvrstk_derivedVars



  ! ------------------------------------------------------------------------ !
  !> Estimate the magnitude of shear terms.
  !!
  !! We use an estimate for mu*(dv/dx) from the deviations and derivative
  !! estimates of the conservative variables:
  !! dm/dx = dv/dx * rho + drho/dx * v
  !! dv/dx = (dm/dx - drho/dx * v)/rho
  !! dv/dx < (max(dm/dx) - max(drho/dx) * max(v))/min(rho)
  !! dv/dx < (max(dm/dx) * min(rho) - max(drho/dx) * max(m) ) / min(rho)**2
  pure function shear_estimate_3d(nvrstk, mean, deviation, grad) &
    &   result(maxshear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_navierstokes_type), intent(in) :: nvrstk

    !> The mean of each state variable.
    real(kind=rk), intent(in) :: mean(:)

    !> Estimation of maximal deviation of each state.
    real(kind=rk), intent(in) :: deviation(:)

    !> Estimation of maximal gradient of each state.
    real(kind=rk), intent(in) :: grad(:)

    !> Resulting estimate for maximal shear.
    real(kind=rk) :: maxshear
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rho_min
    real(kind=rk) :: m_max
    real(kind=rk) :: grad_mag
    ! -------------------------------------------------------------------- !

    rho_min = max(mean(1) - deviation(1), epsilon(rho_min))
    m_max = sqrt(  (abs(mean(2)) + deviation(2))**2 &
      &          + (abs(mean(3)) + deviation(3))**2 &
      &          + (abs(mean(4)) + deviation(4))**2 )
    grad_mag = sqrt(grad(2)**2 + grad(3)**2 + grad(4)**2)
    maxshear = nvrstk%mu / rho_min**2 &
      &        * (grad_mag * rho_min + grad(1) * m_max)

  end function shear_estimate_3d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Estimate the magnitude of shear terms.
  !!
  !! We use an estimate for mu*(dv/dx) from the deviations and derivative
  !! estimates of the conservative variables:
  !! dm/dx = dv/dx * rho + drho/dx * v
  !! dv/dx = (dm/dx - drho/dx * v)/rho
  !! dv/dx < (max(dm/dx) - max(drho/dx) * max(v))/min(rho)
  !! dv/dx < (max(dm/dx) * min(rho) - max(drho/dx) * max(m) ) / min(rho)**2
  pure function shear_estimate_2d(nvrstk, mean, deviation, grad) &
    &   result(maxshear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_navierstokes_type), intent(in) :: nvrstk

    !> The mean of each state variable.
    real(kind=rk), intent(in) :: mean(:)

    !> Estimation of maximal deviation of each state.
    real(kind=rk), intent(in) :: deviation(:)

    !> Estimation of maximal gradient of each state.
    real(kind=rk), intent(in) :: grad(:)

    !> Resulting estimate for maximal shear.
    real(kind=rk) :: maxshear
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rho_min
    real(kind=rk) :: m_max
    real(kind=rk) :: grad_mag
    ! -------------------------------------------------------------------- !

    rho_min = max(mean(1) - deviation(1), epsilon(rho_min))
    m_max = sqrt(  (abs(mean(2)) + deviation(2))**2 &
      &          + (abs(mean(3)) + deviation(3))**2 )
    grad_mag = sqrt(grad(2)**2 + grad(3)**2)
    maxshear = nvrstk%mu / rho_min**2 &
      &        * (grad_mag * rho_min + grad(1) * m_max)

  end function shear_estimate_2d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  subroutine atl_viscindicator_getPoint( fun, varsys, point, time,tree, nPnts, &
      &                                 res                                   )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer(kind=long_k) :: treeId
    integer :: iPnt
    integer :: coord(4)
    integer :: elempos
    integer :: pos
    integer :: level
    logical :: isinviscid
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    do iPnt=1,nPnts
      coord =  tem_CoordOfReal(tree, point(iPnt,:), tree%global%maxLevel)
      treeId = tem_IdOfCoord(coord)
      ! get the position of treeid or position of the parent treeid
      elemPos = abs(tem_PosofId(treeId, tree%treeID))
      level = tem_levelOf( tree%treeID( elemPos ) )
      Pos = fPtr%solverData%levelPointer(elemPos)
      isinviscid = fPtr%solverdata%equationPtr%NavierStokes%inviscous(    &
        &          mean = fPtr%solverdata                                 &
        &                     %statedata_listPtr(level)                   &
        &                     %state(pos, 1, :),                          &
        &          deviation = fPtr%solverdata                            &
        &                          %kerneldata_listPtr(level)             &
        &                          %deviation(pos,:),                     &
        &          grad = fPtr%solverdata                                 &
        &                     %kerneldata_listPtr(level)                  &
        &                     %maxgrad(pos,:)                             )
      if (isinviscid) then
        res(iPnt) = 0.0_rk
      else
        res(iPnt) = 1.0_rk
      end if
    end do

  end subroutine atl_viscindicator_getPoint


  subroutine atl_viscindicator_getElement(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    logical :: isinviscid
    integer :: iElem
    integer :: level
    integer :: firstdof
    integer :: pos
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    res = 0.0_rk

    firstdof = 1
    do iElem=1,nElems
      level = tem_levelOf(tree%treeID(elempos(iElem)))
      pos = fptr%solverData%levelPointer(elemPos(iElem))

      isinviscid = fPtr%solverdata%equationPtr%NavierStokes%inviscous(    &
        &          mean = fPtr%solverdata                     &
        &                     %statedata_listPtr(level)       &
        &                     %state(pos, 1, :),              &
        &          deviation = fPtr%solverdata                &
        &                          %kerneldata_listPtr(level) &
        &                          %deviation(pos,:),         &
        &          grad = fPtr%solverdata                     &
        &                     %kerneldata_listPtr(level)      &
        &                     %maxgrad(pos,:)                 )
      if (isinviscid) then
        res(firstdof) = 0.0_rk
      else
        res(firstdof) = 1.0_rk
      end if

      firstdof = firstdof + nDofs
    end do

  end subroutine atl_viscindicator_getElement



  subroutine atl_shearestimate_getPoint( fun, varsys, point, time,tree, nPnts, &
      &                                  res                                   )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer(kind=long_k) :: treeId
    integer :: iPnt
    integer :: coord(4)
    integer :: elempos
    integer :: pos
    integer :: level
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    do iPnt=1,nPnts
      coord =  tem_CoordOfReal(tree, point(iPnt,:), tree%global%maxLevel)
      treeId = tem_IdOfCoord(coord)
      ! get the position of treeid or position of the parent treeid
      elemPos = abs(tem_PosofId(treeId, tree%treeID))
      level = tem_levelOf( tree%treeID( elemPos ) )
      Pos = fPtr%solverData%levelPointer(elemPos)
      if (fPtr%solverData%kerneldata_listPtr(level)%nDims == 2) then
        res(iPnt) = shear_estimate_2d(                          &
          &  nvrstk = fPtr%solverdata%equationPtr%NavierStokes, &
          &  mean = fPtr%solverdata                             &
          &             %statedata_listPtr(level)               &
          &             %state(pos, 1, :),                      &
          &  deviation = fPtr%solverdata                        &
          &                  %kerneldata_listPtr(level)         &
          &                  %deviation(pos,:),                 &
          &  grad = fPtr%solverdata                             &
          &             %kerneldata_listPtr(level)              &
          &             %maxgrad(pos,:)                         )
      else
        res(iPnt) = shear_estimate_3d(                          &
          &  nvrstk = fPtr%solverdata%equationPtr%NavierStokes, &
          &  mean = fPtr%solverdata                             &
          &             %statedata_listPtr(level)               &
          &             %state(pos, 1, :),                      &
          &  deviation = fPtr%solverdata                        &
          &                  %kerneldata_listPtr(level)         &
          &                  %deviation(pos,:),                 &
          &  grad = fPtr%solverdata                             &
          &             %kerneldata_listPtr(level)              &
          &             %maxgrad(pos,:)                         )
      end if
    end do

  end subroutine atl_shearestimate_getPoint


  subroutine atl_shearestimate_getElement(fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer :: iElem
    integer :: level
    integer :: firstdof
    integer :: pos
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    res = 0.0_rk

    firstdof = 1
    do iElem=1,nElems
      level = tem_levelOf(tree%treeID(elempos(iElem)))
      pos = fptr%solverData%levelPointer(elemPos(iElem))

      if (fPtr%solverData%kerneldata_listPtr(level)%nDims == 2) then
        res(firstdof) = shear_estimate_2d(                      &
          &  nvrstk = fPtr%solverdata%equationPtr%NavierStokes, &
          &  mean = fPtr%solverdata                             &
          &             %statedata_listPtr(level)               &
          &             %state(pos, 1, :),                      &
          &  deviation = fPtr%solverdata                        &
          &                  %kerneldata_listPtr(level)         &
          &                  %deviation(pos,:),                 &
          &  grad = fPtr%solverdata                             &
          &             %kerneldata_listPtr(level)              &
          &             %maxgrad(pos,:)                         )
      else
        res(firstdof) = shear_estimate_3d(                      &
          &  nvrstk = fPtr%solverdata%equationPtr%NavierStokes, &
          &  mean = fPtr%solverdata                             &
          &             %statedata_listPtr(level)               &
          &             %state(pos, 1, :),                      &
          &  deviation = fPtr%solverdata                        &
          &                  %kerneldata_listPtr(level)         &
          &                  %deviation(pos,:),                 &
          &  grad = fPtr%solverdata                             &
          &             %kerneldata_listPtr(level)              &
          &             %maxgrad(pos,:)                         )
      end if

      firstdof = firstdof + nDofs
    end do

  end subroutine atl_shearestimate_getElement

end module atl_eqn_nvrstk_var_module
