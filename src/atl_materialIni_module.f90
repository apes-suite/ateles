! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2018-2020 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
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
! Make sure loglvl is defined to an integer value.
! Usually this should be defined on the command line with -Dloglvl=
! to some value.

module atl_materialIni_module
  use, intrinsic :: iso_c_binding,    only: c_f_pointer
  ! Treelm modules
  use env_module,                     only: rk, labelLen
  use tem_tools_module,               only: tem_horizontalSpacer
  use tem_aux_module,                 only: tem_abort
  use treelmesh_module,               only: treelmesh_type
  use tem_element_module,             only: eT_fluid
  use tem_time_module,                only: tem_time_type
  use tem_spacetime_fun_module,       only: tem_spacetime_for,       &
    &                                       tem_st_fun_listElem_type
  use tem_topology_module,            only: tem_levelOf
  use tem_faceData_module,            only: tem_invFace_map, &
    &                                       tem_left,        &
    &                                       tem_right
  use tem_grow_array_module,          only: grw_intArray_type, &
    &                                       init,              &
    &                                       append,            &
    &                                       destroy
  use tem_dyn_array_module,           only: dyn_intArray_type, &
    &                                       init,              &
    &                                       append,            &
    &                                       destroy,           &
    &                                       PositionOfVal
  use tem_logging_module,             only: logUnit,   &
    &                                       llerror,   &
    &                                       llwarning, &
    &                                       llinfo,    &
    &                                       lldebug
  use tem_comm_module,                only: tem_commPattern_type
  use tem_comm_env_module,            only: tem_comm_env_type
  use tem_varSys_module,              only: tem_varSys_type,              &
    &                                       tem_varSys_op_type,           &
    &                                       tem_varSys_append_derVar,     &
    &                                       tem_varSys_proc_point,        &
    &                                       tem_varSys_proc_element,      &
    &                                       tem_varSys_proc_setParams,    &
    &                                       tem_varSys_proc_getParams,    &
    &                                       tem_varSys_proc_setupIndices, &
    &                                       tem_varSys_proc_getValOfIndex
  use tem_varMap_module,              only: tem_variable_loadMapping,           &
    &                                       tem_possible_variable_type
  use tem_stringKeyValuePair_module,  only: grw_stringKeyValuePairArray_type, &
    &                                       init
  use tem_timer_module,               only: tem_startTimer, &
    &                                       tem_stopTimer

  ! Aotus modules
  use aotus_module,                   only: flu_State, &
    &                                       aot_get_val
  use aot_table_module,               only: aot_table_open, &
    &                                       aot_table_close

  ! Ateles modules
  use atl_timer_module,               only: atl_timerHandles, atl_elemTimers
  use atl_equation_module,            only: atl_equations_type
  use atl_eqn_maxwelldivcorr_var_module, &
    &                                 only: atl_getMaxPropSpeedDivCor
  use atl_eqn_maxwell_var_module,     only: atl_getMaxPropSpeed
  use atl_boundary_module,            only: atl_level_boundary_type
  use atl_parallel_module,            only: atl_init_faceStateBuffer
  use atl_scheme_module,              only: atl_modg_scheme_prp,    &
    &                                       atl_modg_2d_scheme_prp, &
    &                                       atl_modg_1d_scheme_prp
  use atl_reference_element_module,   only: atl_refToPhysCoord
  use atl_cube_elem_module,           only: atl_cube_elem_type
  use atl_scheme_module,              only: atl_scheme_type
  use atl_varSys_module,              only: atl_varSys_solverData_type,  &
    &                                       atl_get_new_varSys_data_ptr
  use atl_materialPrp_module,         only: atl_mixedMat_prp,               &
    &                                       atl_pureConstMat_prp,           &
    &                                       atl_material_type,              &
    &                                       atl_spacetime_fun_pointer_type, &
    &                                       atl_ConstMatIdx,                &
    &                                       atl_VarMatIdx,                  &
    &                                       atl_init_material_type
  use atl_materialFun_module,         only: atl_materialFun_type, &
    &                                       atl_mode_reduction_type

  ! Polynomials modules
  use ply_poly_project_module,        only: ply_get_quadpoints_faces,    &
    &                                       ply_get_quadpoints_faces_2d, &
    &                                       ply_get_quadpoints_faces_1d, &
    &                                       ply_poly_project_type

  implicit none

  private

  public :: atl_init_materialParams
  public :: atl_update_materialParams
  ! Made public solely for unit test access.
  public ::                        &
    & atl_create_materialElemList, &
    & atl_append_newMaterialVars

contains

  ! ****************************************************************************
  !> This routine implements the getElement interface for material variables.
  !!
  !! Material variables basically refer to spacetime functions which are
  !! already present in the varSys. This spacetime function variable is the
  !! only input variable for the material variable. Thus the solely purpose of
  !! this routine is to forward any getElement request to the material variable
  !! to the correspondig spacetime function variable.
  subroutine atl_getMaterialForElement( fun, varsys, elempos, time, tree, &
    &                                   nElems, nDofs, res                )
    ! --------------------------------------------------------------------------

    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the spacetime function to use or the
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
    ! --------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elempos = elempos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = res                                        )

  end subroutine atl_getMaterialForElement
  ! ****************************************************************************


  ! ****************************************************************************
  !> This routine implements the getPoint interface for material variables.
  !!
  !! Material variables basically refer to spacetime functions which are
  !! already present in the varSys. This spacetime function variable is the
  !! only input variable for the material variable. Thus the solely purpose of
  !! this routine is to forward any getPoint request to the material variable
  !! to the correspondig spacetime function variable.
  subroutine atl_getMaterialForPoint(fun, varsys, point, time, tree, nPnts, res)
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the spacetime function to use or the
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

    call varSys%method%val(fun%input_varPos(1))%get_point( varSys = varSys, &
      &                                                    point  = point,  &
      &                                                    time   = time,   &
      &                                                    tree   = tree,   &
      &                                                    nPnts  = nPnts,  &
      &                                                    res    = res     )

  end subroutine atl_getMaterialForPoint
  ! ****************************************************************************


  ! ****************************************************************************
  !> Adds the configured material variables to the variable system.
  !!
  !! All variables in \ref variables are added to the varSys. Their getPoint
  !! and getElement routines are confgured to retrieve their data from the
  !! corresponding spacetime function that was already added to the varSys.
  subroutine atl_append_newMaterialVars( varSys, varSys_data, poss_matVars,    &
    &                                    materialFun, variables                )
    ! --------------------------------------------------------------------------
    !> The variable system to which the souce variables have to be added
    type(tem_varSys_type), intent(inout) :: varSys
    !> Data of the variable System
    type(atl_varSys_solverData_type), intent(inout), target :: varSys_data
    !! index of the eval-source_routine in eval_source.
    type(tem_possible_variable_type) :: poss_matVars
    type(atl_materialFun_type), intent(inout) :: materialFun
    type(grw_stringKeyValuePairArray_type), intent(inout) :: variables
    ! --------------------------------------------------------------------------
    integer :: addedAtPos, iMat, iVar, nComponents, varPos
    integer :: iStFun
    integer :: inputVarPos
    character(labellen) :: varname, inputVarName
    logical :: wasAdded
    type(tem_st_fun_listElem_type), pointer :: stFunList
    procedure(tem_varSys_proc_point), pointer :: get_point
    procedure(tem_varSys_proc_element), pointer :: get_element
    procedure(tem_varSys_proc_setParams), pointer :: set_params
    procedure(tem_varSys_proc_getParams), pointer :: get_params
    procedure(tem_varSys_proc_setupIndices), pointer ::setup_indices
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_ValOfIndex
    ! --------------------------------------------------------------------------

    get_point => atl_getMaterialForPoint
    get_element => atl_getMaterialForElement
    set_params => null()
    get_params => null()
    setup_indices => null()
    get_ValOfIndex => NULL()

    materialFun%nMat = 0
    allocate( materialFun%nScalars( poss_matVars%varname%nVals ) )
    allocate( materialFun%matvar_pos( poss_matVars%varname%nVals ) )
    allocate( materialFun%matParNames( poss_matVars%varname%nVals ) )

    materialFun%is_transient = .false.

    ! loop over all present material variables
    do iVar = 1, poss_matVars%varname%nVals
      varname = ''
      do iMat = 1, variables%nVals

        if ( trim(variables%val(iMat)%key) &
          &  == trim(poss_matVars%varname%val(iVar)) ) then


          ! Here we add the source variable specific spacetime function as the
          ! last input variable. The state variables in the first slots stay
          ! unchanged.
          inputVarName = variables%val(iMat)%value

          varPos = PositionOfVal( me  = varSys%varName,    &
            &                     val = trim(inputVarName) )

          nComponents = poss_matVars%nComponents%val(iVar)

          ! Check whether the input variable is able to provide as much
          ! components as the material expects
          if (nComponents /= varSys%method%val(varPos)%nComponents) then
            call tem_abort( "Source variable for material "          &
              & // trim(variables%val(iMat)%key)                     &
              & // " doesn't have the correct number of components." )
          endif

          ! get the name of the current variable and prefix it with mat_ as it
          ! is a material parameter. Without prefixing, it wouldn't be possible
          ! to define a user defined variable in the variable table that has the
          ! name of the material parameter. This handling of the variable names
          ! is equal to the name handling for source variables.
          varname = 'mat_' // trim(variables%val(iMat)%key)

          EXIT
        end if
      end do

      if (varname /= '') then
        call tem_varSys_append_derVar(                                 &
          & me             = varSys,                                   &
          & varName        = varname,                                  &
          & operType       = 'material',                               &
          & nComponents    = nComponents,                              &
          & input_varname  = (/ inputVarName /),                       &
          & method_data    = atl_get_new_varSys_data_ptr(varSys_data), &
          & get_point      = get_point,                                &
          & get_element    = get_element,                              &
          & set_params     = set_params,                               &
          & get_params     = get_params,                               &
          & setup_indices  = setup_indices,                            &
          & get_valOfIndex = get_valOfIndex,                           &
          & pos            = addedAtPos,                               &
          & wasAdded       = wasAdded                                  )

        if (wasAdded) then
          materialFun%nMat = materialFun%nMat + 1
          materialFun%nScalars( materialFun%nMat ) = nComponents
          materialFun%matvar_pos( materialFun%nMat ) = addedAtPos
          materialFun%matParNames( materialFun%nMat ) = varName
          ! We expect only one input variable, which is the variable containing
          ! the spacetime functions. Thus we hard-code index 1 here.
          inputVarPos = varSys%method%val(addedAtPos)%input_varPos(1)
          ! Now that we have the position of the lua variable, we can cast the
          ! list of spacetime functions for this variable, which is pointed to
          ! by the method data, back into a Fortran type to iterate over it.
          call c_f_pointer(varSys%method%val(inputVarPos)%method_data, &
            &              stFunList)
          do iStFun = 1, stFunList%nVals
            select case(trim(stFunList%val(iStFun)%fun_kind))
            case ('const')
              ! Nothing to do, keep the assumption, that the material is not
              ! transient.
            case ('combined')
              if (trim(stFunList%val(iStFun)%temporal%kind) /= 'const') then
                materialFun%is_transient = .true.
              end if
            case default
              materialFun%is_transient = .true.
            end select
          end do
        else

          write(logUnit(llerror),*) 'Error: variable ' // trim(varname) &
            & // ' was not added to the variable system.'

        end if

      else ! if (varname /= '')

        ! No definition for one of the materials found.
        ! FAIL now.
        !> @todo We might use default settings here instead, if applicable, as
        !!       for example for the conductivity in maxwell equations.
        !!       However, the behavior should be configured by the equation
        !!       system. Maybe this could also be handled when loading the
        !!       materials instead.
        write(logunit(llerror),*) 'FAILED to find a definition for ' &
          & // trim(poss_matVars%varname%val(iVar))
        write(logunit(llerror),*) 'You need to define a variable for this in' &
          & // ' the material table!'
        call tem_abort( 'Not knowing what to do, ABORTING...' )

      end if

    end do

  end subroutine atl_append_newMaterialVars
  ! ****************************************************************************


  ! ****************************************************************************
  !> Reassigns the spacetime function pointers based on the material variable
  !! position and the spacetime function position.
  subroutine atl_reassignStFunPtr( varSys, tree, material_list, iLevel, iDir, &
    &                              iMat                                       )
    ! --------------------------------------------------------------------------
    !> global variable system to which luaVar to be appended
    type(tem_varSys_type), intent(in) :: varSys
    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(in) :: tree
    !> The description of the material properties. The compute lists in the
    !! material description is filled up by calling this subroutine.
    type(atl_material_type), intent(inout) :: &
      & material_list(tree%global%minLevel:tree%global%maxLevel)
    integer, intent(in) :: iLevel
    integer, intent(in) :: iDir
    integer, intent(in) :: iMat
    ! --------------------------------------------------------------------------
    integer :: iFace, iLR, matVarPos, stFunPos, nFaces
    type(tem_st_fun_listElem_type), pointer :: stFunList
    ! --------------------------------------------------------------------------

    nFaces = size( material_list(iLevel)%material_desc       &
          &                             %material_face(iDir) &
          &                             %mat,                &
          &         1                                        )
    do iFace=1,nFaces
      do iLR=1,2
        matVarPos = material_list(iLevel)%material_desc          &
          &                              %material_face(iDir)    &
          &                              %mat(iFace, iLR, iMat)  &
          &                              %matVarPos
        ! There are situations where not all faces have material assigned. This
        ! is due to the fact that the array is oversized to also contain
        ! boundary faces, which are not necessarily initialized yet. Thus we
        ! have to check for valid matVarPos. We assume that stFunPos is set to a
        ! meaningful value when matVarPos is within the valid range.
        if ( (matVarPos > 0) .and. (matVarPos <= varSys%method%nVals) ) then
          stFunPos = material_list(iLevel)%material_desc         &
            &                             %material_face(iDir)   &
            &                             %mat(iFace, iLR, iMat) &
            &                             %stFunPos
          call c_f_pointer( cptr = varSys%method%val(matVarPos)%method_data, &
            &               fptr = stFunList                                 )
          material_list(iLevel)%material_desc         &
            &                  %material_face(iDir)   &
            &                  %mat(iFace, iLR, iMat) &
            &                  %stFunPtr              &
            & => stFunList%val(stFunPos)
        end if
      end do
    end do

  end subroutine atl_reassignStFunPtr
  ! ****************************************************************************

  ! ****************************************************************************
  !> Read the configuration for the material paramters for Maxwell equations
  !! from configuration files and init the material parameter datatype.
  subroutine atl_init_materialParams( equation, tree, varSys_data,           &
    &                                 material_list, mesh_list, scheme_list, &
    &                                 boundary_list, time, conf, proc,       &
    &                                 commPattern, poly_proj_list,           &
    &                                 levelpointer, initMaterial             )
    ! --------------------------------------------------------------------------
    !> Description on the equation system to solve.
    type(atl_Equations_type), intent(inout) :: equation
    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(in) :: tree
    !> Data in variable system
    type(atl_varSys_solverData_type), intent(inout) :: varSys_data
    !> The description of the material properties. The compute lists in the
    !! material description is filled up by calling this subroutine.
    type(atl_material_type), intent(inout) :: &
      & material_list(tree%global%minLevel:tree%global%maxLevel)
    !> Description of the mesh
    type(atl_cube_elem_type), intent(inout) :: &
      & mesh_list(tree%global%minLevel:tree%global%maxLevel)
    !> Information about the scheme
    type(atl_scheme_type), intent(in) :: &
      & scheme_list(tree%global%minLevel:tree%global%maxLevel)
    !> Global description of the boundary conditions.
    type(atl_level_boundary_type), intent(in) :: &
      & boundary_list(tree%global%minLevel:tree%global%maxLevel)
    !> The current simulation time
    type(tem_time_type), intent(in) :: time
    !> Flu binding to configuration file.
    type(flu_State), intent(inout) :: conf
    !> mpi communication environment with mpi communicator
    type(tem_comm_env_type), intent(inout) :: proc
    !> mpi communication pattern type
    type(tem_commPattern_type), intent(inout) :: commPattern
    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout)  :: poly_proj_list(:)
    !> Pointer from treeIDlist entry to level-wise fluid part of total list.
    !! The length of this vector is the total number of cubic elements.
    integer, intent(in) :: levelPointer(:)
    !> Information about expected material parameters and the user defined
    !! variables to the expected material parameters.
    type(atl_init_material_type) :: initMaterial
    ! --------------------------------------------------------------------------
    integer :: iMat, iLevel, iDir, iLR, iError, iProc
    integer :: nFaces, nElemsCons, nElemsVar, nquadpoints, nScalars, nBcs
    integer :: nConstFace, nVarFace
    integer :: tblHandle, matTblHandle
    integer :: projectionPos
    character(len=labelLen), allocatable :: materials(:)
    ! --------------------------------------------------------------------------

    ! initialize just to silence a compiler warning, will be overwritten later
    ! on
    nquadpoints = 0

    call tem_horizontalSpacer(fUnit=logUnit(llerror))
    call aot_table_open( L       = conf,      &
      &                  thandle = tblHandle, &
      &                  key     = "equation" )
    if( tblHandle > 0 ) then
      ! Initilized keyValuePair dict for material
      call init( me = initMaterial%materialDict )

      ! Read the enabled state for mode reduction from the equation.material
      ! table.
      call aot_table_open( L       = conf,         &
        &                  parent  = tblHandle,    &
        &                  thandle = matTblHandle, &
        &                  key     = 'material'    )
      call aot_get_val( L       = conf,                                     &
        &               thandle = matTblHandle,                             &
        &               key     = 'mode_reduction',                         &
        &               val     = equation%material%mode_reduction%enabled, &
        &               default = .false.,                                  &
        &               ErrCode = iError                                    )
      call aot_get_val( L       = conf,                                       &
        &               thandle = matTblHandle,                               &
        &               key     = 'reduction_threshold',                      &
        &               val     = equation%material%mode_reduction%threshold, &
        &               default = 1.0_rk,                                     &
        &               ErrCode = iError                                      )
      call aot_table_close( conf, matTblHandle )
      write(logUnit(1),'(A,L5)') 'Mode reduction for flux computation: ', &
        & equation%material%mode_reduction%enabled
      write(logUnit(1),'(A,F8.2)') 'Mode reduction threshold: ', &
        & equation%material%mode_reduction%threshold

      ! Read the material mapping from the equation table. A mapping is a
      ! reference from an (anonymous) variable to a material parameter.
      call tem_variable_loadMapping(                               &
        &    possVars            = initMaterial%poss_materialVars, &
        &    conf                = conf,                           &
        &    parent              = tblHandle,                      &
        &    key                 = "material",                     &
        &    varSys              = equation%varSys,                &
        &    varDict             = initMaterial%materialDict       )

      call aot_table_close( conf, tblHandle )

    end if
    call atl_append_newMaterialVars(                   &
      & varSys       = equation%varSys,                &
      & varSys_data  = varSys_data,                    &
      & poss_matVars = initMaterial%poss_materialVars, &
      & materialFun  = equation%material,              &
      & variables    = initMaterial%materialDict       )

    ! Material to use for decision whether to reduce order if constant.
    ! Default is to not use mode reduction, for this we set matpos to
    ! a negative value, so it will not match any material.
    equation%material%mode_reduction%matPos = -1
    if (equation%material%mode_reduction%enabled) then
      select case(equation%eq_kind)
      case('euler', 'navier_stokes','euler_2d', 'navier_stokes_2d', &
        &     'filtered_navier_stokes', 'euler_1d', 'filtered_navier_stokes_2d')
        ! Base decision for reduced order on Chi (masking), which will be found
        ! as first material in the flow equations.
        equation%material%mode_reduction%matPos = 1
      end select
    end if

    ! Store the positions of variables acting as data providers for materials
    ! in a sepparate list.
    allocate(materials(equation%material%nMat))
    do iMat=1, equation%material%nMat
      materials(iMat) = equation%varSys                                  &
        &                       %varname                                 &
        &                       %val( equation%material%matvar_pos(iMat) )
    end do

    ! Assing the material variables to the elements
    call atl_assign_elem_varProp(          &
      & tree           = tree,             &
      & mesh_list      = mesh_list,        &
      & affection_list = material_list,    &
      & levelPointer   = levelPointer,     &
      & varSys         = equation%varSys,  &
      & materialFun    = equation%material )

    ! create the element index list differentiating elements with
    ! constant and non-constant material parameters.
    call atl_create_materialElemList(      &
      & tree          = tree,              &
      & levelPointer  = levelPointer,      &
      & material_list = material_list,     &
      & materialFun   = equation%material, &
      & varSys        = equation%varSys,   &
      & materials     = materials          )

    nScalars = sum(equation%material%nScalars)

    do iLevel = tree%global%minlevel, tree%global%maxlevel
      do iProc = 1, mesh_list(iLevel)%descriptor%recvbuffer%nProcs
        call commPattern%initBuf_int(       &
          &    me    = mesh_list(iLevel)    &
          &              %descriptor        &
          &              %recvBuffer        &
          &              %buf_int(iProc),   &
          &    pos   = mesh_list(iLevel)    &
          &              %descriptor        &
          &              %recvBuffer        &
          &              %elemPos(iProc)    &
          &              %val,              &
          &    nVals = mesh_list(iLevel)    &
          &              %descriptor        &
          &              %recvBuffer        &
          &              %nElemsProc(iProc) )
      end do
      do iProc = 1, mesh_list(iLevel)%descriptor%sendbuffer%nProcs
        call commPattern%initBuf_int(       &
          &    me    = mesh_list(iLevel)    &
          &              %descriptor        &
          &              %sendBuffer        &
          &              %buf_int(iProc),   &
          &    pos   = mesh_list(iLevel)    &
          &              %descriptor        &
          &              %sendBuffer        &
          &              %elemPos(iProc)    &
          &              %val,              &
          &    nVals = mesh_list(iLevel)    &
          &              %descriptor        &
          &              %sendBuffer        &
          &              %nElemsProc(iProc) )
      end do

      ProjectionPos = material_list(iLevel)%poly_proj_pos
      nElemsCons = material_list(iLevel)%material_desc                 &
        &                               %computeElems(atl_ConstMatIdx) &
        &                               %nElems
      nElemsVar = material_list(iLevel)%material_desc               &
        &                              %computeElems(atl_VarMatIdx) &
        &                              %nElems
      material_list(iLevel)%material_dat                      &
        &                  %elemMaterialData(atl_ConstMatIdx) &
        &                  %nPointsPerDir &
        & = 1
      material_list(iLevel)%material_dat                  &
        &                  %elemMaterialData(2)           &
        &                  %nPointsPerDir                 &
        & = poly_proj_list(ProjectionPos)%nquadpointsperDir

      select case(scheme_list(iLevel)%scheme)
      case(atl_modg_scheme_prp)
        ! get correct amount of  number of points on the face
        nquadpoints = poly_proj_list(ProjectionPos)%body_3D%nquadpoints
      case(atl_modg_2d_scheme_prp)
        ! get correct amount of  number of points on the face
        nquadpoints = poly_proj_list(ProjectionPos)%body_2D%nquadpoints
      case(atl_modg_1d_scheme_prp)
        ! get correct amount of  number of points on the face
        nquadpoints = poly_proj_list(ProjectionPos)%body_1D%nquadpoints
      case default
        call tem_abort( 'ERROR in atl_init_materialParams: scheme not known' )
      end select

      allocate( material_list(iLevel)%material_dat                           &
        &                            %elemMaterialData(atl_ConstMatIdx)      &
        &                            %materialDat( nElemsCons, 1, nScalars ) )
      allocate( material_list(iLevel)%material_dat                    &
        &                            %elemMaterialData(atl_varMatIdx) &
        &                            %materialDat( nElemsVar,         &
        &                                          nquadpoints,       &
        &                                          nScalars )         )
      allocate( material_list(iLevel)%material_dat%mode_reducable( &
        &         mesh_list(iLevel)%descriptor%nElems)             )
      material_list(iLevel)%material_dat%mode_reducable = .false.

      ! Evaluate all material properties for the fluid elements (i.e. the
      ! element that are added to the computeElems variable of material_desc
      ! on each level)
      projectionPos = material_list(iLevel)%poly_proj_pos
      call atl_evalElemMaterial(                             &
        & mesh           = mesh_list(iLevel),                &
        & scheme         = scheme_list(iLevel),              &
        & material       = material_list(iLevel),            &
        & materialFun    = equation%material,                &
        & time           = time,                             &
        & mode_reduction = equation%material%mode_reduction, &
        & poly_proj      = poly_proj_list(projectionPos),    &
        & proc           = proc,                             &
        & commPattern    = commPattern                       )

      do iDir = 1, equation%nDimensions

        nFaces = mesh_list(iLevel)%faces%faces(iDir)%faceList%faceId%nVals

        ! Allocate space for the face material information for all faces, i.e.
        ! fluid, ghost, halos, etc. in advance
        ! The allocation has to be done before calling atl_assign_face_matPrp
        ! as this call is done levelwise, but currentLevel+1 is being accessed
        ! inside the routine. So we need to allocate the complete array
        ! beforehand.
        allocate( material_list(iLevel)%material_desc                 &
          &                            %material_face( iDir )         &
          &                            %mat( nFaces,                  &
          &                                  2,                       &
          &                                  equation%material%nMat ) )
        material_list(iLevel)%material_desc &
          &                  %material_face( iDir ) &
          &                  %mat &
          &                  %matVarPos = 0
        material_list(iLevel)%material_desc &
          &                  %material_face( iDir ) &
          &                  %mat &
          &                  %stFunPos = 0
      end do
    end do

    ! Now material lists for all levels have been created and we can
    ! inherit information from one to the next as needed.
    do iLevel = tree%global%minlevel, tree%global%maxlevel
      ! Assign the material properties to all faces in a recursive way with
      ! inheritance to cover local refinement)
      ! Note that this includes modifications of the next finer level.
      call atl_assign_face_matPrp(              &
        & mesh_list     = mesh_list,            &
        & material_list = material_list,        &
        & materialFun   = equation%material,    &
        & spatial_dim   = equation%nDimensions, &
        & currentLevel  = iLevel,               &
        & maxLevel      = tree%global%maxlevel, &
        & minLevel      = tree%global%minlevel  )

      ! Init the integer buffers to communicate the material
      ! information.
      call atl_init_faceStateBuffer( nFaceDofs   = 1,                       &
        &                            faces       = mesh_list(iLevel)%faces, &
        &                            nValsState  = 1,                       &
        &                            nValsFlux   = 1,                       &
        &                            boundary    = boundary_list(iLevel),   &
        &                            initRealBuf = .false.,                 &
        &                            commPattern = commPattern              )

      ! Communicate the material properties for the faces
      do iMat = 1, equation%material%nMat
        do iDir = 1, equation%nDimensions
          do iLR = tem_left, tem_right
            call commPattern%exchange_int(                                &
              & send         = mesh_list(iLevel)%faces                    &
              &                                 %faces(iDir)              &
              &                                 %sendBuffer_state(iLR),   &
              & recv         = mesh_list(iLevel)%faces                    &
              &                                 %faces(iDir)              &
              &                                 %recvBuffer_state(iLR),   &
              & state        = material_list(iLevel)%material_desc        &
              &                                     %material_face(iDir)  &
              &                                     %mat(:,:,iMat)        &
              &                                     %stFunPos,            &
              & message_flag = iLevel,                                    &
              & comm         = proc%comm                                  )
            call commPattern%exchange_int(                                &
              & send         = mesh_list(iLevel)%faces                    &
              &                                 %faces(iDir)              &
              &                                 %sendBuffer_state(iLR),   &
              & recv         = mesh_list(iLevel)%faces                    &
              &                                 %faces(iDir)              &
              &                                 %recvBuffer_state(iLR),   &
              & state        = material_list(iLevel)%material_desc        &
              &                                     %material_face(iDir)  &
              &                                     %mat(:,:,iMat)        &
              &                                     %matVarPos,           &
              & message_flag = iLevel,                                    &
              & comm         = proc%comm                                  )
          end do

          ! There are pointers included in the data exchange above. But the data
          ! exchange can be across ranks, so we have to redo the pointers to
          ! make sure that we have the correct adresses for the current rank.
          call atl_reassignStFunPtr( varSys        = equation%varSys, &
            &                        tree          = tree,            &
            &                        material_list = material_list,   &
            &                        iLevel        = iLevel,          &
            &                        iDir          = iDir,            &
            &                        iMat          = iMat             )

        end do
      end do

      ! Create new compute lists, depending on the combinations
      ! of material properties.
      call atl_create_materialComputeList(     &
        & mesh        = mesh_list(iLevel),     &
        & spatial_dim = equation%nDimensions,  &
        & material    = material_list(iLevel), &
        & materialFun = equation%material      )

      ! Evaluate all material properties for the compute faces (left + right
      ! element). Allocate face data in materialDat in advance
      ProjectionPos = material_list(iLevel)%poly_proj_pos
      do iDir = 1, equation%nDimensions
        select case(scheme_list(iLevel)%scheme)
        case(atl_modg_scheme_prp)
          ! detect the correct number of points on the face due to projection
          ! it is assumed that on left and right face thre is the same amount
          ! of points further is the number of points in a spatial direction
          ! for all dimension equal
          nquadpoints = poly_proj_list(ProjectionPos)%body_3D       &
            &                                        %faces(idir,2) &
            &                                        %nquadpoints
        case(atl_modg_2d_scheme_prp)
          nquadpoints = poly_proj_list(ProjectionPos)%body_2D       &
            &                                        %faces(idir,2) &
            &                                        %nquadpoints
        case(atl_modg_1d_scheme_prp)
          nquadpoints = poly_proj_list(ProjectionPos)%body_1D       &
            &                                        %faces(idir,2) &
            &                                        %nquadpoints
        case default
          call tem_abort( 'ERROR in atl_initMaterialParams: not able to ' &
            & // 'get nquadpoints for evaluate material parameters for '  &
            & // 'this scheme, stopping ...'                              )
        end select

        nConstFace = size(                                                &
          & material_list(iLevel)%material_desc                           &
          &                      %computeFace(iDir, atl_pureConstMat_prp) &
          &                      %facePos                                 )
        nVarFace = size(                                             &
          & material_list(iLevel)%material_desc                      &
          &                      %computeFace(iDir,atl_mixedMat_prp) &
          &                      %facePos                            )

        allocate(material_list(iLevel)%material_dat                          &
          &                           %faceMaterialData(iDir, atl_varMatIdx) &
          &                           %leftElemMaterialDat( nVarFace,        &
          &                                                 nquadpoints,     &
          &                                                 nScalars)        )
        allocate(material_list(iLevel)%material_dat                          &
          &                           %faceMaterialData(iDir, atl_varMatIdx) &
          &                           %rightElemMaterialDat( nVarFace,       &
          &                                                  nquadpoints,    &
          &                                                  nScalars)       )
        allocate( &
          & material_list(iLevel)%material_dat                            &
          &                      %faceMaterialData(iDir, atl_constMatIdx) &
          &                      %leftElemMaterialDat( nConstFace,        &
          &                                            1,                 &
          &                                            nScalars)          )
        allocate( &
          & material_list(iLevel)%material_dat                            &
          &                      %faceMaterialData(iDir, atl_constMatIdx) &
          &                      %rightElemMaterialDat( nConstFace,       &
          &                                             1,                &
          &                                             nScalars)         )
      end do

      call atl_evalFaceMaterial(                             &
        & mesh           = mesh_list(iLevel),                &
        & scheme         = scheme_list(iLevel),              &
        & material       = material_list(iLevel),            &
        & materialFun    = equation%material,                &
        & time           = time,                             &
        & spatial_dim    = equation%nDimensions,             &
        & poly_proj      = poly_proj_list(ProjectionPos)     )

      ! Create list of boundary faces, depending on the combinations
      ! of material properties.
      nBCs = boundary_list(iLevel)%nBcs
      allocate( material_list(iLevel)%material_desc              &
        &                            %bnd_faces(atl_constMatIdx) &
        &                            %boundary                   &
        &                            %bnd( nBCs )                )
      allocate( material_list(iLevel)%material_desc            &
        &                            %bnd_faces(atl_varMatIdx) &
        &                            %boundary                 &
        &                            %bnd( nBCs )              )
      call atl_create_materialBoundaryList(    &
        & material    = material_list(iLevel), &
        & materialFun = equation%material,     &
        & mesh        = mesh_list(iLevel),     &
        & boundary    = boundary_list(iLevel)  )
    end do

    ! Check for Maxwell equations
    select case(equation%eq_kind)
    case('maxwell','maxwell_2d')
      ! Get the maximum information propagation speed in the domain, e.g. the
      ! speed of light which soleny depends on the material parameters.
      call atl_getMaxPropSpeed( tree          = tree,              &
        &                       materialFun   = equation%material, &
        &                       material_list = material_list      )

    case('maxwelldivcorrection')
      ! Get the maximum information propagation speed in the domain, e.g. the
      ! speed of light which soleny depends on the material parameters.
      call atl_getMaxPropSpeedDivCor( tree          = tree,              &
        &                             materialFun   = equation%material, &
        &                             material_list = material_list      )

    end select

    deallocate(materials)

    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine atl_init_materialParams
  ! ****************************************************************************


  ! ****************************************************************************
  !> Assign reference to spacetime functions to the affected elements. This
  !! means: The position of the variable in the variable system, which reflects
  !! the spacetime function is determined and assigned to the element_list.
  !!
  !! This assignment is done as follows:
  !! 1. Set each fluid element's spacetime function pointer to null
  !! 2. Iterate over all material variables
  !! 2.1 Iterate over the first input variable's spacetime function list
  !! 2.1.a If it is a global shape spacetime function, have all element's
  !!       material stfun pointer set to this spacetime function.
  !! 2.1.b Take the spacetime function subtree and assign this spacetime
  !!       function to all element's material stfun pointer in this subtree.
  !! 2.2 Iterate over all elements and check whether their stfun pointer is
  !!     associated, which means, that the material is defined for the element.
  !!
  !! If there is one or more elements that are not covered by all materials, the
  !! solver will abort.
  subroutine atl_assign_elem_varProp( tree, mesh_list, affection_list,  &
    &                                 levelPointer, varSys, materialFun )
    ! ---------------------------------------------------------------------------
    !> The global description of the tree.
    type(treelmesh_type), intent(in) :: tree
    !> List of  mesh for different kernels
    type(atl_cube_elem_type), intent(in)  :: mesh_list( &
             &  tree%global%minlevel:tree%global%maxlevel)
    !> The description of the material properties on the element basis.
    type(atl_material_type), intent(inout) :: &
             & affection_list(tree%global%minlevel:tree%global%maxlevel)
    integer, intent(in)  :: levelPointer(:)
    !> global variable system to which luaVar to be appended
    type(tem_varSys_type), intent(in) :: varSys
    !> The list of material variable indizes.
    type(atl_materialFun_type), intent(in) :: materialFun
    ! ---------------------------------------------------------------------------
    integer :: nfluids, elemPos, levelElemPos
    integer :: iLevel, iMat, iElem, iStFun, varPos, inputVarPos
    type(tem_st_fun_listElem_type), pointer :: stFunList
    logical :: allElemsDefined
    ! ---------------------------------------------------------------------------

    do iLevel = tree%global%minlevel,tree%global%maxlevel
      ! Create array for all fluid properties
      nfluids = mesh_list(ilevel)%descriptor%elem%nElems(eT_fluid)
      allocate( affection_list(iLevel)%material_desc%material_elems( &
        & nfluids, size(materialFun%matvar_pos) ) )
      ! We don't have a global material anymore, thus we assign null here. Later
      ! on we check for elements with null, which indicates holes in the user
      ! defined shapes. If we find those holes, we report them as error.
      do iElem = 1, nFluids
        do iMat = 1, size( materialFun%matvar_pos )
          affection_list(iLevel)%material_desc               &
            &                   %material_elems(iElem, iMat) &
            &                   %stFunPtr                    &
            & => null()

        end do
      end do
    end do

    ! Iterate over all the materials and overwrite the material properties
    ! that have been defined before hand.
    do iMat = 1, size(materialFun%matvar_pos)
      ! get the position of the spacetime function in the varSys for the
      ! current material
      varPos = materialFun%matvar_pos(iMat)

      write(logUnit(lldebug),*) "Material at index ", iMat,    &
        & " refers to variable at index ", varPos, " called ", &
        & trim(varSys%varName%val(varPos)), " with ",          &
        & varSys%method%val(varPos)%nComponents, " components"

      ! We expect only one input variable, which is the variable containing the
      ! spacetime functions. Thus we hard-code index 1 here.
      inputVarPos = varSys%method%val(varPos)%input_varPos(1)

      write(logUnit(lldebug),*) "This variable uses the variable at index ", &
        & inputVarPos, " called ", trim(varSys%varName%val(inputVarPos)),    &
        & " as input."

      ! Now that we have the position of the lua variable, we can cast the
      ! list of spacetime functions for this variable, which is pointer to
      ! by the method data, back into a Fortran type to iterate over it.
      call c_f_pointer(varSys%method%val(inputVarPos)%method_data, stFunList)

      do iStFun = 1, stFunList%nVals
        if (stFunList%val(iStFun)%subTree%useGlobalMesh) then
          do iElem = 1, tree%nElems
            ! Make use of the levelpointer to read out level and position in the
            ! levewise list of elements.
            iLevel = tem_levelOf( tree%treeid(iElem) )
            levelElemPos = levelPointer(iElem)
            call assignStFunPointer( affection_list = affection_list, &
              &                      stFunList      = stFunList,      &
              &                      iLevel         = iLevel,         &
              &                      iMat           = iMat,           &
              &                      iStFun         = iStFun,         &
              &                      levelElemPos   = levelElemPos,   &
              &                      inputVarPos    = inputVarPos     )
          end do
        else
          do iElem = 1, stFunList%val(iStFun)%subTree%nElems
            ! We read out the position in current partition of the original
            ! tree.
            elemPos = stFunList%val(iStFun)%subTree%map2global(iElem)
            ! Make use of the levelpointer to read out level and position in the
            ! levewise list of elements.
            iLevel = tem_levelOf( tree%treeid(elemPos) )
            levelElemPos = levelPointer(elemPos)
            call assignStFunPointer( affection_list = affection_list, &
              &                      stFunList      = stFunList,      &
              &                      iLevel         = iLevel,          &
              &                      iMat           = iMat,           &
              &                      iStFun         = iStFun,         &
              &                      levelElemPos   = levelElemPos,   &
              &                      inputVarPos    = inputVarPos     )
          end do
        end if
      end do
      allElemsDefined = .true.
      do iElem = 1, tree%nElems
        ! Make use of the levelpointer to read out level and position in the
        ! levewise list of elements.
        iLevel = tem_levelOf( tree%treeid(iElem) )
        levelElemPos = levelPointer(iElem)
        allElemsDefined = allElemsDefined                                &
          &                 .and. associated(                            &
          &                         affection_list(ilevel)               &
          &                           %material_desc                     &
          &                           %material_elems(levelElemPos,iMat) &
          &                           %stFunPtr                          )

      end do
      if (.not. allElemsDefined) then
        write(logUnit(1),*) 'ERROR: Material ',                    &
              &                 trim(varSys%varName%val(varPos)),  &
              &                 ' is not defined for all elements'
        write(logUnit(1),*) 'Maybe you forgot to define a global const' &
          &                 // ' in your space-time function!'

        call tem_abort()
      end if

  end do

    contains

      ! *************************************************************************
      subroutine assignStFunPointer( affection_list, stFunList, iLevel, iMat, &
        &                            iStFun, levelElemPos, inputVarPos        )
        ! -----------------------------------------------------------------------
        type(atl_material_type), intent(inout) :: &
                 & affection_list(tree%global%minlevel:tree%global%maxlevel)
        type(tem_st_fun_listElem_type), pointer :: stFunList
        integer :: iLevel
        integer :: iMat
        integer :: iStFun
        integer :: levelElemPos
        integer :: inputVarPos
        ! -----------------------------------------------------------------------
        affection_list(iLevel)%material_desc                      &
          &                   %material_elems(levelElemPos, iMat) &
          &                   %stFunPtr                           &
          & => stFunList%val(iStFun)
        ! The following information is used for communication. We can't send
        ! the stFunPtr to another rank, thus we need some other way to
        ! address the spacetime functions. Therefore we also store the
        ! index of the variable as well as the index of the particular
        ! spacetime function within the variable's spacetime function
        ! list.
        affection_list(iLevel)%material_desc                      &
          &                   %material_elems(levelElemPos, iMat) &
          &                   %stFunPos                           &
          & = iStFun
        affection_list(iLevel)%material_desc                      &
          &                   %material_elems(levelElemPos, iMat) &
          &                   %matVarPos                          &
          & = inputVarPos
      end subroutine
      ! *************************************************************************

  end subroutine atl_assign_elem_varProp
  ! ****************************************************************************


  ! ****************************************************************************
  !> Read the configuration for the material paramters for Maxwell equations
  !! from configuration files and init the material parameter datatype.
  subroutine atl_update_materialParams( equation, mesh, scheme, boundary, &
    &                                   material, time, poly_proj, proc,  &
    &                                   commPattern                       )
    ! --------------------------------------------------------------------------
    !> Description on the equation system to solve.
    type(atl_Equations_type), intent(inout) :: equation
    !> Description of the mesh
    type(atl_cube_elem_type), intent(inout) :: mesh
    !> Information about the scheme
    type(atl_scheme_type), intent(in) :: scheme
    !> Global description of the boundary conditions.
    type(atl_level_boundary_type), intent(in)  :: boundary
    !> The description of the material properties. The compute lists in the
    !! material description is filled up by calling this subroutine.
    type(atl_material_type), intent(inout) :: material
    !> The current simulation time
    type(tem_time_type), intent(in) :: time
    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> mpi communication environment with mpi communicator
    type(tem_comm_env_type), intent(inout) :: proc
    !> mpi communication pattern type
    type(tem_commPattern_type), intent(inout) :: commPattern
    ! --------------------------------------------------------------------------

    if (equation%material%is_transient) then
      ! Evaluate all material properties for the fluid elements (i.e. the
      ! element that are added to the computeElems variable of material_desc
      ! on each level)
      call atl_evalElemMaterial(                             &
        & mesh           = mesh,                             &
        & scheme         = scheme,                           &
        & material       = material,                         &
        & materialFun    = equation%material,                &
        & time           = time,                             &
        & mode_reduction = equation%material%mode_reduction, &
        & poly_proj      = poly_proj,                        &
        & time_weights   = .true.,                           &
        & proc           = proc,                             &
        & commPattern    = commPattern                       )

      ! Evaluate all material properties for the compute faces (left + right
      ! element)
      call atl_evalFaceMaterial(                             &
        & mesh           = mesh,          &
        & scheme         = scheme,        &
        & material       = material,      &
        & materialFun    = equation%material,                &
        & time           = time,                             &
        & spatial_dim    = equation%nDimensions,             &
        & poly_proj      = poly_proj,                        &
        & time_weights   = .true.                            )

      ! Create list of boundary faces, depending on the combinations
      ! of material properties.
      call atl_create_materialBoundaryList( &
        & material     = material,           &
        & mesh         = mesh,               &
        & materialFun  = equation%material,  &
        & boundary     = boundary,           &
        & time_weights = .true.              )
    end if

  end subroutine atl_update_materialParams
  ! ****************************************************************************


  ! ****************************************************************************
  !> Subroutine to evaluate the material properties on the nodes
  !! of the compute faces (left and right element's trace).
  subroutine atl_evalFaceMaterial( mesh, scheme, material, materialFun, time, &
    &                              spatial_dim, poly_proj, time_weights       )
    ! ---
    !> Description of the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> Information about the scheme
    type(atl_scheme_type), intent(in) :: scheme
    !> The description of the material properties. The compute lists in the
    !! material description is filled up by calling this subroutine.
    type(atl_material_type), intent(inout) :: material
    !! The indizes of variables in the global varSys that are used as
    !! material's, penalization's or whatever's data sources.
    !!
    !! This data is needed to calculate the number of total material components.
    type(atl_materialFun_type), intent(in) :: materialFun
    !> The current simulation time
    type(tem_time_type), intent(in) :: time
    !> The spatial dimension
    integer, intent(in) :: spatial_dim
    !> Projection method for current level
    type(ply_poly_project_type), intent(inout) :: poly_proj
    logical, optional, intent(in) :: time_weights
    ! --------------------------------------------------------------------------
    integer :: iDir, iFace, iMat, iAlign
    integer :: faceIndex, leftIndex, rightIndex, res_lb, res_ub
    integer :: nMaxComps, compOffset
    type(atl_spacetime_fun_pointer_type) :: rightElemMat, leftElemMat
    real(kind=rk) :: bndBaryCoord(1,3)
    real(kind=rk), allocatable :: bary_matParLeft(:), bary_matParRight(:)
    integer :: nConstFace, nVarFace
    integer :: nquadpoints, elempos
    real(kind=rk), allocatable :: facePoints(:,:), facePhysCoord(:,:)
    logical :: use_timer
    ! --------------------------------------------------------------------------

    use_timer = .false.
    if (present(time_weights)) then
      use_timer = time_weights
    end if

    ! Store the max number of components a single material parameter has. This
    ! number is used to allocate buffers, as every material parameter fits
    ! into these buffers. We could allocate these buffers for every material
    ! parameter explicitly to not waste some memory in some cases, but
    ! allocating and deallocating the buffers is a costly operation. So for
    ! the sake of performance I decided to go for a slightly higher memory
    ! usage.
    nMaxComps = maxval(materialFun%nScalars)

    ! Use this to allocate a temporary array to avoid doing this for each
    ! materia parameter individually.
    allocate( bary_matParLeft( nMaxComps ), bary_matParRight( nMaxComps ) )

    matLoop: do iMat = 1, materialFun%nMat
      dirLoop: do iDir = 1, spatial_dim

        ! Calculate the component offsets for the current material iMat.These
        ! offsets are used for the assembled material array where all material
        ! parameters for one elemen reside next to each other.
        ! This calculation could have been done in the matLoop, but later on in
        ! the dirLoop res_lb and res_ub are potentially overwritten. Thus we
        ! have to renew it's value here.
        res_lb = sum(materialFun%nScalars(1:iMat-1))+1
        res_ub = sum(materialFun%nScalars(1:iMat))

        ! allocate space for the constant material parameters
        nConstFace = size(                                   &
          & material%material_desc                           &
          &         %computeFace(iDir, atl_pureConstMat_prp) &
          &         %facePos                                 )

        ! Iterate over all the compute faces with purely constant material
        ! parameters
        constFaceLoop: do iFace = 1, nConstFace


          ! Get the indices in the original compute face list
          faceIndex = material%material_desc                           &
            &                 %computeFace(iDir, atl_pureConstMat_prp) &
            &                 %facePos(iFace)
          leftIndex = material%material_desc                           &
            &                 %computeFace(iDir, atl_pureConstMat_prp) &
            &                 %leftPos(iFace)
          rightIndex = material%material_desc                           &
            &                  %computeFace(iDir, atl_pureConstMat_prp) &
            &                  %rightPos(iFace)

          ! Get the barycentric coordinate on the face (we use the barycentric
          ! coordinate of the left element for that). We check if one if the
          ! elements is a boundary. We use the one that is not a boundary
          ! element to calculate the barycentric coordinate.
          if ( leftIndex > mesh%descriptor%elem%nElems(eT_fluid)) then
            bndBaryCoord(1,1:3) = mesh%bary_coord(rightIndex,1:3)
            bndBaryCoord(1,iDir) = bndBaryCoord(1,iDir) - 0.5_rk*mesh%length
            elempos = mesh%descriptor%pntTID(rightIndex)
          else
            bndBaryCoord(1,1:3) = mesh%bary_coord(leftIndex,1:3)
            bndBaryCoord(1,iDir) = bndBaryCoord(1,iDir) + 0.5_rk*mesh%length
            elempos = mesh%descriptor%pntTID(leftIndex)
          end if
          if (use_timer) then
            call tem_startTimer( me          = atl_elemTimers, &
              &                  timerHandle = elemPos         )
          end if

          ! Read out the material positions for both elements (the left element
          ! property is stored
          ! as the right property and vice versa)
          leftElemMat = material%material_desc                  &
            &                   %material_face(iDir)            &
            &                   %mat( faceIndex, tem_left, iMat )
          rightElemMat = material%material_desc                   &
            &                    %material_face(iDir)             &
            &                    %mat( faceIndex, tem_right, iMat )

          ! Evaluate the constant at the barycentric face value (left and right
          ! limit)
          material%material_dat                                 &
            &     %faceMaterialData(iDir, atl_constMatIdx)      &
            &     %leftElemMaterialDat(iFace,1:1,res_lb:res_ub) &
            & = tem_spacetime_for(                              &
            &     me    = leftElemMat%stFunPtr,                 &
            &     coord = bndBaryCoord,                         &
            &     time  = time,                                 &
            &     n     = 1,                                    &
            &     nComp = materialFun%nScalars(iMat)            )
          material%material_dat                                  &
            &     %faceMaterialData(iDir, atl_constMatIdx)       &
            &     %rightElemMaterialDat(iFace,1:1,res_lb:res_ub) &
            & = tem_spacetime_for(                               &
            &     me    = rightElemMat%stFunPtr,                 &
            &     coord = bndBaryCoord,                          &
            &     time  = time,                                  &
            &     n     = 1,                                     &
            &     nComp = materialFun%nScalars(iMat)             )

          if (use_timer) then
            call tem_stopTimer( me          = atl_elemTimers, &
              &                 timerHandle = elemPos         )
          end if

        end do constFaceLoop

        ! Detect the correct number of points on the face due to projection;
        ! It is assumed that on left and right face there is the same amount of
        ! points. Further is the number of points in a spatial direction for all
        ! dimension equal
        select case(scheme%scheme)
        case(atl_modg_scheme_prp)
          nquadpoints = poly_proj%body_3D%faces(idir,2)%nquadpoints
        case(atl_modg_2d_scheme_prp)
          nquadpoints = poly_proj%body_2D%faces(idir,2)%nquadpoints
        case(atl_modg_1d_scheme_prp)
          nquadpoints = poly_proj%body_1D%faces(idir,2)%nquadpoints
        case default
          call tem_abort( 'ERROR in atl_evalFaceMaterial: not able to '      &
            & // 'evaluate face nodal material parameters for this scheme, ' &
            & // 'stopping ...'                                              )
        end select


        ! allocate space for the constant material parameters
        nVarFace = size( material%material_desc                       &
          &                      %computeFace(iDir, atl_mixedMat_prp) &
          &                      %facePos                             )
        allocate( facePhysCoord(nquadpoints,3) )


        ! Iterate over all the non-constant faces
        varFaceLoop: do iFace = 1, nVarFace

          ! Get the indices in the original compute face list
          faceIndex = material%material_desc                       &
            &                 %computeFace(iDir, atl_mixedMat_prp) &
            &                 %facePos(iFace)
          leftIndex = material%material_desc                       &
            &                 %computeFace(iDir, atl_mixedMat_prp) &
            &                 %leftPos(iFace)
          rightIndex = material%material_desc                       &
            &                  %computeFace(iDir, atl_mixedMat_prp) &
            &                  %rightPos(iFace)

          ! At least one of them has to be a fluid element, so we check which
          ! one and use the fluid element to shift the reference Chebyshev face
          ! nodes to the physical element
          if( leftIndex.le.mesh%descriptor%elem%nElems(eT_fluid) ) then
            ! Get the barycentric coordinate of the left element
            bndBaryCoord(1,1:3) = mesh%bary_coord(leftIndex,1:3)
            iAlign = 2
            elempos = mesh%descriptor%pntTID(leftIndex)
          else if(rightIndex.le.mesh%descriptor%elem%nElems(eT_fluid)) then
            ! Get the barycentric coordinate of the right element
            bndBaryCoord(1,1:3) = mesh%bary_coord(rightIndex,1:3)
            iAlign = 1
            elempos = mesh%descriptor%pntTID(rightIndex)
          else
            call tem_abort( 'ERROR in atl_evalFaceMaterial: compute '  &
              & // 'face without adjacent fluid element, stopping ...' )
          end if

          if (use_timer) then
            call tem_startTimer( me          = atl_elemTimers, &
              &                  timerHandle = elemPos         )
          end if

          ! get the quadrature points on the boundary faces
          select case(scheme%scheme)
          case(atl_modg_scheme_prp)
            call ply_get_quadpoints_faces( poly_proj = poly_proj, &
              &                            iDir      = iDir,      &
              &                            iAlign    = iAlign,    &
              &                            points    = facePoints )
          case(atl_modg_2d_scheme_prp)
            call ply_get_quadpoints_faces_2d( poly_proj = poly_proj, &
              &                               iDir      = iDir,      &
              &                               iAlign    = iAlign,    &
              &                               points    = facePoints )
          case(atl_modg_1d_scheme_prp)
            call ply_get_quadpoints_faces_1d( poly_proj = poly_proj, &
              &                               iDir      = iDir,      &
              &                               iAlign    = iAlign,    &
              &                               points    = facePoints )
          case default
            call tem_abort( errorMsg = 'ERROR in atl_evalFaceMaterial: not ' &
              & // 'able to get facial quadpoints to evaluate material '     &
              & // 'parameters for this scheme, stopping ...'                )
          end select

          ! Shift the face reference nodes to the physical face. For the left
          ! element, we have to use the reference face points on its right
          ! face (and vice versa).
          call atl_refToPhysCoord( refpoints  = facePoints,   &
            &                      nPoints    = nquadpoints,  &
            &                      baryCoord  = bndBaryCoord, &
            &                      elemLength = mesh%length,  &
            &                      physPoints = facePhysCoord )

          deallocate( facePoints )

          ! Read out the material positions for both elements (the left element
          ! property is stored as the right property and vice versa)
          leftElemMat = material%material_desc                  &
            &                   %material_face(iDir)            &
            &                   %mat( faceIndex, tem_left, iMat )
          rightElemMat = material%material_desc                   &
            &                    %material_face(iDir)             &
            &                    %mat( faceIndex, tem_right, iMat )

          compOffset = sum(materialFun%nScalars(1:iMat-1))
          res_lb = compOffset + 1
          res_ub = sum(materialFun%nScalars(1:iMat))

          ! Evaluate the material parameters at the physical nodes (left and
          ! right limit)
          if ( materialFun%nScalars(iMat) > 1 ) then
            ! Only call the vectorial routine when actually expecting a vector
            material%material_dat                                 &
              &     %faceMaterialData(iDir, atl_varMatIdx)        &
              &     %leftElemMaterialDat(iFace, :, res_lb:res_ub) &
              & = tem_spacetime_for(                              &
              &     me    = leftElemMat%stFunPtr,                 &
              &     coord = facePhysCoord,                        &
              &     time  = time,                                 &
              &     n     = nQuadPoints,                          &
              &     nComp = materialFun%nScalars(iMat)            )
            material%material_dat                                  &
              &     %faceMaterialData(iDir, atl_varMatIdx)         &
              &     %rightElemMaterialDat(iFace, :, res_lb:res_ub) &
              & = tem_spacetime_for(                               &
              &     me    = rightElemMat%stFunPtr,                 &
              &     coord = facePhysCoord,                         &
              &     time  = time,                                  &
              &     n     = nQuadPoints,                           &
              &     nComp = materialFun%nScalars(iMat)             )
          else
            ! Call the scalar routine when expecting only one value
            material%material_dat                          &
              &     %faceMaterialData(iDir, atl_varMatIdx) &
              &     %leftElemMaterialDat(iFace, :, res_lb) &
              & = tem_spacetime_for(                       &
              &     me    = leftElemMat%stFunPtr,          &
              &     coord = facePhysCoord,                 &
              &     time  = time,                          &
              &     n     = nQuadPoints                    )
            material%material_dat                           &
              &     %faceMaterialData(iDir, atl_varMatIdx)  &
              &     %rightElemMaterialDat(iFace, :, res_lb) &
              & = tem_spacetime_for(                        &
              &     me    = rightElemMat%stFunPtr,          &
              &     coord = facePhysCoord,                  &
              &     time  = time,                           &
              &     n     = nQuadPoints                     )
          end if

          if (use_timer) then
            call tem_stopTimer( me          = atl_elemTimers, &
              &                 timerHandle = elemPos         )
          end if

        end do varFaceLoop

        ! Free the temporary memory
        deallocate( facePhysCoord )

      end do dirLoop
    end do matLoop

    deallocate( bary_matParLeft, bary_matParRight )

  end subroutine atl_evalFaceMaterial
  ! ****************************************************************************


  ! ****************************************************************************
  !> Evaluates the material properties for all elements contained
  !! in the computeElems variable of the material_desc datatype.
  subroutine atl_evalElemMaterial( mesh, scheme, material, materialFun, time, &
    &                              poly_proj, mode_reduction, time_weights,   &
    &                              proc, commPattern                          )
    ! --------------------------------------------------------------------------
    !> Description of the mesh
    type(atl_cube_elem_type), intent(inout) :: mesh
    !> Information about the scheme
    type(atl_scheme_type), intent(in) :: scheme
    !> The description of the material properties. The compute lists in the
    !! material description is filled up by calling this subroutine.
    type(atl_material_type), intent(inout) :: material
    !! The indizes of variables in the global varSys that are used as
    !! material's, penalization's or whatever's data sources.
    !!
    !! This data is needed to calculate the number of total material components.
    type(atl_materialFun_type), intent(in) :: materialFun
    !> The current simulation time
    type(tem_time_type), intent(in) :: time
    !> Projection method for current level
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Settings for mode reduction. Used to determine which element can be
    !! calculated with reduced computational effort due to mode reduction.
    type(atl_mode_reduction_type), intent(in) :: mode_reduction
    logical, optional, intent(in) :: time_weights
    !> Communication environment
    type(tem_comm_env_type), intent(inout) :: proc
    !> mpi communication pattern type
    type(tem_commPattern_type), intent(inout) :: commPattern
    ! --------------------------------------------------------------------------
    real(kind=rk) :: bary_coord(1,3)
    real(kind=rk), allocatable :: volumePhysCoord(:,:), volumePoints (:,:)
    integer, allocatable :: commreduce(:)
    integer :: nElemsCons, nElemsVar, compOffset
    integer :: levelIndex, res_lb, res_ub
    integer :: iElem, iMat
    integer :: iFluid
    integer :: elemPos
    integer :: nquadpoints
    real(kind=rk) :: mode_reduction_val, min_mode_reduction_val
    logical :: use_timer
    integer :: nNeighbors
    integer :: iNeighbor
    integer :: neighbor
    logical :: fluidisreducable
    ! --------------------------------------------------------------------------

    use_timer = .false.
    if (present(time_weights)) then
      use_timer = time_weights
    end if

    !> @todo PV 20151027 Write unit tests for this routine as I'm not sure
    !!                   whether the index calculations around nScalars and iMat
    !!                   do what they are supposed to do.

    ! First we handle the constant material parameters
    nElemsCons = material%material_desc%computeElems(atl_constMatIdx)%nElems

    if( nElemsCons > 0 ) then
      constElemLoop: do iElem = 1, nElemsCons
        ! Get the element index in the total list and read out the
        ! position of the material property.
        levelIndex = material%material_desc                 &
          &                  %computeElems(atl_constMatIdx) &
          &                  %totElemIndices(iElem)

        !timer for LB of the global elements
        elempos = mesh%descriptor%pntTID(levelIndex)
        if (use_timer) then
          call tem_startTimer( me          = atl_elemTimers, &
            &                  timerHandle = elemPos         )
        end if

        ! Read out the barycenter of the element
        bary_coord(1,1:3) = mesh%bary_coord(levelIndex,1:3)

        ! We need to loop over the individual material parameters, because
        ! fortran is not capable of accepting an array of pointers as an
        ! argument. Thus we have to call tem_spacetime_for for each stFunPtr
        ! separately. The result of these calls is an array. The problem is to
        ! store the result at the correct indizes. Therefore we have to
        ! calculate the positions based on the number of components each
        ! material parameter has.
        do iMat = 1, materialFun%nMat
          res_lb = sum(materialFun%nScalars(1:iMat-1))+1
          res_ub = sum(materialFun%nScalars(1:iMat))
          !!if (.not. associated (material%material_desc%material_elems(levelindex, iMat) &
          !!        & %stFunPtr )) then
          !!    write(*,*) '++++++++POINTER is not set!!++++++++++'
          !!end if
          material%material_dat                                     &
            &     %elemMaterialData(atl_constMatIdx)                &
            &     %materialDat(iElem,1:1,res_lb:res_ub)             &
            & =  tem_spacetime_for(                                 &
            &      me    = material%material_desc                   &
            &                      %material_elems(levelIndex,iMat) &
            &                      %stFunPtr,                       &
            &      coord = bary_coord,                              &
            &      time  = time,                                    &
            &      n     = 1,                                       &
            &      nComp = materialFun%nScalars(iMat)               )
          if ( iMat == mode_reduction%matPos ) then
            mode_reduction_val = material%material_dat                      &
              &                          %elemMaterialData(atl_constMatIdx) &
              &                          %materialDat( iElem, 1, res_lb )
            material%material_dat%mode_reducable(levelIndex) = &
              & mode_reduction_val >= mode_reduction%threshold
          end if
        end do
        if (use_timer) then
          call tem_stopTimer( me          = atl_elemTimers,  &
            &                 timerHandle = elemPos          )
        end if
      end do constElemLoop

    end if

    ! Second, we handle the non-constant material parameters. This depends on
    ! the scheme we use, so we check for them and call separate routines for
    ! them.
    select case(scheme%scheme)
    case(atl_modg_scheme_prp)
      ! get correct amount of  number of points on the face
      nquadpoints = poly_proj%body_3D%nquadpoints
      allocate(volumePoints(nquadpoints,3))
      volumePoints = poly_proj%body_3d%nodes

    case(atl_modg_2d_scheme_prp)
      ! get correct amount of  number of points on the face
      nquadpoints = poly_proj%body_2D%nquadpoints
      allocate(volumePoints(nquadpoints,3))
      volumePoints = poly_proj%body_2d%nodes

    case(atl_modg_1d_scheme_prp)
      ! get correct amount of  number of points on the face
      nquadpoints = poly_proj%body_1D%nquadpoints
      allocate(volumePoints(nquadpoints,3))
      volumePoints = poly_proj%body_1d%nodes

    case default
      call tem_abort( 'ERROR in atl_evalElemMaterial: not able to evaluate '   &
        & // 'element nodal material parameters for this scheme, stopping ...' )
    end select

    ! For the variable material parameters we use the quadrature points
    ! provided by the FPT
    nElemsVar = material%material_desc%computeElems(atl_varMatIdx)%nElems

    if( nElemsVar > 0 ) then
      ! Allocate some temporary memory for the qudrature nodes on the physical
      ! element and the nodal material parameters at these points.
      allocate( volumePhysCoord(nquadpoints,1:3) )

      ! Now, we iterate over all the elements with variable coefficients and
      ! evaluate the material parameter for them.
      varElemLoop: do iElem = 1, nElemsVar

        ! Get the element position in the total list and recover the material
        ! position.
        levelIndex = material%material_desc               &
          &                  %computeElems(atl_varMatIdx) &
          &                  %totElemIndices(iElem)

        elempos = mesh%descriptor%pntTID(levelIndex)
        if (use_timer) then
          call tem_startTimer( me          = atl_elemTimers, &
            &                  timerHandle = elemPos         )
        end if
        !PV! Alter Code
        !PV!matPos = material%material_desc%material_elems(levelIndex)

        ! shift Chebyshev nodes to physical element coord
        call atl_refToPhysCoord( refPoints  = volumePoints,                   &
          &                      nPoints    = nquadpoints,                    &
          &                      baryCoord  = mesh%bary_coord(levelIndex, :), &
          &                      elemLength = mesh%length,                    &
          &                      physPoints = volumePhysCoord                 )

        do iMat = 1, materialFun%nMat
          compOffset = sum(materialFun%nScalars(1:iMat-1))
          res_lb = compOffset + 1
          res_ub = sum(materialFun%nScalars(1:iMat))

          if( materialFun%nScalars(iMat) > 1 ) then
            ! Only call the vectorial routine when actually expecting a vector
            material%material_dat                                    &
              &     %elemMaterialData(atl_varmatIdx)                 &
              &     %materialDat( iElem, :, res_lb:res_ub )          &
              & = tem_spacetime_for(                                 &
              &     me    = material%material_desc                   &
              &                     %material_elems(levelIndex,iMat) &
              &                     %stFunPtr,                       &
              &     coord = volumePhysCoord,                         &
              &     time  = time,                                    &
              &     n     = nQuadPoints,                             &
              &     nComp = materialFun%nScalars(iMat)               )
          else
            ! Call the scalar routine when only expecting one value
            material%material_dat                                    &
              &     %elemMaterialData(atl_varmatIdx)                 &
              &     %materialDat( iElem, :, res_lb )                 &
              & = tem_spacetime_for(                                 &
              &     me    = material%material_desc                   &
              &                     %material_elems(levelIndex,iMat) &
              &                     %stFunPtr,                       &
              &     coord = volumePhysCoord,                         &
              &     time  = time,                                    &
              &     n     = nQuadPoints                              )
          end if
          if ( iMat == mode_reduction%matPos ) then
            min_mode_reduction_val =                            &
              & minval(material%material_dat                    &
              &                %elemMaterialData(atl_varmatIdx) &
              &                %materialDat( iElem, :, res_lb ))
            material%material_dat%mode_reducable(levelIndex) = &
              & min_mode_reduction_val >= mode_reduction%threshold
          end if

        end do

        if (use_timer) then
          call tem_stopTimer( me          = atl_elemTimers,  &
            &                 timerHandle = elemPos          )
        end if

      end do varElemLoop

      deallocate(volumePhysCoord)
    end if

    ! Communicate the material%material_dat%mode_reducable info
    ! between processes.
    ! Use an integer to reuse existing communication buffer.
    allocate( commreduce(size(material%material_dat%mode_reducable)) )
    commreduce = 0
    where(material%material_dat%mode_reducable) commreduce = 1

    ! Information about neighbors provided by mesh%descriptor
    ! Communicate reducable flag across processes into halo elements.
    call commPattern%exchange_int(         &
      &    send         = mesh%descriptor  &
      &                       %sendBuffer, &
      &    recv         = mesh%descriptor  &
      &                       %recvBuffer, &
      &    state        = commreduce,      &
      &    message_flag = 123,             &
      &    comm         = proc%comm        )

    ! Run over all *fluid* elements, check reducable flag of all neighbors
    ! in the commreduce array and update mode_reducable accordingly.
    nNeighbors = size(mesh%descriptor%neigh(1)%nghElems, 1)
    do iFluid=1,mesh%descriptor%elem%nElems(eT_fluid)

      fluidisreducable = (commreduce(iFluid) == 1)
      do iNeighbor=1,nNeighbors
        neighbor = mesh%descriptor%neigh(1)%nghElems(iNeighbor,iFluid)

        if (neighbor > 0) then
          fluidisreducable = fluidisreducable &
            &                .and. (commreduce(neighbor) == 1)
        end if
      end do
      material%material_dat%mode_reducable(iFluid) = fluidisreducable

    end do

    ! Clean up the temporary memory for the next iteration
    deallocate(commreduce)
    deallocate(volumePoints)

  end subroutine atl_evalElemMaterial
  ! ****************************************************************************


  ! ****************************************************************************
  !> Subroutine to create element index list for constant and non-constant
  !! material parameters.
  !!
  !! To create the two lists, one list with elements with constant material
  !! parameters and one list with elements with variable material parameters,
  !! we have to follow several steps split into separate phases:
  !! Phase 1: Determine variable and constant elements
  !! 1. Loop over all materials
  !! 1.1. Loop over all spacetime functions
  !! 1.1.1 Determine whether the element is variable or not and add it to the
  !!       according list.
  !! Phase 2: Determine the levelwise element lists
  !! 2. Loop over all elements.
  !! 2.1. Determine the element's level.
  !! 2.2. Add them to intermediate lists.
  !! Phase 3: Create the result
  !! 3. Copy intermediate lists into result type.
  subroutine atl_create_materialElemList( tree, levelPointer, material_list, &
    &                                     varSys, materials, materialFun     )
    ! --------------------------------------------------------------------------
    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(in) :: tree
    !> The levelPointer contains the indizes on the levelwise fluid list for all
    !! treeID entries.
    integer, intent(in) :: levelPointer(:)
    !> The description of the material properties. This routine fills the
    !! compute lists in the material description.
    type(atl_material_type), intent(inout) &
      & :: material_list(tree%global%minLevel:tree%global%maxLevel)
    !> global variable system to which lua variables are to be appended
    type(tem_varSys_type), intent(in) :: varSys
    !> The list of the variables that should be used as materials.
    character(len=labelLen), intent(in) :: materials(:)
    !! The indizes of variables in the global varSys that are used as
    !! material's, penalization's or whatever's data sources.
    !!
    !! This data is needed to calculate the number of total material components.
    type(atl_materialFun_type), intent(inout) :: materialFun
    ! --------------------------------------------------------------------------
    integer :: iLevel, iElem, iStFun, iMat, iIndex, level
    integer :: nElemsConst, nElemsVar
    integer :: matPos, inVarPos
    logical :: allElemsAffected
    !> The position of an element according to the sorted array.
    integer :: sortedIndex
    !> Indicates that an element has variable material.
    logical :: elemHasVarMat
    type(dyn_intArray_type) :: variableElements
    type(dyn_intArray_type) ::                                                &
      & variableElementsPerLevel( tree%global%minLevel : tree%global%maxLevel )
    type(dyn_intArray_type) ::                                                &
      & constantElementsPerLevel( tree%global%minLevel : tree%global%maxLevel )
    type(tem_st_fun_listElem_type), pointer :: stFunList
    ! --------------------------------------------------------------------------

    !! Phase 1: Determine variable and constant elements

    !! this variable should be false if no material is there
    allElemsAffected = .false.
    do iMat = 1, materialFun%nMat
      ! get the position of the spacetime function in the varSys for the
      ! current material
      matPos = positionOfVal( me  = varSys%varName, &
        &                     val = materials(iMat) )
      ! The material variable just points to the variable contaning the
      ! spacetime functions. So we have to get this variable's index, which is
      ! the only input variable, so we will always find this variable at index 1
      inVarPos = varSys%method%val(matPos)%input_varPos(1)

      ! With the position of the lua variable, we can cast the list of
      ! spacetime functions for this variable, which is pointed to by
      ! method_data, back into a Fortran type to iterate over it.
      call c_f_pointer(varSys%method%val(inVarPos)%method_data, stFunList)

      ! Loop over all spacetime functions
      do iStFun = 1, stFunList%nVals
        write(logUnit(lldebug),'(A)') 'Kind of space-time function: ' &
          & // trim(stFunList%val(iStFun)%fun_kind)
        ! If the spacetime function is not constant, we need to add it's
        ! elements to the list of variable elements.
        if( stFunList%val(iStFun)%fun_kind /= 'const' ) then

          ! Check whether this material affects the whole mesh. If so, there
          ! is no need to examine further spacetime functions, because all
          ! elements have at least this material as variable material and thus
          ! all elements belong to the variable elements.
          allElemsAffected = stFunList%val(iStFun)%subTree%useGlobalMesh
          if( allElemsAffected ) exit

          do iElem = 1, stFunList%val(iStFun)%subTree%nElems
            call append(                                              &
              & me  = variableElements,                               &
              & val = stFunList%val(iStFun)%subTree%map2global(iElem) )
          end do

        end if

        if( allElemsAffected ) exit

      end do

    end do

    !! Phase 2: Determine the levelwise element lists
    iIndex = 1
    do iElem = 1, tree%nElems
      level = tem_levelOf(tree%treeID(iElem))

      if( allElemsAffected ) then
        elemHasVarMat = .true.
      else
        elemHasVarMat = .false.
        ! Check whether the index of the current element is in variableElements,
        ! which contains the indizes of variable Elements.
        if ( variableElements%nVals >= iIndex ) then
          sortedIndex = variableElements%sorted(iIndex)
          if( sortedIndex /= 0 ) then
            elemHasVarMat = variableElements%val(sortedIndex) == iElem
          end if
        end if
      end if
      if( elemHasVarMat ) then
        ! The element at index iElem is also part of the variableElements list
        ! and thus belongs to the variable elements
        call append( me  = variableElementsPerLevel(level), &
          &          val = levelPointer(iElem)              )
        ! As the current element is part of variableElements, we have to
        ! increase the pointer to the next possible element
        iIndex = iIndex + 1
      else
        ! The element at index iElem is not part of the variableElements list
        ! and thus belongs to the constant elements
        call append( me  = constantElementsPerLevel(level), &
          &          val = levelPointer(iElem)              )
      end if

    end do

    !! Phase 3: Create the result
    do iLevel = tree%global%minLevel, tree%global%maxLevel

      nElemsVar = variableElementsPerLevel(iLevel)%nVals
      nElemsConst = constantElementsPerLevel(iLevel)%nVals

      material_list(iLevel)%material_desc                 &
        &                  %computeElems(atl_ConstMatIdx) &
        &                  %nElems                        &
        & = nElemsConst
      allocate( material_list(iLevel)%material_desc                 &
        &                            %computeElems(atl_ConstMatIdx) &
        &                            %totElemIndices(nElemsConst)   )
      if( nElemsConst > 0 ) then
        material_list(iLevel)%material_desc                 &
          &                  %computeElems(atl_ConstMatIdx) &
          &                  %totElemIndices                &
          & = constantElementsPerLevel(iLevel)%val(:nElemsConst)
      end if

      material_list(iLevel)%material_desc               &
        &                  %computeElems(atl_VarMatIdx) &
        &                  %nElems                      &
        & = nElemsVar
      allocate(                                             &
        & material_list(iLevel)%material_desc               &
        &                      %computeElems(atl_VarMatIdx) &
        &                      %totElemIndices(nElemsVar)   )
      if( nElemsVar > 0 ) then
        material_list(iLevel)%material_desc                  &
          &                  %computeElems(atl_VarMatIdx)    &
          &                  %totElemIndices                 &
          & = variableElementsPerLevel(iLevel)%val(:nElemsVar)
      end if

    end do

    do iLevel = tree%global%minLevel, tree%global%maxLevel
      call destroy( variableElementsPerLevel( iLevel ) )
      call destroy( constantElementsPerLevel( iLevel ) )
    end do

  end subroutine atl_create_materialElemList
  ! ****************************************************************************


  ! ****************************************************************************
  !> Create separate compute list for constant-constant, constant-variable (or
  !! vice versa) and variable-variable material parameter compute faces on this
  !! rank.
  subroutine atl_create_materialBoundaryList( material, materialFun, boundary, &
    &                                         mesh, time_weights )
    ! --------------------------------------------------------------------------
    !> Description of the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The description of the material properties. The compute lists in the
    !! material description is filled up by calling this subroutine.
    type(atl_material_type), intent(inout) :: material
    !! The indizes of variables in the global varSys that are used as
    !! material's, penalization's or whatever's data sources.
    !!
    !! This data is needed to calculate the number of total material components.
    type(atl_materialFun_type), intent(in) :: materialFun
    !> Boundary description for the all the levels.
    type(atl_level_boundary_type), intent(in) :: boundary
    logical, optional, intent(in) :: time_weights
    ! --------------------------------------------------------------------------
    integer :: iBc, iDir, iAlign, iFace, nBCs, iMat
    integer :: facePos, neighPos
    type(atl_spacetime_fun_pointer_type) :: neighMat(materialFun%nMat)
    integer :: matType, elempos
    logical :: constMat(materialFun%nMat)
    logical :: use_timer
    ! --------------------------------------------------------------------------

    use_timer = .false.
    if (present(time_weights)) then
      use_timer = time_weights
    end if

    nBCs = boundary%nBcs
    material%material_desc%bnd_faces(:)%boundary%nBCs = nBCs
    bcLoop: do iBc = 1, nBcs
      dirLoop: do iDir = 1,3
        alignLoop: do iAlign = 1, 2

          ! Init the new datatype which holds the boundary material
          ! information.
          call init( me     = material%material_desc              &
            &                         %bnd_faces(atl_constMatIdx) &
            &                         %boundary                   &
            &                         %bnd(iBc)                   &
            &                         %faces(iDir,iAlign)         &
            &                         %facePos,                   &
            &        length = 16                                  )
          call init( me     = material%material_desc            &
            &                         %bnd_faces(atl_varMatIdx) &
            &                         %boundary                 &
            &                         %bnd(iBc)                 &
            &                         %faces(iDir,iAlign)       &
            &                         %facePos,                 &
            &        length = 16                                )
          call init( me     = material%material_desc              &
            &                         %bnd_faces(atl_constMatIdx) &
            &                         %boundary                   &
            &                         %bnd(iBc)                   &
            &                         %faces(iDir,iAlign)         &
            &                         %neighPos,                  &
            &        length = 16                                  )
          call init( me     = material%material_desc            &
            &                         %bnd_faces(atl_varMatIdx) &
            &                         %boundary                 &
            &                         %bnd(iBc)                 &
            &                         %faces(iDir,iAlign)       &
            &                         %neighPos,                &
            &        length = 16                                )

          faceLoop: do iFace = 1, boundary%bnd(iBc)           &
            &                             %faces(iDir,iAlign) &
            &                             %facePos            &
            &                             %nVals

            ! Get the position of the boundary face and neighboring fluid
            ! element's face
            facePos = boundary%bnd(iBc)%faces(iDir,iAlign)%facePos  &
              &                                           %val(iFace)
            neighPos = boundary%bnd(iBc)%faces(iDir,iAlign)%neighPos &
              &                                            %val(iFace)
            elempos = mesh%descriptor%pntTID(neighPos)
            if (use_timer) then
              call tem_startTimer( me          = atl_elemTimers, &
                &                  timerHandle = elemPos         )
            end if
            ! Now, we get the material property assigned for this face from
            ! its fluid neighbor element.
            neighMat = material%material_desc%material_elems(neighPos, :)

            matLoop: do iMat = 1, materialFun%nMat
              ! Check if this boundary face has the same material function from
              ! both sides
              if( .not. associated( neighMat(iMat)%stFunPtr ) ) then
                call tem_abort( 'ERROR in atl_create_materialBoundaryList: ' &
                  & // 'no material at boundary face, stopping ...'   )
              end if
              constMat(iMat) = neighMat(iMat)%stFunPtr%fun_kind == 'const'
            end do matLoop

            ! Locate the position of the face materials
            if ( all(constMat) ) then
              matType = atl_pureConstMat_prp
            else
              matType = atl_mixedMat_prp
            end if

            ! Store the gathered information.
            call append( me  = material%material_desc      &
              &                        %bnd_faces(matType) &
              &                        %boundary           &
              &                        %bnd(iBc)           &
              &                        %faces(iDir,iAlign) &
              &                        %facePos,           &
              &          val = facePos                     )
            call append( me  = material%material_desc      &
              &                        %bnd_faces(matType) &
              &                        %boundary           &
              &                        %bnd(iBc)           &
              &                        %faces(iDir,iAlign) &
              &                        %neighPos,          &
              &          val = neighPos                    )

            if (use_timer) then
              call tem_stopTimer( me          = atl_elemTimers, &
                &                 timerHandle = elemPos         )
            end if
          end do faceLoop
        end do alignLoop
      end do dirLoop
    end do bcLoop

  end subroutine atl_create_materialBoundaryList
  ! ****************************************************************************


  ! ****************************************************************************
  !> Create separate compute list for constant-constant, constant-variable (or
  !! vice versa) and variable-variable material parameter compute faces on this
  !! rank.
  subroutine atl_create_materialComputeList( mesh, spatial_dim, material, &
    &                                        materialFun )
    ! --------------------------------------------------------------------------
    !> Description of the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The spatial dimension
    integer, intent(in) :: spatial_dim
    !> The description of the material properties. The compute lists in the
    !! material description is filled up by calling this subroutine.
    type(atl_material_type), intent(inout) :: material
    !! The indizes of variables in the global varSys that are used as
    !! material's, penalization's or whatever's data sources.
    !!
    !! This data is needed to calculate the number of total material components.
    type(atl_materialFun_type), intent(in) :: materialFun
    ! --------------------------------------------------------------------------
    type(grw_intArray_type) :: faceLists(2)
    integer :: iDir, iFace, iList, iMat
    integer :: facePos
    logical :: left_isConst, right_isConst
    integer :: compute_type
    integer :: computeIndex, nFaces
    ! --------------------------------------------------------------------------

    ! Iterate over all the levels and the spatial directions and prepare the
    ! separate compute lists
    do iDir = 1,spatial_dim

      ! Init the temporary index lists
      call init( me = faceLists(1), length = 16 )
      call init( me = faceLists(2), length = 16 )

        ! We iterate over all the compute faces and check left and right
      ! material position. Then we check the combination of the adjacent
      ! materials and decide what we do with this face.
      faceLoop: do iFace = 1, size( mesh%faces%faces(iDir)%computeFace &
        &                                                 %facePos     )

        facePos = mesh%faces%faces(iDir)%computeFace%FacePos(iFace)

        ! Check if the left and right material property are constants
        left_isConst = .true.
        right_isConst = .true.
        do iMat = 1, materialFun%nMat
          if (associated(material%material_desc               &
              &                          %material_face(iDir) &
              &                          %mat( FacePos,       &
              &                                tem_left,      &
              &                                iMat )         &
              &                          %stFunPtr)           ) then
            left_isConst = left_isConst                          &
              & .and. 'const' == material%material_desc          &
              &                          %material_face(iDir)    &
              &                          %mat( FacePos,          &
              &                                tem_left,         &
              &                                iMat )            &
              &                          %stFunPtr               &
              &                          %fun_kind
          end if
          if (associated(material%material_desc               &
              &                          %material_face(iDir) &
              &                          %mat( FacePos,       &
              &                                tem_right,     &
              &                                iMat )         &
              &                          %stFunPtr)           ) then
            right_isConst = right_isConst                       &
              & .and. 'const' == material%material_desc         &
              &                          %material_face(iDir)   &
              &                          %mat( FacePos,         &
              &                                tem_right,       &
              &                                iMat )           &
              &                          %stFunPtr              &
              &                          %fun_kind
          end if
        end do !iMat

        ! Check if both are constants or at least one of them has a variable
        ! coefficient.
        if (left_isConst.and.right_isConst) then
          compute_type = atl_pureConstMat_prp
        else
          compute_type = atl_mixedMat_prp
        end if

        ! add this face index to the corresponding list
        call append( me = faceLists(compute_type), val = iFace )

      end do faceLoop

      ! Build up the final compute lists from the temporary arrays.
      ! We build one for the atl_pureConstMat_prp -> iList=1
      ! and another one for atl_mixedMat_prp -> iList=2 .
      do iList = 1, 2
        nFaces = faceLists(iList)%nVals

        if ( .not. allocated( material%material_desc           &
          &                           %computeFace(iDir,iList) &
          &                           %facePos) ) then
          allocate(                            &
            & material%material_desc           &
            &         %computeFace(iDir,iList) &
            &         %facePos(nFaces),        &
            & material%material_desc           &
            &         %computeFace(iDir,iList) &
            &         %leftPos(nFaces),        &
            & material%material_desc           &
            &         %computeFace(iDir,iList) &
            &         %rightPos(nFaces)        )
        end if

        do iFace = 1, nFaces
          ! Get the index in the original compute face list
          computeIndex = faceLists(iList)%val(iFace)
          ! Get the face position from the original face description
          material%material_desc           &
            &     %computeFace(iDir,iList) &
            &     %facePos(iFace)          &
            & = mesh%faces                 &
            &       %faces(iDir)           &
            &       %computeFace           &
            &       %facePos(computeIndex)
          ! Get the left element's position from the original face description
          material%material_desc           &
            &     %computeFace(iDir,iList) &
            &     %leftPos(iFace)          &
            & = mesh%faces                 &
            &       %faces(iDir)           &
            &       %computeFace           &
            &       %leftPos(computeIndex)
          ! Get the right element's position from the original face
          ! description
          material%material_desc           &
            &     %computeFace(iDir,iList) &
            &     %rightPos(iFace)         &
            & = mesh%faces                 &
            &       %faces(iDir)           &
            &       %computeFace           &
            &       %rightPos(computeIndex)
            end do !iFace
      end do !iList

      ! Clean up the temporary arrays on this level
      call destroy( me = faceLists(1) )
      call destroy( me = faceLists(2) )

    end do !iDir

  end subroutine atl_create_materialComputeList
  ! ****************************************************************************


  ! ****************************************************************************
  !> Define material properties for the faces of all fluid elements and inherit
  !! the face material property to all ghost elements on the finer level.
  subroutine atl_assign_face_matPrp( minLevel, maxLevel, mesh_list,            &
    &                                currentLevel, spatial_dim, material_list, &
    &                                materialFun                               )
    ! --------------------------------------------------------------------------
    integer,intent(in) :: minLevel, maxLevel
    !> List of  mesh for different kernels
    type(atl_cube_elem_type), intent(in) :: mesh_list(minlevel:maxlevel)
    !> currentLevel
    integer, intent(in) :: currentLevel
    !> The spatial dimension
    integer, intent(in) :: spatial_dim
    !> The description of the material properties on the element basis.
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxlevel)
    !! The indices of variables in the global varSys that are used as
    !! material's, penalization's or whatever's data sources.
    !!
    !! This data is needed to calculate the number of total material components.
    type(atl_materialFun_type), intent(in) :: materialFun
    ! --------------------------------------------------------------------------
    integer :: nFluids
    integer :: facepos
    integer :: iDir, iFace, iAlign, invAlign, iMat
    integer :: elemPos, childPos(4), lrElemPos(2)
    integer :: nChilds
    integer :: matvarpos
    logical :: failed
    ! --------------------------------------------------------------------------


    failed = .false.
    nfluids = mesh_list(currentlevel)%descriptor%elem%nElems(eT_fluid)

    ! Assign the material properties to all faces
    dirloop: do iDir = 1, spatial_dim

      ! Total number of faces on this level for this direction,
      ! including ghosts.
      ! (Virtual) boundary elements will have indices greater than
      ! this, and we need to set the material of the inner adjacent
      ! element for these.

      faceloop: do iFace = 1, size( mesh_list(currentLevel)%faces       &
        &                                                  %faces(iDir) &
        &                                                  %computeFace &
        &                                                  %facePos )

        ! The position of this face in the complete list of faces on the
        ! level.
        facepos = mesh_list(currentLevel)%faces       &
          &                              %faces(iDir) &
          &                              %computeFace &
          &                              %facePos(iFace)

        ! Get the left and right element positions
        lrElemPos(1) = mesh_list(currentLevel)%faces         &
          &                                   %faces(iDir)   &
          &                                   %computeFace   &
          &                                   %leftPos(iFace)
        lrElemPos(2) = mesh_list(currentLevel)%faces          &
          &                                   %faces(iDir)    &
          &                                   %computeFace    &
          &                                   %rightPos(iFace)

        ! Look up left(1) and right(2) elements for all faces.
        do iAlign=1,2
          invAlign = tem_invFace_map(iAlign)
          if (lrElemPos(iAlign) <= nFluids) then
            elempos = lrelempos(iAlign)
          else if (lrElemPos(invAlign) <= nFluids) then
            ! There is no fluid neighbor on this side.
            ! We try to use the material from of the element
            ! on the other side of the face instead.
            elempos = lrElemPos(invAlign)
          else
            elempos = nFluids+1
            failed = .true.
          end if

          if (elempos <= nFluids) then
            ! Adjacent element is a fluid element, use its material for
            ! this side of the face.
            do iMat=1,materialFun%nMat
              matvarpos = material_list(currentLevel)%material_desc &
                &           %material_face(iDir)                    &
                &           %mat( facepos, iAlign, iMat )           &
                &           %matvarpos
              if (matvarpos == 0) then
                ! Only set the material, if it was not set already.
                ! If there is a refinement, the material might already be
                ! set by the coarser level.
                material_list(currentLevel)%material_desc                  &
                  &                        %material_face(iDir)            &
                  &                        %mat( facepos,                  &
                  &                              iAlign,                   &
                  &                              iMat )                    &
                  & = material_list(currentLevel)%material_desc            &
                  &                              %material_elems( ElemPos, &
                  &                                               iMat     )
              end if
            end do

          end if

        end do

      end do faceloop

      if (failed) then
        write(logunit(1),*) 'Ran into a computed face without fluid neighbor!'
        call tem_abort()
      end if


      ! Inherit the material position for the faces which have children
      ! on the next level (i.e. for all from finer faces).
      !
      ! Fluxes are computed on the finest level, so we only need to set the side
      ! towards the coarser elements to ensure that we have material information
      ! on the finest level for the fine ghost faces.
      if (currentLevel < maxLevel) then
        nChilds = size(mesh_list(currentlevel)%faces%faces(iDir) &
          &                                         %faceDep &
          &                                         %childFacePos, 1)
        do iAlign = 1, 2
          do iFace = 1, size( mesh_list(currentLevel)%faces                 &
            &                                        %faces(iDir)           &
            &                                        %fromFinerFace(iAlign) &
            &                                        %facePos               )

            ! Position of this from finer face in the face list.
            facepos = mesh_list(currentLevel)%faces                 &
              &                              %faces(iDir)           &
              &                              %fromFinerFace(iAlign) &
              &                              %facePos(iFace)

            ! The children's face positions
            childPos(:nChilds) = mesh_list(currentLevel) &
              &                    %faces%faces(iDir)    &
              &                    %facedep%childfacepos(:,facepos)

            ! The parent's neighboring coarse element (on the opposite side of
            ! this from finer face side).
            elemPos = mesh_list(currentLevel)%faces                 &
              &                              %faces(iDir)           &
              &                              %fromFinerFace(iAlign) &
              &                              %elemPosOp(iFace)

            do iMat = 1,  materialFun%nMat
              ! Set the material of the child faces according to
              ! their parents element material.
              material_list(currentLevel+1)                   &
                & %material_desc                              &
                & %material_face(iDir)                        &
                & %mat(childPos(:nChilds), iAlign, iMat)      &
                & = material_list(currentLevel)%material_desc &
                &                              %material_elems(elemPos, iMat)
            end do
          end do
        end do
      end if
    end do dirloop

  end subroutine atl_assign_face_matPrp
  ! ****************************************************************************


end module atl_materialIni_module
