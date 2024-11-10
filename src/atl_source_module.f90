! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012, 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012-2016, 2018, 2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
! Copyright (c) 2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> author: Harald Klimach
!! author: Peter Vitt
!! author: Melven Zoellner
!! author: Jens Zudrop
!!
!! Provides routines to initialize and update source terms (accumulate the RHS
!! of the PDE)
!!
module atl_source_module
  use, intrinsic :: iso_c_binding,   only: c_loc,       &
    &                                      c_f_pointer, &
    &                                      c_ptr
  use env_module,                    only: rk, labellen
  use aotus_module,                  only: flu_state

  ! include treelm modules
  use tem_varSys_module,             only: tem_varSys_type,              &
    &                                      tem_varSys_append_derVar,     &
    &                                      tem_varSys_proc_point,        &
    &                                      tem_varSys_proc_element,      &
    &                                      tem_varSys_proc_setParams,    &
    &                                      tem_varSys_proc_getParams,    &
    &                                      tem_varSys_proc_setupIndices, &
    &                                      tem_varSys_proc_getValOfIndex
  use tem_varMap_module,             only: tem_possible_variable_type, &
    &                                      tem_variable_loadMapping
  use tem_stringKeyValuePair_module, only: grw_stringKeyValuePairArray_type, &
    &                                      append
  use tem_dyn_array_module,          only: PositionOfVal
  use treelmesh_module,              only: treelmesh_type
  use tem_element_module,            only: eT_fluid
  use tem_time_module,               only: tem_time_type
  use tem_aux_module,                only: tem_abort
  use tem_logging_module,            only: logUnit
  use tem_spacetime_fun_module,      only: tem_st_fun_listElem_type
  use tem_grow_array_module,         only: append
  use tem_varSys_module,             only: tem_varSys_solverData_evalElem_type
  use tem_stringKeyValuePair_module, only: tem_stringKeyValuePair_type, &
    &                                      init

  use ply_poly_project_module,       only: ply_poly_project_type, &
    &                                      ply_prj_body_type,     &
    &                                      assignment(=)

  ! include ateles modules
  use atl_equation_module,           only: atl_equations_type
  use atl_scheme_module,             only: atl_modg_scheme_prp,  &
    &                                      atl_scheme_type,      &
    &                                      atl_modg_2d_scheme_prp
  use atl_cube_elem_module,          only: atl_cube_elem_type
  use atl_materialPrp_module,        only: atl_material_type, &
    &                                      atl_ConstMatIdx
  use atl_varSys_module,             only: atl_varSys_solverData_type, &
    &                                      atl_set_stFun_getElement,   &
    &                                      atl_get_new_varSys_data_ptr
  use atl_operator_module,           only: atl_set_opVar_getElement
  use atl_source_types_module,       only: atl_source_type,        &
    &                                      atl_eqn_sourceMap_type, &
    &                                      atl_init_source_type
  use atl_reference_element_module,  only: atl_refToPhysCoord

  implicit none
  private


  public :: atl_update_sourcedata
  public :: atl_allocate_sourceData
  public :: atl_deallocate_sourceData
  public :: atl_initialize_sources
  public :: atl_append_newSourceVars
  public :: atl_fill_sourceIndex
  public :: atl_source_prim2cons


contains


! ******************************************************************************
  subroutine atl_initialize_sources( source, initSource, conf, equation, &
    &                                poly_proj_list, mesh_list, tree,    &
    &                                varSys_data                         )
    ! --------------------------------------------------------------------------
    !> Instance of atl_source_type to be initialized.
    !! This instance will contain the source definitions from lua as well as
    !! the corresponding variables in the global variable system.
    type(atl_source_type), intent(inout) :: source

    !> Initialize source type contains possible source terms and function
    !! pointers to update those source terms
    type(atl_init_source_type), intent(in) :: initSource

    !> lua state
    type(flu_state), intent(inout) :: conf

    !> Description on the equation system to solve.
    type(atl_Equations_type), intent(inout) :: equation

    !> unique list for projection methods
    type(ply_poly_project_type), intent(in) :: poly_proj_list(:)

    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(in) :: tree

    !> Mesh list to access the level descriptors
    type(atl_cube_elem_type), intent(in) :: mesh_list(tree%global%minLevel:)

    type(atl_varSys_solverData_type), intent(in), target :: varSys_data
    ! --------------------------------------------------------------------------
    type(tem_stringKeyValuePair_type) :: predef_src
    ! --------------------------------------------------------------------------

    ! Initilized keyValuePair dict for sources
    call init( me = source%varDict )

    call tem_variable_loadMapping(                     &
      & possVars            = initSource%poss_srcVars, &
      & conf                = conf,                    &
      & key                 = "source",                &
      & varSys              = equation%varSys,         &
      & varDict             = source%varDict           )

    ! In case the equation has a pre-defined/permament source.
    ! The key and the value is appended to the sourceDict here.
    ! Then, this source term would get added if it was added as a
    ! poss_srcVars
    ! @todo NA: This needs to fixed and done elegantly
    select case(equation%eq_kind)
      case('filtered_navier_stokes')
        predef_src%key = 'ransSource'
        predef_src%value = ''
        call append(me = source%varDict, val = predef_src)
      case('filtered_navier_stokes_2d')
        predef_src%key = 'rans2dsource'
        predef_src%value = ''
        call append(me = source%varDict, val = predef_src)
    end select

    ! As we need some information about mesh and projection during source term
    ! evaluation, we have to reference them here.
    call atl_append_newSourceVars( me           = source,                  &
      &                            varSys       = equation%varSys,         &
      &                            varSys_data  = varSys_data,             &
      &                            poss_srcVars = initSource%poss_srcVars, &
      &                            eval_source  = initSource%eval_source   )
    ! convert the sourceVars from prim to cons
    ! At the moment this is only for SpongeLayer which is checked inside this fn
    ! KM: Do conversion from prim 2 cons before setupIndices because for
    ! spatial function values are stored in st_fun setupIndices
    if (equation%hasPrimitiveVariables) then
      call atl_source_prim2cons( equation = equation,        &
        &                        source   = source,          &
        &                        varSys   = equation%varSys  )
    end if
    ! Fill source idx returned from setupIndices of variable to retreive
    ! a value later on in update_sourceTerm using get_valOfIndex.
    ! Source index are not filled for predefined or permanent source terms.
    call atl_fill_sourceIndex( source         = source,                 &
      &                        varSys         = equation%varSys,        &
      &                        nDim           = equation%nDimensions,   &
      &                        tree           = tree,                   &
      &                        poly_proj_list = poly_proj_list,         &
      &                        mesh_list      = mesh_list               )

  end subroutine atl_initialize_sources
! ******************************************************************************


! *******************************************************************************
  subroutine atl_append_newSourceVars( me, varSys, varSys_data, poss_srcVars,  &
    &                                  eval_source )
    ! ---------------------------------------------------------------------------
    !> Instance of atl_source_type to be initialized.
    !! This instance will contain the source definitions from lua as well as
    !! the corresponding variables in the global variable system.
    type(atl_source_type), intent(inout) :: me
    !> The variable system to which the souce variables have to be added
    type(tem_varSys_type), intent(inout) :: varSys
    !> Data for the variable System
    type(atl_varSys_solverData_type), intent(in), target :: varSys_data
    !> The list of possible source variables. This is used to determine the
    !! index of the eval-source_routine in eval_source.
    type(tem_possible_variable_type) :: poss_srcVars
    !> List of function pointers for each possible source to apply the source
    !! term to the state.
    type(atl_eqn_sourceMap_type), allocatable :: eval_source(:)
    ! ---------------------------------------------------------------------------
    integer :: addedSourcePosition, iSource, nComponents, iState, poss_srcVarPos
    integer :: nComp_data
    logical :: wasAdded
    character(labellen), allocatable :: inputVarNames(:)
    character(labellen) :: varname
    procedure(tem_varSys_proc_point), pointer :: get_point => null()
    procedure(tem_varSys_proc_element), pointer :: get_element => null()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => Null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => Null()
    procedure(tem_varSys_proc_setupIndices), pointer ::setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_ValOfIndex => NULL()
    ! ---------------------------------------------------------------------------

    ! List all state variables to use them as input variables for the routines
    ! that apply the source variables to the state.
    ! Additionally, the space-time-function for this source variable is also
    ! needed as an input variable. Thus we allocatate one additional slot for
    ! the name of the space-time-function.
    allocate(inputVarNames(varSys%nStateVars+1))
    do iState = 1, varSys%nStateVars
      write(logUnit(5),*) 'Adding state variable ' &
        & // trim(varSys%varname%val(iState))            &
        & // ' to the list of input variables'
      inputVarNames(iState) = varSys%varname%val(iState)
    end do

    ! Now we have to add the source variables. These are the 'real' source
    ! variables that can make use of the space-time-functions added earlier.
    ! Therefore we first have to allocate the slots in atl_source_type to store
    ! information like affection elements or functions to evaluate the source
    ! term.
    write(logUnit(9),'(A,I1,a)') 'Allocating space for ', me%varDict%nVals, &
      & ' sources'
    allocate( me%method( me%varDict%nVals ) )
    do iSource = 1, me%varDict%nVals

      ! get the name of the current variable
      varname = trim(me%varDict%val(iSource)%key)

      ! The function to apply the source term to the state has already been
      ! added to the eval_source list by the equation systems during their
      ! initialization of the source terms. So we need to find the position of
      ! the variable in the list of possible variables, as this is also the
      ! position of the function pointer in the eval_source list.
      poss_srcVarPos = PositionOfVal( me = poss_srcVars%varName, &
        &                             val = trim(varname)        )
      me%method(iSource)%updateSrc => eval_source(poss_srcVarPos)%do
      nComponents = poss_srcVars%nComponents%val(poss_srcVarPos)

      if (len_trim(me%varDict%val(iSource)%value) == 0 ) then
        ! In case of the predefined sources the input_varname should just contain
        ! the state variables
        ! We use the global atl_varSys_data here. This is possible because the
        ! first access to atl_varSys_data will happen after it is initialized.
        call tem_varSys_append_derVar(                                 &
          & me             = varSys,                                   &
          & varName        = varname,                                  &
          & operType       = 'predefined_source',                      &
          & nComponents    = nComponents,                              &
          & input_varname  = inputVarNames(1:varSys%nStateVars),       &
          & method_data    = atl_get_new_varSys_data_ptr(varSys_data), &
          & get_point      = get_point,                                &
          & get_element    = get_element,                              &
          & set_params     = set_params,                               &
          & get_params     = get_params,                               &
          & setup_indices  = setup_indices,                            &
          & get_valOfIndex = get_valOfIndex,                           &
          & pos            = addedSourcePosition,                      &
          & wasAdded       = wasAdded                                  )

      else
        ! Here we add the source variable specific space-time-function
        ! as the last input variable.
        ! The state variables in the first slots stay unchanged.
        inputVarNames(varSys%nStateVars+1) = me%varDict%val(iSource)%value

        ! We use the global atl_varSys_data here. This is possible because the
        ! first access to atl_varSys_data will happen after it is initialized.
        call tem_varSys_append_derVar(                                 &
          & me             = varSys,                                   &
          & varName        = varname,                                  &
          & operType       = 'source',                                 &
          & nComponents    = nComponents,                              &
          & input_varname  = inputVarNames,                            &
          & method_data    = atl_get_new_varSys_data_ptr(varSys_data), &
          & get_point      = get_point,                                &
          & get_element    = get_element,                              &
          & set_params     = set_params,                               &
          & get_params     = get_params,                               &
          & setup_indices  = setup_indices,                            &
          & get_valOfIndex = get_valOfIndex,                           &
          & pos            = addedSourcePosition,                      &
          & wasAdded       = wasAdded                                  )

      endif

      if (wasAdded) then

        write(logUnit(5),*) 'Appended variable ' // trim(varname) &
          & // ' to the variable system'
        me%method(iSource)%srcTerm_varPos = addedSourcePosition
        ! Take the name of the variable that provides access to the space time
        ! function for this source term, find it's position in the
        ! varSys%varname array as this is also the position in the
        ! varSys%method array and assign it to the source. Thus we can access
        ! the space time variable immediately when evaluating the source.
        ! it makes sense to do this only for the sources which are not pre-defined
        if (len_trim(me%varDict%val(iSource)%value) > 0 ) then

          me%method(iSource)%data_varPos = positionOfVal( &
            & me  = varSys%varname,                       &
            & val = me%varDict%val(iSource)%value )

          nComp_data = varSys%method%val(me%method(iSource)%data_varPos) &
            &                %nComponents
          if ( nComp_data /= nComponents ) then
            write(logUnit(1),'(a)') 'Error: Appending source variable'
            write(logUnit(1),'(a,i0)') 'nComponent of defined variable: "' &
              & // trim(me%varDict%val(iSource)%value)//'"= ', nComp_data
            write(logUnit(1),'(a,i0)') '/= nComponent of expected variable: "' &
              & // trim(varname)//'"= ', nComponents
            call tem_abort()
          end if
        else
          ! permanent source terms
          me%method(iSource)%data_varPos = -1
          me%method(iSource)%isPermanent = .true.
        end if

      else

        write(logUnit(1),*) 'Error: variable ' // trim(varname)          &
          & // ' is not added to the variable system. Adding it to the ' &
          & // ' list for the next turn'

      end if

    end do

  end subroutine atl_append_newSourceVars
! *******************************************************************************


! *******************************************************************************
  !> Create source elements list for given source variable
  !!
  !! This routine was taken from musubi's mus_source_module. This was
  !! necessary, because we can't move the routine to treelm and work with
  !! type extensions, as there are differences between atl_source_elems_type and
  !! mus_source_elems_type, which prevents us from using a shared and extendable
  !! tem_source_op_type.
  !!
  !! If this routine is changed, please also keep the coresponding routine in
  !! mus_source_module up to date.
  subroutine atl_fill_sourceIndex( source, varSys, nDim, tree, poly_proj_list, &
    &                              mesh_list )
    ! ---------------------------------------------------------------------------
    !> Instance of atl_source_type to be initialized.
    !! This instance will contain the source definitions from lua as well as
    !! the corresponding variables in the global variable system.
    type(atl_source_type), intent(inout) :: source

    !> global variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> equation nDimensions
    integer, intent(in) :: nDim

    !> global treelm mesh
    type( treelmesh_type ), intent(in) :: tree

    !> unique list for projection methods
    type(ply_poly_project_type), target, intent(in) :: poly_proj_list(:)

    !> Mesh list to access the level descriptors
    type(atl_cube_elem_type), intent(in) :: mesh_list(tree%global%minLevel:)
    ! ---------------------------------------------------------------------------
    integer :: iElem, iLevel, iSource, nElems, nQuadPnts
    integer :: data_varPos, minLevel, maxLevel
    type(ply_prj_body_type), pointer :: polyProjBody
    real(kind=rk), allocatable, dimension(:,:) :: refPoints, srcPhysPnts
    integer, allocatable :: idx(:)
    ! ---------------------------------------------------------------------------
    write(logUnit(1),*) 'Fill source index'
    ! return if there are not active source terms
    if ( source%varDict%nVals == 0 ) then
      write(logUnit(1),*) 'No active source terms'
      return
    end if

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    do iSource = 1, source%varDict%nVals
      data_varPos = source%method(iSource)%data_varPos
      if ( data_varPos > 0 ) then
        ! Set params for source variables
        call varSys%method%val(data_varPos)%set_params( &
          & varSys   = varSys,                          &
          & instring = 'isSurface = false'              )
      end if

      ! allocate source elems level array
      allocate(source%method(iSource)%elems(minLevel:maxLevel))
    end do


    write(logUnit(3),*) ' Setup indices for source terms'
    ! run over all levels
    do iLevel = minLevel, maxLevel

      ! First get the amount of points to allocate arrays accordingly
      select case(nDim)
      case (1)
        ! set the correct projection type
        polyProjBody => poly_proj_list(source%poly_proj_pos(iLevel))%body_1d
      case (2)
        polyProjBody => poly_proj_list(source%poly_proj_pos(iLevel))%body_2d
      case (3)
        polyProjBody => poly_proj_list(source%poly_proj_pos(iLevel))%body_3d
      end select

      nQuadPnts = polyProjBody%nQuadPoints
      allocate(refPoints(nQuadPnts,3))
      refPoints = polyProjBody%nodes

      allocate(srcPhysPnts(nQuadPnts,3))
      allocate(idx(nQuadPnts))

      ! total number of fluid elements on this level
      nElems = mesh_list(iLevel)%descriptor%elem%nElems(eT_fluid)

      do iElem = 1, nElems

        !  Get the cheb. phys. coordinate for the element
        call atl_refToPhysCoord(                                   &
          & refPoints  = refPoints,                                &
          & nPoints    = nQuadPnts,                                &
          & baryCoord  = mesh_list(iLevel)%bary_coord( iElem, :),  &
          & elemLength = mesh_list(iLevel)%length,                 &
          & physPoints = srcPhysPnts                               )

        do iSource = 1, source%varDict%nVals
          data_varPos = source%method(iSource)%data_varPos
          ! we do it only for the list of sources defined in the lua files
          if ( data_varPos > 0 ) then
            idx = 0
            call varSys%method%val(data_varPos)%setup_indices( &
              &  varSys     = varSys,                          &
              &  point      = srcPhysPnts,                     &
              &  iLevel     = iLevel,                          &
              &  tree       = tree,                            &
              &  nPnts      = nQuadPnts,                       &
              &  idx        = idx                              )

            ! store index only if idx>0 i.e point is active on data variable
            ! shape
            if ( all(idx > 0) ) then
              call append( me  = source%method(iSource)%elems(iLevel) &
                &                                      %posInTotal,   &
                &          val = iElem                                )
              call append( me  = source%method(iSource)%elems(iLevel)%idx, &
                &          val = idx                                       )
            end if
          end if
        end do ! iSource
      end do ! iElem

      ! set number of elements active per source term.
      ! For predefined or permanent sources, set nFluids on this level
      do iSource = 1, source%varDict%nVals
        data_varPos = source%method(iSource)%data_varPos
        if ( data_varPos > 0 ) then
          source%method(iSource)%elems(iLevel)%nElems               &
            & = source%method(iSource)%elems(iLevel)%posInTotal%nVals
        else
          source%method(iSource)%elems(iLevel)%nElems = nElems
        end if
      end do

      ! Deallocate the allocated arrays before the end of ilevel loop
      deallocate(refPoints)
      deallocate(srcPhysPnts)
      deallocate(idx)

    end do ! iLevel

    do iSource = 1, source%varDict%nVals
      write(logUnit(10),*) ' source: (iLevel , nElem) for srcTerm: ' &
        &                // trim(source%varDict%val(iSource)%key)
      do iLevel = minLevel, maxLevel
        write(logUnit(10),*) '(',iLevel,' , ', source%method(iSource)   &
          &                                    %elems(iLevel)%nElems, ')'
      end do
    end do
  end subroutine atl_fill_sourceIndex
! *******************************************************************************

! ******************************************************************************
  !> This routine converts primitive variables in source terms to convervative
  subroutine atl_source_prim2cons(equation, source, varSys)
    ! --------------------------------------------------------------------------
    type(atl_equations_type), intent(in) :: equation
    type(atl_source_type), intent(in) :: source
    type(tem_varSys_type), intent(in) :: varSys
    ! --------------------------------------------------------------------------
    character(len=labelLen) :: vartype
    integer :: iSrc, iFun, data_varPos
    real(kind=rk) :: cons(1,1,varSys%nScalars)
    type(tem_st_fun_listElem_type), pointer :: stFun => NULL()
    logical :: isSpatialSponge
    ! --------------------------------------------------------------------------
    do iSrc = 1, source%varDict%nVals

      data_varPos = source%method(iSrc)%data_varPos

      select case ( trim(source%varDict%val(iSrc)%key) )
      case ('spongelayer')
        call varSys%method%val(data_varPos)%get_params( &
          & varSys    = varSys,                         &
          & instring  = 'vartype',                      &
          & outstring = vartype                         )
        if (trim(vartype) == 'st_fun') then
          call c_f_pointer(                               &
            & varSys%method%val(data_varPos)%method_data, &
            & stFun                                       )

          isSpatialSponge = .true.
          do iFun = 1,size(stFun%val)
            isSpatialSponge = isSpatialSponge                                 &
              & .and. ( stFun%val(iFun)%spatial%kind(1:11) == 'spongelayer' )
          end do
          if (.not. isSpatialSponge) then
            call tem_abort('Error: spongelayer source variable st_fun is ' &
              &          //'not spatial spongelayer')
          end if

          do iFun = 1,size(stFun%val)
            ! write the primitive values in target state to cons
            cons(1,1,:) = stFun%val(iFun)%spatial%spongePlane%targetstate

            call equation%prim2cons( instate  = cons,     &
              &                      nElems   = 1         )

            ! write it back to targetState
            stFun%val(iFun)%spatial%spongePlane%targetstate = cons(1,1,:)
          enddo
        else
          call tem_abort('Error: spongelayer source variable is not st_fun')
        end if
      end select
    enddo
  end subroutine atl_source_prim2cons
! ******************************************************************************


! ******************************************************************************
  !> summary: subroutine to calculate the RHS of the PDE from the sum of all
  !! source terms
  subroutine atl_update_sourcedata( equation, time, mesh, poly_proj,       &
    &                               currentLevel, state, material, source, &
    &                               scheme                                 )
    ! --------------------------------------------------------------------------
    !> The equation with source term data
    type(atl_equations_type), intent(in) :: equation
    !> current time
    type(tem_time_type), intent(in)             :: time
    !> Current level mesh information
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The projection used for update the source terms for.
    type(ply_poly_project_type), intent(inout)  :: poly_proj
    !> The current Level
    integer, intent(in)                         :: currentLevel
    !> state vector (sources may depend on the state)
    real(kind=rk), intent(in)                   :: state(:,:,:)
    !> The material description.
    type(atl_material_type), intent(inout)      :: material
    !> sources for this level
    type(atl_source_type), intent(inout)        :: source
    !> The scheme you update the source terms for.
    type(atl_scheme_type), intent(in)               :: scheme
    ! --------------------------------------------------------------------------
    integer :: iSource
    ! --------------------------------------------------------------------------

    ! loop over all sources
    do iSource = 1, source%varDict%nVals

      ! init sourcedata with zeros
      source%method(iSource)%val = 0.0_rk

      select case(scheme%scheme)
      case(atl_modg_scheme_prp, atl_modg_2d_scheme_prp)

        ! build the resulting source term ( call the do function pointer)
        call source%method(iSource)%updateSrc(                        &
          & varSys       = equation%varSys,                           &
          & time         = time,                                      &
          & mesh         = mesh,                                      &
          & poly_proj    = poly_proj,                                 &
          & currentLevel = currentLevel,                              &
          & state        = state,                                     &
          & material     = material%material_dat                      &
          &                        %elemMaterialData(atl_ConstMatIdx) &
          &                        %materialDat(1,1,:),               &
          & sourcedata   = source%method(iSource)%val                 )

        !> @todo PV 20150925 The call above gives the material value for the
        !!                   first point of the first element to evaluate the
        !!                   source term. As this is not necessarily the correct
        !!                   value, we thought about handing over the whole
        !!                   material_dat. Out of luck, this is not necessary to
        !!                   get material up and running, so I didn't do it yet.

      case default
        write(logUnit(1),*) 'ERROR in update_sourcedata: not able to ' &
          & // 'evaluate source terms for this scheme type, stopping...'
        call tem_abort()
      end select

    end do

  end subroutine atl_update_sourcedata
! *******************************************************************************

! *******************************************************************************
  subroutine atl_allocate_sourceData( source, nDofs, nComponents )
    ! ------------------------------------------------------------------- !
    type(atl_source_type), intent(inout) :: source
    integer, intent(in) :: nDofs
    integer, intent(in) :: nComponents
    ! ------------------------------------------------------------------- !
    integer :: iSource
    ! ------------------------------------------------------------------- !

    ! Allocate the array for storing the StateData
    do iSource = 1, size(source%method,1)
      ! In case of multi-leveles it should only be allocated once
      if (.not. allocated (source%method(iSource)%val) ) then
        allocate( source%method(iSource)%val(                       &
          &         maxVal(source%method(iSource)%elems(:)%nElems), &
          &         nDofs,                                          &
          &         nComponents )                                   )
      end if
    end do

  end subroutine atl_allocate_sourceData
! *******************************************************************************

! *******************************************************************************
  !> Deallocates the array for storing the sourceData for the currentLevel
  !  after it is no longer needed
  subroutine atl_deallocate_sourceData( source )
    !> Levelwise list of sources
    type(atl_source_type), intent(inout) :: source
    ! --------------------------------------------------------------------------
    integer :: iSource
    ! --------------------------------------------------------------------------

    do iSource = 1, size(source%method)
      deallocate(source%method(iSource)%val)
    end do

  end subroutine atl_deallocate_sourceData
! *******************************************************************************


end module atl_source_module
