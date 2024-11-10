! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2018, 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011-2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2013-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2016-2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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

!> author: Jens Zudrop
!! Initialization of Ateles
!!
!! This module collects the various initialization tasks, that need to be done,
!! to set up the simulation.
module atl_initialize_module
  use, intrinsic :: iso_c_binding,   only: c_loc, C_F_POINTER
  use env_module,                    only: rk, pathLen
  use aotus_module,                  only: flu_state, close_config

  use treelmesh_module,              only: treelmesh_type
  use tem_aux_module,                only: tem_open_distconf
  use tem_logging_module,            only: logUnit
  use tem_variable_module,           only: tem_variable_type,                  &
    &                                      tem_variable_load
  use tem_derived_module,            only: tem_varSys_append_luaVar
  use tem_varSys_module,             only: tem_varSys_solverData_evalElem_type
  use tem_meshInfo_module,           only: tem_varSys_append_meshInfoVar
  use tem_stencil_module,            only: tem_stencilHeader_type,   &
    &                                      init,                     &
    &                                      d3q6_cxDir,               &
    &                                      tem_stencil_map_toTreelmDef
  use tem_precice_module,            only: precice_available, &
    &                                      tem_precice_create
  use tem_operation_var_module,      only: tem_opVar_reduction_transient_init

  use tem_spacetime_fun_module,      only: tem_create_subTree_of_st_funList
  use tem_convergence_module,        only: tem_init_convergence

  use ply_dynArray_project_module,   only: dyn_projectionArray_type, init
  use ply_poly_project_module,       only: ply_poly_project_type, &
    &                                      ply_fill_project_list
  use ply_sampled_tracking_module,   only: ply_sampled_track_init

  use atl_bc_header_module,          only: atl_store_bcVarPos
  use atl_bc_state_module,           only: atl_bc_state_set_fromPoint
  use atl_equation_module,           only: atl_equations_type
  use atl_equation_init_module,      only: atl_init_equation
  use atl_godunovFlux_module,        only: atl_flux_initGodunov
  use atl_container_module,          only: atl_init_elem_container, &
    &                                      atl_element_container_type
  use atl_space_basis,               only: atl_init_spacebasis
  use atl_initial_condition_module,  only: atl_load_initial_condition
  use atl_solver_param_module,       only: atl_solver_param_type
  use atl_restart_module,            only: atl_initRestart, &
    &                                      atl_readRestart
  use atl_source_module,             only: atl_initialize_sources, &
    &                                      atl_allocate_sourcedata
  use atl_source_types_module,       only: atl_init_source_type
  use atl_scheme_module,             only: atl_modg_scheme_prp,    &
    &                                      atl_modg_2d_scheme_prp, &
    &                                      atl_modg_1d_scheme_prp, &
    &                                      atl_scheme_type
  use atl_materialIni_module,        only: atl_init_materialParams
  use atl_materialPrp_module,        only: atl_init_material_type
  use atl_load_project_module,       only: atl_load_projection
  use atl_varSys_module,             only: atl_varSys_solverData_type, &
    &                                      atl_init_varSys_solverData, &
    &                                      atl_set_stFun_getElement,   &
    &                                      atl_create_fortranVar,      &
    &                                      atl_varSys_load_user
  use atl_operator_module,           only: atl_set_opVar_getElement
  use atl_boundary_module,           only: atl_fill_BCIndex
  use atl_penalization_module,       only: atl_init_penalization
  use atl_weights_module,            only: atl_initWeights
  use atl_precice_module,            only: atl_init_precice

  implicit none

  private

  public :: atl_initialize


contains


  ! ****************************************************************************
  !> Routine to initialize the complete Ateles solver.
  !!
  !! It reads the configuration from the given Lua script, and
  !! creates the necessary data structures for the computation.
  !!
  !! For a detailed description see [Initialization](../page/init.html)
  subroutine atl_initialize( params, container, equation, tree, varSys_data, &
    &                        timestepLimit, poly_proj_list                   )
    ! --------------------------------------------------------------------------
    type(atl_solver_param_type), intent(inout) :: params
    !> Container of mesh data.
    !! This includes information about the mesh parts, the scheme, but also the
    !! initial kernel states (i.e. initial conditions).
    type(atl_element_container_type), intent(inout), target :: container

    !> Description on the equation system to solve.
    type(atl_Equations_type), intent(inout), target :: equation

    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(inout) :: tree

    !> Data infomation of the variable System
    type(atl_varSys_solverData_type), intent(inout), target :: varSys_data

    !> timesteplimitation due to the precice timestep
    real(kind=rk), intent(inout), optional :: timestepLimit

    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout), allocatable, target &
      & :: poly_proj_list(:)
    ! -----local variables------------------------------------------------------
    !> unique list for projection methods
    type(dyn_projectionArray_type) :: dynprojectArray
    !> position pointing if projection method used for ic
    integer, allocatable :: ic_prj_pos_list(:)
    ! The number of degrees of freedoms to write to/read from
    ! a tracking or restart file
    integer :: nDofsIO
    integer :: scheme_dim
    integer :: stateShape(3)
    type(tem_stencilHeader_type) :: stencil
    !> This structure takes the information about source terms from the equation
    !! system and stores them to evaluate the present source terms in the
    !! configuration file. After adding the active source terms to the variable
    !! system, this structure is not needed anymore.
    type(atl_init_source_type) :: initSource
    type(atl_init_material_type) :: initMaterial
    ! Name of the file, the configuration is read from.
    character(len=PathLen) :: configFile
    ! --------------------------------------------------------------------------
    type(tem_variable_type), allocatable :: userVars(:)
    type(tem_variable_type), allocatable :: atlVars(:)
    ! --------------------------------------------------------------------------
    integer, allocatable :: vError(:)
    type(tem_varSys_solverData_evalElem_type) :: solverData_evalElem
    type(flu_State) :: eqn_conf
    ! --------------------------------------------------------------------------

    solverData_evalElem%solver_bundle = c_loc(varSys_data)
    solverData_evalElem%stFun_setter => atl_set_stfun_getElement
    solverData_evalElem%opVar_setter => atl_set_opVar_getElement

    call atl_bc_state_set_fromPoint(solverData_evalElem)

    configFile = trim(params%general%solver%configFile)
    ! Create the appropriate D3Q6 stencil
    call init( me     = stencil,   &
      &        QQN    = 6,         &
      &        QQ     = 6,         &
      &        useAll = .false.,   &
      &        nDims  = 3,         &
      &        label  = 'd3q6',    &
      &        cxDir  = d3q6_cxDir )
    call tem_stencil_map_toTreelmDef(stencil)

    ! **************************************************************************
    ! EQUATION
    !
    ! Initialize the equation description
    call atl_init_equation( equation     = equation,                      &
      &                     conf         = params%general%solver%conf(1), &
      &                     varSys_Data  = varSys_data,                   &
      &                     initSource   = initSource,                    &
      &                     initMaterial = initMaterial                   )

    if ( (trim(equation%eq_kind) == 'unknown') ) then

      ! If there was no equation system provided in the configuration,
      ! try to get it from the restart file, when available.
      if ( params%general%restart%controller%readRestart ) then

        call tem_open_distconf(                                              &
          & L        = eqn_conf,                                             &
          & filename = trim(params%general%restart%controller%readFileName), &
          & proc     = params%general%proc                                   )

        call atl_init_equation( equation     = equation,            &
          &                     conf         = eqn_conf,            &
          &                     varSys_Data  = varSys_data,         &
          &                     initSource   = initSource,          &
          &                     initMaterial = initMaterial         )

        call close_config(eqn_conf)

        if ( (trim(equation%eq_kind) /= 'unknown') ) then
          write(logunit(1),*) 'Using the equation table from the restart!'
        end if

      end if

    else
      write(logunit(1),*) 'Using the equation table from the configuration.'
    end if

    if ( (trim(equation%eq_kind) == 'unknown') ) &
      &  write(logunit(1),*) 'WARNING: no equation table found !!'


    select case(trim(equation%eq_kind))
    case( 'euler', 'euler_2d', 'euler_1d',                      &
      &   'navier_stokes', 'navier_stokes_2d',                  &
      &   'filtered_navier_stokes', 'filtered_navier_stokes_2d' )
      call atl_flux_initGodunov(equation%euler%isen_coef)
    end select


    ! Init precice if pecice is used
    if (precice_available) then
      write(logUnit(1),*) 'Create precice environment'
      call tem_precice_create( rank      = params%general%proc%rank,     &
        &                      comm_size = params%general%proc%comm_size )
    end if

    ! **************************************************************************
    ! ELEMENT CONTAINER
    !
    ! Initialize the simulation container (including the scheme)
    call atl_init_elem_container(                                        &
      & me              = container,                                     &
      & equation        = equation,                                      &
      & conf            = params%general%solver%conf(1),                 &
      & tree            = tree,                                          &
      & time            = params%general%simControl%now,                 &
      & readRestart     = params%general%restart%controller%readRestart, &
      & proc            = params%general%proc,                           &
      & commPattern     = params%general%commPattern,                    &
      & boundary        = params%boundary                                )


    ! **************************************************************************
    ! PROJECTION
    !
    ! Now the general and all individuall projection tables are loaded and
    ! added to the DA for projection
    call atl_load_projection(                                      &
      & minLevel           = tree%global%minlevel,                 &
      & maxLevel           = tree%global%maxlevel,                 &
      & equation           = equation,                             &
      & poly_proj_pos      = container%cubes%poly_proj_pos,        &
      & dynprojectArray    = dynProjectArray,                      &
      & conf               = params%general%solver%conf(1),        &
      & scheme_list        = container%cubes%scheme_list,          &
      & source_projPos     = container%cubes%source%poly_proj_pos, &
      & boundary_list      = container%cubes%boundary_list,        &
      & boundary_list_stab = container%cubes%boundary_stab_list,   &
      & ic_pos_list        = ic_prj_pos_list,                      &
      & material_list      = container%cubes%material_list         )

    ! After loading all information needed for individuall projection method
    ! we can finally fill up the projection body
    !Figure out the dimension of the scheme to reduce overhead
    select case(container%cubes%scheme_list(tree%global%minLevel)%scheme)
      case (atl_modg_scheme_prp)
         scheme_dim = 3
      case (atl_modg_2d_scheme_prp)
         scheme_dim = 2
      case (atl_modg_1d_scheme_prp)
         scheme_dim = 1
    end select

    call ply_fill_project_list( proj_list        = poly_proj_list,  &
      &                         dyn_projectArray = dynProjectArray, &
      &                         scheme_dim       = scheme_dim       )


    call atl_alloc_temp( minLevel    = tree%global%minlevel,       &
      &                  maxLevel    = tree%global%maxlevel,       &
      &                  equation    = equation,                   &
      &                  proj_list   = poly_proj_list,             &
      &                  scheme_list = container%cubes%scheme_list )

    ! **************************************************************************
    ! VARIABLE SYSTEM
    !
    ! Load the user defined variables from the configuration file and store them
    ! in a temporary list.
    call tem_variable_load( me             = userVars,                      &
      &                     conf           = params%general%solver%conf(1), &
      &                     load_solvervar = atl_varSys_load_user,          &
      &                     vError         = vError                         )

    ! take the temporary list with the lua variables and try to add them to the
    ! variable system.
    call tem_varSys_append_luaVar(                      &
      & luaVar                   = userVars,            &
      & varSys                   = equation%varSys,     &
      & st_funList               = equation%stFunList,  &
      & solverData_evalElem      = solverData_evalElem  )

    ! Create some default ateles variables.
    ! For example: Lot of boundary condition variables requires constant zero as
    ! spacetime-function so create a variable "zero_const" and refer to this
    ! name in boundary condition variables
    call atl_create_fortranVar( me = atlVars )

    ! Append hard coded ateles spacetime function variable to the
    ! variable system.
    ! @TODO KM: Do no add atlVars into stFunList
    call tem_varSys_append_luaVar(                      &
      & luaVar                   = atlVars,             &
      & varSys                   = equation%varSys,     &
      & st_funList               = equation%stFunList,  &
      & solverData_evalElem      = solverData_evalElem  )

    ! Append variables with information on the mesh.
    call tem_varSys_append_meshInfoVar(equation%varSys)

    ! Create a subtree for each space time function based on the defined shape.
    ! This subtree can be used to create a list of elements that are affected
    ! by the space time function. This simplifies the evaluation of the
    ! particular space time function as there is no need to loop over all
    ! elements and check whether they are affected by the space time function.
    ! Instead, the space time function has a list of affected elements to loop
    ! over.
    call tem_create_subTree_of_st_funList( &
      & me      = equation%stFunList,      &
      & tree    = tree,                    &
      & bc_prop = params%boundary          )

    ! Initialize the global instance of atl_varSys_data
    call atl_init_varSys_solverdata(                      &
      & this           = varSys_data,                     &
      & equation       = equation,                        &
      & tree           = tree,                            &
      & schemeList     = container%cubes%scheme_list,     &
      & statedataList  = container%cubes%stateData_list,  &
      & kerneldataList = container%cubes%kernelData_list, &
      & meshList       = container%cubes%mesh_list,       &
      & polyProjList   = poly_proj_list,                  &
      & levelPointer   = container%cubes%levelPointer,    &
      & polyProjPos    = container%cubes%poly_proj_pos    )

    ! **************************************************************************
    ! BOUNDARY CODITIONS
    !
    ! Store position of a variable name defined for each boundary
    ! variable in tem_bc_state_type
    call atl_store_bcVarPos( bc     = container%cubes%bc, &
      &                      varSys = equation%varSys     )

    ! For BC we use the seuptIndices and getValOfIndex routines to access
    ! variables in the varSys. This need to be done after bc are initialized
    ! in init_cube_container and after the projection is
    ! initialized
    call atl_fill_BCIndex( tree           = tree,                          &
      &                    bc             = container%cubes%bc,            &
      &                    boundary_list  = container%cubes%boundary_list, &
      &                    nDim           = equation%nDimensions,          &
      &                    varSys         = equation%varSys,               &
      &                    poly_proj_list = poly_proj_list,                &
      &                    mesh_list      = container%cubes%mesh_list      )

    ! **************************************************************************
    ! SOURCES
    !
    call atl_initialize_sources(                        &
      & source         = container%cubes%source,        &
      & initSource     = initSource,                    &
      & conf           = params%general%solver%conf(1), &
      & equation       = equation,                      &
      & poly_proj_list = poly_proj_list,                &
      & mesh_list      = container%cubes%mesh_list,     &
      & tree           = tree,                          &
      & varSys_data    = varSys_data                    )

    stateShape = shape( container%cubes                                    &
      &                          %statedata_list(container%cubes%minLevel) &
      &                          %state                                    )

    call atl_allocate_sourceData( source      = container%cubes%source, &
      &                           nDofs       = stateShape(2),          &
      &                           nComponents = stateShape(3)           )

    ! **************************************************************************
    ! LB weights measurement setup
    !
    call atl_initWeights( tree )

    ! **************************************************************************
    ! MATERIAL
    !
    call atl_init_materialParams(                             &
      & equation             = equation,                      &
      & tree                 = tree,                          &
      & varSys_data          = varSys_data,                   &
      & material_list        = container%cubes%material_list, &
      & mesh_list            = container%cubes%mesh_list,     &
      & scheme_list          = container%cubes%scheme_list,   &
      & boundary_list        = container%cubes%boundary_list, &
      & time                 = params%general%simControl%now, &
      & conf                 = params%general%solver%conf(1), &
      & proc                 = params%general%proc,           &
      & commPattern          = params%general%commPattern,    &
      & poly_proj_list       = poly_proj_list,                &
      & levelpointer         = container%cubes%levelpointer,  &
      & initMaterial         = initMaterial                   )

    ! **************************************************************************
    ! PENALIZATION
    !
    call atl_init_penalization(                                        &
      & tree                  = tree,                                  &
      & penalizationdata_list = container%cubes%penalizationdata_list, &
      & scheme_list           = container%cubes%scheme_list,           &
      & equation              = equation,                              &
      & mesh_list             = container%cubes%mesh_list              )
!PV!      ,             &
!PV!      & penalization_list     = container%cubes%penalization_list,     &
!PV!      & conf                  = params%general%solver%conf(1),         &
!PV!      & levelpointer          = container%cubes%levelpointer           )


    ! The number of degrees of freedoms to write out in tracking or
    ! restart
    nDofsIO = maxval(container%cubes%scheme_list(:)%nDofs)

    ! Fill the spatial polynomial basis, that is needed for integrations
    ! and point evaluations.
    call atl_init_spacebasis( scheme_list = container%cubes%scheme_list, &
      &                       minlevel    = tree%global%minlevel,        &
      &                       maxlevel    = tree%global%maxlevel         )

    ! **************************************************************************
    ! RESTART
    !
    call atl_initRestart( restart    = params%general%restart, &
      &                   equation   = equation,               &
      &                   solver     = params%general%solver,  &
      &                   mesh       = container%cubes,        &
      &                   tree       = tree,                   &
      &                   nDofsIO    = nDofsIO,                &
      &                   scheme_dim = scheme_dim              )

    ! **************************************************************************
    ! INITIAL CONDITION
    !
    if (params%general%restart%controller%readRestart) then
      ! Actually read the restart data as initial condition
      call atl_readRestart( mesh           = container%cubes,        &
        &                   restart        = params%general%restart, &
        &                   tree           = tree,                   &
        &                   proc           = params%general%proc,    &
        &                   equation       = equation                )
    else
      ! No restart, get the initial condition from the configuration
      call atl_load_initial_condition(                    &
        & cube_container = container%cubes,               &
        & equation       = equation,                      &
        & conf           = params%general%solver%conf(1), &
        & tree           = tree,                          &
        & luaFile        = trim(configFile),              &
        & prj_pos        = ic_prj_pos_list,               &
        & poly_proj_list = poly_proj_list                 )
    end if


    ! **************************************************************************
    ! PRECICE
    !
    !> COUPLING USING EQUIDISTANT AND NON-EQUIDISTANT POINTS
    !! Reminder: For the equidistant and non-equidistant points the call for
    !! the setup_indices is included in the
    !! atl_write_equiPoints/atl_writenonequiPoints routine.
    if (precice_available) then
      call atl_init_precice(                         &
        & params        = params,                    &
        & equation      = equation,                  &
        & tree          = tree,                      &
        & mesh_list     = container%cubes%mesh_list, &
        & timestepLimit = timestepLimit              )
    end if

    ! **************************************************************************
    ! TRACKING
    !
    call ply_sampled_track_init( me      = params%plySampleTrack, &
      &                          mesh    = tree,                  &
      &                          solver  = params%general%solver, &
      &                          varSys  = equation%varSys,       &
      &                          bc      = params%boundary,       &
      &                          stencil = stencil,               &
      &                          proc    = params%general%proc,   &
      &                          ndofs   = ndofsIO,               &
      &                          ndims   = equation%nDimensions   )

    ! **************************************************************************
    ! CONVERGENCE
    !
    if ( params%general%simControl%abortCriteria%steady_state ) then
      call tem_init_convergence( me       = params%general%simControl    &
        &                                                 %abortCriteria &
        &                                                 %convergence,  &
        &                        tree     = tree,                        &
        &                        bc_prop  = params%boundary,             &
        &                        globProc = params%general%proc,         &
        &                        varSys   = equation%varSys,             &
        &                        ndofs    = ndofsIO                      )
    end if

    ! Initialize reduction_transient operation variables
    call tem_opVar_reduction_transient_init(           &
      & varSys         = equation%varSys,              &
      & tree           = tree,                         &
      & redTransVarMap = equation%redTransVarmap,      &
      & nDofs          = nDofsIO,                      &
      & time           = params%general%simControl%now )

  end subroutine atl_initialize
  ! ****************************************************************************


  ! ****************************************************************************
  !> Allocate temporary arrays.
  subroutine atl_alloc_temp( minLevel, maxLevel, proj_list, &
    &                        scheme_list, equation   )
  ! ----------------------------------------------------------------------------
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    type(ply_poly_project_type), intent(inout), allocatable :: proj_list(:)
    type(atl_equations_type), intent(in) :: equation
    type(atl_scheme_type), intent(inout) :: scheme_list(minlevel:maxLevel)
    ! --------------------------------------------------------------------------

    ! allocate the temporary arrays required for the calculation of fluxes
    ! Note : They are preallocated so that the allocation is not needed at the
    ! time they are used ( in the thread parallel region)
    select case(scheme_list(minLevel)%scheme)
      case (atl_modg_scheme_prp)

        allocate( scheme_list(minlevel)%temp_over(      &
          & maxVal(proj_list(:)%body_3d%oversamp_dofs), &
          & equation%temp%nScal,                        &
          & equation%temp%oversamp )                    )
        allocate( scheme_list(minlevel)%temp_nodal(   &
          & maxval(proj_list(:)%body_3d%nquadpoints), &
          & equation%temp%nScal,                      &
          & equation%temp%nodal )                     )
        allocate( scheme_list(minlevel)%temp_modal( &
          & maxVal(proj_list(:)%body_3d%ndofs),     &
          & equation%temp%nScal,                    &
          & equation%temp%modal )                   )

        scheme_list(minlevel)%temp_over = 0.0_rk
        scheme_list(minlevel)%temp_nodal = 0.0_rk
        scheme_list(minlevel)%temp_modal = 0.0_rk

      case (atl_modg_2d_scheme_prp)

        allocate( scheme_list(minlevel)%temp_over(      &
          & maxVal(proj_list(:)%body_2d%oversamp_dofs), &
          & equation%temp%nScal,                        &
          & equation%temp%oversamp )                    )
        allocate( scheme_list(minlevel)%temp_nodal(    &
          & maxval(proj_list(:)%body_2d%nquadpoints),  &
          & equation%temp%nScal,                       &
          & equation%temp%nodal )                      )
        allocate( scheme_list(minlevel)%temp_modal( &
          & maxVal(proj_list(:)%body_2d%ndofs),     &
          & equation%temp%nScal,                    &
          & equation%temp%modal )                   )

        scheme_list(minlevel)%temp_over = 0.0_rk
        scheme_list(minlevel)%temp_nodal = 0.0_rk
        scheme_list(minlevel)%temp_modal = 0.0_rk

      case (atl_modg_1d_scheme_prp)

        allocate( scheme_list(minlevel)%temp_over(       &
          & maxVal(proj_list(:)%body_1d%oversamp_dofs ), &
          & equation%temp%nScal,                         &
          & equation%temp%oversamp )                     )
        allocate( scheme_list(minlevel)%temp_nodal(   &
          & maxval(proj_list(:)%body_1d%nquadpoints), &
          & equation%temp%nScal,                      &
          & equation%temp%nodal )                     )
        allocate( scheme_list(minlevel)%temp_modal( &
          & maxVal(proj_list(:)%body_1d%ndofs),     &
          & equation%temp%nScal,                    &
          & equation%temp%modal )                   )

        scheme_list(minlevel)%temp_over = 0.0_rk
        scheme_list(minlevel)%temp_nodal = 0.0_rk
        scheme_list(minlevel)%temp_modal = 0.0_rk

    end select

  end subroutine atl_alloc_temp
  ! ****************************************************************************

end module atl_initialize_module
