! Copyright (c) 2013, 2015-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013, 2015-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> Helper routines for the euler equation system.
module atl_eqn_euler_hlp_module
  use env_module,                     only: rk, labelLen
  use aotus_module,                   only: flu_State, aot_get_val

  use tem_aux_module,                 only: tem_abort
  use tem_bc_module,                  only: tem_bc_state_type
  use tem_logging_module,             only: logUnit
  use tem_tools_module,               only: upper_to_lower
  use tem_stringKeyValuePair_module,  only: tem_stringKeyValuePair_type
  use tem_stringKeyValuePair_module,  only: grw_stringKeyValuePairArray_type, &
    &                                       init, truncate, append

  use ply_oversample_module,          only: ply_convert2oversample, &
    &                                       ply_convertFromoversample
  use ply_poly_project_module,        only: ply_poly_project_type, &
    &                                       ply_poly_project_m2n,  &
    &                                       ply_poly_project_n2m

  use atl_equation_module,            only: atl_equations_type,     &
    &                                       atl_eqn_var_trafo_type
  use atl_bc_state_module,            only: atl_load_bc_state
  use atl_eqn_euler_module,           only: atl_load_Euler, &
    &                                       atl_euler_type
  use atl_eqn_euler_derive_module,    only: atl_eqn_euler_prim2cons,       &
    &                                       atl_eqn_euler_cons2primTemp,   &
    &                                       atl_eqn_euler_primTemp2Cons,   &
    &                                       atl_eqn_euler_cons2prim,       &
    &                                       atl_eqn_euler_prim2cons_elems, &
    &                                       atl_eqn_euler_cons2prim_elems
  use atl_eqn_euler_2d_derive_module, only: atl_eqn_euler_2d_prim2cons, &
    &                                       atl_eqn_euler_2d_cons2prim, &
    &                                       atl_eqn_euler_2d_cons2primTemp,   &
    &                                       atl_eqn_euler_2d_primTemp2Cons,   &
    &                                       atl_eqn_euler_2d_prim2cons_elems, &
    &                                       atl_eqn_euler_2d_cons2prim_elems
  use atl_eqn_euler_1d_derive_module, only: atl_eqn_euler_1d_prim2cons, &
    &                                       atl_eqn_euler_1d_cons2prim, &
    &                                       atl_eqn_euler_1d_prim2cons_elems, &
    &                                       atl_eqn_euler_1d_cons2prim_elems
  use atl_eqn_euler_var_module,       only: atl_init_euler_vars,        &
    &                                       atl_init_euler_sourceTerms, &
    &                                       atl_init_euler_material
  use atl_eqn_euler_2d_var_module,    only: atl_init_euler_2d_vars,       &
    &                                       atl_init_euler_2d_sourceTerms
  use atl_eqn_euler_1d_var_module,    only: atl_init_euler_1d_vars,       &
    &                                       atl_init_euler_1d_sourceTerms
  use atl_source_types_module,        only: atl_init_source_type
  use atl_materialPrp_module,         only: atl_init_material_type, &
    &                                       atl_material_type
  use atl_varSys_module,              only: atl_varSys_solverData_type

  use atl_godunovFlux_module,         only: atl_GodunovEuler,   &
    &                                       atl_GodunovEuler2D, &
    &                                       atl_GodunovEuler1D
  use atl_hlleFlux_module,            only: atl_HLLEuler, atl_HLLEuler2D, &
    &                                       atl_HLLEuler1D
  use atl_laxFriedrichFlux_module,    only: atl_laxFriedEuler
  use atl_laxFriedrichFlux_2d_module, only: atl_laxFriedEuler_2d
  use atl_laxFriedrichFlux_1d_module, only: atl_laxFriedEuler_1d

  implicit none

  private

  public :: atl_eqn_euler_load_bc
  public :: atl_eqn_euler_init
  public :: atl_eqn_euler_implicit_pen
  public :: atl_getEulerFluxes
  public :: atl_getEulerLinInd


contains


  ! ****************************************************************************
  !> Initialization of the Euler equation.
  !!
  !! This routine sets up the necessary infrastructure for the Euler equations.
  !! It reads the configuration from the given script in conf under the table
  !! provided in thandle and sets function pointers and variables accordingly.
  subroutine atl_eqn_euler_init( conf, thandle, equation, nDimensions, &
    &                            initSource, initMaterial, varSys_data )
    ! --------------------------------------------------------------------------
    !> Handle to the Lua configuration
    type(flu_State), intent(in) :: conf

    !> Handle to the equation table in the Lua script given in conf.
    integer, intent(in) :: thandle

    !> Equation system to set with this routine.
    type(atl_equations_type), intent(inout) :: equation

    !> Number of spatial dimensions, the Euler equations should live on.
    !!
    !! Has to be 1, 2 or 3.
    integer, intent(in) :: nDimensions

    !> Type to be filled with the possible source variables for the equation
    !! system. These source variables are later on used to extract the
    !! corresponding information from the configuration file.
    type(atl_init_source_type), intent(inout) :: initSource

    !> Type to be filled with the possible material variables for the equation
    !! system. These material variables are later on used to extract the
    !! corresponding information from the configuration file.
    type(atl_init_material_type), intent(inout) :: initMaterial

    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), intent(inout) :: varSys_data
    ! --------------------------------------------------------------------------

    equation%isNonlinear = .true.
    equation%nDerivatives = 0
    equation%nDimensions = nDimensions
    ! timestep is change with time since it is nonlinear
    equation%adaptive_timestep = .true.

    select case(nDimensions)
    case(1)
      equation%load_bc => atl_eqn_euler_load_bc

      equation%cons2prim => atl_eqn_euler_1d_cons2prim_elems
      equation%prim2cons => atl_eqn_euler_1d_prim2cons_elems

      call atl_init_euler_1d_vars( equation   = equation,   &
        &                          solverData = varSys_data )

      call atl_init_euler_1d_sourceTerms( initSource%poss_srcVars, &
        &                                 initSource%eval_source   )

    case(2)
      equation%load_bc => atl_eqn_euler_load_bc

      equation%cons2prim => atl_eqn_euler_2d_cons2prim_elems
      equation%prim2cons => atl_eqn_euler_2d_prim2cons_elems

      call atl_init_euler_2d_vars( equation   = equation,   &
        &                          solverData = varSys_data )

      call atl_init_euler_2d_sourceTerms( initSource%poss_srcVars, &
        &                                 initSource%eval_source   )
    case(3)
      equation%load_bc => atl_eqn_euler_load_bc

      equation%cons2prim => atl_eqn_euler_cons2prim_elems
      equation%prim2cons => atl_eqn_euler_prim2cons_elems

      call atl_init_euler_vars( equation   = equation,   &
        &                       solverData = varSys_data )

      call atl_init_euler_sourceTerms( initSource%poss_srcVars, &
        &                              initSource%eval_source   )

    end select

    call atl_load_euler( euler        = equation%Euler,  &
      &                  conf         = conf,            &
      &                  eq_table     = thandle          )

    ! Set the flag, that we require the computation of deviations, if
    ! the adaptive linearization is active.
    equation%requiresDeviation = (equation%euler%linear_limit > 0.0_rk)

    call atl_init_euler_material(                     &
      & possVars    = initMaterial%poss_materialVars, &
      & nDimensions = nDimensions                     )

    ! Getting the numerical flux !
    call atl_getEulerFluxes(euler      = equation%Euler, &
      &                     conf       = conf,           &
      &                     eqn_handle = thandle,        &
      &                     eqn_dim    = nDimensions     )

    ! Getting the indicator to use in linearization !
    call atl_getEulerLinInd(euler      = equation%Euler, &
      &                     conf       = conf,           &
      &                     eqn_handle = thandle,        &
      &                     eqn_dim    = nDimensions     )

  end subroutine atl_eqn_euler_init
  ! ****************************************************************************


  ! ****************************************************************************
  !> Reading boundary conditions for the euler equations.
  !!
  !! Need to set 5 bc_states here, typically the primitive variables.
  !! Vectorial quantities are described either by the normal component and
  !! a tangential definition that has to be the same in all directions,
  !! or in the universal coordinate system.
  !! The normal is defined as pointing inwards.
  !! Internally the tangential definition is duplicated to get the same size
  !! for vectorial quantities irregardless of the coordinate system it is
  !! defined in.
  !!
  !! This routine has to conform to the interface definition
  !! atl_equation_module#eqn_load_bc.
  subroutine atl_eqn_euler_load_bc( equation,                              &
    &                               bc_state, bc_state_gradient,           &
    &                               bc_varDict, bc_varDict_gradient,       &
    &                               bc_normal_vec, bc_normal_vec_gradient, &
    &                               bc_trafo, bc_trafo_gradient,           &
    &                               bc_label, bc_kind, thandle, conf       )
    ! --------------------------------------------------------------------------
    class(atl_equations_type), intent(inout) :: equation
    type(tem_bc_state_type), allocatable, intent(out) :: bc_state(:)
    type(tem_bc_state_type), allocatable, intent(out) :: bc_state_gradient(:)
    !> Dictionary of boundary variables in bc_state
    type(grw_stringKeyValuePairArray_type), intent(out) :: bc_varDict
    !> Dictionary of boundary variables in bc_state_gradient
    type(grw_stringKeyValuePairArray_type), intent(out) :: bc_varDict_gradient
    logical, intent(out) :: bc_normal_vec
    logical, intent(out) :: bc_normal_vec_gradient
    type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo
    type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo_gradient
    character(len=*), intent(in) :: bc_label
    character(len=*), intent(in) :: bc_kind
    integer, intent(in) :: thandle
    type(flu_State) :: conf
    ! --------------------------------------------------------------------------
    integer :: nDims
    integer :: pIndex
    type(tem_stringKeyValuePair_type) :: kvp
    ! --------------------------------------------------------------------------

    nDims = equation%nDimensions
    pIndex = equation%varSys%nScalars
    allocate(bc_state(pIndex))
    allocate(bc_state_gradient(0))
    bc_normal_vec_gradient = .false.

!!VK    allocate(bc_normal_vec_gradient(2))
!!VK    allocate(bc_trafo_gradient(2))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )
    ! Constant zero variable for non-configurable boundary variable
    kvp%value = 'zero_const'

    ! By default we set the function pointer for a conversion,
    ! even if the boundary condition does not use them.
    select case(nDims)
    case(1)
      bc_trafo%from => atl_eqn_euler_1d_prim2cons
      bc_trafo%to => atl_eqn_euler_1d_cons2prim
    case(2)
      bc_trafo%from => atl_eqn_euler_2d_prim2cons
      bc_trafo%to => atl_eqn_euler_2d_cons2prim
    case(3)
      bc_trafo%from => atl_eqn_euler_prim2cons
      bc_trafo%to => atl_eqn_euler_cons2prim
    end select

    select case(bc_kind)
    case('slipwall', 'wall')
      ! This boundary is given in primitive variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .false.
      bc_normal_vec = .true.
      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'dirichlet'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

      if (nDims > 1) then
        ! Extrapolate v_tangential_1
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'neumann'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        if (nDims > 2) then
          ! Extrapolate v_tangential_2
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'neumann'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )
        end if
      end if

      ! Extrapolate pressure
      bc_state(pIndex)%state_name = 'pressure'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

    case('isothermal_wall')
      ! Use a non-slip boundary
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2primTemp
        bc_trafo%from => atl_eqn_euler_2d_primTemp2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2primTemp
        bc_trafo%from => atl_eqn_euler_primTemp2cons
      end select
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )
      ! Prescribe v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'dirichlet'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )
      if (nDims > 1) then
        ! Impose tangential velocity
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'dirichlet'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )
        if (nDims > 2) then
          ! Impose tangential velocity
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'dirichlet'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )
        endif
      endif
      ! Prescribe temperature
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'temperature',    &
        &                     style       = 'dirichlet',      &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varDict     = bc_varDict,       &
        &                     varSys      = equation%varSys   )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition isothermal_wall you have to set'
        write(logUnit(1),*) 'the Temperature, you did not set the value for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if
    case('primitives')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .false.
      bc_normal_vec = .false.

      call atl_load_bc_state( bc          = bc_state(1),     &
        &                     state_name  = 'density',       &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      call atl_load_bc_state( bc          = bc_state(2),     &
        &                     state_name  = 'velocityX',     &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )
      if (nDims > 1) then
        call atl_load_bc_state( bc          = bc_state(3),     &
          &                     state_name  = 'velocityY',     &
          &                     conf        = conf,            &
          &                     bc_handle   = thandle,         &
          &                     varSys      = equation%varSys, &
          &                     varDict     = bc_varDict       )
        if (nDims > 2) then
          call atl_load_bc_state( bc          = bc_state(4),     &
            &                     state_name  = 'velocityZ',     &
            &                     conf        = conf,            &
            &                     bc_handle   = thandle,         &
            &                     varSys      = equation%varSys, &
            &                     varDict     = bc_varDict       )
        end if
      end if
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varSys      = equation%varSys,  &
        &                     varDict     = bc_varDict        )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition primitives you have to set'
        write(logUnit(1),*) 'all primitive variables (density, '
        write(logUnit(1),*) 'velocityX, velocityY, velocityZ'
        write(logUnit(1),*) 'and pressure) this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('conservatives')
      ! Everything is in conservative, so we do not need a transformation.
      bc_trafo%identity = .true.
      bc_normal_vec = .false.

      call atl_load_bc_state( bc          = bc_state(1),      &
        &                     state_name  = 'density',        &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varSys      = equation%varSys,  &
        &                     varDict     = bc_varDict        )

      call atl_load_bc_state( bc          = bc_state(2),      &
        &                     state_name  = 'momentumX',      &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varSys      = equation%varSys,  &
        &                     varDict     = bc_varDict        )
      if (nDims > 1) then
        call atl_load_bc_state( bc          = bc_state(3),    &
          &                     state_name  = 'momentumY',    &
          &                     conf        = conf,           &
          &                     bc_handle   = thandle,        &
        &                     varSys      = equation%varSys,  &
          &                     varDict     = bc_varDict      )
        if (nDims > 2) then
          call atl_load_bc_state( bc          = bc_state(4),  &
            &                     state_name  = 'momentumZ',  &
            &                     conf        = conf,         &
            &                     bc_handle   = thandle,      &
        &                     varSys      = equation%varSys,  &
            &                     varDict     = bc_varDict   )
        end if
      end if
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'energy',         &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varSys      = equation%varSys,  &
        &                     varDict     = bc_varDict        )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition conservatives you have to'
        write(logUnit(1),*) 'set all conservative variables (density, '
        write(logUnit(1),*) 'momentumX, momentumY, momentumZ and energy) '
        write(logUnit(1),*) 'this set is not'
        write(logUnit(1),*) 'complete for ' // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('inflow')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .false.
      bc_normal_vec = .false.

      ! Impose denisty
      call atl_load_bc_state( bc          = bc_state(1),     &
        &                     state_name  = 'density',       &
        &                     style       = 'dirichlet',     &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      ! Impose x velocity
      call atl_load_bc_state( bc          = bc_state(2),     &
        &                     state_name  = 'velocityX',     &
        &                     style       = 'dirichlet',     &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      if (nDims > 1) then
        ! Impose y velocity
        call atl_load_bc_state( bc          = bc_state(3),     &
          &                     state_name  = 'velocityY',     &
          &                     style       = 'dirichlet',     &
          &                     conf        = conf,            &
          &                     bc_handle   = thandle,         &
          &                     varSys      = equation%varSys, &
          &                     varDict     = bc_varDict       )

        if (nDims > 2) then
          ! Impose z velocity
          call atl_load_bc_state( bc          = bc_state(4),     &
            &                     state_name  = 'velocityZ',     &
            &                     style       = 'dirichlet',     &
            &                     conf        = conf,            &
            &                     bc_handle   = thandle,         &
            &                     varSys      = equation%varSys, &
            &                     varDict     = bc_varDict       )
        end if
      end if

      ! Extrapolate pressure
      bc_state(pIndex)%state_name = 'pressure'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition inflow you have to set the'
        write(logUnit(1),*) 'primitive variables density, '
        write(logUnit(1),*) 'velocityX, velocityY, velocityZ'
        write(logUnit(1),*) 'this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('inflow_normal')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! Impose density at inlet
      call atl_load_bc_state( bc          = bc_state(1),     &
        &                     state_name  = 'density',       &
        &                     style       = 'dirichlet',     &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      ! Impose normal velocity
      call atl_load_bc_state( bc          = bc_state(2),     &
        &                     state_name  = 'v_norm',        &
        &                     style       = 'dirichlet',     &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      if (nDims > 1) then
        ! Impose tangential velocity to zero
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'dirichlet'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        if (nDims > 2) then
          ! Copy the tangential definition into the second tangential direction
          ! (we only support normal boundary definitions where this is valid)
          ! Rename the state_name, the first tangential component is always
          ! v_tan and to distinguish the second (superfluos) one we use v_tan2
          ! here.
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'dirichlet'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )
        end if
      end if

      ! Extrapolate pressure
      bc_state(pIndex)%state_name = 'pressure'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition inflow_normal you have to'
        write(logUnit(1),*) 'set the primitive variables density, v_norm and'
        write(logUnit(1),*) 'v_tan this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('supersonic_inflow_normal')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! Impose density at inlet
      call atl_load_bc_state( bc          = bc_state(1),     &
        &                     state_name  = 'density',       &
        &                     style       = 'dirichlet',     &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      ! Impose normal velocity
      call atl_load_bc_state( bc          = bc_state(2),     &
        &                     state_name  = 'v_norm',        &
        &                     style       = 'dirichlet',     &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      if (nDims > 1) then
        ! Impose tangential velocity
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'dirichlet'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        if (nDims > 2) then
          ! Copy the tangential definition into the second tangential direction
          ! (we only support normal boundary definitions where this is valid)
          ! Rename the state_name, the first tangential component is always
          ! v_tan and to distinguish the second (superfluos) one we use v_tan2
          ! here.
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'dirichlet'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )
        end if
      end if

      ! Impose pressure
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     style       = 'dirichlet',      &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varSys      = equation%varSys,  &
        &                     varDict     = bc_varDict        )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition supersonic_inflow_normal'
        write(logUnit(1),*) 'you have to set the primitive variables density,'
        write(logUnit(1),*) 'v_norm, v_tan and pressure.'
        write(logUnit(1),*) 'This set is not complete for ' &
          &                 // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('outflow')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Extrapolate v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'neumann'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

      if (nDims > 1) then
        ! Extrapolate v_tangential_1
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'neumann'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        if (nDims > 2) then
          ! Extrapolate v_tangential_2
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'neumann'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )
        end if
      end if

      ! Impose pressure
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     style       = 'dirichlet',      &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varSys      = equation%varSys,  &
        &                     varDict     = bc_varDict        )

      if (.not. bc_state(pIndex)%isDefined) then
        write(logUnit(1),*) 'For boundary condition outflow you have to set the'
        write(logUnit(1),*) 'pressure!'
        write(logUnit(1),*) 'Something is wrong with that in boundary ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('supersonic_outflow')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Extrapolate v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'neumann'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

      if (nDims > 1) then
        ! Extrapolate v_tangential_1
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'neumann'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        if (nDims > 2) then
          ! Extrapolate v_tangential_2
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'neumann'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )
        end if
      end if

      ! Extrapolate pressure
      bc_state(pIndex)%state_name = 'pressure'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

    case default
      write(logUnit(1),*) 'Unknown boundary kind "' // trim(bc_kind) // '"'
      write(logUnit(1),*) 'for boundary  "' // trim(bc_label) // '".'
      write(logUnit(1),*) 'Available boundary kinds for Euler equations:'
      write(logUnit(1),*) ' * slipwall / wall'
      write(logUnit(1),*) ' * primitives'
      write(logUnit(1),*) ' * conservatives'
      write(logUnit(1),*) ' * inflow and inflow_normal'
      write(logUnit(1),*) ' * supersonic_inflow_normal'
      write(logUnit(1),*) ' * outflow'
      write(logUnit(1),*) ' * supersonic_outflow'
      write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
      call tem_abort()

    end select

    call truncate( me = bc_varDict )
    call truncate( me = bc_varDict_gradient )

    if (size(bc_state) /= bc_varDict%nVals) then
      call tem_abort( 'Nr. of state variables does not match size of varDict' )
    end if

    if (size(bc_state_gradient) /= bc_varDict_gradient%nVals) then
      call tem_abort( 'Nr. of state gradient variables does not match ' &
        &             // 'size of varDict_gradient'                     )
    end if
  end subroutine atl_eqn_euler_load_bc
  ! ****************************************************************************


  ! ****************************************************************************
  ! Getting the numerical flux for Euler equations
  subroutine atl_getEulerFluxes(euler, conf, eqn_handle, eqn_dim)
    !> The equations type to set the numerical flux in.
    type(atl_euler_type), intent(inout) :: euler
    !> Configuration file handle to get the numerical flux setting from.
    type(flu_state) :: conf
    !> Handle to the equation table in the configuration script.
    integer, intent(in) :: eqn_handle
    !> Dimension of the equation to set the flux for.
    integer, intent(in) :: eqn_dim
    ! --------------------------------------------------------------------------
    character(len=labelLen) :: eq_nflux
    integer :: iError
    ! --------------------------------------------------------------------------

    call aot_get_val( L       = conf,           &
      &               thandle = eqn_handle,     &
      &               key     = 'numflux',      &
      &               val     =  eq_nflux,      &
      &               ErrCode = iError,         &
      &               default = 'lax_friedrich' )
    eq_nflux = upper_to_lower(eq_nflux)
    eq_nflux = adjustl(eq_nflux)

    select case(trim(eq_nflux))
    case ('hll')
      write(logunit(2),*) 'Using HLL numerical flux.'
      write(logunit(2),*) 'Warning, this flux ignores materials completely!'
      select case(eqn_dim)
      case(1)
        euler%numflux => atl_HLLEuler1D
      case(2)
        euler%numflux => atl_HLLEuler2D
      case(3)
        euler%numflux => atl_HLLEuler
      end select

    case ('godunov')
      write(logunit(2),*) 'Using Godunov numerical flux.'
      write(logunit(2),*) 'Warning, this flux does not handle materials'
      write(logunit(2),*) 'completely correct!'
      select case(eqn_dim)
      case(1)
        euler%numflux => atl_GodunovEuler1D
      case(2)
        euler%numflux => atl_GodunovEuler2D
      case(3)
        euler%numflux => atl_GodunovEuler
      end select

    case ('lax_friedrich')
      write(logunit(2),*) 'Using Lax Friedrichs numerical flux.'
      select case(eqn_dim)
      case(1)
        euler%numflux => atl_laxFriedEuler_1D
      case(2)
        euler%numflux => atl_laxFriedEuler_2D
      case(3)
        euler%numflux => atl_laxFriedEuler
      end select

    case default
      write(logunit(1),*) 'Unknown numerical flux ', trim(eq_nflux)
      write(logunit(1),*) 'for the Euler equation system.'
      write(logunit(1),*) 'Please choose one of the available:'
      write(logunit(1),*) ' * lax_friedrich (default)'
      write(logunit(1),*) ' * godunov'
      write(logunit(1),*) ' * hll'
      call tem_abort()
    end select
  end subroutine atl_getEulerFluxes
  ! ****************************************************************************


  ! ****************************************************************************
  !> Getting the linearization indicator for Euler equations from the config.
  !!
  !! Set the function pointer to compute  the linearization indicator according
  !! to the setting by the user.
  !! Available indicators are:
  !!
  !! - 'density' to use the maximal relative deviation in density
  !! - 'energy' to use the maximal relative deviation in energy
  !! - 'error' to use an error estimate
  !!
  !! If linear_limit is 0, the indicator is completely deactivated and
  !! euler%linear will always return .false.
  subroutine atl_getEulerLinInd(euler, conf, eqn_handle, eqn_dim)
    !> The equations type to set the numerical flux in.
    type(atl_euler_type), intent(inout) :: euler
    !> Configuration file handle to get the numerical flux setting from.
    type(flu_state) :: conf
    !> Handle to the equation table in the configuration script.
    integer, intent(in) :: eqn_handle
    !> Dimension of the equation to set the flux for.
    integer, intent(in) :: eqn_dim
    ! --------------------------------------------------------------------------
    character(len=labelLen) :: eq_linind
    integer :: iError
    ! --------------------------------------------------------------------------

    if (euler%linear_limit > 0.0_rk) then
      call aot_get_val( L       = conf,                      &
        &               thandle = eqn_handle,                &
        &               key     = 'linearization_indicator', &
        &               val     =  eq_linind,                &
        &               ErrCode = iError,                    &
        &               default = 'error'                    )
      eq_linind = upper_to_lower(eq_linind)
      eq_linind = adjustl(eq_linind)

      select case(trim(eq_linind))
      case ('density')
        write(logunit(2),*) 'Using energy as linearization indicator.'
        euler%linear => linearization_indicator_density

      case ('energy')
        write(logunit(2),*) 'Using energy as linearization indicator.'
        select case(eqn_dim)
        case(1)
          euler%linear => linearization_indicator_energy1d
        case(2)
          euler%linear => linearization_indicator_energy2d
        case(3)
          euler%linear => linearization_indicator_energy3d
        end select

      case ('error')
        write(logunit(2),*) 'Using the error estimate as linearization' &
          &                 //' indicator.'
        select case(eqn_dim)
        case(1)
          euler%linear => linearization_indicator_err1d
        case(2)
          euler%linear => linearization_indicator_err2d
        case(3)
          euler%linear => linearization_indicator_err3d
        end select
      case default
        write(logunit(1),*) 'Unknown linearization_indicator ', trim(eq_linind)
        write(logunit(1),*) 'for the Euler equation system.'
        write(logunit(1),*) 'Please choose one of the available:'
        write(logunit(1),*) ' * error (default)'
        write(logunit(1),*) ' * density'
        write(logunit(1),*) ' * energy'
        call tem_abort()
      end select
    else
      euler%linear => linearization_deactivated
    end if

  end subroutine atl_getEulerLinInd
  ! ****************************************************************************


  ! ------------------------------------------------------------------------ !
  !> Solve the equation system with just the penalization terms to find an
  !! implicit update for the IMEX timestepping procedure.
  subroutine atl_eqn_euler_implicit_pen( material, eqn, weighted_dt, nDims,  &
    &                                    poly_proj, state, timestep_rhs      )
    !> Definition of the material, which directly describes the penalization.
    !!
    !! We expect the mask function Chi to be defined in materialdat(:,:,1),
    !! the obstacle velocity U_o in materialdat(:,:,2:nDims+1) and
    !! the obstacle Temperature T_o in materialdat(:,:,nDims+2).
    type(atl_material_type), intent(in) :: material

    !> Definition of parameters in the Euler equations.
    !!
    !! This has to provide cv, the viscous permeability and the thermal
    !! permeability.
    type(atl_euler_type), intent(in) :: eqn

    !> Timestep which is already weighted by the time integration scheme.
    real(kind=rk), intent(in) :: weighted_dt

    !> Number of dimensions, the equation system is computed in (2 or 3).
    integer, intent(in) :: ndims

    !> Description of the projection for the material.
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> The state variables of the equation system, they will be updated to
    !! the solution of the implicit computation for penalization.
    real(kind=rk), intent(inout) :: state(:,:,:)

    !> Right hand side contribution by the implicit calculation.
    real(kind=rk), intent(out) :: timestep_rhs(:,:,:)
    ! -------------------------------------------------------------------- !
    integer :: iMatElem
    integer :: iElem
    integer :: iPoint
    integer :: iDir
    integer :: nElems
    integer :: nPoints
    integer :: nVars
    real(kind=rk), parameter :: numzero = 8*tiny(weighted_dt)
    real(kind=rk) :: inv_visc_perm
    real(kind=rk) :: inv_thrm_perm
    real(kind=rk) :: viscous_time_weight
    real(kind=rk) :: thermal_time_weight
    real(kind=rk) :: viscous_fact
    real(kind=rk) :: thermal_fact
    real(kind=rk) :: Chi
    real(kind=rk) :: U_o(3)
    real(kind=rk) :: T_o
    real(kind=rk) :: relvel
    real(kind=rk), allocatable :: modalCoeff(:,:), modalCoeff_cur(:,:)
    real(kind=rk), allocatable :: pointVal(:,:), cur(:,:)
    real(kind=rk), allocatable :: velocity(:,:)
    real(kind=rk), allocatable :: velmag(:)
    real(kind=rk), allocatable :: temperature(:)
    ! -------------------------------------------------------------------- !

    select case(nDims)
    case (1)
      nPoints = poly_proj%body_1d%nquadpoints
    case (2)
      nPoints = poly_proj%body_2d%nquadpoints
    case (3)
      nPoints = poly_proj%body_3d%nquadpoints
    end select

    nVars = nDims + 2

    inv_visc_perm = 1.0_rk / eqn%viscous_permeability
    inv_thrm_perm = 1.0_rk / eqn%thermal_permeability

    viscous_time_weight = weighted_dt / eqn%viscous_permeability
    thermal_time_weight = weighted_dt / eqn%thermal_permeability

    allocate(modalCoeff(nPoints,nVars))
    allocate(pointVal(nPoints,nVars))
    allocate(cur(nPoints,nVars-1)) ! Do not need to consider density
    allocate(modalCoeff_cur(nPoints,nVars-1))
    allocate(velocity(nPoints, nDims))
    allocate(temperature(nPoints))
    allocate(velmag(nPoints))

    ! No right hand side contribution for the continuity equation.
    timestep_rhs(:,:,1) = 0.0_rk

    nElems = material%material_desc%computeElems(1)%nElems
    constElems: do iMatElem=1,nElems
      iElem = material%material_desc%computeElems(1)%totElemIndices(iMatElem)

      Chi = material%material_dat%elemMaterialData(1)%materialDat(iMatElem,1,1)
      do iDir=1,nDims
        U_o(iDir) = material%material_dat%elemMaterialData(1)          &
          &                              %materialDat(iMatElem,1,iDir+1)
      end do
      T_o = material%material_dat%elemMaterialData(1) &
        &                        %materialDat(iMatElem,1,nVars)

      if ( abs(Chi) > numzero ) then
        ! Only need to compute the penalization if the constant Chi in this
        ! element is not 0.

        ! Need to convert to nodal as we have to divide by the density.
        call ply_convert2oversample( state       = state(iElem,:,:nVars), &
          &                          ndim        = nDims,                 &
          &                          poly_proj   = poly_proj,             &
          &                          modalCoeffs = modalCoeff             )

        call ply_poly_project_m2n( me         = poly_proj, &
          &                        dim        = nDims,     &
          &                        nVars      = nVars,     &
          &                        nodal_data = pointVal,  &
          &                        modal_data = modalCoeff )

        viscous_fact = Chi*viscous_time_weight
        thermal_fact = Chi*thermal_time_weight

        velmag = 0.0_rk
        relvel = 0.0_rk
        do iDir=1,nDims
          do iPoint=1,nPoints
            velocity(iPoint, iDir) &
              &  = (pointVal(iPoint,iDir+1) + viscous_fact*U_o(iDir)) &
              &    / (pointVal(iPoint,1) + viscous_fact)
            velmag(iPoint) = velmag(iPoint) + velocity(iPoint,iDir)**2
            relvel = relvel + (velocity(iPoint,iDir)-U_o(iDir)) &
              &              * velocity(iPoint,iDir)
          end do
        end do

        do iPoint=1,nPoints
          temperature(iPoint) &
            &  = ( thermal_fact*T_o + pointVal(iPoint,nVars)    &
            &                       - 0.5_rk*pointVal(iPoint,1) &
            &                               *velmag(iPoint)     &
            &                       - viscous_fact * relvel )   &
            &    / ( eqn%cv*pointVal(iPoint,1) + thermal_fact   )
          temperature(iPoint) = max(temperature(iPoint), 0.01*T_o)
           
          cur(iPoint,nDims+1) = Chi*inv_thrm_perm*( T_o - temperature(iPoint))
        end do 
        do iDir=1,nDims
          do iPoint=1,nPoints
            cur(iPoint,iDir) = Chi*inv_visc_perm*( U_o(iDir)               &
              &                                    - velocity(iPoint,iDir) )
            pointVal(iPoint,iDir+1) = pointVal(iPoint,1)*velocity(iPoint,iDir)
          end do
        end do
        do iPoint=1,nPoints
          pointVal(iPoint,nVars) = pointVal(iPoint,1)               &
            &                      * ( 0.5_rk*velmag(iPoint)        &
            &                          + eqn%cv*temperature(iPoint) )
        end do

        ! Convert the updated state back (u_i)
        call ply_poly_project_n2m( me         = poly_proj, &
          &                        dim        = nDims,     &
          &                        nVars      = nVars,     &
          &                        nodal_data = pointVal,  &
          &                        modal_data = modalCoeff )
        call ply_convertFromOversample( modalCoeffs = modalCoeff,      &
          &                             ndim        = nDims,           &
          &                             poly_proj   = poly_proj,       &
          &                             state       = state(iElem,:,:) )

        ! Convert the right hand side back g(u_i)
        call ply_poly_project_n2m( me         = poly_proj,     &
          &                        dim        = nDims,         &
          &                        nVars      = nVars-1,       &
          &                        nodal_data = cur,           &
          &                        modal_data = modalCoeff_cur )
        call ply_convertFromOversample( modalCoeffs = modalCoeff_cur,          &
          &                             ndim        = nDims,                   &
          &                             poly_proj   = poly_proj,               &
          &                             state       = timestep_rhs(iElem,:,2:) )

      else
        timestep_rhs(iElem,:,2:) = 0.0_rk
      end if

    end do constElems


    nElems = material%material_desc%computeElems(2)%nElems
    varElems: do iMatElem=1,nElems
      iElem = material%material_desc%computeElems(2)%totElemIndices(iMatElem)

      call ply_convert2oversample( state       = state(iElem,:,:nVars), &
        &                          ndim        = nDims,                 &
        &                          poly_proj   = poly_proj,             &
        &                          modalCoeffs = modalCoeff             )

      call ply_poly_project_m2n( me         = poly_proj, &
        &                        dim        = nDims,     &
        &                        nVars      = nVars,     &
        &                        nodal_data = pointVal,  &
        &                        modal_data = modalCoeff )

      do iPoint=1,nPoints
        Chi = material%material_dat%elemMaterialData(2) &
          &           %materialDat(iMatElem,iPoint,1)
        do iDir=1,nDims
          U_o(iDir) = material%material_dat%elemMaterialData(2)          &
            &                              %materialDat(iMatElem,iPoint,iDir+1)
        end do
        T_o = material%material_dat%elemMaterialData(2) &
          &                        %materialDat(iMatElem,iPoint,nVars)

        viscous_fact = Chi*viscous_time_weight
        thermal_fact = Chi*thermal_time_weight

        velmag(iPoint) = 0.0_rk
        relvel = 0.0_rk
        do iDir=1,nDims
          velocity(iPoint, iDir) &
            &  = (pointVal(iPoint,iDir+1) + viscous_fact*U_o(iDir)) &
            &    / (pointVal(iPoint,1) + viscous_fact)
          velmag(iPoint) = velmag(iPoint) + velocity(iPoint, iDir)**2
          relvel = relvel + (velocity(iPoint, iDir) - U_o(iDir)) &
                           * velocity(iPoint, iDir)
        end do
          temperature(iPoint) &
            &  = ( thermal_fact*T_o + pointVal(iPoint,nVars)    &
            &                       - 0.5_rk*pointVal(iPoint,1) &
            &                               *velmag(iPoint)     &
            &                       - viscous_fact * relvel)    &
            &    / ( eqn%cv*pointVal(iPoint,1) + thermal_fact   )
          temperature(iPoint) = max(temperature(iPoint), 0.01*T_o)
        cur(iPoint,nDims+1) = Chi*inv_thrm_perm*( T_o - temperature(iPoint))

        do iDir=1,nDims
          cur(iPoint,iDir) = Chi*inv_visc_perm*( U_o(iDir)               &
            &                                    - velocity(iPoint,iDir) )
          pointVal(iPoint,iDir+1) = pointVal(iPoint,1)*velocity(iPoint,iDir)
        end do
        pointVal(iPoint,nVars) = pointVal(iPoint,1)               &
          &                      * ( 0.5_rk*velmag(iPoint)        &
          &                          + eqn%cv*temperature(iPoint) )

      end do

      ! Convert the updated state back (u_i)
      call ply_poly_project_n2m( me         = poly_proj, &
        &                        dim        = nDims,     &
        &                        nVars      = nVars,     &
        &                        nodal_data = pointVal,  &
        &                        modal_data = modalCoeff )
      call ply_convertFromOversample( modalCoeffs = modalCoeff,      &
        &                             ndim        = nDims,           &
        &                             poly_proj   = poly_proj,       &
        &                             state       = state(iElem,:,:) )

      ! Convert the right hand side back g(u_i)
      call ply_poly_project_n2m( me         = poly_proj,     &
        &                        dim        = nDims,         &
        &                        nVars      = nVars-1,       &
        &                        nodal_data = cur,           &
        &                        modal_data = modalCoeff_cur )
      call ply_convertFromOversample( modalCoeffs = modalCoeff_cur,          &
        &                             ndim        = nDims,                   &
        &                             poly_proj   = poly_proj,               &
        &                             state       = timestep_rhs(iElem,:,2:) )

    end do varElems

    deallocate(modalCoeff)
    deallocate(pointVal)
    deallocate(cur)
    deallocate(modalCoeff_cur)
    deallocate(velocity)
    deallocate(temperature)
    deallocate(velmag)

  end subroutine atl_eqn_euler_implicit_pen
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> An indicator that completely deactivates linearization.
  pure function linearization_deactivated(euler, mean, deviation) &
    &        result(islinear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_euler_type), intent(in) :: euler

    !> The mean value of each state
    real(kind=rk), intent(in) :: mean(:)

    !> Maximal deviation of each state
    real(kind=rk), intent(in) :: deviation(:)

    logical :: islinear
    ! -------------------------------------------------------------------- !

    islinear = .false.

  end function linearization_deactivated
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> An indicator to decide whether linearization of fluxes is tolerable
  !! based on the density.
  pure function linearization_indicator_density(euler, mean, deviation) &
    &        result(islinear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_euler_type), intent(in) :: euler

    !> The mean value of each state
    real(kind=rk), intent(in) :: mean(:)

    !> Maximal deviation of each state
    real(kind=rk), intent(in) :: deviation(:)

    logical :: islinear
    ! -------------------------------------------------------------------- !

    islinear = ( deviation(1) < ( euler%linear_limit * mean(1) ) )

  end function linearization_indicator_density
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> An indicator to decide whether linearization of fluxes is tolerable
  !! based on the energy in 3D.
  pure function linearization_indicator_energy3d(euler, mean, deviation) &
    &        result(islinear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_euler_type), intent(in) :: euler

    !> The mean value of each state
    real(kind=rk), intent(in) :: mean(:)

    !> Maximal deviation of each state
    real(kind=rk), intent(in) :: deviation(:)

    logical :: islinear
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    islinear = ( deviation(5) < ( euler%linear_limit &
      &                           * mean(5) )        )

  end function linearization_indicator_energy3d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> An indicator to decide whether linearization of fluxes is tolerable
  !! based on the energy in 2D.
  pure function linearization_indicator_energy2d(euler, mean, deviation) &
    &        result(islinear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_euler_type), intent(in) :: euler

    !> The mean value of each state
    real(kind=rk), intent(in) :: mean(:)

    !> Maximal deviation of each state
    real(kind=rk), intent(in) :: deviation(:)

    logical :: islinear
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    islinear = ( deviation(4) < ( euler%linear_limit &
      &                           * mean(4) )        )

  end function linearization_indicator_energy2d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> An indicator to decide whether linearization of fluxes is tolerable
  !! based on the energy in 2D.
  pure function linearization_indicator_energy1d(euler, mean, deviation) &
    &        result(islinear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_euler_type), intent(in) :: euler

    !> The mean value of each state
    real(kind=rk), intent(in) :: mean(:)

    !> Maximal deviation of each state
    real(kind=rk), intent(in) :: deviation(:)

    logical :: islinear
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    islinear = ( deviation(3) < ( euler%linear_limit &
      &                           * mean(3) )        )

  end function linearization_indicator_energy1d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> An indicator to decide whether linearization of fluxes is tolerable
  !! based on the error estimate.
  pure function linearization_indicator_err3d(euler, mean, deviation) &
    &        result(islinear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_euler_type), intent(in) :: euler

    !> The mean value of each state
    real(kind=rk), intent(in) :: mean(:)

    !> Maximal deviation of each state
    real(kind=rk), intent(in) :: deviation(:)

    logical :: islinear
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: m_dev_mag, m_0_mag, m_max
    real(kind=rk) :: rho_min
    real(kind=rk) :: M_err, E_err
    real(kind=rk) :: gam
    ! -------------------------------------------------------------------- !

    gam = euler%isen_coef

    m_dev_mag = sqrt(deviation(2)**2 + deviation(3)**2 + deviation(4)**2)
    m_0_mag = sqrt(mean(2)**2 + mean(3)**2 + mean(4)**2)
    m_max = sqrt( (abs(mean(2)) + deviation(2))**2   &
      &           + (abs(mean(3)) + deviation(3))**2 &
      &           + (abs(mean(4)) + deviation(4))**2 )
    M_err = mean(1)*m_dev_mag + m_0_mag*deviation(1)
    E_err = mean(1)*deviation(5) + mean(5)*deviation(1)
    rho_min = max(mean(1) - deviation(1), 0.0_rk)

    islinear = ( M_err * ( (gam-1.0_rk)*(m_0_mag + 0.5_rk*m_max)*M_err &
      &                    + gam*rho_min*E_err                       ) &
      &          < (mean(1)*rho_min)**2 * euler%linear_limit           )

  end function linearization_indicator_err3d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> An indicator to decide whether linearization of fluxes is tolerable
  !! based on the error estimate.
  pure function linearization_indicator_err2d(euler, mean, deviation) &
    &        result(islinear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_euler_type), intent(in) :: euler

    !> The mean value of each state
    real(kind=rk), intent(in) :: mean(:)

    !> Maximal deviation of each state
    real(kind=rk), intent(in) :: deviation(:)

    logical :: islinear
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: m_dev_mag, m_0_mag, m_max
    real(kind=rk) :: rho_min
    real(kind=rk) :: M_err, E_err
    real(kind=rk) :: gam
    ! -------------------------------------------------------------------- !

    gam = euler%isen_coef

    m_dev_mag = sqrt(deviation(2)**2 + deviation(3)**2)
    m_0_mag = sqrt(mean(2)**2 + mean(3)**2)
    m_max = sqrt( (abs(mean(2)) + deviation(2))**2   &
      &           + (abs(mean(3)) + deviation(3))**2 )
    M_err = mean(1)*m_dev_mag + m_0_mag*deviation(1)
    E_err = mean(1)*deviation(4) + mean(4)*deviation(1)
    rho_min = max(mean(1) - deviation(1), 0.0_rk)

    islinear = ( M_err * ( (gam-1.0_rk)*(m_0_mag + 0.5_rk*m_max)*M_err &
      &                    + gam*rho_min*E_err                       ) &
      &          < (mean(1)*rho_min)**2 * euler%linear_limit           )

  end function linearization_indicator_err2d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> An indicator to decide whether linearization of fluxes is tolerable
  !! based on the error estimate.
  pure function linearization_indicator_err1d(euler, mean, deviation) &
    &        result(islinear)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_euler_type), intent(in) :: euler

    !> The mean value of each state
    real(kind=rk), intent(in) :: mean(:)

    !> Maximal deviation of each state
    real(kind=rk), intent(in) :: deviation(:)

    logical :: islinear
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: m_dev_mag, m_0_mag, m_max
    real(kind=rk) :: rho_min
    real(kind=rk) :: M_err, E_err
    real(kind=rk) :: gam
    ! -------------------------------------------------------------------- !

    gam = euler%isen_coef

    m_dev_mag = deviation(2)
    m_0_mag = abs(mean(2))
    m_max = abs(mean(2)) + deviation(2)
    M_err = mean(1)*m_dev_mag + m_0_mag*deviation(1)
    E_err = mean(1)*deviation(3) + mean(3)*deviation(1)
    rho_min = max(mean(1) - deviation(1), 0.0_rk)

    islinear = ( M_err * ( (gam-1.0_rk)*(m_0_mag + 0.5_rk*m_max)*M_err &
      &                    + gam*rho_min*E_err                       ) &
      &          < (mean(1)*rho_min)**2 * euler%linear_limit           )

  end function linearization_indicator_err1d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module atl_eqn_euler_hlp_module
