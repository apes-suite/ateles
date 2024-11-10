! Copyright (c) 2013-2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

module atl_eqn_nvrstk_hlp_module
  use env_module, only: rk
  use aotus_module,                   only: flu_State

  use tem_aux_module,                 only: tem_abort
  use tem_bc_module,                  only: tem_bc_state_type
  use tem_logging_module,             only: logUnit
  use tem_stringKeyValuePair_module,  only: tem_stringKeyValuePair_type,      &
    &                                       grw_stringKeyValuePairArray_type, &
    &                                       init, truncate, append

  use atl_equation_module,            only: atl_equations_type,     &
    &                                       atl_eqn_var_trafo_type
  use atl_bc_state_module,            only: atl_load_bc_state
  use atl_eqn_nvrstk_module,          only: atl_load_NavierStokes, &
    &                                       atl_navierstokes_type
  use atl_eqn_nvrstk_var_module,      only: atl_append_nvrstk_derivedVars
  use atl_eqn_euler_derive_module,    only: atl_eqn_euler_cons2primTemp,   &
    &                                       atl_eqn_euler_primTemp2Cons,   &
    &                                       atl_eqn_euler_cons2primVel,    &
    &                                       atl_eqn_euler_primVel2Cons,    &
    &                                       atl_eqn_euler_cons2prim,       &
    &                                       atl_eqn_euler_prim2Cons,       &
    &                                       atl_eqn_euler_cons2prim_grad,  &
    &                                       atl_eqn_euler_prim2Cons_grad,  &
    &                                       atl_eqn_euler_cons2prim_elems, &
    &                                       atl_eqn_euler_prim2Cons_elems
  use atl_eqn_euler_2d_derive_module, only: atl_eqn_euler_2d_cons2primTemp,   &
    &                                       atl_eqn_euler_2d_primTemp2Cons,   &
    &                                       atl_eqn_euler_2d_cons2primVel,    &
    &                                       atl_eqn_euler_2d_primVel2Cons,    &
    &                                       atl_eqn_euler_2d_cons2prim,       &
    &                                       atl_eqn_euler_2d_prim2Cons,       &
    &                                       atl_eqn_euler_2d_cons2prim_grad,  &
    &                                       atl_eqn_euler_2d_prim2Cons_grad,  &
    &                                       atl_eqn_euler_2d_cons2prim_elems, &
    &                                       atl_eqn_euler_2d_prim2Cons_elems
  use atl_eqn_euler_var_module,       only: atl_init_euler_vars,        &
    &                                       atl_init_euler_sourceTerms, &
    &                                       atl_init_euler_material
  use atl_eqn_euler_hlp_module,       only: atl_getEulerFluxes, &
    &                                       atl_getEulerLinInd
  use atl_eqn_euler_2d_var_module,    only: atl_init_euler_2d_vars,       &
    &                                       atl_init_euler_2d_sourceTerms
  use atl_varSys_module,              only: atl_varSys_solverData_type
  use atl_source_types_module,        only: atl_init_source_type
  use atl_materialPrp_module,         only: atl_init_material_type

  implicit none

  private

  public :: atl_eqn_nvrstk_load_bc
  public :: atl_eqn_nvrstk_init


contains

  !> Initialization of the Navier-Stokes equations.
  !!
  !! This routine sets up the necessary infrastructure for the Navier-Stokes
  !! equations.
  !! It reads the configuration from the given script in conf under the table
  !! provided in thandle and sets function pointers and variables accordingly.
  subroutine atl_eqn_nvrstk_init( conf, thandle, equation, nDimensions, &
    &                             initSource, initMaterial, varSys_data )
    ! ---------------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------------

    equation%isNonlinear = .true.
    ! we need 1st spatial derivatives
    equation%nDerivatives = 1
    equation%nDimensions = nDimensions
    ! timesetp is dynamic and changes over simulationtime
    equation%adaptive_timestep = .true.


    select case(nDimensions)
    case(2)
      equation%load_bc => atl_eqn_nvrstk_load_bc
      call atl_init_euler_2d_vars( &
        & equation   = equation,   &
        & solverData = varSys_data )

      call atl_init_euler_2d_sourceTerms( initSource%poss_srcVars, &
        &                                 initSource%eval_source   )

      equation%cons2prim => atl_eqn_euler_2d_cons2prim_elems
      equation%prim2cons => atl_eqn_euler_2d_prim2cons_elems

    case(3)
      equation%load_bc => atl_eqn_nvrstk_load_bc

      call atl_init_euler_vars(    &
        & equation   = equation,   &
        & solverData = varSys_data )

      call atl_init_euler_sourceTerms( initSource%poss_srcVars, &
        &                              initSource%eval_source   )

      equation%cons2prim => atl_eqn_euler_cons2prim_elems
      equation%prim2cons => atl_eqn_euler_prim2cons_elems
    end select

    call atl_load_navierstokes( NavierStokes = equation%NavierStokes, &
      &                         Euler        = equation%Euler,        &
      &                         conf         = conf,                  &
      &                         eq_table     = thandle                )
    equation%requiresDeviation = (equation%euler%linear_limit > 0.0_rk)

    ! Getting the numerical flux !
    call atl_getEulerFluxes(euler      = equation%euler, &
      &                     conf       = conf,           &
      &                     eqn_handle = thandle,        &
      &                     eqn_dim    = nDimensions     )

    ! Getting the indicator to use in linearization !
    call atl_getEulerLinInd(euler      = equation%Euler, &
      &                     conf       = conf,           &
      &                     eqn_handle = thandle,        &
      &                     eqn_dim    = nDimensions     )

    call atl_init_euler_material(                     &
      & possVars    = initMaterial%poss_materialVars, &
      & nDimensions = nDimensions                     )

    equation%requires_gradmax = (equation%NavierStokes%visc_limit > 0.0_rk)
    if (equation%requires_gradmax) then
      select case(nDimensions)
      case(2)
        equation%NavierStokes%inviscous => inviscous_indicator_2d
      case(3)
        equation%NavierStokes%inviscous => inviscous_indicator_3d
      end select
    else
      equation%NavierStokes%inviscous => inviscous_deactivated
    end if

    call atl_append_nvrstk_derivedVars(  &
      &    varSys     = equation%varSys, &
      &    solverData = varSys_data      )

  end subroutine atl_eqn_nvrstk_init


  subroutine atl_eqn_nvrstk_load_bc( equation,                              &
    &                                bc_state, bc_state_gradient,           &
    &                                bc_varDict, bc_varDict_gradient,       &
    &                                bc_normal_vec, bc_normal_vec_gradient, &
    &                                bc_trafo, bc_trafo_gradient,           &
    &                                bc_label, bc_kind, thandle, conf       )
    ! ---------------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------------
    integer :: iVar
    integer :: nDims
    integer :: pIndex
    type(tem_stringKeyValuePair_type) :: kvp
    ! ---------------------------------------------------------------------------

    nDims = equation%nDimensions

    allocate(bc_state(equation%varSys%nScalars))
    allocate(bc_state_gradient(equation%varSys%nScalars))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )
    ! Constant zero variable for non-configurable boundary variable
    kvp%value = 'zero_const'

    pIndex = equation%varSys%nScalars

    ! Check for Navier-Stokes specific boundary conditions.
    select case(bc_kind)
    case('isothermal_wall')
      ! Use a non-slip boundary for walls in the Navier-Stokes equations.
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2primTemp
        bc_trafo%from => atl_eqn_euler_2d_primTemp2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2primTemp
        bc_trafo_gradient%from => atl_eqn_euler_2d_primTemp2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2primTemp
        bc_trafo%from => atl_eqn_euler_primTemp2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2primTemp
        bc_trafo_gradient%from => atl_eqn_euler_primTemp2cons
      end select
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! For the viscous terms ...
      ! ... the gradients on the face
      bc_trafo_gradient%identity = .true.
      bc_normal_vec_gradient = .false.

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      !!VK! density for viscous terms
      !!VKbc_state_gradient(1,1) = bc_state(1)
      !!VKkvp%key = trim(bc_state_gradient(1,1)%state_name)
      !!VKcall append( me = bc_varDict_gradient, val = kvp )

      ! ... it's gradient (extrapolate)
      bc_state_gradient(1) = bc_state(1)
      kvp%key = trim(bc_state_gradient(1)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

      ! Prescribe v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'dirichlet'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

!!VK      ! v_normal for viscous terms
!!VK      bc_state_gradient(2,1)%state_name = 'v_norm'
!!VK      bc_state_gradient(2,1)%style = 'dirichlet_mirror'
!!VK      bc_state_gradient(2,1)%isDefined = .true.
!!VK      kvp%key = trim(bc_state_gradient(2,1)%state_name)
!!VK      call append( me = bc_varDict_gradient, val = kvp )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.
      kvp%key = trim(bc_state_gradient(2)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

      if (nDims > 1) then
        ! Impose tangential velocity
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'dirichlet'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

!!VK        ! v_tan for viscous terms
!!VK        call atl_load_bc_state( bc          = bc_state_gradient(3,1), &
!!VK          &                     state_name  = 'v_tan',                &
!!VK          &                     style       = 'dirichlet_mirror',     &
!!VK          &                     conf        = conf,                   &
!!VK          &                     bc_handle   = thandle,                &
!!VK          &                     varDict     = bc_varDict_gradient,    &
!!VK          &                     varSys      = equation%varSys         )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.
        kvp%key = trim(bc_state_gradient(3)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )

        if (nDims > 2) then
          ! Impose tangential velocity
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'dirichlet'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )

!!VK          ! v_tan for viscous terms
!!VK          call atl_load_bc_state( bc          = bc_state_gradient(4,1), &
!!VK            &                     state_name  = 'v_tan2',               &
!!VK            &                     style       = 'dirichlet_mirror',     &
!!VK            &                     conf        = conf,                   &
!!VK            &                     bc_handle   = thandle,                &
!!VK            &                     varDict     = bc_varDict_gradient,    &
!!VK            &                     varSys      = equation%varSys         )

          ! ... it's gradient (momentum - extrapolate)
          bc_state_gradient(4)%state_name = 'momentum_tan2'
          bc_state_gradient(4)%style = 'neumann'
          bc_state_gradient(4)%isDefined = .true.
          kvp%key = trim(bc_state_gradient(4)%state_name)
          call append( me = bc_varDict_gradient, val = kvp )
        end if
      end if

      ! Prescribe temperature
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'temperature',    &
        &                     style       = 'dirichlet',      &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varDict     = bc_varDict,       &
        &                     varSys      = equation%varSys   )

!!VK      ! Temperature for viscous terms
!!VK      call atl_load_bc_state( bc          = bc_state_gradient(pIndex,1), &
!!VK        &                     state_name  = 'temperature',               &
!!VK        &                     style       = 'dirichlet_mirror',          &
!!VK        &                     conf        = conf,                        &
!!VK        &                     bc_handle   = thandle,                     &
!!VK        &                     varDict     = bc_varDict_gradient,         &
!!VK        &                     varSys      = equation%varSys              )

      ! ... it's gradient (energy - extrapolate)
      bc_state_gradient(pIndex)%state_name = 'energy'
      bc_state_gradient(pIndex)%style = 'neumann'
      bc_state_gradient(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state_gradient(pIndex)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

    ! Adiabatic wall
    case('wall', 'adiabatic_wall')
      ! Use a non-slip boundary for walls in the Navier-Stokes equations.
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2primVel
        bc_trafo%from => atl_eqn_euler_2d_primVel2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2primVel
        bc_trafo_gradient%from => atl_eqn_euler_2d_primVel2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2primVel
        bc_trafo%from => atl_eqn_euler_primVel2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2primVel
        bc_trafo_gradient%from => atl_eqn_euler_primVel2cons
      end select
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! For the viscous terms ...
!!VK      ! ... the face values
!!VK      bc_trafo_gradient(1)%identity = .false.
!!VK      bc_normal_vec_gradient(1) = .true.
      ! ... the gradients on the face
      bc_trafo_gradient%identity = .true.
      bc_normal_vec_gradient = .false.

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

!!VK      ! density for viscous terms
!!VK      bc_state_gradient(1,1) = bc_state(1)
!!VK      kvp%key = trim(bc_state_gradient(1,1)%state_name)
!!VK      call append( me = bc_varDict_gradient, val = kvp )

      ! ... it's gradient (extrapolate)
      bc_state_gradient(1) = bc_state(1)
      kvp%key = trim(bc_state_gradient(1)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

      ! Prescribe v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'dirichlet'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

!!VK      ! v_normal for viscous terms
!!VK      bc_state_gradient(2,1)%state_name = 'v_norm'
!!VK      bc_state_gradient(2,1)%style = 'dirichlet_mirror'
!!VK      bc_state_gradient(2,1)%isDefined = .true.
!!VK      kvp%key = trim(bc_state_gradient(2,1)%state_name)
!!VK      call append( me = bc_varDict_gradient, val = kvp )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.
      kvp%key = trim(bc_state_gradient(2)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

      if (nDims > 1) then
        ! Impose tangential velocity
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'dirichlet'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

!!VK        ! v_tan for viscous terms
!!VK        bc_state_gradient(3,1)%state_name = 'v_tan'
!!VK        bc_state_gradient(3,1)%style = 'dirichlet_mirror'
!!VK        bc_state_gradient(3,1)%isDefined = .true.
!!VK        kvp%key = trim(bc_state_gradient(3,1)%state_name)
!!VK        call append( me = bc_varDict_gradient, val = kvp )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.
        kvp%key = trim(bc_state_gradient(3)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )

        if (nDims > 2) then
          ! Impose tangential velocity
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'dirichlet'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )

!!VK          ! v_tan for viscous terms
!!VK          bc_state_gradient(4,1)%state_name = 'v_tan2'
!!VK          bc_state_gradient(4,1)%style = 'dirichlet_mirror'
!!VK          bc_state_gradient(4,1)%isDefined = .true.
!!VK          kvp%key = trim(bc_state_gradient(4,1)%state_name)
!!VK          call append( me = bc_varDict_gradient, val = kvp )

          ! ... it's gradient (momentum - extrapolate)
          bc_state_gradient(4)%state_name = 'momentum_tan2'
          bc_state_gradient(4)%style = 'neumann'
          bc_state_gradient(4)%isDefined = .true.
          kvp%key = trim(bc_state_gradient(4)%state_name)
          call append( me = bc_varDict_gradient, val = kvp )
        end if
      end if

      ! Extrapolate energy
      bc_state(pIndex)%state_name = 'energy'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

!!VK      ! Energy for viscous terms
!!VK      bc_state_gradient(pIndex,1) = bc_state(pIndex)
!!VK      kvp%key = trim(bc_state_gradient(pIndex,1)%state_name)
!!VK      call append( me = bc_varDict_gradient, val = kvp )

      ! ... it's gradient (energy - extrapolate)
      bc_state_gradient(pIndex) = bc_state(pIndex)
      kvp%key = trim(bc_state_gradient(pIndex)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

    ! Adiabatic slip wall
    case('slipwall', 'adiabatic_slipwall')
      ! Use a non-slip boundary for walls in the Navier-Stokes equations.
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2primVel
        bc_trafo%from => atl_eqn_euler_2d_primVel2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2primVel
        bc_trafo_gradient%from => atl_eqn_euler_2d_primVel2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2primVel
        bc_trafo%from => atl_eqn_euler_primVel2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2primVel
        bc_trafo_gradient%from => atl_eqn_euler_primVel2cons
      end select

      bc_trafo%identity = .false.
      bc_normal_vec = .true.

!!VK      ! For the viscous terms ...
!!VK      ! ... the face values
!!VK      bc_trafo_gradient(1)%identity = .false.
!!VK      bc_normal_vec_gradient(1) = .true.
      ! ... the gradients on the face
      bc_trafo_gradient%identity = .true.
      bc_normal_vec_gradient = .false.

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

!!VK      ! density for viscous terms
!!VK      bc_state_gradient(1,1) = bc_state(1)
!!VK      kvp%key = trim(bc_state_gradient(1,1)%state_name)
!!VK      call append( me = bc_varDict_gradient, val = kvp )

      ! ... it's gradient (extrapolate)
      bc_state_gradient(1) = bc_state(1)
      kvp%key = trim(bc_state_gradient(1)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

      ! Prescribe v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'dirichlet'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

!!VK      ! v_normal for viscous terms
!!VK      bc_state_gradient(2,1)%state_name = 'v_norm'
!!VK      bc_state_gradient(2,1)%style = 'dirichlet_mirror'
!!VK      bc_state_gradient(2,1)%isDefined = .true.
!!VK      kvp%key = trim(bc_state_gradient(2,1)%state_name)
!!VK      call append( me = bc_varDict_gradient, val = kvp )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.
      kvp%key = trim(bc_state_gradient(2)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

      if (nDims > 1) then
        ! Extrapolate tangential velocity
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'neumann'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

!!VK        ! v_tan for viscous terms
!!VK        bc_state_gradient(3,1)%state_name = 'v_tan'
!!VK        bc_state_gradient(3,1)%style = 'neumann'
!!VK        bc_state_gradient(3,1)%isDefined = .true.
!!VK        kvp%key = trim(bc_state_gradient(3,1)%state_name)
!!VK        call append( me = bc_varDict_gradient, val = kvp )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.
        kvp%key = trim(bc_state_gradient(3)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )

        if (nDims > 2) then
          ! Extrapolate tangential velocity
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'neumann'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )

!!VK          ! v_tan for viscous terms
!!VK          bc_state_gradient(4,1)%state_name = 'v_tan2'
!!VK          bc_state_gradient(4,1)%style = 'neumann'
!!VK          bc_state_gradient(4,1)%isDefined = .true.
!!VK          kvp%key = trim(bc_state_gradient(4,1)%state_name)
!!VK          call append( me = bc_varDict_gradient, val = kvp )

          ! ... it's gradient (momentum - extrapolate)
          bc_state_gradient(4)%state_name = 'momentum_tan2'
          bc_state_gradient(4)%style = 'neumann'
          bc_state_gradient(4)%isDefined = .true.
          kvp%key = trim(bc_state_gradient(4)%state_name)
          call append( me = bc_varDict_gradient, val = kvp )
        end if
      end if

      ! Extrapolate energy
      bc_state(pIndex)%state_name = 'energy'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

!!VK      ! Energy for viscous terms
!!VK      bc_state_gradient(pIndex,1) = bc_state(pIndex)
!!VK      kvp%key = trim(bc_state_gradient(pIndex,1)%state_name)
!!VK      call append( me = bc_varDict_gradient, val = kvp )

      ! ... it's gradient (energy - extrapolate)
      bc_state_gradient(pIndex) = bc_state(pIndex)
      kvp%key = trim(bc_state_gradient(pIndex)%state_name)
      call append( me = bc_varDict_gradient, val = kvp )

    case('primitives')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2prim
        bc_trafo%from => atl_eqn_euler_2d_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons_grad
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons_grad
      end select
      bc_trafo%identity = .false.
      bc_trafo_gradient%identity = .false.

      bc_normal_vec_gradient = .false.
      bc_normal_vec = .false.

      ! Impose density
      call atl_load_bc_state( bc          = bc_state(1),    &
        &                     state_name  = 'density',      &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varSys )

      ! Impose velocity x
      call atl_load_bc_state( bc          = bc_state(2),    &
        &                     state_name  = 'velocityX',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varSys )
      if (nDims > 1) then
        ! Impose velocity y
        call atl_load_bc_state( bc          = bc_state(3),    &
          &                     state_name  = 'velocityY',    &
          &                     conf        = conf,           &
          &                     bc_handle   = thandle,        &
          &                     varDict     = bc_varDict,     &
          &                     varSys      = equation%varSys )

        bc_state_gradient(3)%state_name = 'velocityY'
        if (nDims > 2) then
          ! Impose velocity z
          call atl_load_bc_state( bc          = bc_state(4),    &
            &                     state_name  = 'velocityZ',    &
            &                     conf        = conf,           &
            &                     bc_handle   = thandle,        &
            &                     varDict     = bc_varDict,     &
            &                     varSys      = equation%varSys )
          bc_state_gradient(4)%state_name = 'velocityZ'
        end if
      end if

      ! Impose pressure
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varDict     = bc_varDict,       &
        &                     varSys      = equation%varSys   )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition primitives you have to set'
        write(logUnit(1),*) 'all primitive variables (density, velocityX, '//&
          &                 'velocityY, velocityZ'
        write(logUnit(1),*) 'and pressure) this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

      bc_state_gradient(1)%state_name = 'density'
      bc_state_gradient(2)%state_name = 'velocityX'
      bc_state_gradient(pIndex)%state_name = 'pressure'
      bc_state_gradient(:)%style = 'neumann'
      bc_state_gradient(:)%isDefined = .true.
      do iVar=1,pIndex
        kvp%key = trim(bc_state_gradient(iVar)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )
      end do

    case('grad_primitives')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2prim
        bc_trafo%from => atl_eqn_euler_2d_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons_grad
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons_grad
      end select
      bc_trafo%identity = .false.
      bc_trafo_gradient%identity = .false.

      bc_normal_vec_gradient = .false.
      bc_normal_vec = .false.

      ! set the bc style for all gradients to dirchlet
      bc_state_gradient(:)%style = 'dirichlet'
      ! Impose density
      call atl_load_bc_state( bc          = bc_state(1),    &
        &                     state_name  = 'density',      &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varSys )
      ! and the gradient of density
      call atl_load_bc_state( bc          =  bc_state_gradient(1), &
        &                     state_name  = 'grad_density',        &
        &                     conf        = conf,                  &
        &                     bc_handle   = thandle,               &
        &                     varDict     = bc_varDict_gradient,   &
        &                     varSys      = equation%varSys        )

      ! Impose velocity x
      call atl_load_bc_state( bc          = bc_state(2),    &
        &                     state_name  = 'velocityX',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varSys )
      ! and he gradient of velocity x
      call atl_load_bc_state( bc          = bc_state_gradient(2), &
        &                     state_name  = 'grad_velocityX',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict_gradient,  &
        &                     varSys      = equation%varSys       )

      if (nDims > 1) then
        ! Impose velocity y
        call atl_load_bc_state( bc          = bc_state(3),    &
          &                     state_name  = 'velocityY',    &
          &                     conf        = conf,           &
          &                     bc_handle   = thandle,        &
          &                     varDict     = bc_varDict,     &
          &                     varSys      = equation%varSys )
        ! and the gradient of velocity y
        call atl_load_bc_state( bc          = bc_state_gradient(3), &
          &                     state_name  = 'grad_velocityY',     &
          &                     conf        = conf,                 &
          &                     bc_handle   = thandle,              &
          &                     varDict     = bc_varDict_gradient,  &
          &                     varSys      = equation%varSys       )

        if (nDims > 2) then
          ! Impose velocity z
          call atl_load_bc_state( bc          = bc_state(4),    &
            &                     state_name  = 'velocityZ',    &
            &                     conf        = conf,           &
            &                     bc_handle   = thandle,        &
            &                     varDict     = bc_varDict,     &
            &                     varSys      = equation%varSys )
          ! and the gradient of velocity z
          call atl_load_bc_state( bc          = bc_state_gradient(4), &
            &                     state_name  = 'grad_velocityZ',     &
            &                     conf        = conf,                 &
            &                     bc_handle   = thandle,              &
            &                     varDict     = bc_varDict_gradient,  &
            &                     varSys      = equation%varSys       )
        end if
      end if

      ! Impose pressure
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varDict     = bc_varDict,       &
        &                     varSys      = equation%varSys   )
      ! and the gradient for pressue
      call atl_load_bc_state( bc          = bc_state_gradient(pIndex), &
        &                     state_name  = 'grad_pressure',           &
        &                     conf        = conf,                      &
        &                     bc_handle   = thandle,                   &
        &                     varDict     = bc_varDict_gradient,       &
        &                     varSys      = equation%varSys            )


      ! check if all state variables are set
      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition grad_primitives you have to'
        write(logUnit(1),*) 'set all primitive variables (density, velocityX' &
          &                 //'velocityY, velocityZ'
        write(logUnit(1),*) 'and pressure) this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

      ! check if all gradients of variables are set
      if (.not. all(bc_state_gradient(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition grad primitives you have to'
        write(logUnit(1),*) 'all gradients of primitive variables:' &
          &                 //'grad_density, grad_velocityX, grad_velocityY,' &
          &                 // 'grad_velocityZ, grad_pressue'
        write(logUnit(1),*) 'this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

!!VK      ! add the conservative ones to varDict
!!VK      do iVar=1,pIndex
!!VK        kvp%key = trim(bc_state_gradient(iVar)%state_name)
!!VK        call append( me = bc_varDict_gradient, val = kvp )
!!VK      end do

    case('gradients')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2prim
        bc_trafo%from => atl_eqn_euler_2d_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons_grad
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons_grad
      end select
      bc_trafo%identity = .false.
      bc_trafo_gradient%identity = .false.

      bc_normal_vec_gradient = .false.
      bc_normal_vec = .false.

      ! set the all primtive variables to neumann
      bc_state(:)%style = 'neumann'
      ! set the bc style for all gradients to dirchlet
      bc_state_gradient(:)%style = 'dirichlet'

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      !impose the gradient of density
      call atl_load_bc_state( bc          =  bc_state_gradient(1), &
        &                     state_name  = 'grad_density',        &
        &                     conf        = conf,                  &
        &                     bc_handle   = thandle,               &
        &                     varDict     = bc_varDict_gradient,   &
        &                     varSys      = equation%varSys        )

      ! Extrapolate vel X
      bc_state(2)%state_name = 'velocityX'
      bc_state(2)%style = 'neumann'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! impose and he gradient of velocity x
      call atl_load_bc_state( bc          = bc_state_gradient(2), &
        &                     state_name  = 'grad_velocityX',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict_gradient,  &
        &                     varSys      = equation%varSys       )

      if (nDims > 1) then
         ! Extrapolate vel Y
         bc_state(3)%state_name = 'velocityY'
         bc_state(3)%style = 'neumann'
         bc_state(3)%isDefined = .true.
         kvp%key = trim(bc_state(2)%state_name)
         call append( me = bc_varDict, val = kvp )

        ! impose the gradient of velocity y
        call atl_load_bc_state( bc          = bc_state_gradient(3), &
          &                     state_name  = 'grad_velocityY',     &
          &                     conf        = conf,                 &
          &                     bc_handle   = thandle,              &
          &                     varDict     = bc_varDict_gradient,  &
          &                     varSys      = equation%varSys       )

        if (nDims > 2) then
          ! Extrapolate vel Z
          bc_state(4)%state_name = 'velocityZ'
          bc_state(4)%style = 'neumann'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(2)%state_name)
          call append( me = bc_varDict, val = kvp )

          ! impose the gradient of velocity z
          call atl_load_bc_state( bc          = bc_state_gradient(4), &
            &                     state_name  = 'grad_velocityZ',     &
            &                     conf        = conf,                 &
            &                     bc_handle   = thandle,              &
            &                     varDict     = bc_varDict_gradient,  &
            &                     varSys      = equation%varSys       )
        end if
      end if

      ! Extrapolate pressure
      bc_state(pIndex)%state_name = 'pressure'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! impose the gradient for pressue
      call atl_load_bc_state( bc          = bc_state_gradient(pIndex), &
        &                     state_name  = 'grad_pressure',           &
        &                     conf        = conf,                      &
        &                     bc_handle   = thandle,                   &
        &                     varDict     = bc_varDict_gradient,       &
        &                     varSys      = equation%varSys            )

      ! check if all gradients of variables are set
      if (.not. all(bc_state_gradient(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition grad primitives you have to'
        write(logUnit(1),*) 'all gradients of primitive variables:' &
          &                 //'grad_density, grad_velocityX, grad_velocityY,' &
          &                 // 'grad_velocityZ, grad_pressue'
        write(logUnit(1),*) 'this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('inflow_normal')

      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2prim
        bc_trafo%from => atl_eqn_euler_2d_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons_grad
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons_grad
      end select
      bc_trafo%identity = .false.
      bc_normal_vec = .true.
!!VK
!!VK      ! For the viscous parts of the equation
!!VK      ! ... the face values
!!VK      bc_trafo_gradient(1)%identity = .false.
!!VK      bc_normal_vec_gradient(1) = .true.
      ! ... the gradients on the face
      bc_trafo_gradient%identity = .true.
      bc_normal_vec_gradient = .false.

      ! Impose density at inlet
      call atl_load_bc_state( bc          = bc_state(1),    &
        &                     state_name  = 'density',      &
        &                     style       = 'dirichlet',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varSys )

      ! ... it's gradient (density - extrapolate)
      bc_state_gradient(1)%state_name = 'density'
      bc_state_gradient(1)%style = 'neumann'
      bc_state_gradient(1)%isDefined = .true.

      ! Impose normal velocity
      call atl_load_bc_state( bc          = bc_state(2),    &
        &                     state_name  = 'v_norm',       &
        &                     style       = 'dirichlet',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varSys )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.

      if (nDims > 1) then
        ! Impose tangential velocity to zero
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'dirichlet'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.

        if (nDims > 2) then
          ! Impose tangential velocity
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'dirichlet'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )

          ! ... it's gradient (momentum - extrapolate)
          bc_state_gradient(4)%state_name = 'momentum_tan2'
          bc_state_gradient(4)%style = 'neumann'
          bc_state_gradient(4)%isDefined = .true.
        end if
      end if

      ! Extrapolate pressure
      bc_state(pIndex)%state_name = 'pressure'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! ... it's gradient (energy - extrapolate)
      bc_state_gradient(pIndex)%state_name = 'energy'
      bc_state_gradient(pIndex)%style = 'neumann'
      bc_state_gradient(pIndex)%isDefined = .true.

      ! set bc_state_gradient(:,1) to maintain order in bc_varDict_gradient
      ! i.e odd count refer to state_gradient(:,1) and
      ! even count refer to state_gradient(:,2)
      do iVar=1,pIndex
      !!VK  bc_state_gradient(iVar, 1) = bc_state_gradient(iVar,2)
      !!VK  kvp%key = trim(bc_state_gradient(iVar,1)%state_name)
      !!VK  call append( me = bc_varDict_gradient, val = kvp )
        kvp%key = trim(bc_state_gradient(iVar)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )
      end do

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition inflow_normal you have to'
        write(logUnit(1),*) 'set the primitive variables density, v_norm and'
        write(logUnit(1),*) 'v_tan, this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('supersonic_inflow_normal')

      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2prim
        bc_trafo%from => atl_eqn_euler_2d_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons_grad
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons_grad
      end select
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! For the viscous parts of the equation
      ! ... the face values
!!VK      bc_trafo_gradient(1)%identity = .false.
!!VK      bc_normal_vec_gradient(1) = .true.
      ! ... the gradients on the face
      bc_trafo_gradient%identity = .true.
      bc_normal_vec_gradient = .false.

      ! Impose density at inlet
      call atl_load_bc_state( bc          = bc_state(1),    &
        &                     state_name  = 'density',      &
        &                     style       = 'dirichlet',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varSys )

      ! ... it's gradient (density - extrapolate)
      bc_state_gradient(1)%state_name = 'density'
      bc_state_gradient(1)%style = 'neumann'
      bc_state_gradient(1)%isDefined = .true.

      ! Impose normal velocity
      call atl_load_bc_state( bc          = bc_state(2),    &
        &                     state_name  = 'v_norm',       &
        &                     style       = 'dirichlet',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varSys )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.

      if (nDims > 1) then
        ! Impose tangential velocity to zero
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'dirichlet'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.

        if (nDims > 2) then
          ! Impose tangential velocity to zero
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'dirichlet'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )

          ! ... it's gradient (momentum - extrapolate)
          bc_state_gradient(4)%state_name = 'momentum_tan2'
          bc_state_gradient(4)%style = 'neumann'
          bc_state_gradient(4)%isDefined = .true.
        end if
      end if

      ! Impose pressure
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     style       = 'dirichlet',      &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varDict     = bc_varDict,       &
        &                     varSys      = equation%varSys   )

      ! ... it's gradient (energy - extrapolate)
      bc_state_gradient(pIndex)%state_name = 'energy'
      bc_state_gradient(pIndex)%style = 'neumann'
      bc_state_gradient(pIndex)%isDefined = .true.

      ! set bc_state_gradient(:,1) to maintain order in bc_varDict_gradient
      ! i.e odd count refer to state_gradient(:,1) and
      ! even count refer to state_gradient(:,2)
      do iVar=1,pIndex
!!VK        bc_state_gradient(iVar, 1) = bc_state_gradient(iVar,2)
!!VK        kvp%key = trim(bc_state_gradient(iVar,1)%state_name)
!!VK        call append( me = bc_varDict_gradient, val = kvp )
        kvp%key = trim(bc_state_gradient(iVar)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )
      end do

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition inflow_normal you have to'
        write(logUnit(1),*) 'set the primitive variables density, v_norm and'
        write(logUnit(1),*) 'v_tan this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('outflow')

      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2prim
        bc_trafo%from => atl_eqn_euler_2d_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons_grad
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons_grad
      end select
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! For the viscous parts of the equation
      ! ... the face values
!!VK      bc_trafo_gradient(1)%identity = .false.
!!VK      bc_normal_vec_gradient(1) = .true.
      ! ... the gradients on the face
      bc_trafo_gradient%identity = .true.
      bc_normal_vec_gradient = .false.

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! ... it's gradient (density - extrapolate)
      bc_state_gradient(1) = bc_state(1)

      ! Extrapolate v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'neumann'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.

      if (nDims > 1) then
        ! Extrapolate v_tangential_1
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'neumann'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.

        if (nDims > 2) then
          ! Extrapolate v_tangential_1
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'neumann'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )

          ! ... it's gradient (momentum - extrapolate)
          bc_state_gradient(4)%state_name = 'momentum_tan2'
          bc_state_gradient(4)%style = 'neumann'
          bc_state_gradient(4)%isDefined = .true.
        end if
      end if

      ! Impose pressure
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     style       = 'dirichlet',      &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varDict     = bc_varDict,       &
        &                     varSys      = equation%varSys   )

      ! ... it's gradient (energy - extrapolate)
      bc_state_gradient(pIndex)%state_name = 'energy'
      bc_state_gradient(pIndex)%style = 'neumann'
      bc_state_gradient(pIndex)%isDefined = .true.

      do iVar=1,pIndex
!!VK        bc_state_gradient(iVar,1) = bc_state_gradient(iVar,2)
!!VK        kvp%key = trim(bc_state_gradient(iVar,1)%state_name)
!!VK        call append( me = bc_varDict_gradient, val = kvp )
        kvp%key = trim(bc_state_gradient(iVar)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )
      end do

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
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2prim
        bc_trafo%from => atl_eqn_euler_2d_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons_grad
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim_grad
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons_grad
      end select
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! For the viscous parts of the equation
      ! ... the face values
!!VK      bc_trafo_gradient(1)%identity = .false.
!!VK      bc_normal_vec_gradient(1) = .true.
      ! ... the gradients on the face
      bc_trafo_gradient%identity = .true.
      bc_normal_vec_gradient = .false.

      ! Extrapolate density
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! ... it's gradient (density - extrapolate)
      bc_state_gradient(1) = bc_state(1)

      ! Extrapolate v_normal
      bc_state(2)%state_name = 'v_norm'
      bc_state(2)%style = 'neumann'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.

      if (nDims > 1) then
        ! Extrapolate v_tangential_1
        bc_state(3)%state_name = 'v_tan'
        bc_state(3)%style = 'neumann'
        bc_state(3)%isDefined = .true.
        kvp%key = trim(bc_state(3)%state_name)
        call append( me = bc_varDict, val = kvp )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.

        if (nDims > 2) then
          ! Extrapolate v_tangential_2
          bc_state(4)%state_name = 'v_tan2'
          bc_state(4)%style = 'neumann'
          bc_state(4)%isDefined = .true.
          kvp%key = trim(bc_state(4)%state_name)
          call append( me = bc_varDict, val = kvp )

          ! ... it's gradient (momentum - extrapolate)
          bc_state_gradient(4)%state_name = 'momentum_tan2'
          bc_state_gradient(4)%style = 'neumann'
          bc_state_gradient(4)%isDefined = .true.
        end if
      end if

      ! Extrapolate pressure
      bc_state(pIndex)%state_name = 'pressure'
      bc_state(pIndex)%style = 'neumann'
      bc_state(pIndex)%isDefined = .true.
      kvp%key = trim(bc_state(pIndex)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! ... it's gradient (energy - extrapolate)
      bc_state_gradient(pIndex)%state_name = 'energy'
      bc_state_gradient(pIndex)%style = 'neumann'
      bc_state_gradient(pIndex)%isDefined = .true.

      do iVar=1,pIndex
!!VK        bc_state_gradient(iVar,1) = bc_state_gradient(iVar,2)
!!VK        kvp%key = trim(bc_state_gradient(iVar,1)%state_name)
!!VK        call append( me = bc_varDict_gradient, val = kvp )
        kvp%key = trim(bc_state_gradient(iVar)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )
      end do

    case default
      write(logUnit(1),*) 'Unknown boundary kind "' // trim(bc_kind) // '"'
      write(logUnit(1),*) 'for boundary  "' // trim(bc_label) // '".'
      write(logUnit(1),*) 'Available boundary kinds for Navier-Stokes equations:'
      write(logUnit(1),*) ' * wall'
      write(logUnit(1),*) ' * adiabatic_wall'
      write(logUnit(1),*) ' * inflow_normal'
      write(logUnit(1),*) ' * supersonic_inflow_normal'
      write(logUnit(1),*) ' * outflow'
      write(logUnit(1),*) ' * supersonic_outflow'
      write(logUnit(1),*) ' * primitives'
      write(logUnit(1),*) ' * grad_primitives'
      write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
      call tem_abort()
    end select

    call truncate( me = bc_varDict )
    call truncate( me = bc_varDict_gradient )

    if (size(bc_state) /= bc_varDict%nVals) then
      write(logUnit(1),*) 'Nr. of state variables does not match size of '//&
        &                 'varDict'
      call tem_abort()
    end if

    if (size(bc_state_gradient) /= bc_varDict_gradient%nVals) then
      write(*,*) 'bc_state_gradient ', bc_state_gradient(:)%state_name
      write(*,*) 'bc_varDict_gradient ', bc_varDict_gradient%val
      write(logUnit(1),*) 'Nr. of state gradient variables does not match '//&
        &                 'size of varDict_gradient'
      call tem_abort()
    end if

  end subroutine atl_eqn_nvrstk_load_bc

  ! ------------------------------------------------------------------------ !
  !> Estimate the impact of viscous terms in 3D.
  !!
  !! We use an estimate for mu*(dv/dx) from the deviations and derivative
  !! estimates of the conservative variables:
  !! dm/dx = dv/dx * rho + drho/dx * v
  !! dv/dx = (dm/dx - drho/dx * v)/rho
  !! dv/dx < (max(dm/dx) - max(drho/dx) * max(v))/min(rho)
  !! dv/dx < (max(dm/dx) * min(rho) - max(drho/dx) * max(m) ) / min(rho)**2
  pure function inviscous_indicator_3d(nvrstk, mean, deviation, grad) &
    &   result(isinviscous)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_navierstokes_type), intent(in) :: nvrstk

    !> The mean of each state variable.
    real(kind=rk), intent(in) :: mean(:)

    !> Estimation of maximal deviation of each state.
    real(kind=rk), intent(in) :: deviation(:)

    !> Estimation of maximal gradient of each state.
    real(kind=rk), intent(in) :: grad(:)

    !> Resulting indication whether viscous terms can be neglected.
    logical :: isinviscous
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rho_min
    real(kind=rk) :: m_max
    real(kind=rk) :: grad_mag
    ! -------------------------------------------------------------------- !

    rho_min = max(mean(1) - deviation(1), 0.0_rk)
    m_max = sqrt(  (abs(mean(2)) + deviation(2))**2 &
      &          + (abs(mean(3)) + deviation(3))**2 &
      &          + (abs(mean(4)) + deviation(4))**2 )
    grad_mag = sqrt(grad(2)**2 + grad(3)**2 + grad(4)**2)
    isinviscous = nvrstk%mu * (grad_mag * rho_min + grad(1) * m_max) &
      &           < nvrstk%visc_limit*rho_min**2

  end function inviscous_indicator_3d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Estimate the impact of viscous terms in 2D.
  !!
  !! We use an estimate for mu*(dv/dx) from the deviations and derivative
  !! estimates of the conservative variables:
  !! dm/dx = dv/dx * rho + drho/dx * v
  !! dv/dx = (dm/dx - drho/dx * v)/rho
  !! dv/dx < (max(dm/dx) - max(drho/dx) * max(v))/min(rho)
  !! dv/dx < (max(dm/dx) * min(rho) - max(drho/dx) * max(m) ) / min(rho)**2
  pure function inviscous_indicator_2d(nvrstk, mean, deviation, grad) &
    &   result(isinviscous)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_navierstokes_type), intent(in) :: nvrstk

    !> The mean of each state variable.
    real(kind=rk), intent(in) :: mean(:)

    !> Estimation of maximal deviation of each state.
    real(kind=rk), intent(in) :: deviation(:)

    !> Estimation of maximal gradient of each state.
    real(kind=rk), intent(in) :: grad(:)

    !> Resulting indication whether viscous terms can be neglected.
    logical :: isinviscous
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rho_min
    real(kind=rk) :: m_max
    real(kind=rk) :: grad_mag
    ! -------------------------------------------------------------------- !

    rho_min = max(mean(1) - deviation(1), 0.0_rk)
    m_max = sqrt(  (abs(mean(2)) + deviation(2))**2 &
      &          + (abs(mean(3)) + deviation(3))**2 )
    grad_mag = sqrt(grad(2)**2 + grad(3)**2)
    isinviscous = nvrstk%mu * (grad_mag * rho_min + grad(1) * m_max) &
      &           < nvrstk%visc_limit*rho_min**2

  end function inviscous_indicator_2d
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Deactivate adaptive inviscous computations.
  pure function inviscous_deactivated(nvrstk, mean, deviation, grad) &
    &   result(isinviscous)
    ! -------------------------------------------------------------------- !
    !> Description of the equation
    class(atl_navierstokes_type), intent(in) :: nvrstk

    !> The mean of each state variable.
    real(kind=rk), intent(in) :: mean(:)

    !> Estimation of maximal deviation of each state.
    real(kind=rk), intent(in) :: deviation(:)

    !> Estimation of maximal gradient of each state.
    real(kind=rk), intent(in) :: grad(:)

    !> Resulting indication whether viscous terms can be neglected.
    logical :: isinviscous
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    isinviscous = .false.
  end function inviscous_deactivated

end module atl_eqn_nvrstk_hlp_module
