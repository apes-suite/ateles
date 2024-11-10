! Copyright (c) 2014-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2016, 2018, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018-2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> Helper routines for the LinearEuler equation system.
module atl_eqn_LinearEuler_hlp_module
  use env_module,                        only: labelLen
  use aotus_module,                      only: flu_State, aot_get_val

  use tem_aux_module,                    only: tem_abort
  use tem_tools_module,                  only: upper_to_lower
  use tem_bc_module,                     only: tem_bc_state_type
  use tem_logging_module,                only: logUnit
  use tem_time_module,                   only: tem_time_type
  use tem_stringKeyValuePair_module,     only: tem_stringKeyValuePair_type, &
    &                                          init, truncate, append,      &
    &                                          grw_stringKeyValuePairArray_type


  use atl_equation_module,               only: atl_equations_type, &
    &                                          atl_eqn_var_trafo_type
  use atl_bc_state_module,               only: atl_load_bc_state
  use atl_eqn_LinearEuler_module,        only: atl_load_LinearEuler,     &
    &                                          atl_linearEuler_type,     &
    &                                          atl_linEuler_numflux,     &
    &                                          atl_eqn_update_background
  use atl_eqn_LinearEuler_var_module,    only: atl_init_LinearEuler_vars, &
    &                                          atl_init_lineuler_sourceTerms
  use atl_eqn_LinearEuler_2d_var_module, only: atl_init_LinearEuler_2d_vars
  use atl_varSys_module,                 only: atl_varSys_solverData_type
  use atl_source_types_module,           only: atl_init_source_type
  use atl_linearEuler_numFlux_module, only: &
    &   atl_LinearEuler_numFlux_subleft,    &
    &   atl_LinearEuler_numFlux_subright,   &
    &   atl_LinearEuler_numFlux_superleft,  &
    &   atl_LinearEuler_numFlux_superright
  use atl_linearEuler_2d_numFlux_module, only: &
    &   atl_LinearEuler_2d_numFlux_subleft,   &
    &   atl_LinearEuler_2d_numFlux_subright,  &
    &   atl_LinearEuler_2d_numFlux_superleft, &
    &   atl_LinearEuler_2d_numFlux_superright
  use atl_laxFriedrichFlux_module,       only: atl_laxFriedLinearEuler
  use atl_laxFriedrichFlux_2d_module,    only: atl_laxFriedLinearEuler_2d

  implicit none

  private

  public :: atl_eqn_LinearEuler_load_bc
  public :: atl_eqn_linearEuler_init
  public :: atl_getLinearEulerFluxes


contains


  ! ------------------------------------------------------------------------ !
  !> Initialization of the linearized Euler equations.
  !!
  !! This routine sets up the necessary infrastructure for the linearized
  !! Euler equations.
  !! It reads the configuration from the given script in conf under the table
  !! provided in thandle and sets function pointers and variables accordingly.
  subroutine atl_eqn_linearEuler_init( conf, thandle, equation, nDimensions, &
    &                                  varSys_data, initSource               )
    ! -------------------------------------------------------------------- !
    !> Handle to the Lua configuration
    type(flu_State), intent(in) :: conf

    !> Handle to the equation table in the Lua script given in conf.
    integer, intent(in) :: thandle
    !> Equation system to set with this routine.
    type(atl_equations_type), intent(inout) :: equation
    !> Type to be filled with the possible source variables for the equation
    !! system. These source variables are later on used to extract the
    !! corresponding information from the configuration file.
    type(atl_init_source_type), intent(inout) :: initSource
    !> Number of spatial dimensions, the Euler equations should live on.
    !!
    !! Has to be 1, 2 or 3.
    integer, intent(in) :: nDimensions
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), intent(inout) :: varSys_data
    ! -------------------------------------------------------------------- !
    !> local type to initial background for the first calculation of cfl
    !! timestep
    type(tem_time_type) :: init_time
    ! -------------------------------------------------------------------- !

    equation%isNonlinear = .false.
    equation%nDimensions = nDimensions
    equation%nDerivatives = 0
    ! timestep is static and not changing over simulationtime
    equation%adaptive_timestep = .false.

    call atl_load_LinearEuler( linearEuler  = equation%LinearEuler, &
      &                        conf         = conf,                 &
      &                        eq_table     = thandle,              &
      &                        spatial_dim  = nDimensions           )

    select case(nDimensions)
    case(3)
      equation%load_bc => atl_eqn_LinearEuler_load_bc
      call atl_init_LinearEuler_vars( equation   = equation,   &
        &                             methodData = varSys_data )
      call atl_init_lineuler_sourceTerms( initSource%poss_srcVars, &
        &                                 initSource%eval_source   )
    case (2)
      equation%load_bc => atl_eqn_LinearEuler_load_bc
      call atl_init_LinearEuler_2d_vars( equation   = equation,   &
        &                                methodData = varSys_data )
      call atl_init_lineuler_sourceTerms( initSource%poss_srcVars, &
        &                                 initSource%eval_source   )
    end select

    !> ToDo: no source terms implemented yet
    !call atl_init_LinearEuler_sourceTerms( initSource%poss_srcVars, &
    !  &                                    initSource%eval_source   )

    ! NO material for linearized euler
    !call atl_init_linearEuler_material(               &
    !  & possVars    = initMaterial%poss_materialVars, &
    !  & nDimensions = nDimensions                     )

    ! update the background for the first timestep calulation
    ! for that we build up the local tem_time_type inittime which is required
    ! in updating the background
    init_time%iter = 0
    init_time%sim  = 0.0
    call atl_eqn_update_background(         &
      & me          = equation%linearEuler, &
      & time        = init_time,            &
      & nDimensions = nDimensions           )
     write(logUnit(1),*) 'Loaded linearized Euler equation'
     write(logUnit(1),*) ' * isen_coef: ', &
       &                 equation%linearEuler%isen_coef
     write(logUnit(1),*) ' * background density: ', &
       &                 equation%linearEuler%density_0
     write(logUnit(1),*) ' * background velocityX: ', &
       &                 equation%linearEuler%velocity_0(1)
     write(logUnit(1),*) ' * background velocityY: ', &
       &                 equation%linearEuler%velocity_0(2)
     if ( nDimensions == 3) then
       write(logUnit(1),*) ' * background velocityZ: ', &
         &                 equation%linearEuler%velocity_0(3)
     end if
     write(logUnit(1),*) ' * background pressure: ', &
       &                 equation%linearEuler%pressure_0
     write(logUnit(1),*) ' * speed of sound: ', &
       &                 equation%linearEuler%speedofsound

    ! Getting the numerical flux per direction !
    call atl_getLinearEulerFluxes( linearEuler = equation%linearEuler, &
      &                            conf        = conf,                 &
      &                            eqn_handle  = thandle,              &
      &                            eqn_dim     = nDimensions           )

  end subroutine atl_eqn_linearEuler_init
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Reading boundary conditions for the LinearEuler equations.
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
  subroutine atl_eqn_LinearEuler_load_bc( equation, bc_state,                 &
    & bc_state_gradient, bc_varDict, bc_varDict_gradient, bc_normal_vec,      &
    & bc_normal_vec_gradient, bc_trafo, bc_trafo_gradient, bc_label, bc_kind, &
    & thandle, conf                                                           )
    ! -------------------------------------------------------------------- !
    class(atl_equations_type), intent(inout) :: equation
    type(tem_bc_state_type), allocatable, intent(out) :: bc_state(:)
    type(tem_bc_state_type), allocatable, intent(out) :: bc_state_gradient(:)
    !> Dictionary of boundary variables in bc_state
    type(grw_stringKeyValuePairArray_type), intent(out) :: bc_varDict
    !> Dictionary of boundary variables in bc_state_gradient
    type(grw_stringKeyValuePairArray_type), intent(out) :: bc_varDict_gradient
    logical, intent(out) :: bc_normal_vec
    logical, intent(out) :: bc_normal_vec_gradient
    character(len=*), intent(in) :: bc_label
    character(len=*), intent(in) :: bc_kind
    type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo
    type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo_gradient
    integer, intent(in) :: thandle
    type(flu_State) :: conf
    ! -------------------------------------------------------------------- !
    integer :: nDims
    integer :: pIndex
    type(tem_stringKeyValuePair_type) :: kvp
    ! -------------------------------------------------------------------- !

    nDims = equation%nDimensions
    pIndex = equation%varSys%nScalars
    allocate(bc_state(pIndex))
    allocate(bc_state_gradient(0))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )
    ! Constant zero variable for non-configurable boundary variable
    kvp%value = 'zero_const'

    ! For LinearEuler equation primitive and conservative variables are the
    ! same, thus we set the transfamation for all boundary condition
    bc_trafo%identity = .true.

    ! The bc_trafo function pointer %from and % to are intilized to NULL and
    ! since it is not required for LinearEuler equation there are not set to any
    ! function here!
    bc_trafo%from => null()
    bc_trafo%to   => null()

    select case(bc_kind)

    case('slipwall', 'wall')
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

    case('primitives')
      bc_normal_vec = .false.

      call atl_load_bc_state( bc          = bc_state(1),    &
        &                     state_name  = 'density',      &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varsys )

      call atl_load_bc_state( bc          = bc_state(2),    &
        &                     state_name  = 'velocityX',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varsys )

      if (nDims > 1) then
        call atl_load_bc_state( bc          = bc_state(3),    &
          &                     state_name  = 'velocityY',    &
          &                     conf        = conf,           &
          &                     bc_handle   = thandle,        &
          &                     varDict     = bc_varDict,     &
          &                     varSys      = equation%varsys )

        if (nDims > 2) then
          call atl_load_bc_state( bc          = bc_state(4),    &
            &                     state_name  = 'velocityZ',    &
            &                     conf        = conf,           &
            &                     bc_handle   = thandle,        &
            &                     varDict     = bc_varDict,     &
            &                     varSys      = equation%varsys )
        end if
      end if

      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varDict     = bc_varDict,       &
        &                     varSys      = equation%varsys   )

      if (.not. all(bc_state(:)%isDefined)) then
         write(logUnit(1),*) 'For boundary condition primtivies you '//&
           &                 'have to set'
         write(logUnit(1),*) 'all primitive variables (density, velocityX,'
         write(logUnit(1),*) 'velocityY, velocityZ, pressure) this set is not'
         write(logUnit(1),*) ' complete for ' // trim(bc_label) // '!'
         write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
         call tem_abort()
      end if

    case('inflow')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_normal_vec = .false.

      ! Impose denisty
      call atl_load_bc_state( bc          = bc_state(1),    &
        &                     state_name  = 'density',      &
        &                     style       = 'dirichlet',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varsys )

      ! Impose x velocity
      call atl_load_bc_state( bc          = bc_state(2),    &
        &                     state_name  = 'velocityX',    &
        &                     style       = 'dirichlet',    &
        &                     conf        = conf,           &
        &                     bc_handle   = thandle,        &
        &                     varDict     = bc_varDict,     &
        &                     varSys      = equation%varsys )

      if (nDims > 1) then
        ! Impose y velocity
        call atl_load_bc_state( bc          = bc_state(3),    &
          &                     state_name  = 'velocityY',    &
          &                     style       = 'dirichlet',    &
          &                     conf        = conf,           &
          &                     bc_handle   = thandle,        &
          &                     varDict     = bc_varDict,     &
          &                     varSys      = equation%varsys )

        if (nDims > 2) then
        ! Impose z velocity
        call atl_load_bc_state( bc          = bc_state(4),    &
          &                     state_name  = 'velocityZ',    &
          &                     style       = 'dirichlet',    &
          &                     conf        = conf,           &
          &                     bc_handle   = thandle,        &
          &                     varDict     = bc_varDict,     &
          &                     varSys      = equation%varsys )
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
        write(logUnit(1),*) 'primitive variables density, velocityX, '//&
          &                 'velocityY, velocityZ'
        write(logUnit(1),*) 'this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('outflow')
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
        &                     varDict     = bc_varDict,       &
        &                     varSys      = equation%varsys   )

      if (.not. bc_state(pIndex)%isDefined) then
        write(logUnit(1),*) 'For boundary condition outflow you have to set the'
        write(logUnit(1),*) 'pressure!'
        write(logUnit(1),*) 'Something is wrong with that in boundary ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if
    case('zero_gradient')
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
       write(logUnit(1),*) 'Available boundary kinds for LinearEuler equations:'
       write(logUnit(1),*) ' * slipwall / wall'
       write(logUnit(1),*) ' * primitives'
       write(logUnit(1),*) ' * outflow'
       write(logUnit(1),*) ' * inflow'
       write(logUnit(1),*) ' * zero_gradient'
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
      write(logUnit(1),*) 'Nr. of state gradient variables does not match '//&
        &                 'size of varDict_gradient'
      call tem_abort()
    end if

  end subroutine atl_eqn_LinearEuler_load_bc
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! Getting the numerical flux for LinearEuler equations
  ! depending on background state and speedofsound
  subroutine atl_getLinearEulerFluxes(Lineareuler,  conf, eqn_handle, eqn_dim)
    ! -------------------------------------------------------------------- !
    !> The equations type to set the numerical flux in.
    type(atl_linearEuler_type), intent(inout) :: LinearEuler
    !> Configuration file handle to get the numerical flux setting from.
    type(flu_state), intent(in) :: conf
    !> Handle to the equation table in the configuration script.
    integer, intent(in) :: eqn_handle
    !> Dimension of the equation to set the flux for.
    integer, intent(in) :: eqn_dim
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: eq_nflux
    procedure(atl_lineuler_numflux), pointer :: subleft    => NULL()
    procedure(atl_lineuler_numflux), pointer :: superleft  => NULL()
    procedure(atl_lineuler_numflux), pointer :: subright   => NULL()
    procedure(atl_lineuler_numflux), pointer :: superright => NULL()
    integer :: iError
    integer :: iDir
    ! -------------------------------------------------------------------- !

    ! -- possible fluxes: godunov, lax_friedrich
    call aot_get_val( L       = conf,       &
      &               thandle = eqn_handle, &
      &               key     = 'numflux',  &
      &               val     =  eq_nflux,  &
      &               ErrCode = iError,     &
      &               default = 'godunov'   )

    eq_nflux = upper_to_lower(eq_nflux)
    eq_nflux = adjustl(eq_nflux)

    select case(trim(eq_nflux))
    case ('godunov')
      write(logunit(2),*) 'Using godunov numerical flux of linear euler'
      select case(eqn_dim)
      case(1)
        write(*,*) 'Linear Euler for 1d not yet implemented! Stopping'
        call tem_abort()
      case(2)
        subleft    => atl_LinearEuler_2d_numFlux_subleft
        superleft  => atl_LinearEuler_2d_numFlux_superleft
        subright   => atl_LinearEuler_2d_numFlux_subright
        superright => atl_LinearEuler_2d_numFlux_superright
      case(3)
        subleft    => atl_LinearEuler_numFlux_subleft
        superleft  => atl_LinearEuler_numFlux_superleft
        subright   => atl_LinearEuler_numFlux_subright
        superright => atl_LinearEuler_numFlux_superright
      end select

      do iDir = 1, eqn_dim

        ! For the linearized Euler flux, there may be for different states in
        ! the Riemann problem at the interface.
        ! Which one to choose depends on the velocity normal to the interface.
        if ( LinearEuler%velocity_0(iDir) >= 0) then

          ! Velocity to the right (positive), use state from left side.
          if ( LinearEuler%velocity_0(iDir) < LinearEuler%SpeedOfSound ) then
            ! Subsonic
            LinearEuler%dir_proc(iDir)%numflux => subleft
          else
            ! Supersonic
            LinearEuler%dir_proc(iDir)%numflux => superleft
          end if

        else

          ! Velocity to the left (negative), use state from right side.
          if ( abs(LinearEuler%velocity_0(iDir)) &
            &  < LinearEuler%SpeedOfSound        ) then
            ! Subsonic
            LinearEuler%dir_proc(iDir)%numflux => subright
          else
            ! Supersonic
            LinearEuler%dir_proc(iDir)%numflux => superright
          end if

        end if

      end do

    case ('lax_friedrich')
      write(logunit(2),*) 'Using Lax Friedrichs numerical flux.'
      select case(eqn_dim)
      case(1)
        write(*,*) 'Linear Euler for 1d not yet implemented! Stopping'
        call tem_abort()
      case(2)
        LinearEuler%dir_proc(1)%numflux => atl_laxFriedLinearEuler_2D
        LinearEuler%dir_proc(2)%numflux => atl_laxFriedLinearEuler_2D
      case(3)
        LinearEuler%dir_proc(1)%numflux => atl_laxFriedLinearEuler
        LinearEuler%dir_proc(2)%numflux => atl_laxFriedLinearEuler
        LinearEuler%dir_proc(3)%numflux => atl_laxFriedLinearEuler
      end select
    case default
      write(logunit(1),*) 'Unknown numerical flux ', trim(eq_nflux)
      write(logunit(1),*) 'for the Linear Euler equation system.'
      write(logunit(1),*) 'Please choose one of the available:'
      write(logunit(1),*) ' * lax_friedrich (default)'
      write(logunit(1),*) ' * godunov'
      call tem_abort()
    end select
  end subroutine atl_getLinearEulerFluxes
  ! ------------------------------------------------------------------------ !

end module atl_eqn_LinearEuler_hlp_module
