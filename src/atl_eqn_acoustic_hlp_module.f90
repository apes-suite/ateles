! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> Helper routines for the acoustic equation system.
module atl_eqn_acoustic_hlp_module
  use aotus_module,                   only: flu_State

  use tem_aux_module,                 only: tem_abort
  use tem_bc_module,                  only: tem_bc_state_type
  use tem_logging_module,             only: logUnit
  use tem_stringKeyValuePair_module,  only: tem_stringKeyValuePair_type
  use tem_stringKeyValuePair_module,  only: grw_stringKeyValuePairArray_type,  &
    &                                       init, truncate, append


  use atl_equation_module,            only: atl_equations_type,     &
    &                                       atl_eqn_var_trafo_type
  use atl_bc_state_module,            only: atl_load_bc_state

  use atl_eqn_acoustic_module,        only: atl_load_acoustic
  use atl_eqn_acoustic_var_module,    only: atl_init_acoustic_vars,        &
    &                                       atl_init_acoustic_sourceTerms
  use atl_eqn_acoustic_2d_var_module, only: atl_init_acoustic_2d_vars,&
    &                                       atl_init_acoustic_2d_sourceTerms
  use atl_varSys_module,              only: atl_varSys_solverData_type
  use atl_source_types_module,        only: atl_init_source_type

  implicit none

  private

  public :: atl_eqn_acoustic_load_bc
  public :: atl_eqn_acoustic_init


contains


  !> Initialization of the Acoustic equation.
  !!
  !! This routine sets up the necessary infrastructure for the Acoustic
  !! equations.
  !! It reads the configuration from the given script in conf under the table
  !! provided in thandle and sets function pointers and variables accordingly.
  subroutine atl_eqn_acoustic_init( conf, thandle, equation, nDimensions, &
    &                               initSource, varSys_data )
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

    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), intent(inout) :: varSys_data
    ! ---------------------------------------------------------------------------

    equation%isNonlinear = .false.
    equation%nDimensions = nDimensions
    ! timesetp is static and not changing over simulationtime
    equation%adaptive_timestep = .false.
    ! define the dimension of the acoustic datatype, required for different
    ! size of background velocity array
    equation%acoustic%ndims = nDimensions

    equation%load_bc => atl_eqn_acoustic_load_bc

    select case(nDimensions)
    case(2)
      call atl_init_acoustic_2d_vars( equation   = equation,   &
        &                             methodData = varSys_data )

      call atl_load_acoustic( acoustic     = equation%acoustic, &
        &                     conf         = conf,              &
        &                     eq_table     = thandle            )

      call atl_init_acoustic_2d_sourceTerms( initSource%poss_srcVars, &
        &                                    initSource%eval_source   )

    case(3)
      call atl_init_acoustic_vars( equation   = equation,   &
        &                          methodData = varSys_data )

      call atl_load_acoustic( acoustic     = equation%acoustic, &
        &                     conf         = conf,              &
        &                     eq_table     = thandle            )

      call atl_init_acoustic_sourceTerms( initSource%poss_srcVars, &
        &                                 initSource%eval_source   )

    end select

  end subroutine atl_eqn_acoustic_init


  ! *****************************************************************************
  !> Reading boundary conditions for the acoustic equations.
  !!
  !! Need to set 4 bc_states here, typically the primitive variables.
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
  subroutine atl_eqn_acoustic_load_bc( equation,                              &
    &                                  bc_state, bc_state_gradient,           &
    &                                  bc_varDict, bc_varDict_gradient,       &
    &                                  bc_normal_vec, bc_normal_vec_gradient, &
    &                                  bc_trafo, bc_trafo_gradient,           &
    &                                  bc_label, bc_kind, thandle,            &
    &                                  conf                                   )
    ! ---------------------------------------------------------------------------
    class(atl_equations_type), intent(inout) :: equation
    type(tem_bc_state_type), allocatable, intent(out) :: bc_state(:)
    type(tem_bc_state_type), allocatable, intent(out) :: bc_state_gradient(:)
    !> Dictionary of boundary variables in bc_state
    type(grw_stringKeyValuePairArray_type), intent(out) :: bc_varDict
    !> Dictionary of boundary variables in bc_state_gradient
    type(grw_stringKeyValuePairArray_type), intent(out) :: bc_varDict_gradient
    logical, intent(out)                              :: bc_normal_vec
    logical, intent(out)                              :: bc_normal_vec_gradient
    character(len=*), intent(in)                      :: bc_label
    character(len=*), intent(in)                      :: bc_kind
    type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo
    type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo_gradient
    integer, intent(in)                               :: thandle
    type(flu_State)                                   :: conf
    ! ---------------------------------------------------------------------------
    integer :: nDims
    type(tem_stringKeyValuePair_type) :: kvp
    ! ---------------------------------------------------------------------------

    nDims = equation%nDimensions

    ! The acoustic equation requires nDims + 1 variables to be set:
    allocate(bc_state(equation%varSys%nScalars))
    allocate(bc_state_gradient(0))
    bc_normal_vec_gradient = .false.
!!VK    allocate(bc_normal_vec_gradient(2))
!!VK    allocate(bc_trafo_gradient(2))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )
    ! Constant zero variable for non-configurable boundary variable
    kvp%value = 'zero_const'

    ! For acoustic equation primitive and conservative variables are the same,
    ! thus we set the transfamation for all boundary condition
    bc_trafo%identity = .true.

    ! The bc_trafo function pointer %from and % to are intilized to NULL and
    ! since it is not required for acoustic equation there are not set to any
    ! function here!

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

    case('conservatives')
      bc_normal_vec = .false.

      ! Impose density
      call atl_load_bc_state( bc          = bc_state(1),     &
        &                     state_name  = 'density',       &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      ! Impose momentum x
      call atl_load_bc_state( bc          = bc_state(2),     &
        &                     state_name  = 'velocityX',     &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varSys      = equation%varSys, &
        &                     varDict     = bc_varDict       )

      if (nDims > 1) then
        ! Impose momentum y
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

      if (.not. all(bc_state(:)%isDefined)) then
         write(logUnit(1),*) 'For boundary condition conservatives you have to'
         write(logUnit(1),*) 'set all conservative variables (density, '
         write(logUnit(1),*) 'velocityX, velocityY, velocityZ) this set is not'
         write(logUnit(1),*) ' complete for ' // trim(bc_label) // '!'
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

    case default
       write(logUnit(1),*) 'Unknown boundary kind "' // trim(bc_kind) // '"'
       write(logUnit(1),*) 'for boundary  "' // trim(bc_label) // '".'
       write(logUnit(1),*) 'Available boundary kinds for acoustic equations:'
       write(logUnit(1),*) ' * slipwall / wall'
       write(logUnit(1),*) ' * conservatives'
       write(logUnit(1),*) ' * inflow and inflow_normal'
       write(logUnit(1),*) ' * outflow'
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
  end subroutine atl_eqn_acoustic_load_bc
  ! ****************************************************************************

end module atl_eqn_acoustic_hlp_module
