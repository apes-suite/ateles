! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013, 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> Helper routines for the advection equation system in 1D.
module atl_eqn_advection_1d_hlp_module
  use aotus_module,                    only: flu_State

  use tem_aux_module,                  only: tem_abort
  use tem_bc_module,                   only: tem_bc_state_type
  use tem_logging_module,              only: logUnit
  use tem_stringKeyValuePair_module,   only: tem_stringKeyValuePair_type, &
    &                                        init, truncate, append,      &
    &                                        grw_stringKeyValuePairArray_type

  use atl_equation_module,             only: atl_equations_type, &
    &                                        atl_eqn_var_trafo_type
  use atl_bc_state_module,             only: atl_load_bc_state
  use atl_eqn_advection_1d_module,     only: atl_load_advection_1d
  use atl_eqn_advection_1d_var_module, only: atl_init_advection_1d_vars
  use atl_varSys_module,               only: atl_varSys_solverData_type

  implicit none

  private

  public :: atl_eqn_advection_1d_load_bc
  public :: atl_eqn_advection_1d_init

contains

  !> Initialization of the linearized Euler equations.
  !!
  !! This routine sets up the necessary infrastructure for the linearized
  !! Euler equations.
  !! It reads the configuration from the given script in conf under the table
  !! provided in thandle and sets function pointers and variables accordingly.
  subroutine atl_eqn_advection_1d_init( conf, thandle, equation, nDimensions, &
    &                                   varSys_data )
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

    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), intent(inout) :: varSys_data
    ! ---------------------------------------------------------------------------

    equation%nDimensions = nDimensions
    equation%isNonlinear = .false.
    ! timesetp is static and not changing over simulationtime
    equation%adaptive_timestep = .false.
    equation%load_bc => atl_eqn_advection_1d_load_bc

    call atl_load_advection_1d(equation%advection, conf, thandle)

    call atl_init_advection_1d_vars(                            &
      & varSys                = equation%varSys,                &
      & hasPrimitiveVariables = equation%hasPrimitiveVariables, &
      & methodData            = varSys_data                     )

  end subroutine atl_eqn_advection_1d_init


  !> Reading boundary conditions for the advection equations in 1D.
  !!
  !! This routine has to conform to the interface definition
  !! atl_equation_module#eqn_load_bc.
  subroutine atl_eqn_advection_1d_load_bc( equation, bc_state,                &
    & bc_state_gradient, bc_varDict, bc_varDict_gradient, bc_normal_vec,      &
    & bc_normal_vec_gradient, bc_trafo, bc_trafo_gradient, bc_label, bc_kind, &
    & thandle, conf                                                           )
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
    type(tem_stringKeyValuePair_type) :: kvp
    ! ---------------------------------------------------------------------------

    ! 1D advection requires 1 variables to be set:
    allocate(bc_state(1))
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
    ! even if the boundary conditions does not use them.
    bc_trafo%from => null()
    bc_trafo%to => null()

    select case(bc_kind)
    case('inflow')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .true.
      bc_normal_vec = .false.

      ! Impose density
      call atl_load_bc_state( bc          = bc_state(1),     &
        &                     state_name  = 'density',       &
        &                     nComp       = 1,               &
        &                     conf        = conf,            &
        &                     bc_handle   = thandle,         &
        &                     varsys      = equation%varsys, &
        &                     varDict     = bc_varDict       )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition inflow you have to set the'
        write(logUnit(1),*) 'primitive variable density '
        write(logUnit(1),*) 'this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case('outflow')
      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      bc_trafo%identity = .true.

      ! Extrapolate density
      bc_normal_vec = .true.
      bc_state(1)%state_name = 'density'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

    case default
      write(logUnit(1),*) 'Unknown boundary kind "' // trim(bc_kind) // '"'
      write(logUnit(1),*) 'for boundary  "' // trim(bc_label) // '".'
      write(logUnit(1),*) 'Available boundary kinds for advection 1D equations:'
      write(logUnit(1),*) ' * inflow'
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
  end subroutine atl_eqn_advection_1d_load_bc

end module atl_eqn_advection_1d_hlp_module
