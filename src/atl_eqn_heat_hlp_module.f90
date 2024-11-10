! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015-2016 Verena Krupp <verena.krupp@uni-siegen.de>
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

!> Helper routines for the Heat equation system.
module atl_eqn_heat_hlp_module
  use aotus_module,                   only: flu_State

  use tem_aux_module,                 only: tem_abort
  use tem_bc_module,                  only: tem_bc_state_type
  use tem_logging_module,             only: logUnit
  use tem_stringKeyValuePair_module,  only: tem_stringKeyValuePair_type
  use tem_stringKeyValuePair_module,  only: grw_stringKeyValuePairArray_type, &
    &                                       init, truncate, append

  use atl_equation_module,            only: atl_equations_type, &
    &                                       atl_eqn_var_trafo_type
  use atl_eqn_heat_module,            only: atl_load_heat
  use atl_eqn_heat_1d_var_module,     only: atl_init_heat_1d_vars
  use atl_eqn_heat_2d_var_module,     only: atl_init_heat_2d_vars
  use atl_eqn_heat_var_module,        only: atl_init_heat_vars
  use atl_varSys_module,              only: atl_varSys_solverData_type

  implicit none

  private

  public :: atl_eqn_heat_load_bc
  public :: atl_eqn_heat_init


contains


  !> Initialization of the Heat equation.
  !!
  !! This routine sets up the necessary infrastructure for the Heat
  !! equation.
  !! It reads the configuration from the given script in conf under the table
  !! provided in thandle and sets function pointers and variables accordingly.
  subroutine atl_eqn_heat_init( conf, thandle, equation, nDimensions, &
    &                           varSys_data                           )
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

    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), intent(inout) :: varSys_data
    !> The possible material properties.
    ! --------------------------------------------------------------------------

    equation%isNonlinear = .false.
    equation%nDimensions = nDimensions
    ! timesetp is static and not changing over simulationtime
    equation%adaptive_timestep = .false.
    ! we need 1st spatial derivatives of all variables ...
    equation%nDerivatives = 1

    equation%load_bc => atl_eqn_heat_load_bc

    select case(nDimensions)
    case(1)
      call atl_init_heat_1d_vars(             &
        & equation              = equation,   &
        & methodData            = varSys_data )
    case(2)
      call atl_init_heat_2d_vars(             &
        & equation              = equation,   &
        & methodData            = varSys_data )
    case(3)
      call atl_init_heat_vars(                &
        & equation              = equation,   &
        & methodData            = varSys_data )
    end select

    call atl_load_heat( heat         = equation%heat,   &
      &                 conf         = conf,            &
      &                 eq_table     = thandle          )

  end subroutine atl_eqn_heat_init


  !> Reading boundary conditions for the Heat equation.
  !!
  !! This routine has to conform to the interface definition
  !! atl_equation_module#eqn_load_bc.
  subroutine atl_eqn_heat_load_bc( equation,                              &
    &                              bc_state, bc_state_gradient,           &
    &                              bc_varDict, bc_varDict_gradient,       &
    &                              bc_normal_vec, bc_normal_vec_gradient, &
    &                              bc_trafo, bc_trafo_gradient,           &
    &                              bc_label, bc_kind, thandle, conf       )
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
    character(len=*), intent(in) :: bc_label
    character(len=*), intent(in) :: bc_kind
    type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo
    type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo_gradient
    integer, intent(in) :: thandle
    type(flu_State) :: conf
    ! ---------------------------------------------------------------------------
    type(tem_stringKeyValuePair_type) :: kvp


    ! By default we set the function pointer for a conversion,
    ! even if the boundary conditions does not use them.
    bc_trafo%from => null()
    bc_trafo%to => null()

    allocate(bc_state(1))
    allocate(bc_state_gradient(0))
!!VK    allocate(bc_trafo_gradient(0))
!!VK    allocate(bc_normal_vec_gradient(0))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )
    ! Constant zero variable for non-configurable boundary variable
    kvp%value = 'zero_const'

    select case(bc_kind)
    case('const', 'dirichlet')
      ! Prescribe temperature
      bc_state(1)%state_name = 'temp'
      bc_state(1)%style = 'dirichlet'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

    !>@todo case('heat_flux', 'neumann')


    !>@todo case('perfectly_insulated')

    case default
      write(logUnit(1),*) 'Unknown boundary kind "' // trim(bc_kind) // '"'
      write(logUnit(1),*) 'for boundary  "' // trim(bc_label) // '".'
      write(logUnit(1),*) 'Available boundary kinds for Heat equation:'
      write(logUnit(1),*) ' * temp'
      write(logUnit(1),*) ' * heat_flux'
      write(logUnit(1),*) ' * perfectly_insulated'
      write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
      call tem_abort()

    end select

    call truncate( me = bc_varDict )
    call truncate( me = bc_varDict_gradient )

  end subroutine atl_eqn_heat_load_bc

end module atl_eqn_heat_hlp_module
