! Copyright (c) 2013-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013, 2015-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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

module atl_eqn_maxwell_hlp_module
  ! Aotus modules
  use aotus_module,                     only: flu_State

  ! Treelm modules
  use env_module,                       only: rk
  use tem_bc_module,                    only: tem_bc_state_type
  use tem_aux_module,                   only: tem_abort
  use tem_logging_module,               only: logUnit
  use tem_stringKeyValuePair_module,    only: tem_stringKeyValuePair_type,   &
    &                                         init, truncate, append,        &
    &                                         grw_stringKeyValuePairArray_type

  ! Polynomial modules
  use ply_poly_project_module,          only: ply_poly_project_type, &
      &                                       ply_poly_project_n2m,  &
      &                                       ply_poly_project_m2n
  use ply_oversample_module,            only: ply_convert2oversample,    &
    &                                         ply_convertFromoversample, &
    &                                         ply_convert2oversample,    &
    &                                         ply_convertFromoversample

  ! Ateles modules
  use atl_equation_module,              only: atl_equations_type, &
    &                                         atl_eqn_var_trafo_type
  use atl_bc_state_module,              only: atl_load_bc_state
  use atl_eqn_maxwell_var_module,       only: atl_init_maxwell_vars,        &
    &                                         atl_init_maxwell_sourceTerms, &
    &                                         atl_init_maxwell_material
  use atl_eqn_maxwell_2d_var_module,    only: atl_init_maxwell_2d_vars, &
    &                                         atl_init_maxwell_2d_sourceTerms
  use atl_eqn_maxwelldivcorr_var_module,         &
    & only: atl_init_maxwellDivCorr_vars,        &
    &       atl_init_maxwellDivCorr_sourceTerms, &
    &       atl_init_maxwellDivCorr_material

  use atl_eqn_maxwell_derive_module,    only: atl_eqn_maxwell_cons2prim, &
    &                                         atl_eqn_maxwell_prim2cons
  use atl_eqn_maxwell_2d_derive_module, only: atl_eqn_maxwell_2d_cons2prim, &
    &                                         atl_eqn_maxwell_2d_prim2cons
  use atl_eqn_maxwelldivcorr_derive_module,   &
    & only: atl_eqn_maxwelldivcorr_cons2prim, &
    &       atl_eqn_maxwelldivcorr_prim2cons
  use atl_varSys_module,                only: atl_varSys_solverData_type
  use atl_source_types_module,          only: atl_init_source_type
  use atl_materialPrp_module,           only: atl_init_material_type, &
    &                                         atl_material_type

  implicit none

  private

  public :: atl_eqn_maxwell_load_bc
  public :: atl_eqn_maxwellDivCor_load_bc
  public :: atl_eqn_maxwell_init
  public :: atl_eqn_maxwell_implicit_pen


contains


  !> Initialization of the Maxwell equations.
  !!
  !! This routine sets up the necessary infrastructure for the Maxwell
  !! equations.
  !! It reads the configuration from the given script in conf under the table
  !! provided in thandle and sets function pointers and variables accordingly.
  subroutine atl_eqn_maxwell_init( equation, nDimensions, divCorr, initSource, &
    &                              initMaterial, varSys_data )
    ! ---------------------------------------------------------------------------
    !> Equation system to set with this routine.
    type(atl_equations_type), intent(inout) :: equation

    !> Number of spatial dimensions, the Maxwell equations should live on.
    !!
    !! Has to be 2 or 3.
    integer, intent(in) :: nDimensions

    !> Use divergence correction?
    logical, intent(in) :: divCorr

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

    equation%isNonlinear = .false.
    equation%nDimensions = nDimensions
    ! timesetp is static and not changing over simulationtime
    equation%adaptive_timestep = .false.

    if (divCorr) then
      equation%load_bc => atl_eqn_maxwellDivCor_load_bc
      call atl_init_maxwellDivCorr_vars( equation   = equation,   &
        &                                methodData = varSys_data )

      call atl_init_maxwellDivCorr_sourceTerms( initSource%poss_srcVars, &
        &                                       initSource%eval_source   )

      call atl_init_maxwellDivCorr_material( initMaterial%poss_materialVars )
    else
      select case(nDimensions)
      case(2)
        equation%load_bc => eqn_maxwell_2d_load_bc
        call atl_init_maxwell_2d_vars(          &
          & equation              = equation,   &
          & methodData            = varSys_data )

        call atl_init_maxwell_2d_sourceTerms( initSource%poss_srcVars, &
          &                                   initSource%eval_source   )
      case(3)
        equation%load_bc => atl_eqn_maxwell_load_bc

        call atl_init_maxwell_vars( equation    = equation,   &
          &                         methodData  = varSys_data )

        call atl_init_maxwell_sourceTerms( initSource%poss_srcVars, &
          &                                initSource%eval_source   )

      end select

      call atl_init_maxwell_material( initMaterial%poss_materialVars )

    end if

  end subroutine atl_eqn_maxwell_init



  !> Subroutine to load boundary conditions for Maxwell equations from
  !! a Lua configuration file.
  subroutine atl_eqn_maxwell_load_bc( equation,                              &
    &                                 bc_state, bc_state_gradient,           &
    &                                 bc_varDict, bc_varDict_gradient,       &
    &                                 bc_normal_vec, bc_normal_vec_gradient, &
    &                                 bc_trafo, bc_trafo_gradient,           &
    &                                 bc_label, bc_kind, thandle, conf       )
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
    type(tem_stringKeyValuePair_type) :: kvp
    ! --------------------------------------------------------------------------

    allocate(bc_state(equation%varSys%nScalars))
    allocate(bc_state_gradient(0))
!!VK    allocate(bc_normal_vec_gradient(2))
!!VK    allocate(bc_trafo_gradient(2))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )

    ! Constant zero variable for non-configurable boundary variable
    kvp%value = 'zero_const'

    ! By default we set the function pointer for a conversion,
    ! even if the boundary conditions does not use them.
    bc_trafo%from => atl_eqn_maxwell_prim2cons
    bc_trafo%to => atl_eqn_maxwell_cons2prim

    ! Check for Navier-Stokes specific boundary conditions.
    select case(bc_kind)
    case('pec')
      ! Use a non-slip boundary for walls in the Navier-Stokes equations.
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! Prescribe e_normal
      bc_state(1)%state_name = 'e_norm'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe e_tangential_1
      bc_state(2)%state_name = 'e_tan'
      bc_state(2)%style = 'dirichlet'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe v_tangential_2
      bc_state(3)%state_name = 'e_tan2'
      bc_state(3)%style = 'dirichlet'
      bc_state(3)%isDefined = .true.
      kvp%key = trim(bc_state(3)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe b_normal
      bc_state(4)%state_name = 'b_norm'
      bc_state(4)%style = 'dirichlet'
      bc_state(4)%isDefined = .true.
      kvp%key = trim(bc_state(4)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe b_tangential_1
      bc_state(5)%state_name = 'b_tan'
      bc_state(5)%style = 'neumann'
      bc_state(5)%isDefined = .true.
      kvp%key = trim(bc_state(5)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe b_tangential_2
      bc_state(6)%state_name = 'b_tan2'
      bc_state(6)%style = 'neumann'
      bc_state(6)%isDefined = .true.
      kvp%key = trim(bc_state(6)%state_name)
      call append( me = bc_varDict, val = kvp )

    case('conservatives')
      ! Everything is in conservative, so we do not need a transformation.
      bc_trafo%identity = .true.
      bc_normal_vec = .false.
      call atl_load_bc_state( bc          = bc_state(1),          &
        &                     state_name  = 'displacement_fieldX', &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(2),          &
        &                     state_name  = 'displacement_fieldY', &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(3),          &
        &                     state_name  = 'displacement_fieldZ', &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(4),          &
        &                     state_name  = 'magnetic_fieldX',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(5),          &
        &                     state_name  = 'magnetic_fieldY',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(6),          &
        &                     state_name  = 'magnetic_fieldZ',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition conservatives you have to set'
        write(logUnit(1),*) 'all conservative variables (displacement_fieldX,'
        write(logUnit(1),*) 'displacement_fieldY, displacement_fieldZ,'
        write(logUnit(1),*) 'magnetic_fieldX, magnetic_fieldY,'
        write(logUnit(1),*) 'magnetic_fieldZ) this set is not'
        write(logUnit(1),*) ' complete for ' // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case default
      call tem_abort( &
        & 'ERROR in atl_eqn_maxwell_load_bc: unknown boundary condition' )

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

  end subroutine atl_eqn_maxwell_load_bc


  !> Subroutine to load boundary conditions for 2D Maxwell equations from
  !! a Lua configuration file.
  subroutine eqn_maxwell_2d_load_bc( equation,                              &
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
    type(tem_stringKeyValuePair_type) :: kvp
    ! ---------------------------------------------------------------------------

    ! @todo PV pec is setting 7 bc_states and conservatives sets only
    ! 3 bc_states but bc_state is allocated with fixed size of nScalars
    allocate(bc_state(equation%varSys%nScalars))
    allocate(bc_state_gradient(0))
!!VK    allocate(bc_normal_vec_gradient(2))
!!VK    allocate(bc_trafo_gradient(2))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )
    ! Constant zero variable
    kvp%value = 'zero_const'

    ! By default we set the function pointer for a conversion,
    ! even if the boundary conditions does not use them.
    bc_trafo%from => atl_eqn_maxwell_2d_prim2cons
    bc_trafo%to => atl_eqn_maxwell_2d_cons2prim

    ! Check for Navier-Stokes specific boundary conditions.
    select case(bc_kind)
    case('pec')
      ! Use a non-slip boundary for walls in the Navier-Stokes equations.
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! Prescribe e_normal
      bc_state(1)%state_name = 'e_norm'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = bc_state(1)%state_name
      call append( me = bc_varDict, val = kvp )

      ! Prescribe e_tangential_1
      bc_state(2)%state_name = 'e_tan'
      bc_state(2)%style = 'dirichlet'
      bc_state(2)%isDefined = .true.
      kvp%key = bc_state(2)%state_name
      call append( me = bc_varDict, val = kvp )

      ! Prescribe b_z / b_tangential
      bc_state(3)%state_name = 'b_tan'
      bc_state(3)%style = 'neumann'
      bc_state(3)%isDefined = .true.
      kvp%key = bc_state(3)%state_name
      call append( me = bc_varDict, val = kvp )

      ! Prescribe zeros for all PML equations
      bc_state(4)%state_name = 'pmlP_norm'
      bc_state(5)%state_name = 'pmlP_tan'
      bc_state(6)%state_name = 'pmlQ_norm'
      bc_state(7)%state_name = 'pmlQ_tan'
      bc_state(4:7)%style = 'dirichlet'
      bc_state(4:7)%isDefined = .true.
      kvp%key = bc_state(4)%state_name
      call append( me = bc_varDict, val = kvp )
      kvp%key = bc_state(5)%state_name
      call append( me = bc_varDict, val = kvp )
      kvp%key = bc_state(6)%state_name
      call append( me = bc_varDict, val = kvp )
      kvp%key = bc_state(7)%state_name
      call append( me = bc_varDict, val = kvp )

    case('conservatives')
      ! Everything is in conservative, so we do not need a transformation.
      bc_trafo%identity = .true.
      bc_normal_vec = .false.

      call atl_load_bc_state( bc          = bc_state(1),          &
        &                     state_name  = 'displacement_fieldX', &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(2),          &
        &                     state_name  = 'displacement_fieldY', &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(3),          &
        &                     state_name  = 'magnetic_fieldZ',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition conservatives you have to set'
        write(logUnit(1),*) 'all conservative variables (displacement_fieldX,'
        write(logUnit(1),*) 'displacement_fieldY, magnetic_fieldZ) '
        write(logUnit(1),*) 'this set is not'
        write(logUnit(1),*) ' complete for ' // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case default
      call tem_abort(                                                   &
        & 'ERROR in eqn_maxwell_2d_load_bc: unknown boundary condition' )

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

  end subroutine eqn_maxwell_2d_load_bc

  !> Subroutine to load boundary conditions for Maxwell equations with divergence correction from
  !! a Lua configuration file. For the correction in E-field dirichlet bc and for correction in B-field
  !! neumann bc are defined.
  subroutine atl_eqn_maxwellDivCor_load_bc( equation,                     &
    & bc_state, bc_state_gradient, bc_varDict, bc_varDict_gradient,       &
    & bc_normal_vec, bc_normal_vec_gradient, bc_trafo, bc_trafo_gradient, &
    & bc_label, bc_kind, thandle, conf                                    )
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
    !> Unused variable, Needed for interface compatibility in module
    !! atl_equation_init_module.f90.
    character(len=*), intent(in) :: bc_label
    character(len=*), intent(in) :: bc_kind
    !> Unused variable, Needed for interface compatibility in module
    !! atl_equation_init_module.f90.
    integer, intent(in) :: thandle
    !> Unused variable, Needed for interface compatibility in module
    !! atl_equation_init_module.f90.
    type(flu_State) :: conf
    ! ---------------------------------------------------------------------------
    type(tem_stringKeyValuePair_type) :: kvp
    ! ---------------------------------------------------------------------------

    allocate(bc_state(equation%varSys%nScalars))
    allocate(bc_state_gradient(0))
!!VK    allocate(bc_normal_vec_gradient(2))
!!VK    allocate(bc_trafo_gradient(2))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )
    ! Constant zero variable
    kvp%value = 'zero_const'

    ! By default we set the function pointer for a conversion,
    ! even if the boundary conditions does not use them.
    bc_trafo%from => atl_eqn_maxwelldivcorr_prim2cons
    bc_trafo%to => atl_eqn_maxwelldivcorr_cons2prim

    ! Check for Navier-Stokes specific boundary conditions.
    select case(bc_kind)
    case('pec')
      ! Use a non-slip boundary for walls in the Navier-Stokes equations.
      bc_trafo%identity = .false.
      bc_normal_vec = .true.

      ! Prescribe e_normal
      bc_state(1)%state_name = 'e_norm'
      bc_state(1)%style = 'neumann'
      bc_state(1)%isDefined = .true.
      kvp%key = trim(bc_state(1)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe e_tangential_1
      bc_state(2)%state_name = 'e_tan'
      bc_state(2)%style = 'dirichlet'
      bc_state(2)%isDefined = .true.
      kvp%key = trim(bc_state(2)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe v_tangential_2
      bc_state(3)%state_name = 'e_tan'
      bc_state(3)%style = 'dirichlet'
      bc_state(3)%isDefined = .true.
      kvp%key = trim(bc_state(3)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe b_normal
      bc_state(4)%state_name = 'b_norm'
      bc_state(4)%style = 'dirichlet'
      bc_state(4)%isDefined = .true.
      kvp%key = trim(bc_state(4)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe b_tangential_1
      bc_state(5)%state_name = 'b_tan'
      bc_state(5)%style = 'neumann'
      bc_state(5)%isDefined = .true.
      kvp%key = trim(bc_state(5)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe b_tangential_2
      bc_state(6)%state_name = 'b_tan'
      bc_state(6)%style = 'neumann'
      bc_state(6)%isDefined = .true.
      kvp%key = trim(bc_state(6)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe e_normal
      bc_state(7)%state_name = 'e_cor'
      bc_state(7)%style = 'dirichlet'
      bc_state(7)%isDefined = .true.
      kvp%key = trim(bc_state(7)%state_name)
      call append( me = bc_varDict, val = kvp )

      ! Prescribe e_tangential_1
      bc_state(8)%state_name = 'b_cor'
      bc_state(8)%style = 'neumann'
      bc_state(8)%isDefined = .true.
      kvp%key = trim(bc_state(8)%state_name)
      call append( me = bc_varDict, val = kvp )

    case('conservatives')
      ! Everything is in conservative, so we do not need a transformation.
      bc_trafo%identity = .true.
      bc_normal_vec = .false.
      call atl_load_bc_state( bc          = bc_state(1),          &
        &                     state_name  = 'displacement_fieldX', &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(2),          &
        &                     state_name  = 'displacement_fieldY', &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(3),          &
        &                     state_name  = 'displacement_fieldZ', &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(4),          &
        &                     state_name  = 'magnetic_fieldX',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(5),          &
        &                     state_name  = 'magnetic_fieldY',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(6),          &
        &                     state_name  = 'magnetic_fieldZ',     &
        &                     conf        = conf,                 &
        &                     bc_handle   = thandle,              &
        &                     varDict     = bc_varDict,           &
        &                     varSys      = equation%varSys       )

      call atl_load_bc_state( bc          = bc_state(7),           &
        &                     state_name  = 'electric_correction', &
        &                     conf        = conf,                  &
        &                     bc_handle   = thandle,               &
        &                     varDict     = bc_varDict,            &
        &                     varSys      = equation%varSys        )

      call atl_load_bc_state( bc          = bc_state(8),           &
        &                     state_name  = 'magnetic_correction', &
        &                     conf        = conf,                  &
        &                     bc_handle   = thandle,               &
        &                     varDict     = bc_varDict,            &
        &                     varSys      = equation%varSys        )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition conservatives you have to set'
        write(logUnit(1),*) 'all conservative variables (displacement_fieldX,'
        write(logUnit(1),*) 'displacement_fieldY, displacement_fieldZ,'
        write(logUnit(1),*) 'magnetic_fieldX, magnetic_fieldY,'
        write(logUnit(1),*) 'magnetic_fieldZ, electric_correction,'
        write(logUnit(1),*) 'magnetic_correction) this set is not'
        write(logUnit(1),*) ' complete for ' // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

    case default
      call tem_abort(                                                          &
        & 'ERROR in atl_eqn_maxwellDivCor_load_bc: unknown boundary condition' )

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

  end subroutine atl_eqn_maxwellDivCor_load_bc


  ! ------------------------------------------------------------------------ !
  !> Solve the equation system with just the penalization terms to find an
  !! implicit update for the IMEX timestepping procedure.
  subroutine atl_eqn_maxwell_implicit_pen( material, weighted_dt, nDims,  &
    &                                      poly_proj, state, timestep_rhs )
    !> Definition of the material, which directly describes the penalization.
    !!
    !! We expect the conductivity sigma to be defined in materialdat(:,:,3)
    !! and the permittivity epsilon in materialdat(:,:,2).
    type(atl_material_type), intent(in) :: material

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
    integer :: nElems
    integer :: nPoints
    real(kind=rk) :: sigbyeps
    real(kind=rk) :: factor
    real(kind=rk), allocatable :: modalCoeff(:,:), modalCoeff_cur(:,:)
    real(kind=rk), allocatable :: pointVal(:,:), cur(:,:)
    ! -------------------------------------------------------------------- !


    nElems = material%material_desc%computeElems(1)%nElems
    constElems: do iMatElem=1,nElems
      iElem = material%material_desc%computeElems(1)%totElemIndices(iMatElem)
      sigbyeps = material%material_dat%elemMaterialData(1)       &
        &                             %materialDat(iMatElem,1,3) &
        &        / material%material_dat%elemMaterialData(1)     &
        &                               %materialDat(iMatElem,1,2)
      factor = 1 + weighted_dt*sigbyeps
      factor = 1.0_rk / factor

      ! u_i
      state(iElem,:,:nDims) = state(iElem,:,:nDims) * factor

      ! g(u_i)
      timestep_rhs(iElem,:,:nDims) &
        &   = (-1.0_rk) * sigbyeps * state(iElem,:,:nDims)
      timestep_rhs(iElem,:,nDims+1:) = 0.0_rk

    end do constElems

    nElems = material%material_desc%computeElems(2)%nElems
    if (nElems > 0) then
      if (nDims == 2) then
        nPoints = poly_proj%body_2d%nquadpoints
        allocate(modalCoeff(poly_proj%body_2d%oversamp_dofs,nDims))
      else
        nPoints = poly_proj%body_3d%nquadpoints
        allocate(modalCoeff(poly_proj%body_3d%oversamp_dofs,nDims))
      end if
      allocate(pointVal(nPoints,nDims))
      allocate(cur(nPoints,nDims))
      allocate(modalCoeff_cur(nPoints,nDims))
    end if
    varElems: do iMatElem=1,nElems
      iElem = material%material_desc%computeElems(2)%totElemIndices(iMatElem)

      call ply_convert2oversample( state       = state(iElem,:,:nDims), &
        &                          ndim        = nDims,                 &
        &                          poly_proj   = poly_proj,             &
        &                          modalCoeffs = modalCoeff(:,:nDims)   )

      call ply_poly_project_m2n(me         = poly_proj, &
        &                       dim        = nDims,     &
        &                       nVars      = nDims,     &
        &                       nodal_data = pointVal,  &
        &                       modal_data = modalCoeff )

      do iPoint=1,nPoints
        sigbyeps = material%material_dat%elemMaterialData(2)            &
          &                             %materialDat(iMatElem,iPoint,3) &
          &        / material%material_dat%elemMaterialData(2)          &
          &                               %materialDat(iMatElem,iPoint,2)
        factor = 1 + weighted_dt*sigbyeps
        factor = 1.0_rk / factor

        ! compute u_i
        pointVal(iPoint,:nDims) = pointVal(iPoint,:nDims) * factor

        ! compute g(u_i)
        cur(iPoint,:nDims) = (-1.0_rk) * sigbyeps * pointVal(iPoint,:nDims)
      end do

      call ply_poly_project_n2m( me         = poly_proj, &
        &                        dim        = nDims,     &
        &                        nVars      = nDims,     &
        &                        nodal_data = pointVal,  &
        &                        modal_data = modalCoeff )

      ! u_i
      call ply_convertFromOversample( modalCoeffs = modalCoeff(:,:nDims), &
        &                             ndim        = nDims,                &
        &                             poly_proj   = poly_proj,            &
        &                             state       = state(iElem,:,:nDims) )

      ! write g(u_i)
      call ply_poly_project_n2m( me         = poly_proj,     &
        &                        dim        = nDims,         &
        &                        nVars      = nDims,         &
        &                        nodal_data = cur,           &
        &                        modal_data = modalCoeff_cur )
      call ply_convertFromOversample(                   &
        &    modalCoeffs = modalCoeff_cur(:,:nDims),    &
        &    ndim        = nDims,                       &
        &    poly_proj   = poly_proj,                   &
        &    state       = timestep_rhs(iElem,:,:nDims) )
      timestep_rhs(iElem,:,nDims+1:) = 0.0_rk

    end do varElems
    if (nElems > 0) then
      deallocate(modalCoeff)
      deallocate(pointVal)
      deallocate(cur)
      deallocate(modalCoeff_cur)
    end if

  end subroutine atl_eqn_maxwell_implicit_pen
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module atl_eqn_maxwell_hlp_module
