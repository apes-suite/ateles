! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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

module atl_eqn_filnvrStk_hlp_module
  use env_module,                      only: labelLen
  use tem_tools_module,                only: upper_to_lower
  use aotus_module,                    only: flu_State, aot_get_val

  use tem_aux_module,                  only: tem_abort
  use tem_bc_module,                   only: tem_bc_state_type
  use tem_logging_module,              only: logUnit
  use tem_stringKeyValuePair_module,   only: tem_stringKeyValuePair_type
  use tem_stringKeyValuePair_module,   only: grw_stringKeyValuePairArray_type, &
    &                                        init, truncate, append

  use atl_equation_module,             only: atl_equations_type,     &
    &                                        atl_eqn_var_trafo_type
  use atl_bc_state_module,             only: atl_load_bc_state
  use atl_eqn_nvrstk_module,           only: atl_load_filtNS
  use atl_eqn_euler_derive_module,     only: atl_eqn_euler_cons2primTemp,   &
    &                                        atl_eqn_euler_primTemp2Cons,   &
    &                                        atl_eqn_euler_cons2primVel,    &
    &                                        atl_eqn_euler_primVel2Cons,    &
    &                                        atl_eqn_euler_cons2prim,       &
    &                                        atl_eqn_euler_prim2Cons
  use atl_eqn_euler_2d_derive_module,  only: atl_eqn_euler_2d_cons2primTemp,   &
    &                                        atl_eqn_euler_2d_primTemp2Cons,   &
    &                                        atl_eqn_euler_2d_cons2primVel,    &
    &                                        atl_eqn_euler_2d_primVel2Cons,    &
    &                                        atl_eqn_euler_2d_cons2prim,       &
    &                                        atl_eqn_euler_2d_prim2Cons
  use atl_eqn_filNvrStk_derive_module, only: atl_eqn_rans_cons2prim_elems,&
    &                                        atl_eqn_rans_prim2cons_elems,&
    &                                        atl_eqn_rans_2d_prim2cons_elems,&
    &                                        atl_eqn_rans_2d_cons2prim_elems
  use atl_eqn_euler_var_module,        only: atl_init_euler_material
  use atl_varSys_module,               only: atl_varSys_solverData_type
  use atl_source_types_module,         only: atl_init_source_type
  use atl_materialPrp_module,          only: atl_init_material_type
  use atl_eqn_filNvrStk_var_module,    only: atl_init_Rans_vars, &
    &                                        atl_init_rans_2d_vars,&
    &                                        atl_init_filNvrStk_sourceTerms, &
    &                                        atl_init_RANS_closure_coeffs
  use atl_laxFriedrichFlux_module,     only: atl_laxFriedRans,           &
    &                                        atl_laxFriedRans_2D


  implicit none

  private

  public :: atl_eqn_rans_load_bc
  public :: atl_eqn_filtered_nvrstk_init


contains


  !> Initialization of the Navier-Stokes equations.
  !!
  !! This routine sets up the necessary infrastructure for the Navier-Stokes
  !! equations.
  !! It reads the configuration from the given script in conf under the table
  !! provided in thandle and sets function pointers and variables accordingly.
  subroutine atl_eqn_filtered_nvrstk_init( conf, thandle, equation,  &
    &                                      nDimensions, initSource,  &
    &                                      initMaterial, varSys_data )
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
    ! timesetp is dynamic and changing over simulationtime
    equation%adaptive_timestep = .true.
    ! we need 1st spatial derivatives
    equation%nDerivatives = 1
    equation%nDimensions = nDimensions

    ! Load the equation specific parameters from the lua file
    call atl_load_filtNS(                            &
      & NavierStokes     = equation%NavierStokes,    &
      & FiltNavierStokes = equation%FiltNavierStokes,&
      & Euler            = equation%Euler,           &
      & conf             = conf,                     &
      & eq_table         = thandle                   )

    select case(nDimensions)
    case(3)
      select case(trim(equation%FiltNavierStokes%model_type))
      case('rans')
        equation%load_bc => atl_eqn_rans_load_bc


        ! Initialize the Variable system for the equation
        call atl_init_Rans_vars( equation   = equation,   &
          &                      solverData = varSys_data )


       ! call init_Rans_sourceTerms(                                    &
       !   &       possVars    = initSource%poss_srcVars,               &
       !   &       eval_source = initSource%eval_source,                &
       !   &       model_type  = equation%FiltNavierStokes%model_type,  &
       !   &       sourceDict  = initSource%sourceDict                  )

        ! Get the numerical flux !
        call atl_eqn_getFilNvrStkFluxes( equation   = equation,   &
          &                              conf       = conf,       &
          &                              eqn_handle = thandle,    &
          &                              eqn_dim    = nDimensions )

        ! Set up the links to the routine to convert variables
        ! whenever required in the code
        equation%cons2prim => atl_eqn_rans_cons2prim_elems
        equation%prim2cons => atl_eqn_rans_prim2cons_elems

      ! Deactivated due to default case
      !case('les')
      !  write(logUnit(1),*) 'Not yet implemented'
      !  call tem_abort()

      case default
        write(logUnit(1),*) 'This turbulence model is unknown!'
        write(logUnit(1),*) 'Please use turbulence_model="rans".'
        call tem_abort()

      end select

    case(2)

      select case(trim(equation%FiltNavierStokes%model_type))
      case('rans_2d')
        equation%load_bc => atl_eqn_rans_load_bc

        ! Initialize the Variable system for the equation
        call atl_init_rans_2d_vars(                  &
          &                 equation   = equation,   &
          &                 solverData = varSys_data )

        call atl_init_filNvrStk_sourceTerms(                           &
          &       possVars    = initSource%poss_srcVars,               &
          &       eval_source = initSource%eval_source,                &
          &       model_type  = equation%FiltNavierStokes%model_type   )

        ! Get the numerical flux !
        call atl_eqn_getFilNvrStkFluxes( equation   = equation,   &
          &                              conf       = conf,       &
          &                              eqn_handle = thandle,    &
          &                              eqn_dim    = nDimensions )

        ! Set up the links to the routine to convert variables
        ! whenever required in the code
        equation%cons2prim => atl_eqn_rans_2d_cons2prim_elems
        equation%prim2cons => atl_eqn_rans_2d_prim2cons_elems

        ! Fill up the RANS closure coefficients
        call atl_init_RANS_closure_coeffs(equation)

      ! Deactivated due to default case
      !case('les')
      !  write(logUnit(1),*) 'Not yet implemented'
      !  call tem_abort()

      case default
        write(logUnit(1),*) 'This turbulence model is unknown!'
        write(logUnit(1),*) 'Please use turbulence_model="rans_2d".'
        call tem_abort()

      end select

    case default
      write(logUnit(1),*) "Not Implemented yet"
      call tem_abort()

    end select

    call atl_init_euler_material(                     &
      & possVars    = initMaterial%poss_materialVars, &
      & nDimensions = nDimensions                     )

  end subroutine atl_eqn_filtered_nvrstk_init



  ! *****************************************************************************
  ! Getting the numerical flux for Euler equations
  subroutine atl_eqn_getFilNvrStkFluxes(equation, conf, eqn_handle, eqn_dim)
    !> The equations type to set the numerical flux in.
    type(atl_equations_type), intent(inout) :: equation
    !> Configuration file handle to get the numerical flux setting from.
    type(flu_state) :: conf
    !> Handle to the equation table in the configuration script.
    integer, intent(in) :: eqn_handle
    !> Dimension of the equation to set the flux for.
    integer, intent(in) :: eqn_dim

    character(len=labelLen) :: eq_nflux
    integer :: iError

    call aot_get_val(L = conf, thandle = eqn_handle, key = 'numflux', &
      &              val =  eq_nflux, ErrCode = iError,             &
      &              default = 'lax_friedrich')
    eq_nflux = upper_to_lower(eq_nflux)
    eq_nflux = adjustl(eq_nflux)

    select case(trim(equation%FiltNavierStokes%model_type))
    case('rans', 'rans_2d')
      select case(trim(eq_nflux))
        case ('hll')
         !NA! @todo : Not implemented
         !NA! write(logunit(2),*) 'Using HLL numerical flux.'
         !NA! write(logunit(2),*) 'Warning, this flux ignores materials completely!'
         !NA! select case(eqn_dim)
         !NA! case(1)
         !NA!   equation%numflux => atl_HLLEuler1D
         !NA! case(2)
         !NA!   equation%numflux => atl_HLLEuler2D
         !NA! case(3)
         !NA!   equation%numflux => atl_HLLEuler
         !NA! end select

        case ('godunov')
         !NA! @todo : Not implemented
         !NA! write(logunit(2),*) 'Using Godunov numerical flux.'
         !NA! write(logunit(2),*) 'Warning, this flux does not handle materials'
         !NA! write(logunit(2),*) 'completely correct!'
         !NA! select case(eqn_dim)
         !NA! case(1)
         !NA!   equation%numflux => atl_GodunovEuler1D
         !NA! case(2)
         !NA!   equation%numflux => atl_GodunovEuler2D
         !NA! case(3)
         !NA!   equation%numflux => atl_GodunovEuler
         !NA! end select

        case ('lax_friedrich')
          write(logunit(2),*) 'Using Lax Friedrichs numerical flux.'
          select case(eqn_dim)
          case(2)
            equation%euler%numflux => atl_laxFriedRans_2D
          case(3)
            equation%euler%numflux => atl_laxFriedRans
          end select
        case default

          write(logunit(1),*) 'Unknown numerical flux ', trim(eq_nflux)
          write(logunit(1),*) 'for the Euler equation system.'
          write(logunit(1),*) 'Please choose one of the available:'
          write(logunit(1),*) ' * lax_friedrich (default)'
          write(logunit(1),*) ' * godunov'
          write(logunit(1),*) ' * hll'
          write(logunit(1),*) ''
          write(logunit(1),*) 'Aborting!'
          call tem_abort()

      end select
    case('les')
      write(logUnit(1),*) 'Not yet implemented'
      call tem_abort()


    end select


  end subroutine atl_eqn_getFilNvrStkFluxes



  subroutine atl_eqn_rans_load_bc( equation,                              &
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

    ! Navier-Stokes requires 2 + nDims variables to be set:
    ! (pressure is stored in the last position, thus its index is given by:)
    pIndex = 2+nDims

    allocate(bc_state(pIndex))
    allocate(bc_state_gradient(pIndex))
!!VK    allocate(bc_normal_vec_gradient(2))
!!VK    allocate(bc_trafo_gradient(2))

    ! Initialize varDict for current boundary
    call init( me = bc_varDict )
    call init( me = bc_varDict_gradient )
    ! Constant zero variable for non-configurable boundary variable
    kvp%value = 'zero_const'

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
!!VK        call atl_load_bc_state( bc          = bc_state_gradient(3,1), &
!!VK          &                     state_name  = 'v_tan',                &
!!VK          &                     style       = 'dirichlet_mirror',     &
!!VK          &                     conf        = conf,                   &
!!VK          &                     bc_handle   = thandle,                &
!!VK          &                     varsys      = equation%varsys,        &
!!VK          &                     varDict     = bc_varDict_gradient     )

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
!!VK            &                     varsys      = equation%varsys,        &
!!VK            &                     varDict     = bc_varDict_gradient     )

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
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict        )

!!VK      ! Temperature for viscous terms
!!VK      call atl_load_bc_state( bc          = bc_state_gradient(pIndex,1), &
!!VK        &                     state_name  = 'temperature',               &
!!VK        &                     style       = 'dirichlet_mirror',          &
!!VK        &                     conf        = conf,                        &
!!VK        &                     bc_handle   = thandle,                     &
!!VK        &                     varsys      = equation%varsys,             &
!!VK        &                     varDict     = bc_varDict_gradient          )

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

      ! For the viscous terms ...
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

  !!VK    ! v_normal for viscous terms
  !!VK    bc_state_gradient(2,1)%state_name = 'v_norm'
  !!VK    bc_state_gradient(2,1)%style = 'dirichlet_mirror'
  !!VK    bc_state_gradient(2,1)%isDefined = .true.
  !!VK    kvp%key = trim(bc_state_gradient(2,1)%state_name)
  !!VK    call append( me = bc_varDict_gradient, val = kvp )

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
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons
      end select
      bc_trafo%identity = .false.
      bc_trafo_gradient%identity = .true.

      bc_normal_vec_gradient = .false.
      bc_normal_vec = .false.

      ! Impose density
      call atl_load_bc_state( bc          = bc_state(1), &
        &                     state_name  = 'density',   &
        &                     conf        = conf,        &
        &                     bc_handle   = thandle,     &
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict   )

      ! Impose velocity x
      call atl_load_bc_state( bc          = bc_state(2), &
        &                     state_name  = 'velocityX', &
        &                     conf        = conf,        &
        &                     bc_handle   = thandle,     &
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict   )
      if (nDims > 1) then
        ! Impose velocity y
        call atl_load_bc_state( bc          = bc_state(3), &
          &                     state_name  = 'velocityY', &
          &                     conf        = conf,        &
          &                     bc_handle   = thandle,     &
          &                     varsys      = equation%varsys,  &
          &                     varDict     = bc_varDict   )
        bc_state_gradient(3)%state_name = 'momentumY'
        if (nDims > 2) then
          ! Impose velocity z
          call atl_load_bc_state( bc          = bc_state(4), &
            &                     state_name  = 'velocityZ', &
            &                     conf        = conf,        &
            &                     bc_handle   = thandle,     &
            &                     varsys      = equation%varsys,  &
            &                     varDict     = bc_varDict   )
          bc_state_gradient(4)%state_name = 'momentumZ'
        end if
      end if

      ! Impose pressure
      call atl_load_bc_state( bc          = bc_state(pIndex), &
        &                     state_name  = 'pressure',       &
        &                     conf        = conf,             &
        &                     bc_handle   = thandle,          &
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict        )

      if (.not. all(bc_state(:)%isDefined)) then
        write(logUnit(1),*) 'For boundary condition primitives you have to set'
        write(logUnit(1),*) 'all primitive variables (density, velocityX, '
        write(logUnit(1),*) 'velocityY, velocityZ'
        write(logUnit(1),*) 'and pressure) this set is not complete for ' &
          &            // trim(bc_label) // '!'
        write(logUnit(1),*) 'Do not know how to proceed, ABORTING...'
        call tem_abort()
      end if

      bc_state_gradient(1)%state_name = 'density'
      bc_state_gradient(2)%state_name = 'momentumX'
      bc_state_gradient(pIndex)%state_name = 'energy'
      bc_state_gradient(:)%style = 'neumann'
      bc_state_gradient(:)%isDefined = .true.
      do iVar=1,pIndex
!!VK        kvp%key = trim(bc_state_gradient(iVar,1)%state_name)
!!VK        call append( me = bc_varDict_gradient, val = kvp )
        kvp%key = trim(bc_state_gradient(iVar)%state_name)
        call append( me = bc_varDict_gradient, val = kvp )
      end do

    case('inflow_normal')

      ! This boundary is given in primite variables, so we have
      ! to use a conversion.
      select case(nDims)
      case(2)
        bc_trafo%to => atl_eqn_euler_2d_cons2prim
        bc_trafo%from => atl_eqn_euler_2d_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons
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
      call atl_load_bc_state( bc          = bc_state(1), &
        &                     state_name  = 'density',   &
        &                     style       = 'dirichlet', &
        &                     conf        = conf,        &
        &                     bc_handle   = thandle,     &
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict   )

      ! ... it's gradient (density - extrapolate)
      bc_state_gradient(1)%state_name = 'density'
      bc_state_gradient(1)%style = 'neumann'
      bc_state_gradient(1)%isDefined = .true.

      ! Impose normal velocity
      call atl_load_bc_state( bc          = bc_state(2), &
        &                     state_name  = 'v_norm',    &
        &                     style       = 'dirichlet', &
        &                     conf        = conf,        &
        &                     bc_handle   = thandle,     &
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict   )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.

      if (nDims > 1) then
        ! Impose tangential velocity
        call atl_load_bc_state( bc          = bc_state(3), &
          &                     state_name  = 'v_tan',     &
          &                     style       = 'dirichlet', &
          &                     conf        = conf,        &
          &                     bc_handle   = thandle,     &
          &                     varsys      = equation%varsys,  &
          &                     varDict     = bc_varDict   )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.

        if (nDims > 2) then
          ! Impose tangential velocity
          call atl_load_bc_state( bc          = bc_state(4), &
            &                     state_name  = 'v_tan2',    &
            &                     style       = 'dirichlet', &
            &                     conf        = conf,        &
            &                     bc_handle   = thandle,     &
            &                     varsys      = equation%varsys,  &
            &                     varDict     = bc_varDict   )

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
!!VK        bc_state_gradient(iVar, 1) = bc_state_gradient(iVar,2)
!!VK        kvp%key = trim(bc_state_gradient(iVar,1)%state_name)
!!VK        call append( me = bc_varDict_gradient, val = kvp )
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
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons
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
      call atl_load_bc_state( bc          = bc_state(1), &
        &                     state_name  = 'density',   &
        &                     style       = 'dirichlet', &
        &                     conf        = conf,        &
        &                     bc_handle   = thandle,     &
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict   )

      ! ... it's gradient (density - extrapolate)
      bc_state_gradient(1)%state_name = 'density'
      bc_state_gradient(1)%style = 'neumann'
      bc_state_gradient(1)%isDefined = .true.

      ! Impose normal velocity
      call atl_load_bc_state( bc          = bc_state(2), &
        &                     state_name  = 'v_norm',    &
        &                     style       = 'dirichlet', &
        &                     conf        = conf,        &
        &                     bc_handle   = thandle,     &
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict   )

      ! ... it's gradient (momentum - extrapolate)
      bc_state_gradient(2)%state_name = 'momentum_norm'
      bc_state_gradient(2)%style = 'neumann'
      bc_state_gradient(2)%isDefined = .true.

      if (nDims > 1) then
        ! Impose tangential velocity
        call atl_load_bc_state( bc          = bc_state(3), &
          &                     state_name  = 'v_tan',     &
          &                     style       = 'dirichlet', &
          &                     conf        = conf,        &
          &                     bc_handle   = thandle,     &
          &                     varsys      = equation%varsys,  &
          &                     varDict     = bc_varDict   )

        ! ... it's gradient (momentum - extrapolate)
        bc_state_gradient(3)%state_name = 'momentum_tan'
        bc_state_gradient(3)%style = 'neumann'
        bc_state_gradient(3)%isDefined = .true.

        if (nDims > 2) then
          ! Impose tangential velocity
          call atl_load_bc_state( bc          = bc_state(4), &
            &                     state_name  = 'v_tan2',    &
            &                     style       = 'dirichlet', &
            &                     conf        = conf,        &
            &                     bc_handle   = thandle,     &
            &                     varsys      = equation%varsys,  &
            &                     varDict     = bc_varDict   )

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
        &                     varsys      = equation%varsys,  &
        &                     varDict     = bc_varDict        )

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
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons
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
        &                     varSys      = equation%varSys,  &
        &                     varDict     = bc_varDict        )

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
        bc_trafo_gradient%to => atl_eqn_euler_2d_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_2d_prim2cons
      case(3)
        bc_trafo%to => atl_eqn_euler_cons2prim
        bc_trafo%from => atl_eqn_euler_prim2cons
        bc_trafo_gradient%to => atl_eqn_euler_cons2prim
        bc_trafo_gradient%from => atl_eqn_euler_prim2cons
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
  end subroutine atl_eqn_rans_load_bc

end module atl_eqn_filNvrStk_hlp_module
