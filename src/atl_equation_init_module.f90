! Copyright (c) 2013-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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

!> Some helping routines for the equation system.
!!
!! This module mainly implements the initialization of the equation system.
!! It uses the various equation system implementations and calls their
!! initializations accordingly.
!! This design is intented to provide a fairly easy option to introduce new
!! equation systems, with their own implementations.
module atl_equation_init_module
  use env_module,                        only: labelLen

  use aotus_module,                      only: flu_State, &
    &                                          aot_get_val
  use aot_table_module,                  only: aot_table_open, &
    &                                          aot_table_close
  use aot_out_module,                    only: aot_out_type,       &
    &                                          aot_out_open_table, &
    &                                          aot_out_close_table

  use tem_aux_module,                    only: tem_abort
  use tem_tools_module,                  only: upper_to_lower
  use tem_coordinate_module,             only: initCoordinateRotation, &
    &                                          xToXAxes,               &
    &                                          yToXAxes,               &
    &                                          zToXAxes
  use tem_logging_module,                only: logUnit

  use atl_equation_module,               only: atl_equations_type

  use atl_eqn_euler_hlp_module,          only: atl_eqn_euler_init
  !!use atl_eqn_LoclinEuler_hlp_module     only: eqn_LoclinEuler_init
  use atl_eqn_advection_1d_hlp_module,   only: atl_eqn_advection_1d_init
  use atl_eqn_heat_hlp_module,           only: atl_eqn_heat_init
  use atl_eqn_LinearEuler_hlp_module,    only: atl_eqn_LinearEuler_init
  use atl_eqn_nvrstk_hlp_module,         only: atl_eqn_nvrstk_init
  use atl_eqn_filNvrstk_hlp_module,      only: atl_eqn_filtered_nvrstk_init
  use atl_eqn_maxwell_hlp_module,        only: atl_eqn_maxwell_init
  use atl_eqn_bbm_module,                only: atl_load_BBMEM
  use atl_eqn_nerplanck_module,          only: atl_load_nernstPlanck
  use atl_eqn_nerplanck_var_module,      only: atl_init_nerplanck_vars
  use atl_eqn_acoustic_hlp_module,       only: atl_eqn_acoustic_init
  use atl_varSys_module,                 only: atl_varSys_solverData_type
  use atl_source_types_module,           only: atl_init_source_type
  use atl_materialPrp_module,            only: atl_init_material_type
  use atl_materialFun_module,            only: atl_dump_materialFun

  implicit none

  private

  public :: atl_init_equation
  public :: atl_eqn_write


contains


  !*****************************************************************************
  !> This subroutine reads the equation system to solve from the configuration.
  !!
  !! It creates the necessary data structures and sets the required values
  !! according to the configuration.
  subroutine atl_init_equation( equation, conf, varSys_data, initSource, &
    &                           initMaterial                             )
    !--------------------------------------------------------------------------!
    !> Equation system to set with this routine.
    type(atl_Equations_type), intent(out) :: equation

    !> Handle to the Lua configuration
    type(flu_State), intent(in) :: conf

    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), intent(inout) :: varSys_data

    !> Type to be filled with the possible source variables for the equation
    !! system. These source variables are later on used to extract the
    !! corresponding information from the configuration file.
    type(atl_init_source_type), intent(inout) :: initSource

    !> Type to be filled with the possible material variables for the equation
    !! system. These material variables are later on used to extract the
    !! corresponding information from the configuration file.
    type(atl_init_material_type), intent(inout) :: initMaterial
    !--------------------------------------------------------------------------!
    integer                 :: eq_table, iError
    character(len=labelLen) :: eq_name
    !--------------------------------------------------------------------------!

    write(logUnit(1),'(A)') 'Setting up equation parameters ...'

    equation%requiresDeviation = .false.

    ! Try to open the equation table in the Lua script
    call aot_table_open(L=conf, thandle=eq_table, key='equation')
    if (eq_table == 0) then

      ! If there is no equation table set the eq_kind to 'unknown' for
      ! appropiate later processing.
      equation%eq_kind = 'unknown'

    else

      ! Select the equation system based upon the name specified in the
      ! equation table.
      ! All names are turned to lower case, to make the input case
      ! insensitive!
      call aot_get_val( L       = conf,     &
        &               thandle = eq_table, &
        &               key     = 'name',   &
        &               val     = eq_name,  &
        &               ErrCode = iError    )
      equation%eq_kind = upper_to_lower(eq_name)
      equation%eq_kind = adjustl(equation%eq_kind)
      equation%eq_kind = trim(equation%eq_kind)

    end if

    write(logUnit(4),*) 'Looking for equation ', trim(equation%eq_kind)

    select case(trim(equation%eq_kind))

    case('bbmem')
      ! Call a routine to intialize the BBMEM Equation type
      equation%isNonlinear = .false.
      equation%nDimensions = 3
      write(logUnit(1),*) 'Solving the BBMEM equation'
      call atl_load_BBMEM(equation%BBMEM, conf, eq_table)

    case('euler')
      ! call a routine to init euler with respect to
      ! the configuration file
      write(logUnit(1),*) 'Solving the Euler equations (3D)'
      call atl_eqn_euler_init( conf         = conf,         &
        &                      thandle      = eq_table,     &
        &                      equation     = equation,     &
        &                      nDimensions  = 3,            &
        &                      initSource   = initSource,   &
        &                      initMaterial = initMaterial, &
        &                      varSys_data  = varSys_data   )

    case('euler_2d')
      ! call a routine to init euler with respect to
      ! the configuration file
      write(logUnit(1),*) 'Solving the Euler equations (2D)'
      call atl_eqn_euler_init( conf         = conf,         &
        &                      thandle      = eq_table,     &
        &                      equation     = equation,     &
        &                      nDimensions  = 2,            &
        &                      initSource   = initSource,   &
        &                      initMaterial = initMaterial, &
        &                      varSys_data  = varSys_data   )

    case('euler_1d')
      ! call a routine to init euler with respect to
      ! the configuration file
      write(logUnit(1),*) 'Solving the Euler equations (1D)'
      call atl_eqn_euler_init( conf         = conf,         &
        &                      thandle      = eq_table,     &
        &                      equation     = equation,     &
        &                      nDimensions  = 1,            &
        &                      initSource   = initSource,   &
        &                      initMaterial = initMaterial, &
        &                      varSys_data  = varSys_data   )

    case('navier_stokes')
      ! call a routine to init navier stoke with respect to
      ! the configuration file
      write(logUnit(1),*) 'Solving the Navier-Stokes equations (3D)'
      call atl_eqn_nvrstk_init( conf         = conf,         &
        &                       thandle      = eq_table,     &
        &                       equation     = equation,     &
        &                       nDimensions  = 3,            &
        &                       initSource   = initSource,   &
        &                       initMaterial = initMaterial, &
        &                       varSys_data  = varSys_data   )

    case('filtered_navier_stokes')
      ! call a routine to init navier stoke with respect to
      ! the configuration file
      write(logUnit(1),*) 'Solving the filtered Navier-Stokes equations (3D)'
      call atl_eqn_filtered_nvrstk_init( conf         = conf,         &
        &                                thandle      = eq_table,     &
        &                                equation     = equation,     &
        &                                nDimensions  = 3,            &
        &                                initSource   = initSource,   &
        &                                initMaterial = initMaterial, &
        &                                varSys_data  = varSys_data   )

    case('navier_stokes_2d')
      ! call a routine to init navier stoke with respect to
      ! the configuration file
      write(logUnit(1),*) 'Solving the Navier-Stokes equations (2D)'
      call atl_eqn_nvrstk_init( conf         = conf,         &
        &                       thandle      = eq_table,     &
        &                       equation     = equation,     &
        &                       nDimensions  = 2,            &
        &                       initSource   = initSource,   &
        &                       initMaterial = initMaterial, &
        &                       varSys_data  = varSys_data   )

    case('filtered_navier_stokes_2d')
      write(logUnit(1),*) 'Solving the filtered Navier-Stokes equations (2D)'
      call atl_eqn_filtered_nvrstk_init( conf         = conf,         &
        &                                thandle      = eq_table,     &
        &                                equation     = equation,     &
        &                                nDimensions  = 2,            &
        &                                initSource   = initSource,   &
        &                                initMaterial = initMaterial, &
        &                                varSys_data  = varSys_data   )

    case('acoustic_2d')
      write(logUnit(1),*) 'Solving the Acoustic Wave equation (2D)'
      call atl_eqn_acoustic_init( conf         = conf,         &
        &                         thandle      = eq_table,     &
        &                         equation     = equation,     &
        &                         nDimensions  = 2,            &
        &                         initSource   = initSource,   &
        &                         varSys_data  = varSys_data   )

    case('acoustic')
      write(logUnit(1),*) 'Solving the Acoustic Wave equation (3D)'
      call atl_eqn_acoustic_init( conf         = conf,         &
        &                         thandle      = eq_table,     &
        &                         equation     = equation,     &
        &                         nDimensions  = 3,            &
        &                         initSource   = initSource,   &
        &                         varSys_data  = varSys_data   )

    case('lineareuler')
      write(logUnit(1),*) 'Solving the Linearized Euler equations (3D)'
      call atl_eqn_linearEuler_init( conf         = conf,         &
        &                            thandle      = eq_table,     &
        &                            equation     = equation,     &
        &                            nDimensions  = 3,            &
        &                            initSource   = initSource ,  &
        &                            varSys_data  = varSys_data   )
    case('lineareuler_2d')
      write(logUnit(1),*) 'Solving the Linearized Euler equations (2D)'
      call atl_eqn_linearEuler_init( conf         = conf,         &
        &                            thandle      = eq_table,     &
        &                            equation     = equation,     &
        &                            nDimensions  = 2,            &
        &                            initSource   = initSource ,  &
        &                            varSys_data  = varSys_data   )

    case('loclineuler')
      write(logunit(1),*) 'Solving the locally linearized Euler equations (3D)'
      call atl_eqn_euler_init( conf         = conf,          &
        &                      thandle      = eq_table,      &
        &                      equation     = equation,      &
        &                      nDimensions  = 3,             &
        &                      initSource   = initSource,    &
        &                      initMaterial = initMaterial,  &
        &                      varSys_data  = varSys_data    )

      equation%isNonlinear = .true.

    case('loclineuler_1d')
      write(logunit(1),*) 'Solving the locally linearized Euler equations (1D)'
      call atl_eqn_euler_init( conf         = conf,          &
        &                      thandle      = eq_table,      &
        &                      equation     = equation,      &
        &                      nDimensions  = 1,             &
        &                      initSource   = initSource,    &
        &                      initMaterial = initMaterial,  &
        &                      varSys_data  = varSys_data    )

      equation%isNonlinear = .true.

    case('maxwell', 'pec_maxwell', 'pec_maxwell_scatter')
      write(logUnit(1),*) 'Solving the Maxwell equations (3D)'
      call atl_eqn_maxwell_init( equation     = equation,     &
        &                        nDimensions  = 3,            &
        &                        divCorr      = .false.,      &
        &                        initSource   = initSource,   &
        &                        initMaterial = initMaterial, &
        &                        varSys_data  = varSys_data   )
      write(logUnit(1),*) 'Solving the Maxwell equations with PML (3D)'

    case('maxwell_2d', 'pec_maxwell_2d')
      write(logUnit(1),*) 'Solving the Maxwell equations (2D)'
      call atl_eqn_maxwell_init( equation     = equation,     &
        &                        nDimensions  = 2,            &
        &                        divCorr      = .false.,      &
        &                        initSource   = initSource,   &
        &                        initMaterial = initMaterial, &
        &                        varSys_data  = varSys_data   )

    case('maxwelldivcorrection')
      write(logUnit(1),*) 'Solving the Maxwell equations with DivCorr (3D)'
      call atl_eqn_maxwell_init( equation     = equation,     &
        &                        nDimensions  = 3,            &
        &                        divCorr      = .true.,       &
        &                        initSource   = initSource,   &
        &                        initMaterial = initMaterial, &
        &                        varSys_data  = varSys_data   )


    case('nernstplanck')
      write(logUnit(1),*) 'Solving the Nernst-Planck equation'
      equation%isNonlinear = .false.
      equation%nDimensions = 3

      call atl_load_nernstPlanck(equation%nerplanck, conf, eq_table)

      call atl_init_nerplanck_vars(                               &
        & varSys                = equation%varSys,                &
        & hasPrimitiveVariables = equation%hasPrimitiveVariables, &
        & methodData            = varSys_data                     )

    case('advection_1d')
      write(logUnit(1),*) 'Solving the 1D Advection equation'
      call atl_eqn_advection_1d_init( conf         = conf,         &
        &                             thandle      = eq_table,     &
        &                             equation     = equation,     &
        &                             nDimensions  = 1,            &
        &                             varSys_data  = varSys_data   )

    case('heat_1d')
      write(logUnit(1),*) 'Solving the 1D Heat equation'
      call atl_eqn_heat_init( conf         = conf,         &
        &                     thandle      = eq_table,     &
        &                     equation     = equation,     &
        &                     nDimensions  = 1,            &
        &                     varSys_data  = varSys_data   )


    case('heat_2d')
      write(logUnit(1),*) 'Solving the 2D Heat equation'
      call atl_eqn_heat_init( conf         = conf,         &
        &                     thandle      = eq_table,     &
        &                     equation     = equation,     &
        &                     nDimensions  = 2,            &
        &                     varSys_data  = varSys_data   )

    case('heat')
      write(logUnit(1),*) 'Solving the 3D Heat equation'
      call atl_eqn_heat_init( conf         = conf,         &
        &                     thandle      = eq_table,     &
        &                     equation     = equation,     &
        &                     nDimensions  = 3,            &
        &                     varSys_data  = varSys_data   )

    case ('unknown')
      ! No equation table provided, nothing to do...
      equation%isNonlinear = .false.
      equation%nDimensions = 3

    !> @todo All equation specific initialization routines should have a default
    !        case for not supported values for nDimension
    case default
       call tem_abort( 'ERROR in atl_init_equation: the equation type ' &
        & // trim(equation%eq_kind)                                     &
        & // ' from the lua configuration file is not known!'           )

    end select

    ! Close the Lua table with the equation description again.
    call aot_table_close(L = conf, thandle = eq_table)

    ! Build the permutations for the flux calculations across faces in
    ! x directions
    call atl_initCoordinateRotations( equation )

  end subroutine atl_init_equation
  !*****************************************************************************


  !*****************************************************************************
  subroutine atl_initCoordinateRotations( equation )
    !---------------------------------------------------------------------------
    !> Equation system to set with this routine.
    type(atl_Equations_type), intent(inout) :: equation
    !---------------------------------------------------------------------------

    select case(trim(equation%eq_kind))
    case( 'bbmem', 'euler', 'navier_stokes', 'lineareuler',          &
      & 'filtered_navier_stokes', 'maxwell', 'maxwelldivcorrection', &
      & 'nernstplanck', 'acoustic', 'heat','loclineuler'            )
      call initCoordinateRotation( varSys      = equation%varSys,         &
        &                          coordTrans  = xToXAxes,                &
        &                          derivatives = equation%nDerivatives,   &
        &                          rotation    = equation%varRotation(1), &
        &                          dimen       = 3                        )

      call initCoordinateRotation( varSys      = equation%varSys,         &
        &                          coordTrans  = yToXAxes,                &
        &                          derivatives = equation%nDerivatives,   &
        &                          rotation    = equation%varRotation(2), &
        &                          dimen       = 3                        )

      call initCoordinateRotation( varSys      = equation%varSys,         &
        &                          coordTrans  = zToXAxes,                &
        &                          derivatives = equation%nDerivatives,   &
        &                          rotation    = equation%varRotation(3), &
        &                          dimen       = 3                        )
    case( 'euler_2d', 'navier_stokes_2d', 'maxwell_2d', 'acoustic_2d',    &
      & 'heat_2d', 'filtered_navier_stokes_2d', 'lineareuler_2d'          )
      call initCoordinateRotation( varSys      = equation%varSys,         &
        &                          coordTrans  = xToXAxes,                &
        &                          derivatives = equation%nDerivatives,   &
        &                          rotation    = equation%varRotation(1), &
        &                          dimen       = 2                        )

      call initCoordinateRotation( varSys      = equation%varSys,         &
        &                          coordTrans  = yToXAxes,                &
        &                          derivatives = equation%nDerivatives,   &
        &                          rotation    = equation%varRotation(2), &
        &                          dimen       = 2                        )
    case( 'euler_1d', 'advection_1d', 'heat_1d','loclineuler_1d' )
      call initCoordinateRotation( varSys      = equation%varSys,         &
        &                          coordTrans  = xToXAxes,                &
        &                          derivatives = equation%nDerivatives,   &
        &                          rotation    = equation%varRotation(1), &
        &                          dimen       = 1                        )
    case default
      call tem_abort(                                                  &
        & 'ERROR in atl_initCoordinateRotations: the equation system ' &
        & // trim(equation%eq_kind)                                    &
        & // ' is not supported for coordinate transformations!'       )
    end select

    write(logUnit(1),*) ''

  end subroutine atl_initCoordinateRotations
  !*****************************************************************************


  !*****************************************************************************
  subroutine atl_eqn_write(equation, pconf)
    use atl_eqn_euler_module,        only: atl_save_euler
    use atl_eqn_nvrstk_module,       only: atl_save_navierStokes
    use atl_eqn_maxwell_module,      only: atl_save_maxwell
    use atl_eqn_advection_1d_module, only: atl_save_advection_1d
    use atl_eqn_acoustic_module,     only: atl_save_acoustic
    use atl_eqn_LinearEuler_module,  only: atl_save_linearEuler
    use atl_eqn_heat_module,         only: atl_save_heat
    !---------------------------------------------------------------------------
    !> Equation information to write.
    type(atl_Equations_type), intent(in) :: equation

    !> Lua output handle to write the equation information to.
    !!
    !!  Needs to have been opened with aot_out_open.
    type(aot_out_type) :: pconf
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    call aot_out_open_table(put_conf = pconf, tname = 'equation')

    select case(trim(equation%eq_kind))
    case('euler', 'euler_2d', 'euler_1d')
      call atl_save_euler( me       = equation%euler,         &
        &                  eqn_name = trim(equation%eq_kind), &
        &                  conf     = pconf                   )

    case( 'navier_stokes','navier_stokes_2d',                  &
        & 'filtered_navier_stokes', 'filtered_navier_stokes_2d')
      call atl_save_navierStokes( me_euler  = equation%euler,         &
        &                         me_nvrstk = equation%NavierStokes,  &
        &                         eqn_name  = trim(equation%eq_kind), &
        &                         conf      = pconf                   )

    case('advection_1d')
      call atl_save_advection_1d( me       = equation%advection,     &
        &                         eqn_name = trim(equation%eq_kind), &
        &                         conf     = pconf                   )

    case('maxwell','maxwell_2d')
      call atl_save_maxwell( eqn_name = trim(equation%eq_kind), &
        &                    conf     = pconf                   )

    case('lineareuler', 'lineareuler_2d')
      call atl_save_linearEuler( me          = equation%LinearEuler,   &
        &                        eqn_name    = trim(equation%eq_kind), &
        &                        nDimensions = equation%nDimensions,   &
        &                        conf        = pconf                   )

    case('loclineuler', 'loclineuler_1d')
      call atl_save_euler(       me          = equation%euler,         &
        &                        eqn_name    = trim(equation%eq_kind), &
        &                        conf        = pconf                   )

    case('acoustic', 'acoustic_2d')
      call atl_save_acoustic( me       = equation%acoustic,      &
        &                     eqn_name = trim(equation%eq_kind), &
        &                     conf     = pconf                   )

    case('heat_1d','heat_2d','heat_3d')
      call atl_save_heat( me       = equation%heat,          &
        &                 eqn_name = trim(equation%eq_kind), &
        &                 conf     = pconf                   )

    case default
      write(logUnit(1),*) 'Equation system ' // trim(equation%eq_kind) &
        & // ' has no support for its description in restart files yet!'

    end select

    call atl_dump_materialFun( conf     = pconf,             &
      &                        material = equation%material, &
      &                        key      = 'material',        &
      &                        varSys   = equation%varSys    )

    call aot_out_close_table(put_conf = pconf)

  end subroutine atl_eqn_write
  !*****************************************************************************


end module atl_equation_init_module
