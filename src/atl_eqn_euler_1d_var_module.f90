! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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

!> Module to configure the variables of the Euler equations.
module atl_eqn_euler_1d_var_module
  use, intrinsic :: iso_c_binding,    only: c_loc, c_f_pointer, c_ptr
  use env_module,                     only: labellen

  use tem_varSys_module,              only: tem_varSys_type,              &
    &                                       tem_varSys_init,              &
    &                                       tem_varSys_append_stateVar,   &
    &                                       tem_varSys_append_derVar,     &
    &                                       tem_varSys_proc_point,        &
    &                                       tem_varSys_proc_element,      &
    &                                       tem_varSys_proc_setparams,    &
    &                                       tem_varSys_proc_getparams,    &
    &                                       tem_varSys_proc_setupIndices, &
    &                                       tem_varSys_proc_getValOfIndex
  use tem_varMap_module,              only: tem_possible_variable_type, &
    &                                       init, append
  use tem_dyn_array_module,           only: PositionOfVal
  use tem_logging_module,             only: logUnit
  use tem_aux_module,                 only: tem_abort
  use tem_operation_var_module,       only: tem_divideVecByScal_forPoint,  &
    &                                       tem_divideVecByScal_fromIndex, &
    &                                       tem_opVar_setupIndices,        &
    &                                       tem_get_new_varSys_data_ptr

  use atl_varSys_module,              only: atl_varSys_solverData_type,    &
    &                                       atl_varSys_getStateForElement, &
    &                                       atl_varSys_getStateForPoint,   &
    &                                       atl_get_new_varSys_data_ptr,   &
    &                                       atl_varSys_setupStateIndices,  &
    &                                       atl_varSys_getStateValofIndex
  use atl_operator_module,            only: atl_op_divideVecByScal_forElement, &
    &                                       atl_opVar_setupIndices
  use atl_source_types_module,        only: atl_eqn_sourceMap_type
  use atl_equation_module,            only: atl_equations_type
  use atl_eqn_euler_derive_module,    only: atl_pressure_getElement,      &
    &                                       atl_speedOfSound_getPoint,    &
    &                                       atl_speedOfSound_getElement,  &
    &                                       atl_temperature_getPoint,     &
    &                                       atl_temperature_getElement,   &
    &                                       atl_machNumber_getPoint,      &
    &                                       atl_machNumber_getElement,    &
    &                                       atl_KineticEnergy_getPoint,   &
    &                                       atl_kineticEnergy_getElement, &
    &                                       atl_linindicator_getPoint,    &
    &                                       atl_linindicator_getElement
  use atl_eqn_euler_1d_derive_module, only: atl_pressure_1d_getPoint, &
    &                                       atl_pressure_1d_getIndex

  implicit none

  private

  public :: atl_init_euler_1d_vars
  public :: atl_append_euler_1d_consVars
  public :: atl_append_euler_1d_primVars
  public :: atl_append_euler_1d_derivedVars
  public :: atl_init_euler_1d_sourceTerms


contains


! *******************************************************************************
  !> Init the variable system for Euler (inviscid) flow simulations.
  !!
  !! The variable system describes, which variables are to be used and how
  !! they are organized in the memory.
  !! The first few variables up to the sys_mark are those, describing the
  !! state, and thus describe the output for regular restart files.
  !! Here these are the conservative variables density, momentum and energy.
  !! After the mark, there are additional values described that can be deduced
  !! from the state variables.
  subroutine atl_init_euler_1d_vars( equation, solverData )
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type) :: solverData
    ! --------------------------------------------------------------------------

    ! Initialize variable system
    call tem_varSys_init( me         = equation%varSys,            &
      &                   systemName = 'euler_1d_conservative' )

    ! Append the conservative variables
    call atl_append_euler_1d_consVars(equation, solverData)

    ! Append primitive variables
    equation%hasPrimitiveVariables = .true.
    call atl_append_euler_1d_primVars( equation%varSys,  &
      &                                equation%primVar, &
      &                                solverData        )

    ! Append derived quantities (also sets derive routines for primitive vars)
    call atl_append_euler_1d_derivedVars(equation%varSys, solverData)

    ! set the size of the temporary array used for fluxes later
    equation%temp%overSamp = 1
    equation%temp%modal = 0
    equation%temp%nodal = 2
    equation%temp%nScal = equation%varsys%nScalars

  end subroutine atl_init_euler_1d_vars
! *******************************************************************************


! *******************************************************************************
  !> Append conservative variables for Euler equations.
  !!
  !! These are density, momentum and energy here.
  subroutine atl_append_euler_1d_consVars(equation, solverData)
    ! --------------------------------------------------------------------------
    !> The equation type .
    type(atl_equations_type), intent(inout) :: equation
    !> The type to be used as method data. This type has to contain all
    !! information necessary to determine the values of the variables in the
    !! variable system.
    type(atl_varSys_solverData_type), target :: solverData
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    ! --------------------------------------------------------------------------

    allocate(equation%stateVar(3))

    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    setup_indices => atl_varSys_setupStateIndices
    get_valOfIndex => atl_varSys_getStateValofIndex

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'density',                               &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & pos            = equation%stateVar(1),                    &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex                           )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'momentum',                              &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & pos            = equation%stateVar(2),                    &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex                           )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'energy',                                &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & pos            = equation%stateVar(3),                    &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex                           )

    get_element => null()
    get_point => null()
    set_params => null()
    get_params => null()
    setup_Indices => null()
    get_valOfindex => null()

  end subroutine atl_append_euler_1d_consVars
! *******************************************************************************


! *******************************************************************************
  !> Append primitive variables for euler equation
  subroutine atl_append_euler_1d_primVars(varSys, primVar, solverData)
    ! --------------------------------------------------------------------------
    !> The Euler variable system to modify, it has to contain the conservative
    !! variables already.
    type(tem_varSys_type), intent(inout) :: varSys

    !> Indices of the primitive variables in the overall system.
    integer, allocatable, intent(out) :: primVar(:)

    !> the pointer to the data required for the varsys to fulfill all operations
    !! and derivations on the variables
    type(atl_varSys_solverData_type), target :: solverData
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    type(c_ptr) :: solver_bundle
    ! --------------------------------------------------------------------------

    ! Save var position in primVar
    allocate(primVar(3))

    ! First primitive variable is the density.
    ! It is already part of the variable system, thus we have to find it.
    primVar(1) = PositionOfVal(varSys%varname, 'density')

    ! Second primitive variable is the velocity.
    ! To compute it, density and momentum of the state are required.

    get_element => atl_op_divideVecByScal_forElement
    get_point => tem_divideVecByScal_forPoint
    setup_indices => tem_opVar_setupIndices
    get_valOfIndex => tem_divideVecByScal_fromIndex

    ! get c_ptr from solver method_Data and store it tem_varSys_op_data_type
    solver_bundle = atl_get_new_varSys_data_ptr(solverData)

    call tem_varSys_append_derVar(                                   &
      & me             = varSys,                                     &
      & varName        = 'velocity',                                 &
      & operType       = 'divive_vector_by_scalar',                  &
      & nComponents    = 1,                                          &
      & input_varname  = ['momentum', 'density '],                   &
      & method_data    = tem_get_new_varSys_data_ptr(solver_bundle), &
      & get_point      = get_point,                                  &
      & get_element    = get_element,                                &
      & set_params     = set_params,                                 &
      & get_params     = get_params,                                 &
      & setup_indices  = setup_indices,                              &
      & get_valOfIndex = get_valOfIndex,                             &
      & pos            = primVar(2)                                  )

    ! Third primitive variable is the pressure.
    ! To compute it all three state variables are required.
    get_element => atl_pressure_getElement
    get_point => atl_pressure_1d_getPoint
    set_params => null()
    get_params => null()
    setup_Indices => atl_opVar_setupIndices
    get_valOfIndex => atl_pressure_1d_getIndex


    call tem_varSys_append_derVar(                                &
      & me             = varSys,                                  &
      & varName        = 'pressure',                              &
      & nComponents    = 1,                                       &
      & input_varname  = ['density ', 'momentum', 'energy  '],    &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = primVar(3)                               )

    get_element => null()
    get_point => null()
    set_params => null()
    get_params => null()
    setup_Indices => null()
    get_valOfindex => null()

  end subroutine atl_append_euler_1d_primVars
! *******************************************************************************


! *******************************************************************************
  !> summary: append / set data of derived quantities
  subroutine atl_append_euler_1d_derivedVars( varSys, solverData )
    ! ------------------------------------------------------------------------
    !> The Euler variable system to modify. It has to contain the conservative
    !! and primitive variables already.
    type(tem_varSys_type), intent(inout) :: varSys
    !> the pointer to the data required for the varsys to fulfill all operations
    !! and derivations on the variables
    type(atl_varSys_solverData_type), target :: solverData
    ! ------------------------------------------------------------------------
    integer :: nDerivedVars, iVar, nComponents
    character(len=20), allocatable :: derVarName(:)
    character(len=labelLen), allocatable :: invar_name(:)
    character(len=labelLen) :: varname
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    type(c_ptr) :: method_data
    ! ------------------------------------------------------------------------

    wasAdded = .False.

    nDerivedVars = 6
    allocate(derVarName(nDerivedVars))
    derVarName    = [ 'speedOfSound  ', 'temperature   ', &
      &               'mach_number   ', 'mach_vector   ', &
      &               'kinetic_energy', 'linindicator  '  ]

    do iVar = 1, nDerivedVars
      varname = trim(adjustl(derVarName(iVar)))

      !for all varaiables these pointer should be nullified
      nullify(get_point, get_element, set_params, get_params, setup_indices, &
        &     get_valOfIndex)

      select case(varname)

      ! Define the speed of sound as new possible variable to derive.
      case ('speedOfSound')
        get_point => atl_speedOfSound_getPoint
        get_element => atl_speedOfSound_getElement
        setup_Indices => atl_opVar_setupIndices
        method_data = atl_get_new_varSys_data_ptr(solverData)
        nComponents = 1
        allocate(invar_name(2))
        invar_name(1) = 'pressure'
        invar_name(2) = 'density'

      ! Define the linearization indicator as new possible variable to derive.
      case ('linindicator')
        get_point => atl_linindicator_getPoint
        get_element => atl_linindicator_getElement
        setup_indices => atl_opVar_setupIndices
        method_data = atl_get_new_varSys_data_ptr(solverData)
        nComponents = 1
        allocate(invar_name(0))

      case ('temperature')
        get_point => atl_temperature_getPoint
        get_element => atl_temperature_getElement
        setup_Indices => atl_opVar_setupIndices
        method_data = atl_get_new_varSys_data_ptr(solverData)
        nComponents = 1
        allocate(invar_name(2))
        invar_name(1) = 'pressure'
        invar_name(2) = 'density'

      case ('mach_number')
        get_point => atl_machNumber_getPoint
        get_element => atl_machNumber_getElement
        setup_Indices => atl_opVar_setupIndices
        method_data = atl_get_new_varSys_data_ptr(solverData)
        nComponents = 1
        allocate(invar_name(3))
        invar_name(1) = 'density'
        invar_name(2) = 'momentum'
        invar_name(3) = 'speedOfSound'

      case ('mach_vector')
        get_point => tem_divideVecByScal_forPoint
        get_element => atl_op_divideVecByScal_forElement
        get_valOfIndex => tem_divideVecByScal_fromIndex
        setup_Indices => tem_opVar_setupIndices
        ! KM: replace solver method data into treelm method data
        method_data = atl_get_new_varSys_data_ptr(solverData)
        method_data = tem_get_new_varSys_data_ptr(method_data)
        nComponents = 1
        allocate(invar_name(2))
        invar_name(1) = 'velocity'
        invar_name(2) = 'speedOfSound'

      case ('kinetic_energy')
        get_point => atl_KineticEnergy_getPoint
        get_element => atl_kineticEnergy_getElement
        setup_Indices => atl_opVar_setupIndices
        method_data = atl_get_new_varSys_data_ptr(solverData)
        nComponents = 1
        allocate(invar_name(2))
        invar_name(1) = 'density'
        invar_name(2) = 'momentum'

      case default
        write(logUnit(1),*) 'WARNING: Unknown variable: '//trim(varname)
        cycle !go to next variable

      end select

      ! append variable to varSys
      call tem_varSys_append_derVar( me             = varSys,         &
        &                            varName        = varname,        &
        &                            nComponents    = nComponents,    &
        &                            input_varname  = invar_name,     &
        &                            method_data    = method_data,    &
        &                            get_point      = get_point,      &
        &                            get_element    = get_element,    &
        &                            set_params     = set_params,     &
        &                            get_params     = get_params,     &
        &                            setup_indices  = setup_indices,  &
        &                            get_valOfIndex = get_valOfIndex, &
        &                            wasAdded       = wasAdded        )

      if (wasAdded) then
        write(logUnit(10),*) 'Appended variable: ' // trim(varname)
      else
        call tem_abort( 'Error: variable ' // trim(varname) &
          & // ' is not added to variable system'           )
      end if

      deallocate(invar_name)

    end do

    deallocate(derVarName)

  end subroutine atl_append_euler_1d_derivedVars
! *******************************************************************************




! *******************************************************************************
  !> This routine initializes possible source variables and returns the filled
  !! up list of the poss_srcVars
  subroutine atl_init_euler_1d_sourceTerms(possVars, eval_source)
    ! --------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(out)  :: possVars
    type(atl_eqn_sourceMap_type), allocatable, intent(out) :: eval_source(:)
    ! --------------------------------------------------------------------------
    integer :: pos
    ! --------------------------------------------------------------------------

    allocate(eval_source(2))
    call init(me = possVars, length = 2)

    ! Add the spongeLayer Source term
    call append( me          = possVars,      &
      &          varName     = 'spongeLayer', &
      &          nComponents = 4,             &
      &          pos         = pos            )

    ! Add the gravitation Source term
    call append( me          = possVars,      &
      &          varName     = 'gravitation', &
      &          nComponents = 1,             &
      &          pos         = pos            )

    eval_source(pos)%do => null()

  end subroutine atl_init_euler_1d_sourceTerms
! *******************************************************************************

end module atl_eqn_euler_1d_var_module
