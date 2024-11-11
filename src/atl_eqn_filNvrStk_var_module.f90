! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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


module atl_eqn_filNvrStk_var_module
  use, intrinsic :: iso_c_binding,      only: c_loc, c_f_pointer, c_ptr
  use env_module,                       only: rk, labelLen

  use tem_varSys_module,                only: tem_varSys_type,             &
    &                                         tem_varSys_init,             &
    &                                         tem_varSys_append_stateVar,  &
    &                                         tem_varSys_append_derVar,    &
    &                                         tem_varSys_proc_point,       &
    &                                         tem_varSys_proc_element,     &
    &                                         tem_varSys_proc_setparams,    &
    &                                         tem_varSys_proc_getparams,    &
    &                                         tem_varSys_proc_setupIndices, &
    &                                         tem_varSys_proc_getValOfIndex
  use tem_grow_array_module,            only: init, append, truncate
  use tem_dyn_array_module,             only: init, append, truncate, &
    &                                         PositionOfVal
  use tem_logging_module,               only: logUnit
  use tem_varMap_module,                only: tem_possible_variable_type, &
    &                                         init, append
  use tem_aux_module,                   only: tem_abort
  use tem_operation_var_module,         only: tem_divideVecByScal_forPoint,  &
    &                                         tem_divideVecByScal_fromIndex, &
    &                                         tem_opVar_setupIndices,        &
    &                                         tem_get_new_varSys_data_ptr
  use tem_varSys_module,                only: tem_varSys_type
  use tem_time_module,                  only: tem_time_type

  use ply_oversample_module,            only: ply_convert2oversample,   &
    &                                         ply_convertFromoversample
  use ply_poly_project_module,          only: ply_poly_project_type, &
    &                                         assignment(=),         &
    &                                         ply_poly_project_n2m,  &
    &                                         ply_poly_project_m2n
  use ply_leg_diff_module,              only: ply_calcDiff_leg_2d

  use atl_equation_module,              only: atl_equations_type
  use atl_source_types_module,          only: atl_eqn_sourceMap_type
  use atl_eqn_filNvrStk_derive_module,  only: atl_Rans_pressure_getElement, &
    &                                         atl_Rans_pressure_getPoint
  use atl_varSys_module,                only: atl_varSys_solverData_type,    &
    &                                         atl_varSys_data_type,          &
    &                                         atl_varSys_getStateForElement, &
    &                                         atl_varSys_getStateForPoint,   &
    &                                         atl_get_new_varSys_data_ptr,   &
    &                                         atl_varSys_setupStateIndices,  &
    &                                         atl_varSys_getStateValofIndex
  use atl_operator_module,              only: atl_op_divideVecByScal_forElement, &
    &                                         atl_opVar_setupIndices
  use atl_source_types_module,          only: atl_source_op_type
  use atl_cube_elem_module,             only: atl_cube_elem_type

  implicit none

  private

  public :: atl_init_rans_vars
  public :: atl_init_rans_2d_vars
  public :: atl_append_rans_consVars
  public :: atl_append_rans_primVars
  public :: atl_init_filNvrStk_sourceTerms
  public :: atl_get_pointwise_velocity_gradient_2D
  public :: atl_get_pointwise_visc_stress_tensor_2D
  public :: atl_get_lower_bound_turb_disscipation
  public :: atl_init_RANS_closure_coeffs

contains

  subroutine atl_init_RANS_closure_coeffs( equation )
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    ! --------------------------------------------------------------------------
    equation%FiltNavierStokes%rans%turb_prandtl_num =  0.9_rk
    equation%FiltNavierStokes%rans%sig_k = 0.5_rk
    equation%FiltNavierStokes%rans%beta_k = 0.09_rk
    equation%FiltNavierStokes%rans%sig_omg = 0.5_rk
    equation%FiltNavierStokes%rans%alpha_omg = 5.0/9.0
    equation%FiltNavierStokes%rans%beta_omg = 3.0/40.0
    equation%FiltNavierStokes%rans%c_mu = 1.0
    equation%FiltNavierStokes%rans%alpha = 1.0

  end subroutine atl_init_RANS_closure_coeffs

  ! *******************************************************************************


  ! *******************************************************************************
  !> Init the variable system for filtered NavierStokes equation.
  !!
  !! The variable system describes, which variables are to be used and how
  !! they are organized in the memory.
  !! The first few variables up to the sys_mark are those, describing the
  !! state, and thus describe the output for regular restart files.
  !! Here these are the conservative variables density, momentum and energy.
  !! After the mark, there are additional values described that can be deduced
  !! from the state variables.
  subroutine atl_init_Rans_vars( equation, solverData )
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type) :: solverData
    ! --------------------------------------------------------------------------
    integer :: nDim = 3
    ! --------------------------------------------------------------------------

    ! Initialize variable system
    call tem_varSys_init( me         = equation%varSys,      &
      &                   systemName = 'rans'                )

    ! Append the required conservative variables for the RANS model
    call atl_append_rans_consVars(equation, solverData)

    ! Append primitive variables
    equation%hasPrimitiveVariables = .true.
    call atl_append_rans_primVars(equation%varSys, equation%primVar,   &
     &                            solverData, nDim )

    ! Append derived quantities (also sets derive routines for primitive vars)
    ! @todo : check if all the derived vars are relevent
    !call append_euler_derivedVars(equation%varSys, solverData)

    ! for Navier Stokes 3D
    equation%temp%overSamp = 1
    equation%temp%modal = 0
    equation%temp%nodal = 3
    equation%temp%nScal = equation%varsys%nScalars

  end subroutine atl_init_Rans_vars

  ! *******************************************************************************

  subroutine atl_init_rans_2d_vars( equation, solverData )
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type) :: solverData
    ! --------------------------------------------------------------------------
    integer :: nDim = 2
    ! --------------------------------------------------------------------------

    ! Initialize variable system
    call tem_varSys_init( me         = equation%varSys,      &
      &                   systemName = 'rans_2d'             )

    ! Append the required conservative variables for the RANS model
    call atl_append_rans_consVars(equation, solverData)

    ! Append primitive variables
    equation%hasPrimitiveVariables = .true.
    call atl_append_rans_primVars(equation%varSys, equation%primVar,   &
     &                            solverData, nDim                     )

    ! Append derived quantities (also sets derive routines for primitive vars)
    call append_rans_2d_derivedVars(equation%varSys, solverData)

    ! for Rans 2d
    equation%temp%overSamp = 1
    equation%temp%modal = 0
    equation%temp%nodal = 3
    equation%temp%nScal = equation%varsys%nScalars

  end subroutine atl_init_rans_2d_vars

  ! *******************************************************************************
  subroutine atl_append_rans_consVars(equation, solverData)
    ! ---------------------------------------------------------------------------
    !> The equation type.
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), target :: solverData
    ! ---------------------------------------------------------------------------
    procedure(tem_varSys_proc_point), pointer :: get_point => null()
    procedure(tem_varSys_proc_element), pointer :: get_element => null()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    ! ---------------------------------------------------------------------------

    allocate(equation%stateVar(5))

    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    setup_indices => atl_varSys_setupStateIndices
    get_valOfIndex => atl_varSys_getStateValofIndex

    ! Append conservative variables to varSys
    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'density',                               &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(1)                     )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'momentum',                              &
      & nComponents    = equation%nDimensions,                    &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(2)                     )

    call tem_varSys_append_stateVar(                             &
      & me             = equation%varSys,                         &
      & varName        = 'energy',                                &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(3)                     )

    ! Append the Turbulent Kinetic Energy
    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'k',                                     &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(4)                     )

    ! Append the variable for specific dissipation rate
    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'omega',                                 &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(5)                     )

  end subroutine atl_append_rans_consVars

  ! *******************************************************************************

  subroutine atl_append_rans_primVars(varSys, primVar, solverData, nDim)
    ! --------------------------------------------------------------------------
    !> The variable system to append primitive variables to
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Indices of the primitive variables in the overall system.
    integer, allocatable, intent(out)     :: primVar(:)

    !> the pointer to the data required for the varsys to fulfill all operations
    !! and derivations on the variables
    type(atl_varSys_solverData_type), target :: solverData
    integer, intent(in) :: nDim
    ! --------------------------------------------------------------------------
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    type(c_ptr) :: solver_bundle
    ! --------------------------------------------------------------------------

    ! Save var position in primVar
    allocate(primVar(5))

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
      & operType       = 'divide_vector_by_scalar',                  &
      & nComponents    = nDim,                                       &
      & input_varname  = ['momentum', 'density '],                   &
      & method_data    = tem_get_new_varSys_data_ptr(solver_bundle), &
      & get_point      = get_point,                                  &
      & get_element    = get_element,                                &
      & set_params     = set_params,                                 &
      & get_params     = get_params,                                 &
      & setup_indices  = setup_indices,                              &
      & get_valOfIndex = get_valOfIndex,                             &
      & wasAdded       = wasAdded,                                   &
      & pos            = primVar(2)                                  )

    if (wasAdded) then
      write(logUnit(10),*) 'Appended variable: velocity'
    else
      write(logUnit(1),*) 'Error: variable velocity not appended'
      call tem_abort()
    end if

    ! Third primitive variable is the pressure.
    ! To compute it all three state variables are required.
    get_element => atl_Rans_pressure_getElement
    get_point => atl_Rans_pressure_getPoint
    setup_indices => atl_opVar_setupIndices

    call tem_varSys_append_derVar(                                &
      & me             = varSys,                                  &
      & varName        = 'pressure',                              &
      & nComponents    = 1,                                       &
      & input_varname  = ['density ', 'momentum',                 &
      &                   'energy  ', 'k       '],                &
      & method_data    = atl_get_new_varSys_data_ptr(solverData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = primVar(3)                               )

    if (wasAdded) then
      write(logUnit(10),*) 'Appended variable: Pressure'
    else
      write(logUnit(1),*) 'Error: variable Pressure not appended'
      call tem_abort()
    end if

    ! Fourth primitive variable is the turbulent KE.
    ! It is already part of the variable system, thus we have to find it.
    primVar(4) = PositionOfVal(varSys%varname, 'k')


    ! Fifth primitive variable is the specific dissipation rate.
    ! It is already part of the variable system, thus we have to find it.
    primVar(5) = PositionOfVal(varSys%varname, 'omega')


    get_element => null()
    get_point => null()
    set_params => null()
    get_params => null()
    setup_Indices => null()
    get_valOfindex => null()

  end subroutine atl_append_rans_primVars

! *******************************************************************************

  subroutine append_rans_2d_derivedVars( varSys, solverData )
    ! ------------------------------------------------------------------------
    type(tem_varSys_type), intent(inout) :: varSys
    !> the pointer to the data required for the varsys to fulfill all operations
    !! and derivations on the variables
    type(atl_varSys_solverData_type), target :: solverData
    ! ------------------------------------------------------------------------
!    integer :: nDerivedVars, iVar, nComponents
!    character(len=20), allocatable :: derVarName(:)
!    character(len=labelLen), allocatable :: invar_name(:)
!    character(len=labelLen) :: varname
!    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()
    ! ------------------------------------------------------------------------
!    wasAdded = .False.
!
!    !@todo Write derive routines for the primitive variables
!    nDerivedVars = 1
!    allocate(derVarName(nDerivedVars))
!
!
!    derVarName    = [ 'track_source  ']
!
!    do iVar = 1, nDerivedVars
!      varname = trim(adjustl(derVarName(iVar)))
!      select case(varname)
!!
!!      ! Define the speed of sound as new possible variable to derive.
!      case ('track_source')
!        get_point => Null()
!        get_element => Null()
!        nComponents = 6
!        allocate(invar_name(5))
!        invar_name(1) = 'density'
!        invar_name(2) = 'momentum'
!        invar_name(3) = 'energy'
!        invar_name(4) = 'k'
!        invar_name(5) = 'omega'
!
!      case default
!        write(logUnit(1),*) 'WARNING: Unknown variable: '//trim(varname)
!        cycle !go to next variable
!
!      end select
!
!      ! append variable to varSys
!      call tem_varSys_append_derVar(                                &
!        & me             = varSys,                                  &
!        & varName        = varname,                                 &
!        & nComponents    = nComponents,                             &
!        & input_varname  = invar_name,                              &
!        & method_data    = atl_get_new_varSys_data_ptr(methoddata), &
!        & get_point      = get_point,                               &
!        & get_element    = get_element,                             &
!        & wasAdded       = wasAdded                                 )
!
!      call outputAppendResult( wasAdded, varname )
!
!      deallocate(invar_name)
!
!    end do
!
!    deallocate(derVarName)

  end subroutine append_rans_2d_derivedVars
 ! *******************************************************************************


  !> Init source terms for flow simulations.
  !> This routine initializes possible source variables and returns the filled
  !! up list of the poss_srcVars
  subroutine atl_init_filNvrStk_sourceTerms(possVars, eval_source, model_type)
    ! --------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(inout)  :: possVars
    type(atl_eqn_sourceMap_type), allocatable, intent(out) :: eval_source(:)
    character(len=32), intent(in) :: model_type
    ! --------------------------------------------------------------------------
    integer :: pos
    ! --------------------------------------------------------------------------

    allocate(eval_source(1))
    call init(me = possVars%varName, length = 1)


    select case (trim(model_type))

    case('rans')
!NA!      ! Add the arbitrary Source term
!NA!      call append( me          = possVars,     &
!NA!        &          varName     = 'ransSource', &
!NA!        &          nComponents = 7,            &
!NA!        &          pos         = pos           )
!NA!      eval_source(pos)%do => eval_source_rans

    case('rans_2d')

      ! Add the source term in the
      call append( me          = possVars,       &
        &          varName     = 'rans2dsource', &
        &          nComponents = 6,              &
        &          pos         = pos             )
      eval_source(pos)%do => eval_source_rans2D

    end select


    call truncate(possVars%varName)
    call truncate(possVars%nComponents)

  end subroutine atl_init_filNvrStk_sourceTerms

  subroutine eval_source_rans2D( fun, varSys, time, mesh, poly_proj,       &
    &                            currentLevel, state, material, sourcedata )
    ! ---------------------------------------------------------------------------
    !> Description of method to update source
    class(atl_source_op_type), intent(in) :: fun
    !> variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time
    !> Current level mesh information
    type(atl_cube_elem_type), intent(in) :: mesh
    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> current level
    integer, intent(in) :: currentLevel
    !> The state in modal space.
    !! This is needed for several source terms that have to be applied to the
    !! current state
    real(kind=rk), intent(in) :: state(:,:,:)
    !> Material description for the background on the current level
    real(kind=rk), intent(in) :: material(:)
    !> The source data to update. When all source terms are added to this
    !! buffer, it is applied to the state.
    real(kind=rk), intent(inout) :: sourcedata(:,:,:)
    ! ---------------------------------------------------------------------------
    !procedure(compute_source_interface) , pointer:: evaluate_source
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:), modal_grad(:,:,:)
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:), pointVal_grad(:,:,:)
    real(kind=rk), allocatable :: rhsNodal(:,:), rhsModal(:,:), rhs(:,:)
    real(kind=rk) :: rey_stress_tensor(2,2)
    integer :: iElem, iDir
    integer :: iPoint
    real(kind=rk) :: beta_k, omega_r, alp_omg, beta_omg, sig_omg
    real(kind=rk) :: mu_turb, mu, c_mu, alpha,  source(2), omega, k_bar
    real(kind=rk) :: velGrad(2,2), ViscStressTensor(2,2)
    type(atl_varSys_data_type), pointer :: fPtr
    type(atl_equations_type), pointer :: eqn
    real(kind=rk) :: d_omega1, d_omega2
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(fun%srcTerm_varPos)%method_Data, fPtr )
    eqn => fPtr%solverData%equationPtr

    beta_k   = eqn%FiltNavierStokes%rans%beta_k
    beta_omg = eqn%FiltNavierStokes%rans%beta_omg
    alp_omg  = eqn%FiltNavierStokes%rans%alpha_omg
    c_mu     = eqn%FiltNavierStokes%rans%c_mu
    alpha    = eqn%FiltNavierStokes%rans%alpha
    mu       = eqn%NavierStokes%mu
    sig_omg  = eqn%FiltNavierStokes%rans%sig_omg


    allocate( modalCoeffs (poly_proj%Body_2D%OverSamp_dofs ,varSys%nScalars ))
    allocate( modal_grad (poly_proj%Body_2D%OverSamp_dofs,varSys%nScalars,2 ))
    allocate( pointVal (poly_proj%Body_2D%nQuadPoints,varSys%nScalars )     )
    allocate( pointVal_grad (poly_proj%Body_2D%nQuadPoints,varSys%nScalars,2 ))
    allocate( rhsNodal (poly_proj%Body_2D%nQuadPoints, varSys%nScalars )    )
    allocate( rhsModal (poly_proj%Body_2D%OverSamp_dofs, varSys%nScalars )  )
    allocate( rhs (poly_proj%body_2d%nDofs, varSys%nScalars ))

    rhsNodal = 0.0_rk
    do iElem = 1, fun%elems(currentLevel)%nElems

      ! --> modal space
      call ply_convert2oversample(state       = state(iElem,:,:),  &
        &                         poly_proj   = poly_proj,         &
        &                         nDim        = 2,                 &
        &                         modalCoeffs = modalCoeffs        )
      ! --> oversamp modal space

      ! Calculate the gradient of modal Coeffs
      call ply_calcDiff_leg_2d( legCoeffs     = modalCoeffs,               &
        &                       legCoeffsDiff = modal_grad,                &
        &                       maxPolyDegree = poly_proj%oversamp_degree, &
        &                       nVars         = varsys%nScalars,           &
        &                       elemLength    = mesh%length                )

      ! Now, we transform the modal representation of this element to nodal
      ! space by making use of fast polynomial transformations (FPT)
      call ply_poly_project_m2n(me = poly_proj,                  &
       &                       dim = 2 ,                         &
       &                       nVars = varSys%nScalars,          &
       &                       nodal_data= pointVal,             &
       &                       modal_data= modalCoeffs           )
      ! --> oversamp nodal space

      ! Now, transform modal gradient representation of this element to
      ! nodal space by making use of fast polynomial transformations (FPT)
      do iDir = 1, 2
        call ply_poly_project_m2n(me   = poly_proj,               &
          &                 dim        = 2,                       &
          &                 nVars      = varSys%nScalars,         &
          &                 nodal_data = pointVal_grad(:,:,iDir), &
          &                 modal_data = modal_grad(:,:,iDir)     )
      end do

      do iPoint = 1, poly_proj%Body_2D%nQuadPoints


        call atl_get_pointwise_velocity_gradient_2D( &
          & pointVal_grad(iPoint,:,:),               &
          & pointVal(iPoint,:),                      &
          & velGrad                                  )

        call atl_get_pointwise_visc_stress_tensor_2D(         &
          &          velGrad = velGrad,                       &
          &          S       = ViscStressTensor               )

        ! get the value of \omega_r
        call atl_get_lower_bound_turb_disscipation(                      &
          &          S       = ViscStressTensor,                         &
          &        c_mu      = c_mu,                                     &
          &        omega     = pointVal(iPoint, 6 )/ pointVal(iPoint,1), &
          &        omega_r   = omega_r                                   )

        ! Calculate the limited turbulent eddy viscosity
        mu_turb = alpha * max(pointVal(iPoint,5),0.0_rk)* exp(-omega_r)

        !Careful - This routine returns the reynolds tensor divided by
        ! the limited "k" (\bar{k})
        call get_rans_reynolds_tensor_2D(                      &
          &    state             = pointVal(iPoint,:),         &
          &    alpha             = alpha,                      &
          &    omega_r           = omega_r,                    &
          &    velGrad           = velGrad,                    &
          &    rey_stress_tensor = rey_stress_tensor           )

        omega = exp( pointVal(iPoint,6)/ pointVal(iPoint,1) )
        k_bar = max(pointVal(iPoint,5)/pointVal(iPoint,1),0.0_rk)
        d_omega1 = (pointVal_grad(iPoint,6,1) - pointVal(iPoint,6)         &
          &           *pointval_grad(iPoint,1,1)/pointVal(iPoint,1)    )   &
          &              *pointVal(iPoint,1)
        d_omega2 = (pointVal_grad(iPoint,6,2) - pointVal(iPoint,6)         &
          &           *pointval_grad(iPoint,1,2)/pointVal(iPoint,1)    )   &
          &              / pointVal(iPoint,1)

        ! Note: The first four comonents of the source term = 0
        ! This is the source corresponding to the fifth var i.e turbulent KE
        ! the multiplication with k_bar as the rey_stress_tensor used here is
        ! divided by k_bar
        source(1) = rey_stress_tensor(1,1)*velGrad(1,1)*k_bar &
          &       + rey_stress_tensor(1,2)*velGrad(1,2)*k_bar &
          &       + rey_stress_tensor(2,1)*velGrad(2,1)*k_bar &
          &       + rey_stress_tensor(2,2)*velGrad(2,2)*k_bar &
          &       - beta_k*max(pointVal(iPoint,5),0.0_rk)     &
                          *exp(omega_r)

        ! This is the source corresponding to the sixth conservative var
        ! i.e specific dissipation rate
        source(2) = (                                         &
          &   rey_stress_tensor(1,1)*velGrad(1,1)             &
          & + rey_stress_tensor(1,2)*velGrad(1,2)             &
          & + rey_stress_tensor(2,1)*velGrad(2,1)             &
          & + rey_stress_tensor(2,2)*velGrad(2,2) )           &
          &   *alp_omg*omega                                  &
          & - beta_omg*pointVal(iPoint,1)*exp(omega_r)        &
          & + (mu + sig_omg*mu_turb)                          &
          & * ( d_omega1*d_omega1                             &
          & +   d_omega2*d_omega2                             )

          !write(*,*) "**>", iPoint, source, pointVal(iPoint,:)
          ! Debug statement
          rhsNodal(iPoint,5:6) = rhsNodal(iPoint,5:6) +  source(1:2)

      end do

      ! Convert everything back to modal space
      call ply_poly_project_n2m( me         = poly_proj,       &
        &                        dim        = 2,               &
        &                        nVars      = varSys%nScalars, &
        &                        nodal_data = rhsNodal(:,:),   &
        &                        modal_data = rhsModal(:,:)    )

      ! --> oversamp modal space
      call ply_convertFromoversample( modalCoeffs = rhsModal(:,:), &
        &                             poly_proj   = poly_proj,     &
        &                             nDim        = 2,             &
        &                             state       = rhs(:,:)       )
      ! --> oversampled modal space

      ! Add rhsmodal to the source representation
      sourcedata(iElem,:,:) = sourcedata(iElem,:,:) + rhs(:,:)


    end do ! Elems Loop

    deallocate( modalCoeffs )
    deallocate( modal_grad )
    deallocate( pointVal )
    deallocate( pointVal_grad )
    deallocate( rhsNodal )
    deallocate( rhsModal )
    deallocate( rhs )


  end subroutine eval_source_rans2D

! *******************************************************************************

  ! Calculates the velocity gradient from the available state gradient
  subroutine atl_get_pointwise_velocity_gradient_2D( state_grad, state, res )
    ! ---------------------------------------------------------------------------
    real(kind=rk), intent(in) :: state_grad(:,:)
    ! The gradient of the state
    real(kind=rk), intent(in) :: state(:)
    ! Reynolds_stress_tensor
    real(kind=rk), intent(inout) :: res(:,:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: dens_sq
    ! ---------------------------------------------------------------------------
    dens_sq = state(1)*state(1)
    ! fill up the four velocity gradients.
    res(1,1) = state_grad(2,1)/ state(1) - state(2) * state_grad(1,1)/ dens_sq
    res(2,1) = state_grad(3,1)/ state(1) - state(3) * state_grad(1,1)/ dens_sq
    res(1,2) = state_grad(2,2)/ state(1) - state(2) * state_grad(1,2)/ dens_sq
    res(2,2) = state_grad(3,2)/ state(1) - state(3) * state_grad(1,2)/ dens_sq


  end subroutine atl_get_pointwise_velocity_gradient_2D

! *******************************************************************************

  ! Takes the pointwise velocity gradient as the input and evaluates
  ! the four components of viscous stress tensor
  subroutine atl_get_pointwise_visc_stress_tensor_2D( velGrad, S )
    ! ---------------------------------------------------------------------------
    ! The gradient of the state
    real(kind=rk), intent(in) :: velGrad(2,2)
    ! Viscous stress tensor - the output
    real(kind=rk), intent(inout) :: S(2,2)
    ! ---------------------------------------------------------------------------

    S(1,1) = 4.0* velGrad(1,1)/ 3.0 - 2.0* velGrad(2,2)/ 3.0
    S(1,2) = velGrad(1,2) + velGrad(2,1)
    S(2,1) = S(1,2)
    S(2,2) = 4.0 * velGrad(2,2)/ 3.0 - 2.0* velGrad(1,1)/ 3.0

  end subroutine atl_get_pointwise_visc_stress_tensor_2D

! *******************************************************************************

  ! takes pointwise input and evaluates the lower bound of the turbulent
  ! dissipation to ensure the positivity of normal turbulent stresses.
  ! The fulfilement of the inequality solved in this routine is the same as
  ! one proposed by Wilcox for the k-\omega turbulence model and also used
  ! in the Hartmann paper. Here we just solve this quadratic inequality
  subroutine atl_get_lower_bound_turb_disscipation( S, c_mu, omega, Omega_r )
    ! Viscous stress tensor
    real(kind=rk), intent(in) :: S(2,2)
    ! The constant c_mu
    real(kind=rk), intent(in) :: c_mu
    ! The state omega
    real(kind=rk), intent(in) :: omega
    ! The lower bound on omega
    real(kind=rk), intent(out) :: omega_r
    ! ---------------------------------------------------------------------------
    ! The constants needed
    real(kind=rk) :: C3,C4
    ! The omegas
    real(kind=rk) :: omega_r0
    real(kind=rk) :: omg1, omg2, omg3
    ! ---------------------------------------------------------------------------

    if (S(1,1) > 0.0_rk) then
      omg1 =   3.0 *c_mu * (S(1,1) ) / 2.0
    else
      omg1 = 0.0_rk
    end if

    if (S(2,2) > 0.0_rk) then
      omg2 =   3.0 *c_mu * (S(2,2) ) / 2.0
    else
      omg2 = 0.0_rk
    end if

    C3 = - 3.0 *c_mu * ( S(1,1) + S(2,2))/2.0
    C4 =   9.0 *c_mu * c_mu * (S(1,1)* S(2,2) - S(1,2)*S(1,2))/ 4.0

    if ( (C3 * C3 - 4.0*C4) > 0 ) then
      omg3 = - (C3 + sqrt( C3 * C3 - 4.0*C4) )/ 2.0
    else
      omg3 = 0.0_rk
    end if

    omega_r0 = max(omg1, omg2, omg3, tiny(omega))
    omega_r0 = log(omega_r0)

    omega_r = max(omega, omega_r0)


  end subroutine atl_get_lower_bound_turb_disscipation


! *******************************************************************************
  ! This routine calculates the reynolds tensor and returns the result divided
  ! by \bar{k} = \frac{\bar{\mu_t}}{\alpha \rho e^{\tilde{\omega_r}}}
  subroutine get_rans_reynolds_tensor_2D( state, alpha, omega_r, velGrad, &
    &                                     rey_stress_tensor )
    ! ---------------------------------------------------------------------------
    !! current state
    real(kind=rk), intent(in) :: state(:)
    ! \alpha for the evaluation of limited \mu_t (turbulent eddy viscosity)
    real(kind=rk), intent(in) :: alpha
    ! lower bound of \omega = \omega_r for the evaluation of limited
    ! \mu_t (turbulent eddy viscosity)
    real(kind=rk), intent(in) :: omega_r
    ! The velocity gradient
    real(kind=rk), intent(in) :: velGrad(2,2)
    ! Reynolds_stress_tensor
    real(kind=rk), intent(inout) :: rey_stress_tensor(:,:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: fact
    ! ---------------------------------------------------------------------------

    fact = alpha*state(1)*exp(-omega_r)


    rey_stress_tensor(1,1) =  (4.0*velGrad(1,1)  - 2.0*velGrad(2,2) )  &
      &                          *fact / 3.0   - 2.0 * state(1) / 3.0

    rey_stress_tensor(1,2) = fact*(velGrad(1,2) + velGrad(2,1))

    rey_stress_tensor(2,1) = rey_stress_tensor(1,2)

    rey_stress_tensor(2,2) = fact*( 4.0*velGrad(2,2) - 2.0*velGrad(1,1)) / 3.0 &
      &                        - 2.0 * state(1) / 3.0

  end subroutine get_rans_reynolds_tensor_2D


! *******************************************************************************

end module atl_eqn_filNvrStk_var_module
