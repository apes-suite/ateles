! Copyright (c) 2013-2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Jana Gericke <jana.gericke@student.uni-siegen.de>
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

!> summary: module to configure information about the variables of the maxwell
!!equations with divergence correction
module atl_eqn_maxwelldivcorr_var_module
  use, intrinsic :: iso_c_binding,  only: c_loc, c_ptr
  use env_module,                   only: rk

  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use tem_varSys_module,            only: tem_varSys_type,              &
    &                                     tem_varSys_init,              &
    &                                     tem_varSys_append_stateVar,   &
    &                                     tem_varSys_proc_element,      &
    &                                     tem_varSys_proc_point,        &
    &                                     tem_varSys_proc_element,      &
    &                                     tem_varSys_proc_setparams,    &
    &                                     tem_varSys_proc_getparams,    &
    &                                     tem_varSys_proc_setupIndices, &
    &                                     tem_varSys_proc_getValOfIndex
  use tem_varMap_module,            only: tem_possible_variable_type, &
    &                                     init,                       &
    &                                     append
  use tem_grow_array_module,        only: init,   &
    &                                     append, &
    &                                     truncate
  use tem_dyn_array_module,         only: init,    &
    &                                     append,  &
    &                                     truncate
  use tem_logging_module,           only: logUnit
  use tem_aux_module,               only: tem_abort

  use ply_poly_project_module,      only: ply_poly_project_type, &
    &                                     assignment(=)

  use atl_eqn_maxwell_var_module,   only: atl_append_maxwell_vars, &
    &                                     atl_eval_source_currentDensity
  use atl_equation_module,          only: atl_equations_type
  use atl_varSys_module,            only: atl_varSys_solverData_type,    &
    &                                     atl_varSys_getStateForElement, &
    &                                     atl_varSys_getStateForPoint,   &
    &                                     atl_get_new_varSys_data_ptr,   &
    &                                     atl_varSys_setupStateIndices,  &
    &                                     atl_varSys_getStateValofIndex
  use atl_source_types_module,      only: atl_eqn_sourceMap_type, &
    &                                     atl_source_op_type
  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_equation_source_module,   only: atl_equation_evaluate_source_modal, &
    &                                     atl_compute_source_interface
  use atl_materialFun_module,       only: atl_materialFun_type
  use atl_materialPrp_module,       only: atl_material_type, &
    &                                     atl_ConstMatIdx,   &
    &                                     atl_VarMatIdx

  implicit none

  private

  public ::                                &
    & atl_init_maxwellDivCorr_vars,        &
    & atl_append_maxwellDivCorr_vars,      &
    & atl_init_maxwellDivCorr_sourceTerms, &
    & atl_init_maxwellDivCorr_material,    &
    & atl_getMaxPropSpeedDivCor


contains


! ******************************************************************************
  !> summary: init the variables for maxwell equation
  subroutine atl_init_maxwellDivCorr_vars(equation, methodData)
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type) :: methodData
    ! --------------------------------------------------------------------------
    ! initialize variable system
    call tem_varSys_init( me = equation%varSys, systemName = 'maxwellDivCorr' )

    ! Append conservative Variables to variable system
    call atl_append_maxwellDivCorr_vars(equation, methodData)

    equation%hasPrimitiveVariables = .false.

    ! Set values fro allocating temp flux arrays
    equation%temp%overSamp = 1
    equation%temp%modal = 0
    equation%temp%nodal = 2
    equation%temp%nScal = equation%varsys%nScalars

  end subroutine atl_init_maxwellDivCorr_vars
! ******************************************************************************


! ******************************************************************************
  !> summary: append the variables for electrodynamic simulations that include
  !! divergence cleaning.
  subroutine atl_append_maxwellDivCorr_vars(equation, methodData)
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), target :: methodData
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer &
      & :: get_valOfIndex => NULL()
    ! --------------------------------------------------------------------------
    logical :: wasAdded
    ! --------------------------------------------------------------------------

    allocate(equation%stateVar(4))
    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    setup_indices => atl_varSys_setupStateIndices
    get_valOfIndex => atl_varSys_getStateValofIndex

    ! Append maxwell variables displacement field and magentic field
    ! before correction variables
    call atl_append_maxwell_vars(equation, methodData)

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'electric_correction',                   &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(3),                    &
      & wasAdded       = wasAdded                                 )

    if (wasAdded) then
      write(logUnit(10),*) 'Appended variable: electric_correction'
    else
      write(logUnit(1),*) 'Error: variable electric_correction not appended'
      call tem_abort()
    end if

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'magnetic_correction',                   &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(4),                    &
      & wasAdded       = wasAdded                                 )
    if (wasAdded) then
      write(logUnit(10),*) 'Appended variable: magnetic_correction'
    else
      write(logUnit(1),*) 'Error: variable magnetic_correction not appended'
      call tem_abort()
    end if

  end subroutine atl_append_maxwellDivCorr_vars
! ******************************************************************************


! ******************************************************************************
  !> summary: evaluate "currentDensity" source
  subroutine eval_charge(rhs, source, state, constants)
    ! --------------------------------------------------------------------------
    !> The Right Hand side to be updated
    real(kind=rk), intent(inout) :: rhs(:,:)
    !> The source data to be used
    real(kind=rk), intent(in) :: source(:,:)
    !> The state in the modal form
    real(kind=rk), intent(in) :: state(:,:)
    !> the constants required for the evaluation of source
    real(kind = rk ), intent(in) :: constants(:)
    ! --------------------------------------------------------------------------
    rhs = 0.0_rk
    ! Compute RHS using the modal values of source and state

    ! current density is relevant in Ampere's law:
    ! magnetic field
    rhs(:,7) = source(:,1)

  end subroutine eval_charge
! ******************************************************************************


! ******************************************************************************
  !> summary: evaluate "charge" source
  subroutine eval_source_charge( fun, varSys, time, mesh, poly_proj,       &
    &                            currentLevel, state, material, sourcedata )
    !---------------------------------------------------------------------------

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

    !> Material description for the complete domain. Used for evaluation of some
    !! source terms.
    real(kind=rk), intent(in) :: material(:)

    !> The source data to update. When all source terms are added to this
    !! buffer, it is applied to the state.
    real(kind=rk), intent(inout) :: sourcedata(:,:,:)
    ! --------------------------------------------------------------------------
    procedure(atl_compute_source_interface) , pointer:: evaluate_source
    ! --------------------------------------------------------------------------

    ! Set the function pointer for the evaluation of arbitray source
    evaluate_source => eval_charge

    ! Call the common function for updating the sourceData
    call atl_equation_evaluate_source_modal(    &
      & fun          = fun,                     &
      & varSys       = varSys,                  &
      & currentLevel = currentLevel,            &
      & nDim         = 3,                       &
      & time         = time,                    &
      & eval_rhs     = evaluate_source,         &
      & state        = state,                   &
      & poly_proj    = poly_proj,               &
      & polyProjBody = poly_proj%body_3d,       &
      & sourceData   = sourceData               )

  end subroutine eval_source_charge
! ******************************************************************************


! ******************************************************************************
  !> summary: init source terms for electrodynamic simulations with divergence
  !! cleaning
  subroutine atl_init_maxwellDivCorr_sourceTerms(possVars, eval_source)
    ! --------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(out)  :: possVars
    type(atl_eqn_sourceMap_type), allocatable, intent(out) :: eval_source(:)
    ! --------------------------------------------------------------------------
    integer :: pos
    ! --------------------------------------------------------------------------

    allocate(eval_source(2))
    call init(me = possVars, length = 2)


    call append( me          = possVars,         &
      &          varName     = 'currentDensity', &
      &          nComponents = 3,                &
      &          pos         = pos               )
    eval_source(pos)%do => atl_eval_source_currentDensity

    call append( me          = possVars, &
      &          varName     = 'charge', &
      &          nComponents = 1,        &
      &          pos         = pos       )
    eval_source(pos)%do => eval_source_charge

  end subroutine atl_init_maxwellDivCorr_sourceTerms
! ******************************************************************************


! ******************************************************************************
  !> Adds the properties of the expected source terms to the list of possible
  !! variables to extract these expected variables later on from the
  !! configuration file.
  subroutine atl_init_maxwellDivCorr_material( possVars )
    ! --------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(out)  :: possVars
    ! --------------------------------------------------------------------------

    call init(me = possVars%varName, length = 4)

    call append( me = possVars, varName = 'permeability', nComponents = 1 )
    call append( me = possVars, varName = 'permittivity', nComponents = 1 )
    call append( me = possVars, varName = 'gam', nComponents = 1 )
    call append( me = possVars, varName = 'chi', nComponents = 1 )

  end subroutine atl_init_maxwellDivCorr_material
! ******************************************************************************


! ******************************************************************************
  !> Determines maximum propagation speed for Maxwell equation
  !! with divergence cleaning (hyperbolic), i.e. the speed of light depends only
  !! on material parameters.
  subroutine atl_getMaxPropSpeedDivCor(tree, materialFun, material_list)
    ! --------------------------------------------------------------------------
    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(in) :: tree
    !> Information about the material parameters. Used to figure out the
    !! order of the material paramters as well as the number of components.
    type(atl_materialFun_type), intent(in) :: materialFun
    !> The description of the material properties. The compute lists in the
    !! material description is filled up by calling this subroutine.
    type(atl_material_type), intent(inout) :: material_list( &
      & tree%global%minLevel:tree%global%maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iMat
    integer :: permeaPos = 0, permittPos = 0, gamPos = 0, chiPos = 0
    ! --------------------------------------------------------------------------

    do iMat = 1, materialFun%nMat
      if( materialFun%matParNames(iMat) == 'mat_permeability' ) then
        permeaPos = sum(materialFun%nScalars(1:iMat))
      else if( materialFun%matParNames(iMat) == 'mat_permittivity' ) then
        permittPos = sum(materialFun%nScalars(1:iMat))
      else if( materialFun%matParNames(iMat) == 'mat_gam' ) then
        gamPos = sum(materialFun%nScalars(1:iMat))
      else if( materialFun%matParNames(iMat) == 'mat_chi' ) then
        chiPos = sum(materialFun%nScalars(1:iMat))
      end if
    end do

    if( permeaPos == 0 ) then
      call tem_abort( "Permeability not found. Can't calculate propagation" &
        & // " speed, stopping..."                                          )
    end if
    if( permittPos == 0 ) then
      call tem_abort( "Permittivity not found. Can't calculate propagation" &
        & // " speed, stopping..."                                          )
    end if
    if( gamPos == 0 ) then
      call tem_abort( "Gam not found. Can't calculate propagation" &
        & // " speed, stopping..."                                 )
    end if
    if( chiPos == 0 ) then
      call tem_abort( "Chi not found. Can't calculate propagation" &
        & // " speed, stopping..."                                 )
    end if

    do iLevel = tree%global%minLevel, tree%global%maxLevel

      ! Speed of light is given by: sol = 1.0 / sqrt(mu*epsi)
      ! The wave speed for Maxwell with hyperbolic divergence cleaning are c,
      ! gamma*c and chi*c

      ! Get the fastest waves for the constant material parameters.
      ! ... c
      material_list(iLevel)%maxPropSpeed =                                &
        & maxval(                                                         &
        &   sqrt(                                                         &
        &     1.0                                                         &
        &     / ( material_list(iLevel)%material_dat                      &
        &                              %elemMaterialData(atl_constMatIdx) &
        &                              %materialDat(:,:,permeaPos)        &
        &       * material_list(iLevel)%material_dat                      &
        &                              %elemMaterialData(atl_constMatIdx) &
        &                              %materialDat(:,:,permittPos) ) ) )
      ! ... c*gamma
      material_list(iLevel)%maxPropSpeed =                                  &
        & max(                                                              &
        &   material_list(iLevel)%maxPropSpeed,                             &
        &   maxval(                                                         &
        &     sqrt(                                                         &
        &       1.0                                                         &
        &       / ( material_list(iLevel)%material_dat                      &
        &                                %elemMaterialData(atl_constMatIdx) &
        &                                %materialDat(:,:,permeaPos)        &
        &         * material_list(iLevel)%material_dat                      &
        &                                %elemMaterialData(atl_constMatIdx) &
        &                                %materialDat(:,:,permittPos) )     &
        &         * material_list(iLevel)%material_dat                      &
        &                                %elemMaterialData(atl_constMatIdx) &
        &                                %materialDat(:,:,gamPos) ) ) )
      ! ... c*chi
      material_list(iLevel)%maxPropSpeed =                                  &
        & max(                                                              &
        &   material_list(iLevel)%maxPropSpeed,                             &
        &   maxval(                                                         &
        &     sqrt(                                                         &
        &       1.0                                                         &
        &       / ( material_list(iLevel)%material_dat                      &
        &                                %elemMaterialData(atl_constMatIdx) &
        &                                %materialDat(:,:,permeaPos)        &
        &         * material_list(iLevel)%material_dat                      &
        &                                %elemMaterialData(atl_constMatIdx) &
        &                                %materialDat(:,:,permittPos) )     &
        &         * material_list(iLevel)%material_dat                      &
        &                                %elemMaterialData(atl_constMatIdx) &
        &                                %materialDat(:,:,chiPos) ) ) )

      ! And now, get the fastest waves for the non-constant material parameters.
      ! ... c
      material_list(iLevel)%maxPropSpeed =                                  &
        & max(                                                              &
        &   material_list(iLevel)%maxPropSpeed,                             &
        &   maxval(                                                         &
        &     sqrt(                                                         &
        &       1.0                                                         &
        &       / ( material_list(iLevel)%material_dat                      &
        &                              %elemMaterialData(atl_varMatIdx)     &
        &                              %materialDat(:,:,permeaPos)          &
        &         * material_list(iLevel)%material_dat                      &
        &                                %elemMaterialData(atl_varMatIdx)   &
        &                                %materialDat(:,:,permittPos) ) ) ) )
      ! ... c*gamma
      material_list(iLevel)%maxPropSpeed =                                &
        & max(                                                            &
        &   material_list(iLevel)%maxPropSpeed,                           &
        &   maxval(                                                       &
        &     sqrt(                                                       &
        &       1.0                                                       &
        &       / ( material_list(iLevel)%material_dat                    &
        &                                %elemMaterialData(atl_varMatIdx) &
        &                                %materialDat(:,:,permeaPos)      &
        &         * material_list(iLevel)%material_dat                    &
        &                                %elemMaterialData(atl_varMatIdx) &
        &                                %materialDat(:,:,permittPos) )   &
        &         * material_list(iLevel)%material_dat                    &
        &                                %elemMaterialData(atl_varMatIdx) &
        &                                %materialDat(:,:,gamPos) ) ) )
      ! ... c*chi
      material_list(iLevel)%maxPropSpeed =                                &
        & max(                                                            &
        &   material_list(iLevel)%maxPropSpeed,                           &
        &   maxval(                                                       &
        &     sqrt(                                                       &
        &       1.0                                                       &
        &       / ( material_list(iLevel)%material_dat                    &
        &                                %elemMaterialData(atl_varMatIdx) &
        &                                %materialDat(:,:,permeaPos)      &
        &       *   material_list(iLevel)%material_dat                    &
        &                                %elemMaterialData(atl_varMatIdx) &
        &                                %materialDat(:,:,permittPos) )   &
        &       *   material_list(iLevel)%material_dat                    &
        &                                %elemMaterialData(atl_varMatIdx) &
        &                                %materialDat(:,:,chiPos) ) ) )

    end do

  end subroutine atl_getMaxPropSpeedDivCor
! ******************************************************************************


end module atl_eqn_maxwelldivcorr_var_module
