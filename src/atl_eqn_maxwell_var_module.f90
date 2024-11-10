! Copyright (c) 2013-2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2015 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!! equations
module atl_eqn_maxwell_var_module
  use, intrinsic :: iso_c_binding,    only: c_loc, c_ptr
  use env_module,                     only: rk

  use tem_time_module,                only: tem_time_type
  use treelmesh_module,               only: treelmesh_type
  use tem_varSys_module,              only: tem_varSys_type,              &
    &                                       tem_varSys_init,              &
    &                                       tem_varSys_append_stateVar,   &
    &                                       tem_varSys_proc_point,        &
    &                                       tem_varSys_proc_element,      &
    &                                       tem_varSys_proc_setparams,    &
    &                                       tem_varSys_proc_getparams,    &
    &                                       tem_varSys_proc_setupIndices, &
    &                                       tem_varSys_proc_getValOfIndex
  use tem_varMap_module,              only: tem_possible_variable_type, &
    &                                       init,                       &
    &                                       append
  use tem_aux_module,                 only: tem_abort

  use ply_poly_project_module,        only: ply_poly_project_type, &
    &                                       assignment(=)

  use atl_equation_module,            only: atl_equations_type
  use atl_materialFun_module,         only: atl_materialFun_type
  use atl_materialPrp_module,         only: atl_material_type, &
    &                                       atl_ConstMatIdx,   &
    &                                       atl_VarMatIdx
  use atl_varSys_module,              only: atl_varSys_solverData_type,    &
    &                                       atl_varSys_getStateForElement, &
    &                                       atl_varSys_getStateForPoint,   &
    &                                       atl_get_new_varSys_data_ptr,   &
    &                                       atl_varSys_setupStateIndices,  &
    &                                       atl_varSys_getStateValofIndex
  use atl_cube_elem_module,           only: atl_cube_elem_type
  use atl_source_types_module,        only: atl_eqn_sourceMap_type, &
    &                                       atl_source_op_type
  use atl_equation_source_module,     only: atl_equation_evaluate_source_modal,&
    &                                       atl_compute_source_interface

  implicit none

  private

  public ::                           &
    & atl_init_maxwell_vars,          &
    & atl_append_maxwell_vars,        &
    & atl_init_maxwell_sourceTerms,   &
    & atl_eval_source_currentDensity, &
    & atl_init_maxwell_material,      &
    & atl_getMaxPropSpeed


contains


! ******************************************************************************
  !> summary: init the variables for maxwell equation
  subroutine atl_init_maxwell_vars(equation, methodData)
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout)         :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type) :: methodData
    ! --------------------------------------------------------------------------
    ! initialize variable system
    call tem_varSys_init( me = equation%varSys, systemName = 'maxwell' )

    allocate(equation%stateVar(2))

    ! Append conservative Variables to variable system
    call atl_append_maxwell_vars(equation, methodData)

    equation%hasPrimitiveVariables = .false.

    ! Set values fro allocating temp flux arrays
    equation%temp%overSamp = 1
    equation%temp%modal = 0
    equation%temp%nodal = 2
    equation%temp%nScal = equation%varsys%nScalars

  end subroutine atl_init_maxwell_vars
! ******************************************************************************


! ******************************************************************************
  !> summary: append the variables for electrodynamic simulations
  subroutine atl_append_maxwell_vars(equation, methodData)
    ! --------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout)         :: equation
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

    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    setup_indices => atl_varSys_setupStateIndices
    get_valOfIndex => atl_varSys_getStateValofIndex

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'displacement_field',                    &
      & nComponents    = 3,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(1)                     )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'magnetic_field',                        &
      & nComponents    = 3,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex,                          &
      & pos            = equation%stateVar(2)                     )

  end subroutine atl_append_maxwell_vars
! ******************************************************************************


! ******************************************************************************
  !> summary: evaluate "currentDensity" source
  subroutine eval_currentDensity(rhs, source, state, constants)
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
    rhs(:,1:3) = - source(:,1:3)

  end subroutine eval_currentDensity
! ******************************************************************************


! ******************************************************************************
  subroutine atl_eval_source_currentDensity( fun, varSys, time, mesh,        &
    &                                        poly_proj, currentLevel, state, &
    &                                        material, sourcedata            )
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
    evaluate_source => eval_currentDensity

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

  end subroutine atl_eval_source_currentDensity
! ******************************************************************************


! ******************************************************************************
  !> Adds the properties of the expected source terms to the list of possible
  !! variables to extract these expected variables later on from the
  !! configuration file.
  subroutine atl_init_maxwell_material( possVars )
    ! --------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(out)  :: possVars
    ! --------------------------------------------------------------------------

    call init( me = possVars, length = 3 )

    call append( me          = possVars,       &
      &          varName     = 'permeability', &
      &          nComponents = 1               )
    call append( me          = possVars,       &
      &          varName     = 'permittivity', &
      &          nComponents = 1               )
    call append( me          = possVars,       &
      &          varName     = 'conductivity', &
      &          nComponents = 1               )

  end subroutine atl_init_maxwell_material
! ******************************************************************************


! ******************************************************************************
  !> summary: init source terms for electrodynamic simulations.
  subroutine atl_init_maxwell_sourceTerms(possVars, eval_source)
    ! --------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(out)  :: possVars
    type(atl_eqn_sourceMap_type), allocatable, intent(out) :: eval_source(:)
    ! --------------------------------------------------------------------------
    integer :: pos
    ! --------------------------------------------------------------------------

    allocate(eval_source(1))
    call init( me = possVars, length = 1 )

    call append( me          = possVars,          &
      &          varName     = 'current_density', &
      &          nComponents = 3,                 &
      &          pos         = pos                )
    eval_source(pos)%do => atl_eval_source_currentDensity

  end subroutine atl_init_maxwell_sourceTerms
! ******************************************************************************


! ******************************************************************************
  !> Determines maximum propagation speed, i.e. the speed of light depends only
  !! on material parameters.
  subroutine atl_getMaxPropSpeed(tree, materialFun, material_list)
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
    integer :: iLevel, iMat, permeaPos = 0, permittPos = 0
    ! --------------------------------------------------------------------------

    do iMat = 1, materialFun%nMat
      if( materialFun%matParNames(iMat) == 'mat_permeability' ) then
        permeaPos = sum(materialFun%nScalars(1:iMat))
      else if( materialFun%matParNames(iMat) == 'mat_permittivity' ) then
        permittPos = sum(materialFun%nScalars(1:iMat))
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

    do iLevel = tree%global%minLevel, tree%global%maxLevel

      ! Speed of light is given by: sol = 1.0 / sqrt(mu*epsi)

      ! Get the smallest value for the constant material parameters.
      material_list(iLevel)%maxPropSpeed =                            &
        & sqrt( 1.0 / minval(                                         &
        &   material_list(iLevel)%material_dat                        &
        &                        %elemMaterialData(atl_constMatIdx)   &
        &                        %materialDat(:,:,permeaPos)          &
        &   * material_list(iLevel)%material_dat                      &
        &                          %elemMaterialData(atl_constMatIdx) &
        &                          %materialDat(:,:,permittPos) ) )

      ! And now, get the smallest value for the non-constant material
      ! parameters.
      material_list(iLevel)%maxPropSpeed =                            &
        & max( material_list(iLevel)%maxPropSpeed,                    &
        &   sqrt( 1.0 / minval(                                       &
        &     material_list(iLevel)%material_dat                      &
        &                          %elemMaterialData(atl_varMatIdx)   &
        &                          %materialDat(:,:,permeaPos)        &
        &     * material_list(iLevel)%material_dat                    &
        &                            %elemMaterialData(atl_varMatIdx) &
        &                            %materialDat(:,:,permittPos) ) ) )
    end do

  end subroutine atl_getMaxPropSpeed
! ******************************************************************************


end module atl_eqn_maxwell_var_module
