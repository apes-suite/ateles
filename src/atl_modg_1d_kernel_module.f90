! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt and Harald Klimach for
! University of Siegen.
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
!
! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the number of degrees of freedom for Q polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * A variable to store the number of degrees of freedom for a P tensor product
!   polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.

! The x ansatz degree is turned into the degree of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x ansatz degree is turned into the degree of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
!> author: Jens Zudrop
!! Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
!!
!!@todo HK: Split this module into three: a common one, a Q part and a P part.
module atl_modg_1d_kernel_module
  use env_module,               only: rk, long_k

  use tem_aux_module,           only: tem_abort
  use treelmesh_module,         only: treelmesh_type
  use tem_element_module,       only: eT_fluid
  use tem_param_module,         only: q__E, q__W, qAxis
  use tem_time_module,          only: tem_time_type
  use tem_element_module,       only: eT_fluid
  use tem_faceData_module,      only: tem_dirToFace_map, tem_left, tem_right
  use tem_logging_module,       only: logUnit
  use tem_comm_module,          only: tem_commPattern_type
  use tem_comm_env_module,      only: tem_comm_env_type

  use ply_modg_basis_module,    only: ply_faceValLeftBndAns,       &
    &                                 ply_faceValLeftBndTest,      &
    &                                 ply_faceValRightBndTest,     &
    &                                 ply_faceValLeftBndgradTest,  &
    &                                 ply_faceValRightBndgradTest, &
    &                                 ply_scalProdDualLeg,         &
    &                                 ply_scalProdDualLegDiff
  use ply_dof_module,           only: Q_space
  use ply_poly_project_module,  only: ply_poly_project_type, assignment(=)
  use ply_leg_diff_module,      only: ply_calcDiff_leg_1d
  use atl_modg_1d_scheme_module,only: atl_modg_1d_scheme_type
  use atl_equation_module,      only: atl_equations_type
  use atl_kerneldata_module,    only: atl_kerneldata_type, &
    &                                 atl_init_kerneldata, &
    &                                 atl_init_statedata,  &
    &                                 atl_statedata_type
  use atl_cube_elem_module,     only: atl_cube_elem_type
  use atl_parallel_module,      only: atl_init_parallel_module
  use atl_scheme_module,        only: atl_scheme_type
  use atl_elemental_time_integration_module,                      &
    &                           only: atl_elemental_timestep_vec, &
    &                                 atl_timestep_type
  use atl_source_types_module,  only: atl_source_type
  use atl_facedata_module,      only: atl_facedata_type, &
    &                                 atl_facerep_type,  &
    &                                 atl_elemfaceToNormal_prp
  use atl_boundary_module,      only: atl_level_boundary_type, &
    &                                 atl_get_numBndElems
  use atl_penalization_module,  only: atl_penalizationData_type
  use atl_materialPrp_module,   only: atl_material_type
  use atl_materialIni_module,   only: atl_update_materialParams

  implicit none

  private

  public :: atl_init_modg_1d_kernel,        &
    & atl_preprocess_modg_1d_kernel,        &
    & atl_modg_1d_ensure_pos_face,          &
    & atl_modg_1d_invMassMatrix,            &
    & atl_modg_1d_modalVolToModalFace,      &
    & atl_modg_1d_project_testFunc,         &
    & atl_modg_1d_project_PhysFlux_testFunc

contains

  subroutine atl_preprocess_modg_1d_kernel(                                &
    &          currentlevel, minLevel, maxLevel,                           &
    &          equation, statedata, mesh_list, boundary_list, scheme_list, &
    &          poly_proj_material, material_list, commPattern, proc        )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: currentLevel
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    type(atl_equations_type), intent(inout) :: equation
    type(atl_statedata_type), intent(inout) :: statedata
    type(atl_level_boundary_type), intent(in) :: &
      & boundary_list(minlevel:maxLevel)
    type(atl_cube_elem_type), intent(inout) :: mesh_list(minlevel:maxlevel)
    type(atl_scheme_type), intent(inout) :: scheme_list(minlevel:maxlevel)
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxlevel)
    type(ply_poly_project_type), intent(inout) :: poly_proj_material
    !> mpi communication environment with mpi communicator
    type(tem_comm_env_type), intent(inout) :: proc
    !> mpi communication pattern type
    type(tem_commPattern_type), intent(inout) :: commPattern
    ! --------------------------------------------------------------------------

    call atl_update_materialParams( equation    = equation,                    &
      &                             mesh        = mesh_list(currentLevel),     &
      &                             scheme      = scheme_list(currentlevel),   &
      &                             boundary    = boundary_list(currentlevel), &
      &                             material    = material_list(currentlevel), &
      &                             time        = statedata%local_time,        &
      &                             poly_proj   = poly_proj_material,          &
      &                             proc        = proc,                        &
      &                             commPattern = commPattern                  )


  end subroutine atl_preprocess_modg_1d_kernel


  !> Initiate the MODG kernel for cubic elements on all levels.
  subroutine atl_init_modg_1d_kernel(                                     &
    &          kerneldata_list, statedata_list, statedata_stab_list,      &
    &          mesh_list, scheme_list, boundary_list, boundary_stab_list, &
    &          equation, tree, time, commPattern,                         &
    &          need_deviation                                             )
    ! --------------------------------------------------------------------------
    !> Array of kerneldata_types, one for each level. They will be initialized.
    type(atl_kerneldata_type), allocatable :: kerneldata_list(:)
    !> Array of statedata_types, one for each level.
    type(atl_statedata_type), allocatable :: statedata_list(:)
    type(atl_statedata_type), allocatable :: statedata_stab_list(:,:)
    !> List of cubic meshes. One entry per level.
    type(atl_cube_elem_type), allocatable :: mesh_list(:)
    !> List of schemes on the levels.
    type(atl_scheme_type), allocatable :: scheme_list(:)
    !> The boundary description for the faces on the current level.
    type(atl_level_boundary_type), allocatable :: boundary_list(:)
    !> The boundary description for the faces on the current level.
    type(atl_level_boundary_type), allocatable :: boundary_stab_list(:)
    !> The equation type.
    type(atl_equations_type) :: equation
    !> The tree of the simulation domain
    type(treelmesh_type), intent(in) :: tree
    !> current time
    type(tem_time_type), intent(in) :: time
    !> mpi communication pattern type
    type(tem_commPattern_type), intent(in) :: commPattern
    !> Flag to indicate whether maximal polynomial deviations should be
    !! computed for the state.
    logical, intent(in) :: need_deviation
    ! --------------------------------------------------------------------------
    integer :: iList, iStab, nTotal
    integer :: nScalars, nDer, nDof, nScalarsState_Comm, nScalarsFlux_Comm
    logical :: stab_reqNeigh
    integer :: nBndStabElems(tree%global%minLevel:tree%global%maxLevel,1:3)
    ! --------------------------------------------------------------------------

    ! the number of scalar variables of our equation system
    nScalars = equation%varSys%nScalars

    ! If the stabilization requires further neighbor information,
    ! we init the buffer for it
    stab_reqNeigh = .false.
    do iStab = 1, size(scheme_list(tree%global%minLevel)%stabilization)
      if( scheme_list(tree%global%minLevel)%stabilization(iStab)%reqNeigh ) then
        stab_reqNeigh = .true.
      end if
    end do

    ! init the kerneldata on each level.
    do iList = tree%global%minLevel, tree%global%maxLevel
      ! The number of derived quantities is the number of polynomial degrees of
      ! freedoms per scalar variable.

      select case(scheme_list(iList)%modg_1d%basisType)
      case(Q_space)
        nDer = (scheme_list(iList)%modg_1d%maxPolyDegree+1)**1
        nDof = nDer
      end select

      call atl_init_kerneldata(                                             &
        &    kerneldata     = kerneldata_list(iList),                       &
        &    statedata      = statedata_list(iList),                        &
        &    nTotal         = mesh_list(iList)%descriptor%elem              &
        &                                                %nElems(eT_fluid), &
        &    nVars          = nScalars,                                     &
        &    nDofs          = nDof,                                         &
        &    nDervQuant     = nDer,                                         &
        &    time           = time,                                         &
        &    maxPolyDegree  = scheme_list(iList)%modg_1d%maxPolyDegree,     &
        &    nDims          = 1,                                            &
        &    poly_space     = scheme_list(iList)%modg_1d%basisType,         &
        &    need_deviation = need_deviation,                               &
        &    need_maxgrad   = equation%requires_gradmax                     )

      if(stab_reqNeigh) then
        nBndStabElems = atl_get_numBndElems(      &
          & minLevel      = tree%global%minLevel, &
          & maxLevel      = tree%global%maxLevel, &
          & boundary_list = boundary_stab_list    )
        nTotal = mesh_list(iList)%faces%dimByDimDesc(1)%nElems &
          &      + nBndStabElems(iList,1)
        call atl_init_statedata(                              &
          &    statedata      = statedata_stab_list(iList,1), &
          &    nTotal         = nTotal,                       &
          &    nVars          = nScalars,                     &
          &    nDofs          = nDof,                         &
          &    time           = time                          )
      end if
    end do

    ! Init the parallel module here, as well. DG requires only face values, so
    ! we init only the face buffers and do not init the cell buffers.
    ! ... for the state, we need the state and all its derivatives
    !     (in all directions)
    nScalarsState_Comm = nScalars*(1+equation%nDerivatives*1)
    ! ... for the flux, we need the numerical flux and the stabilization flux
    nScalarsFlux_Comm = nScalars*(1+equation%nDerivatives)
    call atl_init_parallel_module(                   &
      & commPattern          = commPattern,          &
      & scheme               = scheme_list,          &
      & nValsElem            = nScalars,             &
      & nValsStateFace       = nScalarsState_Comm,   &
      & nValsFluxFace        = nScalarsFlux_Comm,    &
      & cube                 = mesh_list,            &
      & boundary             = boundary_list,        &
      & createCellBuffer     = .false.,              &
      & createFaceBuffer     = .true.,               &
      & createStabFaceBuffer = .false.,              &
      & createStabElemBuffer = stab_reqNeigh,        &
      & nBndStabElems        = nBndStabElems,        &
      & minLevel             = tree%global%minLevel, &
      & maxLevel             = tree%global%maxLevel  )

  end subroutine atl_init_modg_1d_kernel



  !> Subroutine to project modal representations of physical flux, numerical
  !! flux and source terms onto test functions.
  subroutine atl_modg_1d_project_testFunc( mesh, equation, kerneldata,        &
    &                                      facedata, sourcedata,              &
    &                                      penalizationdata, usePenalization, &
    &                                      level, scheme                      )
    ! --------------------------------------------------------------------------
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The representation on the face + representation of the flux.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Volumetric data for the penalization
    type(atl_penalizationData_type), intent(in) :: penalizationdata
    !> Flag to indicate whether the penalization data has to be considered, or
    !! if it is taken care of somewhere else (imex)
    logical, intent(in) :: usePenalization
    !> Modal representation of the sources.
    type(atl_source_type), intent(inout) :: sourcedata
    !> The level you compute for
    integer, intent(in) :: level
    !> The parameters of the MODG scheme
    type(atl_modg_1d_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------

    select case(scheme%basisType)
    case(Q_space)
      call modg_1d_project_testFunc_Q( mesh             = mesh,              &
        &                              equation          = equation,         &
        &                              kerneldata        = kerneldata,       &
        &                              facedata          = facedata,         &
        &                              penalizationdata  = penalizationdata, &
        &                              usePenalization   = usePenalization,  &
        &                              sourcedata        = sourcedata,       &
        &                              Currentlevel      = level,            &
        &                              scheme            = scheme            )
    end select

  end subroutine atl_modg_1d_project_testFunc

  !> Subroutine to project modal representations of physical flux, numerical flux
  !! and source terms onto test functions.
  subroutine atl_modg_1d_project_PhysFlux_testFunc( equation, kerneldata, &
    &                                               scheme, iElem, nDofs, &
    &                                               state_data            )
    ! --------------------------------------------------------------------------
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The parameters of the MODG scheme
    type(atl_modg_1d_scheme_type), intent(in) :: scheme
    !> The total degrees of freedom
    integer, intent(in) :: nDofs
    !> The element index
    integer, intent(in) :: iElem
    !> The physical fluxes that needs to be projected
    real(kind=rk), intent(in)  :: state_data(nDofs,equation%varSys%nScalars)
    ! --------------------------------------------------------------------------

    select case(scheme%basisType)
    case(Q_space)
      ! Projection of the physical flux
      ! ... x direction
      call modg_1d_project_physFluxX_Q(             &
        & nScalars      = equation%varSys%nScalars, &
        & maxPolyDegree = scheme%maxPolyDegree,     &
        & iElem         = iElem,                    &
        & state_der     = state_data,               &
        & state         = kerneldata%state_der      )
    end select

  end subroutine atl_modg_1d_project_PhysFlux_testFunc


  !> Subroutine to project modal representations of physical flux, numerical flux
  !! and source terms onto test functions.
  subroutine modg_1d_project_testFunc_Q( mesh, equation, kerneldata, facedata, &
                            & penalizationdata, usePenalization, sourcedata,   &
                            & Currentlevel, scheme                             )
    ! --------------------------------------------------------------------------
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The representation on the face + representation of the flux.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Volumetric data for the penalization
    type(atl_penalizationData_type), intent(in) :: penalizationdata
    !> Flag to indicate whether the penalization data has to be considered, or
    !! if it is taken care of somewhere else (imex)
    logical, intent(in) :: usePenalization
    !> Modal representation of the sources.
    type(atl_source_type), intent(inout) :: sourcedata
    !> The level you compute for
    integer, intent(in) :: Currentlevel
    !> The parameters of the MODG scheme
    type(atl_modg_1d_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------

    ! Iterate over all the elements and do the following:
    ! 1. Project physical fluxes (3 directions) to test functions (stiffness terms).
    !    Attention: physical x flux: state_der(:,:,2,:)
    ! 2. Project numerical fluxes (4 faces) to test functions.
    ! 3. Project source terms to test functions.
    ! Attention: write the projections to state_der(:,:,1,:) because inverse
    ! of mass matrix will be applied to these entries.

    ! Before, we start the projections, we initialize the entries (where we sum up) with
    ! zeros.


    ! Projection of the numerical flux
    ! ... x faces
    call modg_1d_project_numFluxX_Q(                              &
      & numFluxLeftFace  = facedata%faceFlux(1)%dat,              &
      & numFluxRightFace = facedata%faceFlux(1)%dat,              &
      & nScalars         = equation%varSys%nScalars,              &
      & maxPolyDegree    = scheme%maxPolyDegree,                  &
      & nElems_fluid     = mesh%descriptor%elem%nElems(eT_fluid), &
      & projection       = kerneldata%state_der                   )

    select case(equation%eq_kind)
      case('heat_1d')
        ! Projection of the appropriate  numerical flux onto the derivative of
        ! the test functions
        ! ... x faces
        call modg_1d_project_numFlux_diffTestX_Q(                     &
          & numFluxLeftFace  = facedata%faceFlux(1)%dat,              &
          & numFluxRightFace = facedata%faceFlux(1)%dat,              &
          & nScalars         = equation%varSys%nScalars,              &
          & maxPolyDegree    = scheme%maxPolyDegree,                  &
          & nElems_fluid     = mesh%descriptor%elem%nElems(eT_fluid), &
          & projection       = kerneldata%state_der                   )

    end select

    ! Projection of the source terms onto the test functions.
    call modg_1d_project_source_Q( sourcedata = sourcedata,            &
      &                         nScalars = equation%varSys%nScalars,   &
      &                         mesh = mesh,                           &
      &                         maxPolyDegree = scheme%maxPolyDegree,  &
      &                         kerneldata = kerneldata,               &
      &                         currentLevel = currentLevel            )

    ! Project the penalization term to test function
    ! Projection of the penalization term if is not computed somewhere else
    if (penalizationdata%isActive .and. usePenalization) then
      call modg_1d_project_penalization(                               &
        &                 nScalars         = equation%varSys%nScalars, &
        &                 mesh             = mesh,                     &
        &                 maxPolyDegree    = scheme%maxPolyDegree,     &
        &                 penalizationdata = penalizationdata,         &
        &                 kerneldata       = kerneldata                )
    end if


  end subroutine modg_1d_project_testFunc_Q


  !> Projection of the source terms (in modal representation) to the test functions.
  subroutine modg_1d_project_source_Q( sourcedata, nScalars, mesh, &
    &                                  maxPolyDegree, kerneldata,  &
    &                                  currentLevel )
    ! ---------------------------------------------------------------------------
    !> The number scalar variables in the equation system.
    integer, intent(in) :: nScalars
    !> The modal representation of the source
    type(atl_source_type), intent(in) :: sourcedata
    !> The maximal polynomial degree of the modg scheme
    integer, intent(in) :: maxPolyDegree
    !> The cubical elements.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The data of the kernel. This one is updated by the projection.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The current level
    integer, intent(in) :: currentLevel
    ! ---------------------------------------------------------------------------
    integer(kind=long_k) :: elemPos
    integer :: iElem, xTestFunc, testPos, xAnsFuncMin, xAnsFunc, &
      &        ansPos, varPos, iSource, nSourceElems
    real(kind=rk) :: jacobiDet, xScalProd
    ! ---------------------------------------------------------------------------

    jacobiDet = (0.5_rk*mesh%length)**1

    do iSource = 1, size(sourcedata%method)

      nSourceElems = sourcedata%method(iSource)%elems(currentLevel)%nElems

      do varPos = 1, nScalars
        do iElem = 1, nSourceElems

          ! Position of the current element
          if (sourcedata%method(iSource)%isPermanent) then
            elempos = iElem
          else
            elempos = sourcedata%method(iSource)%elems(currentLevel) &
                                %posInTotal%val(iElem)
          end if

          ! Now, we loop over all the test functions for this element and calculate
          ! the projection of the source terms onto this test functions.
          do xTestFunc = 1, maxPolyDegree+1

            ! Get the position of the test functions in the serialized list
            ! of dofs
  testpos = xtestfunc                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

            ! Loop over relevant ans functions
            xAnsFuncMin = xTestFunc-2
            if( xAnsFuncMin < 1 ) then
              xAnsFuncMin = xTestFunc
            end if
            do xAnsFunc = xAnsFuncMin,xTestFunc,2

              ! get position of ansatz functions in the serialized list
              ! of dofs.
  anspos = xansfunc                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

              ! project the current ansatz function onto the test function
              xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
              kernelData%state_der(elemPos, testPos,  varPos) &
                  & = kernelData%state_der(elemPos, testPos,  varPos) &
                  & + xScalProd &
                  & * jacobiDet &
                  & * sourcedata%method(iSource)%val(iElem, ansPos, varPos)
            end do ! x ansatz functions

          end do ! x test functions

        end do ! elem loop
      end do ! variable loop
    end do ! source loop


  end subroutine modg_1d_project_source_Q


  !> Projection of the penalization terms (in modal representation) to the test functions.
  subroutine modg_1d_project_penalization( nScalars, mesh, maxPolyDegree, &
                                           & kerneldata, penalizationdata )
    ! ---------------------------------------------------------------------------
    !> The number scalar variables in the equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree of the modg scheme
    integer, intent(in) :: maxPolyDegree
    !> The cubical elements.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The data of the kernel. This one is updated by the projection.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> Volumetric data for the penalization
    type(atl_penalizationData_type), intent(in) :: penalizationdata
    ! ---------------------------------------------------------------------------
    integer :: iElem, xTestFunc, testPos, &
             & xAnsFuncMin, xAnsFunc,     &
             & ansPos
    real(kind=rk) :: jacobiDet, xScalProd
    ! ---------------------------------------------------------------------------

    jacobiDet = (0.5_rk*mesh%length)

    do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

      ! Now, we loop over all the test functions for this element and calculate
      ! the projection of the term onto this test functions.
      do testpos=1, maxPolyDegree+1

        ! Get the position of the test functions in the serialized list
        ! of dofs and store it in xTestFunc
  xtestfunc = testpos                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        ! Loop over relevant ans functions
        xAnsFuncMin = xTestFunc-2
        if( xAnsFuncMin < 1 ) then
          xAnsFuncMin = xTestFunc
        end if
        do xAnsFunc = xAnsFuncMin,xTestFunc,2

            ! get position of ansatz functions in the serialized list
            ! of dofs.
  anspos = xansfunc                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

            ! project the current ansatz function onto the test function
            xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
            kernelData%state_der(iElem, testPos,  1:nScalars) &
                         & = kernelData%state_der(iElem, testPos, 1:nScalars) &
                         & + xScalProd &
                         & * jacobiDet &
                         & * penalizationdata%penalization_data(iElem, ansPos,1:nScalars)
        end do ! x ansatz functions

      end do ! test functions

    end do ! elem loop

  end subroutine modg_1d_project_penalization

  !> Projection of the numerical flux in x direction onto the testfunctions.
  subroutine modg_1d_project_numFluxX_Q( numFluxLeftFace, numFluxRightFace, &
    &                                    nScalars, maxPolyDegree,           &
    &                                    nElems_fluid,  projection          )
    ! --------------------------------------------------------------------------
    !> The numerical flux on the left face in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: numFluxLeftFace(:,:,:,:)
    !> The numerical flux on the right face in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: numFluxRightFace(:,:,:,:)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The element index
    integer, intent(in) :: nElems_fluid
    !> The numerical flux projected onto the test functions.
    real(kind=rk), intent(inout) :: projection(:,:,:)
    ! --------------------------------------------------------------------------
    integer :: xTestFunc
    integer :: testPos, ansPos
    real(kind=rk) :: outerNormalLeft, outerNormalRight
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iElem, iVar
    ! --------------------------------------------------------------------------


    ! Jacobi determinant for the projections of the numerical fluxes onto the test functions
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)

    ! Loop over all the test functions and project the numerical flux to them.
    do xTestFunc = 1,min(maxPolyDegree+1,2)

      faceValLeft = ply_faceValLeftBndTest(xTestFunc) * outerNormalLeft
      faceValRight = ply_faceValRightBndTest(xTestFunc) * outerNormalRight
      ! position of the test functions in the kerneldata
  testpos = xtestfunc                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

      ! the position of the modal coefficeint of this ansatz functions
  anspos = 1                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)


      ! buffer the result in kernel data and take care of the outer surface unit normal
      ! ... for the left face
      elementLoop: do iElem = 1,nElems_fluid
        do iVar = 1,nScalars
          projection(iElem,testPos,iVar) = projection(iElem,testPos,iVar) &
                          & - faceValLeft &
                          & * numFluxLeftFace(iElem,ansPos,iVar,1)
          ! ... for the right face
          projection(iElem,testPos,iVar) = projection(iElem,testPos,iVar) &
                          & - faceValRight &
                          & * numFluxRightFace(iElem,ansPos,iVar,2)

        end do

      end do elementLoop
    end do

  end subroutine modg_1d_project_numFluxX_Q


  !> Projection of the numerical flux in x direction onto the differentiated
  !  testfunctions.
  subroutine modg_1d_project_numFlux_diffTestX_Q( &
    &          numFluxLeftFace, numFluxRightFace, nScalars, maxPolyDegree, &
    &          nElems_fluid,  projection )
    ! --------------------------------------------------------------------------
    !> The numerical flux on the left face in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: numFluxLeftFace(:,:,:,:)
    !> The numerical flux on the right face in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: numFluxRightFace(:,:,:,:)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The element index
    integer, intent(in) :: nElems_fluid
    !> The numerical flux projected onto the test functions.
    real(kind=rk), intent(inout) :: projection(:,:,:)
    ! --------------------------------------------------------------------------
    integer :: xTestFunc
    integer :: testPos, ansPos
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iElem, iVar
    ! --------------------------------------------------------------------------


    ! Jacobi determinant for the projections of the numerical fluxes
    ! onto the test functions

    ! Loop over all the test functions and project the numerical flux to them.
    do xTestFunc = 1,maxPolyDegree+1

      faceValLeft = ply_faceValLeftBndgradTest(xTestFunc)
      faceValRight = ply_faceValRightBndgradTest(xTestFunc)

      ! position of the test functions in the kerneldata
  testpos = xtestfunc                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

      ! the position of the modal coefficeint of this ansatz functions
  anspos = 1                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)


      ! buffer the result in kernel data and take care of the outer surface unit normal
      ! ... for the left face
      elementLoop: do iElem = 1,nElems_fluid
        do iVar = 1,nScalars
          projection(iElem,testPos,iVar) = projection(iElem,testPos,iVar) &
                          & - faceValLeft &
                          & * numFluxLeftFace(iElem,ansPos,nScalars+iVar,1)
          ! ... for the right face
          projection(iElem,testPos,iVar) = projection(iElem,testPos,iVar) &
                          & - faceValRight &
                          & * numFluxRightFace(iElem,ansPos,nscalars+iVar,2)
        end do
      end do elementLoop

    end do
  end subroutine modg_1d_project_numFlux_diffTestX_Q


  !> Projection of the physical flux in x direction onto the testfunctions.
  subroutine modg_1d_project_physFluxX_Q( nScalars, maxPolyDegree, state, &
    &                                     iElem, state_der )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> The element index
    integer, intent(in) :: iElem
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der((maxPolyDegree+1),nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: scalProdX(maxPolyDegree)
    integer :: ansPos(4), testPos
    integer :: iAnsX
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(4)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test
    ! functions.
    ! This is the stiffness term!

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    ! for x direction (x test function differentiated)
    ! get the relevant indices for the ansatz function
    do iAnsX=1,maxPolyDegree
      scalProdX(iAnsX) = ply_scalProdDualLegDiff(iAnsX, iAnsX+1)
    end do

    ! now, project onto all test functions
  testpos = 1                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
    do iVar=var_lb,var_ub
      state(iElem,testPos,iVar) = 0.0_rk
    end do

    ! Need to add one term
    do iAnsX = 1, maxPolyDegree
      ! the position of the current test functions
  testpos = iansx+1                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
      scalProd(1) = scalProdX(iAnsX)
  anspos(1) = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

      do iVar=var_lb,var_ub
        state(iElem,testPos,iVar) &
          &  = state_der(ansPos(1),iVar) * scalProd(1)
      end do

    end do

  end subroutine modg_1d_project_physFluxX_Q



  !> Lift the value on the face to a positive value if necessary.
  !!
  !! In some equations, some variables may not be allowed to be negative,
  !! for example the density and energy in flow simulations.
  !! With this routine we enforce this limitation for the projected mean
  !! state on the element surfaces.
  subroutine atl_modg_1d_ensure_pos_face(                                   &
    &          nelems_fluid, volState, faceRep, nScalars, ensure_positivity )
    ! -------------------------------------------------------------------- !
    !> Number of fluid elements
    integer, intent(in) :: nElems_fluid

    !> Volumetric modal states for each element
    real(kind=rk), intent(in) :: volState(:,:,:)

    !> Modal representation on the face that might need limiting.
    !!
    !! All means below 0 of variables where positivity is to be ensured will
    !! be increased to a small value above 0.
    type(atl_faceRep_type), intent(inout) :: faceRep(:)

    !> Number of variables
    integer, intent(in) :: nScalars

    !> Indication for each variable, whether positivity is to be ensured
    !! or not.
    !!
    !! Only variables where this is set to .true. will be modified by this
    !! routine.
    logical, intent(in) :: ensure_positivity(:)
    ! --------------------------------------------------------------------------
    integer :: iVar
    real(kind=rk) :: epsfact
    ! --------------------------------------------------------------------------

    epsfact = epsilon(volState(1,1,1))

    do iVar=1,nScalars
      if (ensure_positivity(iVar)) then

        ! We assume the integral mean of the element is positive for all
        ! elements. The integral mean on the surface will now be lifted to
        ! a small fraction of that element mean if it is below that threshold.
        facerep(1)%dat(1:nElems_fluid,1,iVar,1)             &
          & = max( facerep(1)%dat(1:nElems_fluid,1,iVar,1), &
          &        epsfact * volState(1:nElems_fluid,1,iVar)   )
        facerep(1)%dat(1:nElems_fluid,1,iVar,2)             &
          & = max( facerep(1)%dat(1:nElems_fluid,1,iVar,2), &
          &        epsfact * volState(1:nElems_fluid,1,iVar)   )

      end if
    end do

  end subroutine atl_modg_1d_ensure_pos_face


  !> Projects modal representation of each cell to its faces, i.e.
  !! this subroutine creates a modal representation on the faces.
  subroutine atl_modg_1d_modalVolToModalFace( mesh , statedata, facedata, &
    &                                         nScalars, modg_1d, equation )
    ! --------------------------------------------------------------------------
    !> The elements we apply the projection for.
    type(atl_cube_elem_type),  intent(in) :: mesh
    !> Volumetric, modal states for each element.
    type(atl_statedata_type), intent(in) :: statedata
    !> Modal representation on the face (will be updated by this routine for all
    !! fluid elements in mesh).
    type(atl_facedata_type), intent(inout) :: facedata
    !> The number of scalars varaibales in your equation system.
    integer, intent(in) :: nScalars
    !> The parameters of your modg scheme.
    type(atl_modg_1d_scheme_type), intent(in) :: modg_1d
    !> The equation you solve.
    type(atl_equations_type) :: equation
    ! --------------------------------------------------------------------------
    integer :: spaceDir, iDir
    ! --------------------------------------------------------------------------

    ! Iterate over all the fluid elements and project its modal representations
    ! to the faces of the element.

    select case(modg_1d%basisType)
      case(Q_space) ! Q tensor product ansatz functions

        do iDir = 1, size(facedata%faceRep)
          facedata%faceRep(iDir)                              &
            &     %dat(1:mesh%descriptor                      &
            &                %elem                            &
            &                %nElems(eT_fluid),:,:,:) = 0.0_rk
        end do


        ! Now, iterate over all the faces and project to this face.
        ! ... x faces
        ! Project to the face in the current direction.
        spaceDir = qAxis(1)
        call modg_1d_volToFace_Q(                                  &
          & volState      = statedata%state,                       &
          & maxPolyDegree = modg_1d%maxPolyDegree,                 &
          & faceDir       = 1,                                     &
          & nScalars      = nScalars,                              &
          & nElems        = mesh%descriptor%elem%nElems(eT_fluid), &
          & faceState     = facedata%faceRep(spaceDir)%dat         )
        ! ... x
        ! Project to the face in the current direction.
        spaceDir = qAxis(4)
        call modg_1d_volToFace_Q(                                  &
          & volState      = statedata%state,                       &
          & maxPolyDegree = modg_1d%maxPolyDegree,                 &
          & faceDir       = 4,                                     &
          & nScalars      = nScalars,                              &
          & nElems        = mesh%descriptor%elem%nElems(eT_fluid), &
          & faceState     = facedata%faceRep(spaceDir)%dat         )


        ! Now, iterate over all the faces and project the gradients to this face.
        if(equation%nDerivatives == 1 ) then
          spaceDir = qAxis(1)

          call modg_1d_VolToFace_grad_Q(                             &
            & volState      = statedata%state,                       &
            & maxPolyDegree = modg_1d%maxPolyDegree,                 &
            & faceDir       = 1,                                     &
            & nScalars      = nScalars,                              &
            & nElems        = mesh%descriptor%elem%nElems(eT_fluid), &
            & elemLength    = mesh%length,                           &
            & faceState     = facedata%faceRep(spaceDir)%dat         )

          spaceDir = qAxis(4)

          call modg_1d_VolToFace_grad_Q(                             &
            & volState      = statedata%state,                       &
            & maxPolyDegree = modg_1d%maxPolyDegree,                 &
            & faceDir       = 4,                                     &
            & nScalars      = nScalars,                              &
            & nElems        = mesh%descriptor%elem%nElems(eT_fluid), &
            & elemLength    = mesh%length,                           &
            & faceState     = facedata%faceRep(spaceDir)%dat         )
          end if
      case default
        call tem_abort( 'ERROR in atl_modg_1d_modalVolToModalFace:' &
          & // ' Unknown tensor product'                            )
    end select
  end subroutine atl_modg_1d_modalVolToModalFace

  subroutine modg_1d_VolToFace_grad_Q( volstate, maxPolyDegree, faceDir,       &
    &                         nScalars, nElems,elemLength, faceState           )
    ! --------------------------------------------------------------------------
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the number
    !! of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink, e.g.
    !! \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The number of elements
    integer, intent(in) :: nElems
    !> Length of elements
    real(kind=rk), intent(in) :: elemLength
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: iElem, lb, ub
    integer :: pos, iAnsX,leftOrRight
    real(kind=rk) :: faceVal
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    real(kind=rk), allocatable :: state_gradient(:,:)
    ! --------------------------------------------------------------------------

    ! Allocate memory for the gradient buffer
    allocate( state_gradient(maxPolyDegree+1,nScalars) )
    allocate( modalCoeffs(maxPolyDegree+1,nScalars) )

    lb = nScalars + 1
    ub = lb + nScalars -1

    leftOrRight = tem_dirToFace_map(faceDir)

    do iElem = 1, nElems

      modalCoeffs(:,:) = volState(iElem,:,:)

      ! Now we calculate the gradient of the modal representation required
      call ply_calcDiff_leg_1d( legCoeffs     = modalCoeffs,    &
        &                       legcoeffsDiff = state_gradient, &
        &                       maxPolyDegree = maxPolyDegree,  &
        &                       elemLength    = elemLength      )

        ! now, project to a certain face.
      if( faceDir.eq.q__E ) then
        do iAnsX = 1, maxPolyDegree+1
          ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          faceState(iElem,1,lb:ub,leftOrRight) &
            &  = faceState(iElem,1,lb:ub,leftOrRight) &
            &  + state_gradient(pos,1:nScalars)
        end do
      elseif(faceDir.eq.q__W) then
        do iAnsX = 1, maxPolyDegree+1
          ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          ! get the face value of the ansatz function with fixed coordinate
          faceVal = ply_faceValLeftBndAns(iAnsX)
          faceState(iElem,1,lb:ub,leftOrRight) &
            &  = faceState(iElem,1,lb:ub,leftOrRight) &
            &  + faceVal*state_gradient(pos,1:nScalars)
        end do
      else
        write(logUnit(1),*) 'ERROR in modg_1d_volToFace: unsupported       &
          &                                 face direction, stopping...'
        call tem_abort()
      end if

    end do  ! elems loop

  end subroutine modg_1d_VolToFace_grad_Q

  !> Projects derivative of modal representation of each cell to its faces, i.e.
  !! this subroutine creates a modal representation on the faces.
!NA!  subroutine modg_1d_modalVolToModalFace_der( mesh , statedata, facedata,      &
!NA!    &                                          nScalars, modg_1d, poly_proj    )
!NA!    ! --------------------------------------------------------------------------
!NA!    !> The elements we apply the projection for.
!NA!    type(atl_cube_elem_type),  intent(in) :: mesh
!NA!    !> Volumetric, modal states for each element.
!NA!    type(atl_statedata_type), intent(in) :: statedata
!NA!    !> Modal representation on the face (will be updated by this routine for all
!NA!    !! fluid elements in mesh).
!NA!    type(atl_facedata_type), intent(inout) :: facedata
!NA!    !> The number of scalars varaibales in your equation system.
!NA!    integer, intent(in) :: nScalars
!NA!    !> The parameters of your modg scheme.
!NA!    type(atl_modg_1d_scheme_type), intent(in) :: modg_1d
!NA!    !> Parameters for projection
!NA!    type(ply_poly_project_type), intent(inout) :: poly_proj
!NA!    ! --------------------------------------------------------------------------
!NA!    integer :: spaceDir, leftOrRightFace, iDir
!NA!    integer :: dofs , nquadpointsPerDir, iDegX, iElem, iDof
!NA!    real(kind=rk), allocatable :: modalCoeffs(:,:)
!NA!    real(kind=rk), allocatable :: modalCoeffsDiff(:,:)
!NA!    real(kind=rk), allocatable :: state_der(:,:,:)
!NA!    ! --------------------------------------------------------------------------
!NA!
!NA!    dofs = poly_proj%body_1d%ndofs
!NA!    nquadpointsPerDir = poly_proj%nquadpointsPerDir
!NA!
!NA!    allocate( modalCoeffs(dofs,1) )
!NA!    allocate( modalCoeffsDiff(dofs,1) )
!NA!    allocate( state_der(mesh%descriptor%nElems_fluid,dofs,1) )
!NA!
!NA!    elemLoop: do iElem = 1, mesh%descriptor%nElems_fluid
!NA!
!NA!      ! get the modal coefficients of the current cell (for all variables)
!NA!      ! -->modal space
!NA!      modalCoeffs(:,:) = 0.0_rk
!NA!      do iDegX = 1, dofs
!NA!        idof = 1 + (iDegX-1)
!NA!        modalCoeffs(idof,:) = statedata%state(iElem,idof,:)
!NA!      end do
!NA!
!NA!      ! Now we calculate the gradient of the modal representation required
!NA!      call ply_calcDiff_leg_1d( legCoeffs     = modalCoeffs,             &
!NA!        &                       legcoeffsDiff = modalCoeffsDiff,         &
!NA!        &                       maxPolyDegree = poly_proj%maxPolyDegree, &
!NA!        &                       nVars         = 1                        )
!NA!
!NA!        state_der(iElem,:,1) = modalCoeffsDiff(:,1)
!NA!    end do elemLoop
!NA!
!NA!    ! Iterate over all the fluid elements and project its modal representations
!NA!    ! of the state derivatives to the faces of the element.
!NA!    select case(modg_1d%basisType)
!NA!      case(Q_space) ! Q tensor product ansatz functions
!NA!
!NA!! @todo commented out!
!NA!write(*,*) 'Commented out der_dat in 1D modg kernel, ... STOPPING'
!NA!stop
!NA!!JZ:        do iDir = 1, size(facedata%faceRep)
!NA!!JZ:          facedata%faceRep(iDir)%der_dat(1:mesh%descriptor%nElems_fluid,:,:,:,1) &
!NA!!JZ:            &                                                    = 0.0_rk
!NA!!JZ:        end do
!NA!
!NA!
!NA!        ! Now, iterate over all the faces and project to this face.
!NA!        ! ... x faces
!NA!        ! Project to the face in the current direction.
!NA!        spaceDir = qAxis(1)
!NA!! @todo commented out!
!NA!write(*,*) 'Commented out der_dat in 1D modg kernel, ... STOPPING'
!NA!stop
!NA!!JZ:        call modg_1d_volToFace_Q( volState = state_der,                       &
!NA!!JZ:              & maxPolyDegree = modg_1d%maxPolyDegree,         &
!NA!!JZ:              & faceDir = 1,                                   &
!NA!!JZ:              & nScalars = nScalars,                           &
!NA!!JZ:              & nElems = mesh%descriptor%nElems_fluid,         &
!NA!!JZ:              & faceState = facedata%faceRep(spaceDir)%der_dat(:,:,:,:,1) )
!NA!        ! ... x faces
!NA!        ! Project to the face in the current direction.
!NA!        spaceDir = qAxis(4)
!NA!! @todo commented out!
!NA!write(*,*) 'Commented out der_dat in 1D modg kernel, ... STOPPING'
!NA!stop
!NA!!JZ:        call modg_1d_volToFace_Q( volState = state_der,                       &
!NA!!JZ:              & maxPolyDegree = modg_1d%maxPolyDegree,         &
!NA!!JZ:              & faceDir = 4,                                   &
!NA!!JZ:              & nScalars = nScalars,                           &
!NA!!JZ:              & nElems = mesh%descriptor%nElems_fluid,         &
!NA!!JZ:              & faceState = facedata%faceRep(spaceDir)%der_dat(:,:,:,:,1) )
!NA!
!NA!
!NA!      case default
!NA!        write(logUnit(1),*) 'ERROR in modg_1d_modalVolToModalFace_der:         &
!NA!          &                      Unknown tensor product, stopping ...'
!NA!        call tem_abort()
!NA!    end select
!NA!
!NA!  end subroutine modg_1d_modalVolToModalFace_der
!NA!


  !> Project modal representation of an element to one of its faces for Q space.
  !!
  !! Project modal representation of an element onto one of its faces. Therefore,
  !! this function returns the modal representation of the solution on the face. This
  !! function can project onto an arbitrary face direction.
  subroutine modg_1d_volToFace_Q( volState, maxPolyDegree, faceDir, nScalars, nELems, faceState )
    ! --------------------------------------------------------------------------
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the number
    !! of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink, e.g.
    !! \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iVar
    integer :: pos, iAnsX,leftOrRight
    real(kind=rk) :: faceVal
    ! --------------------------------------------------------------------------

    leftOrRight = tem_dirToFace_map(faceDir)

    ! now, project to a certain face.
    if( faceDir.eq.q__E ) then
      do iAnsX = 1, maxPolyDegree+1
        ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        do iElem=1,nElems
          do iVar = 1, nScalars
            faceState(iElem,1,iVar,leftOrRight) &
              &  = faceState(iElem,1,iVar,leftOrRight) &
              &  + volState(iElem,pos,iVar)
          end do
        end do
      end do
    elseif(faceDir.eq.q__W) then
      do iAnsX = 1, maxPolyDegree+1
        ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        ! get the face value of the ansatz function with fixed coordinate
        faceVal = ply_faceValLeftBndAns(iAnsX)
        do iElem=1,nElems
          do iVar = 1, nScalars
            faceState(iElem,1,iVar,leftOrRight) &
              &  = faceState(iElem,1,iVar,leftOrRight) &
              &  + faceVal*volState(iElem,pos,iVar)
          end do
        end do
      end do
    else
      write(logUnit(1),*) 'ERROR in modg_1d_volToFace: unsupported face direction, stopping...'
      call tem_abort()
    end if

  end subroutine modg_1d_volToFace_Q


  !> Applies the inverse of the mass matrix for a 3D scheme.
  subroutine atl_modg_1d_invMassMatrix( mesh, kerneldata, statedata,        &
    &                                   elementalTimestep, timestep, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The data of the kernel.
    type(atl_kerneldata_type) :: kerneldata
    !> THe state if the equation
    type(atl_statedata_type) :: statedata
    !> The elemental timestepping routine, because of performance, this
    !! is constant.
    procedure(atl_elemental_timestep_vec), pointer, intent(in) &
      & :: elementalTimestep
    !> The timestepping data.
    type(atl_timestep_type) :: timestep
    !> Parameters of the modal dg scheme
    type(atl_modg_1d_scheme_type), intent(in) :: scheme
    ! -------------------------------------------------------------------------!
    integer :: nElems

    nElems = mesh%descriptor%elem%nElems(eT_fluid)

    select case(scheme%basisType)
      case(Q_space)
        call modg_1d_invMassMatrix_Q( mesh, kerneldata, scheme, nElems )
    end select

    ! Since this is the last part of the MOdal Discontinuous Galerkin (MODG)
    ! algorithm we can make the final timestep, now.
    call elementalTimestep( timestep, statedata%state, kerneldata)

  end subroutine atl_modg_1d_invMassMatrix

  !> Applies the inverse of the mass matrix for a 3D scheme.
  subroutine modg_1d_invMassMatrix_Q( mesh, kerneldata, scheme, nElems )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The data of the kernel.
    type(atl_kerneldata_type) :: kerneldata
    !> Parameters of the modal dg scheme
    type(atl_modg_1d_scheme_type), intent(in) :: scheme
    integer, intent(in) :: nElems
    ! -------------------------------------------------------------------------!
    ! Loop indices for ansatz functions
    integer :: iAnsX
    ! Positions for the given ansatz functions
    integer :: ansLow, ansPrevLow
    ! Inverse of the determintant of the jacobian of the mapping from reference
    ! to physical element.
    real(kind=rk) :: inv_jacobiDet
    ! -------------------------------------------------------------------------!

    inv_jacobiDet = (2.0_rk/mesh%length)

    ! apply the 1D inverse of the mass matrix
    do iAnsX = 3, scheme%maxPolyDegree+1
  anslow = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansprevlow = iansx-2                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
      kerneldata%state_der(:nElems, ansLow,  :) &
        &        = kerneldata%state_der(:nElems, ansLow,  :) &
        &        + kerneldata%state_der(:nElems, ansPrevLow,  :)
    end do


    do iAnsX = 1, scheme%maxPolyDegree+1
  anslow = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
      kerneldata%state_der(:nElems, ansLow,  :) &
        &    = kerneldata%state_der(:nElems, ansLow,  :) &
        &      * 0.5_rk*(2*iAnsX-1)*inv_jacobiDet
    end do

  end subroutine modg_1d_invMassMatrix_Q


end module atl_modg_1d_kernel_module
