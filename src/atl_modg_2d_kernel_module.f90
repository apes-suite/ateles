! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2015, 2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2019 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
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
module atl_modg_2d_kernel_module
  use env_module,                       only: rk
  ! treelm
  use tem_aux_module,                   only: tem_abort
  use treelmesh_module,                 only: treelmesh_type
  use tem_element_module,               only: eT_fluid
  use tem_param_module,                 only: q__E, &
    &                                         q__W, &
    &                                         q__N, &
    &                                         q__S, &
    &                                         qAxis
  use tem_time_module,                  only: tem_time_type
  use tem_faceData_module,              only: tem_dirToFace_map, &
    &                                         tem_left,  &
    &                                         tem_right
  use tem_logging_module,               only: logUnit
  use tem_comm_module,                  only: tem_commPattern_type
  use tem_comm_env_module,              only: tem_comm_env_type
  use tem_element_module,               only: eT_fluid
  ! polynomials
  use ply_modg_basis_module,            only: ply_faceValLeftBndTest,      &
    &                                         ply_faceValRightBndTest,     &
    &                                         ply_scalProdDualLeg,         &
    &                                         ply_scalProdDualLegDiff,     &
    &                                         ply_faceValLeftBndTestGrad,  &
    &                                         ply_faceValRightBndTestGrad, &
    &                                         ply_faceValLeftBndgradTest,  &
    &                                         ply_faceValRightBndgradTest
  use ply_dof_module,                   only: ply_change_poly_space, &
    &                                         Q_space,               &
    &                                         P_space
  use ply_poly_project_module,          only: ply_poly_project_type, &
    &                                         assignment(=),         &
    &                                         ply_poly_project_m2n,  &
    &                                         ply_poly_project_n2m
  ! ateles
  use atl_modg_2d_scheme_module,        only: atl_modg_2d_scheme_type
  use atl_equation_module,              only: atl_equations_type
  use atl_kerneldata_module,            only: atl_kerneldata_type, &
    &                                         atl_init_kerneldata, &
    &                                         atl_init_statedata,  &
    &                                         atl_statedata_type
  use atl_cube_elem_module,             only: atl_cube_elem_type
  use atl_parallel_module,              only: atl_init_parallel_module
  use atl_scheme_module,                only: atl_scheme_type
  use atl_elemental_time_integration_module,                              &
    &                                   only: atl_elemental_timestep_vec, &
    &                                         atl_timestep_type
  use atl_source_types_module,          only: atl_source_type
  use atl_facedata_module,              only: atl_facedata_type,      &
    &                                         atl_faceRep_type,       &
    &                                         atl_elemfaceToNormal_prp
  use atl_boundary_module,              only: atl_level_boundary_type, &
    &                                         atl_get_numBndElems
  use atl_physFluxNvrstk_2d_module,     only: atl_mult_nu11_NavierStokes_2d, &
    &                                         atl_mult_nu21_NavierStokes_2d, &
    &                                         atl_mult_nu12_NavierStokes_2d, &
    &                                         atl_mult_nu22_NavierStokes_2d
  use atl_physFluxFilNvrStk_module,     only: atl_mult_nu11_Rans_2d, &
    &                                         atl_mult_nu21_Rans_2d, &
    &                                         atl_mult_nu12_Rans_2d, &
    &                                         atl_mult_nu22_Rans_2d
  use atl_materialPrp_module,           only: atl_material_type
  use atl_materialIni_module,           only: atl_update_materialParams
  use atl_penalization_module,          only: atl_penalizationData_type
  use atl_volToFace_module,             only: atl_modg_2d_volToFace_Q,     &
    &                                         atl_modg_2d_volToFace_P,     &
    &                                         atl_modg_2d_volToFace_grad_Q

  implicit none
  private

  public :: atl_init_modg_2d_kernel,    &
    & atl_modg_2d_invMassMatrix,        &
    & atl_preprocess_modg_2d_kernel,    &
    & atl_modg_2d_modalVolToModalFace,  &
    & atl_modg_2d_ensure_pos_facemean,  &
    & atl_modg_2d_project_source,       &
    & atl_modg_2d_project_NumFlux


contains


  !> Initiate the MODG kernel for cubic elements on all levels.
  subroutine atl_init_modg_2d_kernel(                                     &
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
    integer :: iList, iStab, nTotal, iDir
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

      select case(scheme_list(iList)%modg_2d%basisType)
      case(Q_space)
        nDer = (scheme_list(iList)%modg_2d%maxPolyDegree+1)**2
        nDof = nDer
      case(P_space)
        nDer = ( (scheme_list(iList)%modg_2d%maxPolyDegree+1)     &
         &     * (scheme_list(iList)%modg_2d%maxPolyDegree+2) ) / 2
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
        &    maxPolyDegree  = scheme_list(iList)%modg_2d%maxPolyDegree,     &
        &    nDims          = 2,                                            &
        &    poly_space     = scheme_list(iList)%modg_2d%basisType,         &
        &    need_deviation = need_deviation,                               &
        &    need_maxgrad   = equation%requires_gradmax                     )

      if(stab_reqNeigh) then
        nBndStabElems = atl_get_numBndElems(      &
          & minLevel      = tree%global%minLevel, &
          & maxLevel      = tree%global%maxLevel, &
          & boundary_list = boundary_stab_list    )
        do iDir = 1,2
          nTotal = mesh_list(iList)%faces%dimByDimDesc(iDir)%nElems &
            &      + nBndStabElems(iList,iDir)
          call atl_init_statedata(                              &
            & statedata      = statedata_stab_list(iList,iDir), &
            & nTotal         = nTotal,                          &
            & nVars          = nScalars,                        &
            & nDofs          = nDof,                            &
            & time           = time                             )
        end do
      end if
    end do


    ! Init the parallel module here, as well. DG requires only face values, so we
    ! init only the face buffers and do not init the cell buffers.
    ! ... for the state, we need the state and all its derivatives (in all directions)
    nScalarsState_Comm = nScalars*(1+equation%nDerivatives*2)
    ! ... for the flux, we need the numerical flux and the stabilization flux
    nScalarsFlux_Comm = nScalars*(1+equation%nDerivatives)
    call atl_init_parallel_module(                   &
      & commPattern          = commPattern,          &
      & scheme               = scheme_list,          &
      & nValsElem            = nScalars,             &
      & nValsStateFace       = nScalarsState_Comm,   &
      & nValsFluxFace        =  nScalarsFlux_Comm,   &
      & cube                 = mesh_list,            &
      & boundary             = boundary_list,        &
      & createCellBuffer     = .false.,              &
      & createFaceBuffer     = .true.,               &
      & createStabFaceBuffer = .false.,              &
      & createStabElemBuffer = stab_reqNeigh,        &
      & nBndStabElems        = nBndStabElems,        &
      & minLevel             = tree%global%minLevel, &
      & maxLevel             = tree%global%maxLevel  )

  end subroutine atl_init_modg_2d_kernel


  !> Subroutine to execute the preprocessing for the MODG kernels.
  !! Currently this includes: Convert external source terms to modal
  !! representation.
  subroutine atl_preprocess_modg_2d_kernel(                                &
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
    type(tem_comm_env_type) :: proc
    !> mpi communication pattern type
    type(tem_commPattern_type) :: commPattern
    ! --------------------------------------------------------------------------

    call atl_update_materialParams( equation    = equation,                    &
      &                             mesh        = mesh_list(currentLevel),     &
      &                             scheme      = scheme_list(currentLevel),   &
      &                             boundary    = boundary_list(currentLevel), &
      &                             material    = material_list(currentLevel), &
      &                             time        = statedata%local_time,        &
      &                             poly_proj   = poly_proj_material,          &
      &                             proc        = proc,                        &
      &                             commPattern = commPattern                  )


  end subroutine atl_preprocess_modg_2d_kernel


  !> Subroutine to project modal representations of numerical flux
  !! and source terms onto test functions.
  subroutine atl_modg_2d_project_NumFlux( mesh, equation, kerneldata, &
                   & facedata, penalizationdata, usePenalization,     &
                   & scheme, poly_proj, dl_prod, dl_prodDiff          )
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
    !> The parameters of the MODG scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    !> Projection for the current level
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Precomputed scalar products of the test and ansatz function
    real(kind=rk) , intent(in) :: dl_prod(2, scheme%maxPolyDegree+1)
    real(kind=rk) , intent(in) :: dl_prodDiff(2, scheme%maxPolyDegree+1)
    ! --------------------------------------------------------------------------
    integer :: nScalars, nElems_fluid
    real(kind=rk), allocatable ::state_der_Q(:,:,:)
    ! --------------------------------------------------------------------------

    nScalars = equation%varSys%nScalars
    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid)

    ! Projection of the numerical flux
    select case(scheme%basisType)
    case(Q_space)
      ! ... x faces
      call modg_2d_project_numFluxX_Q(              &
        & numFlux       = facedata%faceFlux(1)%dat, &
        & nScalars      = equation%varSys%nScalars, &
        & maxPolyDegree = scheme%maxPolyDegree,     &
        & length        = mesh%length,              &
        & nElems_fluid  = nElems_fluid,             &
        & dl_prod       = dl_prod,                  &
        & projection    = kerneldata%state_der      )
      ! ... y faces
      call modg_2d_project_numFluxY_Q(              &
        & numFlux       = facedata%faceFlux(2)%dat, &
        & nScalars      = equation%varSys%nScalars, &
        & maxPolyDegree = scheme%maxPolyDegree,     &
        & length        = mesh%length,              &
        & nElems_fluid  = nElems_fluid,             &
        & dl_prod       = dl_prod,                  &
        & projection    = kerneldata%state_der      )

    case(P_space)

      ! In order to use the same routines for the numflux projection for
      ! both Q- and P-space we need to change the polynomial space before
      ! the projetion and change back to original polynomial space
      ! after the projection.
      ! Note: faceFlux%dat is a 1D state so the indexing is the same
      !       for Q-space and P-space. No not need to change anything.

      allocate(state_der_Q(nElems_fluid,(scheme%maxPolyDegree+1)**2, &
        &                    nScalars))

      call ply_change_poly_space( inspace    = P_space,              &
        &                         instate    = kerneldata%state_der, &
        &                         outstate   = state_der_Q,          &
        &                         maxPolyDeg = scheme%maxPolyDegree, &
        &                         nElems     = nElems_fluid,         &
        &                         nVars      = nScalars,             &
        &                         nDims      = 2                     )

      ! ... x faces
      call modg_2d_project_numFluxX_Q(              &
        & numFlux       = facedata%faceFlux(1)%dat, &
        & nScalars      = equation%varSys%nScalars, &
        & maxPolyDegree = scheme%maxPolyDegree,     &
        & length        = mesh%length,              &
        & nElems_fluid  = nElems_fluid,             &
        & dl_prod       = dl_prod,                  &
        & projection    = state_der_Q               )
      ! ... y faces
      call modg_2d_project_numFluxY_Q(              &
        & numFlux       = facedata%faceFlux(2)%dat, &
        & nScalars      = equation%varSys%nScalars, &
        & maxPolyDegree = scheme%maxPolyDegree,     &
        & length        = mesh%length,              &
        & nElems_fluid  = nElems_fluid,             &
        & dl_prod       = dl_prod,                  &
        & projection    = state_der_Q               )

      call ply_change_poly_space( inspace    = Q_space,              &
        &                         instate    = state_der_Q,          &
        &                         outstate   = kerneldata%state_der, &
        &                         maxPolyDeg = scheme%maxPolyDegree, &
        &                         nElems     = nElems_fluid,         &
        &                         nVars      = nScalars,             &
        &                         nDims      = 2                     )
      deallocate(state_der_Q)

    end select

    if(equation%nDerivatives == 1) then

      select case(equation%eq_kind)
      case('heat_2d')
        select case(scheme%basisType)
        case(Q_space)
          ! Projection of the appropriate  numerical flux onto the derivative of
          ! the test functions
          ! ... x faces
          call modg_2d_project_numFluxX_diffTestX_Q(        &
            & numFluxLeftFace  = facedata%faceFlux(1)%dat,  &
            & numFluxRightFace = facedata%faceFlux(1)%dat,  &
            & nScalars         = equation%varSys%nScalars,  &
            & maxPolyDegree    = scheme%maxPolyDegree,      &
            & length           = mesh%length,               &
            & nElems_fluid     = nElems_fluid,              &
            & dl_prod          = dl_prodDiff,               &
            & projection       = kerneldata%state_der       )
          ! ... y faces
          call modg_2d_project_numFluxY_diffTestY_Q(       &
            & numFluxLeftFace  = facedata%faceFlux(2)%dat, &
            & numFluxRightFace = facedata%faceFlux(2)%dat, &
            & nScalars         = equation%varSys%nScalars, &
            & maxPolyDegree    = scheme%maxPolyDegree,     &
            & length           = mesh%length,              &
            & nElems_fluid     = nElems_fluid,             &
            & dl_prod          = dl_prodDiff,              &
            & projection       = kerneldata%state_der      )
        case(P_space)
          allocate(state_der_Q(nElems_fluid,(scheme%maxPolyDegree+1)**2, &
            &                    nScalars))

          call ply_change_poly_space( inspace    = P_space,              &
            &                         instate    = kerneldata%state_der, &
            &                         outstate   = state_der_Q,          &
            &                         maxPolyDeg = scheme%maxPolyDegree, &
            &                         nElems     = nElems_fluid,         &
            &                         nVars      = nScalars,             &
            &                         nDims      = 2                     )
          ! Projection of the appropriate  numerical flux onto the derivative of
          ! the test functions
          ! ... x facesa
          call modg_2d_project_numFluxX_diffTestX_Q(        &
            & numFluxLeftFace  = facedata%faceFlux(1)%dat,  &
            & numFluxRightFace = facedata%faceFlux(1)%dat,  &
            & nScalars         = equation%varSys%nScalars,  &
            & maxPolyDegree    = scheme%maxPolyDegree,      &
            & length           = mesh%length,               &
            & nElems_fluid     = nElems_fluid,              &
            & dl_prod          = dl_prodDiff,               &
            & projection       = state_der_Q                )
          ! ... y faces
          call modg_2d_project_numFluxY_diffTestY_Q(       &
            & numFluxLeftFace  = facedata%faceFlux(2)%dat, &
            & numFluxRightFace = facedata%faceFlux(2)%dat, &
            & nScalars         = equation%varSys%nScalars, &
            & maxPolyDegree    = scheme%maxPolyDegree,     &
            & length           = mesh%length,              &
            & nElems_fluid     = nElems_fluid,             &
            & dl_prod          = dl_prodDiff,              &
            & projection       = state_der_Q               )

          call ply_change_poly_space( inspace    = Q_space,              &
            &                         instate    = state_der_Q,          &
            &                         outstate   = kerneldata%state_der, &
            &                         maxPolyDeg = scheme%maxPolyDegree, &
            &                         nElems     = nElems_fluid,         &
            &                         nVars      = nScalars,             &
            &                         nDims      = 2                     )

          deallocate(state_der_Q)
        end select

      case('navier_stokes_2d', 'filtered_navier_stokes_2d')
        select case(scheme%basisType)
        case(Q_space)
          ! Projection of the numerical flux
          ! ... x faces
          call modg_2d_project_stabViscNumFluxX_Q(       &
            & numFlux       = facedata%faceFlux(1)%dat,  &
            & faceState     = facedata%faceRep(1)%dat,   &
            & equation      = equation,                  &
            & maxPolyDegree = scheme%maxPolyDegree,      &
            & length        = mesh%length,               &
            & nElems_fluid  = nElems_fluid,              &
            & projection    = kerneldata%state_der,      &
            & poly_proj     = poly_proj                  )
          ! ... y faces
          call modg_2d_project_stabViscNumFluxY_Q(       &
            & numFlux       = facedata%faceFlux(2)%dat,  &
            & faceState     = facedata%faceRep(2)%dat,   &
            & equation      = equation,                  &
            & maxPolyDegree = scheme%maxPolyDegree,      &
            & length        = mesh%length,               &
            & nElems_fluid  = nElems_fluid,              &
            & projection    = kerneldata%state_der,      &
            & poly_proj     = poly_proj                  )
        case(P_space)
          allocate(state_der_Q(nElems_fluid,(scheme%maxPolyDegree+1)**2, &
            &                    nScalars))

          call ply_change_poly_space( inspace    = P_space,              &
            &                         instate    = kerneldata%state_der, &
            &                         outstate   = state_der_Q,          &
            &                         maxPolyDeg = scheme%maxPolyDegree, &
            &                         nElems     = nElems_fluid,         &
            &                         nVars      = nScalars,             &
            &                         nDims      = 2                     )
          ! Projection of the numerical flux
          ! ... x faces
          call modg_2d_project_stabViscNumFluxX_Q(      &
            & numFlux       = facedata%faceFlux(1)%dat, &
            & faceState     = facedata%faceRep(1)%dat,  &
            & equation      = equation,                 &
            & maxPolyDegree = scheme%maxPolyDegree,     &
            & length        = mesh%length,              &
            & nElems_fluid  = nElems_fluid,             &
            & projection    = state_der_Q,              &
            & poly_proj     = poly_proj                 )
          ! ... y faces
          call modg_2d_project_stabViscNumFluxY_Q(      &
            & numFlux       = facedata%faceFlux(2)%dat, &
            & faceState     = facedata%faceRep(2)%dat,  &
            & equation      = equation,                 &
            & maxPolyDegree = scheme%maxPolyDegree,     &
            & length        = mesh%length,              &
            & nElems_fluid  = nElems_fluid,             &
            & projection    = state_der_Q,              &
            & poly_proj     = poly_proj                 )

          call ply_change_poly_space( inspace    = Q_space,              &
            &                         instate    = state_der_Q,          &
            &                         outstate   = kerneldata%state_der, &
            &                         maxPolyDeg = scheme%maxPolyDegree, &
            &                         nElems     = nElems_fluid,         &
            &                         nVars      = nScalars,             &
            &                         nDims      = 2                     )

        deallocate(state_der_Q)
        end select

     case default
       write(logUnit(1),*) 'ERROR in atlmodg_2d_project_NumFlux:'
       write(logUnit(1),*) 'Unknown equation',equation%eq_kind
       write(logUnit(1),*) 'projections of viscous stabilization,'
       write(logUnit(1),*) ' stopping ... '
       call tem_abort()
      end select

    end if

    ! Projection of the penalization term if is not computed somewhere else
    if (penalizationdata%isActive .and. usePenalization) then
      select case(scheme%basisType)
      case(Q_space)
        call modg_2d_project_penalization_Q(nScalars = equation%varSys%nScalars, &
          &                         mesh = mesh, &
          &                         maxPolyDegree = scheme%maxPolyDegree, &
          &                         penalizationdata = penalizationdata, &
          &                         kerneldata = kerneldata )
      case(P_space)
         write(logUnit(1),*) 'ERROR in atl_modg_2d_project_NumFlux: '
         write(logUnit(1),*) 'Penelization not yet implemented for P_space'
         call tem_abort()
     end select
    end if

  end subroutine atl_modg_2d_project_NumFlux


  !> Projection of the penalization terms (in modal representation) to the test
  !! functions.
  subroutine modg_2d_project_penalization_Q( nScalars, mesh, maxPolyDegree, &
    &                                        kerneldata, penalizationdata )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    integer :: iElem, xTestFunc, yTestFunc, testPos, &
             & xAnsFuncMin, xAnsFunc, yAnsFuncMin, yAnsFunc,  &
             & ansPos
    real(kind=rk) :: jacobiDet, xScalProd, yScalProd
    integer :: mpd1, mpd1_square
    ! --------------------------------------------------------------------------

    jacobiDet = (0.5_rk*mesh%length)**2
    mpd1 = maxPolyDegree+1
    mpd1_square = mpd1**2

    !NA!do iElem = 1, size(kerneldata%state_der,1)
    do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

      ! Now, we loop over all the test functions for this element and calculate
      ! the projection of the source terms onto this test functions.
      do testpos=1,mpd1_square
        yTestFunc = (testpos-1)/mpd1 + 1
        xTestFunc = testpos - (yTestFunc-1)*mpd1

        ! Loop over relevant ans functions
        xAnsFuncMin = xTestFunc-2
        if( xAnsFuncMin < 1 ) then
          xAnsFuncMin = xTestFunc
        end if
        do xAnsFunc = xAnsFuncMin,xTestFunc,2
          ! Loop over relevant ans functions
          yAnsFuncMin = yTestFunc-2
          if( yAnsFuncMin < 1 ) then
            yAnsFuncMin = yTestFunc
          end if
          do yAnsFunc = yAnsFuncMin,yTestFunc,2

            ! get position of ansatz functions in the serialized list
            ! of dofs.
  anspos = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

            ! project the current ansatz function onto the test function
            xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
            yScalProd = ply_scalProdDualLeg(yAnsFunc, yTestFunc)
            kernelData%state_der(iElem, testPos, 1:nScalars)                   &
              & = kernelData%state_der(iElem, testPos, 1:nScalars)             &
              &   + xScalProd * yScalProd                                      &
              &   * jacobiDet                                                  &
              &   * penalizationdata%penalization_data(iElem, ansPos,1:nScalars)
          end do ! y ansatz functions
        end do ! x ansatz functions

      end do ! test functions

    end do ! elem loop

  end subroutine modg_2d_project_penalization_Q


  !> Projection of the source terms (in modal representation) to
  !! the test functions.
  subroutine atl_modg_2d_project_source( sourcedata, nScalars, mesh,      &
    &                                    scheme, kerneldata, currentLevel )
    ! --------------------------------------------------------------------------

    !> The modal representation of the source
    type(atl_source_type), intent(in) :: sourcedata
    !> The number scalar variables in the equation system.
    integer, intent(in) :: nScalars
    !> The cubical elements.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The parameters of the MODG scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    !> The data of the kernel. This one is updated by the projection.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The current level
    integer, intent(in) :: currentLevel

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Projection of the source terms onto the test functions.
    select case(scheme%basisType)
    case(Q_space)
      call modg_2d_project_source_Q(            &
        & sourcedata    = sourcedata,           &
        & nScalars      = nScalars,             &
        & mesh          = mesh,                 &
        & maxPolyDegree = scheme%maxPolyDegree, &
        & kerneldata    = kerneldata,           &
        & currentLevel  = currentLevel          )
    case(P_space)
      call modg_2d_project_source_P(            &
        & sourcedata    = sourcedata,           &
        & nScalars      = nScalars,             &
        & mesh          = mesh,                 &
        & maxPolyDegree = scheme%maxPolyDegree, &
        & kerneldata    = kerneldata,           &
        & currentLevel  = currentLevel          )
    end select

  end subroutine atl_modg_2d_project_source


  !> Projection of the source terms (in modal representation) to
  !! the test functions.
  subroutine modg_2d_project_source_Q( nScalars, sourcedata, maxPolyDegree, &
    &                                  mesh, kerneldata, currentLevel       )
    ! --------------------------------------------------------------------------

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

    ! --------------------------------------------------------------------------
    integer :: iElem, elemPos, xTestFunc, yTestFunc, testPos, &
      & xAnsFuncMin, xAnsFunc, yAnsFuncMin, yAnsFunc,  &
      & ansPos, varPos, iSource, nSourceElems
    real(kind=rk) :: jacobiDet, xScalProd, yScalProd
    integer :: mpd1, mpd1_square
    ! --------------------------------------------------------------------------
    jacobiDet = (0.5_rk*mesh%length)**2
    mpd1 = maxPolyDegree+1
    mpd1_square = mpd1**2

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

          ! Now, we loop over all the test functions for this element and
          ! calculate the projection of the source terms onto this
          ! test functions.
          do testpos=1,mpd1_square
            yTestFunc = (testpos-1)/mpd1 + 1
            xTestFunc = testpos - (yTestFunc-1)*mpd1

            ! Loop over relevant ans functions
            xAnsFuncMin = xTestFunc-2
            if( xAnsFuncMin < 1 ) then
              xAnsFuncMin = xTestFunc
            end if
            do xAnsFunc = xAnsFuncMin,xTestFunc,2
              ! Loop over relevant ans functions
              yAnsFuncMin = yTestFunc-2
              if( yAnsFuncMin < 1 ) then
                yAnsFuncMin = yTestFunc
              end if
              do yAnsFunc = yAnsFuncMin,yTestFunc,2

                ! get position of ansatz functions in the serialized list
                ! of dofs.
  anspos = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

                ! project the current ansatz function onto the test function
                xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
                yScalProd = ply_scalProdDualLeg(yAnsFunc, yTestFunc)
                kernelData%state_der(elemPos, testPos,  varPos)               &
                  & = kernelData%state_der(elemPos, testPos,  varPos)         &
                  &   + xScalProd * yScalProd                                 &
                  &     * jacobiDet                                           &
                  &     * sourcedata%method(iSource)%val(iElem, ansPos, varPos)
              end do ! y ansatz functions
            end do ! x ansatz functions

          end do ! test functions

        end do ! elem loop
      end do ! variable loop
    end do ! source loop

  end subroutine modg_2d_project_source_Q


  !> Projection of the source terms (in modal representation) to
  !! the test functions.
  subroutine modg_2d_project_source_P( nScalars, sourcedata, maxPolyDegree, &
    &                                  mesh, kerneldata, currentLevel       )
    ! --------------------------------------------------------------------------

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

    ! --------------------------------------------------------------------------
    integer :: iElem, elemPos, xTestFunc, yTestFunc, testPos, &
      & xAnsFuncMin, xAnsFunc, yAnsFuncMin, yAnsFunc,  &
      & ansPos, varPos, iSource, nSourceElems, testPosMax
    real(kind=rk) :: jacobiDet, xScalProd, yScalProd
    ! --------------------------------------------------------------------------
    jacobiDet = (0.5_rk*mesh%length)**2

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

          ! Now, we loop over all the test functions for this element and
          ! calculate the projection of the source terms onto this
          ! test functions.
          xTestFunc = 1
          yTestFunc = 1
  testposmax = ((maxpolydegree)+1)*((maxpolydegree)+2)/2
          do testpos=1, testPosMax

            ! Loop over relevant ans functions
            xAnsFuncMin = xTestFunc-2
            if( xAnsFuncMin < 1 ) then
              xAnsFuncMin = xTestFunc
            end if
            do xAnsFunc = xAnsFuncMin,xTestFunc,2
              ! Loop over relevant ans functions
              yAnsFuncMin = yTestFunc-2
              if( yAnsFuncMin < 1 ) then
                yAnsFuncMin = yTestFunc
              end if
              do yAnsFunc = yAnsFuncMin,yTestFunc,2

                ! get position of ansatz functions in the serialized list
                ! of dofs.
  ! integer divisions are no mistake here.
  anspos = ((((xansfunc - 1) + (yansfunc - 1))            &
    &   * (((xansfunc - 1) + (yansfunc - 1)) + 1)) / 2 + 1) &
    & + (yansfunc - 1)

                ! project the current ansatz function onto the test function
                xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
                yScalProd = ply_scalProdDualLeg(yAnsFunc, yTestFunc)
                kernelData%state_der(elemPos, testPos,  varPos)               &
                  & = kernelData%state_der(elemPos, testPos,  varPos)         &
                  &   + xScalProd * yScalProd                                 &
                  &     * jacobiDet                                           &
                  &     * sourcedata%method(iSource)%val(iElem, ansPos, varPos)
              end do ! y ansatz functions
            end do ! x ansatz functions

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (xtestfunc .ne. 1) then
    ! next item
    xtestfunc = xtestfunc - 1
    ytestfunc = ytestfunc + 1
  else
    ! next layer
    xtestfunc = ytestfunc + 1
    ytestfunc = 1
  end if
          end do ! test functions
        end do ! elem loop
      end do ! variable loop
    end do ! source loop

  end subroutine modg_2d_project_source_P


  !> Projection of the numerical flux in x direction onto the testfunctions.
  subroutine modg_2d_project_stabViscNumFluxX_Q( numFlux, faceState, &
    &                                 equation, maxPolyDegree, length, &
    &                                 nElems_fluid, projection, poly_proj )
    ! --------------------------------------------------------------------------
    !> The numerical flux on the faces in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: numFlux(:,:,:,:)
    !> The state on the faces in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    !> The equation system under consideration
    type(atl_equations_type), intent(in) :: equation
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The element index
    integer, intent(in) :: nElems_fluid
    !> The numerical flux projected onto the test functions.
    real(kind=rk), intent(inout) :: projection(:,:,:)
    !> Projection for the current level
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! --------------------------------------------------------------------------
    integer :: iVar, iElem, iPoint, testPos, iVP
    integer :: d, e, f, g, h, i
    integer :: nScalars, nPoints, nPVars, nOversamp
    real(kind=rk), allocatable :: flux_left(:,:), pointVal_flux_left(:,:), &
      & pointVal_left(:,:), state_left(:,:), &
      & p_a_left(:), p_b_left(:), &
      & nodal_a_left(:,:), nodal_b_left(:,:), &
      & modalA_left(:,:), modalB_left(:,:)
    real(kind=rk), allocatable :: flux_right(:,:), pointVal_flux_right(:,:), &
      & pointVal_right(:,:), state_right(:,:), &
      & p_a_right(:), p_b_right(:), &
      & nodal_a_right(:,:), nodal_b_right(:,:), &
      & modalA_right(:,:), modalB_right(:,:)
    real(kind=rk) :: velocity_left(2), velocity_right(2), jacobiDet, &
      & testX_val_left, testX_grad_val_left, &
      & testX_val_right, testX_grad_val_right, &
      & outerNormalLeft, outerNormalRight
    ! --------------------------------------------------------------------------

    nScalars = equation%varSys%nScalars
    nOversamp = poly_proj%body_1D%oversamp_dofs
    nPoints = poly_proj%body_1D%nquadpoints

    nPVars = (maxPolyDegree+1)*equation%varSys%nScalars

    allocate( flux_left(nOversamp, nScalars), flux_right(nOversamp, nScalars) )
    allocate( state_left(nOversamp, nScalars), &
      & state_right(nOversamp, nScalars)       )
    allocate( pointVal_flux_left(nPoints, nScalars), &
      & pointVal_flux_right(nPoints, nScalars)       )
    allocate( pointVal_left(nPoints, nScalars), &
      & pointVal_right(nPoints, nScalars)       )
    allocate( nodal_a_left(nPoints, nScalars), &
      & nodal_b_left(nPoints, nScalars),       &
      & nodal_a_right(nPoints, nScalars),      &
      & nodal_b_right(nPoints, nScalars)       )
    allocate( modalA_left(nOversamp,nScalars), &
      & modalB_left(nOversamp, nScalars),      &
      & modalA_right(nOversamp,nScalars),      &
      & modalB_right(nOversamp, nScalars)      )
    allocate( p_a_left(nScalars), &
      & p_b_left(nScalars),       &
      & p_a_right(nScalars),      &
      & p_b_right(nScalars)       )

    ! The outer unit normals
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)

    ! The 1D Jacobi determinant
    jacobiDet = (0.5_rk*length)

    elementLoop: do iElem = 1,nElems_fluid

      do f=lbound(state_left,2),ubound(state_left,2)
        state_left(:,f) = 0.0_rk
        do g=lbound(state_right,2),ubound(state_right,2)
          state_right(:,g) = 0.0_rk
          do h=lbound(flux_left,2),ubound(flux_left,2)
            flux_left(:,h) = 0.0_rk
            do i=lbound(flux_right,2),ubound(flux_right,2)
              flux_right(:,i) = 0.0_rk
            end do
          end do
        end do
      end do

      do iVP = 1,nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)
        ! Get u
        ! ... for left face
        state_left(iPoint,iVar) = faceState(iElem,iPoint,iVar,1)
        ! ... for right face
        state_right(iPoint,iVar) = faceState(iElem,iPoint,iVar,2)
        ! Caluclate (u^* - u)
        ! ... for the left face
        flux_left(iPoint,iVar) = numFlux(iElem,iPoint,iVar+nScalars,1)  &
          & - faceState(iElem,iPoint,iVar,1)
        ! ... for the right face
        flux_right(iPoint,iVar) = numFlux(iElem,iPoint,iVar+nScalars,2) &
          & - faceState(iElem,iPoint,iVar,2)
      end do

      ! Transform u to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me         = poly_proj,     &
        &                       dim        = 1,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = pointVal_left, &
        &                       modal_data = state_left     )
      ! ... for the right face
      call ply_poly_project_m2n(me         = poly_proj,      &
        &                       dim        = 1,              &
        &                       nVars      = nScalars,       &
        &                       nodal_data = pointVal_right, &
        &                       modal_data = state_right     )

      ! Transform (u^* - u) to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me         = poly_proj,          &
        &                       dim        = 1,                  &
        &                       nVars      = nScalars,           &
        &                       nodal_data = pointVal_flux_left, &
        &                       modal_data = flux_left           )
      ! ... for the right face
      call ply_poly_project_m2n(me         = poly_proj,           &
        &                       dim        = 1,                   &
        &                       nVars      = nScalars,            &
        &                       nodal_data = pointVal_flux_right, &
        &                       modal_data = flux_right           )

      ! Loop over all the points
      pointLoop: do iPoint = 1, nPoints

        ! Caculate velocity at this point
        ! ... for the left face
        velocity_left(1:2) = pointVal_left(iPoint,2:3)/pointVal_left(iPoint,1)
        ! ... for the right face
        velocity_right(1:2) = &
          & pointVal_right(iPoint,2:3)/pointVal_right(iPoint,1)

        select case(equation%eq_kind)
        case('navier_stokes_2d')
          ! Build matrix-vector product of nu_11 and values at current point
          ! ... for the left face
          nodal_a_left(iPoint,:) = atl_mult_nu11_NavierStokes_2d( &
            & density   = pointVal_left(iPoint,1),                &
            & velocity  = velocity_left,                          &
            & totEnergy = pointVal_left(iPoint,4),                &
            & inVec     = pointVal_flux_left(iPoint,:),           &
            & mu        = equation%NavierStokes%mu,               &
            & lambda    = equation%NavierStokes%lambda,           &
            & thermCond = equation%NavierStokes%therm_cond,       &
            & heatCap   = equation%euler%cv                       )
          ! ... for the right face
          nodal_a_right(iPoint,:) = atl_mult_nu11_NavierStokes_2d( &
            & density   = pointVal_right(iPoint,1),                &
            & velocity  = velocity_right,                          &
            & totEnergy = pointVal_right(iPoint,4),                &
            & inVec     = pointVal_flux_right(iPoint,:),           &
            & mu        = equation%NavierStokes%mu,                &
            & lambda    = equation%NavierStokes%lambda,            &
            & thermCond = equation%NavierStokes%therm_cond,        &
            & heatCap   = equation%euler%cv                        )

          ! Build matrix-vector product of nu_21 and values at current point
          ! ... for the left face
          nodal_b_left(iPoint,:) = atl_mult_nu21_NavierStokes_2d( &
            & density  = pointVal_left(iPoint,1),                 &
            & velocity = velocity_left,                           &
            & inVec    = pointVal_flux_left(iPoint,:),            &
            & mu       = equation%NavierStokes%mu,                &
            & lambda   = equation%NavierStokes%lambda             )
          ! ... for the right face
          nodal_b_right(iPoint,:) = atl_mult_nu21_NavierStokes_2d( &
            & density   = pointVal_right(iPoint,1),                &
            & velocity  = velocity_right,                          &
            & inVec     = pointVal_flux_right(iPoint,:),           &
            & mu        = equation%NavierStokes%mu,                &
            & lambda    = equation%NavierStokes%lambda             )

       case('filtered_navier_stokes_2d')
          ! Build matrix-vector product of nu_11 and values at current point
          ! ... for the left face
          nodal_a_left(iPoint,:) = atl_mult_nu11_Rans_2d(     &
            & velocity    = velocity_left,                    &
            & state       = pointVal_left(iPoint,:),          &
            & inVec       = pointVal_flux_left(iPoint,:),     &
            & isenCoeff   = equation%euler%isen_coef,         &
            & mu          = equation%NavierStokes%mu,         &
            & lambda      = equation%NavierStokes%lambda,     &
            & thermCond   = equation%NavierStokes%therm_cond, &
            & rans_params = equation%FiltNavierStokes%rans,   &
            & heatCap     = equation%euler%cv                 )
          ! ... for the right face
          nodal_a_right(iPoint,:) = atl_mult_nu11_Rans_2d(    &
            & velocity    = velocity_right,                   &
            & state       = pointVal_right(iPoint,:),         &
            & inVec       = pointVal_flux_right(iPoint,:),    &
            & isenCoeff   = equation%euler%isen_coef,         &
            & mu          = equation%NavierStokes%mu,         &
            & lambda      = equation%NavierStokes%lambda,     &
            & thermCond   = equation%NavierStokes%therm_cond, &
            & rans_params = equation%FiltNavierStokes%rans,   &
            & heatCap     = equation%euler%cv                 )

          ! Build matrix-vector product of nu_21 and values at current point
          ! ... for the left face
          nodal_b_left(iPoint,:) = atl_mult_nu21_Rans_2d(  &
            & velocity    = velocity_left,                 &
            & state       = pointVal_left(iPoint,:),       &
            & inVec       = pointVal_flux_left(iPoint,:),  &
            & mu          = equation%NavierStokes%mu,      &
            & lambda      = equation%NavierStokes%lambda,  &
            & rans_params = equation%FiltNavierStokes%rans )
          ! ... for the right face
          nodal_b_right(iPoint,:) = atl_mult_nu21_Rans_2d( &
            & velocity    = velocity_right,                &
            & state       = pointVal_right(iPoint,:),      &
            & inVec       = pointVal_flux_right(iPoint,:), &
            & mu          = equation%NavierStokes%mu,      &
            & lambda      = equation%NavierStokes%lambda,  &
            & rans_params = equation%FiltNavierStokes%rans )

        end select

      end do pointLoop

      ! Transform nodal_a and nodal_b back to modal space for projections
      ! ... for the left face
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 1,            &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_a_left, &
        &                       modal_data = modalA_left   )
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 1,            &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_b_left, &
        &                       modal_data = modalB_left   )
      ! ... for the right face
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 1,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_a_right, &
        &                       modal_data = modalA_right   )
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 1,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_b_right, &
        &                       modal_data = modalB_right   )

      ! Project onto all the test functions
      xTestLoop: do d = 1, maxPolyDegree+1
        ! Evaluate test function in x direction at face coordinate
        ! ... for the left face
        testX_val_left = ply_faceValLeftBndTest(d)
        testX_grad_val_left = ply_faceValLeftBndTestGrad(d)
        ! ... for the right face
        testX_val_right = ply_faceValRightBndTest(d)
        testX_grad_val_right = ply_faceValRightBndTestGrad(d)

        ! Now we compute the projections. Since the test functions
        ! change when e > 2 we moved the first two projections of
        ! of the projection loop.

        ! Project onto the first test function in y direction
        e = 1
        ! ... for the left face
        p_a_left(:) = jacobiDet * testX_grad_val_left &
          & * ply_scalProdDualLeg(e,e)                &
          & * modalA_left(e,:)
        testPos = d + (e-1)*(maxPolyDegree+1)
        projection(iElem,testPos,:) = projection(iElem,testPos,:) &
          & - outerNormalLeft * ( p_a_left(:) )
        ! ... for the right face
        p_a_right(:) = jacobiDet * testX_grad_val_right &
          & * ply_scalProdDualLeg(e,e)                  &
          & * modalA_right(e,:)
        testPos = d + (e-1)*(maxPolyDegree+1)
        projection(iElem,testPos,:) = projection(iElem,testPos,:) &
          & - outerNormalRight * ( p_a_right(:) )

        if( maxPolyDegree+1 > 1) then
          ! Project onto the second test function in y direction
          e = 2
          ! ... for the left face
          p_a_left(:) = jacobiDet * testX_grad_val_left &
            & * ply_scalProdDualLeg(e,e)                &
            & * modalA_left(e,:)
          p_b_left(:) = testX_val_left * ply_scalProdDualLegDiff(e-1,e) &
            & * modalB_left(e-1,:)
          testPos = d + (e-1)*(maxPolyDegree+1)
          projection(iElem,testPos,:) = projection(iElem,testPos,:) &
            & - outerNormalLeft * ( p_a_left(:) + p_b_left(:) )
          ! ... for the right face
          p_a_right(:) = jacobiDet * testX_grad_val_right &
            & * ply_scalProdDualLeg(e,e) * modalA_right(e,:)
          p_b_right(:) = testX_val_right* ply_scalProdDualLegDiff(e-1,e) &
            & * modalB_right(e-1,:)
          testPos = d + (e-1)*(maxPolyDegree+1)
          projection(iElem,testPos,:) = projection(iElem,testPos,:) &
            & - outerNormalRight* ( p_a_right(:) + p_b_right(:) )

          ! Now, we project onto all the other test functions in y direction
          yTestLoop: do e = 3, maxPolyDegree+1
            ! ... for the left face
            p_a_left(:) = jacobiDet * testX_grad_val_left   &
              & * ( ply_scalProdDualLeg(e,e) * modalA_left(e,:) &
              & + ply_scalProdDualLeg(e-2,e) * modalA_left(e-2,:) )
            p_b_left(:) = testX_val_left * ply_scalProdDualLegDiff(e-1,e) &
              & * modalB_left(e-1,:)
            testPos = d + (e-1)*(maxPolyDegree+1)
            projection(iElem,testPos,:) = projection(iElem,testPos,:) &
              & - outerNormalLeft * ( p_a_left(:) + p_b_left(:) )
            ! ... for the right face
            p_a_right(:) = jacobiDet * testX_grad_val_right  &
              & * ( ply_scalProdDualLeg(e,e) * modalA_right(e,:) &
              & + ply_scalProdDualLeg(e-2,e) * modalA_right(e-2,:) )
            p_b_right(:) = testX_val_right * ply_scalProdDualLegDiff(e-1,e) &
              & * modalB_right(e-1,:)
            testPos = d + (e-1)*(maxPolyDegree+1)
            projection(iElem,testPos,:) = projection(iElem,testPos,:) &
              & - outerNormalRight * ( p_a_right(:) + p_b_right(:) )
          end do yTestLoop
        end if
      end do xTestLoop

    end do elementLoop

  end subroutine modg_2d_project_stabViscNumFluxX_Q


  !> Projection of the numerical flux in y direction onto the testfunctions.
  subroutine modg_2d_project_stabViscNumFluxY_Q( numFlux, faceState, &
    &                                 equation, maxPolyDegree, length, &
    &                                 nElems_fluid, projection, poly_proj )
    ! --------------------------------------------------------------------------
    !> The numerical flux on the faces in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: numFlux(:,:,:,:)
    !> The state on the faces in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    !> The equation system under consideration
    type(atl_equations_type), intent(in) :: equation
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The element index
    integer, intent(in) :: nElems_fluid
    !> The numerical flux projected onto the test functions.
    real(kind=rk), intent(inout) :: projection(:,:,:)
    !> Projection for the current level
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! --------------------------------------------------------------------------
    integer :: iVar, iElem, iPoint, testPos, iVP
    integer :: d, e, f, g, h, i
    integer :: nScalars, nPoints, nPVars, nOversamp
    real(kind=rk), allocatable :: flux_left(:,:), pointVal_flux_left(:,:), &
      & pointVal_left(:,:), state_left(:,:),                               &
      & p_a_left(:), p_b_left(:),                                          &
      & nodal_a_left(:,:), nodal_b_left(:,:),                              &
      & modalA_left(:,:), modalB_left(:,:)
    real(kind=rk), allocatable :: flux_right(:,:), pointVal_flux_right(:,:), &
      & pointVal_right(:,:), state_right(:,:),                               &
      & p_a_right(:), p_b_right(:),                                          &
      & nodal_a_right(:,:), nodal_b_right(:,:),                              &
      & modalA_right(:,:), modalB_right(:,:)
    real(kind=rk) :: velocity_left(2), velocity_right(2), jacobiDet, &
      & testY_val_left, testY_grad_val_left,                         &
      & testY_val_right, testY_grad_val_right,                       &
      & outerNormalLeft, outerNormalRight
    ! --------------------------------------------------------------------------

    nScalars = equation%varSys%nScalars
    nOversamp = poly_proj%body_1D%oversamp_dofs
    nPoints = poly_proj%body_1D%nquadpoints

    nPVars = (maxPolyDegree+1)*equation%varSys%nScalars

    allocate( flux_left(nOversamp, nScalars), flux_right(nOversamp, nScalars) )
    allocate( state_left(nOversamp, nScalars), &
      & state_right(nOversamp, nScalars)       )
    allocate( pointVal_flux_left(nPoints, nScalars), &
      & pointVal_flux_right(nPoints, nScalars)       )
    allocate( pointVal_left(nPoints, nScalars), &
      & pointVal_right(nPoints, nScalars)       )
    allocate( nodal_a_left(nPoints, nScalars), &
      & nodal_b_left(nPoints, nScalars),       &
      & nodal_a_right(nPoints, nScalars),      &
      & nodal_b_right(nPoints, nScalars)       )
    allocate( modalA_left(nOversamp,nScalars), &
      & modalB_left(nOversamp, nScalars),      &
      & modalA_right(nOversamp,nScalars),      &
      & modalB_right(nOversamp, nScalars)      )
    allocate( p_a_left(nScalars), p_b_left(nScalars), &
      & p_a_right(nScalars), p_b_right(nScalars)      )

    ! The outer unit normals
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)

    ! The 1D Jacobi determinant
    jacobiDet = (0.5_rk*length)

    elementLoop: do iElem = 1,nElems_fluid
      do f=lbound(state_left,2),ubound(state_left,2)
        state_left(:,f) = 0.0_rk
        do g=lbound(state_right,2),ubound(state_right,2)
          state_right(:,g) = 0.0_rk
          do h=lbound(flux_left,2),ubound(flux_left,2)
            flux_left(:,h) = 0.0_rk
            do i=lbound(flux_right,2),ubound(flux_right,2)
              flux_right(:,i) = 0.0_rk
            end do
          end do
        end do
      end do

      do iVP = 1,nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)

        ! Get u
        ! ... for left face
        state_left(iPoint,iVar) = faceState(iElem,iPoint,iVar,1)
        ! ... for right face
        state_right(iPoint,iVar) = faceState(iElem,iPoint,iVar,2)
        ! Caluclate (u^* - u)
        ! ... for the left face
        flux_left(iPoint,iVar) = numFlux(iElem,iPoint,iVar+nScalars,1)  &
          & - faceState(iElem,iPoint,iVar,1)
        ! ... for the right face
        flux_right(iPoint,iVar) = numFlux(iElem,iPoint,iVar+nScalars,2) &
          & - faceState(iElem,iPoint,iVar,2)
      end do

      ! Transform  u to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me         = poly_proj,     &
        &                       dim        = 1,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = pointVal_left, &
        &                       modal_data = state_left     )
      ! ... for the right face
      call ply_poly_project_m2n(me         = poly_proj,      &
        &                       dim        = 1 ,             &
        &                       nVars      = nScalars,       &
        &                       nodal_data = pointVal_right, &
        &                       modal_data = state_right     )

      ! Transform (u^* - u) to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me         = poly_proj,          &
        &                       dim        = 1 ,                 &
        &                       nVars      = nScalars,           &
        &                       nodal_data = pointVal_flux_left, &
        &                       modal_data = flux_left           )
      ! ... for the right face
      call ply_poly_project_m2n(me         = poly_proj,           &
        &                       dim        = 1 ,                  &
        &                       nVars      = nScalars,            &
        &                       nodal_data = pointVal_flux_right, &
        &                       modal_data = flux_right           )

      ! Loop over all the points
      pointLoop: do iPoint = 1, nPoints

        ! Caculate velocity at this point
        ! ... for the left face
        velocity_left(1:2) = pointVal_left(iPoint,2:3)/pointVal_left(iPoint,1)
        ! ... for the right face
        velocity_right(1:2) = pointVal_right(iPoint,2:3) &
          & / pointVal_right(iPoint,1)

        ! Build matrix-vector product of nu_12 and values at current point
        ! ... for the left face
        select case(equation%eq_kind)
        case('navier_stokes_2d')
          nodal_a_left(iPoint,:) = atl_mult_nu12_NavierStokes_2d( &
            & density  = pointVal_left(iPoint,1),                 &
            & velocity = velocity_left,                           &
            & inVec    = pointVal_flux_left(iPoint,:),            &
            & mu       = equation%NavierStokes%mu,                &
            & lambda   = equation%NavierStokes%lambda             )
          ! ... for the right face
          nodal_a_right(iPoint,:) = atl_mult_nu12_NavierStokes_2d( &
            & density  = pointVal_right(iPoint,1),                 &
            & velocity = velocity_right,                           &
            & inVec    = pointVal_flux_right(iPoint,:),            &
            & mu       = equation%NavierStokes%mu,                 &
            & lambda   = equation%NavierStokes%lambda              )

          ! Build matrix-vector product of nu_22 and values at current point
          ! ... for the left face
          nodal_b_left(iPoint,:) = atl_mult_nu22_NavierStokes_2d( &
            & density   = pointVal_left(iPoint,1),                &
            & velocity  = velocity_left,                          &
            & totEnergy = pointVal_left(iPoint,4),                &
            & inVec     = pointVal_flux_left(iPoint,:),           &
            & mu        = equation%NavierStokes%mu,               &
            & lambda    = equation%NavierStokes%lambda,           &
            & thermCond = equation%NavierStokes%therm_cond,       &
            & heatCap   = equation%euler%cv                       )
          ! ... for the right face
          nodal_b_right(iPoint,:) = atl_mult_nu22_NavierStokes_2d( &
            & density   = pointVal_right(iPoint,1),                &
            & velocity  = velocity_right,                          &
            & totEnergy = pointVal_right(iPoint,4),                &
            & inVec     = pointVal_flux_right(iPoint,:),           &
            & mu        = equation%NavierStokes%mu,                &
            & lambda    = equation%NavierStokes%lambda,            &
            & thermCond = equation%NavierStokes%therm_cond,        &
            & heatCap   = equation%euler%cv                        )
        case('filtered_navier_stokes_2d')
          nodal_a_left(iPoint,:) = atl_mult_nu12_Rans_2d(  &
            & velocity    = velocity_left,                 &
            & state       = pointVal_left(iPoint,:),       &
            & inVec       = pointVal_flux_left(iPoint,:),  &
            & mu          = equation%NavierStokes%mu,      &
            & lambda      = equation%NavierStokes%lambda,  &
            & rans_params = equation%FiltNavierStokes%rans )
          ! ... for the right face
          nodal_a_right(iPoint,:) = atl_mult_nu12_Rans_2d( &
            & velocity    = velocity_right,                &
            & state       = pointVal_right(iPoint,:),      &
            & inVec       = pointVal_flux_right(iPoint,:), &
            & mu          = equation%NavierStokes%mu,      &
            & lambda      = equation%NavierStokes%lambda,  &
            & rans_params = equation%FiltNavierStokes%rans )

          ! Build matrix-vector product of nu_22 and values at current point
          ! ... for the left face
          nodal_b_left(iPoint,:) = atl_mult_nu22_Rans_2d(     &
            & velocity    = velocity_left,                    &
            & state       = pointVal_left(iPoint,:),          &
            & inVec       = pointVal_flux_left(iPoint,:),     &
            & isenCoeff   = equation%euler%isen_coef,         &
            & mu          = equation%NavierStokes%mu,         &
            & lambda      = equation%NavierStokes%lambda,     &
            & thermCond   = equation%NavierStokes%therm_cond, &
            & rans_params = equation%FiltNavierStokes%rans,   &
            & heatCap     = equation%euler%cv                 )

          ! ... for the right face
          nodal_b_right(iPoint,:) = atl_mult_nu22_Rans_2d(    &
            & velocity    = velocity_right,                   &
            & state       = pointVal_right(iPoint,:),         &
            & inVec       = pointVal_flux_right(iPoint,:),    &
            & isenCoeff   = equation%euler%isen_coef,         &
            & mu          = equation%NavierStokes%mu,         &
            & lambda      = equation%NavierStokes%lambda,     &
            & thermCond   = equation%NavierStokes%therm_cond, &
            & rans_params = equation%FiltNavierStokes%rans,   &
            & heatCap     = equation%euler%cv                 )
        end select

      end do pointLoop


      ! Transform nodal_a and nodal_b back to modal space for projections
      ! ... for the left face
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 1,            &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_a_left, &
        &                       modal_data = modalA_left   )
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 1,            &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_b_left, &
        &                       modal_data = modalB_left   )

      ! ... for the right face
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 1,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_a_right, &
        &                       modal_data = modalA_right   )
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 1,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_b_right, &
        &                       modal_data = modalB_right   )

      ! Project onto all the test functions
      yTestLoop: do e = 1, maxPolyDegree+1
        ! Evaluate test function in x direction at face coordinate
        ! ... for the left face
        testY_val_left = ply_faceValLeftBndTest(e)
        testY_grad_val_left = ply_faceValLeftBndTestGrad(e)
        ! ... for the right face
        testY_val_right = ply_faceValRightBndTest(e)
        testY_grad_val_right = ply_faceValRightBndTestGrad(e)

        ! Now we compute the projections. Since the test functions
        ! change when e > 2 we moved the first two projections of
        ! of the projection loop.

        ! Project onto the first test function in y direction
        d = 1
        ! ... for the left face
        p_b_left(:) = jacobiDet * testY_grad_val_left &
          & * ply_scalProdDualLeg(d,d)                &
          & * modalB_left(d,:)
        testPos = d + (e-1)*(maxPolyDegree+1)
        projection(iElem,testPos,:) = projection(iElem,testPos,:) &
          & - outerNormalLeft * ( p_b_left(:) )
        ! ... for the right face
        p_b_right(:) = jacobiDet * testY_grad_val_right &
          & * ply_scalProdDualLeg(d,d)                  &
          & * modalb_right(d,:)
        testPos = d + (e-1)*(maxPolyDegree+1)
        projection(iElem,testPos,:) = projection(iElem,testPos,:) &
          & - outerNormalRight * ( p_b_right(:) )

        if( maxPolyDegree+1 > 1) then
          ! Project onto the second test function in y direction
          d = 2
          ! ... for the left face
          p_b_left(:) = jacobiDet * testY_grad_val_left &
            & * ply_scalProdDualLeg(d,d)                &
            & * modalb_left(d,:)
          p_a_left(:) = testY_val_left * ply_scalProdDualLegDiff(d-1,d) &
            & * modalA_left(d-1,:)
          testPos = d + (e-1)*(maxPolyDegree+1)
          projection(iElem,testPos,:) = projection(iElem,testPos,:) &
            & - outerNormalLeft * ( p_a_left(:) + p_b_left(:) )
          ! ... for the right face
          p_b_right(:) = jacobiDet * testY_grad_val_right &
            & * ply_scalProdDualLeg(d,d) * modalB_right(d,:)
          p_a_right(:) = testY_val_right* ply_scalProdDualLegDiff(d-1,d) &
            & * modalA_right(d-1,:)
          testPos = d + (e-1)*(maxPolyDegree+1)
          projection(iElem,testPos,:) = projection(iElem,testPos,:) &
            & - outerNormalRight* ( p_a_right(:) + p_b_right(:) )

          ! Now, we project onto all the other test functions in y direction
          xTestLoop: do d = 3, maxPolyDegree+1
            ! ... for the left face
            p_b_left(:) = jacobiDet * testY_grad_val_left     &
              & * ( ply_scalProdDualLeg(d,d) * modalB_left(d,:)   &
              & + ply_scalProdDualLeg(d-2,d) * modalB_left(d-2,:) )
            p_a_left(:) = testY_val_left * ply_scalProdDualLegDiff(d-1,d) &
              & * modalA_left(d-1,:)
            testPos = d + (e-1)*(maxPolyDegree+1)
            projection(iElem,testPos,:) = projection(iElem,testPos,:) &
              & - outerNormalLeft * ( p_a_left(:) + p_b_left(:) )
            ! ... for the right face
            p_b_right(:) = jacobiDet * testY_grad_val_right    &
              & * ( ply_scalProdDualLeg(d,d) * modalB_right(d,:)   &
              & + ply_scalProdDualLeg(d-2,d) * modalB_right(d-2,:) )
            p_a_right(:) = testY_val_right * ply_scalProdDualLegDiff(d-1,d) &
              & * modalA_right(d-1,:)
            testPos = d + (e-1)*(maxPolyDegree+1)
            projection(iElem,testPos,:) = projection(iElem,testPos,:) &
              & - outerNormalRight * ( p_a_right(:) + p_b_right(:) )
          end do xTestLoop
        end if
      end do yTestLoop

    end do elementLoop

  end subroutine modg_2d_project_stabViscNumFluxY_Q


  !> Projection of the numerical flux in x direction onto the testfunctions
  !! for Q_space.
  subroutine modg_2d_project_numFluxX_Q( numFlux, nScalars, maxPolyDegree, &
    &                                    length, nElems_fluid, dl_prod,    &
    &                                    projection                        )
    ! --------------------------------------------------------------------------
    !> The numerical flux on the faces in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: numFlux(:,:,:,:)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The element index
    integer, intent(in) :: nElems_fluid
    !> The numerical flux projected onto the test functions.
    real(kind=rk), intent(inout) :: projection(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    ! --------------------------------------------------------------------------
    integer :: xTestFunc, yTestFunc
    integer :: yAnsFunc
    integer :: testPos, ansPos
    integer :: yAnsFuncMin
    real(kind=rk) :: yScalProd
    real(kind=rk) :: outerNormalLeft, outerNormalRight
    real(kind=rk) :: jacobiDetFaceProj
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iElem, iXY
    integer :: min2mpd, nTests
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for the projections of the numerical fluxes onto the
    ! test functions
    jacobiDetFaceProj = (0.5_rk*length)**1
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)

    min2mpd = min(maxPolyDegree+1,2)
    nTests = min2mpd*(maxPolyDegree+1)

    ! Loop over all the test functions and project the numerical flux to them.
    do iXY=1,nTests
      yTestFunc = (iXY-1)/min2mpd + 1
      xTestFunc = iXY - (yTestFunc-1)*min2mpd
      testPos = xTestFunc + (yTestFunc-1)*(maxPolyDegree+1)

      faceValLeft = ply_faceValLeftBndTest(xTestFunc) * outerNormalLeft
      faceValRight = ply_faceValRightBndTest(xTestFunc) * outerNormalRight

      ! get the relevant indices of ansatz functions for the projection
      yAnsFuncMin = 1
      if( yTestFunc <= 2 ) then
        yAnsFuncMin = 2
      end if
      do yAnsFunc = yAnsFuncMin,2
        ! calculate the projection of the ansatz and test function
        yScalProd = dl_prod(yAnsFunc, yTestFunc) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos = ytestfunc+(yansfunc-2)*2                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        elementLoop: do iElem = 1,nElems_fluid
          ! ... for the left face
          projection(iElem,testPos,1:nScalars) =   &
            & projection(iElem,testPos,1:nScalars) &
            & - faceValLeft * yScalProd            &
            & * numFlux(iElem,ansPos,1:nScalars,1)
          ! ... for the right face
          projection(iElem,testPos,1:nScalars) =   &
            & projection(iElem,testPos,1:nScalars) &
            & - faceValRight * yScalProd           &
            & * numFlux(iElem,ansPos,1:nScalars,2)
        end do elementLoop
      end do

    end do

  end subroutine modg_2d_project_numFluxX_Q

  !> Projection of the numerical flux in x direction onto the testfunctions.
  subroutine modg_2d_project_numFluxX_diffTestX_Q( numFluxLeftFace,    &
    & numFluxRightFace, nScalars, maxPolyDegree, length, nElems_fluid, &
    & dl_prod, projection                                              )
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
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The element index
    integer, intent(in) :: nElems_fluid
    !> The numerical flux projected onto the test functions.
    real(kind=rk), intent(inout) :: projection(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    ! --------------------------------------------------------------------------
    integer :: xTestFunc, yTestFunc
    integer :: yAnsFunc
    integer :: testPos, ansPos
    integer :: yAnsFuncMin
    real(kind=rk) :: yScalProd
    real(kind=rk) :: jacobiDetFaceProj
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iElem, iXY
    integer :: nTests, iVar
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for the projections of the numerical fluxes onto the
    ! test functions
    jacobiDetFaceProj = (0.5_rk*length)**1
    !outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    !outerNormalRight = atl_elemfaceToNormal_prp(tem_right)

    nTests = (maxPolyDegree)**2

    ! Loop over all the test functions and project the numerical flux to them.
    do iXY=1,nTests
      yTestFunc = (iXY-1)/maxPolyDegree + 2
      xTestFunc = iXY - (yTestFunc-1)*maxPolyDegree + maxPolyDegree +1
      testPos = xTestFunc + (yTestFunc-1)*(maxPolyDegree+1)

      faceValLeft  = ply_faceValLeftBndgradTest(xTestFunc)
      faceValRight = ply_faceValRightBndgradTest(xTestFunc)

      ! get the relevant indices of ansatz functions for the projection
      yAnsFuncMin = 1
      if( yTestFunc <= 5 ) then
        yAnsFuncMin = 2
      end if

      do yAnsFunc = yAnsFuncMin,2

        ! calculate the projection of the ansatz and test function
        yScalProd = dl_prod(yAnsFunc, yTestFunc) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
        ansPos = yTestFunc-1  + (yAnsFunc-2)*4

        ! buffer the result in kernel data and take care of the outer surface
        ! ... unit normal for the left face
        elementLoop: do iElem = 1,nElems_fluid
          do iVar = 1,nScalars
            ! ... for the left face
            projection(iElem,testPos,iVar) = projection(iElem,testPos,iVar) &
              & - faceValLeft * yScalProd                                   &
              & * numFluxLeftFace(iElem,ansPos,nScalars+iVar,1)
            ! ... for the right face
            projection(iElem,testPos,iVar) = projection(iElem,testPos,iVar) &
              & - faceValRight * yScalProd                                  &
              & * numFluxRightFace(iElem,ansPos,nScalars+iVar,2)
          end do
        end do elementLoop
      end do
    end do

  end subroutine modg_2d_project_numFluxX_diffTestX_Q


  !> Projection of the numerical flux in y direction onto the testfunctions
  !! for Q_space.
  subroutine modg_2d_project_numFluxY_Q( numFlux, nScalars, maxPolyDegree, &
    &                                    length, nElems_fluid, dl_prod,    &
    &                                    projection                        )
    ! --------------------------------------------------------------------------
    !> The numerical flux on the faces in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: numFlux(:,:,:,:)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The element index
    integer, intent(in) :: nElems_fluid
    !> The numerical flux projected onto the test functions.
    real(kind=rk), intent(inout) :: projection(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    ! --------------------------------------------------------------------------
    integer :: xTestFunc, yTestFunc
    integer :: xAnsFunc
    integer :: testPos, ansPos
    integer :: xAnsFuncMin
    real(kind=rk) :: xScalProd
    real(kind=rk) :: outerNormalLeft, outerNormalRight
    real(kind=rk) :: jacobiDetFaceProj
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iElem
    integer :: min2mpd, nTests
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for the projections of the numerical fluxes onto the
    ! test functions
    jacobiDetFaceProj = (0.5_rk*length)**1
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)

    min2mpd = min(maxPolyDegree+1,2)
    nTests = min2mpd*(maxPolyDegree+1)

    ! Loop over all the test functions and project the numerical flux to them.
    do testpos = 1,nTests
      yTestFunc = (testpos-1)/(maxPolyDegree+1) + 1
      xTestFunc = testpos - (yTestFunc-1)*(maxPolyDegree+1)

      faceValLeft = ply_faceValLeftBndTest(yTestFunc) * outerNormalLeft
      faceValRight = ply_faceValRightBndTest(yTestFunc) * outerNormalRight

      ! get the relevant indices of ansatz functions for the projection
      xAnsFuncMin = 1
      if( xTestFunc <= 2 ) then
        xAnsFuncMin = 2
      end if
      do xAnsFunc = xAnsFuncMin,2
        xScalProd = dl_prod(xAnsFunc, xTestFunc) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos = xtestfunc+(xansfunc-2)*2                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ! buffer the result in kernel data and take care of the outer
        !  surface unit normal
        elementLoop: do iElem = 1,nElems_fluid
          ! ... for the left face
          projection(iElem,testPos,1:nScalars) =   &
            & projection(iElem,testPos,1:nScalars) &
            & - faceValLeft * xScalProd            &
            & * numFlux(iElem,ansPos,1:nScalars,1)
          ! ... for the right face
          projection(iElem,testPos,1:nScalars) =   &
            & projection(iElem,testPos,1:nScalars) &
            & - faceValRight * xScalProd           &
            & * numFlux(iElem,ansPos,1:nScalars,2)
        end do elementLoop
      end do
    end do

  end subroutine modg_2d_project_numFluxY_Q


  !> Projection of the numerical flux in y direction onto the testfunctions.
  subroutine modg_2d_project_numFluxY_diffTestY_Q( numFluxLeftFace,    &
    & numFluxRightFace, nScalars, maxPolyDegree, length, nElems_fluid, &
    & dl_prod, projection                                              )
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
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The element index
    integer, intent(in) :: nElems_fluid
    !> The numerical flux projected onto the test functions.
    real(kind=rk), intent(inout) :: projection(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    ! --------------------------------------------------------------------------
    integer :: xTestFunc, yTestFunc, iVar
    integer :: xAnsFunc
    integer :: testPos, ansPos
    integer :: xAnsFuncMin
    real(kind=rk) :: xScalProd
    real(kind=rk) :: jacobiDetFaceProj
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iElem, iXY
    integer :: nTests
    ! --------------------------------------------------------------------------


    ! Jacobi determinant for the projections of the numerical fluxes onto
    !  the test functions
    jacobiDetFaceProj = (0.5_rk*length)**1
    !outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    !outerNormalRight = atl_elemfaceToNormal_prp(tem_right)

    nTests = (maxPolyDegree)**2

    ! Loop over all the test functions and project the numerical flux to them.
    do iXY = 1,nTests

      xTestFunc = (iXY-1)/(maxPolyDegree) + 2
      yTestFunc = iXY - (xTestFunc-1)*maxPolyDegree + maxPolyDegree + 1
      testPos = xTestFunc + (yTestFunc-1)*(maxPolyDegree+1)

      faceValLeft = ply_faceValLeftBndgradTest(yTestFunc)
      faceValRight = ply_faceValRightBndgradTest(yTestFunc)

      ! get the relevant indices of ansatz functions for the projection
      xAnsFuncMin = 1
      if( xTestFunc <= 5 ) then
        xAnsFuncMin = 2
      end if

      do xAnsFunc = xAnsFuncMin,2

        xScalProd = dl_prod(xAnsFunc, xTestFunc) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
        ansPos = (xTestFunc -1) + (xAnsFunc-2)*4

        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        elementLoop: do iElem = 1,nElems_fluid
            do iVar = 1,nScalars
              ! ... for the left face
              projection(iElem,testPos,iVar) = projection(iElem,testPos,iVar) &
                & - faceValLeft * xScalProd                                   &
                & * numFluxLeftFace(iElem,ansPos,nScalars+iVar,1)
              ! ... for the right face
              projection(iElem,testPos,iVar) = projection(iElem,testPos,iVar) &
                & - faceValRight * xScalProd                                  &
                & * numFluxRightFace(iElem,ansPos,nScalars+iVar,2)
          end do
        end do elementLoop
      end do
    end do
  end subroutine modg_2d_project_numFluxY_diffTestY_Q


  !> Projects modal representation of each cell to its faces, i.e.
  !! this subroutine creates a modal representation on the faces.
  subroutine atl_modg_2d_modalVolToModalFace( mesh , statedata, facedata, &
    &                                         equation, modg_2d           )
    ! --------------------------------------------------------------------------
    !> The elements we apply the projection for.
    type(atl_cube_elem_type),  intent(in) :: mesh
    !> Volumetric, modal states for each element.
    type(atl_statedata_type), intent(in) :: statedata
    !> Modal representation on the face (will be updated by this routine for all
    !! fluid elements in mesh).
    type(atl_facedata_type), intent(inout) :: facedata
    !> The equation system under consideration
    type(atl_equations_type), intent(in) :: equation
    !> The parameters of your modg scheme.
    type(atl_modg_2d_scheme_type), intent(in) :: modg_2d
    ! --------------------------------------------------------------------------
    integer :: iDir, spaceDir, nScalars, nElems_fluid
    real(kind=rk), allocatable :: volState_Q(:,:,:)
    ! --------------------------------------------------------------------------

    nScalars = equation%varSys%nScalars
    nElems_fluid = mesh%descriptor%elem%nElems(eT_fluid)

    ! Iterate over all the fluid elements and project its modal representations
    ! to the faces of the element.
    select case(modg_2d%basisType)
      case(Q_space) ! Q tensor product ansatz functions

        do iDir = 1, equation%nDimensions
          facedata%faceRep(iDir) &
            &     %dat(1:mesh%descriptor%elem%nElems(eT_fluid),:,:,:) = 0.0_rk
        end do

        ! Now, iterate over all the faces and project to this face.
        ! ... x and y faces
        do iDir = 1, 2
          ! Project to the face in the current direction.
          spaceDir = qAxis(iDir)
          call atl_modg_2d_volToFace_Q(                      &
            & volState      = statedata%state,               &
            & maxPolyDegree = modg_2d%maxPolyDegree,         &
            & faceDir       = iDir,                          &
            & nScalars      = nScalars,                      &
            & nElems        = nElems_fluid,                  &
            & faceState     = facedata%faceRep(spaceDir)%dat )
        end do
        ! ... x and y faces
        do iDir = 4,5
          ! Project to the face in the current direction.
          spaceDir = qAxis(iDir)
          call atl_modg_2d_volToFace_Q(                      &
            & volState      = statedata%state,               &
            & maxPolyDegree = modg_2d%maxPolyDegree,         &
            & faceDir       = iDir,                          &
            & nScalars      = nScalars,                      &
            & nElems        = nElems_fluid,                  &
            & faceState     = facedata%faceRep(spaceDir)%dat )
        end do


        ! Now, iterate over all the faces and project the gradients to this face.
        if(equation%nDerivatives == 1 ) then
          ! ... x faces
          do iDir = 1, 2
            ! Project to the face in the current direction.
            spaceDir = qAxis(iDir)
            ! Project derivatives to the faces
            call atl_modg_2d_volToFace_grad_Q(                 &
              & volState      = statedata%state,               &
              & maxPolyDegree = modg_2d%maxPolyDegree,         &
              & faceDir       = iDir,                          &
              & nScalars      = nScalars,                      &
              & nElems        = nElems_fluid,                  &
              & elemLength    = mesh%length,                   &
              & faceState     = facedata%faceRep(spaceDir)%dat )
          end do
          ! ... y faces
          do iDir = 4,5
            ! Project to the face in the current direction.
            spaceDir = qAxis(iDir)
            call atl_modg_2d_volToFace_grad_Q(                 &
              & volState      = statedata%state,               &
              & maxPolyDegree = modg_2d%maxPolyDegree,         &
              & faceDir       = iDir,                          &
              & nScalars      = nScalars,                      &
              & nElems        = nElems_fluid,                  &
              & elemLength    = mesh%length,                   &
              & faceState     = facedata%faceRep(spaceDir)%dat )
          end do
        end if

      case(P_space) ! P tensor product ansatz functions

        do iDir = 1, equation%nDimensions
          facedata%faceRep(iDir) &
            &     %dat(1:mesh%descriptor%elem%nElems(eT_fluid),:,:,:) = 0.0_rk
        end do

        ! Now, iterate over all the faces and project to this face.
        ! ... x and y faces
        do iDir = 1, 2
          ! Project to the face in the current direction.
          spaceDir = qAxis(iDir)
          call atl_modg_2d_volToFace_P(                      &
            & volState      = statedata%state,               &
            & maxPolyDegree = modg_2d%maxPolyDegree,         &
            & faceDir       = iDir,                          &
            & nScalars      = nScalars,                      &
            & nElems        = nElems_fluid,                  &
            & faceState     = facedata%faceRep(spaceDir)%dat )
        end do
        ! ... x and y faces
        do iDir = 4,5
          ! Project to the face in the current direction.
          spaceDir = qAxis(iDir)
          call atl_modg_2d_volToFace_P(                      &
            & volState      = statedata%state,               &
            & maxPolyDegree = modg_2d%maxPolyDegree,         &
            & faceDir       = iDir,                          &
            & nScalars      = nScalars,                      &
            & nElems        = nElems_fluid,                  &
            & faceState     = facedata%faceRep(spaceDir)%dat )
        end do


        ! Now, iterate over all the faces and project the gradients to this face.
        if(equation%nDerivatives == 1 ) then

          allocate(volstate_Q(nElems_fluid,(modg_2d%maxPolyDegree+1)**2, &
            &                 nScalars))

          call ply_change_poly_space( inspace    = P_space,               &
            &                         instate    = statedata%state,       &
            &                         outstate   = volstate_Q,            &
            &                         maxPolyDeg = modg_2d%maxPolyDegree, &
            &                         nElems     = nElems_fluid,          &
            &                         nVars      = nScalars,              &
            &                         nDims      = 2                      )

          ! ... x faces
          do iDir = 1, 2
            ! Project to the face in the current direction.
            spaceDir = qAxis(iDir)
            ! Project derivatives to the faces
            call atl_modg_2d_volToFace_grad_Q(                 &
              & volState      = volstate_Q,                    &
              & maxPolyDegree = modg_2d%maxPolyDegree,         &
              & faceDir       = iDir,                          &
              & nScalars      = nScalars,                      &
              & nElems        = nElems_fluid,                  &
              & elemLength    = mesh%length,                   &
              & faceState     = facedata%faceRep(spaceDir)%dat )
          end do
          ! ... y faces
          do iDir = 4,5
            ! Project to the face in the current direction.
            spaceDir = qAxis(iDir)
            call atl_modg_2d_volToFace_grad_Q(                 &
              & volState      = volstate_Q,                    &
              & maxPolyDegree = modg_2d%maxPolyDegree,         &
              & faceDir       = iDir,                          &
              & nScalars      = nScalars,                      &
              & nElems        = nElems_fluid,                  &
              & elemLength    = mesh%length,                   &
              & faceState     = facedata%faceRep(spaceDir)%dat )
          end do
          deallocate(volstate_Q)
        end if

      case default
        call tem_abort( 'ERROR in atl_modg_2d_modalVolToModalFace:' &
          & // 'Unknown tensor product, stopping ...'               )
    end select

  end subroutine atl_modg_2d_modalVolToModalFace


  !> Lift the mean on the face to a positive value if necessary.
  !!
  !! In some equations, some variables may not be allowed to be negative,
  !! for example the density and energy in flow simulations.
  !! With this routine we enforce this limitation for the projected mean
  !! state on the element surfaces.
  subroutine atl_modg_2d_ensure_pos_facemean( nelems_fluid, volState, faceRep, &
    &                                         nScalars, ensure_positivity      )
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
    integer :: iVar, iDir
    real(kind=rk) :: epsfact
    ! --------------------------------------------------------------------------

    epsfact = epsilon(volState(1,1,1))

    do iVar=1,nScalars
      if (ensure_positivity(iVar)) then

        ! We assume the integral mean of the element is positive for all
        ! elements. The integral mean on the surface will now be lifted to
        ! a small fraction of that element mean if it is below that threshold.
        do iDir=1,2
          facerep(iDir)%dat(1:nElems_fluid,1,iVar,1)             &
            & = max( facerep(iDir)%dat(1:nElems_fluid,1,iVar,1), &
            &        epsfact * volState(1:nElems_fluid,1,iVar)   )
          facerep(iDir)%dat(1:nElems_fluid,1,iVar,2)             &
            & = max( facerep(iDir)%dat(1:nElems_fluid,1,iVar,2), &
            &        epsfact * volState(1:nElems_fluid,1,iVar)   )
        end do

      end if
    end do

  end subroutine atl_modg_2d_ensure_pos_facemean


  !> Applies the inverse of the mass matrix for a 2D scheme.
  subroutine atl_modg_2d_invMassMatrix( mesh, equation, kerneldata, statedata, &
    &                                   elementalTimestep, timestep, scheme    )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The equation you solve.
    type(atl_equations_type) :: equation
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
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    ! -------------------------------------------------------------------------!
    integer :: nElems

    nElems = mesh%descriptor%elem%nElems(eT_fluid)

    select case(scheme%basisType)
      case(Q_space)
        call modg_2d_invMassMatrix_Q( mesh, equation, kerneldata, scheme, &
          &                           nElems                              )
      case(P_space)
        call modg_2d_invMassMatrix_P( mesh, kerneldata, scheme, nElems )
    end select

    ! Since this is the last part of the MOdal Discontinuous Galerkin (MODG)
    ! algorithm we can make the final timestep, now.
    call elementalTimestep( timestep, statedata%state,kerneldata )

  end subroutine atl_modg_2d_invMassMatrix

  !> Applies the inverse of the mass matrix for a 2D scheme in Q_space.
  subroutine modg_2d_invMassMatrix_Q( mesh, equation, kerneldata, scheme, &
    &                                 nElems                              )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The equation you solve.
    type(atl_equations_type) :: equation
    !> The data of the kernel.
    type(atl_kerneldata_type) :: kerneldata
    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    integer, intent(in) :: nElems
    ! -------------------------------------------------------------------------!
    ! Loop indices for ansatz functions
    integer :: iAnsX, iAnsY
    ! Positions for the given ansatz functions
    integer :: ansLow, ansPrevLow
    ! Inverse of the determintant of the jacobian of the mapping from reference
    ! to physical element.
    real(kind=rk) :: inv_jacobiDet
    integer :: mpd1_square
    ! -------------------------------------------------------------------------!

    inv_jacobiDet = (2.0_rk/mesh%length)**2
    mpd1_square = (scheme%maxPolyDegree+1)**2

    ! apply the 1D inverse of the mass matrix
    do iAnsY = 1, scheme%maxPolyDegree+1
      do iAnsX = 3, scheme%maxPolyDegree+1
  anslow = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansprevlow = iansx-2                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
        kerneldata%state_der(:nElems, ansLow, :equation%varSys%nScalars) &
          & = kerneldata%state_der(:nElems,                              &
          &                        ansLow,                               &
          &                        :equation%varSys%nScalars)            &
          & + kerneldata%state_der(:nElems,                              &
          &                        ansPrevLow,                           &
          &                        :equation%varSys%nScalars)
      end do
    end do


    do AnsLow=1,mpd1_square
      iAnsX = mod(AnsLow-1, scheme%maxPolyDegree+1) + 1
      kerneldata%state_der(:nElems, ansLow, :equation%varSys%nScalars) &
        & = kerneldata%state_der(:nElems, ansLow, :equation%varSys%nScalars) &
        &   * 0.5_rk*(2*iAnsX-1)
    end do


    ! apply the 2D inverse of the mass matrix
    do iAnsX = 1, scheme%maxPolyDegree+1
      do iAnsY = 3, scheme%maxPolyDegree+1
  anslow = 1                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
        ansLow = ansLow + iAnsX - 1
  ansprevlow = 1                                      &
    &      + ( ( iansy-2-1)                             &
    &      + (1-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
        ansPrevLow = ansPrevLow + iAnsX - 1
        kerneldata%state_der(:nElems, ansLow, :equation%varSys%nScalars) &
          & = kerneldata%state_der(:nElems,                              &
          &                        ansLow,                               &
          &                        :equation%varSys%nScalars)            &
          & + kerneldata%state_der(:nElems,                              &
          &                        ansPrevLow,                           &
          &                        :equation%varSys%nScalars)
      end do
    end do


    do AnsLow=1,mpd1_square
      iAnsY = (AnsLow-1) / (scheme%maxPolyDegree+1) + 1
      kerneldata%state_der(:nElems, ansLow, :equation%varSys%nScalars) &
        & = kerneldata%state_der(:nElems, ansLow, :equation%varSys%nScalars) &
        &   * 0.5_rk*(2*iAnsY-1) * inv_jacobiDet
    end do

  end subroutine modg_2d_invMassMatrix_Q


  !> Applies the inverse of the mass matrix for a 2D scheme in P_space.
  subroutine modg_2d_invMassMatrix_P( mesh, kerneldata, scheme, nElems )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The data of the kernel.
    type(atl_kerneldata_type) :: kerneldata
    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    integer, intent(in) :: nElems
    ! -------------------------------------------------------------------------!
    ! Loop indices for ansatz functions
    integer :: iAnsX, iAnsY
    ! Positions for the given ansatz functions
    integer :: ansLow, ansPrevLow
    ! Inverse of the determintant of the jacobian of the mapping from reference
    ! to physical element.
    real(kind=rk) :: inv_jacobiDet
    ! -------------------------------------------------------------------------!

    inv_jacobiDet = (2.0_rk/mesh%length)**2

  ! apply the 1D inverse of the mass matrix
    do iAnsY = 1, scheme%maxPolyDegree+1
      do iAnsX = 3, scheme%maxPolyDegree+1 - (iAnsY-1)
  ! integer divisions are no mistake here.
  anslow = ((((iansx - 1) + (iansy - 1))            &
    &   * (((iansx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
  ! integer divisions are no mistake here.
  ansprevlow = ((((iansx-2 - 1) + (iansy - 1))            &
    &   * (((iansx-2 - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
        kerneldata%state_der(:nElems, ansLow, :) &
          &        = kerneldata%state_der(:nElems, ansLow, :) &
          &        + kerneldata%state_der(:nElems, ansPrevLow, :)
      end do
    end do


    do iAnsY = 1, scheme%maxPolyDegree+1
      do iAnsX = 1, scheme%maxPolyDegree+1 - (iAnsY-1)
  ! integer divisions are no mistake here.
  anslow = ((((iansx - 1) + (iansy - 1))            &
    &   * (((iansx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
        kerneldata%state_der(:nElems, ansLow,  :) &
          &    = kerneldata%state_der(:nElems, ansLow,  :) &
          &      * 0.5_rk*(2*iAnsX-1)
      end do
    end do


    ! apply the 2D inverse of the mass matrix
    do iAnsX = 1, scheme%maxPolyDegree+1
      do iAnsY = 3, scheme%maxPolyDegree+1 - (iAnsX-1)
  ! integer divisions are no mistake here.
  anslow = ((((iansx - 1) + (iansy - 1))            &
    &   * (((iansx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
  ! integer divisions are no mistake here.
  ansprevlow = ((((iansx - 1) + (iansy-2 - 1))            &
    &   * (((iansx - 1) + (iansy-2 - 1)) + 1)) / 2 + 1) &
    & + (iansy-2 - 1)
        kerneldata%state_der(:nElems, ansLow,  :) &
          &    = kerneldata%state_der(:nElems, ansLow,  :) &
          &    + kerneldata%state_der(:nElems, ansPrevLow,  :)
      end do
    end do


    do iAnsY = 1, scheme%maxPolyDegree+1
      do iAnsX = 1, scheme%maxPolyDegree+1 - (iAnsY-1)
  ! integer divisions are no mistake here.
  anslow = ((((iansx - 1) + (iansy - 1))            &
    &   * (((iansx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
        kerneldata%state_der(:nElems, ansLow,  :) &
          &    = kerneldata%state_der(:nElems, ansLow,  :) &
          &      * 0.5_rk*(2*iAnsY-1) * inv_jacobiDet
      end do
    end do

  end subroutine modg_2d_invMassMatrix_P


end module atl_modg_2d_kernel_module
