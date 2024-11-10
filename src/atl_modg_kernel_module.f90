! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
!!
!!@todo HK: Split this module into three: a common one, a Q part and a P part.
module atl_modg_kernel_module
  use env_module,                 only: rk, long_k, eps
  use tem_compileconf_module,     only: vlen

  use tem_aux_module,             only: tem_abort
  use treelmesh_module,           only: treelmesh_type
  use tem_element_module,         only: eT_fluid
  use tem_param_module,           only: qAxis
  use tem_time_module,            only: tem_time_type
  use tem_faceData_module,        only: tem_left, &
    &                                   tem_right
  use tem_logging_module,         only: logUnit
  use tem_comm_module,            only: tem_commPattern_type
  use tem_comm_env_module,        only: tem_comm_env_type
  use tem_element_module,         only: eT_fluid

  use ply_modg_basis_module,      only: ply_faceValLeftBndTest,      &
    &                                   ply_faceValRightBndTest,     &
    &                                   ply_faceValLeftBndGradTest,  &
    &                                   ply_faceValRightBndGradTest, &
    &                                   ply_scalProdDualLeg,         &
    &                                   ply_scalProdDualLegDiff,     &
    &                                   ply_faceValLeftBndTestGrad,  &
    &                                   ply_faceValRightBndTestGrad
  use ply_dof_module,             only: Q_space, P_space, &
    &                                   ply_change_poly_space
  use ply_poly_project_module,    only: ply_poly_project_type, assignment(=), &
    &                                   ply_poly_project_m2n, &
    &                                   ply_poly_project_n2m

  use atl_equation_module,        only: atl_equations_type
  use atl_cube_elem_module,       only: atl_cube_elem_type
  use atl_kerneldata_module,      only: atl_kerneldata_type, &
    &                                   atl_init_kerneldata, &
    &                                   atl_init_statedata,  &
    &                                   atl_statedata_type
  use atl_parallel_module,        only: atl_init_parallel_module
  use atl_scheme_module,          only: atl_scheme_type, atl_modg_scheme_type
  use atl_elemental_time_integration_module,                        &
    &                             only: atl_elemental_timestep_vec, &
    &                                   atl_timestep_type
  use atl_source_types_module,    only: atl_source_type

  use atl_facedata_module,        only: atl_facedata_type, &
    &                                   atl_faceRep_type,  &
    &                                   atl_elemfaceToNormal_prp
  use atl_boundary_module,        only: atl_level_boundary_type, &
    &                                   atl_get_numBndElems
  use atl_materialPrp_module,     only: atl_material_type
  use atl_penalization_module,    only: atl_penalizationData_type
  use atl_physFluxNvrstk_module,  only: atl_mult_nu11_NavierStokes, &
    &                                   atl_mult_nu21_NavierStokes, &
    &                                   atl_mult_nu31_NavierStokes, &
    &                                   atl_mult_nu12_NavierStokes, &
    &                                   atl_mult_nu22_NavierStokes, &
    &                                   atl_mult_nu32_NavierStokes, &
    &                                   atl_mult_nu13_NavierStokes, &
    &                                   atl_mult_nu23_NavierStokes, &
    &                                   atl_mult_nu33_NavierStokes
  use atl_materialIni_module,     only: atl_update_materialParams
  use atl_volToFace_module,       only: atl_modg_volToFace_grad_Q, &
    &                                   atl_modg_volToFace_P

  implicit none

  private

  public :: atl_init_modg_kernel, &
    & atl_modg_invMassMatrix,     &
    & atl_preprocess_modg_kernel, &
    & atl_modg_modalVolToModalFace

  public :: atl_modg_ensure_pos_facemean

  public :: atl_modg_scaledTransposedInvMassMatrix_Q, &
    & atl_modg_scaledTransposedInvMassMatrix_P,       &
    & atl_modg_scaledTransposedProject_physFlux_Q,    &
    & atl_modg_scaledTransposedProject_physFlux_P

  public :: atl_modg_project_source, &
    &       atl_modg_project_NumFlux

  ! only for the utests
  public :: atl_modg_kernel_utests


contains


  !> Initiate the MODG kernel for cubic elements on all levels.
  subroutine atl_init_modg_kernel(                                        &
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
    integer :: iList, iStab, iDir
    integer :: nScalars, nScalarsState_Comm, nScalarsFlux_Comm, nDer, nDof
    integer :: nTotal
    logical :: stab_reqNeigh
    integer :: nBndStabElems(tree%global%minLevel:tree%global%maxLevel,1:3)
    ! --------------------------------------------------------------------------

    ! the number of scalar variables of our equation system
    nScalars = equation%varSys%nScalars

    ! If the stabilization requires further neighbor information,
    ! we init the buffer for it
    stab_reqNeigh = .false.
    do iStab = 1, size(scheme_list(tree%global%minLevel)%stabilization)
      if ( scheme_list(tree%global%minLevel)%stabilization(iStab) &
        &                                   %reqNeigh             ) then
        stab_reqNeigh = .true.
      end if
    end do

    ! init the kerneldata on each level.
    do iList = tree%global%minLevel, tree%global%maxLevel
      ! The number of derived quantities is the number of polynomial degrees of
      ! freedoms per scalar variable.

      select case(scheme_list(iList)%modg%basisType)
      case(Q_space)
        nDer = (scheme_list(iList)%modg%maxPolyDegree+1)**3
        nDof = nDer
      case(P_space)
  nder = (((scheme_list(ilist)%modg%maxpolydegree) + 1) &
    &   * ((scheme_list(ilist)%modg%maxpolydegree) + 2) &
    &   * ((scheme_list(ilist)%modg%maxpolydegree) + 3)) &
    & / 6
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
        &    maxPolyDegree  = scheme_list(iList)%modg%maxPolyDegree,        &
        &    nDims          = 3,                                            &
        &    poly_space     = scheme_list(iList)%modg%basisType,            &
        &    need_deviation = need_deviation,                               &
        &    need_maxgrad   = equation%requires_gradmax                     )

      if(stab_reqNeigh) then
        nBndStabElems = atl_get_numBndElems(      &
          & minLevel      = tree%global%minLevel, &
          & maxLevel      = tree%global%maxLevel, &
          & boundary_list = boundary_stab_list    )
        do iDir = 1,3
          nTotal = mesh_list(iList)%faces%dimByDimDesc(iDir)%nElems &
            &    + nBndStabElems(iList,iDir)
          call atl_init_statedata(                                 &
            &    statedata      = statedata_stab_list(iList,iDir), &
            &    nTotal         = nTotal,                          &
            &    nVars          = nScalars,                        &
            &    nDofs          = nDof,                            &
            &    time           = time                             )
        end do

      end if
    end do

    ! Init the parallel module here, as well. DG requires only face values, so
    ! we init only the face buffers and do not init the cell buffers.
    ! ... for the state, we need the state and all its derivatives
    !     (in all directions)
    nScalarsState_Comm = nScalars*(1+equation%nDerivatives*3)
    ! ... for the flux, we need the numerical flux and the stabilization flux
    nScalarsFlux_Comm = nScalars*(1+equation%nDerivatives)
    call atl_init_parallel_module(                      &
      &    commPattern          = commPattern,          &
      &    scheme               = scheme_list,          &
      &    nValsElem            = nScalars,             &
      &    nValsStateFace       = nScalarsState_Comm,   &
      &    nValsFluxFace        = nScalarsFlux_Comm,    &
      &    cube                 = mesh_list,            &
      &    boundary             = boundary_list,        &
      &    createCellBuffer     = .false.,              &
      &    createFaceBuffer     = .true.,               &
      &    createStabFaceBuffer = .false.,              &
      &    createStabElemBuffer = stab_reqNeigh,        &
      &    nBndStabElems        = nBndStabElems,        &
      &    minLevel             = tree%global%minLevel, &
      &    maxLevel             = tree%global%maxLevel  )

  end subroutine atl_init_modg_kernel



  !> Subroutine to execute the preprocessing for the MODG kernels.
  !! Currently this includes: Convert external source terms to modal
  !! representation.
  subroutine atl_preprocess_modg_kernel( equation, statedata, mesh, scheme, &
    &                                    poly_proj_material, boundary,      &
  &                                      material, proc, commPattern        )
    ! --------------------------------------------------------------------------
    type(atl_equations_type), intent(inout) :: equation
    type(atl_statedata_type), intent(inout) :: statedata
    type(atl_cube_elem_type), intent(inout) :: mesh
    type(atl_scheme_type), intent(inout) :: scheme
    type(atl_material_type), intent(inout) :: material
    !> Global description of the boundary conditions.
    type(atl_level_boundary_type), intent(in) :: boundary
    type(ply_poly_project_type), intent(inout) :: poly_proj_material
    !> mpi communication environment with mpi communicator
    type(tem_comm_env_type) :: proc
    !> mpi communication pattern type
    type(tem_commPattern_type) :: commPattern
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    call atl_update_materialParams( equation    = equation,             &
      &                             mesh        = mesh,                 &
      &                             scheme      = scheme,               &
      &                             boundary    = boundary,             &
      &                             material    = material,             &
      &                             time        = statedata%local_time, &
      &                             poly_proj   = poly_proj_material,   &
      &                             proc        = proc,                 &
      &                             commPattern = commPattern           )

  end subroutine atl_preprocess_modg_kernel




  !> Subroutine to project modal representations of physical flux, numerical
  !! flux and source terms onto test functions.
  subroutine atl_modg_project_NumFlux( mesh, equation, kerneldata, facedata, &
    &                              scheme, poly_proj, dl_prod, dl_prodDiff,  &
    &                              dirvec, penalizationdata, usePenalization )
    ! --------------------------------------------------------------------------
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The representation on the face + representation of the flux.
    type(atl_facedata_type), intent(inout) :: facedata
    !> The parameters of the MODG scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    !> Projection for the current level
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> stored scalar products of the testfunction and anstaz function
    real(kind=rk), intent(in) :: dl_prod(2, scheme%maxPolyDegree+1)
    real(kind=rk), intent(in) :: dl_prodDiff(2, scheme%maxPolyDegree+1)
    !> vector for direction indicators
    integer, intent(in) :: dirVec(3,3)
    !> Volumetric data for the penalization
    type(atl_penalizationData_type), intent(in) :: penalizationdata
    !> Flag to indicate, whether we need to take care of the penalization
    !! terms here or not.
    logical, intent(in) :: usePenalization
    ! --------------------------------------------------------------------------
    integer :: nScalars, nFaces, leftOrRight, nElems_fluid
    !> The direction
    integer :: iDir
    type(atl_faceRep_type), allocatable :: faceReP_Q(:)
    type(atl_faceRep_type), allocatable :: faceFlux_Q(:)
    real(kind=rk), allocatable ::state_der_Q(:,:,:)
    ! --------------------------------------------------------------------------
    nScalars = equation%varSys%nScalars
    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid)


    ! Projection of the numerical flux
    do iDir = 1, 3
      select case(scheme%basisType)
      case(Q_space)
        call modg_project_numFlux_Q(                                         &
          &    nTotalFaces   = size(facedata%faceFlux(iDir)%dat,1),          &
          &    nFaceDofs     = size(facedata%faceFlux(iDir)%dat,2),          &
          &    nScalars      = nScalars,                                     &
          &    faceFlux      = facedata%faceFlux(iDir)%dat(:,:,:nScalars,:), &
          &    maxPolyDegree = scheme%maxPolyDegree,                         &
          &    length        = mesh%length,                                  &
          &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid),        &
          &    dl_prod       = dl_prod,                                      &
          &    projection    = kerneldata%state_der,                         &
          &    dirVec        = dirVec(:,iDir)                                )

      case(P_space)
        call modg_project_numFlux_P(                                         &
          &    nTotalFaces   = size(facedata%faceFlux(iDir)%dat,1),          &
          &    nFaceDofs     = size(facedata%faceFlux(iDir)%dat,2),          &
          &    nScalars      = nScalars,                                     &
          &    faceFlux      = facedata%faceFlux(iDir)%dat(:,:,:nScalars,:), &
          &    maxPolyDegree = scheme%maxPolyDegree,                         &
          &    length        = mesh%length,                                  &
          &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid),        &
          &    dl_prod       = dl_prod,                                      &
          &    projection    = kerneldata%state_der,                         &
          &    dirVec        = dirVec(:,iDir)                                )

      end select
    end do

    if (equation%nDerivatives == 1) then
      select case(equation%eq_kind)
      case('heat')
        select case(scheme%basisType)
        case(Q_space)
        do iDir = 1, 3
          ! Projection of the appropriate  numerical flux onto the derivative
          ! of the test functions
          call modg_project_numFlux_diffTest_Q(                             &
            &    nScalars      = nScalars,                                  &
            &    faceFlux      = facedata%faceFlux(iDir)                    &
            &                            %dat(:,:,1+nScalars:2*nScalars,:), &
            &    maxPolyDegree = scheme%maxPolyDegree,                      &
            &    length        = mesh%length,                               &
            &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid),     &
            &    dl_prod       = dl_prodDiff,                               &
            &    projection    = kerneldata%state_der,                      &
            &    dirVec        = dirVec(:,iDir)                             )
        end do

        case(P_space)
          allocate(faceFlux_Q(3))
          allocate(state_der_Q(nElems_fluid,(scheme%maxPolyDegree+1)**3, &
            &                    nScalars))
          do iDir = 1,3
            nFaces = size(facedata%faceflux(iDir)%dat,1)
            allocate(faceFlux_Q(iDir)%dat(nFaces,(scheme%maxPolyDegree+1)**2, &
              &                           2*nScalars,2))
            do leftOrRight = 1,2
              call ply_change_poly_space( inspace    = P_space,                 &
                &                         instate    = facedata%faceflux(iDir)  &
                &                                      %dat(:,:,:,leftOrRight), &
                &                         outstate   = faceFlux_Q(iDir)         &
                &                                      %dat(:,:,:,leftORRight), &
                &                         maxPolyDeg = scheme%maxPolyDegree,    &
                &                         nElems     = nElems_fluid,            &
                &                         nVars      = 2*nScalars,              &
                &                         nDims      = 2                        )
            end do
          end do

          call ply_change_poly_space( inspace    = P_space,              &
            &                         instate    = kerneldata%state_der, &
            &                         outstate   = state_der_Q,          &
            &                         maxPolyDeg = scheme%maxPolyDegree, &
            &                         nElems     = nElems_fluid,         &
            &                         nVars      = nScalars,             &
            &                         nDims      = 3                     )
          ! Projection of the appropriate  numerical flux onto the derivative
          ! of the test functions
          do iDir = 1,3
            call modg_project_numFlux_diffTest_Q(                         &
              &    nScalars      = nScalars,                              &
              &    faceFlux      = faceFlux_Q(iDir)                       &
              &                    %dat(:,:,1+nScalars:2*nScalars,:),     &
              &    maxPolyDegree = scheme%maxPolyDegree,                  &
              &    length        = mesh%length,                           &
              &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid), &
              &    dl_prod       = dl_prodDiff,                           &
              &    projection    = state_der_Q,                           &
              &    dirVec        = dirVec(:,iDir)                         )
          end do

          do iDir = 1,3
            do leftOrRight = 1,2
              call ply_change_poly_space( inspace    = Q_space,                 &
                &                         instate    = faceFlux_Q(iDir)         &
                &                                      %dat(:,:,:,leftORRight), &
                &                         outstate   = facedata%faceflux(iDir)  &
                &                                      %dat(:,:,:,leftOrRight), &
                &                         maxPolyDeg = scheme%maxPolyDegree,    &
                &                         nElems     = nElems_fluid,            &
                &                         nVars      = 2*nScalars,              &
                &                         nDims      = 2                        )
            end do
          end do
          call ply_change_poly_space( inspace    = Q_space,              &
            &                         instate    = state_der_Q,          &
            &                         outstate   = kerneldata%state_der, &
            &                         maxPolyDeg = scheme%maxPolyDegree, &
            &                         nElems     = nElems_fluid,         &
            &                         nVars      = nScalars,             &
            &                         nDims      = 3                     )

          deallocate(faceFlux_Q)
          deallocate(state_der_Q)

        end select

      case('navier_stokes', 'filtered_navier_stokes')
        ! Projection of the numerical flux
        ! ... x faces
        select case(scheme%basisType)
        case(Q_space)
          call modg_project_stabViscNumFluxX_Q(                         &
            &    numFlux       = facedata%faceFlux(1)%dat,              &
            &    faceState     = facedata%faceRep(1)%dat,               &
            &    equation      = equation,                              &
            &    maxPolyDegree = scheme%maxPolyDegree,                  &
            &    length        = mesh%length,                           &
            &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid), &
            &    projection    = kerneldata%state_der,                  &
            &    poly_proj     = poly_proj                              )
          ! ... y faces
          call modg_project_stabViscNumFluxY_Q(                         &
            &    numFlux       = facedata%faceFlux(2)%dat,              &
            &    faceState     = facedata%faceRep(2)%dat,               &
            &    equation      = equation,                              &
            &    maxPolyDegree = scheme%maxPolyDegree,                  &
            &    length        = mesh%length,                           &
            &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid), &
            &    projection    = kerneldata%state_der,                  &
            &    poly_proj     = poly_proj                              )

          ! ... z faces
          call modg_project_stabViscNumFluxZ_Q(                         &
            &    numFlux       = facedata%faceFlux(3)%dat,              &
            &    faceState     = facedata%faceRep(3)%dat,               &
            &    equation      = equation,                              &
            &    maxPolyDegree = scheme%maxPolyDegree,                  &
            &    length        = mesh%length,                           &
            &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid), &
            &    projection    = kerneldata%state_der,                  &
            &    poly_proj     = poly_proj                              )
        case(P_space)

          ! For P_space we need to change the polynomial space if we
          ! want to use the same routines that are used for the
          ! Q_space case.

          allocate(faceFlux_Q(3))
          allocate(faceRep_Q(3))
          allocate(state_der_Q(nElems_fluid,(scheme%maxPolyDegree+1)**3, &
            &                    nScalars))

          do iDir = 1,3
            nFaces = size(facedata%faceflux(iDir)%dat,1)
            allocate(faceFlux_Q(iDir)%dat(nFaces,(scheme%maxPolyDegree+1)**2, &
              &                             2*nScalars,2))
            allocate(faceRep_Q(iDir)%dat(nFaces,(scheme%maxPolyDegree+1)**2, &
              &                            2*nScalars,2))
            do leftOrRight = 1,2
              call ply_change_poly_space( inspace    = P_space,                 &
                &                         instate    = facedata%faceflux(iDir)  &
                &                                      %dat(:,:,:,leftOrRight), &
                &                         outstate   = faceFlux_Q(iDir)         &
                &                                      %dat(:,:,:,leftORRight), &
                &                         maxPolyDeg = scheme%maxPolyDegree,    &
                &                         nElems     = nElems_fluid,            &
                &                         nVars      = 2*nScalars,              &
                &                         nDims      = 2                        )

              call ply_change_poly_space( inspace    = P_space,                 &
                &                         instate    = facedata%faceRep(iDir)   &
                &                                      %dat(:,:,:,leftOrRight), &
                &                         outstate   = faceRep_Q(iDir)          &
                &                                      %dat(:,:,:,leftOrRight), &
                &                         maxPolyDeg = scheme%maxPolyDegree,    &
                &                         nElems     = nElems_fluid,            &
                &                         nVars      = 2*nScalars,              &
                &                         nDims      = 2                        )
            end do
          end do

          call ply_change_poly_space( inspace    = P_space,              &
            &                         instate    = kerneldata%state_der, &
            &                         outstate   = state_der_Q,          &
            &                         maxPolyDeg = scheme%maxPolyDegree, &
            &                         nElems     = nElems_fluid,         &
            &                         nVars      = nScalars,             &
            &                         nDims      = 3                     )

          call modg_project_stabViscNumFluxX_Q(                         &
            &    numFlux       = faceFlux_Q(1)%dat,                     &
            &    faceState     = faceRep_Q(1)%dat,                      &
            &    equation      = equation,                              &
            &    maxPolyDegree = scheme%maxPolyDegree,                  &
            &    length        = mesh%length,                           &
            &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid), &
            &    projection    = state_der_Q,                           &
            &    poly_proj     = poly_proj                              )
          ! ... y faces
          call modg_project_stabViscNumFluxY_Q(                         &
            &    numFlux       = faceFlux_Q(2)%dat,                     &
            &    faceState     = faceRep_Q(2)%dat,                      &
            &    equation      = equation,                              &
            &    maxPolyDegree = scheme%maxPolyDegree,                  &
            &    length        = mesh%length,                           &
            &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid), &
            &    projection    = state_der_Q,                           &
            &    poly_proj     = poly_proj                              )

          ! ... z faces
          call modg_project_stabViscNumFluxZ_Q(                         &
            &    numFlux       = faceFlux_Q(3)%dat,                     &
            &    faceState     = faceRep_Q(3)%dat,                      &
            &    equation      = equation,                              &
            &    maxPolyDegree = scheme%maxPolyDegree,                  &
            &    length        = mesh%length,                           &
            &    nElems_fluid  = mesh%descriptor%elem%nElems(eT_fluid), &
            &    projection    = state_der_Q,                           &
            &    poly_proj     = poly_proj                              )

          do iDir = 1,3
            do leftOrRight = 1,2
              call ply_change_poly_space( inspace    = Q_space,                 &
                &                         instate    = faceFlux_Q(iDir)         &
                &                                      %dat(:,:,:,leftORRight), &
                &                         outstate   = facedata%faceflux(iDir)  &
                &                                      %dat(:,:,:,leftOrRight), &
                &                         maxPolyDeg = scheme%maxPolyDegree,    &
                &                         nElems     = nElems_fluid,            &
                &                         nVars      = 2*nScalars,              &
                &                         nDims      = 2                        )

              call ply_change_poly_space( inspace    = Q_space,                 &
                &                         instate    = faceRep_Q(iDir)          &
                &                                      %dat(:,:,:,leftOrRight), &
                &                         outstate   = facedata%faceRep(iDir)   &
                &                                      %dat(:,:,:,leftOrRight), &
                &                         maxPolyDeg = scheme%maxPolyDegree,    &
                &                         nElems     = nElems_fluid,            &
                &                         nVars      = 2*nScalars,              &
                &                         nDims      = 2                        )
            end do
          end do
          call ply_change_poly_space( inspace    = Q_space,              &
            &                         instate    = state_der_Q,          &
            &                         outstate   = kerneldata%state_der, &
            &                         maxPolyDeg = scheme%maxPolyDegree, &
            &                         nElems     = nElems_fluid,         &
            &                         nVars      = nScalars,             &
            &                         nDims      = 3                     )

          deallocate(faceRep_Q)
          deallocate(faceFlux_Q)
          deallocate(state_der_Q)
        end select

      case default
        write(logUnit(1),*) 'ERROR in atl_modg_project_NumFlux:'
        write(logUnit(1),*) 'Unknown equation',equation%eq_kind
        write(logUnit(1),*) 'projections of viscous stabilization,'
        write(logUnit(1),*) ' stopping ... '
        call tem_abort()
      end select
    end if

    !> TODO NA - maybe move this call out of this routine ?
    ! Projection of the penalization terms if not computed elsewhere.
    if (penalizationdata%isActive .and. usePenalization) then
      select case(scheme%basisType)
      case(Q_space)
       call modg_project_penalization_Q(                   &
         &    mesh             = mesh,                     &
         &    maxPolyDegree    = scheme%maxPolyDegree,     &
         &    penalizationdata = penalizationdata,         &
         &    kerneldata       = kerneldata                )
      case(P_space)
        write(logUnit(1),*) 'ERROR in atl_modg_project_NumFlux: '
        write(logUnit(1),*) 'Penelization not yet implemented for P_space'
        call tem_abort
      end select
    end if
  end subroutine atl_modg_project_NumFlux


  !> Projection of the numerical flux in x direction onto the testfunctions.
  subroutine modg_project_stabViscNumFluxX_Q(                       &
    &          numFlux, faceState, equation, maxPolyDegree, length, &
    &          nElems_fluid, projection, poly_proj                  )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    integer :: iVar, iElem, iPoint, testPos, iVP
    integer :: iTestX, iTestY, iTestZ
    integer :: iOversamp, iOrig
    integer :: iDof, iDofY, iDofX
    integer :: nScalars, nPoints, nPVars, nOversamp

    real(kind=rk) :: p_a_left(equation%varSys%nScalars)
    real(kind=rk) :: p_b_left(equation%varSys%nScalars)
    real(kind=rk) :: p_c_left(equation%varSys%nScalars)

    real(kind=rk) :: p_a_right(equation%varSys%nScalars)
    real(kind=rk) :: p_b_right(equation%varSys%nScalars)
    real(kind=rk) :: p_c_right(equation%varSys%nScalars)

    real(kind=rk), allocatable :: flux_left(:,:)
    real(kind=rk), allocatable :: pointVal_flux_left(:,:)
    real(kind=rk), allocatable :: pointVal_left(:,:)
    real(kind=rk), allocatable :: state_left(:,:)
    real(kind=rk), allocatable :: nodal_a_left(:,:)
    real(kind=rk), allocatable :: nodal_b_left(:,:)
    real(kind=rk), allocatable :: nodal_c_left(:,:)
    real(kind=rk), allocatable :: modalA_left(:,:)
    real(kind=rk), allocatable :: modalB_left(:,:)
    real(kind=rk), allocatable :: modalC_left(:,:)

    real(kind=rk), allocatable :: flux_right(:,:)
    real(kind=rk), allocatable :: pointVal_flux_right(:,:)
    real(kind=rk), allocatable :: pointVal_right(:,:)
    real(kind=rk), allocatable :: state_right(:,:)
    real(kind=rk), allocatable :: nodal_a_right(:,:)
    real(kind=rk), allocatable :: nodal_b_right(:,:)
    real(kind=rk), allocatable :: nodal_c_right(:,:)
    real(kind=rk), allocatable :: modalA_right(:,:)
    real(kind=rk), allocatable :: modalB_right(:,:)
    real(kind=rk), allocatable :: modalC_right(:,:)

    real(kind=rk) :: velocity_left(3)
    real(kind=rk) :: velocity_right(3)
    real(kind=rk) :: jacobiDet
    real(kind=rk) :: testX_val_left, testX_grad_val_left
    real(kind=rk) :: testX_val_right, testX_grad_val_right
    real(kind=rk) :: outerNormalLeft, outerNormalRight
    ! -------------------------------------------------------------------- !

    nScalars = equation%varSys%nScalars
    nPoints = poly_proj%body_2D%nquadpoints
    nOversamp = poly_proj%body_2D%oversamp_dofs

    nPVars = (poly_proj%min_Degree+1)**2*equation%varSys%nScalars

    allocate( flux_left(nOversamp, nScalars), flux_right(nOversamp, nScalars) )
    allocate( state_left(nOversamp, nScalars) )
    allocate( state_right(nOversamp, nScalars) )
    allocate( pointVal_flux_left(nPoints, nScalars) )
    allocate( pointVal_flux_right(nPoints, nScalars) )
    allocate( pointVal_left(nPoints, nScalars) )
    allocate( pointVal_right(nPoints, nScalars) )
    allocate( nodal_a_left(nPoints, nScalars) )
    allocate( nodal_b_left(nPoints, nScalars) )
    allocate( nodal_c_left(nPoints, nScalars) )
    allocate( nodal_a_right(nPoints, nScalars) )
    allocate( nodal_b_right(nPoints, nScalars) )
    allocate( nodal_c_right(nPoints, nScalars) )
    allocate( modalA_left(nOversamp,nScalars) )
    allocate( modalB_left(nOversamp, nScalars) )
    allocate( modalC_left(nOversamp,nScalars) )
    allocate( modalA_right(nOversamp,nScalars) )
    allocate( modalB_right(nOversamp, nScalars) )
    allocate( modalC_right(nOversamp, nScalars) )

    ! The 1D Jacobi determinant
    jacobiDet = 0.5_rk * length

    ! The outer unit normals scaled by jacobiDet as this is a common factor
    outerNormalLeft = jacobiDet * atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = jacobiDet * atl_elemfaceToNormal_prp(tem_right)


    elementLoop: do iElem = 1,nElems_fluid

      state_left = 0.0_rk
      state_right = 0.0_rk
      flux_left = 0.0_rk
      flux_right = 0.0_rk

      do iVP = 1,nPVars
        iVar = (iVP-1)/((poly_proj%min_degree+1)**2) + 1
        iDof = iVP - (iVar-1)*((poly_proj%min_degree+1)**2)
        iDofY = (iDof-1) / (poly_proj%min_Degree+1)
        iDofX = mod(iDof-1,poly_proj%min_degree+1)
        iOversamp = iDofY*(poly_proj%oversamp_degree+1) + iDofX + 1
        iOrig = iDofY*(poly_proj%maxPolyDegree+1) + iDofX + 1

        ! Get u
        ! ... for left face
        state_left(iOversamp,iVar) = faceState(iElem,iOrig,iVar,1)
        ! ... for right face
        state_right(iOversamp,iVar) = faceState(iElem,iOrig,iVar,2)

        ! Caluclate (u^* - u)
        ! ... for the left face
        flux_left(iOversamp,iVar) = numFlux(iElem,iOrig,iVar+nScalars,1)  &
          &                       - state_left(iOversamp,iVar)
        ! ... for the right face
        flux_right(iOversamp,iVar) = numFlux(iElem,iOrig,iVar+nScalars,2) &
          &                        - state_right(iOversamp,iVar)
      end do

      ! Transform u to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me = poly_proj,              &
        &                       dim = 2 ,                    &
        &                       nVars = nScalars,            &
        &                       nodal_data = pointVal_left,  &
        &                       modal_data = state_left      )
      ! ... for the right face
      call ply_poly_project_m2n(me = poly_proj,              &
        &                       dim = 2 ,                    &
        &                       nVars = nScalars,            &
        &                       nodal_data = pointVal_right, &
        &                       modal_data = state_right     )

      ! Transform (u^* - u) to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me = poly_proj,                  &
        &                       dim = 2 ,                        &
        &                       nVars = nScalars,                &
        &                       nodal_data = pointVal_flux_left, &
        &                       modal_data = flux_left           )
      ! ... for the right face
      call ply_poly_project_m2n(me = poly_proj,                   &
        &                       dim = 2 ,                         &
        &                       nVars = nScalars,                 &
        &                       nodal_data = pointVal_flux_right, &
        &                       modal_data = flux_right           )

      ! Loop over all the points
      pointLoop: do iPoint = 1, nPoints

        ! Caculate velocity at this point
        ! ... for the left face
        velocity_left(1:3) = pointVal_left(iPoint,2:4) &
          &                  / pointVal_left(iPoint,1)
        ! ... for the right face
        velocity_right(1:3) = pointVal_right(iPoint,2:4) &
          &                   / pointVal_right(iPoint,1)

        ! Build matrix-vector product of nu_11 and values at current point
        ! ... for the left face
        nodal_a_left(iPoint,:)                                 &
          &  = atl_mult_nu11_NavierStokes(                     &
          &      density   = pointVal_left(iPoint,1),          &
          &      velocity  = velocity_left,                    &
          &      totEnergy = pointVal_left(iPoint,5),          &
          &      inVec     = pointVal_flux_left(iPoint,:),     &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda,     &
          &      thermCond = equation%NavierStokes%therm_cond, &
          &      heatCap   = equation%euler%cv                 )

        ! ... for the right face
        nodal_a_right(iPoint,:)                                &
          &  = atl_mult_nu11_NavierStokes(                     &
          &      density   = pointVal_right(iPoint,1),         &
          &      velocity  = velocity_right,                   &
          &      totEnergy = pointVal_right(iPoint,5),         &
          &      inVec     = pointVal_flux_right(iPoint,:),    &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda,     &
          &      thermCond = equation%NavierStokes%therm_cond, &
          &      heatCap   = equation%euler%cv                 )

        ! Build matrix-vector product of nu_21 and values at current point
        ! ... for the left face
        nodal_b_left(iPoint,:)                                 &
          &  = atl_mult_nu21_NavierStokes(                     &
          &      density   = pointVal_left(iPoint,1),          &
          &      velocity  = velocity_left,                    &
          &      inVec     = pointVal_flux_left(iPoint,:),     &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda      )

        ! ... for the right face
        nodal_b_right(iPoint,:)                                &
          &  = atl_mult_nu21_NavierStokes(                     &
          &      density   = pointVal_right(iPoint,1),         &
          &      velocity  = velocity_right,                   &
          &      inVec     = pointVal_flux_right(iPoint,:),    &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda      )

        ! Build matrix-vector product of nu_31 and values at current point
        ! ... for the left face
        nodal_c_left(iPoint,:)                                 &
          &  = atl_mult_nu31_NavierStokes(                     &
          &      density   = pointVal_left(iPoint,1),          &
          &      velocity  = velocity_left,                    &
          &      inVec     = pointVal_flux_left(iPoint,:),     &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda      )

        ! ... for the right face
        nodal_c_right(iPoint,:)                                &
          &  = atl_mult_nu31_NavierStokes(                     &
          &      density   = pointVal_right(iPoint,1),         &
          &      velocity  = velocity_right,                   &
          &      inVec     = pointVal_flux_right(iPoint,:),    &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda      )

      end do pointLoop

      ! Transform nodal_a and nodal_b back to modal space for projections
      ! ... for the left face
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 2,            &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_a_left, &
        &                       modal_data = modalA_left   )
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 2,            &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_b_left, &
        &                       modal_data = modalB_left   )
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 2,            &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_c_left, &
        &                       modal_data = modalC_left   )

      ! ... for the right face
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 2,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_a_right, &
        &                       modal_data = modalA_right   )
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 2,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_b_right, &
        &                       modal_data = modalB_right   )
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 2,             &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_c_right, &
        &                       modal_data = modalC_right   )

      ! Project onto all the test functions
      xTestLoop: do iTestX = 1, maxPolyDegree+1
        ! Evaluate test function in x direction at face coordinate
        ! ... for the left face
        testX_val_left = ply_faceValLeftBndTest(iTestX)
        testX_grad_val_left = ply_faceValLeftBndTestGrad(iTestX)
        ! ... for the right face
        testX_val_right = ply_faceValRightBndTest(iTestX)
        testX_grad_val_right = ply_faceValRightBndTestGrad(iTestX)

        ! Now we compute the projections. Since the test functions
        ! change when iTestY > 2 we moved the first two projections of
        ! of the projection loop.

        ! Project onto the first test function in y direction
        iTestY = 1
        iTestZ = 1
        testpos = iTestX
        ! ... for the left face
        p_a_left = jacobiDet * testX_grad_val_left  &
          &        * ply_scalProdDualLeg(iTestY,iTestY) &
          &        * ply_scalProdDualLeg(iTestZ,iTestZ) &
          &        * modalA_left(iTestY+(iTestZ-1)*(maxPolyDegree+1),:)
        projection(iElem,testPos,:) = projection(iElem,testPos,:)     &
          &                         - outerNormalLeft * ( p_a_left(:) )

        ! ... for the right face
        p_a_right = jacobiDet * testX_grad_val_right              &
          &         * ply_scalProdDualLeg(iTestY,iTestY) &
          &         * ply_scalProdDualLeg(iTestZ,iTestZ) &
          &         * modalA_right(iTestY+(iTestZ-1)*(maxPolyDegree+1),:)
        projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
          &                         - outerNormalRight * ( p_a_right(:) )

        if (maxPolyDegree+1 > 1) then

          ! Continue iTestY=1...
          iTestZ = 2
          testPos = iTestX &
            &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
            &       * (maxPolyDegree+1)

          ! ... for the left face
          p_a_left = jacobiDet * testX_grad_val_left      &
            &        * ply_scalProdDualLeg(iTestY,iTestY) &
            &        * ply_scalProdDualLeg(iTestZ,iTestZ) &
            &        * modalA_left(iTestY+(iTestZ-1)*(maxPolyDegree+1),:)
          p_c_left = testX_val_left                             &
            &        * ply_scalProdDualLeg(iTestY,iTestY)       &
            &        * ply_scalProdDualLegDiff(iTestZ-1,iTestZ) &
            &        * modalC_left(iTestY+(iTestZ-2)*(maxPolyDegree+1),:)
          projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
            &                         - outerNormalLeft * ( p_a_left(:)   &
            &                                               + p_c_left(:) )

          ! ... for the right face
          p_a_right = jacobiDet * testX_grad_val_right &
            &         * ply_scalProdDualLeg(iTestY,iTestY) &
            &         * ply_scalProdDualLeg(iTestZ,iTestZ) &
            &         * modalA_right(iTestY+(iTestZ-1)*(maxPolyDegree+1),:)
          p_c_right = testX_val_right &
            &         * ply_scalProdDualLeg(iTestY,iTestY) &
            &         * ply_scalProdDualLegDiff(iTestZ-1,iTestZ) &
            &         * modalC_right(iTestY+(iTestZ-2)*(maxPolyDegree+1),:)
          projection(iElem,testPos,:) = projection(iElem,testPos,:)         &
            &                         - outerNormalRight * ( p_a_right(:)   &
            &                                                + p_c_right(:) )

          do iTestZ = 3, maxPolyDegree+1
            testPos = iTestX &
              &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
              &       * (maxPolyDegree+1)

            ! ... for the left face
            p_a_left = jacobiDet * testX_grad_val_left                        &
              &        * ( ply_scalProdDualLeg(iTestY,iTestY)                 &
              &            * ply_scalProdDualLeg(iTestZ,iTestZ)               &
              &            * modalA_left(iTestY+(iTestZ-1)*(maxPolyDegree+1), &
              &                          :)                                   &
              &            + ply_scalProdDualLeg(iTestY,iTestY)               &
              &              * ply_scalProdDualLeg(iTestZ-2,iTestZ)           &
              &              * modalA_left(iTestY                             &
              &                            +(iTestZ-3)*(maxPolyDegree+1),:)   &
              &          )
            p_c_left = testX_val_left                                     &
              &        * ( ply_scalProdDualLeg(iTestY,iTestY)             &
              &            * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)     &
              &            * modalC_left(iTestY                           &
              &                          +(iTestZ-2)*(maxPolyDegree+1),:) &
              &          )
            projection(iElem,testPos,:) = projection(iElem,testPos,:)     &
              &                         - outerNormalLeft * ( p_a_left(:) &
              &                                             + p_c_left(:) )

            ! ... for the right face
            p_a_right = jacobiDet * testX_grad_val_right                    &
              &         * ( ply_scalProdDualLeg(iTestY,iTestY)              &
              &             * ply_scalProdDualLeg(iTestZ,iTestZ)            &
              &             * modalA_right(iTestY                           &
              &                            +(iTestZ-1)*(maxPolyDegree+1),:) &
              &             + ply_scalProdDualLeg(iTestY,iTestY)            &
              &               * ply_scalProdDualLeg(iTestZ-2,iTestZ)        &
              &             * modalA_right(iTestY                           &
              &                            +(iTestZ-3)*(maxPolyDegree+1),:) &
              &           )
            p_c_right = testX_val_right                                     &
              &         * ( ply_scalProdDualLeg(iTestY,iTestY)              &
              &             * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)      &
              &             * modalC_right(iTestY                           &
              &                            +(iTestZ-2)*(maxPolyDegree+1),:) &
              &           )
            projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
              &                         - outerNormalRight * ( p_a_right(:) &
              &                                              + p_c_right(:) )
          end do

          ! Next iTestY iteration
          iTestY = 2

          iTestZ = 1
          testPos = iTestX &
            &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
            &       * (maxPolyDegree+1)

          ! ... for the left face
          p_a_left = jacobiDet * testX_grad_val_left      &
            &        * ply_scalProdDualLeg(iTestY,iTestY) &
            &        * ply_scalProdDualLeg(iTestZ,iTestZ) &
            &        * modalA_left(iTestY+(iTestZ-1)*(maxPolyDegree+1),:)

          p_b_left(:) = testX_val_left                             &
            &           * ply_scalProdDualLegDiff(iTestY-1,iTestY) &
            &           * ply_scalProdDualLeg(iTestZ,iTestZ)       &
            &           * modalB_left(iTestY-1+(iTestZ-1)*(maxPolyDegree+1),:)
          projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
            &                         - outerNormalLeft * ( p_a_left(:)   &
            &                                               + p_b_left(:) )

          ! ... for the right face
          p_a_right = jacobiDet * testX_grad_val_right     &
            &         * ply_scalProdDualLeg(iTestY,iTestY) &
            &         * ply_scalProdDualLeg(iTestZ,iTestZ) &
            &         * modalA_right(iTestY+(iTestZ-1)*(maxPolyDegree+1),:)
          p_b_right = testX_val_right                            &
            &         * ply_scalProdDualLegDiff(iTestY-1,iTestY) &
            &         * ply_scalProdDualLeg(iTestZ,iTestZ)       &
            &         * modalB_right(iTestY-1+(iTestZ-1)*(maxPolyDegree+1),:)
          projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
            &                         - outerNormalRight * ( p_a_right(:) &
            &                                              + p_b_right(:) )

          iTestZ = 2
          testPos = iTestX &
            &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
            &       * (maxPolyDegree+1)

          ! ... for the left face
          p_a_left = jacobiDet * testX_grad_val_left      &
            &        * ply_scalProdDualLeg(iTestY,iTestY) &
            &        * ply_scalProdDualLeg(iTestZ,iTestZ) &
            &        * modalA_left(iTestY+(iTestZ-1)*(maxPolyDegree+1),:)
          p_b_left = testX_val_left                             &
            &        * ply_scalProdDualLegDiff(iTestY-1,iTestY) &
            &        * ply_scalProdDualLeg(iTestZ,iTestZ)       &
            &        * modalB_left(iTestY-1+(iTestZ-1)*(maxPolyDegree+1),:)
          p_c_left = testX_val_left                             &
            &        * ply_scalProdDualLeg(iTestY,iTestY)       &
            &        * ply_scalProdDualLegDiff(iTestZ-1,iTestZ) &
            &        * modalC_left(iTestY+(iTestZ-2)*(maxPolyDegree+1),:)
          projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
            &                         - outerNormalLeft * ( p_a_left(:)   &
            &                                               + p_b_left(:) &
            &                                               + p_c_left(:) )

          ! ... for the right face
          p_a_right = jacobiDet * testX_grad_val_right     &
            &         * ply_scalProdDualLeg(iTestY,iTestY) &
            &         * ply_scalProdDualLeg(iTestZ,iTestZ) &
            &         * modalA_right(iTestY+(iTestZ-1)*(maxPolyDegree+1),:)
          p_b_right = testX_val_right                            &
            &         * ply_scalProdDualLegDiff(iTestY-1,iTestY) &
            &         * ply_scalProdDualLeg(iTestZ,iTestZ)       &
            &         * modalB_right(iTestY-1+(iTestZ-1)*(maxPolyDegree+1),:)
          p_c_right = testX_val_right                            &
            &         * ply_scalProdDualLeg(iTestY,iTestY)       &
            &         * ply_scalProdDualLegDiff(iTestZ-1,iTestZ) &
            &         * modalC_right(iTestY+(iTestZ-2)*(maxPolyDegree+1),:)
          projection(iElem,testPos,:) = projection(iElem,testPos,:)         &
            &                         - outerNormalRight * ( p_a_right(:)   &
            &                                                + p_b_right(:) &
            &                                                + p_c_right(:) )
          do iTestZ = 3, maxPolyDegree+1
            testPos = iTestX &
              &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
              &       * (maxPolyDegree+1)

            ! ... for the left face
            p_a_left = jacobiDet * testX_grad_val_left                      &
              &        * ( ply_scalProdDualLeg(iTestY,iTestY)               &
              &            * ply_scalProdDualLeg(iTestZ,iTestZ)             &
              &            * modalA_left(iTestY                             &
              &                          +(iTestZ-1)*(maxPolyDegree+1),:)   &
              &            + ply_scalProdDualLeg(iTestY,iTestY)             &
              &              * ply_scalProdDualLeg(iTestZ-2,iTestZ)         &
              &              * modalA_left(iTestY                           &
              &                            +(iTestZ-3)*(maxPolyDegree+1),:) &
              &          )
            p_b_left = testX_val_left                                       &
              &        * ( ply_scalProdDualLegDiff(iTestY-1,iTestY)         &
              &            * ply_scalProdDualLeg(iTestZ,iTestZ)             &
              &            * modalB_left(iTestY-1                           &
              &                          +(iTestZ-1)*(maxPolyDegree+1),:)   &
              &            + ply_scalProdDualLegDiff(iTestY-1,iTestY)       &
              &              * ply_scalProdDualLeg(iTestZ-2,iTestZ)         &
              &              * modalB_left(iTestY-1                         &
              &                            +(iTestZ-3)*(maxPolyDegree+1),:) &
              &          )
            p_c_left = testX_val_left                                     &
              &        * ( ply_scalProdDualLeg(iTestY,iTestY)             &
              &            * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)     &
              &            * modalC_left(iTestY                           &
              &                          +(iTestZ-2)*(maxPolyDegree+1),:) &
              &          )
            projection(iElem,testPos,:) = projection(iElem,testPos,:)     &
              &                         - outerNormalLeft * ( p_a_left(:) &
              &                                             + p_b_left(:) &
              &                                             + p_c_left(:) )

            ! ... for the right face
            p_a_right = jacobiDet * testX_grad_val_right                      &
              &         * ( ply_scalProdDualLeg(iTestY,iTestY)                &
              &             * ply_scalProdDualLeg(iTestZ,iTestZ)              &
              &             * modalA_right(iTestY                             &
              &                            +(iTestZ-1)*(maxPolyDegree+1),:)   &
              &             + ply_scalProdDualLeg(iTestY,iTestY)              &
              &               * ply_scalProdDualLeg(iTestZ-2,iTestZ)          &
              &               * modalA_right(iTestY                           &
              &                              +(iTestZ-3)*(maxPolyDegree+1),:) &
              &          )
            p_b_right = testX_val_right                                       &
              &         * ( ply_scalProdDualLegDiff(iTestY-1,iTestY)          &
              &             * ply_scalProdDualLeg(iTestZ,iTestZ)              &
              &             * modalB_right(iTestY-1                           &
              &                            +(iTestZ-1)*(maxPolyDegree+1),:)   &
              &             + ply_scalProdDualLegDiff(iTestY-1,iTestY)        &
              &               * ply_scalProdDualLeg(iTestZ-2,iTestZ)          &
              &               * modalB_right(iTestY-1                         &
              &                              +(iTestZ-3)*(maxPolyDegree+1),:) &
              &           )
            p_c_right = testX_val_right                                     &
              &         * ( ply_scalProdDualLeg(iTestY,iTestY)              &
              &             * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)      &
              &             * modalC_right(iTestY                           &
              &                            +(iTestZ-2)*(maxPolyDegree+1),:) &
              &           )
            projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
              &                         - outerNormalRight * ( p_a_right(:) &
              &                                              + p_b_right(:) &
              &                                              + p_c_right(:) )
          end do


          ! Now, we project onto all the other test functions in y direction
          yTestLoop: do iTestY = 3, maxPolyDegree+1
            iTestZ = 1
            testPos = iTestX                                        &
              &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
              &       * (maxPolyDegree+1)

            ! ... for the left face
            p_a_left = jacobiDet * testX_grad_val_left                      &
              &        * ( ply_scalProdDualLeg(iTestY,iTestY)               &
              &            * ply_scalProdDualLeg(iTestZ,iTestZ)             &
              &            * modalA_left(iTestY                             &
              &                          +(iTestZ-1)*(maxPolyDegree+1),:)   &
              &            + ply_scalProdDualLeg(iTestY-2,iTestY)           &
              &              * ply_scalProdDualLeg(iTestZ,iTestZ)           &
              &              * modalA_left(iTestY-2                         &
              &                            +(iTestZ-1)*(maxPolyDegree+1),:) &
              &          )
            p_b_left = testX_val_left                             &
              &        * ply_scalProdDualLegDiff(iTestY-1,iTestY) &
              &        * ply_scalProdDualLeg(iTestZ,iTestZ)       &
              &        * modalB_left(iTestY-1+(iTestZ-1)*(maxPolyDegree+1),:)
            projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
              &                         - outerNormalLeft * ( p_a_left(:)   &
              &                                               + p_b_left(:) )

            ! ... for the right face
            p_a_right = jacobiDet * testX_grad_val_right                      &
              &         * ( ply_scalProdDualLeg(iTestY,iTestY)                &
              &             * ply_scalProdDualLeg(iTestZ,iTestZ)              &
              &             * modalA_right(iTestY                             &
              &                            +(iTestZ-1)*(maxPolyDegree+1),:)   &
              &             + ply_scalProdDualLeg(iTestY-2,iTestY)            &
              &               * ply_scalProdDualLeg(iTestZ,iTestZ)            &
              &               * modalA_right(iTestY-2                         &
              &                              +(iTestZ-1)*(maxPolyDegree+1),:) &
              &           )
            p_b_right = testX_val_right                            &
              &         * ply_scalProdDualLegDiff(iTestY-1,iTestY) &
              &         * ply_scalProdDualLeg(iTestZ,iTestZ)       &
              &         * modalB_right(iTestY-1+(iTestZ-1)*(maxPolyDegree+1),:)
            projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
              &                         - outerNormalRight * ( p_a_right(:) &
              &                                              + p_b_right(:) )

            iTestZ = 2
            testPos = iTestX                                        &
              &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
              &       * (maxPolyDegree+1)

            ! ... for the left face
            p_a_left = jacobiDet * testX_grad_val_left                      &
              &        * ( ply_scalProdDualLeg(iTestY,iTestY)               &
              &            * ply_scalProdDualLeg(iTestZ,iTestZ)             &
              &            * modalA_left(iTestY                             &
              &                          +(iTestZ-1)*(maxPolyDegree+1),:)   &
              &            + ply_scalProdDualLeg(iTestY-2,iTestY)           &
              &              * ply_scalProdDualLeg(iTestZ,iTestZ)           &
              &              * modalA_left(iTestY-2                         &
              &                            +(iTestZ-1)*(maxPolyDegree+1),:) &
              &          )
            p_b_left = testX_val_left                             &
              &        * ply_scalProdDualLegDiff(iTestY-1,iTestY) &
              &        * ply_scalProdDualLeg(iTestZ,iTestZ)       &
              &        * modalB_left(iTestY-1+(iTestZ-1)*(maxPolyDegree+1),:)
            p_c_left = testX_val_left                                       &
              &        * ( ply_scalProdDualLeg(iTestY,iTestY)               &
              &            * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)       &
              &            * modalC_left(iTestY                             &
              &                          +(iTestZ-2)*(maxPolyDegree+1),:)   &
              &            + ply_scalProdDualLeg(iTestY-2,iTestY)           &
              &              * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)     &
              &              * modalC_left(iTestY-2                         &
              &                            +(iTestZ-2)*(maxPolyDegree+1),:) &
              &          )
            projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
              &                         - outerNormalLeft * ( p_a_left(:)   &
              &                                               + p_b_left(:) &
              &                                               + p_c_left(:) )

            ! ... for the right face
            p_a_right = jacobiDet * testX_grad_val_right                      &
              &         * ( ply_scalProdDualLeg(iTestY,iTestY)                &
              &             * ply_scalProdDualLeg(iTestZ,iTestZ)              &
              &             * modalA_right(iTestY                             &
              &                            +(iTestZ-1)*(maxPolyDegree+1),:)   &
              &             + ply_scalProdDualLeg(iTestY-2,iTestY)            &
              &               * ply_scalProdDualLeg(iTestZ,iTestZ)            &
              &               * modalA_right(iTestY-2                         &
              &                              +(iTestZ-1)*(maxPolyDegree+1),:) &
              &           )
            p_b_right = testX_val_right                            &
              &         * ply_scalProdDualLegDiff(iTestY-1,iTestY) &
              &         * ply_scalProdDualLeg(iTestZ,iTestZ)       &
              &         * modalB_right(iTestY-1+(iTestZ-1)*(maxPolyDegree+1),:)
            p_c_right = testX_val_right                                     &
              &         * ( ply_scalProdDualLeg(iTestY,iTestY)              &
              &             * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)      &
              &             * modalC_right(iTestY                           &
              &                            +(iTestZ-2)*(maxPolyDegree+1),:) &
              &             + ply_scalProdDualLeg(iTestY-2,iTestY)          &
              &               * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)    &
              &               * modalC_right(iTestY-2                       &
              &                              +(iTestZ-2)*(maxPolyDegree+1), &
              &                              :)                             &
              &           )
            projection(iElem,testPos,:) = projection(iElem,testPos,:)         &
              &                         - outerNormalRight * ( p_a_right(:)   &
              &                                                + p_b_right(:) &
              &                                                + p_c_right(:) )

            zTestLoop: do iTestZ = 3, maxPolyDegree+1
              testPos = iTestX                                        &
                &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
                &       * (maxPolyDegree+1)

              ! ... for the left face
              p_a_left = jacobiDet * testX_grad_val_left                      &
                &        * ( ply_scalProdDualLeg(iTestY,iTestY)               &
                &            * ply_scalProdDualLeg(iTestZ,iTestZ)             &
                &            * modalA_left(iTestY                             &
                &                          +(iTestZ-1)*(maxPolyDegree+1),:)   &
                &            + ply_scalProdDualLeg(iTestY,iTestY)             &
                &              * ply_scalProdDualLeg(iTestZ-2,iTestZ)         &
                &              * modalA_left(iTestY                           &
                &                            +(iTestZ-3)*(maxPolyDegree+1),:) &
                &            + ply_scalProdDualLeg(iTestY-2,iTestY)           &
                &              * ply_scalProdDualLeg(iTestZ-2,iTestZ)         &
                &              * modalA_left(iTestY-2                         &
                &                            +(iTestZ-3)*(maxPolyDegree+1),:) &
                &            + ply_scalProdDualLeg(iTestY-2,iTestY)           &
                &              * ply_scalProdDualLeg(iTestZ,iTestZ)           &
                &              * modalA_left(iTestY-2                         &
                &                            +(iTestZ-1)*(maxPolyDegree+1),:) &
                &         )
              p_b_left = testX_val_left                                       &
                &        * ( ply_scalProdDualLegDiff(iTestY-1,iTestY)         &
                &            * ply_scalProdDualLeg(iTestZ,iTestZ)             &
                &            * modalB_left(iTestY-1                           &
                &                          +(iTestZ-1)*(maxPolyDegree+1),:)   &
                &            + ply_scalProdDualLegDiff(iTestY-1,iTestY)       &
                &              * ply_scalProdDualLeg(iTestZ-2,iTestZ)         &
                &              * modalB_left(iTestY-1                         &
                &                            +(iTestZ-3)*(maxPolyDegree+1),:) &
                &          )
              p_c_left = testX_val_left                                       &
                &        * ( ply_scalProdDualLeg(iTestY,iTestY)               &
                &            * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)       &
                &            * modalC_left(iTestY                             &
                &                          +(iTestZ-2)*(maxPolyDegree+1),:)   &
                &            + ply_scalProdDualLeg(iTestY-2,iTestY)           &
                &              * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)     &
                &              * modalC_left(iTestY-2                         &
                &                            +(iTestZ-2)*(maxPolyDegree+1),:) &
                &          )
              projection(iElem,testPos,:) = projection(iElem,testPos,:)     &
                &                         - outerNormalLeft * ( p_a_left(:) &
                &                                             + p_b_left(:) &
                &                                             + p_c_left(:) )

              ! ... for the right face
              p_a_right = jacobiDet * testX_grad_val_right                     &
                &         * ( ply_scalProdDualLeg(iTestY,iTestY)               &
                &             * ply_scalProdDualLeg(iTestZ,iTestZ)             &
                &             * modalA_right(iTestY                            &
                &                            +(iTestZ-1)*(maxPolyDegree+1),:)  &
                &             + ply_scalProdDualLeg(iTestY,iTestY)             &
                &               * ply_scalProdDualLeg(iTestZ-2,iTestZ)         &
                &               * modalA_right(iTestY                          &
                &                              +(iTestZ-3)*(maxPolyDegree+1),  &
                &                              :)                              &
                &             + ply_scalProdDualLeg(iTestY-2,iTestY)           &
                &               * ply_scalProdDualLeg(iTestZ-2,iTestZ)         &
                &               * modalA_right(iTestY-2                        &
                &                             +(iTestZ-3)*(maxPolyDegree+1),:) &
                &             + ply_scalProdDualLeg(iTestY-2,iTestY)           &
                &               * ply_scalProdDualLeg(iTestZ,iTestZ)           &
                &               * modalA_right(iTestY-2                        &
                &                              +(iTestZ-1)*(maxPolyDegree+1),  &
                &                              :)                              &
                &           )
              p_b_right = testX_val_right                                     &
                &         * ( ply_scalProdDualLegDiff(iTestY-1,iTestY)        &
                &             * ply_scalProdDualLeg(iTestZ,iTestZ)            &
                &             * modalB_right(iTestY-1                         &
                &                            +(iTestZ-1)*(maxPolyDegree+1),:) &
                &             + ply_scalProdDualLegDiff(iTestY-1,iTestY)      &
                &             * ply_scalProdDualLeg(iTestZ-2,iTestZ)          &
                &             * modalB_right(iTestY-1                         &
                &                            +(iTestZ-3)*(maxPolyDegree+1),:) )
              p_c_right = testX_val_right                                     &
                &         * ( ply_scalProdDualLeg(iTestY,iTestY)              &
                &             * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)      &
                &             * modalC_right(iTestY                           &
                &                            +(iTestZ-2)*(maxPolyDegree+1),:) &
                &             + ply_scalProdDualLeg(iTestY-2,iTestY)          &
                &               * ply_scalProdDualLegDiff(iTestZ-1,iTestZ)    &
                &               * modalC_right(iTestY-2                       &
                &                              +(iTestZ-2)*(maxPolyDegree+1), &
                &                              :)                             &
                &           )
              projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
                &                         - outerNormalRight * ( p_a_right(:) &
                &                                              + p_b_right(:) &
                &                                              + p_c_right(:) )
            end do zTestLoop
          end do yTestLoop
        end if
      end do xTestLoop


    end do elementLoop


  end subroutine modg_project_stabViscNumFluxX_Q


  !> Projection of the numerical flux in y direction onto the testfunctions.
  subroutine modg_project_stabViscNumFluxY_Q(                       &
    &          numFlux, faceState, equation, maxPolyDegree, length, &
    &          nElems_fluid, projection, poly_proj                  )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    integer :: iVar, iElem, iPoint, testPos, iVP
    integer :: iTestX, iTestY, iTestZ
    integer :: iTestFace
    integer :: iOversamp, iOrig
    integer :: iDof, iDofY, iDofX
    integer :: nScalars, nPoints, nPVars, nOversamp

    real(kind=rk) :: p_a_left(equation%varSys%nScalars)
    real(kind=rk) :: p_b_left(equation%varSys%nScalars)
    real(kind=rk) :: p_c_left(equation%varSys%nScalars)

    real(kind=rk) :: p_a_right(equation%varSys%nScalars)
    real(kind=rk) :: p_b_right(equation%varSys%nScalars)
    real(kind=rk) :: p_c_right(equation%varSys%nScalars)

    real(kind=rk), allocatable :: flux_left(:,:)
    real(kind=rk), allocatable :: pointVal_flux_left(:,:)
    real(kind=rk), allocatable :: pointVal_left(:,:)
    real(kind=rk), allocatable :: state_left(:,:)
    real(kind=rk), allocatable :: nodal_a_left(:,:)
    real(kind=rk), allocatable :: nodal_b_left(:,:)
    real(kind=rk), allocatable :: nodal_c_left(:,:)
    real(kind=rk), allocatable :: modalA_left(:,:)
    real(kind=rk), allocatable :: modalB_left(:,:)
    real(kind=rk), allocatable :: modalC_left(:,:)

    real(kind=rk), allocatable :: flux_right(:,:)
    real(kind=rk), allocatable :: pointVal_flux_right(:,:)
    real(kind=rk), allocatable :: pointVal_right(:,:)
    real(kind=rk), allocatable :: state_right(:,:)
    real(kind=rk), allocatable :: nodal_a_right(:,:)
    real(kind=rk), allocatable :: nodal_b_right(:,:)
    real(kind=rk), allocatable :: nodal_c_right(:,:)
    real(kind=rk), allocatable :: modalA_right(:,:)
    real(kind=rk), allocatable :: modalB_right(:,:)
    real(kind=rk), allocatable :: modalC_right(:,:)

    real(kind=rk) :: velocity_left(3)
    real(kind=rk) :: velocity_right(3)
    real(kind=rk) :: jacobiDet
    real(kind=rk) :: testY_val_left, testY_grad_val_left
    real(kind=rk) :: testY_val_right, testY_grad_val_right
    real(kind=rk) :: outerNormalLeft, outerNormalRight
    real(kind=rk) :: legprod
    real(kind=rk) :: legprod_square
    real(kind=rk) :: tmp_left(equation%varSys%nScalars)
    real(kind=rk) :: tmp_right(equation%varSys%nScalars)
    ! -------------------------------------------------------------------- !

    nScalars = equation%varSys%nScalars
    nPoints = poly_proj%body_2D%nquadpoints
    nOversamp = poly_proj%body_2D%oversamp_dofs

    nPVars = (poly_proj%min_degree+1)**2*equation%varSys%nScalars

    allocate( flux_left(nOversamp, nScalars), flux_right(nOversamp, nScalars) )
    allocate( state_left(nOversamp, nScalars) )
    allocate( state_right(nOversamp, nScalars) )
    allocate( pointVal_flux_left(nPoints, nScalars) )
    allocate( pointVal_flux_right(nPoints, nScalars) )
    allocate( pointVal_left(nPoints, nScalars) )
    allocate( pointVal_right(nPoints, nScalars) )
    allocate( nodal_a_left(nPoints, nScalars) )
    allocate( nodal_b_left(nPoints, nScalars) )
    allocate( nodal_c_left(nPoints, nScalars) )
    allocate( nodal_a_right(nPoints, nScalars) )
    allocate( nodal_b_right(nPoints, nScalars) )
    allocate( nodal_c_right(nPoints, nScalars) )
    allocate( modalA_left(nOversamp,nScalars) )
    allocate( modalB_left(nOversamp, nScalars) )
    allocate( modalC_left(nOversamp,nScalars) )
    allocate( modalA_right(nOversamp,nScalars) )
    allocate( modalB_right(nOversamp, nScalars) )
    allocate( modalC_right(nOversamp, nScalars) )

    ! The 1D Jacobi determinant
    jacobiDet = 0.5_rk * length

    ! The outer unit normals scaled by jacobiDet as this is a common factor
    outerNormalLeft = jacobiDet * atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = jacobiDet * atl_elemfaceToNormal_prp(tem_right)


    elementLoop: do iElem = 1,nElems_fluid

      state_left = 0.0_rk
      state_right = 0.0_rk
      flux_left = 0.0_rk
      flux_right = 0.0_rk

      do iVP = 1,nPVars
        iVar = (iVP-1)/((poly_proj%min_degree+1)**2) + 1
        iDof = iVP - (iVar-1)*((poly_proj%min_degree+1)**2)
        iDofY = (iDof-1) / (poly_proj%min_Degree+1)
        iDofX = mod(iDof-1,poly_proj%min_degree+1)
        iOversamp = iDofY*(poly_proj%oversamp_degree+1) + iDofX + 1
        iOrig = iDofY*(poly_proj%maxPolyDegree+1) + iDofX + 1

        ! Get u
        ! ... for left face
        state_left(iOversamp,iVar) = faceState(iElem,iOrig,iVar,1)
        ! ... for right face
        state_right(iOversamp,iVar) = faceState(iElem,iOrig,iVar,2)

        ! Caluclate (u^* - u)
        ! ... for the left face
        flux_left(iOversamp,iVar) = numFlux(iElem,iOrig,iVar+nScalars,1)  &
          &                       - state_left(iOversamp,iVar)
        ! ... for the right face
        flux_right(iOversamp,iVar) = numFlux(iElem,iOrig,iVar+nScalars,2) &
          &                        - state_right(iOversamp,iVar)
      end do

      ! Transform  u to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me = poly_proj,             &
        &                       dim = 2 ,                   &
        &                       nVars = nScalars,           &
        &                       nodal_data = pointVal_left, &
        &                       modal_data = state_left     )
      ! ... for the right face
      call ply_poly_project_m2n(me = poly_proj,              &
        &                       dim = 2 ,                    &
        &                       nVars = nScalars,            &
        &                       nodal_data = pointVal_right, &
        &                       modal_data = state_right     )

      ! Transform (u^* - u) to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me = poly_proj,                  &
        &                       dim = 2 ,                        &
        &                       nVars = nScalars,                &
        &                       nodal_data = pointVal_flux_left, &
        &                       modal_data = flux_left           )
      ! ... for the right face
      call ply_poly_project_m2n(me = poly_proj,                   &
        &                       dim = 2 ,                         &
        &                       nVars = nScalars,                 &
        &                       nodal_data = pointVal_flux_right, &
        &                       modal_data = flux_right           )

      ! Loop over all the points
      pointLoop: do iPoint = 1, nPoints

        ! Caculate velocity at this point
        ! ... for the left face
        velocity_left = pointVal_left(iPoint,2:4) &
          &             / pointVal_left(iPoint,1)
        ! ... for the right face
        velocity_right = pointVal_right(iPoint,2:4) &
          &              / pointVal_right(iPoint,1)

        ! Build matrix-vector product of nu_12 and values at current point
        ! ... for the left face
        nodal_a_left(iPoint,:)                            &
          &  = atl_mult_nu12_NavierStokes(                &
          &      density  = pointVal_left(iPoint,1),      &
          &      velocity = velocity_left,                &
          &      inVec    = pointVal_flux_left(iPoint,:), &
          &      mu       = equation%NavierStokes%mu,     &
          &      lambda   = equation%NavierStokes%lambda  )

        ! ... for the right face
        nodal_a_right(iPoint,:)                            &
          &  = atl_mult_nu12_NavierStokes(                 &
          &      density  = pointVal_right(iPoint,1),      &
          &      velocity = velocity_right,                &
          &      inVec    = pointVal_flux_right(iPoint,:), &
          &      mu       = equation%NavierStokes%mu,      &
          &      lambda   = equation%NavierStokes%lambda   )

        ! Build matrix-vector product of nu_22 and values at current point
        ! ... for the left face
        nodal_b_left(iPoint,:) &
          &  = atl_mult_nu22_NavierStokes(                     &
          &      density   = pointVal_left(iPoint,1),          &
          &      velocity  = velocity_left,                    &
          &      totEnergy = pointVal_left(iPoint,5),          &
          &      inVec     = pointVal_flux_left(iPoint,:),     &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda,     &
          &      thermCond = equation%NavierStokes%therm_cond, &
          &      heatCap   = equation%euler%cv                 )

        ! ... for the right face
        nodal_b_right(iPoint,:) &
          &  = atl_mult_nu22_NavierStokes(                     &
          &      density   = pointVal_right(iPoint,1),         &
          &      velocity  = velocity_right,                   &
          &      totEnergy = pointVal_right(iPoint,5),         &
          &      inVec     = pointVal_flux_right(iPoint,:),    &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda,     &
          &      thermCond = equation%NavierStokes%therm_cond, &
          &      heatCap   = equation%euler%cv                 )

        ! Build matrix-vector product of nu_22 and values at current point
        ! ... for the left face
        nodal_c_left(iPoint,:) &
          &  = atl_mult_nu32_NavierStokes(                &
          &      density  = pointVal_left(iPoint,1),      &
          &      velocity = velocity_left,                &
          &      inVec    = pointVal_flux_left(iPoint,:), &
          &      mu       = equation%NavierStokes%mu,     &
          &      lambda   = equation%NavierStokes%lambda  )

        ! ... for the right face
        nodal_c_right(iPoint,:) &
          &  = atl_mult_nu32_NavierStokes(                 &
          &      density  = pointVal_right(iPoint,1),      &
          &      velocity = velocity_right,                &
          &      inVec    = pointVal_flux_right(iPoint,:), &
          &      mu       = equation%NavierStokes%mu,      &
          &      lambda   = equation%NavierStokes%lambda   )

      end do pointLoop

      ! Transform nodal_a and nodal_b back to modal space for projections
      ! ... for the left face
      call ply_poly_project_n2m(me = poly_proj,            &
        &                       dim = 2 ,                  &
        &                       nVars = nScalars,          &
        &                       nodal_data = nodal_a_left, &
        &                       modal_data = modalA_left   )
      call ply_poly_project_n2m(me = poly_proj,            &
        &                       dim = 2 ,                  &
        &                       nVars = nScalars,          &
        &                       nodal_data = nodal_b_left, &
        &                       modal_data = modalB_left   )
      call ply_poly_project_n2m(me = poly_proj,            &
        &                       dim = 2 ,                  &
        &                       nVars = nScalars,          &
        &                       nodal_data = nodal_c_left, &
        &                       modal_data = modalC_left   )

      ! ... for the right face
      call ply_poly_project_n2m(me = poly_proj,             &
        &                       dim = 2 ,                   &
        &                       nVars = nScalars,           &
        &                       nodal_data = nodal_a_right, &
        &                       modal_data = modalA_right   )
      call ply_poly_project_n2m(me = poly_proj,             &
        &                       dim = 2 ,                   &
        &                       nVars = nScalars,           &
        &                       nodal_data = nodal_b_right, &
        &                       modal_data = modalB_right   )
      call ply_poly_project_n2m(me = poly_proj,             &
        &                       dim = 2 ,                   &
        &                       nVars = nScalars,           &
        &                       nodal_data = nodal_c_right, &
        &                       modal_data = modalC_right   )

      ! Project onto all the test functions
      yTestLoop: do iTestY = 1, maxPolyDegree+1
        ! Evaluate test function in y direction at face coordinate
        ! ... for the left face
        testY_val_left = ply_faceValLeftBndTest(iTestY)
        testY_grad_val_left = ply_faceValLeftBndTestGrad(iTestY)
        ! ... for the right face
        testY_val_right = ply_faceValRightBndTest(iTestY)
        testY_grad_val_right = ply_faceValRightBndTestGrad(iTestY)

        ! Now we compute the projections. Since the test functions
        ! change when iTestY > 2 we moved the first two projections of
        ! of the projection loop.

        xTestLoop: do iTestX = 1, maxPolyDegree+1
          zTestLoop: do iTestZ = 1, maxPolyDegree+1

            testPos = iTestX                                        &
              &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
              &       * (maxPolyDegree+1)

            legprod_square = ply_scalProdDualLeg(iTestZ, iTestZ)
            legprod = ply_scalProdDualLeg(iTestX, iTestX) &
              &       * legprod_square
            iTestFace = iTestX+(iTestZ-1)*(maxPolyDegree+1)

            tmp_left = legprod * modalB_left(iTestFace,:)
            tmp_right = legprod * modalB_right(iTestFace,:)

            if (iTestX > 2) then
              legprod = ply_scalProdDualLeg(iTestX-2, iTestX) &
                &       * legprod_square
              iTestFace = iTestX+(iTestZ-1)*(maxPolyDegree+1) - 2

              tmp_left = tmp_left + legprod * modalB_left(iTestFace,:)
              tmp_right = tmp_right + legprod * modalB_right(iTestFace,:)
            end if
            if (iTestZ > 2) then
              legprod = ply_scalProdDualLeg(iTestX, iTestX) &
                &       * ply_scalProdDualLeg(iTestZ-2, iTestZ)
              iTestFace = iTestX+(iTestZ-3)*(maxPolyDegree+1)

              tmp_left = tmp_left &
                &      + legprod * modalB_left(iTestFace,:)
              tmp_right = tmp_right &
                &       + legprod * modalB_right(iTestFace,:)
              if (iTestX > 2) then
                legprod = ply_scalProdDualLeg(iTestX-2,iTestX) &
                   &      * ply_scalProdDualLeg(iTestZ-2,iTestZ)
                iTestFace = iTestX+(iTestZ-3)*(maxPolyDegree+1) - 2

                tmp_left = tmp_left &
                  &      + legprod * modalB_left(iTestFace,:)
                tmp_right = tmp_right &
                  &       + legprod * modalB_right(iTestFace,:)
              end if
            end if

            p_b_left = jacobiDet * testY_grad_val_left * tmp_left
            p_b_right = jacobiDet * testY_grad_val_right * tmp_right

            if (iTestX > 1) then
              legprod = ply_scalProdDualLegDiff(iTestX-1,iTestX) &
                &       * legprod_square
              iTestFace = iTestX+(iTestZ-1)*(maxPolyDegree+1) - 1

              tmp_left = legprod * modalA_left(iTestFace,:)
              tmp_right = legprod * modalA_right(iTestFace,:)

              if (iTestZ > 2) then
                legprod = ply_scalProdDualLegDiff(iTestX-1, iTestX) &
                  &       * ply_scalProdDualLeg(iTestZ-2,iTestZ)
                iTestFace = iTestX+(iTestZ-3)*(maxPolyDegree+1) - 1

                tmp_left = tmp_left + legprod * modalA_left(iTestFace,:)
                tmp_right = tmp_right + legprod * modalA_right(iTestFace,:)
              end if

              p_a_left = testY_val_left * tmp_left
              p_a_right = testY_val_right * tmp_right

            else

              p_a_left = 0.0_rk
              p_a_right = 0.0_rk

            end if

            if (iTestZ > 1) then
              legprod = ply_scalProdDualLeg(iTestX, iTestX) &
                &       * ply_scalProdDualLegDiff(iTestZ-1, iTestZ)
              iTestFace = iTestX+(iTestZ-2)*(maxPolyDegree+1)

              tmp_left = legprod * modalC_left(iTestFace,:)
              tmp_right = legprod * modalC_right(iTestFace,:)

              if (iTestX > 2) then
                legprod = ply_scalProdDualLeg(iTestX-2, iTestX) &
                  &       * ply_scalProdDualLegDiff(iTestZ-1, iTestZ)
                iTestFace = iTestX+(iTestZ-2)*(maxPolyDegree+1) - 2

                tmp_left = tmp_left + legprod * modalC_left(iTestFace,:)
                tmp_right = tmp_right + legprod * modalC_right(iTestFace,:)
              end if

              p_c_left = testY_val_left * tmp_left
              p_c_right = testY_val_right * tmp_right

            else

              p_c_left = 0.0_rk
              p_c_right = 0.0_rk

            end if

            projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
              &                         - outerNormalLeft * ( p_a_left(:)   &
              &                                             + p_b_left(:)   &
              &                                             + p_c_left(:) ) &
              &                         - outerNormalRight * ( p_a_right(:) &
              &                                              + p_b_right(:) &
              &                                              + p_c_right(:) )

          end do zTestLoop
        end do xTestLoop
      end do yTestLoop


    end do elementLoop


  end subroutine modg_project_stabViscNumFluxY_Q


  !> Projection of the numerical flux in y direction onto the testfunctions.
  subroutine modg_project_stabViscNumFluxZ_Q(                       &
    &          numFlux, faceState, equation, maxPolyDegree, length, &
    &          nElems_fluid, projection, poly_proj                  )
    ! -------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    integer :: iVar, iElem, iPoint, testPos, iVP
    integer :: iTestX, iTestY, iTestZ
    integer :: iTestFace
    integer :: iOversamp, iOrig
    integer :: iDof, iDofY, iDofX
    integer :: nScalars, nPoints, nPVars, nOversamp

    real(kind=rk) :: p_a_left(equation%varSys%nScalars)
    real(kind=rk) :: p_b_left(equation%varSys%nScalars)
    real(kind=rk) :: p_c_left(equation%varSys%nScalars)

    real(kind=rk) :: p_a_right(equation%varSys%nScalars)
    real(kind=rk) :: p_b_right(equation%varSys%nScalars)
    real(kind=rk) :: p_c_right(equation%varSys%nScalars)

    real(kind=rk), allocatable :: flux_left(:,:)
    real(kind=rk), allocatable :: pointVal_flux_left(:,:)
    real(kind=rk), allocatable :: pointVal_left(:,:)
    real(kind=rk), allocatable :: state_left(:,:)
    real(kind=rk), allocatable :: nodal_a_left(:,:)
    real(kind=rk), allocatable :: nodal_b_left(:,:)
    real(kind=rk), allocatable :: nodal_c_left(:,:)
    real(kind=rk), allocatable :: modalA_left(:,:)
    real(kind=rk), allocatable :: modalB_left(:,:)
    real(kind=rk), allocatable :: modalC_left(:,:)

    real(kind=rk), allocatable :: flux_right(:,:)
    real(kind=rk), allocatable :: pointVal_flux_right(:,:)
    real(kind=rk), allocatable :: pointVal_right(:,:)
    real(kind=rk), allocatable :: state_right(:,:)
    real(kind=rk), allocatable :: nodal_a_right(:,:)
    real(kind=rk), allocatable :: nodal_b_right(:,:)
    real(kind=rk), allocatable :: nodal_c_right(:,:)
    real(kind=rk), allocatable :: modalA_right(:,:)
    real(kind=rk), allocatable :: modalB_right(:,:)
    real(kind=rk), allocatable :: modalC_right(:,:)

    real(kind=rk) :: velocity_left(3)
    real(kind=rk) :: velocity_right(3)
    real(kind=rk) :: jacobiDet
    real(kind=rk) :: testZ_val_left, testZ_grad_val_left
    real(kind=rk) :: testZ_val_right, testZ_grad_val_right
    real(kind=rk) :: outerNormalLeft, outerNormalRight
    real(kind=rk) :: legprod
    real(kind=rk) :: legprod_square
    real(kind=rk) :: tmp_left(equation%varSys%nScalars)
    real(kind=rk) :: tmp_right(equation%varSys%nScalars)
    ! -------------------------------------------------------------------- !

    nScalars = equation%varSys%nScalars
    nPoints = poly_proj%body_2D%nquadpoints
    nOversamp = poly_proj%body_2D%oversamp_dofs

    nPVars = (poly_proj%min_degree+1)**2*equation%varSys%nScalars

    allocate( flux_left(nOversamp, nScalars), flux_right(nOversamp, nScalars) )
    allocate( state_left(nOversamp, nScalars) )
    allocate( state_right(nOversamp, nScalars) )
    allocate( pointVal_flux_left(nPoints, nScalars) )
    allocate( pointVal_flux_right(nPoints, nScalars) )
    allocate( pointVal_left(nPoints, nScalars) )
    allocate( pointVal_right(nPoints, nScalars) )
    allocate( nodal_a_left(nPoints, nScalars) )
    allocate( nodal_b_left(nPoints, nScalars) )
    allocate( nodal_c_left(nPoints, nScalars) )
    allocate( nodal_a_right(nPoints, nScalars) )
    allocate( nodal_b_right(nPoints, nScalars) )
    allocate( nodal_c_right(nPoints, nScalars) )
    allocate( modalA_left(nOversamp,nScalars) )
    allocate( modalB_left(nOversamp, nScalars) )
    allocate( modalC_left(nOversamp,nScalars) )
    allocate( modalA_right(nOversamp,nScalars) )
    allocate( modalB_right(nOversamp, nScalars) )
    allocate( modalC_right(nOversamp, nScalars) )

    ! The 1D Jacobi determinant
    jacobiDet = 0.5_rk * length

    ! The outer unit normals scaled by jacobiDet as this is a common factor
    outerNormalLeft = jacobiDet * atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = jacobiDet * atl_elemfaceToNormal_prp(tem_right)


    elementLoop: do iElem = 1,nElems_fluid

      state_left = 0.0_rk
      state_right = 0.0_rk
      flux_left = 0.0_rk
      flux_right = 0.0_rk

      do iVP = 1,nPVars
        iVar = (iVP-1)/((poly_proj%min_degree+1)**2) + 1
        iDof = iVP - (iVar-1)*((poly_proj%min_degree+1)**2)
        iDofY = (iDof-1) / (poly_proj%min_Degree+1)
        iDofX = mod(iDof-1,poly_proj%min_degree+1)
        iOversamp = iDofY*(poly_proj%oversamp_degree+1) + iDofX + 1
        iOrig = iDofY*(poly_proj%maxPolyDegree+1) + iDofX + 1

        ! Get u
        ! ... for left face
        state_left(iOversamp,iVar) = faceState(iElem,iOrig,iVar,1)
        ! ... for right face
        state_right(iOversamp,iVar) = faceState(iElem,iOrig,iVar,2)

        ! Caluclate (u^* - u)
        ! ... for the left face
        flux_left(iOversamp,iVar) = numFlux(iElem,iOrig,iVar+nScalars,1)  &
          &                       - state_left(iOversamp,iVar)
        ! ... for the right face
        flux_right(iOversamp,iVar) = numFlux(iElem,iOrig,iVar+nScalars,2) &
          &                        - state_right(iOversamp,iVar)
      end do

      ! Transform  u to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me = poly_proj,             &
        &                       dim = 2 ,                   &
        &                       nVars = nScalars,           &
        &                       nodal_data = pointVal_left, &
        &                       modal_data = state_left     )
      ! ... for the right face
      call ply_poly_project_m2n(me = poly_proj,              &
        &                       dim = 2 ,                    &
        &                       nVars = nScalars,            &
        &                       nodal_data = pointVal_right, &
        &                       modal_data = state_right     )

      ! Transform (u^* - u) to nodal representation
      ! ... for the left face
      call ply_poly_project_m2n(me = poly_proj,                  &
        &                       dim = 2 ,                        &
        &                       nVars = nScalars,                &
        &                       nodal_data = pointVal_flux_left, &
        &                       modal_data = flux_left           )
      ! ... for the right face
      call ply_poly_project_m2n(me = poly_proj,                   &
        &                       dim = 2 ,                         &
        &                       nVars = nScalars,                 &
        &                       nodal_data = pointVal_flux_right, &
        &                       modal_data = flux_right           )

      ! Loop over all the points
      pointLoop: do iPoint = 1, nPoints

        ! Caculate velocity at this point
        ! ... for the left face
        velocity_left(1:3) = pointVal_left(iPoint,2:4) &
          &                  / pointVal_left(iPoint,1)
        ! ... for the right face
        velocity_right(1:3) = pointVal_right(iPoint,2:4) &
          &                   / pointVal_right(iPoint,1)

        ! Build matrix-vector product of nu_13 and values at current point
        ! ... for the left face
        nodal_a_left(iPoint,:)                             &
          &  = atl_mult_nu13_NavierStokes(                 &
          &      density   = pointVal_left(iPoint,1),      &
          &      velocity  = velocity_left,                &
          &      inVec     = pointVal_flux_left(iPoint,:), &
          &      mu        = equation%NavierStokes%mu,     &
          &      lambda    = equation%NavierStokes%lambda  )

        ! ... for the right face
        nodal_a_right(iPoint,:)                             &
          &  = atl_mult_nu13_NavierStokes(                  &
          &      density   = pointVal_right(iPoint,1),      &
          &      velocity  = velocity_right,                &
          &      inVec     = pointVal_flux_right(iPoint,:), &
          &      mu        = equation%NavierStokes%mu,      &
          &      lambda    = equation%NavierStokes%lambda   )

        ! Build matrix-vector product of nu_23 and values at current point
        ! ... for the left face
        nodal_b_left(iPoint,:)                             &
          &  = atl_mult_nu23_NavierStokes(                 &
          &      density   = pointVal_left(iPoint,1),      &
          &      velocity  = velocity_left,                &
          &      inVec     = pointVal_flux_left(iPoint,:), &
          &      mu        = equation%NavierStokes%mu,     &
          &      lambda    = equation%NavierStokes%lambda  )

        ! ... for the right face
        nodal_b_right(iPoint,:)                             &
          &  = atl_mult_nu23_NavierStokes(                  &
          &      density   = pointVal_right(iPoint,1),      &
          &      velocity  = velocity_right,                &
          &      inVec     = pointVal_flux_right(iPoint,:), &
          &      mu        = equation%NavierStokes%mu,      &
          &      lambda    = equation%NavierStokes%lambda   )

        ! Build matrix-vector product of nu_33 and values at current point
        ! ... for the left face
        nodal_c_left(iPoint,:)                                 &
          &  = atl_mult_nu33_NavierStokes(                     &
          &      density   = pointVal_left(iPoint,1),          &
          &      velocity  = velocity_left,                    &
          &      totEnergy = pointVal_left(iPoint,5),          &
          &      inVec     = pointVal_flux_left(iPoint,:),     &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda,     &
          &      thermCond = equation%NavierStokes%therm_cond, &
          &      heatCap   = equation%euler%cv                 )

        ! ... for the right face
        nodal_c_right(iPoint,:)                                &
          &  = atl_mult_nu33_NavierStokes(                     &
          &      density   = pointVal_right(iPoint,1),         &
          &      velocity  = velocity_right,                   &
          &      totEnergy = pointVal_right(iPoint,5),         &
          &      inVec     = pointVal_flux_right(iPoint,:),    &
          &      mu        = equation%NavierStokes%mu,         &
          &      lambda    = equation%NavierStokes%lambda,     &
          &      thermCond = equation%NavierStokes%therm_cond, &
          &      heatCap   = equation%euler%cv                 )

      end do pointLoop


      ! Transform nodal_a and nodal_b back to modal space for projections
      ! ... for the left face
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 2 ,           &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_a_left, &
        &                       modal_data = modalA_left   )
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 2 ,           &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_b_left, &
        &                       modal_data = modalB_left   )
      call ply_poly_project_n2m(me         = poly_proj,    &
        &                       dim        = 2 ,           &
        &                       nVars      = nScalars,     &
        &                       nodal_data = nodal_c_left, &
        &                       modal_data = modalC_left   )

      ! ... for the right face
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 2 ,            &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_a_right, &
        &                       modal_data = modalA_right   )
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 2 ,            &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_b_right, &
        &                       modal_data = modalB_right   )
      call ply_poly_project_n2m(me         = poly_proj,     &
        &                       dim        = 2 ,            &
        &                       nVars      = nScalars,      &
        &                       nodal_data = nodal_c_right, &
        &                       modal_data = modalC_right   )

      ! Project onto all the test functions
      zTestLoop: do iTestZ = 1, maxPolyDegree+1
        ! Evaluate test function in x direction at face coordinate
        ! ... for the left face
        testZ_val_left = ply_faceValLeftBndTest(iTestZ)
        testZ_grad_val_left = ply_faceValLeftBndTestGrad(iTestZ)
        ! ... for the right face
        testZ_val_right = ply_faceValRightBndTest(iTestZ)
        testZ_grad_val_right = ply_faceValRightBndTestGrad(iTestZ)

        xTestLoop: do iTestX = 1, maxPolyDegree+1
          yTestLoop: do iTestY = 1, maxPolyDegree+1

            testPos = iTestX &
              &     + ( (iTestY-1) + (iTestZ-1)*(maxPolyDegree+1) ) &
              &       * (maxPolyDegree+1)

            legprod_square = ply_scalProdDualLeg(iTestY,iTestY)
            legprod = ply_scalProdDualLeg(iTestX,iTestX) &
              &       * legprod_square
            iTestFace = iTestX+(iTestY-1)*(maxPolyDegree+1)

            tmp_left = legprod * modalC_left(iTestFace,:)
            tmp_right = legprod * modalC_right(iTestFace,:)

            if (iTestX > 2) then
              legprod = ply_scalProdDualLeg(iTestX-2,iTestX) &
                &       * legprod_square
              iTestFace = iTestX+(iTestY-1)*(maxPolyDegree+1) - 2

              tmp_left = tmp_left + legprod * modalC_left(iTestFace,:)
              tmp_right = tmp_right + legprod * modalC_right(iTestFace,:)
            end if

            if (iTestY > 2) then
              legprod = ply_scalProdDualLeg(iTestX,iTestX) &
                &       * ply_scalProdDualLeg(iTestY-2,iTestY)
              iTestFace = iTestX+(iTestY-3)*(maxPolyDegree+1)

              tmp_left = tmp_left + legprod * modalC_left(iTestFace,:)
              tmp_right = tmp_right + legprod * modalC_right(iTestFace,:)

              if (iTestX > 2) then
                legprod = ply_scalProdDualLeg(iTestX-2,iTestX) &
                  &       * ply_scalProdDualLeg(iTestY-2,iTestY)
                iTestFace = iTestX+(iTestY-3)*(maxPolyDegree+1) - 2

                tmp_left = tmp_left + legprod * modalC_left(iTestFace,:)
                tmp_right = tmp_right + legprod * modalC_right(iTestFace,:)
              end if

            end if

            p_c_left = jacobiDet * testZ_grad_val_left * tmp_left
            p_c_right = jacobiDet * testZ_grad_val_right * tmp_right

            if (iTestX > 1) then
              legprod = ply_scalProdDualLegDiff(iTestX-1,iTestX) &
                &       * legprod_square
              iTestFace = iTestX+(iTestY-1)*(maxPolyDegree+1) - 1

              tmp_left = legprod * modalA_left(iTestFace,:)
              tmp_right = legprod * modalA_right(iTestFace,:)

              if (iTestY > 2) then
                legprod = ply_scalProdDualLegDiff(iTestX-1,iTestX) &
                  &       * ply_scalProdDualLeg(iTestY-2,iTestY)
                iTestFace = iTestX+(iTestY-3)*(maxPolyDegree+1) - 1

                tmp_left = tmp_left + legprod * modalA_left(iTestFace,:)
                tmp_right = tmp_right + legprod * modalA_right(iTestFace,:)
              end if

              p_a_left = testZ_val_left * tmp_left
              p_a_right = testZ_val_right * tmp_right

            else

              p_a_left = 0.0_rk
              p_a_right = 0.0_rk

            end if

            if (iTestY > 1) then
              legprod = ply_scalProdDualLeg(iTestX,iTestX) &
                &       * ply_scalProdDualLegDiff(iTestY-1,iTestY)
              iTestFace = iTestX+(iTestY-2)*(maxPolyDegree+1)

              tmp_left = legprod * modalB_left(iTestFace,:)
              tmp_right = legprod * modalB_right(iTestFace,:)

              if (iTestX > 2) then
                legprod = ply_scalProdDualLeg(iTestX-2,iTestX) &
                  &       * ply_scalProdDualLegDiff(iTestY-1,iTestY)
                iTestFace = iTestX+(iTestY-2)*(maxPolyDegree+1) - 2

                tmp_left = tmp_left + legprod * modalB_left(iTestFace,:)
                tmp_right = tmp_right + legprod * modalB_right(iTestFace,:)
              end if

              p_b_left = testZ_val_left * tmp_left
              p_b_right = testZ_val_right * tmp_right

            else

              p_b_left = 0.0_rk
              p_b_right = 0.0_rk

            end if

            projection(iElem,testPos,:) = projection(iElem,testPos,:)       &
              &                         - outerNormalLeft * ( p_a_left(:)   &
              &                                             + p_b_left(:)   &
              &                                             + p_c_left(:) ) &
              &                         - outerNormalRight * ( p_a_right(:) &
              &                                              + p_b_right(:) &
              &                                              + p_c_right(:) )

          end do yTestLoop
        end do xTestLoop
      end do zTestLoop

    end do elementLoop

  end subroutine modg_project_stabViscNumFluxZ_Q

  !> Projection of the penalization terms (in modal representation) to the
  !! test functions.
  subroutine modg_project_penalization_Q( mesh, maxPolyDegree, &
                                        & kerneldata, penalizationdata)
    ! --------------------------------------------------------------------------
    !> The maximal polynomial degree of the modg scheme
    integer, intent(in) :: maxPolyDegree
    !> The cubical elements.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The data of the kernel. This one is updated by the projection.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> Volumetric data for the penalization
    type(atl_penalizationData_type), intent(in) :: penalizationdata
    ! --------------------------------------------------------------------------
    integer :: iElem, xTestFunc, yTestFunc, zTestFunc, testPos, &
             & xAnsFuncMin, xAnsFunc, yAnsFuncMin, yAnsFunc, zAnsFuncMin, &
             & zAnsFunc, ansPos
    real(kind=rk) :: jacobiDet, xScalProd, yScalProd, zScalProd
    integer :: nDoFs
    ! --------------------------------------------------------------------------

    jacobiDet = (0.5_rk*mesh%length)**3
    nDoFs = (maxPolyDegree+1)**3


    do iElem = 1, kerneldata%nTotal

      ! Now, we loop over all the test functions for this element and calculate
      ! the projection of the penalization terms onto this test functions.
      do testpos = 1,nDoFs
        xTestFunc = mod(testpos-1, maxpolydegree+1) + 1
        zTestFunc = (testpos-1)/(maxPolyDegree+1)**2 + 1
        yTestFunc = (testpos-1 - (zTestFunc-1)*(maxPolyDegree+1)**2) &
          &       / (maxPolyDegree+1) + 1
        ! Loop over relevant ans functions
        xAnsFuncMin = xTestFunc-2
        if (xAnsFuncMin < 1) then
          xAnsFuncMin = xTestFunc
        end if
        do xAnsFunc = xAnsFuncMin,xTestFunc,2
          ! Loop over relevant ans functions
          yAnsFuncMin = yTestFunc-2
          if( yAnsFuncMin < 1 ) then
            yAnsFuncMin = yTestFunc
          end if
          do yAnsFunc = yAnsFuncMin,yTestFunc,2
            ! Loop over relevant ans functions
            zAnsFuncMin = zTestFunc-2
            if( zAnsFuncMin < 1 ) then
              zAnsFuncMin = zTestFunc
            end if
            do zAnsFunc = zAnsFuncMin,zTestFunc,2

              ! get position of ansatz functions in the serialized list
              ! of dofs.
  anspos = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (zansfunc-1)*(maxpolydegree+1))*(maxpolydegree+1)

              ! project the current ansatz function onto the test function
              xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
              yScalProd = ply_scalProdDualLeg(yAnsFunc, yTestFunc)
              zScalProd = ply_scalProdDualLeg(zAnsFunc, zTestFunc)
              kernelData%state_der(iElem, testPos,  :) &
                               & = kernelData%state_der(iElem, testPos, :) &
                               & + xScalProd * yScalProd * zScalProd &
                               & * jacobiDet &
                               & * penalizationdata% &
                               &   penalization_data(iElem, ansPos,:)
            end do ! z ansatz functions
          end do ! y ansatz functions
        end do ! x ansatz functions

      end do
    end do ! elem loop


  end subroutine modg_project_penalization_Q


  !> Projection of the source terms (in modal representation) to the
  !! test functions.
  subroutine atl_modg_project_source( nScalars, sourcedata, scheme,  &
    &                                 mesh, kerneldata, currentLevel )
    ! --------------------------------------------------------------------------
    !> The number scalar variables in the equation system.
    integer, intent(in) :: nScalars
    !> The modal representation of the source
    type(atl_source_type), intent(in) :: sourcedata
    !> The parameters of the MODG scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    !> The cubical elements.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The data of the kernel. This one is updated by the projection.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The current level
    integer, intent(in) :: currentLevel
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Projection of the source terms onto the test functions.
    select case(scheme%basisType)
    case(Q_space)
      call atl_modg_project_source_Q(           &
        & sourcedata    = sourcedata,           &
        & nScalars      = nScalars,             &
        & mesh          = mesh,                 &
        & maxPolyDegree = scheme%maxPolyDegree, &
        & kerneldata    = kerneldata,           &
        & currentLevel  = currentLevel          )
    case(P_space)
      call atl_modg_project_source_P(           &
        & sourcedata    = sourcedata,           &
        & nScalars      = nScalars,             &
        & mesh          = mesh,                 &
        & maxPolyDegree = scheme%maxPolyDegree, &
        & kerneldata    = kerneldata,           &
        & currentLevel  = currentLevel          )
    end select

  end subroutine atl_modg_project_source


  !> Projection of the source terms (in modal representation) to the
  !! test functions.
  subroutine atl_modg_project_source_Q( sourcedata, nScalars, mesh, &
    &                                   maxPolyDegree, kerneldata,  &
    &                                   currentLevel                )
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
    integer :: iElem, elemPos, xTestFunc, yTestFunc, zTestFunc, testPos, &
             & xAnsFuncMin, xAnsFunc, yAnsFuncMin, yAnsFunc, zAnsFuncMin, &
             & zAnsFunc, ansPos, varPos, iSource, nSourceElems
    real(kind=rk) :: jacobiDet, xScalProd, yScalProd, zScalProd
    integer :: nDoFs
    ! --------------------------------------------------------------------------

    jacobiDet = (0.5_rk*mesh%length)**3
    nDoFs = (maxPolyDegree+1)**3

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
          do testpos = 1,nDoFs
            xTestFunc = mod(testpos-1, maxpolydegree+1) + 1
            zTestFunc = (testpos-1)/(maxPolyDegree+1)**2 + 1
            yTestFunc = (testpos-1 - (zTestFunc-1)*(maxPolyDegree+1)**2) &
              &       / (maxPolyDegree+1) + 1
            ! Loop over relevant ans functions
            xAnsFuncMin = xTestFunc-2
            if (xAnsFuncMin < 1) then
              xAnsFuncMin = xTestFunc
            end if
            do xAnsFunc = xAnsFuncMin,xTestFunc,2
              ! Loop over relevant ans functions
              yAnsFuncMin = yTestFunc-2
              if( yAnsFuncMin < 1 ) then
                yAnsFuncMin = yTestFunc
              end if
              do yAnsFunc = yAnsFuncMin,yTestFunc,2
                ! Loop over relevant ans functions
                zAnsFuncMin = zTestFunc-2
                if( zAnsFuncMin < 1 ) then
                  zAnsFuncMin = zTestFunc
                end if
                do zAnsFunc = zAnsFuncMin,zTestFunc,2

                  ! get position of ansatz functions in the serialized list
                  ! of dofs.
  anspos = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (zansfunc-1)*(maxpolydegree+1))*(maxpolydegree+1)

                  ! project the current ansatz function onto the test function
                  xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
                  yScalProd = ply_scalProdDualLeg(yAnsFunc, yTestFunc)
                  zScalProd = ply_scalProdDualLeg(zAnsFunc, zTestFunc)
                  kernelData%state_der(elemPos, testPos,  varPos)             &
                    & = kernelData%state_der(elemPos, testPos,  varPos)       &
                    &   + xScalProd * yScalProd * zScalProd                   &
                    &   * jacobiDet                                           &
                    &   * sourcedata%method(iSource)%val(iElem, ansPos, varPos)
                end do ! z ansatz functions
              end do ! y ansatz functions
            end do ! x ansatz functions

          end do

        end do ! elem loop
      end do ! variable loop
    end do ! source loop


  end subroutine atl_modg_project_source_Q


  !> Projection of the source terms (in modal representation) to the
  !! test functions.
  subroutine atl_modg_project_source_P( sourcedata, nScalars, mesh, &
    &                                   maxPolyDegree, kerneldata,  &
    &                                   currentLevel                )
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
    integer(kind=long_k) :: elemPos
    integer :: iElem, xTestFunc, yTestFunc, zTestFunc, testPos, xAnsFuncMin, &
      &        xAnsFunc, yAnsFuncMin, yAnsFunc, zAnsFuncMin, zAnsFunc,       &
      &        ansPos, varPos, iSource, nSourceElems, testPosMax
    real(kind=rk) :: jacobiDet, xScalProd, yScalProd, zScalProd
    ! --------------------------------------------------------------------------

    jacobiDet = (0.5_rk*mesh%length)**3

    do iSource = 1, size(sourcedata%method)

      nSourceElems = sourcedata%method(iSource)%elems(currentLevel)%nElems

      do varPos = 1, nScalars
        do iElem = 1, nSourceElems

          ! Position of the current element
          elempos = sourcedata%method(iSource)%elems(currentLevel) &
            &                 %posInTotal%val(iElem)

          ! Now, we loop over all the test functions for this element and
          ! calculate the projection of the source terms onto this
          ! test functions.
          xTestFunc = 1
          yTestFunc = 1
          zTestFunc = 1
  testposmax = (((maxpolydegree) + 1) &
    &   * ((maxpolydegree) + 2) &
    &   * ((maxpolydegree) + 3)) &
    & / 6
          do testPos = 1, testPosMax

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
                ! Loop over relevant ans functions
                zAnsFuncMin = zTestFunc-2
                if( zAnsFuncMin < 1 ) then
                  zAnsFuncMin = zTestFunc
                end if
                do zAnsFunc = zAnsFuncMin,zTestFunc,2

                  ! get position of ansatz functions in the serialized list
                  ! of dofs.
  ! integer divisions are no mistake here.
  anspos = (((xansfunc + yansfunc + zansfunc - 3) &
    &     * (xansfunc + yansfunc + zansfunc - 2) &
    &     * (xansfunc + yansfunc + zansfunc - 1)) &
    &   / 6 + 1)             &
    & + ((zansfunc-1) * (xansfunc + yansfunc + zansfunc -2) &
    &   - ((zansfunc-2) * (zansfunc-1)) / 2) &
    & + (yansfunc-1)

                  ! project the current ansatz function onto the test function
                  xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
                  yScalProd = ply_scalProdDualLeg(yAnsFunc, yTestFunc)
                  zScalProd = ply_scalProdDualLeg(zAnsFunc, zTestFunc)
                  kernelData%state_der(elemPos, testPos,  varPos) &
                      & = kernelData%state_der(elemPos, testPos,  varPos) &
                      & + xScalProd * yScalProd * zScalProd &
                      & * jacobiDet &
                      & * sourcedata%method(iSource)%val(iElem, ansPos, varPos)
                end do ! z ansatz functions
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
  elseif (ytestfunc .ne. 1) then
    ! next block
    xtestfunc = ytestfunc - 1
    ytestfunc = 1
    ztestfunc = ztestfunc + 1
  else
    ! next layer
    xtestfunc = ztestfunc + 1
    ytestfunc = 1
    ztestfunc = 1
  end if
          end do ! test functions
        end do ! elem loop
      end do ! variable loop
    end do ! source loop

  end subroutine atl_modg_project_source_P


  !> Projection of the numerical flux onto the testfunctions for Q_space.
  subroutine modg_project_numFlux_Q( nTotalFaces, &
    &                                nFaceDofs, nScalars, faceFlux, &
    &                                maxPolyDegree, length, &
    &                                nElems_fluid, dl_prod, projection, dirVec )
    ! --------------------------------------------------------------------------
    !> dimensions
    integer, intent(in) :: nTotalFaces, nFaceDofs, nScalars
    !> The numerical flux on the left face in modal representations.
    !! Dimension is nTotalFaces, (maxPolyDegree+1)^2 , nScalars, 2
    real(kind=rk), intent(in) :: faceFlux(:,:,:,:)
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
    !> direction (x,y,z)
    integer, intent(in) :: dirVec(3)
    ! --------------------------------------------------------------------------
    integer :: testFunc1, testFunc2, testFunc3, testPosFunc(3)
    integer :: testPos, ansPos(4)
    real(kind=rk) :: yScalProd, zScalProd(4)
    real(kind=rk) :: outerNormalLeft, outerNormalRight
    real(kind=rk) :: jacobiDetFaceProj
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iVar, min2mpd
    integer :: jk
    integer :: maxpd_m1
    ! --------------------------------------------------------------------------


    ! Jacobi determinant for the projections of the numerical fluxes onto the
    ! test functions
    jacobiDetFaceProj = (0.5_rk*length)**2
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)
    min2mpd = min(maxPolyDegree+1,2)
    maxpd_m1 = max(maxPolyDegree-1,0)


    ! Loop over all the test functions and project the numerical flux to them.
    do testFunc1 = 1,min2mpd
      faceValLeft  = ply_faceValLeftBndTest(testFunc1) * outerNormalLeft
      faceValRight = ply_faceValRightBndTest(testFunc1) * outerNormalRight

      do jk=1,min2mpd**2
        testFunc2 = mod(jk-1,min2mpd) + 1
        testFunc3 = (jk-1)/min2mpd + 1
        ! Just a single term to sum

        ! position of the test functions in the kerneldata
        testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  testpos = testposfunc(dirvec(1))                                      &
    &      + ( ( testposfunc(dirvec(2))-1)                             &
    &      + (testposfunc(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! calculate the projection of the ansatz and test function
        zScalProd(1) = dl_prod(2, testFunc3) &
          &            * dl_prod(2, testFunc2) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos(1) = testfunc2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        do iVar=1,nScalars
          projection(:nElems_fluid,testPos,iVar) &
            &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
            ! ... for the left face
            &      * ( faceValLeft &
            &          * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
            ! ... for the right face
            &        + faceValRight &
            &          * faceFlux(:nElems_fluid,ansPos(1),iVar,2) )
        end do
      end do

      do jk=1,min2mpd*maxpd_m1
        testFunc2 = mod(jk-1,min2mpd) + 1
        testFunc3 = (jk-1)/min2mpd + 3
        ! Two terms to sum

        ! position of the test functions in the kerneldata
        testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  testpos = testposfunc(dirvec(1))                                      &
    &      + ( ( testposfunc(dirvec(2))-1)                             &
    &      + (testposfunc(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! calculate the projection of the ansatz and test function
        yScalProd = dl_prod(2, testFunc2) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos(1) = testfunc2                                      &
    &      + ( ( testfunc3-2-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(2) = testfunc2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! calculate the projection of the ansatz and test function
        zScalProd(1) = dl_prod(1, testFunc3) * yScalProd
        zScalProd(2) = dl_prod(2, testFunc3) * yScalProd

        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        do iVar=1,nScalars
          projection(:nElems_fluid,testPos,iVar) &
            &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
            ! ... for the left face
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
            &  - zScalProd(2) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) )
        end do

      end do

      do jk=1,min2mpd*maxpd_m1
        testFunc3 = mod(jk-1,min2mpd) + 1
        testFunc2 = (jk-1)/min2mpd + 3
        ! Two terms to sum

        ! position of the test functions in the kerneldata
        testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  testpos = testposfunc(dirvec(1))                                      &
    &      + ( ( testposfunc(dirvec(2))-1)                             &
    &      + (testposfunc(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ! calculate the projection of the ansatz and test function
        zScalProd(1) = dl_prod(2, testFunc3) &
          &            * dl_prod(1, testFunc2) * jacobiDetFaceProj
        zScalProd(2) = dl_prod(2, testFunc3) &
          &            * dl_prod(2, testFunc2) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos(1) = testfunc2-2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(2) = testfunc2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        do iVar=1,nScalars
          projection(:nElems_fluid,testPos,iVar) &
            &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
            ! ... for the left face
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
            &  - zScalProd(2) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) )
        end do

      end do

      do jk=1,maxpd_m1**2
        testFunc2 = mod(jk-1,maxpd_m1) + 3
        testFunc3 = (jk-1)/maxpd_m1 + 3
        ! Four terms to sum

        ! position of the test functions in the kerneldata
        testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  testpos = testposfunc(dirvec(1))                                      &
    &      + ( ( testposfunc(dirvec(2))-1)                             &
    &      + (testposfunc(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ! calculate the projection of the ansatz and test function
        zScalProd(1) = dl_prod(1, testFunc3) * dl_prod(1, testFunc2) &
          &                                  * jacobiDetFaceProj
        zScalProd(2) = dl_prod(1, testFunc3) * dl_prod(2, testFunc2) &
          &                                  * jacobiDetFaceProj
        zScalProd(3) = dl_prod(2, testFunc3) * dl_prod(1, testFunc2) &
          &                                  * jacobiDetFaceProj
        zScalProd(4) = dl_prod(2, testFunc3) * dl_prod(2, testFunc2) &
          &                                  * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos(1) = testfunc2-2                                      &
    &      + ( ( testfunc3-2-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(2) = testfunc2                                      &
    &      + ( ( testfunc3-2-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(3) = testfunc2-2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(4) = testfunc2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        do iVar=1,nScalars
          projection(:nElems_fluid,testPos,iVar) &
            &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
            ! ... for the left face
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
            &  - zScalProd(2) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) ) &
            &  - zScalProd(3) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(3),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(3),iVar,2) ) &
            &  - zScalProd(4) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(4),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(4),iVar,2) )
        end do

      end do

    end do
  end subroutine modg_project_numFlux_Q



  !> Projection of the numerical flux onto the differentiated testfunctions.
  subroutine modg_project_numFlux_diffTest_Q( nScalars, faceFlux,      &
    & maxPolyDegree, length, nElems_fluid, dl_prod, projection, dirVec )
    ! --------------------------------------------------------------------------
    !> dimensions
    integer, intent(in) :: nScalars
    !> The numerical flux on the left face in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(inout) :: faceFlux(:,:,:,:)
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
    !> direction (x,y,z)
    integer, intent(in) :: dirVec(3)
    ! --------------------------------------------------------------------------
    integer :: testFunc1, testFunc2, testFunc3, testPosFunc(3)
    integer :: testPos, ansPos(4)
    real(kind=rk) :: yScalProd, zScalProd(4)
    real(kind=rk) :: jacobiDetFaceProj
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iVar, min2mpd
    integer :: jk
    integer :: maxpd_m1
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for the projections of the numerical fluxes onto the
    ! test functions

    jacobiDetFaceProj = (0.5_rk*length)**2
    zscalProd = 0.0_rk

    min2mpd = min(maxPolyDegree+1,5)
    maxpd_m1 = max(maxPolyDegree-4,0)

    do testFunc1 = 2,maxPolyDegree+1
      faceValLeft  = ply_faceValLeftBndGradTest(testFunc1)
      faceValRight = ply_faceValRightBndGradTest(testFunc1)

      do jk=1,(min2mpd)**2
        testFunc2 = mod(jk-1,min2mpd) + 1
        testFunc3 = (jk-1)/min2mpd + 1
        ! Just a single term to sum

        ! position of the test functions in the kerneldata
        testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  testpos = testposfunc(dirvec(1))                                      &
    &      + ( ( testposfunc(dirvec(2))-1)                             &
    &      + (testposfunc(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ! calculate the projection of the ansatz and test function
        zScalProd(1) = dl_prod(2, testFunc3) &
        &            * dl_prod(2, testFunc2) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos(1) = testfunc2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        do iVar=1,nScalars

          projection(:nElems_fluid,testPos,iVar) &
            &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
            ! ... for the left face
            &      * ( faceValLeft &
            &          * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
            ! ... for the right face
            &        + faceValRight &
            &          * faceFlux(:nElems_fluid,ansPos(1),iVar,2) )
        end do

      end do

      ! for the testfunct3 greater than 5 and testfunc2 less than 5
      do jk=1,maxpd_m1*5
        testFunc2 = mod(jk-1,min2mpd) + 1
        testFunc3 = (jk-1)/min2mpd + 6
        ! Two terms to sum

        ! position of the test functions in the kerneldata
        testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  testpos = testposfunc(dirvec(1))                                      &
    &      + ( ( testposfunc(dirvec(2))-1)                             &
    &      + (testposfunc(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ! calculate the projection of the ansatz and test function
        yScalProd = dl_prod(2, testFunc2) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos(1) = testfunc2                                      &
    &      + ( ( testfunc3-2-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(2) = testfunc2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! calculate the projection of the ansatz and test function
        zScalProd(1) = dl_prod(1, testFunc3) * yScalProd
        zScalProd(2) = dl_prod(2, testFunc3) * yScalProd

       ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        do iVar=1,nScalars
          projection(:nElems_fluid,testPos,iVar) &
            &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
            ! ... for the left face
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
            &  - zScalProd(2) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) )
        end do
      end do

      ! for the testfunct2 greater than 5 and testfunc3 less than 5
      do jk=1,maxpd_m1*5
        testFunc3 = mod(jk-1,min2mpd) + 1
        testFunc2 = (jk-1)/min2mpd + 6
        ! Two terms to sum

        ! position of the test functions in the kerneldata
        testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  testpos = testposfunc(dirvec(1))                                      &
    &      + ( ( testposfunc(dirvec(2))-1)                             &
    &      + (testposfunc(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ! calculate the projection of the ansatz and test function
        zScalProd(1) = dl_prod(2, testFunc3) &
          &            * dl_prod(1, testFunc2) * jacobiDetFaceProj
        zScalProd(2) = dl_prod(2, testFunc3) &
          &            * dl_prod(2, testFunc2) * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos(1) = testfunc2-2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(2) = testfunc2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        do iVar=1,nScalars
          projection(:nElems_fluid,testPos,iVar) &
            &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
            ! ... for the left face
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
            &  - zScalProd(2) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) )
        end do

      end do

      do jk=1,(maxpd_m1)**2
        testFunc2 = mod(jk-1,maxpd_m1) + 6
        testFunc3 = (jk-1)/maxpd_m1 + 6
        ! Four terms to sum
        ! position of the test functions in the kerneldata
        testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  testpos = testposfunc(dirvec(1))                                      &
    &      + ( ( testposfunc(dirvec(2))-1)                             &
    &      + (testposfunc(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ! calculate the projection of the ansatz and test function
        zScalProd(1) = dl_prod(1, testFunc3) * dl_prod(1, testFunc2) &
          &                                  * jacobiDetFaceProj
        zScalProd(2) = dl_prod(1, testFunc3) * dl_prod(2, testFunc2) &
          &                                  * jacobiDetFaceProj
        zScalProd(3) = dl_prod(2, testFunc3) * dl_prod(1, testFunc2) &
          &                                  * jacobiDetFaceProj
        zScalProd(4) = dl_prod(2, testFunc3) * dl_prod(2, testFunc2) &
          &                                  * jacobiDetFaceProj

        ! the position of the modal coefficeint of this ansatz functions
  anspos(1) = testfunc2-2                                      &
    &      + ( ( testfunc3-2-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(2) = testfunc2                                      &
    &      + ( ( testfunc3-2-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(3) = testfunc2-2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(4) = testfunc2                                      &
    &      + ( ( testfunc3-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ! buffer the result in kernel data and take care of the outer
        ! surface unit normal
        do iVar=1,nScalars
          projection(:nElems_fluid,testPos,iVar) &
            &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
            ! ... for the left face
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
            &  - zScalProd(2) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) ) &
            &  - zScalProd(3) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(3),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(3),iVar,2) ) &
            &  - zScalProd(4) &
            &    * ( faceValLeft &
            &        * faceFlux(:nElems_fluid,ansPos(4),iVar,1) &
            ! ... for the right face
            &      + faceValRight &
            &        * faceFlux(:nElems_fluid,ansPos(4),iVar,2) )
        end do
      end do

    end do
  end subroutine modg_project_numFlux_difftest_Q


  !> Projection of the numerical flux onto the testfunctions for P_space.
  subroutine modg_project_numFlux_P( nTotalFaces, nFaceDofs, nScalars,    &
    & faceFlux, maxPolyDegree, length, nElems_fluid, dl_prod, projection, &
    & dirVec                                                              )
    ! --------------------------------------------------------------------------
    !> dimensions
    integer, intent(in) :: nTotalFaces, nFaceDofs, nScalars
    !> The numerical flux on the left face in modal representations.
    !! Dimension is (maxPolyDegree+1)^2 , nScalars
    real(kind=rk), intent(in) :: faceFlux(nTotalFaces,nFaceDofs,nScalars,2)
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
    !> direction (x,y,z)
    integer, intent(in) :: dirVec(3)
    ! --------------------------------------------------------------------------
    integer :: testFunc1, testFunc2, testFunc3, testPosFunc(3)
    integer :: testPos, ansPos(4)
    real(kind=rk) :: yScalProd, zScalProd(4)
    real(kind=rk) :: outerNormalLeft, outerNormalRight
    real(kind=rk) :: jacobiDetFaceProj
    real(kind=rk) :: faceValLeft, faceValRight
    integer :: iVar
    integer :: min2mpd
    ! --------------------------------------------------------------------------


    ! Jacobi determinant for the projections of the numerical fluxes onto the
    ! test functions
    jacobiDetFaceProj = (0.5_rk*length)**2
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)
    min2mpd = min(maxPolyDegree+1,2)


    ! Loop over all the test functions and project the numerical flux to them.
    do testFunc1 = 1,min2mpd
      faceValLeft = ply_faceValLeftBndTest(testFunc1) * outerNormalLeft
      faceValRight = ply_faceValRightBndTest(testFunc1) * outerNormalRight

      do testFunc2 = 1,min(2, maxPolyDegree+1 - (testFunc1-1))
        do testFunc3 = 1,min(2, maxPolyDegree+1 - (testFunc1-1) - (testfunc2-1))
          ! Just a single term to sum

          ! position of the test functions in the kerneldata
          testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  ! integer divisions are no mistake here.
  testpos = (((testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 3) &
    &     * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 2) &
    &     * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((testposfunc(dirvec(3))-1) * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) -2) &
    &   - ((testposfunc(dirvec(3))-2) * (testposfunc(dirvec(3))-1)) / 2) &
    & + (testposfunc(dirvec(2))-1)

          ! calculate the projection of the ansatz and test function
          zScalProd(1) = dl_prod(2, testFunc3) &
            &            * dl_prod(2, testFunc2) * jacobiDetFaceProj

          ! the position of the modal coefficeint of this ansatz functions
  ! integer divisions are no mistake here.
  anspos(1) = ((((testfunc2 - 1) + (testfunc3 - 1))            &
    &   * (((testfunc2 - 1) + (testfunc3 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3 - 1)

          ! buffer the result in kernel data and take care of the outer
          ! surface unit normal
          do iVar=1, nScalars
              projection(:nElems_fluid,testPos,iVar) &
                &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
                ! ... for the left face
                &      * ( faceValLeft &
                &          * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
                ! ... for the right face
                &        + faceValRight &
                &          * faceFlux(:nElems_fluid,ansPos(1),iVar,2) )
          end do

        end do
      end do

      do testFunc2 = 1,min(2, maxPolyDegree+1 - (testFunc1-1))
        do testFunc3 = 3,maxPolyDegree+1 - (testFunc1-1) - (testfunc2-1)
          ! Two terms to sum

          ! position of the test functions in the kerneldata
          testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  ! integer divisions are no mistake here.
  testpos = (((testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 3) &
    &     * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 2) &
    &     * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((testposfunc(dirvec(3))-1) * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) -2) &
    &   - ((testposfunc(dirvec(3))-2) * (testposfunc(dirvec(3))-1)) / 2) &
    & + (testposfunc(dirvec(2))-1)

          ! calculate the projection of the ansatz and test function
          yScalProd = dl_prod(2, testFunc2) * jacobiDetFaceProj

          ! the position of the modal coefficeint of this ansatz functions
  ! integer divisions are no mistake here.
  anspos(1) = ((((testfunc2 - 1) + (testfunc3-2 - 1))            &
    &   * (((testfunc2 - 1) + (testfunc3-2 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3-2 - 1)
  ! integer divisions are no mistake here.
  anspos(2) = ((((testfunc2 - 1) + (testfunc3 - 1))            &
    &   * (((testfunc2 - 1) + (testfunc3 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3 - 1)

          ! calculate the projection of the ansatz and test function
          zScalProd(1) = dl_prod(1, testFunc3) * yScalProd
          zScalProd(2) = dl_prod(2, testFunc3) * yScalProd

          ! buffer the result in kernel data and take care of the outer
          ! surface unit normal
          do iVar=1, nScalars
            projection(:nElems_fluid,testPos,iVar) &
              &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
              ! ... for the left face
              &    * ( faceValLeft &
              &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
              ! ... for the right face
              &      + faceValRight &
              &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
              &  - zScalProd(2) &
              &    * ( faceValLeft &
              &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
              ! ... for the right face
              &      + faceValRight &
              &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) )
          end do

        end do
      end do

      do testFunc2 = 3,maxPolyDegree+1 - (testFunc1-1)
        do testFunc3 = 1,min(2,maxPolyDegree+1 - (testFunc1-1) - (testfunc2-1))
          ! Two terms to sum

          ! position of the test functions in the kerneldata
          testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  ! integer divisions are no mistake here.
  testpos = (((testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 3) &
    &     * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 2) &
    &     * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((testposfunc(dirvec(3))-1) * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) -2) &
    &   - ((testposfunc(dirvec(3))-2) * (testposfunc(dirvec(3))-1)) / 2) &
    & + (testposfunc(dirvec(2))-1)


          ! calculate the projection of the ansatz and test function
          zScalProd(1) = dl_prod(2, testFunc3) &
            &            * dl_prod(1, testFunc2) * jacobiDetFaceProj
          zScalProd(2) = dl_prod(2, testFunc3) &
            &            * dl_prod(2, testFunc2) * jacobiDetFaceProj

          ! the position of the modal coefficeint of this ansatz functions
  ! integer divisions are no mistake here.
  anspos(1) = ((((testfunc2-2 - 1) + (testfunc3 - 1))            &
    &   * (((testfunc2-2 - 1) + (testfunc3 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3 - 1)
  ! integer divisions are no mistake here.
  anspos(2) = ((((testfunc2 - 1) + (testfunc3 - 1))            &
    &   * (((testfunc2 - 1) + (testfunc3 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3 - 1)

          ! buffer the result in kernel data and take care of the outer
          ! surface unit normal
          do iVar=1, nScalars
            projection(:nElems_fluid,testPos,iVar) &
              &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
              ! ... for the left face
              &    * ( faceValLeft &
              &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
              ! ... for the right face
              &      + faceValRight &
              &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
              &  - zScalProd(2) &
              &    * ( faceValLeft &
              &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
              ! ... for the right face
              &      + faceValRight &
              &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) )
          end do

        end do
      end do

      do testFunc2 = 3,maxPolyDegree+1 - (testFunc1-1)
        do testFunc3 = 3,maxPolyDegree+1 - (testFunc1-1) - (testfunc2-1)
          ! Four terms to sum

          ! position of the test functions in the kerneldata
          testPosFunc(:) = (/testFunc1, testFunc2, testFunc3/)
  ! integer divisions are no mistake here.
  testpos = (((testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 3) &
    &     * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 2) &
    &     * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((testposfunc(dirvec(3))-1) * (testposfunc(dirvec(1)) + testposfunc(dirvec(2)) + testposfunc(dirvec(3)) -2) &
    &   - ((testposfunc(dirvec(3))-2) * (testposfunc(dirvec(3))-1)) / 2) &
    & + (testposfunc(dirvec(2))-1)


          ! calculate the projection of the ansatz and test function
          zScalProd(1) = dl_prod(1, testFunc3) * dl_prod(1, testFunc2) &
            &                                  * jacobiDetFaceProj
          zScalProd(2) = dl_prod(1, testFunc3) * dl_prod(2, testFunc2) &
            &                                  * jacobiDetFaceProj
          zScalProd(3) = dl_prod(2, testFunc3) * dl_prod(1, testFunc2) &
            &                                  * jacobiDetFaceProj
          zScalProd(4) = dl_prod(2, testFunc3) * dl_prod(2, testFunc2) &
            &                                  * jacobiDetFaceProj

          ! the position of the modal coefficeint of this ansatz functions
  ! integer divisions are no mistake here.
  anspos(1) = ((((testfunc2-2 - 1) + (testfunc3-2 - 1))            &
    &   * (((testfunc2-2 - 1) + (testfunc3-2 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3-2 - 1)
  ! integer divisions are no mistake here.
  anspos(2) = ((((testfunc2 - 1) + (testfunc3-2 - 1))            &
    &   * (((testfunc2 - 1) + (testfunc3-2 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3-2 - 1)
  ! integer divisions are no mistake here.
  anspos(3) = ((((testfunc2-2 - 1) + (testfunc3 - 1))            &
    &   * (((testfunc2-2 - 1) + (testfunc3 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3 - 1)
  ! integer divisions are no mistake here.
  anspos(4) = ((((testfunc2 - 1) + (testfunc3 - 1))            &
    &   * (((testfunc2 - 1) + (testfunc3 - 1)) + 1)) / 2 + 1) &
    & + (testfunc3 - 1)

          ! buffer the result in kernel data and take care of the outer
          ! surface unit normal
          do iVar=1, nScalars
            projection(:nElems_fluid,testPos,iVar) &
              &  = projection(:nElems_fluid,testPos,iVar) - zScalProd(1) &
              ! ... for the left face
              &    * ( faceValLeft &
              &        * faceFlux(:nElems_fluid,ansPos(1),iVar,1) &
              ! ... for the right face
              &      + faceValRight &
              &        * faceFlux(:nElems_fluid,ansPos(1),iVar,2) ) &
              &  - zScalProd(2) &
              &    * ( faceValLeft &
              &        * faceFlux(:nElems_fluid,ansPos(2),iVar,1) &
              ! ... for the right face
              &      + faceValRight &
              &        * faceFlux(:nElems_fluid,ansPos(2),iVar,2) ) &
              &  - zScalProd(3) &
              &    * ( faceValLeft &
              &        * faceFlux(:nElems_fluid,ansPos(3),iVar,1) &
              ! ... for the right face
              &      + faceValRight &
              &        * faceFlux(:nElems_fluid,ansPos(3),iVar,2) ) &
              &  - zScalProd(4) &
              &    * ( faceValLeft &
              &        * faceFlux(:nElems_fluid,ansPos(4),iVar,1) &
              ! ... for the right face
              &      + faceValRight &
              &        * faceFlux(:nElems_fluid,ansPos(4),iVar,2) )
          end do

        end do
      end do

    end do

  end subroutine modg_project_numFlux_P


  !> Lift the mean on the face to a positive value if necessary.
  !!
  !! In some equations, some variables may not be allowed to be negative,
  !! for example the density and energy in flow simulations.
  !! With this routine we enforce this limitation for the projected mean
  !! state on the element surfaces.
  subroutine atl_modg_ensure_pos_facemean(                                  &
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
    integer :: iVar, iDir
    real(kind=rk) :: epsfact
    ! --------------------------------------------------------------------------

    epsfact = epsilon(volState(1,1,1))

    do iVar=1,nScalars
      if (ensure_positivity(iVar)) then

        ! We assume the integral mean of the element is positive for all
        ! elements. The integral mean on the surface will now be lifted to
        ! a small fraction of that element mean if it is below that threshold.
        do iDir=1,3
          facerep(iDir)%dat(1:nElems_fluid,1,iVar,1)             &
            & = max( facerep(iDir)%dat(1:nElems_fluid,1,iVar,1), &
            &        epsfact * volState(1:nElems_fluid,1,iVar)   )
          facerep(iDir)%dat(1:nElems_fluid,1,iVar,2)             &
            & = max( facerep(iDir)%dat(1:nElems_fluid,1,iVar,2), &
            &        epsfact * volState(1:nElems_fluid,1,iVar)   )
        end do

      end if
    end do

  end subroutine atl_modg_ensure_pos_facemean


  !> Projects modal representation of each cell to its faces, i.e.
  !! this subroutine creates a modal representation on the faces.
  subroutine atl_modg_modalVolToModalFace( nElems_fluid, length , volState, &
    &                                      faceRep, nScalars, nDerivatives, &
    &                                      modg                             )
    ! --------------------------------------------------------------------------
    !> Volumetric, modal states for each element.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> Modal representation on the face (will be updated by this routine for all
    !! fluid elements in mesh).
    type(atl_faceRep_type), intent(inout) :: faceRep(:)
    !> number of variables
    integer, intent(in) :: nScalars
    !> number of derivatives
    integer, intent(in) :: nDerivatives
    !> number of fluid elements
    integer, intent(in) :: nElems_fluid
    !> length of cubic element
    real(kind=rk), intent(in) :: length
    !> The parameters of your modg scheme.
    type(atl_modg_scheme_type), intent(in) :: modg
    ! --------------------------------------------------------------------------
    integer :: iDir, spaceDir
    integer :: nFaces
    real(kind=rk), allocatable :: volState_Q(:,:,:)
    real(kind=rk), allocatable :: faceState_Q(:,:,:,:)
    ! --------------------------------------------------------------------------

    do iDir = 1,3
      faceRep(iDir)%dat(1:nElems_fluid,:,:,:) = 0.0_rk
    end do

    ! Iterate over all the fluid elements and project its modal representations
    ! to the faces of the element.
    select case(modg%basisType)
      case(Q_space) ! Q tensor product ansatz functions

        ! Project the state on the faces in each direction.
        call modg_volToFace_Q_x(                &
          & volState      = volstate,           &
          & maxPolyDegree = modg%maxPolyDegree, &
          & nScalars      = nScalars,           &
          & nElems        = nElems_fluid,       &
          & faceState     = faceRep(1)%dat      )
        call modg_volToFace_Q_y(                &
          & volState      = volstate,           &
          & maxPolyDegree = modg%maxPolyDegree, &
          & nScalars      = nScalars,           &
          & nElems        = nElems_fluid,       &
          & faceState     = faceRep(2)%dat      )
        call modg_volToFace_Q_z(                &
          & volState      = volstate,           &
          & maxPolyDegree = modg%maxPolyDegree, &
          & nScalars      = nScalars,           &
          & nElems        = nElems_fluid,       &
          & faceState     = faceRep(3)%dat      )

        if (nDerivatives == 1) then
          do iDir = 1, 6
            ! Project to the face in the current direction.
            spaceDir = qAxis(iDir)
            nFaces = size(faceRep(spaceDir)%dat,1)

            ! Project derivatives to the faces
              call atl_modg_volToFace_grad_Q(           &
                & volState      = volstate,             &
                & maxPolyDegree = modg%maxPolyDegree,   &
                & faceDir       = iDir,                 &
                & nScalars      = nScalars,             &
                & nElems        = nElems_fluid,         &
                & elemLength    = length,               &
                & faceState     = faceRep(spaceDir)%dat )

          end do
        end if

      case(P_space) ! P-tensor product ansatz functions

        ! Now, iterate over all the faces and project to this face.
        do iDir = 1, 6
          ! Project to the face in the current direction.
          spaceDir = qAxis(iDir)

          call atl_modg_volToFace_P(                                  &
            & nTotalElems   = size(volstate,1),                       &
            & nTotalFaces   = size(faceRep(spaceDir)%dat,1),          &
            & nDofs         = size(volstate,2),                       &
            & nFaceDofs     = size(faceRep(spaceDir)%dat,2),          &
            & nScalars      = nScalars,                               &
            & volState      = volstate,                               &
            & maxPolyDegree = modg%maxPolyDegree,                     &
            & faceDir       = iDir,                                   &
            & nElems        = nElems_fluid,                           &
            & faceState     = faceRep(spaceDir)%dat(:,:,1:nScalars,:) )
        end do

        if (nDerivatives == 1) then

          ! For the projection of the derivatives to the faces in P_space
          ! we can use the same subroutine that we use for the projection
          ! in Q_space but we need to change the polynomial space to a
          ! Q_space representation and change it back to P_space afterwards.

          allocate(volstate_Q(nElems_fluid,(modg%maxPolyDegree+1)**3,nScalars))

          call ply_change_poly_space( inspace    = P_space,            &
            &                         instate    = volstate,           &
            &                         outstate   = volstate_Q,         &
            &                         maxPolyDeg = modg%maxPolyDegree, &
            &                         nElems     = nElems_fluid,       &
            &                         nVars      = nScalars,           &
            &                         nDims      = 3                   )

          do iDir = 1, 6
            ! Project to the face in the current direction.
            spaceDir = qAxis(iDir)

            nFaces = size(faceRep(spaceDir)%dat,1)

            allocate(faceState_Q(nFaces,(modg%maxPolyDegree+1)**2,4*nScalars,2))

            ! Change space for the left and for the right face seperatly
            call ply_change_poly_space( inspace    = P_space,              &
              &                         instate    = faceRep(spaceDir)     &
              &                                       %dat(:,:,:,1),       &
              &                         outstate   = faceState_Q(:,:,:,1), &
              &                         maxPolyDeg = modg%maxPolyDegree,   &
              &                         nElems     = nFaces,               &
              &                         nVars      = 4*nScalars,           &
              &                         nDims      = 2                     )

            call ply_change_poly_space( inspace    = P_space,              &
              &                         instate    = faceRep(spaceDir)     &
              &                                       %dat(:,:,:,2),       &
              &                         outstate   = faceState_Q(:,:,:,2), &
              &                         maxPolyDeg = modg%maxPolyDegree,   &
              &                         nElems     = nFaces,               &
              &                         nVars      = 4*nScalars,           &
              &                         nDims      = 2                     )

            ! Project derivatives to the faces
            call atl_modg_volToFace_grad_Q(         &
              & volState      = volstate_Q,         &
              & maxPolyDegree = modg%maxPolyDegree, &
              & faceDir       = iDir,               &
              & nScalars      = nScalars,           &
              & nElems        = nElems_fluid,       &
              & elemLength    = length,             &
              & faceState     = faceState_Q         )

            call ply_change_poly_space( inspace    = Q_space,              &
              &                         instate    = faceState_Q(:,:,:,1), &
              &                         outstate   = faceRep(spaceDir)     &
              &                                       %dat(:,:,:,1),       &
              &                         maxPolyDeg = modg%maxPolyDegree,   &
              &                         nElems     = nFaces,               &
              &                         nVars      = 4*nScalars,           &
              &                         nDims      = 2                     )

            call ply_change_poly_space( inspace    = Q_space,              &
              &                         instate    = faceState_Q(:,:,:,2), &
              &                         outstate   = faceRep(spaceDir)     &
              &                                       %dat(:,:,:,2),       &
              &                         maxPolyDeg = modg%maxPolyDegree,   &
              &                         nElems     = nFaces,               &
              &                         nVars      = 4*nScalars,           &
              &                         nDims      = 2                     )


            deallocate(faceState_Q)
          end do

          deallocate(volstate_Q)

        end if

      case default
        write(logUnit(1),*) 'ERROR in atl_modg_modalVolToModalFace:'
        write(logUnit(1),*) 'Unknown tensor product, stopping ...'
        call tem_abort()
    end select


  end subroutine atl_modg_modalVolToModalFace


  !> Project modal representation of an element to its two faces in X.
  subroutine modg_volToFace_Q_x( nScalars, volState, maxPolyDegree, nElems, &
    &                            faceState )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: nScalars
    !> The modal representation in the volume. First dimension is the number of
    !! volumetric numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the faces in X direction
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: pos, facePos, iAnsX
    integer :: iVar
    integer :: mpd1, mpd1_square
    integer :: iElem
    integer :: iVEF, ipv
    integer :: nIndeps
    ! --------------------------------------------------------------------------

    mpd1 = maxpolydegree+1
    mpd1_square = mpd1**2
    nIndeps = nScalars*nElems*mpd1_square

    ! Split off first ansatz function, to avoid initialization with 0.
    iAnsX = 1

    !$NEC ivdep
    do iVEF=1,nIndeps
      iVar = (iVEF-1) / (nElems*mpd1_square) + 1
      ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
      facepos = (ipv-1) / nElems + 1
      iElem = mod(ipv-1, nElems) + 1
      ! get position of the current ansatz function
      pos = (facepos-1)*mpd1 + iAnsX

      faceState(iElem,facePos,iVar,1) = volState(iElem, pos, iVar)
      faceState(iElem,facePos,iVar,2) = volState(iElem, pos, iVar)
    end do

    ! Legendre polynomial on the right is just the sum of all modes,
    ! but on the left we alternatingly need to add and subtract the modes.
    ! Thus, we split the loop to allow the correct operation to be used for
    ! each mode.
    do iAnsX=3,maxPolyDegree+1,2

      !$NEC ivdep
      do iVEF=1,nIndeps
        iVar = (iVEF-1) / (nElems*mpd1_square) + 1
        ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
        facepos = (ipv-1) / nElems + 1
        iElem = mod(ipv-1, nElems) + 1
        ! get position of the current ansatz function
        pos = (facepos-1)*mpd1 + iAnsX

        faceState(iElem,facePos,iVar,1) = faceState(iElem,facePos,iVar,1) &
          &                             + volState(iElem, pos, iVar)
        faceState(iElem,facePos,iVar,2) = faceState(iElem,facePos,iVar,2) &
          &                             + volState(iElem, pos, iVar)
      end do

    end do
    do iAnsX=2,maxPolyDegree+1,2

      !$NEC ivdep
      do iVEF=1,nIndeps
        iVar = (iVEF-1) / (nElems*mpd1_square) + 1
        ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
        facepos = (ipv-1) / nElems + 1
        iElem = mod(ipv-1, nElems) + 1
        ! get position of the current ansatz function
        pos = (facepos-1)*mpd1 + iAnsX

        faceState(iElem,facePos,iVar,1) = faceState(iElem,facePos,iVar,1) &
          &                             - volState(iElem, pos, iVar)
        faceState(iElem,facePos,iVar,2) = faceState(iElem,facePos,iVar,2) &
          &                             + volState(iElem, pos, iVar)
      end do

    end do

  end subroutine modg_voltoface_Q_x


  !> Project modal representation of an element to its two faces in Y.
  subroutine modg_volToFace_Q_y( nScalars, volState, maxPolyDegree, nElems, &
    &                            faceState )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: nScalars
    !> The modal representation in the volume. First dimension is the number of
    !! volumetric numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the faces in X direction
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: pos, facePos
    integer :: iAnsX, iAnsY, iAnsZ
    integer :: iVar
    integer :: mpd1, mpd1_square
    integer :: iElem
    integer :: iVEF, ipv
    integer :: nIndeps
    ! --------------------------------------------------------------------------

    mpd1 = maxpolydegree+1
    mpd1_square = mpd1**2
    nIndeps = nScalars*nElems*mpd1_square

    ! Split off first ansatz function, to avoid initialization with 0.
    iAnsY = 1

    !$NEC ivdep
    do iVEF=1,nIndeps
      iVar = (iVEF-1) / (nElems*mpd1_square) + 1
      ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
      facepos = (ipv-1) / nElems + 1
      iElem = mod(ipv-1, nElems) + 1
      iAnsX = mod(facepos-1, mpd1) + 1
      iAnsZ = (facepos-1)/mpd1 + 1
      ! get position of the current ansatz function
      pos = iAnsX + (iAnsY-1)*mpd1 + (iAnsZ-1)*mpd1_square

      faceState(iElem,facePos,iVar,1) = volState(iElem, pos, iVar)
      faceState(iElem,facePos,iVar,2) = volState(iElem, pos, iVar)
    end do

    ! Legendre polynomial on the right is just the sum of all modes,
    ! but on the left we alternatingly need to add and subtract the modes.
    ! Thus, we split the loop to allow the correct operation to be used for
    ! each mode.
    do iAnsY=3,maxPolyDegree+1,2

      !$NEC ivdep
      do iVEF=1,nIndeps
        iVar = (iVEF-1) / (nElems*mpd1_square) + 1
        ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
        facepos = (ipv-1) / nElems + 1
        iElem = mod(ipv-1, nElems) + 1
        iAnsX = mod(facepos-1, mpd1) + 1
        iAnsZ = (facepos-1)/mpd1 + 1
        ! get position of the current ansatz function
        pos = iAnsX + (iAnsY-1)*mpd1 + (iAnsZ-1)*mpd1_square

        faceState(iElem,facePos,iVar,1) = faceState(iElem,facePos,iVar,1) &
          &                             + volState(iElem, pos, iVar)
        faceState(iElem,facePos,iVar,2) = faceState(iElem,facePos,iVar,2) &
          &                             + volState(iElem, pos, iVar)
      end do

    end do
    do iAnsY=2,maxPolyDegree+1,2

      !$NEC ivdep
      do iVEF=1,nIndeps
        iVar = (iVEF-1) / (nElems*mpd1_square) + 1
        ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
        facepos = (ipv-1) / nElems + 1
        iElem = mod(ipv-1, nElems) + 1
        iAnsX = mod(facepos-1, mpd1) + 1
        iAnsZ = (facepos-1)/mpd1 + 1
        ! get position of the current ansatz function
        pos = iAnsX + (iAnsY-1)*mpd1 + (iAnsZ-1)*mpd1_square

        faceState(iElem,facePos,iVar,1) = faceState(iElem,facePos,iVar,1) &
          &                             - volState(iElem, pos, iVar)
        faceState(iElem,facePos,iVar,2) = faceState(iElem,facePos,iVar,2) &
          &                             + volState(iElem, pos, iVar)
      end do

    end do

  end subroutine modg_voltoface_Q_y


  !> Project modal representation of an element to its two faces in Y.
  subroutine modg_volToFace_Q_z( nScalars, volState, maxPolyDegree, nElems, &
    &                            faceState )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: nScalars
    !> The modal representation in the volume. First dimension is the number of
    !! volumetric numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the faces in X direction
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: pos, facePos
    integer :: iAnsZ
    integer :: iVar
    integer :: mpd1, mpd1_square
    integer :: iElem
    integer :: iVEF, ipv
    integer :: nIndeps
    ! --------------------------------------------------------------------------

    mpd1 = maxpolydegree+1
    mpd1_square = mpd1**2
    nIndeps = nScalars*nElems*mpd1_square

    ! Split off first ansatz function, to avoid initialization with 0.
    iAnsZ = 1

    !$NEC ivdep
    do iVEF=1,nIndeps
      iVar = (iVEF-1) / (nElems*mpd1_square) + 1
      ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
      facepos = (ipv-1) / nElems + 1
      iElem = mod(ipv-1, nElems) + 1
      ! get position of the current ansatz function
      pos = facepos + (iAnsZ-1)*mpd1_square

      faceState(iElem,facePos,iVar,1) = volState(iElem, pos, iVar)
      faceState(iElem,facePos,iVar,2) = volState(iElem, pos, iVar)
    end do

    ! Legendre polynomial on the right is just the sum of all modes,
    ! but on the left we alternatingly need to add and subtract the modes.
    ! Thus, we split the loop to allow the correct operation to be used for
    ! each mode.
    do iAnsZ=3,maxPolyDegree+1,2

      !$NEC ivdep
      do iVEF=1,nIndeps
        iVar = (iVEF-1) / (nElems*mpd1_square) + 1
        ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
        facepos = (ipv-1) / nElems + 1
        iElem = mod(ipv-1, nElems) + 1
        ! get position of the current ansatz function
        pos = facepos + (iAnsZ-1)*mpd1_square

        faceState(iElem,facePos,iVar,1) = faceState(iElem,facePos,iVar,1) &
          &                             + volState(iElem, pos, iVar)
        faceState(iElem,facePos,iVar,2) = faceState(iElem,facePos,iVar,2) &
          &                             + volState(iElem, pos, iVar)
      end do

    end do
    do iAnsZ=2,maxPolyDegree+1,2

      !$NEC ivdep
      do iVEF=1,nIndeps
        iVar = (iVEF-1) / (nElems*mpd1_square) + 1
        ipv = iVEF - (iVar - 1) * (nElems*mpd1_square)
        facepos = (ipv-1) / nElems + 1
        iElem = mod(ipv-1, nElems) + 1
        ! get position of the current ansatz function
        pos = facepos + (iAnsZ-1)*mpd1_square

        faceState(iElem,facePos,iVar,1) = faceState(iElem,facePos,iVar,1) &
          &                             - volState(iElem, pos, iVar)
        faceState(iElem,facePos,iVar,2) = faceState(iElem,facePos,iVar,2) &
          &                             + volState(iElem, pos, iVar)
      end do

    end do

  end subroutine modg_voltoface_Q_z


  !> Applies the inverse of the mass matrix for a 3D scheme.
  subroutine atl_modg_invMassMatrix( mesh, kerneldata, statedata,        &
    &                                elementalTimestep, timestep, scheme )
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
    type(atl_modg_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------
    integer :: nElems

    nElems = mesh%descriptor%elem%nElems(eT_fluid)

    select case(scheme%basisType)
      case(Q_space)
        call modg_invMassMatrix_Q( mesh, kerneldata, &
          &                        scheme, nElems )
      case(P_space)
        call modg_invMassMatrix_P( mesh, kerneldata, scheme, nElems )
    end select

    ! Since this is the last part of the MOdal Discontinuous Galerkin (MODG)
    ! algorithm we can make the final timestep, now.
    call elementalTimestep( timestep, statedata%state, kerneldata )

  end subroutine atl_modg_invMassMatrix

  !> Applies the inverse of the mass matrix for a 3D scheme.
  subroutine modg_invMassMatrix_Q( mesh, kerneldata,  &
    &                              scheme, nElems )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The data of the kernel.
    type(atl_kerneldata_type) :: kerneldata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    integer, intent(in) :: nElems
    ! --------------------------------------------------------------------------
    ! Loop indices for ansatz functions
    integer :: iAnsX, iAnsY, iAnsZ
    integer :: iAnsYZ, iAnsXZ, iAnsXY
    ! Inverse of the determintant of the jacobian of the mapping from reference
    ! to physical element.
    real(kind=rk) :: inv_jacobiDet
    real(kind=rk) :: tmp(vlen)
    integer :: posCoeff
    integer :: yz_offset, xz_offset
    integer :: iVar, iElem, iStrip
    integer :: i
    integer :: nVars, nEntries
    integer :: mp1, square_mp1
    ! --------------------------------------------------------------------------

    nVars = kerneldata%nVars
    mp1 = scheme%maxPolyDegree+1
    square_mp1 = mp1**2

    inv_jacobiDet = (2.0_rk/mesh%length)**3

    ! apply the 1D inverse of the mass matrix in X
    if (mp1 > 1) then
      do iAnsYZ=1,square_mp1
        yz_offset = (iAnsYZ-1)*mp1

        do iStrip=0,nElems-1,vlen
          nEntries = min(vlen, nElems - iStrip)
          do iVar=1,nVars

            do i=1,nEntries
              tmp(i) = kerneldata%state_der(iStrip + i, 1 + yz_offset,  iVar)
              kerneldata%state_der(iStrip+i, 1+yz_offset,  iVar) &
                  &  = tmp(i) * 0.5_rk
            end do
            do iAnsX=3,mp1,2
              do i=1,nEntries
                iElem = iStrip + i
                tmp(i) = kerneldata%state_der(iElem, iAnsX+yz_offset,  iVar) &
                  &      + tmp(i)
                kerneldata%state_der(iElem, iAnsX+yz_offset,  iVar) &
                  &  = tmp(i) * 0.5_rk*(2*iAnsX-1)
              end do
            end do

            do i=1,nEntries
              tmp(i) = kerneldata%state_der(iStrip + i, 2 + yz_offset,  iVar)
              kerneldata%state_der(iStrip+i, 2+yz_offset,  iVar) &
                  &  = tmp(i) * 1.5_rk
            end do
            do iAnsX=4,mp1,2
              do i=1,nEntries
                iElem = iStrip + i
                tmp(i) = kerneldata%state_der(iElem, iAnsX+yz_offset,  iVar) &
                  &      + tmp(i)
                kerneldata%state_der(iElem, iAnsX+yz_offset,  iVar) &
                  &  = tmp(i) * 0.5_rk*(2*iAnsX-1)
              end do
            end do

          end do
        end do

      end do


      ! apply the 2D inverse of the mass matrix (in Y)
      do iAnsXZ=1,square_mp1
        iAnsX = mod(iAnsXZ-1,mp1) + 1
        xz_offset = (iAnsXZ - iAnsX)*mp1 + iAnsX

        do iStrip=0,nElems-1,vlen
          nEntries = min(vlen, nElems - iStrip)
          do iVar=1,nVars

            do i=1,nEntries
              tmp(i) = kerneldata%state_der(iStrip+i, xz_offset,  iVar)
              kerneldata%state_der(iStrip+i, xz_offset,  iVar) &
                &  = tmp(i) * 0.5_rk
            end do
            do iAnsY=3,mp1,2
              posCoeff = xz_offset + (iAnsY-1)*mp1
              do i=1,nEntries
                iElem = iStrip + i
                tmp(i) = kerneldata%state_der(iElem, posCoeff,  iVar) + tmp(i)
                kerneldata%state_der(iElem, posCoeff,  iVar) &
                  &  = tmp(i) * 0.5_rk*(2*iAnsY-1)
              end do
            end do

            do i=1,nEntries
              tmp(i) = kerneldata%state_der(iStrip+i, mp1 + xz_offset,  iVar)
              kerneldata%state_der(iStrip+i, mp1 + xz_offset,  iVar) &
                &  = tmp(i) * 1.5_rk
            end do
            do iAnsY=4,mp1,2
              posCoeff = xz_offset + (iAnsY-1)*mp1
              do i=1,nEntries
                iElem = iStrip + i
                tmp(i) = kerneldata%state_der(iElem, posCoeff,  iVar) + tmp(i)
                kerneldata%state_der(iElem, posCoeff,  iVar) &
                  &  = tmp(i) * 0.5_rk*(2*iAnsY-1)
              end do
            end do

          end do
        end do

      end do


      ! apply the 3D inverse of the mass matrix (in Z)
      do iAnsXY=1,square_mp1

        do iStrip=0,nElems-1,vlen
          nEntries = min(vlen, nElems - iStrip)
          do iVar=1,nVars

            do i=1,nEntries
              tmp(i) = kerneldata%state_der(iStrip+i, iAnsXY,  iVar)
              kerneldata%state_der(iStrip+i, iAnsXY,  iVar) &
                &  = tmp(i) * 0.5_rk * inv_jacobiDet
            end do
            do iAnsZ=3,mp1,2
              posCoeff = iAnsXY + (iAnsZ-1)*square_mp1
              do i=1,nEntries
                iElem = iStrip + i
                tmp(i) = kerneldata%state_der(iElem, posCoeff,  iVar) + tmp(i)
                kerneldata%state_der(iElem, posCoeff,  iVar) &
                  &  = tmp(i) * 0.5_rk*(2*iAnsZ-1) * inv_jacobiDet
              end do
            end do

            do i=1,nEntries
              tmp(i) = kerneldata%state_der(iStrip+i, iAnsXY+square_mp1,  iVar)
              kerneldata%state_der(iStrip+i, iAnsXY+square_mp1,  iVar) &
                &  = tmp(i) * 1.5_rk * inv_jacobiDet
            end do
            do iAnsZ=4,mp1,2
              posCoeff = iAnsXY + (iAnsZ-1)*square_mp1
              do i=1,nEntries
                iElem = iStrip + i
                tmp(i) = kerneldata%state_der(iElem, posCoeff,  iVar) + tmp(i)
                kerneldata%state_der(iElem, posCoeff,  iVar) &
                  &  = tmp(i) * 0.5_rk*(2*iAnsZ-1) * inv_jacobiDet
              end do
            end do

          end do
        end do

      end do

    else

      ! Single degree of freedom.
      do iVar=1,nVars
        kerneldata%state_der(:, 1,  iVar) &
          &  = 0.125_rk * kerneldata%state_der(:, 1, iVar)
      end do

    end if

  end subroutine modg_invMassMatrix_Q


  !> Applies the inverse of the mass matrix for a 3D scheme.
  subroutine modg_invMassMatrix_P( mesh, kerneldata, scheme, nElems )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The data of the kernel.
    type(atl_kerneldata_type) :: kerneldata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    integer, intent(in) :: nElems
    ! --------------------------------------------------------------------------
    ! Loop indices for ansatz functions
    integer :: iAnsX, iAnsY, iAnsZ
    ! Positions for the given ansatz functions
    integer :: ansLow, ansPrevLow
    ! Inverse of the determintant of the jacobian of the mapping from reference
    ! to physical element.
    real(kind=rk) :: inv_jacobiDet
    ! --------------------------------------------------------------------------

    inv_jacobiDet = (2.0_rk/mesh%length)**3

    ! apply the 1D inverse of the mass matrix
    do iAnsZ = 1, scheme%maxPolyDegree+1
      do iAnsY = 1, scheme%maxPolyDegree+1 - (iAnsZ-1)
        do iAnsX = 3, scheme%maxPolyDegree+1 - (iAnsZ-1) - (iAnsY-1)
  ! integer divisions are no mistake here.
  anslow = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
  ! integer divisions are no mistake here.
  ansprevlow = (((iansx-2 + iansy + iansz - 3) &
    &     * (iansx-2 + iansy + iansz - 2) &
    &     * (iansx-2 + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx-2 + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
          kerneldata%state_der(:nElems, ansLow, :) &
            &        = kerneldata%state_der(:nElems, ansLow, :) &
            &        + kerneldata%state_der(:nElems, ansPrevLow, :)
        end do
      end do
    end do


    do iAnsZ = 1, scheme%maxPolyDegree+1
      do iAnsY = 1, scheme%maxPolyDegree+1 - (iAnsZ-1)
        do iAnsX = 1, scheme%maxPolyDegree+1 - (iAnsZ-1) - (iAnsY-1)
  ! integer divisions are no mistake here.
  anslow = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
          kerneldata%state_der(:nElems, ansLow,  :) &
            &    = kerneldata%state_der(:nElems, ansLow,  :) &
            &      * 0.5_rk*(2*iAnsX-1)
        end do
      end do
    end do


    ! apply the 2D inverse of the mass matrix
    do iAnsX = 1, scheme%maxPolyDegree+1
      do iAnsZ = 1, scheme%maxPolyDegree+1 - (iAnsX-1)
        do iAnsY = 3, scheme%maxPolyDegree+1 - (iAnsX-1) - (iAnsZ-1)
  ! integer divisions are no mistake here.
  anslow = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
  ! integer divisions are no mistake here.
  ansprevlow = (((iansx + iansy-2 + iansz - 3) &
    &     * (iansx + iansy-2 + iansz - 2) &
    &     * (iansx + iansy-2 + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy-2 + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-2-1)
          kerneldata%state_der(:nElems, ansLow,  :) &
            &    = kerneldata%state_der(:nElems, ansLow,  :) &
            &    + kerneldata%state_der(:nElems, ansPrevLow,  :)
        end do
      end do
    end do


    do iAnsZ = 1, scheme%maxPolyDegree+1
      do iAnsY = 1, scheme%maxPolyDegree+1 - (iAnsZ-1)
        do iAnsX = 1, scheme%maxPolyDegree+1 - (iAnsZ-1) - (iAnsY-1)
  ! integer divisions are no mistake here.
  anslow = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
          kerneldata%state_der(:nElems, ansLow,  :) &
            &    = kerneldata%state_der(:nElems, ansLow,  :) &
            &      * 0.5_rk*(2*iAnsY-1)
        end do
      end do
    end do


    ! apply the 3D inverse of the mass matrix
    do iAnsX = 1, scheme%maxPolyDegree+1
      do iAnsY = 1, scheme%maxPolyDegree+1 - (iAnsX-1)
        do iAnsZ = 3, scheme%maxPolyDegree+1 - (iAnsX-1) - (iAnsY-1)
  ! integer divisions are no mistake here.
  anslow = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
  ! integer divisions are no mistake here.
  ansprevlow = (((iansx + iansy + iansz-2 - 3) &
    &     * (iansx + iansy + iansz-2 - 2) &
    &     * (iansx + iansy + iansz-2 - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-2-1) * (iansx + iansy + iansz-2 -2) &
    &   - ((iansz-2-2) * (iansz-2-1)) / 2) &
    & + (iansy-1)
          kerneldata%state_der(:nElems, ansLow,:) &
            &        = kerneldata%state_der(:nElems, ansLow,  :) &
            &        + kerneldata%state_der(:nElems, ansPrevLow,  :)
        end do
      end do
    end do


    ! Apply also inverse of determinant of the jacobain of the mapping from
    ! reference to physical element.
    do iAnsZ = 1, scheme%maxPolyDegree+1
      do iAnsY = 1, scheme%maxPolyDegree+1 - (iAnsZ-1)
        do iAnsX = 1, scheme%maxPolyDegree+1 - (iAnsZ-1) - (iAnsY-1)
  ! integer divisions are no mistake here.
  anslow = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
          kerneldata%state_der(:nElems, ansLow,  :) &
            &    = kerneldata%state_der(:nElems, ansLow,  :) &
            &      * 0.5_rk*(2*iAnsZ-1) * inv_jacobiDet
        end do
      end do
    end do


  end subroutine modg_invMassMatrix_P


  !> Applies a scaled transposed inverse of the mass matrix for a 3D scheme.
  !! The result is a transformation of the polynomial basis from
  !! the ansatz- to the test-polynomials.
  !! This is useful for a local predictor
  subroutine atl_modg_scaledTransposedInvMassMatrix_Q( nTotalElems, nDofs, &
    &                                                  nScalars,           &
    &                                                  maxPolyDegree,      &
    &                                                  nElems, state       )
    ! --------------------------------------------------------------------------
    !> defines the dimensions of the state array
    integer, intent(in) :: nTotalElems, nDofs, nScalars
    !> polynomial degree of the modal dg scheme
    integer, intent(in) :: maxPolyDegree
    !> number of elements to work with
    integer, intent(in) :: nElems
    !> the state to transform
    real(kind=rk), intent(inout) :: state(nTotalElems,nDofs,nScalars)
    ! --------------------------------------------------------------------------
    ! Loop indices for ansatz functions
    integer :: iAnsX, iAnsY, iAnsZ
    integer :: ij, ik, jk
    ! Positions for the given ansatz functions
    integer :: ans, ansNext
    integer :: mpd1, mpd1_square
    ! --------------------------------------------------------------------------

    mpd1 = maxPolyDegree+1
    mpd1_square = mpd1**2

    ! apply the 1D inverse of the mass matrix
    do ij=1,mpd1_square
      do iAnsZ = maxPolyDegree-1, 1, -1
        ans = (iAnsZ-1)*mpd1_square + ij
        ansNext = ans + 2*mpd1_square
        state(:nElems, ans, :) &
          &        = state(:nElems, ans, :) &
          &        + state(:nElems, ansNext, :)
      end do
    end do

    ! apply the 2D inverse of the mass matrix
    do ik=1,mpd1_square
      iAnsX = mod(ik-1, mpd1) + 1
      ! iAnsZ = (ik-1)/mpd1 + 1
      iAnsZ = (ik-1)/mpd1 * mpd1
      do iAnsY = maxPolyDegree-1, 1, -1
        ans = (iAnsZ+(iAnsY-1))*mpd1 + iAnsX
        ansNext = ans + 2*mpd1
        state(:nElems, ans, :) &
          &    = state(:nElems, ans, :) &
          &    + state(:nElems, ansNext, :)
      end do
    end do

    ! apply the 3D inverse of the mass matrix
    do jk=1,mpd1_square
      do iAnsX = maxPolyDegree-1, 1, -1
        ans = (jk-1)*mpd1 + iAnsX
        ansNext = ans + 2
        state(:nElems, ans, :) &
          &        = state(:nElems, ans, :) &
          &        + state(:nElems, ansNext, :)
      end do
    end do


  end subroutine atl_modg_scaledTransposedInvMassMatrix_Q


  !> Applies a scaled transposed inverse of the mass matrix for a 3D scheme.
  !! The result is a transformation of the polynomial basis from the
  !! ansatz- to the test-polynomials.
  !! This is useful for a local predictor
  subroutine atl_modg_scaledTransposedInvMassMatrix_P( nTotalElems, nDofs, &
    &                                                  nScalars,           &
    &                                                  maxPolyDegree,      &
    &                                                  nElems, state       )
    ! --------------------------------------------------------------------------
    !> defines the dimensions of the state array
    integer, intent(in) :: nTotalElems, nDofs, nScalars
    !> polynomial degree of the modal dg scheme
    integer, intent(in) :: maxPolyDegree
    !> number of elements to work with
    integer, intent(in) :: nElems
    !> the state to transform
    real(kind=rk), intent(inout) :: state(nTotalElems,nDofs,nScalars)
    ! --------------------------------------------------------------------------
    ! Loop indices for ansatz functions
    integer :: iAnsX, iAnsY, iAnsZ
    ! Positions for the given ansatz functions
    integer :: ans, ansNext
    ! --------------------------------------------------------------------------


    ! apply the 1D inverse of the mass matrix
    do iAnsX = 1, maxPolyDegree+1, 1
      do iAnsY = 1, maxPolyDegree+1 - (iAnsX-1), 1
        do iAnsZ = maxPolyDegree-1 - (iAnsX-1) - (iAnsY-1), 1, -1
  ! integer divisions are no mistake here.
  ans = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
  ! integer divisions are no mistake here.
  ansnext = (((iansx + iansy + iansz+2 - 3) &
    &     * (iansx + iansy + iansz+2 - 2) &
    &     * (iansx + iansy + iansz+2 - 1)) &
    &   / 6 + 1)             &
    & + ((iansz+2-1) * (iansx + iansy + iansz+2 -2) &
    &   - ((iansz+2-2) * (iansz+2-1)) / 2) &
    & + (iansy-1)
          state(:nElems, ans, :) &
            &        = state(:nElems, ans, :) &
            &        + state(:nElems, ansNext, :)
        end do
      end do
    end do


    ! apply the 2D inverse of the mass matrix
    do iAnsX = 1, maxPolyDegree+1, 1
      do iAnsZ = 1, maxPolyDegree+1 - (iAnsX-1), 1
        do iAnsY = maxPolyDegree-1 - (iAnsX-1) - (iAnsZ-1), 1, -1
  ! integer divisions are no mistake here.
  ans = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
  ! integer divisions are no mistake here.
  ansnext = (((iansx + iansy+2 + iansz - 3) &
    &     * (iansx + iansy+2 + iansz - 2) &
    &     * (iansx + iansy+2 + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy+2 + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy+2-1)
          state(:nElems, ans, :) &
            &    = state(:nElems, ans, :) &
            &    + state(:nElems, ansNext, :)
        end do
      end do
    end do


    ! apply the 3D inverse of the mass matrix
    do iAnsZ = 1, maxPolyDegree+1, 1
      do iAnsY = 1, maxPolyDegree+1 - (iAnsZ-1), 1
        do iAnsX = maxPolyDegree-1 - (iAnsZ-1) - (iAnsY-1), 1, -1
  ! integer divisions are no mistake here.
  ans = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
  ! integer divisions are no mistake here.
  ansnext = (((iansx+2 + iansy + iansz - 3) &
    &     * (iansx+2 + iansy + iansz - 2) &
    &     * (iansx+2 + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx+2 + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
          state(:nElems, ans, :) &
            &        = state(:nElems, ans, :) &
            &        + state(:nElems, ansNext, :)
        end do
      end do
    end do


  end subroutine atl_modg_scaledTransposedInvMassMatrix_P


  !> Projection of the physical flux for the local predictor:
  !! This function expects the physical flux transformed in the basis of
  !! test-polynomials and projects it onto the ansatz-functions.
  !! Thus, the result is directly given in the basis of ansatz-polynomials
  !! (no invMassMatrix-call needed afterwards)
  !!
  !! This function implements the multiplication with the transposed stiffness
  !! matrix and some useful scaling factors
  subroutine atl_modg_scaledTransposedProject_physFlux_Q( nScalars, &
    & maxPolyDegree, length, nElems, state, iDir, dirVec            )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The element index to project for
    integer, intent(in) :: nElems
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:,:)
    integer, intent(in) :: iDir
    !> ordering of xyz for current direction
    integer, intent(in) :: dirVec(3)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: scaledJacobiDetStiffProj, scalProd1(maxPolyDegree)
    integer :: testPos
    integer :: iTest2, iTest1, iTest3, iTestVec(3)
    integer :: iAnsVec(3)
    integer :: iVar
    integer :: ansPos(4)
    real(kind=rk) :: scalProd
    integer :: ij
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes in test basis
    ! onto the ansatz functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    !jacobiDetStiffProj = (0.5_rk*length)**2
    ! Here we apply some additional scaling with invJacobiDet = (2/length)**3,
    ! so the final scaling is (- because it's on the rhs)
    scaledJacobiDetStiffProj = - (2/length)


    ! we directly scale the result with the inverse of the diagonal matrix
    ! <ansatz_i,ansatz_j>
    do iTest1=1,maxPolyDegree
      scalProd1(iTest1) = (2.0_rk*iTest1-1)*scaledJacobiDetStiffProj
    end do


    ! unrolled loop

    do iTest3 = maxPolyDegree, maxPolyDegree+1
      do ij=1,2*maxPolyDegree
        iTest2 = (ij-1)/maxPolyDegree + maxPolyDegree
        iTest1 = mod(ij-1, maxPolyDegree) + 1

        ! one entry

        iTestVec = (/iTest1, iTest2, iTest3/)
  testpos = itestvec(dirvec(1))                                      &
    &      + ( ( itestvec(dirvec(2))-1)                             &
    &      + (itestvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        iAnsVec = (/iTest1+1, iTest2, iTest3/)
  anspos(1) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        scalProd = scalProd1(iTest1)

        do iVar=1,nScalars
          state(:nElems,testPos,1,iVar) &
            & = state(:nElems,testPos,1,iVar) &
            & + state(:nElems,ansPos(1),iDir+1,iVar) * scalProd
        end do

      end do

      do ij=1,maxPolyDegree*(maxPolyDegree-1)
        iTest2 = (ij-1)/maxPolyDegree + 1
        iTest1 = mod(ij-1, maxPolyDegree) + 1

        ! two entries

        iTestVec = (/iTest1, iTest2, iTest3/)
  testpos = itestvec(dirvec(1))                                      &
    &      + ( ( itestvec(dirvec(2))-1)                             &
    &      + (itestvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAnsVec = (/iTest1+1, iTest2, iTest3/)
  anspos(1) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        iAnsVec = (/iTest1+1, iTest2+2, iTest3/)
  anspos(2) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        scalProd = scalProd1(iTest1)

        do iVar=1,nScalars
          state(:nElems,testPos,1,iVar) &
            & = state(:nElems,testPos,1,iVar) &
            & + ( state(:nElems,ansPos(1),iDir+1,iVar) &
            &   - state(:nElems,ansPos(2),iDir+1,iVar) ) * scalProd
        end do

      end do
    end do

    do iTest3 = 1, maxPolyDegree-1
      do ij=1,2*maxPolyDegree
        iTest2 = (ij-1)/maxPolyDegree + maxPolyDegree
        iTest1 = mod(ij-1, maxPolyDegree) + 1

        ! two entries

        iTestVec = (/iTest1, iTest2, iTest3/)
  testpos = itestvec(dirvec(1))                                      &
    &      + ( ( itestvec(dirvec(2))-1)                             &
    &      + (itestvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAnsVec = (/iTest1+1, iTest2, iTest3/)
  anspos(1) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        iAnsVec = (/iTest1+1, iTest2, iTest3+2/)
  anspos(2) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        scalProd = scalProd1(iTest1)

        do iVar=1,nScalars
          state(:nElems,testPos,1,iVar) &
            & = state(:nElems,testPos,1,iVar) &
            & + ( state(:nElems,ansPos(1),iDir+1,iVar) &
            &   - state(:nElems,ansPos(2),iDir+1,iVar) ) * scalProd
        end do

      end do


      do ij=1,maxPolyDegree*(maxPolyDegree-1)
        iTest2 = (ij-1)/maxPolyDegree + 1
        iTest1 = mod(ij-1, maxPolyDegree) + 1

        ! four entries

        iTestVec = (/iTest1, iTest2, iTest3/)
  testpos = itestvec(dirvec(1))                                      &
    &      + ( ( itestvec(dirvec(2))-1)                             &
    &      + (itestvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAnsVec = (/iTest1+1, iTest2, iTest3/)
  anspos(1) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        iAnsVec = (/iTest1+1, iTest2+2, iTest3/)
  anspos(2) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        iAnsVec = (/iTest1+1, iTest2+2, iTest3+2/)
  anspos(3) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        iAnsVec = (/iTest1+1, iTest2, iTest3+2/)
  anspos(4) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        scalProd = scalProd1(iTest1)


        do iVar=1,nScalars
          state(:nElems,testPos,1,iVar) &
            & = state(:nElems,testPos,1,iVar) &
            & + ( state(:nElems,ansPos(1),iDir+1,iVar) &
            &   - state(:nElems,ansPos(2),iDir+1,iVar) &
            &   + state(:nElems,ansPos(3),iDir+1,iVar) &
            &   - state(:nElems,ansPos(4),iDir+1,iVar) ) * scalProd
        end do

      end do
    end do

  end subroutine atl_modg_scaledTransposedProject_physFlux_Q


  !> Projection of the physical flux for the local predictor:
  !! This function expects the physical flux transformed in the basis of
  !! test-polynomials and projects it onto the ansatz-functions.
  !! Thus, the result is directly given in the basis of ansatz-polynomials
  !! (no invMassMatrix-call needed afterwards)
  !!
  !! This function implements the multiplication with the transposed stiffness
  !! matrix and some useful scaling factors
  subroutine atl_modg_scaledTransposedProject_physFlux_P( nScalars, &
    & maxPolyDegree, length, nElems, state, iDir, dirVec            )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The element index to project for
    integer, intent(in) :: nElems
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:,:)
    integer, intent(in) :: iDir
    !> ordering of xyz for current direction
    integer, intent(in) :: dirVec(3)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: scaledJacobiDetStiffProj, scalProd1(maxPolyDegree)
    integer :: testPos
    integer :: iTest2, iTest1, iTest3, iTestVec(3)
    integer :: iAnsVec(3)
    integer :: iVar
    integer :: ansPos(4)
    real(kind=rk) :: scalProd
    integer :: maxDeg2, maxDeg1
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes in test basis
    ! onto the ansatz functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    !jacobiDetStiffProj = (0.5_rk*length)**2
    ! Here we apply some additional scaling with invJacobiDet = (2/length)**3,
    ! so the final scaling is (- because it's on the rhs)
    scaledJacobiDetStiffProj = - (2/length)


    ! we directly scale the result with the inverse of the diagonal matrix
    ! <ansatz_i,ansatz_j>
    do iTest1=1,maxPolyDegree
      scalProd1(iTest1) = (2.0_rk*iTest1-1)*scaledJacobiDetStiffProj
    end do


    ! unrolled loop

    do iTest3 = 1, maxPolyDegree+1,1
      maxDeg2 = max(1,maxPolyDegree - (iTest3-1))

      do iTest2 = 1, maxDeg2+1, 1
        maxDeg1 = maxPolyDegree - (iTest3-1) - (iTest2-1)

        do iTest1 = max(1,maxDeg1-1), maxDeg1, 1

          ! one entry

          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)

          iAnsVec = (/iTest1+1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd =   scalProd1(iTest1)

          do iVar=1,nScalars
            state(:nElems,testPos,1,iVar) &
              & = state(:nElems,testPos,1,iVar) &
              & + state(:nElems,ansPos(1),iDir+1,iVar) * scalProd
          end do

        end do
      end do
    end do


    do iTest3 = maxPolyDegree, maxPolyDegree+1
      maxDeg2 = max(1,maxPolyDegree - (iTest3-1))

      do iTest2 = 1, maxDeg2-1, 1
        maxDeg1 = maxPolyDegree - (iTest3-1) - (iTest2+1)

        do iTest1 = max(1, maxDeg1-1), maxDeg1, 1

          ! two entries

          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)


          iAnsVec = (/iTest1+1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)

          iAnsVec = (/iTest1+1, iTest2+2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(2) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd = scalProd1(iTest1)

          do iVar=1,nScalars
            state(:nElems,testPos,1,iVar) &
              & = state(:nElems,testPos,1,iVar) &
              & + ( state(:nElems,ansPos(1),iDir+1,iVar) &
              &   - state(:nElems,ansPos(2),iDir+1,iVar) ) * scalProd
          end do

        end do
      end do
    end do

    do iTest3 = 1, maxPolyDegree-1, 1
      maxDeg2 = max(1,maxPolyDegree - (iTest3-1))

      do iTest2 = max(1,maxDeg2), maxDeg2+1, 1
        maxDeg1 = maxPolyDegree - (iTest3+1) - (iTest2-1)

        do iTest1 = max(1,maxDeg1-1), maxDeg1, 1

          ! two entries

          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)


          iAnsVec = (/iTest1+1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)

          iAnsVec = (/iTest1+1, iTest2, iTest3+2/)
  ! integer divisions are no mistake here.
  anspos(2) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd =   scalProd1(iTest1)


          do iVar=1,nScalars
            state(:nElems,testPos,1,iVar) &
              & = state(:nElems,testPos,1,iVar) &
              & + ( state(:nElems,ansPos(1),iDir+1,iVar) &
              &   - state(:nElems,ansPos(2),iDir+1,iVar) ) * scalProd
          end do

        end do
      end do


      do iTest2 = 1, maxDeg2-1, 1
        maxDeg1 = maxPolyDegree - (iTest3-1) - (iTest2-1) - 2

        do iTest1 = max(1,maxDeg1-1), maxDeg1, 1

          ! three entries


          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)


          iAnsVec = (/iTest1+1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)

          iAnsVec = (/iTest1+1, iTest2, iTest3+2/)
  ! integer divisions are no mistake here.
  anspos(2) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)

          iAnsVec = (/iTest1+1, iTest2+2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(3) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd = scalProd1(iTest1)



          do iVar=1,nScalars
            state(:nElems,testPos,1,iVar) &
              & = state(:nElems,testPos,1,iVar) &
              & + ( state(:nElems,ansPos(1),iDir+1,iVar) &
              &   - state(:nElems,ansPos(2),iDir+1,iVar) &
              &   - state(:nElems,ansPos(3),iDir+1,iVar) ) * scalProd
          end do

        end do
      end do


      do iTest2 = 1, maxDeg2-1, 1
        maxDeg1 = maxPolyDegree - (iTest3+1) - (iTest2+1)

        do iTest1 = 1, maxDeg1, 1

          ! four entries

          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)


          iAnsVec = (/iTest1+1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)

          iAnsVec = (/iTest1+1, iTest2, iTest3+2/)
  ! integer divisions are no mistake here.
  anspos(2) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)

          iAnsVec = (/iTest1+1, iTest2+2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(3) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)

          iAnsVec = (/iTest1+1, iTest2+2, iTest3+2/)
  ! integer divisions are no mistake here.
  anspos(4) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd = scalProd1(iTest1)



          do iVar=1,nScalars
            state(:nElems,testPos,1,iVar) &
              & = state(:nElems,testPos,1,iVar) &
              & + ( state(:nElems,ansPos(1),iDir+1,iVar) &
              &   - state(:nElems,ansPos(2),iDir+1,iVar) &
              &   - state(:nElems,ansPos(3),iDir+1,iVar) &
              &   + state(:nElems,ansPos(4),iDir+1,iVar) ) * scalProd
          end do

        end do
      end do
    end do

  end subroutine atl_modg_scaledTransposedProject_physFlux_P


  !> Subroutine to test the various routines of this module.
  subroutine atl_modg_kernel_utests(passed)
    use atl_equation_module, only: atl_equations_type
    use atl_eqn_euler_derive_module, only: atl_eqn_euler_cons2prim_elems, &
      &                                    atl_eqn_euler_prim2cons_elems
    use atl_eqn_euler_var_module, only: atl_init_euler_vars
    use atl_varsys_module, only: atl_varsys_solverdata_type
    ! -------------------------------------------------------------------- !
    logical, intent(out) :: passed
    ! -------------------------------------------------------------------- !
    type(atl_equations_type) :: eqn
    type(atl_varSys_solverData_type) :: varsys_data
    integer :: iDir
    integer :: order
    real(kind=rk), allocatable :: resmom(:,:)
    ! -------------------------------------------------------------------- !

    eqn%eq_kind = 'navier_stokes'
    eqn%isNonlinear = .true.
    eqn%nDerivatives =  1
    eqn%adaptive_timestep = .true.

    eqn%cons2prim => atl_eqn_euler_cons2prim_elems
    eqn%prim2cons => atl_eqn_euler_prim2cons_elems

    call atl_init_euler_vars(    &
      & equation   = eqn,        &
      & solverData = varSys_data )

    eqn%euler%isen_coef = 1.4_rk
    eqn%euler%r = 1.0_rk
    eqn%euler%cv = eqn%euler%r / (eqn%euler%isen_coef - 1.0_rk)
    eqn%euler%porosity = 0.0_rk
    eqn%euler%viscous_permeability = 0.0_rk
    eqn%euler%thermal_permeability = 0.0_rk

    eqn%navierstokes%mu = 1.0e-3_rk
    eqn%navierStokes%lambda = 2.0_rk * eqn%navierstokes%mu / 3.0_rk
    eqn%navierstokes%therm_cond = 0.7_rk
    eqn%navierstokes%ip_param = 0.5_rk

    order = 5
    allocate(resmom(order**3,3))
    do iDir = 1,3
      write(*,*) '==========================================================='
      call test_project_stabViscNumFlux( polydegree       = order-1,       &
        &                                oversamplefactor = 1.0_rk,        &
        &                                dir              = iDir,          &
        &                                equation         = eqn,           &
        &                                length           = 1.0_rk,        &
        &                                rotated_mom      = resmom(:,iDir) )
      write(*,*) '==========================================================='
      write(*,*) ''
    end do

    passed = ( (maxval(abs(resmom(:,2) - resmom(:,1)))      &
      &         < epsilon(1.0_rk))                          &
      &       .and. (maxval(abs(resmom(:,3) - resmom(:,1))) &
      &              < epsilon(1.0_rk))                     )

  end subroutine atl_modg_kernel_utests


  subroutine test_project_stabViscNumFlux(polydegree, oversamplefactor, dir, &
    &                                     equation, length, rotated_mom      )
    use ply_dof_module, only: q_space
    use ply_poly_project_module, only: ply_poly_project_fillbody
    use ply_dynArray_project_module, only: ply_prj_init_type
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: polydegree
    real(kind=rk), intent(in) :: oversamplefactor
    integer, intent(in) :: dir
    type(atl_equations_type), intent(in) :: equation
    real(kind=rk), intent(in) :: length
    real(kind=rk), intent(out) :: rotated_mom(:)
    ! -------------------------------------------------------------------- !
    integer :: nScalars
    integer :: idof, jdof, kdof
    integer :: signfact
    integer :: rotdof
    integer :: dofstart
    real(kind=rk), allocatable :: numFlux(:,:,:,:)
    real(kind=rk), allocatable :: faceState(:,:,:,:)
    real(kind=rk), allocatable :: projection(:,:,:)
    type(ply_poly_project_type) :: poly_proj
    type(ply_prj_init_type) :: proj_init
    ! -------------------------------------------------------------------- !

    nScalars = equation%varSys%nScalars

    ! Just use l2p here to not depend on fftw for this test.
    proj_init%header%kind = 'l2p'
    proj_init%maxpolydegree = polydegree
    proj_init%basisType = q_space
    proj_init%header%l2p_header%factor = oversamplefactor
    proj_init%header%l2p_header%nodes_header%nodes_kind = 'gauss-legendre'

    call ply_poly_project_fillbody( me         = poly_proj, &
      &                             proj_init  = proj_init, &
      &                             scheme_dim = 3          )

    ! Flux for:
    ! 1 element,
    ! ndofs degrees of freedom,
    ! 2*nScalars variables,
    ! 2 sides
    allocate(numflux(1, poly_proj%body_2D%ndofs, 2*nScalars, 2))

    ! State on the face for:
    ! 1 element,
    ! ndofs degrees of freedom,
    ! 2*nScalars variables,
    ! 2 sides
    allocate(facestate(1, poly_proj%body_2D%ndofs, nScalars, 2))

    ! Constant state with a velocity in the current direction.
    facestate = 0.0_rk
    facestate(:, 1, [1,dir+1,5], :) = 1.0_rk


    signfact = 1
    do jdof=0,polydegree
      do idof=0,polydegree
        dofstart = jdof*(polydegree+1) + iDof + 1
        numflux(:,dofstart,nScalars+1:2*nScalars,:) = signfact * 1.0_rk/dofstart
        signfact = -signfact
      end do
    end do

    numflux(:,:,nScalars+1:2*nScalars,:) &
      &  = numflux(:,:,nScalars+1:2*nScalars,:) + facestate

    ! Projection onto testfunctions for:
    ! 1 element,
    ! ndofs degrees of freedom,
    ! nScalars variables
    allocate(projection(1, poly_proj%body_3D%ndofs, nScalars))

    projection = 0.0_rk

    select case(dir)
    case(1)
      call modg_project_stabViscNumFluxX_Q( &
        &    numFlux       = numFlux,       &
        &    faceState     = faceState,     &
        &    equation      = equation,      &
        &    maxpolydegree = polydegree,    &
        &    length        = length,        &
        &    nElems_fluid  = 1,             &
        &    projection    = projection,    &
        &    poly_proj     = poly_proj      )
      rotated_mom = projection(1,:,2)

    case(2)
      call modg_project_stabViscNumFluxY_Q( &
        &    numFlux       = numFlux,       &
        &    faceState     = faceState,     &
        &    equation      = equation,      &
        &    maxpolydegree = polydegree,    &
        &    length        = length,        &
        &    nElems_fluid  = 1,             &
        &    projection    = projection,    &
        &    poly_proj     = poly_proj      )
      rotdof = 0
      do iDof=0,polydegree
        do kDof=0,polydegree
          do jDof=0,polydegree
            rotdof = rotdof + 1
            dofstart = 1 + kdof + ( idof*(polydegree+1) &
              &                     + jDof ) * (polydegree+1)
            rotated_mom(rotdof) = projection(1,dofstart,3)
          end do
        end do
      end do

    case(3)
      call modg_project_stabViscNumFluxZ_Q( &
        &    numFlux       = numFlux,       &
        &    faceState     = faceState,     &
        &    equation      = equation,      &
        &    maxpolydegree = polydegree,    &
        &    length        = length,        &
        &    nElems_fluid  = 1,             &
        &    projection    = projection,    &
        &    poly_proj     = poly_proj      )
      rotdof = 0
      do jDof=0,polydegree
        do iDof=0,polydegree
          do kDof=0,polydegree
            rotdof = rotdof + 1
            dofstart = 1 + idof + ( kdof*(polydegree+1) &
              &                     + jDof ) * (polydegree+1)
            rotated_mom(rotdof) = projection(1,dofstart,4)
          end do
        end do
      end do

    end select

    do kdof=0,polydegree
      write(*,*) '----------------------------------------------------------'
      do jdof=0,polydegree
        dofstart = (kdof*(polydegree+1) + jdof)*(polydegree+1) + 1
        write(*,*) rotated_mom(dofstart:dofstart+polydegree)
      end do
    end do
    write(*,*) '----------------------------------------------------------'

  end subroutine test_project_stabViscNumFlux


end module atl_modg_kernel_module
