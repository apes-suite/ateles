! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
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

module atl_modg_maxwell_kernel_module
  use env_module,               only: rk
  use tem_aux_module,           only: tem_abort

  use ply_dof_module,           only: Q_space, &
    &                                 P_space
  use ply_poly_project_module,  only: ply_poly_project_type, assignment(=), &
   &                                  ply_poly_project_n2m, ply_poly_project_m2n
  use ply_oversample_module,    only: ply_convert2oversample, &
   &                                  ply_convertFromoversample

  use atl_equation_module,      only: atl_equations_type
  use atl_maxwell_flux_module,  only: atl_maxwell_flux
  use atl_modg_scheme_module,   only: atl_modg_scheme_type
  use atl_facedata_module,      only: atl_facedata_type
  use atl_materialPrp_module,   only: atl_material_type
  use atl_penalization_module,  only: atl_penalizationData_type
  use atl_scheme_module,        only: atl_scheme_type

  implicit none

  private

  public :: atl_modg_maxwell_numFlux,  &
    & atl_modg_maxwell_physFlux_const, &
    & atl_modg_maxwell_physFlux_NonConst


contains


  ! ************************************************************************ !
  !> Calculate the numerical flux for Maxwell equation and MODG scheme
  subroutine atl_modg_maxwell_numFlux( equation, facedata, scheme, poly_proj, &
    &                                  material )
    ! -------------------------------------------------------------------- !
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(inout) :: scheme
    !> Data for projection method
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! -------------------------------------------------------------------- !
    integer :: iDir, nFaceDofs
    real(kind=rk), allocatable :: modalCoeffs(:,:,:)
    real(kind=rk), allocatable :: pntVal(:,:,:)
    real(kind=rk), allocatable :: nodalNumFlux(:,:)
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    integer :: nquadpoints, oversamp_dofs
    ! -------------------------------------------------------------------- !

    ! Numerical flux for faces in all 3 spatial face directions
    select case(scheme%basisType)
      case(Q_space)
        nFaceDofs = (scheme%maxPolyDegree+1)**2
      case(P_space)
  nfacedofs = ((scheme%maxpolydegree)+1)*((scheme%maxpolydegree)+2)/2
    end select

    ! correct amount of points due to projection
    nquadpoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs


    ! Calculate the numerical fluxes for the faces with constant material
    ! paramters
    do iDir = 1,3
      if(size(material%material_desc%computeFace(iDir,1)%leftPos) > 0) then
        call atl_maxwell_flux(                                   &
          & nTotalFaces    = size(facedata%faceRep(iDir)%dat,1), &
          & nSides         = size(material%material_desc         &
          &                               %computeFace(iDir,1)   &
          &                               %leftPos),             &
          & nFaceDofs      = nFaceDofs,                          &
          & faceRep        = facedata%faceRep(iDir)%dat,         &
          & faceFlux       = facedata%faceFlux(iDir)%dat,        &
          & leftPos        = material%material_desc              &
          &                          %computeFace(iDir,1)        &
          &                          %leftPos,                   &
          & rightPos       = material%material_desc              &
          &                          %computeFace(iDir,1)        &
          &                          %rightPos,                  &
          & var            = equation%varRotation(iDir)          &
          &                          %varTransformIndices(:),    &
          & material_left  = material%material_dat               &
          &                          %faceMaterialData(iDir,1)   &
          &                          %leftElemMaterialDat,       &
          & material_right = material%material_dat               &
          &                          %faceMaterialData(iDir,1)   &
          &                          %rightElemMaterialDat       )
      end if
    end do


    ! Calculate the numerical fluxes for the faces with non-constant material
    ! paramters
    do iDir = 1,3
      if(size(material%material_desc%computeFace(iDir,2)%leftPos) > 0) then
        ! Variable material parameters work only for Q_space polynomials
        if(.not.(scheme%basisType.eq.Q_space)) then
          call tem_abort( 'Error in atl_modg_maxwell_numFlux: Variable'  &
            & // ' material parameters work only for Q polynomial space' )
        end if

        allocate( modalCoeffs(oversamp_dofs,equation%varSys%nScalars,2) )
        allocate( pntVal(nQuadPoints,equation%varSys%nScalars,2) )
        allocate( nodalNumFlux(nQuadPoints,equation%varSys%nScalars) )
        allocate( numFluxBuffer(oversamp_dofs,equation%varSys%nScalars) )

        call atl_maxwell_flux(                                   &
          & nTotalFaces    = size(facedata%faceRep(iDir)%dat,1), &
          & nSides         = size(material%material_desc         &
          &                               %computeFace(iDir,2)   &
          &                               %leftPos),             &
          & nFaceDofs      = nFaceDofs,                          &
          & faceRep        = facedata%faceRep(iDir)%dat,         &
          & faceFlux       = facedata%faceFlux(iDir)%dat,        &
          & leftPos        = material%material_desc              &
          &                          %computeFace(iDir,2)        &
          &                          %leftPos,                   &
          & rightPos       = material%material_desc              &
          &                          %computeFace(iDir,2)        &
          &                          %rightPos,                  &
          & var            = equation%varRotation(iDir)          &
          &                          %varTransformIndices(:),    &
          & material_left  = material%material_dat               &
          &                          %faceMaterialData(iDir,2)   &
          &                          %leftElemMaterialDat,       &
          & material_right = material%material_dat               &
          &                          %faceMaterialData(iDir,2)   &
          &                          %rightElemMaterialDat,      &
          & poly_proj      = poly_proj,                          &
          & modalCoeffs    = modalCoeffs,                        &
          & pntVal         = pntVal,                             &
          & nodalNumFlux   = nodalNumFlux,                       &
          & numFluxBuffer  = numFluxBuffer                       )

        deallocate( modalCoeffs )
        deallocate( pntVal )
        deallocate( nodalNumFlux )
        deallocate( numFluxBuffer )

      end if
    end do


  end subroutine atl_modg_maxwell_numFlux
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_modg_maxwell_physFlux_const ( equation, res, state, iElem,   &
    & iDir, penalizationData, poly_proj, material, nodal_data,nodal_gradData, &
    & nodal_res, elemLength,  scheme_min, scheme_current                      )
    ! -------------------------------------------------------------------- !
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> To store the resulting phy flux in modal form
    real(kind=rk), intent(inout) :: res(:,:)
    !> The state of the equation
    real(kind=rk), intent(in), optional :: state(:,:)
    !> The current Element
    integer, intent(in) :: iElem
    !> The current Direction
    integer, intent(in) :: iDir
    !> The penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> Poly project
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    !> The state in nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    real(kind=rk), intent(in), optional :: nodal_GradData(:,:,:)
    real(kind=rk), intent(inout) :: nodal_res(:,:)
    !> Length of the element
    real(kind=rk), intent(in) :: ElemLength
    !> The scheme information
    type(atl_scheme_type), intent(inout) :: scheme_min
    type(atl_scheme_type), intent(inout) :: scheme_current
    ! -------------------------------------------------------------------- !
    ! Rotation indices for physical flux calculation in y and z direction
    integer :: rot(6)
    integer :: nDofs, nScalars
    real(kind=rk) :: inv_mu, inv_epsi
    ! -------------------------------------------------------------------- !

    ! get the rotation for the physical flux calculation in y and z direction.
    rot = equation%varRotation(iDir)%varTransformIndices(1:6)
    nDofs = size(state,1)
    nScalars = size(state,2)

    inv_mu = 1.0_rk / material%material_dat%elemMaterialData(1)  &
      &                                    %materialDat(iElem,1,1)
    inv_epsi = 1.0_rk / material%material_dat%elemMaterialData(1)  &
      &                                      %materialDat(iElem,1,2)

    call compute_physFlux( nDofs     = nDofs,         &
      &                    nScalars  = nScalars,      &
      &                    state_der = res(:ndofs,:), &
      &                    state     = state,         &
      &                    rot       = rot,           &
      &                    inv_mu    = inv_mu,        &
      &                    inv_epsi  = inv_epsi       )

  end subroutine atl_modg_maxwell_physFlux_const
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_modg_maxwell_physFlux_NonConst( equation, res, state, iElem,  &
    & iDir, penalizationData, poly_proj, material, nodal_data, nodal_gradData, &
    & nodal_res, elemLength,  scheme_min, scheme_current                       )
    ! -------------------------------------------------------------------- !
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> To store the resulting phy flux in modal form
    real(kind=rk), intent(inout) :: res(:,:)
    !> The state of the equation
    real(kind=rk), intent(in), optional :: state(:,:)
    !> The current Element
    integer, intent(in) :: iElem
    !> The current Direction
    integer, intent(in) :: iDir
    !> The penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> Poly project
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    !> The state in nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    real(kind=rk), intent(in), optional :: nodal_GradData(:,:,:)
    real(kind=rk), intent(inout) :: nodal_res(:,:)
    !> Length of the element
    real(kind=rk), intent(in) :: ElemLength
    !> The scheme information
    type(atl_scheme_type), intent(inout) :: scheme_min
    type(atl_scheme_type), intent(inout) :: scheme_current
    ! -------------------------------------------------------------------- !
    ! Rotation indices for physical flux calculation in y and z direction
    integer :: rot(6)
    integer :: nDofs, nScalars
    ! -------------------------------------------------------------------- !

    ! get the rotation for the physical flux calculation in y and z direction.
    rot = equation%varRotation(iDir)%varTransformIndices(1:6)
    nDofs = size(state,1)
    nScalars = size(state,2)


    ! Variable material parameters work only for Q_space polynomials
    if(.not.(scheme_current%modg%basisType.eq.Q_space)) then
      call tem_abort( 'Error in modg_maxwell_physFlux: Variable material' &
        & // ' parameters work only for Q poly space'                     )
    end if

    call compute_physFlux_nonConst(                                      &
      &  nDofs         = nDofs,                                          &
      &  nScalars      = nScalars,                                       &
      &  state_der     = res(:ndofs,:),                                  &
      &  state         = state,                                          &
      &  rot           = rot,                                            &
      &  nelems        = material%material_desc                          &
      &                          %computeElems(2)%nElems ,               &
      &  material      = material%material_dat                           &
      &                          %elemMaterialData(2)%materialDat,       &
      &  poly_proj     = poly_proj,                                      &
      &  modalCoeffs   = scheme_min%temp_over,                           &
      &  nodalPhysFlux = scheme_min%temp_nodal,                          &
      &  iElem         = iElem                                           )

  end subroutine atl_modg_maxwell_physFlux_NonConst
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Compute the physical flux in x direction.
  !!
  !! For other directions a properly defined variable permutation can be used.
  !! This routine covers only constant material parameters.
  subroutine compute_physFlux( nDofs, nScalars, state_der, state, rot, &
    &                          inv_mu, inv_epsi                        )
    ! -------------------------------------------------------------------- !
    !> dimensions
    integer, intent(in) :: nDofs, nScalars
    real(kind=rk), intent(inout) :: state_der(nDofs,nScalars)
    !> State to compute the fluxes from.
    real(kind=rk), intent(in) :: state(nDofs,nScalars)
    !> Rotationing to index the variables.
    integer, intent(in) :: rot(6)
    real(kind=rk), intent(in) :: inv_mu, inv_epsi
    ! -------------------------------------------------------------------- !
    integer :: iDoF, iter
    ! -------------------------------------------------------------------- !

    do iter=lbound(state_der,2),ubound(state_der,2)
      state_der(:,iter) = 0.0_rk
    end do
    do iDoF=1,nDoFs
      state_der(iDoF,rot(1)) = 0.0_rk
      state_der(iDoF,rot(2)) = state(iDoF,rot(6)) * inv_mu
      state_der(iDoF,rot(3)) = - state(iDoF,rot(5)) * inv_mu
      state_der(iDoF,rot(4)) = 0.0_rk
      state_der(iDoF,rot(5)) = - state(iDoF,rot(3)) * inv_epsi
      state_der(iDoF,rot(6)) = state(iDoF,rot(2)) * inv_epsi
    end do

  end subroutine compute_physFlux
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Function for physical flux of the Maxwell equations in terms of D and B.
  function atl_physFluxMaxwell(state, material) result(flux)
    ! -------------------------------------------------------------------- !
    !> State to compute the fluxes from (D,B).
    real(kind=rk), intent(in) :: state(6)
    !> Material parameters (mu, epsilon) the flux calculation
    real(kind=rk), intent(in) :: material(2)
    !> The resulting flux in x direction
    real(kind=rk) :: flux(6)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: inv_mu, inv_epsi
    ! -------------------------------------------------------------------- !

    ! Get magnetic permeability
    inv_mu = 1.0_rk / material(1)
    ! Get electric permitivity
    inv_epsi = 1.0_rk / material(2)

    flux(1) = 0.0_rk
    flux(2) = state(6) * inv_mu
    flux(3) = state(5) * (-1.0_rk) * inv_mu
    flux(4) = 0.0_rk
    flux(5) = state(3) * (-1.0_rk) * inv_epsi
    flux(6) = state(2) * inv_epsi

  end function atl_physFluxMaxwell
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Compute the physical flux in x direction.
  !!
  !! For other directions a properly defined variable permutation can be used.
  !! This routine covers non-constant material parameters.
  subroutine compute_physFlux_nonConst( nDofs, nScalars, nElems, state_der, &
    &                                   state, rot, material, poly_proj,    &
    &                                   modalCoeffs, nodalPhysFlux, iElem   )
    ! -------------------------------------------------------------------- !
    !> dimensions
    integer, intent(in) :: nDofs, nScalars
    !> Array to store the fluxes in.
    !real(kind=rk), intent(inout) :: state_der(nTotalElems,nDofs,nScalars)
    real(kind=rk), intent(inout) :: state_der(nDofs,nScalars)
    !> State to compute the fluxes from.
    !real(kind=rk), intent(in) :: state(nTotalElems,nDofs,nScalars)
    real(kind=rk), intent(in) :: state(nDofs,nScalars)
    !> Rotationing to index the variables.
    !integer, intent(in) :: rotY(6), rotZ(6)
    integer, intent(in) :: rot(6)
    !> Number of elements.
    integer, intent(in) :: nElems
    !> Current element index
    integer, intent(in) :: iElem
    !> Material parameters (mu, epsilon) for all elements
    real(kind=rk), intent(in) :: material(nElems,nDofs,2)
    !> Data for projection method
    type(ply_poly_project_type) :: poly_proj
    !> Working array for modal coefficients of the current element in the loop.
    real(kind=rk), intent(inout) &
      & :: modalCoeffs(poly_proj%body_3D%oversamp_dofs,size(state,2),1)
    !> Working array for nodal representation of the physical flux along the
    !! 3 spatial directions.
    real(kind=rk), intent(inout) &
      & :: nodalPhysFlux(poly_proj%body_3D%nquadpoints,size(state,2),2)
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !

    nquadpoints = poly_proj%body_3D%nquadpoints

    ! --> modal space
    call ply_convert2oversample( state       = state(:,:),        &
      &                          poly_proj   = poly_proj,         &
      &                          nDim        = 3,                 &
      &                          modalCoeffs = modalCoeffs(:,:,1) )
    ! --> oversampling modal space

    ! Now, we transform the modal representation of this element to nodal
    ! space by making use of fast polynomial transformations (FPT)
    call ply_poly_project_m2n( me         = poly_proj,            &
      &                        dim        = 3 ,                   &
      &                        nVars      = nScalars,             &
      &                        nodal_data = nodalPhysFlux(:,:,2), &
      &                        modal_data = modalCoeffs(:,:,1)    )
    ! --> oversampling nodal space

    ! Calculate the physical flux point by point within this cell - x direction
    do iPoint = 1,nQuadPoints
      nodalPhysFlux(iPoint,rot,1) = atl_physFluxMaxwell( &
        & state    = nodalPhysFlux(iPoint,rot,2),        &
        & material = material(iElem,iPoint,:)            )
    end do

    ! Transform the nodal physical flux back to modal space
    ! --> oversampled nodal space
    call ply_poly_project_n2m( me         = poly_proj,            &
      &                        dim        = 3 ,                   &
      &                        nVars      = nScalars,             &
      &                        nodal_data = nodalPhysFlux(:,:,1), &
      &                        modal_data = modalCoeffs(:,:,1)    )
    ! --> oversampled modal space
    call ply_convertFromOversample( modalCoeffs = modalCoeffs(:,:,1), &
      &                             poly_proj   = poly_proj,          &
      &                             nDim        = 3,                  &
      &                             state       = state_der(:,:)      )
    ! --> modal space


  end subroutine compute_physFlux_nonConst
  ! ************************************************************************ !

end module atl_modg_maxwell_kernel_module
