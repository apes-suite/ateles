! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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

!> Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme for the Reynolds Avereaged Navier-Stokes equation. This scheme is a spectral scheme for linear, convection dominated
!! partial differential equation systems.
module atl_modg_2d_filNvrStk_kernel_module
  use env_module,                   only: rk
  use tem_faceData_module,          only: tem_faceIterator_type
  use atl_equation_module,          only: atl_equations_type
  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_scheme_module,            only: atl_scheme_type
  use atl_modg_2d_scheme_module,    only: atl_modg_2d_scheme_type
  use atl_facedata_module,          only: atl_facedata_type
  use ply_poly_project_module,      only: ply_poly_project_type, &
    &                                     assignment(=),         &
    &                                     ply_poly_project_m2n,  &
    &                                     ply_poly_project_n2m
  use atl_materialPrp_module,       only: atl_material_type
  use ply_oversample_module,        only: ply_convert2oversample, &
    &                                     ply_convertFromOversample
  use atl_penalization_module,      only: atl_penalizationData_type
  use atl_physFluxFilNvrStk_module, only: atl_viscPhysFluxRans_2d, &
    &                                      atl_physFluxRans_2d
  use atl_modg_2d_navierstokes_kernel_module,                &
    &                               only: atl_get_penaltyIP_2d
  use atl_numFlux_filNvrStk_module, only: atl_viscRans_2d

  implicit none
  private

  public :: atl_modg_2d_filNvrStk_numFlux
  public :: atl_modg_2d_filNvrStk_physFlux_const
  public :: atl_modg_2d_filNvrStk_physFlux_NonConst


contains


  !> Calculate the physical flux for the MODG scheme and
  !! Navier-Stokes equation (with constant penalizations).
  subroutine atl_modg_2d_filNvrStk_physFlux_const( equation, res, state, &
    & iElem, iDir, penalizationData, poly_proj, material, nodal_data,    &
    & nodal_gradData, nodal_res, elemLength, scheme_min, scheme_current  )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------!
    integer :: iPoint
    integer :: rot(6), derRot(2)
    ! --------------------------------------------------------------------------!

    ! get the rotation for the physical flux calculation in y direction
    rot = equation%varRotation(iDir)%varTransformIndices(1:equation%varsys%nScalars)
    ! get rotation for the derivatives
    derRot = equation%varRotation(iDir)%derTransformIndices(2:3) &
      & - equation%varRotation(iDir)%derTransformIndices(1)


    ! Calculate the physical flux point by point within this cell - x direction
    do iPoint = 1, poly_proj%body_2D%nQuadPoints
      scheme_min%temp_nodal(iPoint,rot,1) = atl_physFluxRans_2d(       &
        & state        = nodal_data(iPoint,rot),                       &
        & isenCoeff    = equation%euler%isen_coef,                     &
        & penalty_char = material%material_dat%elemMaterialData(1)     &
        &                                     %materialDat(iElem,1,1), &
        & porosity     = equation%euler%porosity                       )
    end do

    ! Calculate viscous physical flux point by point within this cell - x direction
    do iPoint = 1, poly_proj%body_2D%nQuadPoints
      scheme_min%temp_nodal(iPoint,rot,2) = atl_viscPhysFluxRans_2d( &
        & state          = nodal_data(iPoint,rot),                   &
        & state_gradient = nodal_GradData(iPoint,rot, DerRot),       &
        & isenCoeff      = equation%euler%isen_coef,                 &
        & mu             = equation%NavierStokes%mu,                 &
        & lambda         = equation%NavierStokes%lambda,             &
        & thermCond      = equation%NavierStokes%therm_cond,         &
        & rans_params    = equation%FiltNavierStokes%rans,           &
        & heatCap        = equation%euler%cv                         )
    end do

    ! Add up the nodal data
    nodal_res(:,:) = scheme_min%temp_nodal(:,:,1) - scheme_min%temp_nodal(:,:,2)

  end subroutine atl_modg_2d_filNvrStk_physFlux_const


  !> Calculate the physical flux for the MODG scheme and
  !! Navier-Stokes equation (with non-constant penalizations).
  subroutine atl_modg_2d_filNvrStk_physFlux_NonConst(  equation, res, state, &
    & iElem,iDir, penalizationData, poly_proj, material, nodal_data,         &
    & nodal_GradData, nodal_res, elemLength, scheme_min, scheme_current      )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------!
    integer :: iPoint
    integer :: nquadpoints
    integer :: rot(6), derRot(2)
    real(kind=rk) :: penalization( poly_proj%body_2D%nquadpoints, 4 )
    ! --------------------------------------------------------------------------!
    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_2D%nquadpoints
    penalization = material%material_dat%elemMaterialData(2)    &
      &                                 %materialDat(iElem, :, :)


    ! get the rotation for the physical flux calculation in y direction
    rot = equation%varRotation(iDir)                             &
      &           %varTransformIndices(1:equation%varsys%nScalars)
    ! get rotation for the derivatives
    derRot = equation%varRotation(iDir)%derTransformIndices(2:3) &
      & - equation%varRotation(iDir)%derTransformIndices(1)


    do iPoint = 1, nQuadPoints
      scheme_min%temp_nodal(iPoint,rot,1) = atl_physFluxRans_2d( &
        & state        = nodal_data(iPoint,rot),                 &
        & isenCoeff    = equation%euler%isen_coef,               &
        & penalty_char = penalization(iPoint, 1),                &
        & porosity     = equation%euler%porosity                 )
    end do

    do iPoint = 1, nQuadPoints
      scheme_min%temp_nodal(iPoint,:,2) = atl_viscPhysFluxRans_2d( &
        & state          = nodal_data(iPoint,rot),                 &
        & state_gradient = nodal_GradData(iPoint,rot,derRot),      &
        & isenCoeff      = equation%euler%isen_coef,               &
        & mu             = equation%NavierStokes%mu,               &
        & lambda         = equation%NavierStokes%lambda,           &
        & thermCond      = equation%NavierStokes%therm_cond,       &
        & rans_params    = equation%FiltNavierStokes%rans,         &
        & heatCap        = equation%euler%cv                       )
    end do

    ! Add up the nodal data
    nodal_res(:,:) = scheme_min%temp_nodal(:,:,1) - scheme_min%temp_nodal(:,:,2)

  end subroutine atl_modg_2d_filNvrStk_physFlux_nonconst


  !> Calculate the numerical flux for Navier-Stokes equation and MODG scheme
  subroutine atl_modg_2d_filNvrStk_numFlux( mesh, equation, facedata, scheme, &
    &                                       poly_proj, material )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    !> Parameter for used projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! --------------------------------------------------------------------------
    integer :: iDir
    ! --------------------------------------------------------------------------

    ! Numerical flux for faces in all 2 spatial face directions (x and y dir)
    do iDir = 1,2
      ! convective part of the Euler equations (F^*)
      ! ... fluxes for constant penalization parameters
      call modg_2d_rans_oneDim_numFlux_const( equation = equation ,        &
        & nSides         = size(material%material_desc%computeFace(iDir,1) &
        &                                             %leftPos),           &
        & faceRep        = facedata%faceRep(iDir)%dat,                     &
        & faceFlux       = facedata%faceFlux(iDir)%dat,                    &
        & leftPos        = material%material_desc%computeFace(iDir,1)      &
        &                                        %leftPos,                 &
        & rightPos       = material%material_desc%computeFace(iDir,1)      &
        &                                        %rightPos,                &
        & poly_proj      = poly_proj,                                      &
        & varRotation    = equation%varRotation(iDir)                      &
        &                          %varTransformIndices(1:6),              &
        & material_left  = material%material_dat%faceMaterialData(iDir,1)  &
        &                                       %leftElemMaterialDat,      &
        & material_right = material%material_dat%faceMaterialData(iDir,1)  &
        &                                       %rightElemMaterialDat      )

      ! ... fluxes for non-constant penalization parameters
      call modg_2d_rans_oneDim_numFlux_nonconst(                           &
        & equation       = equation,                                       &
        & nSides         = size(material%material_desc%computeFace(iDir,2) &
        &                                             %leftPos),           &
        & faceRep        = facedata%faceRep(iDir)%dat,                     &
        & faceFlux       = facedata%faceFlux(iDir)%dat,                    &
        & leftPos        = material%material_desc%computeFace(iDir,2)      &
        &                                        %leftPos,                 &
        & rightPos       = material%material_desc%computeFace(iDir,2)      &
        &                                        %rightPos,                &
        & poly_proj      = poly_proj,                                      &
        & varRotation    = equation%varRotation(iDir)                      &
        &                          %varTransformIndices(1:6),              &
        & material_left  = material%material_dat%faceMaterialData(iDir,2)  &
        &                                       %leftElemMaterialDat,      &
        & material_right = material%material_dat%faceMaterialData(iDir,2)  &
        &                                       %rightElemMaterialDat      )

      ! viscous numerical flux of the Navier-Stokes equation (sigma^*)
      call modg_2d_viscRans_oneDim_numFlux(                &
        & equation   = equation,                           &
        & facedata   = facedata,                           &
        & scheme     = scheme ,                            &
        & faces      = mesh%faces%faces(iDir)%computeFace, &
        & faceDir    = iDir,                               &
        & poly_proj  = poly_proj,                          &
        & elemLength = mesh%length                         )

      ! stabilization viscous numerical flux of the Navier-Stokes equation (u^*)
      call modg_2d_stabViscRans_oneDim_numFlux(           &
        & equation  = equation,                           &
        & facedata  = facedata,                           &
        & scheme    = scheme,                             &
        & faces     = mesh%faces%faces(iDir)%computeFace, &
        & faceDir   = iDir,                               &
        & poly_proj = poly_proj                           )

    end do

  end subroutine atl_modg_2d_filNvrStk_numFlux


  !> Numerical flux calculation for Rans 2D equation across the faces in a
  !! single spatial direction (with constant penalization parameters).
  subroutine modg_2d_rans_oneDim_numFlux_const( equation, nSides, faceRep, &
    & faceFlux, leftPos, rightPos, poly_proj, varRotation, material_left,  &
    & material_right                                                       )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The number of faces to compute the flux for
    integer, intent(in) :: nSides
    !> The state on the face.
    real(kind=rk), intent(in) :: faceRep(:,:,:,:)
    !> The fluxes on the face.
    real(kind=rk), intent(inout) :: faceFlux(:,:,:,:)
    !> The positions of the faces to calculate the fluxes for (for elements
    !! left and right of the face).
    integer, intent(in) :: leftPos(:), rightPos(:)
    !> Parameter for used projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! The rotation indices for the flux calculation
    integer, intent(in) :: varRotation(6)
    !> The penalization material left and right of the face
    real(kind=rk), intent(in) :: material_left(:,:,:), material_right(:,:,:)
    ! --------------------------------------------------------------------------!
    ! Modal coefficients for elements left and right of the face:
    ! First dimension is the number of modal coefficients on the face, second
    ! is the number of variables.
    real(kind=rk), allocatable :: leftModalCoeffs(:,:), &
                               &  rightModalCoeffs(:,:)
    ! Loop var for the faces.
    integer :: iside
    ! Element positions of the left and right element of the face.
    integer :: left_neighbor, right_neighbor
    ! Nodal representation on the face (for left and right neighbor)
    real(kind=rk), allocatable :: pointValLeft(:,:), pointValRight(:,:)
    ! Nodal representation of the numerical flux
    real(kind=rk), allocatable :: nodalNumFlux(:,:)
    ! Modal representation of the flux on the face
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    ! Loop over variables (due to single variable FPTs)
    integer :: nquadpoints, oversamp_dofs
    ! --------------------------------------------------------------------------

    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_1D%nquadpoints
    oversamp_dofs = poly_proj%body_1D%oversamp_dofs

    allocate( leftModalCoeffs(oversamp_dofs, equation%varSys%nScalars))
    allocate( rightModalCoeffs(oversamp_dofs, equation%varSys%nScalars))

    allocate(numFluxBuffer(oversamp_dofs, equation%varSys%nScalars))
    allocate(nodalNumFlux(nQuadPoints, equation%varSys%nScalars))

    allocate( pointValLeft(nQuadPoints, equation%varSys%nScalars), &
      & pointValRight(nQuadPoints, equation%varSys%nScalars) )

    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, nSides
      ! Get the fluid neighbors for this face.
      left_neighbor = leftPos(iside)
      right_neighbor = rightPos(iside)

      ! for the left element, we have to access the right face values
      ! and for the right, we have to acess the left face values.
      ! --> modal space
      call ply_convert2oversample(state       = faceRep(left_neighbor,:,:,2), &
        &                         poly_proj   = poly_proj,                    &
        &                         nDim        = 1,                            &
        &                         modalCoeffs = leftModalCoeffs,              &
        &                         nScalars = equation%varSys%nScalars         )
      call ply_convert2oversample(state       = faceRep(right_neighbor,:,:,1), &
        &                         poly_proj   = poly_proj,                     &
        &                         nDim        = 1,                             &
        &                         modalCoeffs = rightModalCoeffs,              &
        &                         nScalars    = equation%varSys%nScalars       )
      ! --> oversamp modal space

      ! transform the 1D modal representation to nodal surface points
      call ply_poly_project_m2n(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = pointValLeft,             &
        &                       modal_data = leftModalCoeffs           )
      call ply_poly_project_m2n(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = pointValRight,            &
        &                       modal_data = rightModalCoeffs          )
      ! --> oversamp nodal space

      ! Use the numerical flux set by the equation system.
      ! Note, the rotation of the input state, the output is not rotated back
      ! here yet, as this is not allowed for output variables.
      call equation%Euler%numflux(                       &
        & state_left     = pointValLeft(:,varRotation),  &
        & state_right    = pointValRight(:,varRotation), &
        & material_left  = material_left(iSide,:,:),     &
        & material_right = material_right(iSide,:,:),    &
        & nPoints        = nQuadPoints,                  &
        & flux           = nodalNumFlux                  )

      ! transform back to modal space (facial polynomial)
      call ply_poly_project_n2m(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = nodalNumFlux,             &
        &                       modal_data = numFluxBuffer             )

      ! --> oversamp modal space
      call ply_convertFromOversample(                   &
        & modalCoeffs = numFluxBuffer,                  &
        & poly_proj   = poly_proj,                      &
        & nDim        = 1,                              &
        & state       = faceFlux(left_neighbor,:,:, 2), &
        & nScalars    = equation%varSys%nScalars        )

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      faceFlux(right_neighbor,:,varRotation,1)                               &
        &             = faceFlux(left_neighbor,:,:equation%varSys%nScalars,2 )
      faceFlux(left_neighbor,:,:equation%varSys%nScalars,2)                  &
        &             = faceFlux(right_neighbor,:,:equation%varSys%nScalars,1)

    end do FaceLoop

  end subroutine modg_2d_rans_oneDim_numFlux_const


  !> Numerical flux calculation for Rans equation across the faces in a single
  !! spatial direction (with non-constant penalization parameters).
  subroutine modg_2d_rans_oneDim_numFlux_nonconst( equation, nSides, faceRep, &
    & faceFlux, leftPos, rightPos, poly_proj, varRotation, material_left,     &
    & material_right                                                          )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The number of faces to compute the flux for
    integer, intent(in) :: nSides
    !> The state on the face.
    real(kind=rk), intent(in) :: faceRep(:,:,:,:)
    !> The fluxes on the face.
    real(kind=rk), intent(inout) :: faceFlux(:,:,:,:)
    !> The positions of the faces to calculate the fluxes for (for elements
    !! left and right of the face).
    integer, intent(in) :: leftPos(:), rightPos(:)
    !> Parameter for used projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! The rotation indices for the flux calculation
    integer, intent(in) :: varRotation(6)
    !> The penalization material left and right of the face
    real(kind=rk), intent(in) :: material_left(:,:,:), material_right(:,:,:)
    ! --------------------------------------------------------------------------!
    ! Modal coefficients for elements left and right of the face:
    ! First dimension is the number of modal coefficients on the face, second
    ! is the number of variables.
    real(kind=rk), allocatable :: leftModalCoeff(:,:), rightModalCoeff(:,:)
    ! Loop var for the faces.
    integer :: iside
    ! Element positions of the left and right element of the face.
    integer :: left_neighbor, right_neighbor
    ! Nodal representation on the face (for left and right neighbor)
    real(kind=rk), allocatable :: pointValLeft(:,:), pointValRight(:,:)
    ! Nodal representation of the numerical flux
    real(kind=rk), allocatable :: nodalNumFlux(:,:)
    ! Modal representation of the flux on the face
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    ! Loop over variables (due to single variable FPTs)
    integer :: nquadpoints, oversamp_dofs
    ! --------------------------------------------------------------------------

    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_1D%nquadpoints
    oversamp_dofs = poly_proj%body_1D%oversamp_dofs

    allocate( leftModalCoeff(oversamp_dofs, equation%varSys%nScalars))
    allocate( rightModalCoeff(oversamp_dofs, equation%varSys%nScalars))

    allocate(numFluxBuffer(oversamp_dofs, equation%varSys%nScalars))
    allocate(nodalNumFlux(nQuadPoints, equation%varSys%nScalars))

    allocate( pointValLeft(nQuadPoints, equation%varSys%nScalars), &
      & pointValRight(nQuadPoints, equation%varSys%nScalars) )

    ! The permutation indices we apply to enable numerical flux calculation across faces in x
    ! direction.


    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, nSides
      ! Get the fluid neighbors for this face.
      left_neighbor = leftPos(iside)
      right_neighbor = rightPos(iside)

      ! for the left element, we have to access the right face values
      ! and for the right, we have to acess the left face values.
      ! --> modal space
      call ply_convert2oversample(state       = faceRep(left_neighbor,:,:,2), &
        &                         poly_proj   = poly_proj,                    &
        &                         nDim        = 1,                            &
        &                         modalCoeffs = leftModalCoeff,               &
        &                         nScalars    = equation%varSys%nScalars      )
      call ply_convert2oversample(state       = faceRep(right_neighbor,:,:,1), &
        &                         poly_proj   = poly_proj,                     &
        &                         nDim        = 1,                             &
        &                         modalCoeffs = rightModalCoeff,               &
        &                         nScalars    = equation%varSys%nScalars       )
      ! --> oversamp modal space

      ! transform the 1D modal representation to nodal surface points
      call ply_poly_project_m2n(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = pointValLeft,             &
        &                       modal_data = leftModalCoeff            )
      call ply_poly_project_m2n(me         = poly_proj,                &
        &                       dim        = 1 ,                       &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = pointValRight,            &
        &                       modal_data = rightModalCoeff           )
      ! --> oversamp nodal space`

      ! Use the numerical flux set by the equation system.
      ! Note, the rotation of the input state, the output is not rotated back
      ! here yet, as this is not allowed for output variables.
      call equation%Euler%numflux(                       &
        & state_left     = pointValLeft(:,varRotation),  &
        & state_right    = pointValRight(:,varRotation), &
        & material_left  = material_left(iSide,:,:),     &
        & material_right = material_right(iSide,:,:),    &
        & nPoints        = nQuadPoints,                  &
        & flux           = nodalNumFlux                  )

      ! transform back to modal space (facial polynomial)
      call ply_poly_project_n2m( me         = poly_proj,                &
        &                        dim        = 1,                        &
        &                        nVars      = equation%varSys%nScalars, &
        &                        nodal_data = nodalNumFlux,             &
        &                        modal_data = numFluxBuffer             )

      ! --> oversamp modal space
      call ply_convertFromOversample(                     &
        & modalCoeffs = numFluxBuffer,                    &
        & poly_proj   = poly_proj,                        &
        & nDim        = 1,                                &
        & state       = faceFlux(left_neighbor, :, :, 2), &
        & nScalars    = equation%varSys%nScalars          )
      ! --> modal space

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      faceFlux(right_neighbor,:,varRotation,1)                               &
        &             = faceFlux(left_neighbor,:,:equation%varSys%nScalars,2 )
      faceFlux(left_neighbor,:,:equation%varSys%nScalars,2)                  &
        &             = faceFlux(right_neighbor,:,:equation%varSys%nScalars,1)

    end do FaceLoop

  end subroutine modg_2d_rans_oneDim_numFlux_nonconst


  !> Numerical flux calculation for viscous part of the RANS equation across the faces in a single
  !! spatial direction.
  subroutine modg_2d_viscRans_oneDim_numFlux( equation, facedata, scheme, &
                                            & faces, faceDir, poly_proj,  &
                                            & elemLength                  )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face state if the equation
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    !> The faces to calculate the fluxes for.
    type(tem_faceIterator_type), intent(in) :: faces
    !> The spatial direction of the faces you calc the fluxes for, use the following:
    !! 1 --> x direction. \n
    !! 2 --> y direction. \n
    !! 3 --> z direction.
    integer, intent(in) :: faceDir
    !> Parameter for used projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The length of an element
    real(kind=rk), intent(in) :: elemLength
    ! --------------------------------------------------------------------------!
    ! Loop vars
    integer :: iPoint
    ! Modal coefficients for elements left and right of the face:
    ! First dimension is the number of modal coefficients on the face, second
    ! is the number of variables.
    real(kind=rk), allocatable :: leftModalCoeff(:,:), &
     & rightModalCoeff(:,:),                           &
     & leftModalCoeff_gradX(:,:),                      &
     & rightModalCoeff_gradX(:,:),                     &
     & leftModalCoeff_gradY(:,:),                      &
     & rightModalCoeff_gradY(:,:)
    real(kind=rk) :: flux(6)
    ! Loop var for the faces.
    integer :: iside
    ! Element positions of the left and right element of the face.
    integer :: left_neighbor, right_neighbor
    ! The rotation indices for the flux calculation
    integer :: varRotation(6), gradRot(2)
    ! Nodal representation on the face (for left and right neighbor)
    real(kind=rk), allocatable :: pointValLeft(:,:), pointValRight(:,:)
    ! Nodal representations of gradients on the face (for left and right
    ! neighbor)
    real(kind=rk), allocatable :: pointValLeft_grad(:,:,:), &
      & pointValRight_grad(:,:,:)
    ! Nodal representation of the numerical flux
    real(kind=rk), allocatable :: nodalNumFlux(:,:)
    ! Modal representation of the flux on the face
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    ! Loop over variables (due to single variable FPTs)
    integer :: iVar
    integer :: nquadpoints, oversamp_dofs
    integer :: iVP, nPVars, nScalars
    real(kind=rk) :: penaltyIP
    ! --------------------------------------------------------------------------

    penaltyIP = atl_get_penaltyIP_2d(scheme%maxPolyDegree, elemLength, &
      & equation%navierstokes%ip_param)

    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_1D%nquadpoints
    oversamp_dofs = poly_proj%body_1D%oversamp_dofs

    allocate( leftModalCoeff(oversamp_dofs,  equation%varSys%nScalars) , &
      & rightModalCoeff(oversamp_dofs, equation%varSys%nScalars) ,       &
      & leftModalCoeff_gradX(oversamp_dofs,  equation%varSys%nScalars) , &
      & rightModalCoeff_gradX(oversamp_dofs, equation%varSys%nScalars) , &
      & leftModalCoeff_gradY(oversamp_dofs,  equation%varSys%nScalars) , &
      & rightModalCoeff_gradY(oversamp_dofs, equation%varSys%nScalars)   )

    allocate(numFluxBuffer(oversamp_dofs, equation%varSys%nScalars))
    allocate(nodalNumFlux(nQuadPoints, equation%varSys%nScalars))

    allocate( pointValLeft(nQuadPoints, equation%varSys%nScalars),  &
      & pointValRight(nQuadPoints, equation%varSys%nScalars),       &
      & pointValLeft_grad(nQuadPoints, equation%varSys%nScalars,2), &
      & pointValRight_grad(nQuadPoints, equation%varSys%nScalars,2) )

    ! The permutation indices we apply to enable numerical flux calculation
    ! across faces in x direction.
    varRotation(:) = equation%varRotation(faceDir)%varTransformIndices(1:6)
    gradRot = equation%varRotation(faceDir)%derTransformIndices(2:3) &
      & - equation%varRotation(faceDir)%derTransformIndices(1)


    nScalars = equation%varSys%nScalars
    nPVars = (scheme%maxPolyDegree+1)*equation%varSys%nScalars


    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, size(faces%leftPos)
      ! Get the fluid neighbors for this face.
      left_neighbor = faces%leftPos(iside)
      right_neighbor = faces%rightPos(iside)

      ! for the left element, we have to access the right face values
      ! and for the right, we have to acess the left face values.
      ! --> modal space
      leftModalCoeff(:,:) = 0.0_rk
      rightModalCoeff(:,:) = 0.0_rk
      leftModalCoeff_gradX(:,:) = 0.0_rk
      rightModalCoeff_gradX(:,:) = 0.0_rk
      leftModalCoeff_gradY(:,:) = 0.0_rk
      rightModalCoeff_gradY(:,:) = 0.0_rk
      do iVP = 1,nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)

        ! Modal coeffs of the state (left and right)
        leftModalCoeff(iPoint,iVar) =                                   &
          & facedata%faceRep(faceDir)%dat(left_neighbor, iPoint, iVar, 2)
        rightModalCoeff(iPoint,iVar) =                                   &
          & facedata%faceRep(faceDir)%dat(right_neighbor, iPoint, iVar, 1)

        ! Modal coeffs of the gradient in x direction (left and right)
        leftModalCoeff_gradX(iPoint,iVar) =                       &
          & facedata%faceRep(faceDir)                             &
          &         %dat(left_neighbor, iPoint, iVar + nScalars, 2)
        rightModalCoeff_gradX(iPoint,iVar) =                       &
          & facedata%faceRep(faceDir)                              &
          &         %dat(right_neighbor, iPoint, iVar + nScalars, 1)

        ! Modal coeffs of the gradient in y direction (left and right)
        leftModalCoeff_gradY(iPoint,iVar) =                           &
          & facedata%faceRep(faceDir)                                 &
          &         %dat(left_neighbor, iPoint, iVar + 2 * nScalars, 2)
        rightModalCoeff_gradY(iPoint,iVar) =                           &
          & facedata%faceRep(faceDir)                                  &
          &         %dat(right_neighbor, iPoint, iVar + 2 * nScalars, 1)

      end do
      ! --> oversamp modal space

      ! transform the 1D modal representation to nodal surface points
      ! State left
      call ply_poly_project_m2n(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = pointValLeft,             &
        &                       modal_data = leftModalCoeff            )
      ! Gradient of state - left
      call ply_poly_project_m2n(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = pointValLeft_grad(:,:,1), &
        &                       modal_data = leftModalCoeff_gradX      )
      call ply_poly_project_m2n(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = pointValLeft_grad(:,:,2), &
        &                       modal_data = leftModalCoeff_gradY      )

      ! State right
      call ply_poly_project_m2n(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = pointValRight,            &
        &                       modal_data = rightModalCoeff           )
      ! Gradient of state - right
      call ply_poly_project_m2n(me         = poly_proj,                 &
        &                       dim        = 1,                         &
        &                       nVars      = equation%varSys%nScalars,  &
        &                       nodal_data = pointValRight_grad(:,:,1), &
        &                       modal_data = rightModalCoeff_gradX      )
      call ply_poly_project_m2n(me         = poly_proj,                 &
        &                       dim        = 1 ,                        &
        &                       nVars      = equation%varSys%nScalars,  &
        &                       nodal_data = pointValRight_grad(:,:,2), &
        &                       modal_data = rightModalCoeff_gradY      )

      ! --> oversamp nodal space`

      ! for each of the surface points calculate the numerical flux
      do iPoint = 1, nQuadPoints

        call atl_viscRans_2d( left = pointValLeft(iPoint,varRotation),   &
          & left_grad   = pointValLeft_grad(iPoint,varRotation,gradRot), &
          & right       = pointValRight(iPoint,varRotation),             &
          & right_grad  = pointValRight_grad(iPoint,varRotation,gradRot),&
          & isen_coeff  = equation%euler%isen_coef,                      &
          & mu          = equation%NavierStokes%mu,                      &
          & lambda      = equation%NavierStokes%lambda,                  &
          & thermCond   = equation%NavierStokes%therm_cond,              &
          & heatCap     = equation%euler%cv,                             &
          & penaltyIP   = penaltyIP,                                     &
          & rans_params = equation%FiltNavierStokes%rans,                &
          & flux        = flux                                           )
        nodalNumFlux(iPoint,varRotation) = flux

      end do

      ! transform back to modal space (facial polynomial)
      call ply_poly_project_n2m(me         = poly_proj,                &
        &                       dim        = 1,                        &
        &                       nVars      = equation%varSys%nScalars, &
        &                       nodal_data = nodalNumFlux,             &
        &                       modal_data = numFluxBuffer             )
      ! --> oversamp modal space

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      do iVP = 1,nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)
        facedata%faceFlux(faceDir)%dat(left_neighbor,iPoint,iVar, 2)       &
          & = facedata%faceFlux(faceDir)%dat(left_neighbor,iPoint,iVar, 2) &
          &   - numFluxBuffer(iPoint,iVar)
        facedata%faceFlux(faceDir)%dat(right_neighbor,iPoint,iVar, 1)       &
          & = facedata%faceFlux(faceDir)%dat(right_neighbor,iPoint,iVar, 1) &
          &   - numFluxBuffer(iPoint,iVar)
      end do
      ! --> modal space

    end do FaceLoop

  end subroutine modg_2d_viscRans_oneDim_numFlux


  !> Numerical flux calculation for stab-viscous part of the RANS equation
  !!across the faces in a single spatial direction.
  subroutine modg_2d_stabViscRans_oneDim_numFlux( equation, facedata, scheme, &
                                                 & faces, faceDir, poly_proj  )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face state if the equation
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    !> The faces to calculate the fluxes for.
    type(tem_faceIterator_type), intent(in) :: faces
    !> The spatial direction of the faces you calc the fluxes for, use the following:
    !! 1 --> x direction. \n
    !! 2 --> y direction. \n
    !! 3 --> z direction.
    integer, intent(in) :: faceDir
    !> Parameter for used projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! --------------------------------------------------------------------------!
    ! Loop vars
    integer :: iPoint
    ! Loop var for the faces.
    integer :: iside
    ! Element positions of the left and right element of the face.
    integer :: left_neighbor, right_neighbor
    ! Modal representation of the flux on the face
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    ! Loop over variables (due to single variable FPTs)
    integer :: iVar
    integer :: oversamp_dofs
    integer :: iVP, nPVars, nScalars
    ! --------------------------------------------------------------------------

    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    oversamp_dofs = poly_proj%body_1D%oversamp_dofs

    allocate(numFluxBuffer(oversamp_dofs,6))

    nScalars = equation%varSys%nScalars
    nPVars = (scheme%maxPolyDegree+1)*equation%varSys%nScalars


    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, size(faces%leftPos)
      ! Get the fluid neighbors for this face.
      left_neighbor = faces%leftPos(iside)
      right_neighbor = faces%rightPos(iside)

      ! Calc the flux in modal way
      do iVP = 1,nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)

        numFluxBuffer(iPoint,iVar) &
          & = ( facedata%faceRep(faceDir) &
          &             %dat(left_neighbor, iPoint, iVar, 2)    &
          &     + facedata%faceRep(faceDir) &
          &               %dat(right_neighbor, iPoint, iVar, 1) ) &
          &   / 2.0_rk
      end do
      ! --> oversamp modal space

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      do iVP = 1,nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)
        facedata%faceFlux(faceDir)                               &
          &     %dat(left_neighbor, iPoint, iVar + nScalars, 2)  &
          & = numFluxBuffer(iPoint,iVar)
        facedata%faceFlux(faceDir)                               &
          &     %dat(right_neighbor, iPoint, iVar + nScalars, 1) &
          & = numFluxBuffer(iPoint,iVar)
      end do
      ! --> modal space

    end do FaceLoop

  end subroutine modg_2d_stabViscRans_oneDim_numFlux

end module atl_modg_2d_filNvrStk_kernel_module
