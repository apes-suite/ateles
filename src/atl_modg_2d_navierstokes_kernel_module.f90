! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> author: Jens Zudrop
!! Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme for the compressible Navier-Stokes equation. This scheme is a spectral scheme for linear, convection dominated
!! partial differential equation systems.
module atl_modg_2d_navierstokes_kernel_module
  use env_module,               only: rk
  use treelmesh_module,         only: treelmesh_type
  use tem_faceData_module,      only: tem_faceIterator_type
  use treelmesh_module,         only: treelmesh_type

  use ply_poly_project_module,     only: ply_poly_project_type, assignment(=), &
                                       & ply_poly_project_m2n, &
                                       & ply_poly_project_n2m
  use ply_oversample_module,       only: ply_convertFromOversample

  use atl_equation_module,         only: atl_equations_type
  use atl_cube_elem_module,        only: atl_cube_elem_type
  use atl_scheme_module,           only: atl_scheme_type
  use atl_modg_2d_scheme_module,   only: atl_modg_2d_scheme_type
  use atl_facedata_module,         only: atl_facedata_type
  use atl_modg_2d_euler_kernel_module, only: atl_modg_2d_euler_oneDim_numFlux_const, &
                                           & atl_modg_2d_euler_oneDim_numFlux_nonconst
  use atl_physFluxEuler_2d_module, only: atl_physFluxEuler_2d
  use atl_physFluxNvrStk_2d_module, only: atl_viscPhysFluxNavierStokes_2d
  use atl_viscNumFlux_Nvrstk_2d_module, only: atl_viscNavierStokes_2d
  use atl_materialPrp_module,      only: atl_material_type
  use atl_penalization_module,     only: atl_penalizationData_type


  implicit none

  private

  public :: atl_modg_2d_navierstokes_numFlux,         &
    & atl_modg_2d_navierstokes_physFlux_const,        &
    & atl_modg_2d_navierstokes_physFlux_NonConst,     &
    & atl_modg_2d_navierstokes_penalization_const,    &
    & atl_modg_2d_navierstokes_penalization_NonConst, &
    & atl_get_penaltyIP_2d


contains


  !> Calculate the physical flux for the MODG scheme and
  !! Navier-Stokes equation (with constant penalizations).
  subroutine atl_modg_2d_navierstokes_physFlux_const( equation, res, state, &
    & iElem, iDir, penalizationData, poly_proj, material, nodal_data,       &
    & nodal_gradData, nodal_res, elemLength,  scheme_min, scheme_current    )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> To store the resulting phy flux in modal form
    real(kind=rk), intent(inout)     :: res(:,:)
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
    integer :: iPoint
    integer :: nquadpoints
    integer :: rot(4), derRot(2)
    ! -------------------------------------------------------------------- !
    nodal_res = 0.0_rk
    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_2D%nquadpoints

    ! get the rotation for the physical flux calculation in y direction
    rot = equation%varRotation(iDir)%varTransformIndices(1:4)
    ! get rotation for the derivatives
    derRot = equation%varRotation(iDir)%derTransformIndices(2:3) &
      & - equation%varRotation(iDir)%derTransformIndices(1)


    ! Calculate the physical flux point by point within this cell - x direction
    do iPoint = 1, nQuadPoints
      scheme_min%temp_nodal(iPoint,rot,1) = atl_physFluxEuler_2d(             &
        & state        = nodal_data(iPoint,rot),                              &
        & isenCoeff    = equation%euler%isen_coef,                            &
        & penalty_char = material%material_dat%elemMaterialData(1)            &
        &                                     %materialDat(iElem, 1, 1),      &
        & U_o          = material%material_dat%elemMaterialData(1)            &
        &                                     %materialDat(iElem, 1, iDir+1), &
        & porosity     = equation%euler%porosity                              )
    end do

    ! Calculate viscous physical flux point by point within this cell - x
    ! direction
    do iPoint = 1, nQuadPoints
      scheme_min%temp_nodal(iPoint,rot,2) = atl_viscPhysFluxNavierStokes_2d( &
        & state          = nodal_data(iPoint,rot),                           &
        & state_gradient = nodal_GradData(iPoint,rot, DerRot),               &
        & mu             = equation%NavierStokes%mu,                         &
        & lambda         = equation%NavierStokes%lambda,                     &
        & thermCond      = equation%NavierStokes%therm_cond,                 &
        & heatCap        = equation%euler%cv                                 )
    end do

    ! Add up the nodal data
    nodal_res(:,:) = scheme_min%temp_nodal(:,:,1) - scheme_min%temp_nodal(:,:,2)

  end subroutine atl_modg_2d_navierstokes_physFlux_const


  ! Calculate the penalization terms (for density, momentum, energy)
  ! The penalization terms are calculated in the sammer manner as
  ! the physical fluxes, i.e. in a nodal-modal approach
  subroutine atl_modg_2d_navierstokes_penalization_const( equation, poly_proj, &
    &            nodal_data, scheme_min, penalizationData, iElem, material     )
    ! -------------------------------------------------------------------- !
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> Poly project
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The state in nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    !> The scheme information
    type(atl_scheme_type), intent(inout) :: scheme_min
    !> The penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> The current Element
    integer, intent(in) :: iElem
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! -------------------------------------------------------------------- !
    integer :: nquadpoints, iPoint, glob_elem_ind
    real(kind=rk) :: temperature, pressure
    real(kind=rk) :: penalization(4)
    ! -------------------------------------------------------------------- !
    ! Calculate the penalization terms (for density, momentum, energy)
    ! The penalization terms are calculated in the samme way as the physical
    ! fluxes, i.e. in a nodal-modal approach
    ! Penalization for momentum and energy in pointwise manner
    glob_elem_ind = material%material_desc%computeElems(1)%totElemIndices(iElem)
    nquadpoints = poly_proj%body_2D%nquadpoints
    !> @todo PV 20150820 Get the correct penalization data here
    penalization = material%material_dat%elemMaterialData(1) &
      &                                 %materialDat(iElem, 1, :)

    do iPoint = 1,nQuadPoints
      pressure = (equation%euler%isen_coef - 1.0_rk) &
        & * ( nodal_data(iPoint,4) &
        &   - 0.5_rk * sum (nodal_data(iPoint,2:3)**2) &
        &   / nodal_data(iPoint,1) )
      temperature =  pressure / ( nodal_data(iPoint,1) * equation%euler%R )

      ! ... momentum
      scheme_min%temp_nodal(iPoint,1:2,1) = (-1.0_rk) * penalization(1) &
        & * ( nodal_data(iPoint,2:3) / nodal_data(iPoint,1)    &
        &   - penalization(2:3) )                    &
        & / equation%euler%viscous_permeability
      ! ... energy
      scheme_min%temp_nodal(iPoint,3,1) = (-1.0_rk) *penalization(1) &
        & * ( temperature - penalization(4) )      &
        & / equation%euler%thermal_permeability
    end do

    ! Transform penalizations back to modal space
    call ply_poly_project_n2m( me         = poly_proj,                    &
      &                        dim        = 2,                            &
      &                        nVars      = 3,                            &
      &                        nodal_data = scheme_min%temp_nodal(:,:,1), &
      &                        modal_data = scheme_min%temp_over(:,:,1)   )

    ! ... no penalization for density

    penalizationdata%penalization_data(glob_elem_ind,:,1) = 0.0_rk

    ! --> oversamp modal space for penalty terms (momentum + energy)
    call ply_convertFromOversample(                                           &
      & modalCoeffs = scheme_min%temp_over(:,:,1),                            &
      & poly_proj   = poly_proj,                                              &
      & nDim        = 2,                                                      &
      & state       = penalizationdata%penalization_data(glob_elem_ind,:,2:4) )

  end subroutine atl_modg_2d_navierstokes_penalization_const


  !> Calculate the physical flux for the MODG scheme and
  !! Navier-Stokes equation (with non-constant penalizations).
  subroutine atl_modg_2d_navierstokes_physFlux_NonConst( equation, res, state, &
    & iElem,iDir, penalizationData, poly_proj, material, nodal_data,           &
    & nodal_GradData, nodal_res, elemLength, scheme_min, scheme_current        )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> To store the resulting phy flux in modal form
    real(kind=rk), intent(inout)     :: res(:,:)
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
    integer :: iPoint
    integer :: nquadpoints
    integer :: rot(4), derRot(2)
    ! -------------------------------------------------------------------- !
    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_2D%nquadpoints


    ! get the rotation for the physical flux calculation in y direction
    rot = equation%varRotation(iDir)%varTransformIndices(1:4)
    ! get rotation for the derivatives
    derRot = equation%varRotation(iDir)%derTransformIndices(2:3) &
      & - equation%varRotation(iDir)%derTransformIndices(1)
    !ind =  minloc(abs(elems - iElem), 1)


    do iPoint = 1, nQuadPoints
      scheme_min%temp_nodal(iPoint,rot,1) = atl_physFluxEuler_2d(   &
        & state        = nodal_data(iPoint,rot),                    &
        & isenCoeff    = equation%euler%isen_coef,                  &
        & penalty_char = material%material_dat                      &
        &                        %elemMaterialData(2)               &
        &                        %materialDat(iElem,iPoint,1),      &
        & U_o          = material%material_dat                      &
        &                        %elemMaterialData(2)               &
        &                        %materialDat(iElem,iPoint,iDir+1), &
        & porosity     = equation%euler%porosity                    )
    end do

    do iPoint = 1, nQuadPoints
      scheme_min%temp_nodal(iPoint,rot,2) = atl_viscPhysFluxNavierStokes_2d( &
        & state          = nodal_data(iPoint,rot),                           &
        & state_gradient = nodal_GradData(iPoint,rot,derRot),                &
        & mu             = equation%NavierStokes%mu,                         &
        & lambda         = equation%NavierStokes%lambda,                     &
        & thermCond      = equation%NavierStokes%therm_cond,                 &
        & heatCap        = equation%euler%cv                                 )
    end do

    ! Add up the nodal data
    nodal_res(:,:) = scheme_min%temp_nodal(:,:,1) - scheme_min%temp_nodal(:,:,2)

  end subroutine atl_modg_2d_navierstokes_physFlux_nonconst



  ! Calculate the penalization terms (for density, momentum, energy)
  ! The penalization terms are calculated in the sammer manner as
  ! the physical fluxes, i.e. in a nodal-modal approach
  subroutine atl_modg_2d_navierstokes_penalization_Nonconst( equation,     &
    & poly_proj, nodal_data, scheme_min, penalizationData, iElem, material )
    ! --------------------------------------------------------------------------!
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> Poly project
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The state in nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    !> The scheme information
    type(atl_scheme_type), intent(inout) :: scheme_min
    !> The penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> The current Element
    integer, intent(in) :: iElem
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! --------------------------------------------------------------------------!
    integer :: nquadpoints, iPoint, glob_elem_ind
    real(kind=rk) :: temperature, pressure
    real(kind=rk) :: penalization(poly_proj%body_2D%nquadpoints,4)
    ! --------------------------------------------------------------------------!
    ! Calculate the penalization terms (for density, momentum, energy)
    ! The penalization terms are calculated in the samme way as the physical fluxes,
    ! i.e. in a nodal-modal approach
    ! Penalization for momentum and energy in pointwise manner
    nquadpoints   = poly_proj%body_2D%nquadpoints
    !> @todo PV 20150820 Get the correct penalization data here
    penalization  = material%material_dat%elemMaterialData(2)    &
      &                                  %materialDat(iElem, :, :)
    glob_elem_ind = material%material_desc%computeElems(2)%totElemIndices(iElem)


    do iPoint = 1,nQuadPoints
      pressure = (equation%euler%isen_coef - 1.0_rk) &
        & * ( nodal_data(iPoint,4) &
        &   - 0.5_rk * sum(nodal_data(iPoint,2:3)**2) / nodal_data(iPoint,1) )
      temperature =  pressure / ( nodal_data(iPoint,1) * equation%euler%R )

      ! ... momentum
      scheme_min%temp_nodal(iPoint,1:2,1) = (-1.0_rk)     &
        & *penalization(iPoint,1)                         &
        & * ( nodal_data(iPoint,2:3)/nodal_data(iPoint,1) &
        & - penalization(iPoint,2:3) )                    &
        & / equation%euler%viscous_permeability
      ! ... energy
      scheme_min%temp_nodal(iPoint,3,1) = (-1.0_rk)  &
        & * penalization(iPoint,1)                   &
        & * ( temperature - penalization(iPoint,4) ) &
        & / equation%euler%thermal_permeability
    end do

    ! Transform penalizations back to modal space
    call ply_poly_project_n2m( me         = poly_proj,                    &
     &                         dim        = 2,                            &
     &                         nVars      = 3,                            &
     &                         nodal_data = scheme_min%temp_nodal(:,:,1), &
     &                         modal_data = scheme_min%temp_over(:,:,1)   )

    ! ... no penalization for density
    penalizationdata%penalization_data(glob_elem_ind,:,1) = 0.0_rk

    ! --> oversamp modal space for penalty terms (momentum + energy)
    call ply_convertFromOversample(                                           &
      & modalCoeffs = scheme_min%temp_over(:,:,1),                            &
      & poly_proj   = poly_proj,                                              &
      & nDim        = 2,                                                      &
      & state       = penalizationdata%penalization_data(glob_elem_ind,:,2:4) )

  end subroutine atl_modg_2d_navierstokes_penalization_Nonconst


  !> Calculate the numerical flux for Navier-Stokes equation and MODG scheme
  subroutine atl_modg_2d_navierstokes_numFlux( mesh, equation, facedata,   &
    &                                          scheme, poly_proj, material )
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
      call atl_modg_2d_euler_oneDim_numFlux_const(                         &
        & equation       = equation,                                       &
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
        &                          %varTransformIndices(1:4),              &
        & material_left  = material%material_dat%faceMaterialData(iDir,1)  &
        &                                       %leftElemMaterialDat,      &
        & material_right = material%material_dat%faceMaterialData(iDir,1)  &
        &                                       %rightElemMaterialDat      )

      ! ... fluxes for non-constant penalization parameters
      call atl_modg_2d_euler_oneDim_numFlux_nonconst(                      &
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
        &                          %varTransformIndices(1:4),              &
        & material_left  = material%material_dat%faceMaterialData(iDir,2)  &
        &                                       %leftElemMaterialDat,      &
        & material_right = material%material_dat%faceMaterialData(iDir,2)  &
        &                                       %rightElemMaterialDat      )

      ! viscous numerical flux of the Navier-Stokes equation (sigma^*)
      call modg_2d_viscNavierStokes_oneDim_numFlux(        &
        & equation   = equation,                           &
        & facedata   = facedata,                           &
        & scheme     = scheme,                             &
        & faces      = mesh%faces%faces(iDir)%computeFace, &
        & faceDir    = iDir,                               &
        & poly_proj  = poly_proj,                          &
        & elemLength = mesh%length                         )

      ! stabilization viscous numerical flux of the Navier-Stokes equation (u^*)
      call modg_2d_stabViscNavierStokes_oneDim_numFlux(   &
        & equation  = equation,                           &
        & facedata  = facedata,                           &
        & scheme    = scheme,                             &
        & faces     = mesh%faces%faces(iDir)%computeFace, &
        & faceDir   = iDir,                               &
        & poly_proj = poly_proj                           )

    end do

  end subroutine atl_modg_2d_navierstokes_numFlux


  !> Numerical flux calculation for viscous part of the Navier-Stokes equation
  !! across the faces in a single spatial direction.
  subroutine modg_2d_viscNavierStokes_oneDim_numFlux( equation, facedata,     &
    &                                                 scheme, faces, faceDir, &
    &                                                 poly_proj, elemLength   )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation

    !> The face state if the equation
    type(atl_facedata_type), intent(inout) :: facedata

    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme

    !> The faces to calculate the fluxes for.
    type(tem_faceIterator_type), intent(in) :: faces

    !> The spatial direction of the faces you calc the fluxes for, use the
    !! following:
    !! * 1 --> x direction.
    !! * 2 --> y direction.
    !! * 3 --> z direction.
    integer, intent(in) :: faceDir

    !> Parameter for used projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> The length of an element
    real(kind=rk), intent(in) :: elemLength
    ! -------------------------------------------------------------------------!
    ! Loop vars
    integer :: iPoint

    ! Modal coefficients for elements left and right of the face:
    ! First dimension is the number of modal coefficients on the face, second
    ! is the number of variables.
    real(kind=rk), allocatable :: leftModalCoeff(:,:),        &
      &                           rightModalCoeff(:,:),       &
      &                           leftModalCoeff_gradX(:,:),  &
      &                           rightModalCoeff_gradX(:,:), &
      &                           leftModalCoeff_gradY(:,:),  &
      &                           rightModalCoeff_gradY(:,:)

    real(kind=rk) :: flux(4)

    ! Loop var for the faces.
    integer :: iside

    ! Element positions of the left and right element of the face.
    integer :: left_neighbor, right_neighbor

    ! The rotation indices for the flux calculation
    integer :: varRotation(4), gradRot(2)

    ! Nodal representation on the face (for left and right neighbor)
    real(kind=rk), allocatable :: pointValLeft(:,:), pointValRight(:,:)

    ! Nodal representations of gradients on the face (for left and right
    ! neighbor)
    real(kind=rk), allocatable :: pointValLeft_grad(:,:,:), &
      &                           pointValRight_grad(:,:,:)

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
    ! method. oversamp_dof is used for the oversampling loop
    nquadpoints = poly_proj%body_1D%nquadpoints
    oversamp_dofs = poly_proj%body_1D%oversamp_dofs
    nScalars = equation%varSys%nScalars
    nPVars = (scheme%maxPolyDegree + 1) * nScalars

    allocate( leftModalCoeff(oversamp_dofs, nScalars),        &
      &       rightModalCoeff(oversamp_dofs, nScalars),       &
      &       leftModalCoeff_gradX(oversamp_dofs, nScalars),  &
      &       rightModalCoeff_gradX(oversamp_dofs, nScalars), &
      &       leftModalCoeff_gradY(oversamp_dofs, nScalars),  &
      &       rightModalCoeff_gradY(oversamp_dofs, nScalars)  )

    ! Initialize oversampled modal coefficients with 0 to make sure that higher
    ! modes are properly filled.
    ! (Only possible if transformation is not 3D, as the modes would otherwise
    !  be used and modified by the transformation.)
    leftModalCoeff       = 0.0_rk
    leftModalCoeff_gradX = 0.0_rk
    leftModalCoeff_gradY = 0.0_rk

    rightModalCoeff       = 0.0_rk
    rightModalCoeff_gradX = 0.0_rk
    rightModalCoeff_gradY = 0.0_rk

    allocate( numFluxBuffer(oversamp_dofs, nScalars) )
    allocate( nodalNumFlux(nQuadPoints, nScalars) )

    allocate( pointValLeft(nQuadPoints, nScalars),         &
      &       pointValRight(nQuadPoints, nScalars),        &
      &       pointValLeft_grad(nQuadPoints, nScalars, 2), &
      &       pointValRight_grad(nQuadPoints, nScalars, 2) )

    ! The permutation indices we apply to enable numerical flux calculation
    ! across faces in x direction.
    varRotation = equation%varRotation(faceDir)%varTransformIndices(1:nScalars)
    gradRot = equation%varRotation(faceDir)%derTransformIndices(2:3) &
      &       - equation%varRotation(faceDir)%derTransformIndices(1)


    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, size(faces%leftPos)
      ! Get the fluid neighbors for this face.
      left_neighbor = faces%leftPos(iside)
      right_neighbor = faces%rightPos(iside)

      ! for the left element, we have to access the right face values
      ! and for the right, we have to acess the left face values.

      ! --> modal space
      do iVP = 1, nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)

        ! Modal coeffs of the state (left and right)
        leftModalCoeff(iPoint,iVar)                                       &
          & = facedata%faceRep(faceDir)%dat(left_neighbor, iPoint, iVar, 2)
        rightModalCoeff(iPoint,iVar)                                       &
          & = facedata%faceRep(faceDir)%dat(right_neighbor, iPoint, iVar, 1)

        ! Modal coeffs of the gradient in x direction (left and right)
        leftModalCoeff_gradX(iPoint, iVar)                        &
          & = facedata%faceRep(faceDir)                           &
          &           %dat(left_neighbor, iPoint, iVar+nScalars, 2)
        rightModalCoeff_gradX(iPoint, iVar)                        &
          & = facedata%faceRep(faceDir)                            &
          &           %dat(right_neighbor, iPoint, iVar+nScalars, 1)

        ! Modal coeffs of the gradient in y direction (left and right)
        leftModalCoeff_gradY(iPoint, iVar)                          &
          & = facedata%faceRep(faceDir)                             &
          &           %dat(left_neighbor, iPoint, iVar+2*nScalars, 2)
        rightModalCoeff_gradY(iPoint, iVar)                          &
          & = facedata%faceRep(faceDir)                              &
          &           %dat(right_neighbor, iPoint, iVar+2*nScalars, 1)

      end do
      ! --> oversamp modal space

      ! transform the 1D modal representation to nodal surface points
      ! State left
      call ply_poly_project_m2n( me         = poly_proj,     &
        &                        dim        = 1,             &
        &                        nVars      = nScalars,      &
        &                        nodal_data = pointValLeft,  &
        &                        modal_data = leftModalCoeff )
      ! Gradient of state - left
      call ply_poly_project_m2n( me         = poly_proj,                &
        &                        dim        = 1,                        &
        &                        nVars      = nScalars,                 &
        &                        nodal_data = pointValLeft_grad(:,:,1), &
        &                        modal_data = leftModalCoeff_gradX      )
      call ply_poly_project_m2n( me         = poly_proj,                &
        &                        dim        = 1,                        &
        &                        nVars      = nScalars,                 &
        &                        nodal_data = pointValLeft_grad(:,:,2), &
        &                        modal_data = leftModalCoeff_gradY      )

      ! State right
      call ply_poly_project_m2n( me         = poly_proj,      &
        &                        dim        = 1,              &
        &                        nVars      = nScalars,       &
        &                        nodal_data = pointValRight,  &
        &                        modal_data = rightModalCoeff )
      ! Gradient of state - right
      call ply_poly_project_m2n( me         = poly_proj,                 &
        &                        dim        = 1,                         &
        &                        nVars      = nScalars,                  &
        &                        nodal_data = pointValRight_grad(:,:,1), &
        &                        modal_data = rightModalCoeff_gradX      )
      call ply_poly_project_m2n( me         = poly_proj,                 &
        &                        dim        = 1 ,                        &
        &                        nVars      = nScalars,                  &
        &                        nodal_data = pointValRight_grad(:,:,2), &
        &                        modal_data = rightModalCoeff_gradY      )

      ! --> oversamp nodal space`

      ! for each of the surface points calculate the numerical flux
      do iPoint = 1, nQuadPoints

        call atl_viscNavierStokes_2d(                                    &
          & left       = pointValLeft(iPoint,varRotation),               &
          & left_grad  = pointValLeft_grad(iPoint,varRotation,gradRot),  &
          & right      = pointValRight(iPoint,varRotation),              &
          & right_grad = pointValRight_grad(iPoint,varRotation,gradRot), &
          & mu         = equation%NavierStokes%mu,                       &
          & lambda     = equation%NavierStokes%lambda,                   &
          & thermCond  = equation%NavierStokes%therm_cond,               &
          & heatCap    = equation%euler%cv,                              &
          & penaltyIP  = penaltyIP,                                      &
          & flux       = flux                                            )

        nodalNumFlux(iPoint,varRotation) = flux

      end do

      ! transform back to modal space (facial polynomial)
      call ply_poly_project_n2m( me         = poly_proj,    &
        &                        dim        = 1,            &
        &                        nVars      = nScalars,     &
        &                        nodal_data = nodalNumFlux, &
        &                        modal_data = numFluxBuffer )
      ! --> oversamp modal space

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      do iVP = 1,nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)

        facedata%faceFlux(faceDir)%dat(left_neighbor, iPoint, iVar, 2)       &
          & = facedata%faceFlux(faceDir)%dat(left_neighbor, iPoint, iVar, 2) &
          &   - numFluxBuffer(iPoint, iVar)
        facedata%faceFlux(faceDir)%dat(right_neighbor, iPoint, iVar, 1)       &
          & = facedata%faceFlux(faceDir)%dat(right_neighbor, iPoint, iVar, 1) &
          &   - numFluxBuffer(iPoint,iVar)
      end do
      ! --> modal space


    end do FaceLoop


  end subroutine modg_2d_viscNavierStokes_oneDim_numFlux


  !> Return the penalty parameter for the IP discretizations of higher order
  !! equations.
  function atl_get_penaltyIP_2d(maxPolyDegree, elemLength, ip_param) &
    & result(penaltyIP)
    ! -------------------------------------------------------------------------!
    !> The maximal polynomial degree of the discretization
    !! (starting from 0 for FVM)
    integer, intent(in) :: maxPolyDegree
    !> The length of an element.
    real(kind=rk), intent(in) :: elemLength
    !> The Interior Penalty paramter
    !! (should be large enough to ensure stability)
    real(kind=rk), intent(in) :: ip_param
    !> The resulting penalty parameter
    real(kind=rk) :: penaltyIP
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    penaltyIP = ip_param                                       &
      &         * real((maxPolyDegree+1)*(maxPolyDegree+2),rk) &
      &         / 2.0_rk                                       &
      &         / elemLength

  end function atl_get_penaltyIP_2d


  !> Numerical flux calculation for stab-viscous part of the Navier-Stokes
  !! equation across the faces in a single
  !! spatial direction.
  subroutine modg_2d_stabViscNavierStokes_oneDim_numFlux( equation, facedata, &
    &                                                     scheme, faces,      &
    &                                                     faceDir, poly_proj  )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face state if the equation
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    !> The faces to calculate the fluxes for.
    type(tem_faceIterator_type), intent(in) :: faces
    !> The spatial direction of the faces you calc the fluxes for, use the
    !! following:
    !! 1 --> x direction. \n
    !! 2 --> y direction. \n
    !! 3 --> z direction.
    integer, intent(in) :: faceDir
    !> Parameter for used projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! -------------------------------------------------------------------- !
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
    integer :: nquadpoints, ndofs, oversamp_dofs
    integer :: iVP, nPVars, nScalars
    ! --------------------------------------------------------------------------

     ! get correct amount of quadrature points and degree due to projection
     ! method. oversamp_dof and oversamp_degree is used for the oversampling
     ! loop
     nquadpoints = poly_proj%body_1D%nquadpoints
     ndofs = poly_proj%body_1D%ndofs
     oversamp_dofs = poly_proj%body_1D%oversamp_dofs

     allocate(numFluxBuffer(oversamp_dofs,4))

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

         numFluxBuffer(iPoint,iVar) =                                          &
           & ( facedata%faceRep(faceDir)%dat(left_neighbor,iPoint,iVar,2)      &
           &   + facedata%faceRep(faceDir)%dat(right_neighbor,iPoint,iVar,1) ) &
           & / 2.0_rk
       end do
       ! --> oversamp modal space

       ! Store the modal coefficients of the numerical flux. For the left
       ! element we have calculated the flux on the right face and vice versa.
       do iVP = 1,nPVars
         iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
         iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)
         facedata%faceFlux(faceDir)%dat(left_neighbor,iPoint,iVar+nScalars, 2) &
           & = numFluxBuffer(iPoint,iVar)
         facedata%faceFlux(faceDir)%dat(right_neighbor,   &
           &                            iPoint,           &
           &                            iVar+nScalars, 1) &
           & = numFluxBuffer(iPoint,iVar)
       end do
       ! --> modal space

     end do FaceLoop

  end subroutine modg_2d_stabViscNavierStokes_oneDim_numFlux

end module atl_modg_2d_navierstokes_kernel_module
