! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> author: Jens Zudrop
!! Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme for the Euler equation. This scheme is a spectral scheme for linear,
!! purley hyperbolic partial differential equation systems.
module atl_modg_euler_kernel_module
  use env_module,               only: rk

  use ply_poly_project_module,  only: ply_poly_project_type, &
    &                                 assignment(=),         &
    &                                 ply_poly_project_n2m,  &
    &                                 ply_poly_project_m2n

  use atl_equation_module,      only: atl_equations_type
  use atl_materialPrp_module,   only: atl_material_type, &
    &                                 atl_constMatIdx
  use atl_scheme_module,        only: atl_scheme_type
  use atl_physFluxEuler_module, only: atl_physFluxEuler, &
    &                                 atl_physFluxEuler_vec
  use atl_facedata_module,      only: atl_facedata_type
  use ply_oversample_module,    only: ply_convert2oversample,   &
   &                                  ply_convertFromoversample
  use atl_penalization_module,  only: atl_penalizationData_type

  implicit none

  private

  public :: atl_modg_euler_numflux,           &
    & atl_modg_euler_oneDim_numFlux_const,    &
    & atl_modg_euler_oneDim_numFlux_nonconst, &
    & atl_modg_euler_physFlux_const,          &
    & atl_modg_euler_physFlux_NonConst,       &
    & atl_modg_euler_penalization_NonConst,   &
    & atl_modg_euler_penalization_const

contains


  !> Calculate the physical flux for the MODG scheme and
  !! Euler equation.
  subroutine atl_modg_euler_physFlux_const ( equation, res, state, iElem,      &
    & iDir, penalizationData, poly_proj, material, nodal_data, nodal_gradData, &
    & nodal_res, elemLength,  scheme_min, scheme_current                       )
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
    integer :: rot(5)
    real(kind=rk), allocatable :: local_penalty(:)
    real(kind=rk), allocatable :: local_Uo(:)
    ! -------------------------------------------------------------------- !

    ! get the rotation for the physical flux calculation in y and z direction.
    rot = equation%varRotation(iDir)%varTransformIndices(1:5)

    allocate( local_penalty(poly_proj%body_3D%nquadpoints) )
    local_penalty = material%material_dat%elemMaterialData(1)    &
      &                                  %materialDat(iElem, 1, 1)
    allocate( local_Uo(poly_proj%body_3D%nquadpoints) )
    local_Uo = material%material_dat%elemMaterialData(1)           &
      &                             %materialDat(iElem, 1, iDir + 1)

    ! this routine needs an array of penalty, for nonconst case
    call atl_physFluxEuler_vec(                        &
       & state        = nodal_data(:,:),               &
       & isenCoeff    = equation%euler%isen_coef,      &
       & penalty_char = local_penalty,                 &
       & U_o          = local_Uo,                      &
       & porosity     = equation%euler%porosity,       &
       & nPoints      = poly_proj%body_3D%nquadpoints, &
       & rot          = rot,                           &
       & physFlux     = nodal_res(:,:)                 )
    deallocate( local_penalty )
    deallocate( local_Uo )

  end subroutine atl_modg_euler_physFlux_const


  ! Calculate the penalization terms (for density, momentum, energy)
  ! The penalization terms are calculated in the sammer manner as
  ! the physical fluxes, i.e. in a nodal-modal approach
  subroutine atl_modg_euler_penalization_const(equation, poly_proj, &
    & nodal_data, scheme_min, penalizationData, iElem, material     )
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
    integer :: nquadpoints, iPoint
    real(kind=rk) :: temperature, pressure, inv_dens
    integer ::  glob_elem_ind
    real(kind=rk) :: penalization(5)
    ! -------------------------------------------------------------------- !
    nquadpoints = poly_proj%body_3D%nquadpoints
    glob_elem_ind = material%material_desc%computeElems(1)%totElemIndices(iElem)
    !> @todo PV 20150820 Get the correct penalization data here
    penalization = material%material_dat                      &
      &                    %elemMaterialData(atl_ConstMatIdx) &
      &                    %materialDat(iElem,1,:)

    ! Penalization for momentum and energy in pointwise manner
    do iPoint = 1,nQuadPoints

      inv_dens = 1.0 / nodal_data(iPoint,1)
      pressure = (equation%euler%isen_coef - 1.0_rk)             &
        & * (nodal_data (iPoint,5)                               &
        &   - 0.5_rk * sum(nodal_data(iPoint,2:4)**2) * inv_dens )
      temperature =  pressure*inv_dens / (equation%euler%R )

      ! ... momentum
      ! Here just rewriting the nodal_data to store the nodal penalty values
      scheme_min%temp_nodal(iPoint,1:3,1) = (-1.0_rk) * penalization(1) &
        & * (nodal_data(iPoint,2:4) / nodal_data(iPoint,1)              &
        &   - penalization(2:4) )                                       &
        & / equation%euler%viscous_permeability
      ! ... energy
      scheme_min%temp_nodal(iPoint,4,1) = (-1.0_rk) * penalization(1) &
        & * ( temperature - penalization(5) )                         &
        & / equation%euler%thermal_permeability
    end do

    ! Transform penalizations back to modal space
    call ply_poly_project_n2m(me         = poly_proj,                    &
      &                       dim        = 3,                            &
      &                       nVars      = 4,                            &
      &                       nodal_data = scheme_min%temp_nodal(:,:,1), &
      &                       modal_data = scheme_min%temp_over(:,:,1)  )


    call ply_convertFromOversample(                                            &
      & state       = penalizationdata%penalization_data(glob_elem_ind,:,2:5), &
      & poly_proj   = poly_proj,                                               &
      & nDim        = 3,                                                       &
      & modalCoeffs = scheme_min%temp_over(:,:4,1)                             )

    ! ... no penalization for density
    penalizationdata%penalization_data(glob_elem_ind,:,1) = 0.0_rk

  end subroutine atl_modg_euler_penalization_const


  !> Calculate the physical flux for the MODG scheme and
  !! Euler equation.
  subroutine atl_modg_euler_physFlux_NonConst( equation, res, state, iElem,    &
    & iDir, penalizationData, poly_proj, material, nodal_data, nodal_gradData, &
    & nodal_res, elemLength,  scheme_min, scheme_current                       )
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
    integer :: rot(5)
    ! -------------------------------------------------------------------- !

    ! get the rotation for the physical flux calculation in y and z direction.
    rot = equation%varRotation(iDir)%varTransformIndices(1:5)

    call atl_physFluxEuler_vec(                                                &
       & state        = nodal_data(:,:),                                       &
       & isenCoeff    = equation%euler%isen_coef,                              &
       & penalty_char = material%material_dat                                  &
       &                        %elemMaterialData(2)                           &
       &                        %materialDat( iElem,                           &
       &                                      1:poly_proj%body_3D%nquadpoints, &
       &                                      1 ),                             &
       & U_o          = material%material_dat                                  &
       &                        %elemMaterialData(2)                           &
       &                        %materialDat( iElem,                           &
       &                                      1:poly_proj%body_3D%nquadpoints, &
       &                                      iDir + 1 ),                      &
       & porosity     = equation%euler%porosity,                               &
       & nPoints      = poly_proj%body_3D%nquadpoints,                         &
       & rot          = rot,                                                   &
       & physFlux     = nodal_res(:,:)                                         )

  end subroutine atl_modg_euler_physFlux_NonConst

  ! Calculate the penalization terms (for density, momentum, energy)
  ! The penalization terms are calculated in the sammer manner as
  ! the physical fluxes, i.e. in a nodal-modal approach
  subroutine atl_modg_euler_penalization_NonConst(equation, poly_proj, &
    & nodal_data, scheme_min, penalizationData, iElem, material        )
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
    integer :: iPoint
    real(kind=rk) :: temperature, pressure, inv_dens
    real(kind=rk) :: penalization(5)
    integer :: glob_elem_ind
    ! -------------------------------------------------------------------- !
    glob_elem_ind = material%material_desc%computeElems(2)%totElemIndices(iElem)

    ! Penalization for momentum and energy in pointwise manner
    do iPoint = 1, poly_proj%body_3D%nquadpoints
      !> @todo PV 20150820 Get the correct penalization data here
      penalization = material%material_dat%elemMaterialData(2) &
        &                    %materialDat(iElem,iPoint,:)

      inv_dens = 1.0 / nodal_data(iPoint,1)
      pressure = (equation%euler%isen_coef - 1.0_rk )            &
        & * ( nodal_data(iPoint,5)                               &
        &   - 0.5_rk * sum(nodal_data(iPoint,2:4)**2) * inv_dens )
      temperature =  pressure * inv_dens / (equation%euler%R)

      ! ... momentum
      ! Here just rewriting the nodal_data to store the nodal penalty values
      scheme_min%temp_nodal( iPoint, 1:3, 1 ) =                   &
        & (-1.0_rk) * penalization(1)                             &
        & * ( nodal_data( iPoint, 2:4 ) / nodal_data( iPoint, 1 ) &
        &   - penalization( 2:4 ) )                               &
        & / equation%euler%viscous_permeability
      ! ... energy
      scheme_min%temp_nodal( iPoint, 4, 1 ) = &
        & (-1.0_rk) * penalization(1)         &
        & * ( temperature - penalization(5) ) &
        & / equation%euler%thermal_permeability
    end do

    ! Transform penalizations back to modal space
    call ply_poly_project_n2m( me         = poly_proj,                    &
      &                        dim        = 3 ,                           &
      &                        nVars      = 4,                            &
      &                        nodal_data = scheme_min%temp_nodal(:,:,1), &
      &                        modal_data = scheme_min%temp_over(:,:,1)   )

    call ply_convertFromOversample(                                           &
      & state       = penalizationdata%penalization_data(glob_elem_ind,:,2:), &
      & poly_proj   = poly_proj,                                              &
      & nDim        = 3,                                                      &
      & modalCoeffs = scheme_min%temp_over(:,:4,1)                            )

    ! ... no penalization for density
    penalizationdata%penalization_data(glob_elem_ind,:,1) = 0.0_rk

  end subroutine atl_modg_euler_penalization_NonConst


  !> Calculate the numerical flux for Euler equation and MODG scheme
  subroutine atl_modg_euler_numFlux( equation, facedata, poly_proj, material )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameter for projection.
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! --------------------------------------------------------------------------
    integer :: iDir
    ! --------------------------------------------------------------------------

    ! Numerical flux for faces in all 3 spatial face directions
    do iDir = 1,3
      call atl_modg_euler_oneDim_numFlux_const(               &
        & equation       = equation ,                         &
        & nSides         = size(material%material_desc        &
        &                               %computeFace(iDir,1)  &
        &                               %leftPos),            &
        & faceRep        = facedata%faceRep(iDir)%dat,        &
        & faceFlux       = facedata%faceFlux(iDir)%dat,       &
        & leftPos        = material%material_desc             &
        &                          %computeFace(iDir,1)       &
        &                          %leftPos,                  &
        & rightPos       = material%material_desc             &
        &                          %computeFace(iDir,1)       &
        &                          %rightPos,                 &
        & poly_proj      = poly_proj,                         &
        & varRotation    = equation%varRotation(iDir)         &
        &                          %varTransformIndices(1:5), &
        & material_left  = material%material_dat              &
        &                          %faceMaterialData(iDir,1)  &
        &                          %leftElemMaterialDat,      &
        & material_right = material%material_dat              &
        &                          %faceMaterialData(iDir,1)  &
        &                          %rightElemMaterialDat      )


      call atl_modg_euler_oneDim_numFlux_nonconst(            &
        & equation       = equation,                          &
        & nSides         = size(material%material_desc        &
        &                               %computeFace(iDir,2)  &
        &                               %leftPos),            &
        & faceRep        = facedata%faceRep(iDir)%dat,        &
        & faceFlux       = facedata%faceFlux(iDir)%dat,       &
        & leftPos        = material%material_desc             &
        &                          %computeFace(iDir,2)       &
        &                          %leftPos,                  &
        & rightPos       = material%material_desc             &
        &                          %computeFace(iDir,2)       &
        &                          %rightPos,                 &
        & poly_proj      = poly_proj,                         &
        & varRotation    = equation%varRotation(iDir)         &
        &                          %varTransformIndices(1:5), &
        & material_left  = material%material_dat              &
        &                          %faceMaterialData(iDir,2)  &
        &                          %leftElemMaterialDat,      &
        & material_right = material%material_dat              &
        &                          %faceMaterialData(iDir,2)  &
        &                          %rightElemMaterialDat      )
    end do

  end subroutine atl_modg_euler_numFlux

  !> Numerical flux calculation for Euler equation across the faces in a single
  !! spatial direction.
  subroutine atl_modg_euler_oneDim_numFlux_const( equation, nSides, faceRep, &
    & faceFlux,  leftPos, rightPos, poly_proj, varRotation, material_left,   &
    & material_right                                                         )
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
    integer, intent(in) :: varRotation(equation%varSys%nScalars)
    !> The penalization material left and right of the face
    real(kind=rk), intent(in) :: material_left(:,:,:), material_right(:,:,:)
    ! -------------------------------------------------------------------- !
    ! Modal coefficients for elements left and right of the face:
    ! First dimension is the number of modal coefficients on the face, second
    ! is the number of variables.
    real(kind=rk), allocatable :: leftModalCoeff(:,:), &
                               &  rightModalCoeff(:,:)
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
    nquadpoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    allocate( leftModalCoeff(oversamp_dofs, equation%varSys%nScalars))
    allocate( rightModalCoeff(oversamp_dofs, equation%varSys%nScalars))
    allocate( numFluxBuffer(oversamp_dofs, equation%varSys%nScalars))
    allocate( nodalNumFlux(nQuadPoints, equation%varSys%nScalars))
    allocate( pointValLeft(nQuadPoints, equation%varSys%nScalars), &
      &       pointValRight(nQuadPoints, equation%varSys%nScalars) )

    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, nSides
      ! Get the fluid neighbors for this face.
      left_neighbor = leftPos(iside)
      right_neighbor = rightPos(iside)

      ! for the left element, we have to access the right face values
      ! and for the right, we have to acess the left face values.
      ! --> modal space
      if (equation%euler%ensure_positivity) then
        call ply_convert2oversample(                                        &
          & state             = faceRep(left_neighbor,:,:,2),               &
          & poly_proj         = poly_proj,                                  &
          & nDim              = 2,                                          &
          & modalCoeffs       = leftModalCoeff,                             &
          & nScalars          = equation%varSys%nScalars,                   &
          & ensure_positivity = [.true., .false., .false., .false., .true.] )
        call ply_convert2oversample(                                        &
          & state             = faceRep(right_neighbor,:,:,1),              &
          & poly_proj         = poly_proj,                                  &
          & nDim              = 2,                                          &
          & modalCoeffs       = rightModalCoeff,                            &
          & nScalars          = equation%varSys%nScalars,                   &
          & ensure_positivity = [.true., .false., .false., .false., .true.] )
      else
        call ply_convert2oversample(                    &
          & state       = faceRep(left_neighbor,:,:,2), &
          & poly_proj   = poly_proj,                    &
          & nDim        = 2,                            &
          & modalCoeffs = leftModalCoeff,               &
          & nScalars    = equation%varSys%nScalars      )
        call ply_convert2oversample(                     &
          & state       = faceRep(right_neighbor,:,:,1), &
          & poly_proj   = poly_proj,                     &
          & nDim        = 2,                             &
          & modalCoeffs = rightModalCoeff,               &
          & nScalars    = equation%varSys%nScalars       )
      end if
      ! --> oversampled modal space

      ! transform the 2D modal representation to nodal surface points
      call ply_poly_project_m2n( me         = poly_proj,                &
        &                        dim        = 2,                        &
        &                        nVars      = equation%varSys%nScalars, &
        &                        nodal_data = pointValLeft,             &
        &                        modal_data = leftModalCoeff            )
      call ply_poly_project_m2n( me         = poly_proj,                &
        &                        dim        = 2 ,                       &
        &                        nVars      = equation%varSys%nScalars, &
        &                        nodal_data = pointValRight,            &
        &                        modal_data = rightModalCoeff           )
      ! --> oversampling nodal space
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
        &                        dim        = 2,                        &
        &                        nVars      = equation%varSys%nScalars, &
        &                        nodal_data = nodalNumFlux,             &
        &                        modal_data = numFluxBuffer             )

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      !oversampling loop
      call ply_convertFromOversample(                                           &
        & modalCoeffs = numFluxBuffer,                                          &
        & poly_proj   = poly_proj,                                              &
        & nDim        = 2,                                                      &
        & state       = faceFlux(left_neighbor,:,1:equation%varSys%nScalars,2), &
        & nScalars    = equation%varSys%nScalars                                )

      ! Set the left flux to the right flux and rotate the state back according
      ! to the varRotation.
      ! We are basically using the one flux copy as a buffer, and replace it
      ! finally by the rotated state. Doing the rotation here is potentially
      ! cheaper, as only the not oversampled state needs to be copied.
      ! In case equation%nDerivatives = 1 the size of faceFlux(:,:,2*nScalars,2)
      ! so the first nScalar values need to be swapped
      faceFlux(right_neighbor,:,varRotation,1)                   &
        & = faceFlux(left_neighbor,:,1:equation%varSys%nScalars,2)
      faceFlux(left_neighbor,:,1:equation%varSys%nScalars,2)      &
        & = faceFlux(right_neighbor,:,1:equation%varSys%nScalars,1)

    end do FaceLoop

  end subroutine atl_modg_euler_oneDim_numFlux_const


  !> Numerical flux calculation for Euler equation across the faces in a single
  !! spatial direction.
  subroutine atl_modg_euler_oneDim_numFlux_nonconst( equation, nSides, &
    & faceRep, faceFlux, leftPos, rightPos, poly_proj, varRotation,    &
    & material_left, material_right                                    )
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
    integer, intent(in) :: varRotation(equation%varSys%nScalars)
    !> The penalization material left and right of the face
    real(kind=rk), intent(in) :: material_left(:,:,:), material_right(:,:,:)
    ! -------------------------------------------------------------------- !
    ! Modal coefficients for elements left and right of the face:
    ! First dimension is the number of modal coefficients on the face, second
    ! is the number of variables.
    real(kind=rk), allocatable :: leftModalCoeff(:,:), &
                               &  rightModalCoeff(:,:)
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
    nquadpoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    allocate( leftModalCoeff(oversamp_dofs, equation%varSys%nScalars))
    allocate( rightModalCoeff(oversamp_dofs, equation%varSys%nScalars))
    allocate( numFluxBuffer(oversamp_dofs,equation%varSys%nScalars))
    allocate( nodalNumFlux(nQuadPoints, equation%varSys%nScalars))
    allocate( pointValLeft(nQuadPoints, equation%varSys%nScalars), &
      &       pointValRight(nQuadPoints, equation%varSys%nScalars) )


    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, nSides
      ! Get the fluid neighbors for this face.
      left_neighbor = leftPos(iside)
      right_neighbor = rightPos(iside)

      ! for the left element, we have to access the right face values
      ! and for the right, we have to acess the left face values.
      ! --> modal space
      if (equation%euler%ensure_positivity) then
        call ply_convert2oversample(                                        &
          & state             = faceRep(left_neighbor,:,:,2),               &
          & poly_proj         = poly_proj,                                  &
          & nDim              = 2,                                          &
          & modalCoeffs       = leftModalCoeff,                             &
          & nScalars          = equation%varSys%nScalars,                   &
          & ensure_positivity = [.true., .false., .false., .false., .true.] )
        call ply_convert2oversample(                                        &
          & state             = faceRep(right_neighbor,:,:,1),              &
          & poly_proj         = poly_proj,                                  &
          & nDim              = 2,                                          &
          & modalCoeffs       = rightModalCoeff,                            &
          & nScalars          = equation%varSys%nScalars,                   &
          & ensure_positivity = [.true., .false., .false., .false., .true.] )
      else
        call ply_convert2oversample(                    &
          & state       = faceRep(left_neighbor,:,:,2), &
          & poly_proj   = poly_proj,                    &
          & nDim        = 2,                            &
          & modalCoeffs = leftModalCoeff,               &
          & nScalars    = equation%varSys%nScalars      )
        call ply_convert2oversample(                     &
          & state       = faceRep(right_neighbor,:,:,1), &
          & poly_proj   = poly_proj,                     &
          & nDim        = 2,                             &
          & modalCoeffs = rightModalCoeff,               &
          & nScalars    = equation%varSys%nScalars       )
      end if
      ! --> oversampled modal space

      ! transform the 2D modal representation to nodal surface points
      call ply_poly_project_m2n( me         = poly_proj,                &
        &                        dim        = 2,                        &
        &                        nVars      = equation%varSys%nScalars, &
        &                        nodal_data = pointValLeft,             &
        &                        modal_data = leftModalCoeff            )
      call ply_poly_project_m2n( me         = poly_proj,                &
        &                        dim        = 2,                        &
        &                        nVars      = equation%varSys%nScalars, &
        &                        nodal_data = pointValRight,            &
        &                        modal_data = rightModalCoeff           )
      ! --> oversampling nodal space

      ! Use the numerical flux set by the equation system.
      call equation%Euler%numflux(                       &
        & state_left     = pointValLeft(:,varRotation),  &
        & state_right    = pointValRight(:,varRotation), &
        & material_left  = material_left(iSide,:,:),     &
        & material_right = material_right(iSide,:,:),    &
        & nPoints        = nQuadPoints,                  &
        & flux           = nodalNumFlux                  )

      ! transform back to modal space (facial polynomial)
      call ply_poly_project_n2m( me         = poly_proj,                &
        &                        dim        = 2,                        &
        &                        nVars      = equation%varSys%nScalars, &
        &                        nodal_data = nodalNumFlux,             &
        &                        modal_data = numFluxBuffer             )

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      !oversampling loop
      call ply_convertFromOversample(                                           &
        & modalCoeffs = numFluxBuffer,                                          &
        & poly_proj   = poly_proj,                                              &
        & nDim        = 2,                                                      &
        & state       = faceFlux(left_neighbor,:,1:equation%varSys%nScalars,2), &
        & nScalars    = equation%varSys%nScalars                                )

      ! Set the left flux to the right flux and rotate the state back according
      ! to the varRotation.
      ! We are basically using the one flux copy as a buffer, and replace it
      ! finally by the rotated state. Doing the rotation here is potentially
      ! cheaper, as only the not oversampled state needs to be copied.

      faceFlux(right_neighbor,:,varRotation,1) =               &
        & faceFlux(left_neighbor,:,1:equation%varSys%nScalars,2)
      faceFlux(left_neighbor,:,1:equation%varSys%nScalars,2) =  &
        & faceFlux(right_neighbor,:,1:equation%varSys%nScalars,1)

    end do FaceLoop

  end subroutine atl_modg_euler_oneDim_numFlux_nonconst

end module atl_modg_euler_kernel_module
