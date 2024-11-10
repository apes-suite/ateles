! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
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
!! scheme for the compressible filtered Navier-Stokes equation. This scheme is
!! a spectral scheme for linear, convection dominated partial differential
!! equation systems.
module atl_modg_filNvrStk_kernel_module
  use env_module,                  only: rk
  use tem_aux_module,              only: tem_abort
  use tem_logging_module,          only: logUnit
  use atl_equation_module,         only: atl_equations_type
  use atl_cube_elem_module,        only: atl_cube_elem_type
  use atl_scheme_module,           only: atl_scheme_type
  use atl_modg_scheme_module,      only: atl_modg_scheme_type
  use atl_facedata_module,         only: atl_facedata_type
  use ply_poly_project_module,     only: ply_poly_project_type, assignment(=)
  use atl_modg_euler_kernel_module, only: atl_modg_euler_oneDim_numFlux_const, &
                                        & atl_modg_euler_oneDim_numFlux_nonconst
  use atl_physFluxFilNvrStk_module, only: atl_viscPhysFluxRans, &
                                          atl_physFluxRans
  use atl_materialPrp_module,      only: atl_material_type
  use atl_penalization_module,     only: atl_penalizationData_type
  use atl_modg_navierStokes_kernel_module,              &
    & only: atl_modg_viscNavierStokes_oneDim_numFlux,   &
    &       atl_modg_stabViscNavierStokes_oneDim_numFlux


  implicit none
  private

  public :: atl_modg_filNvrStk_numFlux,  &
    & atl_modg_filNvrStk_physFlux_const, &
    & atl_modg_filNvrStk_physFlux_NonConst

contains

  !> Calculate the physical flux for the MODG scheme and
  !! Navier-Stokes equation (with constant penalizations).
  subroutine atl_modg_filNvrStk_physFlux_const (equation, res, state, iElem,  &
    & iDir,penalizationData, poly_proj, material, nodal_data, nodal_gradData, &
    & nodal_res, elemLength, scheme_min, scheme_current                       )
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
    integer :: rot(7), derRot(3)
    ! -------------------------------------------------------------------- !

    ! get the rotation for the physical flux calculation in y direction
    rot = equation%varRotation(iDir)                              &
      &           %varTransformIndices(1:equation%varsys%nScalars )

    ! get rotation for the derivatives
    derRot = equation%varRotation(iDir)%derTransformIndices(2:4) &
      & - equation%varRotation(iDir)%derTransformIndices(1)


    ! Calculate the physical flux point by point within this cell
    do iPoint = 1, poly_proj%body_3D%nquadpoints
      scheme_min%temp_nodal(iPoint,rot,1) = atl_physFluxRans(            &
        & state        = nodal_data(iPoint,rot),                         &
        & isenCoeff    = equation%euler%isen_coef,                       &
        & penalty_char = material%material_dat%elemMaterialData(1)       &
        &                                     %materialDat(iElem, 1, 1), &
        & porosity     = equation%euler%porosity                         )
    end do

    ! Calculate viscous physical flux point by point within this cell
    do iPoint = 1, poly_proj%body_3D%nquadpoints
      scheme_min%temp_nodal(iPoint,rot,2) = atl_viscPhysFluxRans( &
        & state          = nodal_data(iPoint,rot),                &
        & state_gradient = nodal_GradData(iPoint,rot, DerRot),    &
        & isenCoeff      = equation%euler%isen_coef,              &
        & mu             = equation%NavierStokes%mu,              &
        & lambda         = equation%NavierStokes%lambda,          &
        & thermCond      = equation%NavierStokes%therm_cond,      &
        & heatCap        = equation%euler%cv                      )
    end do

    ! Add up the nodal data
    nodal_res(:,:) = scheme_min%temp_nodal(:,:,1) &
      & - scheme_min%temp_nodal(:,:,2)

  end subroutine atl_modg_filNvrStk_physFlux_const


  !> Calculate the physical flux for the MODG scheme and
  !! Navier-Stokes equation (with non-constant penalizations).
  subroutine atl_modg_filNvrStk_physFlux_NonConst( equation, res, state, &
    & iElem, iDir, penalizationData, poly_proj, material, nodal_data,    &
    & nodal_GradData, nodal_res, elemLength, scheme_min, scheme_current  )
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
    integer :: rot(7), derRot(3)
    real(kind=rk) :: penalization(poly_proj%body_3D%nquadpoints)
    ! -------------------------------------------------------------------- !
    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    penalization = material%material_dat%elemMaterialData(2)    &
      &                                 %materialDat(iElem, :, 1)

    ! get the rotation for the physical flux calculation
    rot = equation%varRotation(iDir)                             &
                  %varTransformIndices(1:equation%varSys%nScalars)
    ! get rotation for the derivatives
    derRot = equation%varRotation(iDir)%derTransformIndices(2:4) &
      & - equation%varRotation(iDir)%derTransformIndices(1)


    do iPoint = 1, poly_proj%body_3D%nquadpoints
      scheme_min%temp_nodal(iPoint,rot,1) = atl_physFluxRans( &
        & state        = nodal_data(iPoint,rot),              &
        & isenCoeff    = equation%euler%isen_coef,            &
        & penalty_char = penalization(iPoint),                &
        & porosity     = equation%euler%porosity              )
    end do

    do iPoint = 1, poly_proj%body_3D%nquadpoints
      scheme_min%temp_nodal(iPoint,:,2) = atl_viscphysFluxRans( &
        & state          = nodal_data(iPoint,rot),              &
        & state_gradient = nodal_GradData(iPoint,rot,derRot),   &
        & isenCoeff      = equation%euler%isen_coef,            &
        & mu             = equation%NavierStokes%mu,            &
        & lambda         = equation%NavierStokes%lambda,        &
        & thermCond      = equation%NavierStokes%therm_cond,    &
        & heatCap        = equation%euler%cv                    )
    end do

    ! Add up the nodal data
    nodal_res(:,:) = scheme_min%temp_nodal(:,:,1) - scheme_min%temp_nodal(:,:,2)

  end subroutine atl_modg_filNvrStk_physFlux_nonconst


  !> Calculate the numerical flux for Navier-Stokes equation and MODG scheme
  subroutine atl_modg_filNvrStk_numFlux( mesh, equation, facedata, scheme, &
    &                                poly_proj, material )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    !> Parameter for used projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! --------------------------------------------------------------------------
    integer :: iDir
    ! --------------------------------------------------------------------------

    select case(trim(equation%FiltNavierStokes%model_type))
      case('rans')
        ! Numerical flux for faces in all 3 spatial face directions (x,y,z)
        do iDir = 1,3

          ! convective part of the Euler equations (F^*)
          ! ... fluxes for constant penalization parameters
          call atl_modg_euler_oneDim_numFlux_const(                            &
             & equation       = equation ,                                     &
             & nSides         = size(material%material_desc                    &
             &                       %computeFace(iDir,1)%leftPos),            &
             & faceRep        = facedata%faceRep(iDir)%dat,                    &
             & faceFlux       = facedata%faceFlux(iDir)%dat,                   &
             & leftPos        = material%material_desc                         &
             &                   %computeFace(iDir,1)%leftPos,                 &
             & rightPos       = material%material_desc                         &
             &                    %computeFace(iDir,1)%rightPos,               &
             & poly_proj      = poly_proj,                                     &
             & varRotation    = equation%varRotation(iDir)%varTransformIndices &
             &                       (1:equation%varSys%nScalars),             &
             & material_left  = material%material_dat%faceMaterialData(iDir,1) &
             &                   %leftElemMaterialDat,                         &
             & material_right = material%material_dat%faceMaterialData(iDir,1) &
             &                          %rightElemMaterialDat                  )


          ! ... fluxes for non-constant penalization parameters
          call atl_modg_euler_oneDim_numFlux_nonconst(                         &
             & equation       = equation ,                                     &
             & nSides         = size(material%material_desc                    &
             &                       %computeFace(iDir,2)%leftPos),            &
             & faceRep        = facedata%faceRep(iDir)%dat,                    &
             & faceFlux       = facedata%faceFlux(iDir)%dat,                   &
             & leftPos        = material%material_desc                         &
             &                           %computeFace(iDir,2)%leftPos,         &
             & rightPos       = material%material_desc                         &
             &                          %computeFace(iDir,2)%rightPos,         &
             & poly_proj      = poly_proj,                                     &
             & varRotation    = equation%varRotation(iDir)                     &
             &                  %varTransformIndices(1:equation%varSys         &
             &                                                 %nScalars),     &
             & material_left  = material%material_dat%faceMaterialData(iDir,2) &
             &                          %leftElemMaterialDat,                  &
             & material_right = material%material_dat%faceMaterialData(iDir,2)%&
             &                          rightElemMaterialDat                   )

          ! viscous numerical flux of the Navier-Stokes equation (sigma^*)
          call atl_modg_viscNavierStokes_oneDim_numFlux( equation = equation , &
             & facedata = facedata,                                            &
             & scheme = scheme ,                                               &
             & faces = mesh%faces%faces(iDir)%computeFace,                     &
             & faceDir = iDir, poly_proj = poly_proj, elemLength = mesh%length )

          ! stabilization viscous numerical flux of the Navier-Stokes
          ! equation (u^*)
          call atl_modg_stabViscNavierStokes_oneDim_numFlux( &
            & equation = equation,                           &
             & facedata = facedata,                          &
             & scheme = scheme ,                             &
             & faces = mesh%faces%faces(iDir)%computeFace,   &
             & faceDir = iDir, poly_proj = poly_proj         )

        end do

      case('les')
        write(logUnit(1),*) 'Not yet implemented'
        call tem_abort()


      case default
        write(logUnit(1),*) 'Turbulence model not defined ... stopping'
        call tem_abort()

    end select

  end subroutine atl_modg_filNvrStk_numFlux


end module atl_modg_filNvrStk_kernel_module
