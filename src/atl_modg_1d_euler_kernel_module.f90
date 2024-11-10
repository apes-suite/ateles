! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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
!! scheme for the 1D Euler equation.
!! This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
module atl_modg_1d_euler_kernel_module
  use env_module,                   only: rk
  use treelmesh_module,             only: treelmesh_type
  use tem_faceData_module,          only: tem_faceIterator_type
  use treelmesh_module,             only: treelmesh_type

  use atl_equation_module,          only: atl_equations_type
  use atl_materialPrp_module,       only: atl_material_type, &
    &                                     atl_faceMaterialData_type
  use atl_scheme_module,            only: atl_scheme_type
  use atl_physFluxEuler_1d_module,  only: atl_physFluxEuler_1d
  use atl_facedata_module,          only: atl_facedata_type
  use atl_penalization_module,      only: atl_penalizationData_type
  use ply_poly_project_module,      only: ply_poly_project_type, &
    &                                     assignment(=), &
    &                                     ply_poly_project_n2m
  use ply_oversample_module,        only: ply_convertFromoversample

  implicit none

  private

  public :: atl_modg_1d_euler_numflux
  public :: atl_modg_1d_euler_physFlux_const
  public :: atl_modg_1d_euler_physFlux_nonConst
  public :: atl_modg_1d_euler_penalization_const
  public :: atl_modg_1d_euler_penalization_nonConst


contains


  ! ------------------------------------------------------------------------ !
  !> Calculate the physical flux for the MODG scheme and
  !! Euler equation with constant characteristic (mask function) in the
  !! element.
  subroutine atl_modg_1d_euler_physFlux_const( equation, poly_proj,         &
    &                                          res, pointVal, penalty_char, &
    &                                          U_o                          )
    ! -------------------------------------------------------------------- !
    !> The equation you solve.
    type(atl_equations_type) :: equation

    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> The physical flux result to be stored in
    real(kind=rk), intent(inout) :: res(:,:)

    !> Nodal representation of the polynomial with in each cell.
    real(kind=rk), intent(in) :: pointVal(:,:)

    !> Characteristic (mask) function of the material penalization.
    !!
    !! This should be 0 everywhere where there is no material, and
    !! 1 where there is material.
    real(kind=rk), intent(in) :: penalty_char
    !> Obstacle velocity
    real(kind=rk), intent(in) :: U_o
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    integer :: nquadpoints
    real(kind=rk) :: penalty_scaling
    ! -------------------------------------------------------------------- !

    penalty_scaling = (1.0_rk/equation%euler%porosity - 1.0_rk)*penalty_char

    ! get correct amount of quadrature points and degree
    nQuadpoints = poly_proj%body_1d%nquadpoints


    ! Calculate the physical flux point by point within this cell - x direction
    do iPoint = 1, nQuadPoints
      Res(iPoint,:) = atl_physFluxEuler_1d(           &
        & state           = pointVal(iPoint,:),       &
        & isenCoeff       = equation%euler%isen_coef, &
        & U_o             = U_o,                      &
        & penalty_scaling = penalty_scaling           )
    end do

  end subroutine atl_modg_1d_euler_physFlux_const
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Calculate the physical flux for the MODG scheme and
  !! Euler equation with variable characteristic (mask function) in the
  !! element.
  subroutine atl_modg_1d_euler_physFlux_nonConst( equation, poly_proj,         &
    &                                             res, pointVal, penalty_char, &
    &                                             U_o                          )
    ! -------------------------------------------------------------------- !
    !> The equation you solve.
    type(atl_equations_type) :: equation

    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> The physical flux result to be stored in
    real(kind=rk), intent(inout) :: res(:,:)

    !> Nodal representation of the polynomial with in each cell.
    real(kind=rk), intent(in) :: pointVal(:,:)

    !> Characteristic (mask) function of the material penalization.
    !!
    !! This should be 0 everywhere where there is no material, and
    !! 1 where there is material.
    real(kind=rk), intent(in) :: penalty_char(:)
    !! obstacle velocity
    real(kind=rk), intent(in) :: U_o(:)
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    integer :: nQuadpoints
    real(kind=rk) :: porosity_fact
    real(kind=rk) :: penalty_scaling
    ! -------------------------------------------------------------------- !

    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nQuadpoints = poly_proj%body_1d%nquadpoints

    porosity_fact = 1.0_rk/equation%euler%porosity - 1.0_rk

    ! Calculate the physical flux point by point within this cell - x direction
    do iPoint = 1, nQuadPoints
      penalty_scaling =  porosity_fact*penalty_char(iPoint)
      Res(iPoint,:) = atl_physFluxEuler_1d(           &
        & state           = pointVal(iPoint,:),       &
        & isenCoeff       = equation%euler%isen_coef, &
        & U_o             = U_o(iPoint),              &
        & penalty_scaling = penalty_scaling           )
    end do

  end subroutine atl_modg_1d_euler_physFlux_nonConst
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Calculate the numerical flux for Euler equation and MODG scheme
  subroutine atl_modg_1d_euler_numFlux( equation, facedata, material )
    ! -------------------------------------------------------------------- !
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation

    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata

    !> Material description for the faces on the current level
    type(atl_material_type), intent(in) :: material
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ! Numerical flux for faces in all 1 spatial face directions (x dir)
    ! Constant material
    call modg_1d_euler_oneDim_numFlux( equation = equation,                  &
      &                                facedata = facedata,                  &
      &                                faces    = material%material_desc     &
      &                                                   %computeFace(1,1), &
      &                                material = material%material_dat      &
      &                                           %faceMaterialData(1,1)     )

    ! Variable material
    call modg_1d_euler_oneDim_numFlux( equation = equation,                  &
      &                                facedata = facedata,                  &
      &                                faces    = material%material_desc     &
      &                                                   %computeFace(1,2), &
      &                                material = material%material_dat      &
      &                                           %faceMaterialData(1,2)     )

  end subroutine atl_modg_1d_euler_numFlux
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Numerical flux calculation for Euler equation across the faces in a single
  !! spatial direction.
  subroutine modg_1d_euler_oneDim_numFlux( equation, facedata, faces, &
    &                                      material )
    ! -------------------------------------------------------------------- !
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation

    !> The face state if the equation
    type(atl_facedata_type), intent(inout) :: facedata

    !> The faces to calculate the fluxes for.
    type(tem_faceIterator_type), intent(in) :: faces

    !> The penalization material left and right of the face
    type(atl_facematerialData_type), intent(in) :: material
    ! -------------------------------------------------------------------- !
    ! Nodal representation of the numerical flux
    real(kind=rk) :: nodalNumFlux(1, equation%varSys%nScalars)
    ! Loop over variables (due to single variable FPTs)
    integer :: nSides
    integer :: iSide
    integer :: left_neighbor, right_neighbor
    ! -------------------------------------------------------------------- !

    !!nSides = size(faces%leftPos)


    ! Use the numerical flux set by the equation system.
    ! Note, the rotation of the input state, the output is not rotated back
    ! here yet, as this is not allowed for output variables.

    nSides = size(faces%leftPos)
    do iSide=1,nSides
      left_neighbor = faces%leftpos(iSide)
      right_neighbor = faces%rightpos(iSide)
      call equation%Euler%numflux( &
        &  state_left     = facedata%faceRep(1)%dat(left_neighbor,1:1,:,2),  &
        &  state_right    = facedata%faceRep(1)%dat(right_neighbor,1:1,:,1), &
        &  material_left  = material%leftElemMaterialDat(iSide,:,:),         &
        &  material_right = material%rightElemMaterialDat(iSide,:,:),        &
        &  nPoints        = 1,                                               &
        &  flux           = nodalNumFlux                                     )

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      facedata%faceFlux(1)%dat(left_neighbor,1,:,2)  = nodalnumFlux(1,:)
      facedata%faceFlux(1)%dat(right_neighbor,1,:,1) = nodalnumFlux(1,:)
    end do


  end subroutine modg_1d_euler_oneDim_numFlux
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! Calculate the penalization terms (for density, momentum, energy)
  ! The penalization terms are calculated in the sammer manner as
  ! the physical fluxes, i.e. in a nodal-modal approach
  subroutine atl_modg_1d_euler_penalization_Const( equation, poly_proj, &
    & nodal_data, scheme_min, penalizationData, iElem, material         )
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
    type(atl_material_type), intent(in) :: material
    ! -------------------------------------------------------------------- !
    integer :: iPoint, glob_elem_ind
    real(kind=rk) :: temperature, pressure
    real(kind=rk) :: penalization( 3 )
    ! -------------------------------------------------------------------- !
    glob_elem_ind = material%material_desc%computeElems(1)%totElemIndices(iElem)
    penalization = material%material_dat%elemMaterialData(1) &
      &                                 %materialDat(iElem,1,:)

    ! Calculate the penalization terms (for density, momentum, energy)
    ! The penalization terms are calculated in the sammer manner as
    ! the physical fluxes,
    ! i.e. in a nodal-modal approach

    ! Penalization for momentum and energy in pointwise manner

    do iPoint = 1, poly_proj%body_1D%nquadpoints

      pressure = (equation%euler%isen_coef-1.0_rk)   &
        &       * ( nodal_data(iPoint,3)             &
        &           - 0.5_rk*nodal_data(iPoint,2)**2 &
        &                   / nodal_data(iPoint,1)   )
      temperature =  pressure / ( nodal_data(iPoint,1) * equation%euler%R )

      ! ... momentum
      scheme_min%temp_nodal(iPoint,1,1)                              &
        &  = (-1.0_rk) * penalization(1)                             &
        &              * ( nodal_data(iPoint,2)/nodal_data(iPoint,1) &
        &                  - penalization(2) )                       &
        &              / equation%euler%viscous_permeability
      ! ... energy
      scheme_min%temp_nodal(iPoint,2,1)                    &
        &  = (-1.0_rk) * penalization(1)                   &
        &              * ( temperature - penalization(3) ) &
        &              / equation%euler%thermal_permeability

    end do

    ! Transform vector of penalization data here of size (nScalars-1)
    ! back to modal space
    call ply_poly_project_n2m( me         = poly_proj,                    &
      &                        dim        = 1,                            &
      &                        nVars      = equation%varSys%nScalars-1,   &
      &                        nodal_data = scheme_min%temp_nodal(:,:2,1), &
      &                        modal_data = scheme_min%temp_over(:,:2,1)   )

    ! ... no penalization for density
    penalizationdata%penalization_data(glob_elem_ind,:,1) = 0.0_rk

    ! --> oversamp modal space for penalty terms (momentum + energy)
    call ply_convertFromOversample(                               &
      &    state       = penalizationdata                         &
      &                  %penalization_data(glob_elem_ind,:,2:3), &
      &    poly_proj   = poly_proj,                               &
      &    nDim        = 1,                                       &
      &    modalCoeffs = scheme_min%temp_over(:,:2,1)              )

  end subroutine atl_modg_1d_euler_penalization_Const
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! Calculate the penalization terms (for density, momentum, energy)
  ! The penalization terms are calculated in the sammer manner as
  ! the physical fluxes, i.e. in a nodal-modal approach
  subroutine atl_modg_1d_euler_penalization_NonConst(equation, poly_proj, &
    & nodal_data, scheme_min, penalizationData, iElem, material           )
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
    type(atl_material_type), intent(in) :: material
    ! -------------------------------------------------------------------- !
    integer :: iPoint, glob_elem_ind
    real(kind=rk) :: temperature, pressure
    real(kind=rk) :: penalization( poly_proj%body_1D%nquadpoints, 3 )
    ! -------------------------------------------------------------------- !
    glob_elem_ind = material%material_desc%computeElems(2)%totElemIndices(iElem)
    penalization = material%material_dat%elemMaterialData(2) &
      &                                 %materialDat(iElem,:,:)

    ! Calculate the penalization terms (for density, momentum, energy)
    ! The penalization terms are calculated in the sammer manner as
    ! the physical fluxes,
    ! i.e. in a nodal-modal approach

    ! Penalization for momentum and energy in pointwise manner

    do iPoint = 1, poly_proj%body_1D%nquadpoints

      pressure = (equation%euler%isen_coef-1.0_rk)   &
        &       * ( nodal_data(iPoint,3)             &
        &           - 0.5_rk*nodal_data(iPoint,2)**2 &
        &                   / nodal_data(iPoint,1)   )
      temperature =  pressure / ( nodal_data(iPoint,1) * equation%euler%R )

      ! ... momentum
      scheme_min%temp_nodal(iPoint,1,1)                              &
        &  = (-1.0_rk) * penalization(iPoint,1)                      &
        &              * ( nodal_data(iPoint,2)/nodal_data(iPoint,1) &
        &                  - penalization(iPoint,2) )                &
        &              / equation%euler%viscous_permeability
      ! ... energy
      scheme_min%temp_nodal(iPoint,2,1)                           &
        &  = (-1.0_rk) * penalization(iPoint,1)                   &
        &              * ( temperature - penalization(iPoint,3) ) &
        &              / equation%euler%thermal_permeability

    end do

    ! Transform vector of penalization data here of size (nScalars-1)
    ! back to modal space
    call ply_poly_project_n2m( me         = poly_proj,                     &
      &                        dim        = 1,                             &
      &                        nVars      = equation%varSys%nScalars-1,    &
      &                        nodal_data = scheme_min%temp_nodal(:,:2,1), &
      &                        modal_data = scheme_min%temp_over(:,:2,1)   )

    ! ... no penalization for density
    penalizationdata%penalization_data(glob_elem_ind,:,1) = 0.0_rk

    ! --> oversamp modal space for penalty terms (momentum + energy)
    call ply_convertFromOversample(                               &
      &    state       = penalizationdata                         &
      &                  %penalization_data(glob_elem_ind,:,2:3), &
      &    poly_proj   = poly_proj,                               &
      &    nDim        = 1,                                       &
      &    modalCoeffs = scheme_min%temp_over(:,:2,1)              )

  end subroutine atl_modg_1d_euler_penalization_NonConst

end module atl_modg_1d_euler_kernel_module
