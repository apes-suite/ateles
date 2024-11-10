! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014, 2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
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

!> Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme for the Heat equation. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
module atl_modg_heat_kernel_module
  use env_module,                      only: rk

  use atl_equation_module,             only: atl_equations_type
  use atl_cube_elem_module,            only: atl_cube_elem_type
  use atl_scheme_module,               only: atl_modg_scheme_type, &
    &                                        atl_scheme_type
  use atl_facedata_module,             only: atl_facedata_type
  use atl_numFluxHeat_module,          only: atl_modg_heat_numFlux_sipg
  use atl_penalization_module,         only: atl_penalizationData_type
  use atl_materialPrp_module,          only: atl_material_type

  use ply_poly_project_module,         only: ply_poly_project_type,assignment(=)
  use ply_leg_diff_module,             only: ply_calcDiff_leg_normal
  use ply_dof_module,                  only: ply_change_poly_space, &
    &                                        P_space,               &
    &                                        Q_space


  implicit none

  private

  public :: atl_modg_heat_numflux, atl_modg_heat_physFlux

contains

  !> Calculate the physical flux for the MODG scheme and
  !! Heat equation.
  subroutine atl_modg_heat_physFlux( equation, res, state, iElem, iDir,     &
    &                                penalizationData, poly_proj, material, &
    &                                nodal_data,nodal_gradData, nodal_res,  &
    &                                elemLength, scheme_min, scheme_current )

    ! --------------------------------------------------------------------------
    !> The equation system we are working with
    type(atl_equations_type), intent(in) :: equation
    !> The result in the modal form
    real(kind=rk), intent(inout)     :: res(:,:)
    !> The state in the modal form
    real(kind=rk), intent(in), optional :: state(:,:)
    !> The current element index
    integer, intent(in) :: iElem
    !> The current direction
    integer, intent(in) :: iDir
    !> The Penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> The projection datatype for the projection information
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The material information
    type(atl_material_type), intent(inout) :: material
    !> The state data in the nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    real(kind=rk), intent(in), optional :: nodal_GradData(:,:,:)
    !> The result in the nodal form
    real(kind=rk), intent(inout)     :: nodal_res(:,:)
    !> The length of the current element
    real(kind=rk), intent(in) :: ElemLength
    !> The scheme information of the min level (This is needed for the temp
    ! buffer array for evaluating the physical fluxes )
    type(atl_scheme_type), intent(inout) :: scheme_min
    !> Information about the current level
    type(atl_scheme_type), intent(inout) :: scheme_current
    ! --------------------------------------------------------------------------!
    integer :: nScal, nDofs
    real(kind=rk), allocatable ::  temp_modal_P(:,:,:,:)
    real(kind=rk), allocatable ::  temp_modal_Q(:,:,:,:)
    real(kind=rk) ::  therm_diff
    ! --------------------------------------------------------------------------!

    therm_diff = equation%heat%k
    nScal = equation%varSys%nScalars
    nDofs = poly_proj%body_3d%ndofs

    ! get the modal coefficients of the current cell
    ! ATTENTION: have to be duplicated as the FPT modifies the input vector.
    scheme_min%temp_modal(:ndofs,:nScal,1) = state(:,:)

    select case (scheme_min%modg%basisType)
    case (Q_space)
      call ply_calcDiff_leg_normal(                                     &
        &       legCoeffs     = scheme_min%temp_modal(:ndofs,:nScal,1), &
        &       legcoeffsDiff = scheme_min%temp_modal(:ndofs,:nScal,2), &
        &       mPd           = poly_proj%maxPolyDegree,                &
        &       elemLength    = elemLength,                             &
        &       nVars         = 1,                                      &
        &       iDir          = iDir                                    )

    case (P_space)
      allocate(temp_modal_P(1,nDofs,nScal,2))
      allocate(temp_modal_Q(1,(poly_proj%maxPolyDegree+1)**3,nScal,2))

      temp_modal_P(1,:ndofs,:nscal,:2) = scheme_min%temp_modal(:,:,:)

      call ply_change_poly_space( inspace    = P_space,                 &
        &                         instate    = temp_modal_P(:,:,:,1),   &
        &                         outstate   = temp_modal_Q(:,:,:,1),   &
        &                         maxPolyDeg = poly_proj%maxPolyDegree, &
        &                         nElems     = 1,                       &
        &                         nVars      = nScal,                   &
        &                         nDims      = 3                        )

      call ply_calcDiff_leg_normal(                      &
        &       legCoeffs     = temp_modal_Q(1,:,:,1),   &
        &       legcoeffsDiff = temp_modal_Q(1,:,:,2),   &
        &       mPd           = poly_proj%maxPolyDegree, &
        &       elemLength    = elemLength,              &
        &       nVars         = 1,                       &
        &       iDir          = iDir                     )

      call ply_change_poly_space( inspace    = Q_space,                 &
        &                         instate    = temp_modal_Q(:,:,:,2),   &
        &                         outstate   = temp_modal_P(:,:,:,2),   &
        &                         maxPolyDeg = poly_proj%maxPolyDegree, &
        &                         nElems     = 1,                       &
        &                         nVars      = nScal,                   &
        &                         nDims      = 3                        )

      scheme_min%temp_modal(:ndofs,:nScal,2) = temp_modal_P(1,:,:,2)

      deallocate(temp_modal_P)
      deallocate(temp_modal_Q)

    end select

    res(:,1) = -therm_diff*scheme_min%temp_modal(:,1,2)

  end subroutine atl_modg_heat_physFlux

  !> Calculate the numerical flux for Heat equation and MODG scheme
  subroutine atl_modg_heat_numFlux( mesh, equation, facedata, scheme, &
    &                               poly_proj                         )
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
    ! --------------------------------------------------------------------------
    integer :: iDir
    ! --------------------------------------------------------------------------

    ! Numerical flux for faces in all spatial face directions

    !%OMP PRIVATE(iDir) &
    do iDir = 1,3
      call atl_modg_heat_numFlux_sipg(                     &
        & equation   = equation ,                          &
        & facedata   = facedata,                           &
        & faces      = mesh%faces%faces(iDir)%computeFace, &
        & faceDir    = iDir,                               &
        & dofs       = poly_proj%body_2d%ndofs,            &
        & elem_len   = mesh%length,                        &
        & maxPolyDeg = scheme%maxPolyDegree                )
    end do

  end subroutine atl_modg_heat_numFlux


end module atl_modg_heat_kernel_module
