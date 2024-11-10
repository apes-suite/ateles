! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> Module for routines and datatypes of Modal Discontinuous Galerkin (MODG)
!! scheme for the LinearEuler equation. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
module atl_modg_2d_LinearEuler_kernel_module
  use env_module,               only: rk

  use ply_poly_project_module,  only: ply_poly_project_type, assignment(=)
  use ply_dof_module,           only: Q_space, P_space

  use atl_equation_module,      only: atl_equations_type
  use atl_facedata_module,      only: atl_facedata_type
  use atl_cube_elem_module,     only: atl_cube_elem_type
  use atl_scheme_module,        only: atl_scheme_type
  use atl_modg_2d_scheme_module,only: atl_modg_2d_scheme_type
  use atl_LinearEuler_2d_physFlux_module, only: atl_LinearEuler_2d_physFlux
  use atl_penalization_module,          only: atl_penalizationData_type
  use atl_materialPrp_module,           only: atl_material_type

  implicit none

  private

  public :: atl_modg_2d_LinearEuler_numflux, atl_modg_2d_LinearEuler_physFlux


contains


  ! ****************************************************************************
  !> Calculate the physical flux for the MODG scheme and
  !! Linearized euler equation.
  subroutine atl_modg_2d_LinearEuler_physFlux( equation, res, state, iElem,    &
    & iDir, penalizationData, poly_proj, material, nodal_data, nodal_gradData, &
    & nodal_res, elemLength,  scheme_min, scheme_current                       )

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
    ! Loop var for all the dof in an element
    integer :: iDof, nDofs
    ! Rotation indices for physical flux calculation
    integer :: rot(4)
    ! --------------------------------------------------------------------------!

    ! get the rotation for the physical flux calculation
    rot = equation%varRotation(iDir)%varTransformIndices(1:4)
    nDofs = poly_proj%body_2d%ndofs

    !This subroutine is being called inside a parallel region
    dofLoop: do iDof = 1, ndofs

        ! Calculate the physical flux point by point within this cell
        res(iDof,rot) = atl_LinearEuler_2d_physFlux(                        &
          &                         state = state(iDof,rot),                &
          &                         LinearEuler = equation%LinearEuler,     &
          &                         idir = iDir                             )

    end do dofLoop


  end subroutine atl_modg_2d_LinearEuler_physFlux
  ! ****************************************************************************


  ! ****************************************************************************
  !> Calculate the numerical flux for LinearEuler equation and MODG scheme
  subroutine atl_modg_2d_LinearEuler_numFlux( mesh, equation, facedata, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------
    integer :: iDir, nFaceDofs
    ! --------------------------------------------------------------------------

    ! Numerical flux for faces in all 2 spatial face directions
    select case(scheme%basisType)
      case(Q_space)
        nFaceDofs = (scheme%maxPolyDegree+1)
      case(P_space)
  nfacedofs = ((scheme%maxpolydegree)+1)
    end select


    ! Calculate the numerical fluxes for the faces in all 2 spatial face
    ! directions
    do iDir = 1,2
      call equation%LinearEuler%dir_proc(iDir)%numFlux(               &
        &  nSides = size(mesh%faces%faces(iDir)%computeFace%leftPos), &
        &  nFaceDofs = nFaceDofs,                                     &
        &  faceRep  = facedata%faceRep(iDir)%dat,                     &
        &  faceFlux = facedata%faceFlux(iDir)%dat,                    &
        &  leftPos  = mesh%faces%faces(iDir)%computeFace%leftPos,     &
        &  rightPos = mesh%faces%faces(iDir)%computeFace%rightPos,    &
        &  var = equation%varRotation(iDir)%varTransformIndices(1:4), &
        &  LinearEuler = equation%LinearEuler ,                       &
        &  iDir = iDir                                                )
    end do


  end subroutine atl_modg_2d_LinearEuler_numFlux
  ! ****************************************************************************

end module atl_modg_2d_LinearEuler_kernel_module
