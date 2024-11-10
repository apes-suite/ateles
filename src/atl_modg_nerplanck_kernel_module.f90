! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012-2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
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
!> Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
module atl_modg_nerplanck_kernel_module
  use env_module,               only: rk, long_k
  use tem_element_module,       only: eT_fluid
  use tem_element_module,       only: eT_fluid

  use atl_equation_module,       only: atl_equations_type
  use atl_kerneldata_module,     only: atl_kerneldata_type, atl_statedata_type
  use atl_cube_elem_module,      only: atl_cube_elem_type
  use atl_nerplanck_flux_module, only: atl_nerplanck_numflux_solve,      &
    &                                  atl_nerplanck_physflux_solve,     &
    &                                  atl_nerplanck_numflux_preprocess, &
    &                                  atl_nerplanck_physflux_preprocess
  use atl_scheme_module,         only: atl_modg_scheme_type
  use ply_modg_basis_module,     only: ply_faceValLeftBndAns, &
    &                                  ply_faceValLeftBndTest, &
    &                                  ply_faceValRightBndTest, &
    &                                  ply_scalProdDualLeg,     &
    &                                  ply_scalProdDualLegDiff

  implicit none
  private

  public :: atl_preprocess_modg_nerplanck_kernel

  !> interface for preprocessing the data for the kernel
  interface atl_preprocess_modg_nerplanck_kernel
    module procedure preprocess_modg_nerplanck_kernel
  end interface atl_preprocess_modg_nerplanck_kernel

contains

  subroutine preprocess_modg_nerplanck_kernel( equation, kerneldata, statedata, &
                                   & mesh, scheme)
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The equation you solve.
    type(atl_equations_type) :: equation
    !> The data of the kernel.
    type(atl_kerneldata_type) :: kerneldata
    !> THe state if the equation
    type(atl_statedata_type) :: statedata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type) :: scheme
    ! --------------------------------------------------------------------------

    ! Initialize the workspace for the modg kernel with zeros. The following routines
    ! will add the influence of their projections to this array.
    kerneldata%state_der = 0.0_rk

    ! Loop over all the elements and calculate the projection of
    ! the cell internal physical flux. Therefore, we do the following for each cell:
    ! First, we loop over the modal dofs and create a modal representation of the
    ! physical flux.
    ! Second, we loop over the test functions and project the modal representation
    ! of the physical flux onto them.

    ! add physical flux calculation.
    call modg_nerplanck_physFlux_pre( mesh, equation, statedata, scheme )

    ! Loop over all the faces and calculate the projection of the
    ! numerical flux. Therefore, we do the following for each faces:
    ! First, we loop over all modal dofs and calculate the modal representation
    ! of the flux on this face.
    ! Second, we loop over all test functions and project the modal representation
    ! of the numerical flux for this face onto the test functions.

    ! ... numerical flux across faces in x direction
    call modg_nerplanck_x_numFlux_pre( mesh, equation, statedata, scheme )

    ! ... numerical flux across faces in y direction
    call modg_nerplanck_y_numFlux_pre( mesh, equation, statedata, scheme )

    ! ... numerical flux across faces in z direction
    call modg_nerplanck_z_numFlux_pre( mesh, equation, statedata, scheme )

    ! Finally, we loop over all the cells and apply the cell local inverse mass matrix.
    ! All of the previous steps and this one include the correct Determinant of the
    ! Jacobians of the mappings from reference element/face to the physical element/face.
    ! It also includes the elemental timestepping of the time_integration_module.
    call modg_invMassMatrix_pre( mesh, scheme )

  end subroutine preprocess_modg_nerplanck_kernel


  !> Numerical flux calculation for Nernst-Planck equation across the faces in X direction.
  subroutine modg_nerplanck_x_numFlux_pre( mesh, equation, statedata, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The equation you solve.
    type(atl_equations_type) :: equation
    !> THe state if the equation
    type(atl_statedata_type) :: statedata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------!
    ! Loop vars
    integer :: xAnsFunc, yAnsFunc, zAnsFunc, xTestFunc, yTestFunc, zTestFunc
    ! Positions for the modal coefficients
    integer :: modalCoeff, testPos, ansPos
    ! Modal coefficients for elements left and right of the face
    real(kind=rk) :: leftModalCoeff(4), rightModalCoeff(4), flux(3)
    ! Modal representation of the flux on the face
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    ! The values of the x ansatz and test functions at the face from left and right
    ! limit.
    real(kind=rk) :: faceValLeft, faceValRight
    ! Loop var for the faces.
    integer :: xside
    ! Element positions of the left and right element of the face.
    integer(kind=long_k) :: left_neighbor, right_neighbor
    ! Speed of light, i.e. maximum information propagation speed.
    real(kind=rk) :: diffusivity
    ! Scalar products between y and z ansatz and test functions.
    real(kind=rk) :: yScalProd, zScalProd
    ! Determinant of the jacobian of the mapping from reference face to
    ! the current face.
    real(kind=rk) :: jacobiDet
    ! The minimum indices for the relevant ansatz functions.
    integer :: zAnsFuncMin, yAnsFuncMin
    ! --------------------------------------------------------------------------

    allocate(numFluxBuffer((scheme%maxPolyDegree+1)**3,3))

    ! Calculate the speed of light
    diffusivity = equation%nerplanck%D

    ! Build the determinant of the Jacobian of the mapping from reference
    ! face to the physical face. Since we are considering cubes it is
    ! the same everywhere.
    jacobiDet = (mesh%length/2.0_rk)**2.0_rk

    ! Loop over all fluid the faces in x direction
    XFaceLoop: do xside = 1, size(mesh%faces%faces(1)%computeFace%leftPos)
      ! Get the fluid neighbors for this face.
      left_neighbor = mesh%faces%faces(1)%computeFace%leftPos(xside)
      right_neighbor = mesh%faces%faces(1)%computeFace%rightPos(xside)

      ! Loop over all the ansatz functions. Here, we use simple tensor products
      ! of one-dimensional ansatz functions.
      do xAnsFunc = 1, scheme%maxPolyDegree+1
        faceValRight = ply_faceValLeftBndAns(xAnsFunc)
        do yAnsFunc = 1, (scheme%maxPolyDegree+1)
          do zAnsFunc = 1, (scheme%maxPolyDegree+1)
            ! Get the modal coefficients of the left and right ansatz functions and
            ! calculate the modal flux with them.
  modalcoeff = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (zansfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
            leftModalCoeff(:) = statedata%state(left_neighbor, modalCoeff, :)
            rightModalCoeff(:) = statedata%state(right_neighbor, modalCoeff, :) * faceValRight
            call atl_nerplanck_numflux_preprocess( &
              & left        = leftModalCoeff,      &
              & right       = rightModalCoeff,     &
              & diffusivity = diffusivity,         &
              & flux        = flux                 )
            ! Buffer the modal values of the numerical flux
            numFluxBuffer( modalCoeff, :) = flux
          end do
        end do
      end do

      ! Loop over all the test functions and project the numerical flux to them.
      do xTestFunc = 1,min(scheme%maxPolyDegree+1,2)
        faceValLeft = ply_faceValRightBndTest(xTestFunc)
        faceValRight = ply_faceValLeftBndTest(xTestFunc)
        do yTestFunc = 1,scheme%maxPolyDegree+1
          do zTestFunc = 1,scheme%maxPolyDegree+1

            ! position of the test functions in the kerneldata
  testpos = xtestfunc                                      &
    &      + ( ( ytestfunc-1)                             &
    &      + (ztestfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)

            ! get the relevant indices of ansatz functions for the projection
            do xAnsFunc = 1, scheme%maxPolyDegree+1
              yAnsFuncMin = yTestFunc-2
              if( yAnsFuncMin < 1 ) then
                yAnsFuncMin = yTestFunc
              end if
              do yAnsFunc = yAnsFuncMin,yTestFunc,2
                zAnsFuncMin = zTestFunc-2
                if( zAnsFuncMin < 1 ) then
                  zAnsFuncMin = zTestFunc
                end if
                do zAnsFunc = zAnsFuncMin,zTestFunc,2

                  ! the position of the modal coefficeint of this ansatz functions
  anspos = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (zansfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)

                  ! calculate the projection of the ansatz and test function
                  yScalProd = ply_scalProdDualLeg(yAnsFunc, yTestFunc)
                  zScalProd = ply_scalProdDualLeg(zAnsFunc, zTestFunc)

!NA!                  ! buffer the result in kernel data and take care of the outer surface unit normal
!NA!                  ! ... for the element left of the face
!NA!                  outerNormalLeft = 1.0_rk
!NA!                  projectedFlux = outerNormalLeft * numFluxBuffer(ansPos,:)
!NA!                  kernelData%state_der(left_neighbor, testPos, 1, 2:4) &
!NA!                            & = kernelData%state_der(left_neighbor, testPos, 1, 2:4) &
!NA!                            & - faceValLeft * yScalProd * zScalProd &
!NA!                            & * projectedFlux * jacobiDet
!NA!                  ! ... for the element right of the face
!NA!                  outerNormalRight = -1.0_rk
!NA!                  projectedFlux = outerNormalRight * numFluxBuffer(ansPos,:)
!NA!                  kernelData%state_der(right_neighbor, testPos,1,2:4) &
!NA!                            & = kernelData%state_der(right_neighbor, testPos, 1, 2:4) &
!NA!                            & - faceValRight * yScalProd * zScalProd &
!NA!                            & * projectedFlux * jacobiDet

                end do
              end do
            end do

          end do
        end do
      end do

    end do XFaceLoop

    !GhostFaceLoop: do xside = 1, mesh%faces(1)%ghostfromfiner%leftTreeId%nVals


  end subroutine modg_nerplanck_x_numFlux_pre

  !> Numerical flux calculation for Nernst-Planck equation across the faces in Y direction.
  subroutine modg_nerplanck_y_numFlux_pre( mesh, equation, statedata, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The equation you solve.
    type(atl_equations_type) :: equation
    !> THe state if the equation
    type(atl_statedata_type) :: statedata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------!
    ! Loop vars
    integer :: xAnsFunc, yAnsFunc, zAnsFunc, xTestFunc, yTestFunc, zTestFunc
    ! Positions for the modal coefficients
    integer :: modalCoeff, testPos, ansPos
    ! Modal coefficients for elements left and right of the face
    real(kind=rk) :: leftModalCoeff(4), rightModalCoeff(4), flux(3)
    ! Modal representation of the flux on the face
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    ! The values of the x ansatz and test functions at the face from left and right
    ! limit.
    real(kind=rk) :: faceValLeft, faceValRight
    ! Loop var for the faces.
    integer :: yside
    ! Element positions of the left and right element of the face.
    integer(kind=long_k) :: left_neighbor, right_neighbor
    ! Speed of light, i.e. maximum information propagation speed.
    real(kind=rk) :: diffusivity
    ! Scalar products between x and z ansatz and test functions.
    real(kind=rk) :: xScalProd, zScalProd
    ! Determinant of the jacobian of the mapping from reference face to
    ! the current face.
    real(kind=rk) :: jacobiDet
    ! The minimum indices for the relevant ansatz functions.
    integer :: zAnsFuncMin, xAnsFuncMin
    ! Rotation for the flux
    integer :: varRotation(4)
    ! --------------------------------------------------------------------------

    allocate(numFluxBuffer((scheme%maxPolyDegree+1)**3,3))

    ! Calculate the speed of light
    diffusivity = equation%nerplanck%D

    ! Build the determinant of the Jacobian of the mapping from reference
    ! face to the physical face. Since we are considering cubes it is
    ! the same everywhere.
    jacobiDet = (mesh%length/2.0_rk)**2.0_rk

    ! Loop over fluid the faces in y direction
    varRotation = equation%varRotation(2)%varTransformIndices(:)
    yFaceLoop: do yside = 1,size(mesh%faces%faces(2)%computeFace%leftPos)
      ! Get the fluid neighbors for this face.
      left_neighbor = mesh%faces%faces(2)%computeFace%leftPos(yside)
      right_neighbor = mesh%faces%faces(2)%computeFace%rightPos(yside)


      ! VK: old face descriptor
      ! Get the neighbors for this face.
      !left_neighbor  = mesh%side%ydir%left_elemid(yside)
      !right_neighbor = mesh%side%ydir%right_elemid(yside)

      ! Loop over all the ansatz functions. Here, we use simple tensor products
      ! of one-dimensional ansatz functions.
      do yAnsFunc = 1, scheme%maxPolyDegree+1
        faceValRight = ply_faceValLeftBndAns(yAnsFunc)
        do xAnsFunc = 1, (scheme%maxPolyDegree+1)
          do zAnsFunc = 1, (scheme%maxPolyDegree+1)
            ! Get the modal coefficients of the left and right ansatz functions and
            ! calculate the modal flux with them.
  modalcoeff = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (zansfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
            ! ... keep in mind that we have to rotate the input and output of the flux
            leftModalCoeff(:) = statedata%state(left_neighbor, &
                                               & modalCoeff, varRotation)
            rightModalCoeff(:) = statedata%state(right_neighbor, &
                                               &  modalCoeff, varRotation) * faceValRight
            call atl_nerplanck_numflux_preprocess( &
              & left        = leftModalCoeff,      &
              & right       = rightModalCoeff,     &
              & diffusivity = diffusivity,         &
              & flux        = flux                 )
            ! Buffer the modal values of the numerical flux
            numFluxBuffer( modalCoeff, :) = flux
          end do
        end do
      end do

      ! Loop over all the test functions and project the numerical flux to them.
      do yTestFunc = 1,min(scheme%maxPolyDegree+1,2)
        faceValLeft = ply_faceValRightBndTest(yTestFunc)
        faceValRight = ply_faceValLeftBndTest(yTestFunc)
        do xTestFunc = 1,scheme%maxPolyDegree+1
          do zTestFunc = 1,scheme%maxPolyDegree+1

            ! position of the test functions in the kerneldata
  testpos = xtestfunc                                      &
    &      + ( ( ytestfunc-1)                             &
    &      + (ztestfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)

            ! get the relevant indices of ansatz functions for the projection
            do yAnsFunc = 1, scheme%maxPolyDegree+1
              xAnsFuncMin = xTestFunc-2
              if( xAnsFuncMin < 1 ) then
                xAnsFuncMin = xTestFunc
              end if
              do xAnsFunc = xAnsFuncMin,xTestFunc,2
                zAnsFuncMin = zTestFunc-2
                if( zAnsFuncMin < 1 ) then
                  zAnsFuncMin = zTestFunc
                end if
                do zAnsFunc = zAnsFuncMin,zTestFunc,2

                  ! the position of the modal coefficeint of this ansatz functions
  anspos = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (zansfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)

                  ! calculate the projection of the ansatz and test function
                  xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
                  zScalProd = ply_scalProdDualLeg(zAnsFunc, zTestFunc)

!NA!                  ! buffer the result in kernel data and take care of the outer surface unit normal
!NA!                  ! ... for the element left of the face
!NA!                  outerNormalLeft = 1.0_rk
!NA!                  projectedFlux = outerNormalLeft * numFluxBuffer(ansPos,:)
!NA!                  kernelData%state_der(left_neighbor, testPos, 1, 2:4) &
!NA!                            & = kernelData%state_der(left_neighbor, testPos, 1, 2:4) &
!NA!                            & - faceValLeft * xScalProd * zScalProd &
!NA!                            & * projectedFlux * jacobiDet
!NA!                  ! ... for the element right of the face
!NA!                  outerNormalRight = -1.0_rk
!NA!                  projectedFlux = outerNormalRight * numFluxBuffer(ansPos,:)
!NA!                  kernelData%state_der(right_neighbor, testPos,1,2:4) &
!NA!                            & = kernelData%state_der(right_neighbor, testPos, 1, 2:4) &
!NA!                            & - faceValRight * xScalProd * zScalProd &
!NA!                            & * projectedFlux * jacobiDet

                end do
              end do
            end do

          end do
        end do
      end do

    end do YFaceLoop


  end subroutine modg_nerplanck_y_numFlux_pre

  !> Numerical flux calculation for Nernst-Planck equation across the faces in Z direction.
  subroutine modg_nerplanck_z_numFlux_pre( mesh, equation, statedata, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The equation you solve.
    type(atl_equations_type) :: equation
    !> THe state if the equation
    type(atl_statedata_type) :: statedata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------!
    ! Loop vars
    integer :: xAnsFunc, yAnsFunc, zAnsFunc, xTestFunc, yTestFunc, zTestFunc
    ! Positions for the modal coefficients
    integer :: modalCoeff, testPos, ansPos
    ! Modal coefficients for elements left and right of the face
    real(kind=rk) :: leftModalCoeff(4), rightModalCoeff(4), flux(3)
    ! Modal representation of the flux on the face
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    ! The values of the x ansatz and test functions at the face from left and right
    ! limit.
    real(kind=rk) :: faceValLeft, faceValRight
    ! Loop var for the faces.
    integer :: zside
    ! Element positions of the left and right element of the face.
    integer(kind=long_k) :: left_neighbor, right_neighbor
    ! Speed of light, i.e. maximum information propagation speed.
    real(kind=rk) :: diffusivity
    ! Scalar products between x and z ansatz and test functions.
    real(kind=rk) :: xScalProd, yScalProd
    ! Determinant of the jacobian of the mapping from reference face to
    ! the current face.
    real(kind=rk) :: jacobiDet
    ! The minimum indices for the relevant ansatz functions.
    integer :: xAnsFuncMin, yAnsFuncMin
    ! Rotation for the flux
    integer :: varRotation(4)
    ! --------------------------------------------------------------------------

    allocate(numFluxBuffer((scheme%maxPolyDegree+1)**3,3))

    ! Calculate the speed of light
    diffusivity = equation%nerplanck%D

    ! Build the determinant of the Jacobian of the mapping from reference
    ! face to the physical face. Since we are considering cubes it is
    ! the same everywhere.
    jacobiDet = (mesh%length/2.0_rk)**2.0_rk

    ! Loop over all fluid faces in z direction
    varRotation = equation%varRotation(3)%varTransformIndices(:)
    zFaceLoop: do zside = 1,size(mesh%faces%faces(3)%computeFace%leftPos)
      ! Get the fluid neighbors for this face.
      left_neighbor = mesh%faces%faces(3)%computeFace%leftPos(zside)
      right_neighbor = mesh%faces%faces(3)%computeFace%rightPos(zside)


      ! Get the neighbors for this face.
      !left_neighbor  = mesh%side%zdir%left_elemid(zside)
      !right_neighbor = mesh%side%zdir%right_elemid(zside)

      ! Loop over all the ansatz functions. Here, we use simple tensor products
      ! of one-dimensional ansatz functions.
      do zAnsFunc = 1, scheme%maxPolyDegree+1
        faceValRight = ply_faceValLeftBndAns(zAnsFunc)
        do xAnsFunc = 1, (scheme%maxPolyDegree+1)
          do yAnsFunc = 1, (scheme%maxPolyDegree+1)
            ! Get the modal coefficients of the left and right ansatz functions and
            ! calculate the modal flux with them.
  modalcoeff = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (zansfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
            ! ... keep in mind that we have to rotate the input and output of the flux
            leftModalCoeff(:) = statedata%state(left_neighbor, &
                                               & modalCoeff, varRotation)
            rightModalCoeff(:) = statedata%state(right_neighbor, &
                                               &  modalCoeff, varRotation) * faceValRight
            call atl_nerplanck_numflux_preprocess( &
              & left        = leftModalCoeff,      &
              & right       = rightModalCoeff,     &
              & diffusivity = diffusivity,         &
              & flux        = flux                 )

            ! Buffer the modal values of the numerical flux
            numFluxBuffer( modalCoeff, :) = flux
          end do
        end do
      end do

      ! Loop over all the test functions and project the numerical flux to them.
      do zTestFunc = 1,min(scheme%maxPolyDegree+1,2)
        faceValLeft = ply_faceValRightBndTest(zTestFunc)
        faceValRight = ply_faceValLeftBndTest(zTestFunc)
        do xTestFunc = 1,scheme%maxPolyDegree+1
          do yTestFunc = 1,scheme%maxPolyDegree+1

            ! position of the test functions in the kerneldata
  testpos = xtestfunc                                      &
    &      + ( ( ytestfunc-1)                             &
    &      + (ztestfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)

            ! get the relevant indices of ansatz functions for the projection
            do zAnsFunc = 1, scheme%maxPolyDegree+1
              xAnsFuncMin = xTestFunc-2
              if( xAnsFuncMin < 1 ) then
                xAnsFuncMin = xTestFunc
              end if
              do xAnsFunc = xAnsFuncMin,xTestFunc,2
                yAnsFuncMin = yTestFunc-2
                if( yAnsFuncMin < 1 ) then
                  yAnsFuncMin = yTestFunc
                end if
                do yAnsFunc = yAnsFuncMin,yTestFunc,2

                  ! the position of the modal coefficeint of this ansatz functions
  anspos = xansfunc                                      &
    &      + ( ( yansfunc-1)                             &
    &      + (zansfunc-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)

                  ! calculate the projection of the ansatz and test function
                  xScalProd = ply_scalProdDualLeg(xAnsFunc, xTestFunc)
                  yScalProd = ply_scalProdDualLeg(yAnsFunc, yTestFunc)

!NA!                  ! buffer the result in kernel data and take care of the outer surface unit normal
!NA!                  ! ... for the element left of the face
!NA!                  outerNormalLeft = 1.0_rk
!NA!                  projectedFlux = outerNormalLeft * numFluxBuffer(ansPos,:)
!NA!                  kernelData%state_der(left_neighbor, testPos, 1, 2:4) &
!NA!                                  & = kernelData%state_der(left_neighbor, testPos, 1, 2:4) &
!NA!                                  & - faceValLeft * xScalProd * yScalProd &
!NA!                                  & * projectedFlux * jacobiDet
!NA!                  ! ... for the element right of the face
!NA!                  outerNormalRight = -1.0_rk
!NA!                  projectedFlux = outerNormalRight * numFluxBuffer(ansPos,:)
!NA!                  kernelData%state_der(right_neighbor, testPos,1, 2:4) &
!NA!                                  & = kernelData%state_der(right_neighbor, testPos, 1, 2:4) &
!NA!                                  & - faceValRight * xScalProd * yScalProd &
!NA!                                  & * projectedFlux * jacobiDet

                end do
              end do
            end do

          end do
        end do
      end do

    end do ZFaceLoop
  end subroutine modg_nerplanck_z_numFlux_pre


  !> Calculate the projection of the physical flux for the MODG scheme and
  !! Nernst-Planck equation.
  subroutine modg_nerplanck_physFlux_pre( mesh, equation, statedata, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The equation you solve.
    type(atl_equations_type) :: equation
    !> THe state if the equation
    type(atl_statedata_type) :: statedata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------!
    ! Loop var for all the fluid elements
    integer :: iElem
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars for ansatz and test functions
    integer :: iAnsX, iAnsY, iAnsZ, iTestX, iTestY, iTestZ
    ! Position of ansatz and test functions in the linearized list of coefficients
    integer :: ansPos, testPos
    ! Determinant of the jacobian of the mapping of reference element to physical
    ! element.
    real(kind=rk) :: jacobiDet
    ! Scalar products for the volume integral
    real(kind=rk) :: scalProdX, scalProdY, scalProdZ
    ! Modal represenation of the physical flux in x, y and z direction
    real(kind=rk), allocatable :: physFlux(:,:)
    ! Buffer for the flux
    real(kind=rk) :: flux
    ! Material parameters of the Nernst-Planck equations
    real(kind=rk) :: diffusivity
    ! Lower loop bounds for the projections onto the test functions
    integer :: iAnsXMin, iAnsYMin, iAnsZMin
    ! --------------------------------------------------------------------------!

    allocate( modalCoeffs((scheme%maxPolyDegree+1)**3,4) )
    allocate( physFlux((scheme%maxPolyDegree+1)**3,1) )

    ! get the Nernst-Planck material parameters
    diffusivity = equation%nerplanck%D

    ! We have cubiic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDet = (mesh%length/2.0_rk)**2

    ! Loop over all the elements and calculate the projection
    elemLoop: do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

      ! get the modal coefficients of the current cell (for all variables
      ! of the equation, therefore we use ":" for the third index).
      modalCoeffs = statedata%state(iElem,:,:)

      ! calculate the physical flux in modal representation, so we loop over all
      ! of them
      do iAnsX = 1, scheme%maxPolyDegree+1
        do iAnsY = 1, scheme%maxPolyDegree+1
          do iAnsZ = 1, scheme%maxPolyDegree+1
            ! position of the modal coefficient
  anspos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
            call atl_nerplanck_physflux_preprocess(  &
              & state       = modalCoeffs(ansPos,:), &
              & diffusivity = diffusivity,           &
              & flux        = flux                   )
             physFlux(ansPos,:) = flux
          end do
        end do
      end do

      ! now, project onto all test functions
      do iTestX = 1, scheme%maxPolyDegree+1
        do iTestY = 1, scheme%maxPolyDegree+1
          do iTestZ = 1, scheme%maxPolyDegree+1
            ! the position of the current test functions
  testpos = itestx                                      &
    &      + ( ( itesty-1)                             &
    &      + (itestz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)

            ! for x direction (x test function differentiated)
            ! get the relevant indices for the ansatz function
            if(iTestX > 1) then
              iAnsX = iTestX - 1
              scalProdX = ply_scalProdDualLegDiff(iAnsX, iTestX)
              iAnsYMin = iTestY-2
              if(iAnsYMin < 1) then
                iAnsYMin = iTestY
              end if
              do iAnsY = iAnsYMin, iTestY, 2
                scalProdY = ply_scalProdDualLeg(iAnsY, iTestY)
                iAnsZMin = iTestZ-2
                if(iAnsZMin < 1) then
                  iAnsZMin = iTestZ
                end if
                do iAnsZ = iAnsZMin, iTestZ, 2
                  scalProdZ = ply_scalProdDualLeg(iAnsZ, iTestZ)
                  ! now, we calculate the projection
  anspos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA!                  kernelData%state_der(iElem, testPos,1,1) &
!NA!                                  & = kernelData%state_der(iElem, testPos, 1, 1) &
!NA!                                  & + physFlux(ansPos,1) * jacobiDet * scalProdX &
!NA!                                  & * scalProdY * scalProdZ
                end do
              end do
            end if


            !  for y direction (y test function differentiated)
            if(iTestY > 1) then
              iAnsY = iTestY - 1
              scalProdY = ply_scalProdDualLegDiff(iAnsY, iTestY)
              iAnsXMin = iTestX-2
              if(iAnsXMin < 1) then
                iAnsXMin = iTestX
              end if
              do iAnsX = iAnsXMin, iTestX, 2
                scalProdX = ply_scalProdDualLeg(iAnsX, iTestX)
                iAnsZMin = iTestZ-2
                if(iAnsZMin < 1) then
                  iAnsZMin = iTestZ
                end if
                do iAnsZ = iAnsZMin, iTestZ, 2
                  scalProdZ = ply_scalProdDualLeg(iAnsZ, iTestZ)
                  ! now, we calculate the projection
  anspos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA!                  kernelData%state_der(iElem, testPos,1,1) &
!NA!                                  & = kernelData%state_der(iElem, testPos, 1, 1) &
!NA!                                  & + physFlux(ansPos,1) * jacobiDet * scalProdX &
!NA!                                  & * scalProdY * scalProdZ
                end do
              end do
            end if

            ! for z direction (z test function differentiated)
            if(iTestZ > 1) then
              iAnsZ = iTestZ - 1
              scalProdZ = ply_scalProdDualLegDiff(iAnsZ, iTestZ)
              iAnsYMin = iTestY-2
              if(iAnsYMin < 1) then
                iAnsYMin = iTestY
              end if
              do iAnsY = iAnsYMin, iTestY, 2
                scalProdY = ply_scalProdDualLeg(iAnsY, iTestY)
                iAnsXMin = iTestX-2
                if(iAnsXMin < 1) then
                  iAnsXMin = iTestX
                end if
                do iAnsX = iAnsXMin, iTestX, 2
                  scalProdX = ply_scalProdDualLeg(iAnsX, iTestX)
                  ! now, we calculate the projection
  anspos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA!                  kernelData%state_der(iElem, testPos,1,1) &
!NA!                                  & = kernelData%state_der(iElem, testPos, 1, 1) &
!NA!                                  & + physFlux(ansPos,1) * jacobiDet * scalProdX &
!NA!                                  & * scalProdY * scalProdZ
                end do
              end do
            end if

          end do
        end do
      end do

    end do elemLoop

  end subroutine modg_nerplanck_physFlux_pre


  !> Applies the inverse of the mass matrix for a 3D scheme.
  subroutine modg_invMassMatrix_pre( mesh, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------!
    ! Loop indices for ansatz functions
    integer :: iAnsX, iAnsY, iAnsZ
    ! Positions for the given ansatz functions
    integer :: ansLow, ansPrevLow, ansUp, ansPrevUp
    ! Loop var for the elements
    integer :: iElem
    ! Inverse of the determintant of the jacobian of the mapping from reference
    ! to physical element.
    real(kind=rk) :: inv_jacobiDet
    ! --------------------------------------------------------------------------!

    inv_jacobiDet = (2.0_rk/mesh%length)**3.0_rk

    ! Loop over all the elements and apply the inverse mass matrix
    elemLoop: do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

      ! apply the 1D inverse of the mass matrix
      do iAnsX = 3, scheme%maxPolyDegree+1
        do iAnsY = 1, scheme%maxPolyDegree+1
          do iAnsZ = 1, scheme%maxPolyDegree+1
  anslow = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansprevlow = iansx-2                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA!            kerneldata%state_der(iElem, ansLow,1,:) = kerneldata%state_der(iElem, ansLow,1,:) &
!NA!                                                    & + kerneldata%state_der(iElem, ansPrevLow,1,:)
          end do
        end do
      end do
      do iAnsX = 1, scheme%maxPolyDegree+1
        do iAnsY = 1, scheme%maxPolyDegree+1
          do iAnsZ = 1, scheme%maxPolyDegree+1
  anslow = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA!            kerneldata%state_der(iElem, ansLow,1,:) = kerneldata%state_der(iElem, ansLow,1,:) &
!NA!                                                    & * (2.0_rk*iAnsX-1.0_rk)/2.0_rk
          end do
        end do
      end do

      ! apply the 2D inverse of the mass matrix
      do iAnsY = 3, scheme%maxPolyDegree+1
        do iAnsZ = 1, scheme%maxPolyDegree+1
  anslow = 1                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansup = scheme%maxpolydegree+1                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansprevlow = 1                                      &
    &      + ( ( iansy-2-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansprevup = scheme%maxpolydegree+1                                      &
    &      + ( ( iansy-2-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA!          kerneldata%state_der(iElem, ansLow:ansUp,1,:) &
!NA!                                   & = kerneldata%state_der(iElem, ansLow:ansUp,1,:) &
!NA!                                   & + kerneldata%state_der(iElem, ansPrevLow:ansPrevUp,1,:)
        end do
      end do
      do iAnsY = 1, scheme%maxPolyDegree+1
        do iAnsZ = 1, scheme%maxPolyDegree+1
  anslow = 1                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansup = scheme%maxpolydegree+1                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA          kerneldata%state_der(iElem, ansLow:ansUp,1,:) &
!NA                                   & = kerneldata%state_der(iElem, ansLow:ansUp,1,:) &
!NA                                   & * (2.0_rk*iAnsY-1.0_rk)/2.0_rk
        end do
      end do

      ! apply the 3D inverse of the mass matrix
      do iAnsZ = 3, scheme%maxPolyDegree+1
  anslow = 1                                      &
    &      + ( ( 1-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansup = scheme%maxpolydegree+1                                      &
    &      + ( ( scheme%maxpolydegree+1-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansprevlow = 1                                      &
    &      + ( ( 1-1)                             &
    &      + (iansz-2-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansprevup = scheme%maxpolydegree+1                                      &
    &      + ( ( scheme%maxpolydegree+1-1)                             &
    &      + (iansz-2-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA!        kerneldata%state_der(iElem, ansLow:ansUp,1,:) &
!NA!                                   & = kerneldata%state_der(iElem, ansLow:ansUp,1,:) &
!NA!                                   & + kerneldata%state_der(iElem, ansPrevLow:ansPrevUp,1,:)
      end do
      do iAnsZ = 1, scheme%maxPolyDegree+1
  anslow = 1                                      &
    &      + ( ( 1-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
  ansup = scheme%maxpolydegree+1                                      &
    &      + ( ( scheme%maxpolydegree+1-1)                             &
    &      + (iansz-1)*(scheme%maxpolydegree+1))*(scheme%maxpolydegree+1)
!NA!        kerneldata%state_der(iElem, ansLow:ansUp,1,:) &
!NA!                                   & = kerneldata%state_der(iElem, ansLow:ansUp,1,:) &
!NA!                                   & * (2.0_rk*iAnsZ-1.0_rk)/2.0_rk
      end do
    end do elemLoop

  end subroutine modg_invMassMatrix_pre
end module atl_modg_nerplanck_kernel_module
