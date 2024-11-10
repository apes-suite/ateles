! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2014, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!! scheme for the Heat equation.

module atl_modg_1d_heat_kernel_module
  use env_module,                only: rk
  use atl_equation_module,       only: atl_equations_type
  use atl_cube_elem_module,      only: atl_cube_elem_type
  use atl_modg_1d_scheme_module, only: atl_modg_1d_scheme_type
  use ply_poly_project_module,   only: ply_poly_project_type, assignment(=)
  use atl_facedata_module,       only: atl_facedata_type, atl_elemfaceToNormal_prp
  use tem_faceData_module,       only: tem_left, tem_right,    &
    &                                  tem_faceIterator_type
  use ply_leg_diff_module,       only: ply_calcDiff_leg_1d

  implicit none

  private

  public :: atl_modg_1d_heat_numflux, atl_modg_1d_heat_physFlux


contains


  !> Calculate the physical flux for the MODG scheme and
  !! Heat equation.
  subroutine atl_modg_1d_heat_physFlux( mesh, equation, state, res, &
    &                                   poly_proj, modalCoeffs,     &
    &                                   modalCoeffsDiff             )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type) :: mesh
    !> The equation you solve.
    type(atl_equations_type) :: equation
    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj
    real(kind=rk), intent(in) :: state(:,:)
    real(kind=rk), intent(inout) :: res(:,:)
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)
    real(kind=rk), intent(inout) :: modalCoeffsDiff(:,:)

    ! --------------------------------------------------------------------------!
    integer :: dof, dofs
    real(kind=rk) ::  therm_diff
    ! ---------------------------------------------------------------------------

    dofs = poly_proj%body_1d%ndofs


    therm_diff = equation%heat%k

    ! get the modal coefficients of the current cell (for all variables)
    ! -->modal space

    do dof=lbound(modalCoeffs,2),ubound(modalCoeffs,2)
      modalCoeffs(:,dof) = state(:,dof)
    end do

    ! Now we calculate the gradient of the modal representation required
    ! for calculating the Physical Fluxes
    call ply_calcDiff_leg_1d( legCoeffs     = modalCoeffs,             &
      &                       legcoeffsDiff = modalCoeffsDiff,         &
      &                       maxPolyDegree = poly_proj%maxPolyDegree, &
      &                       eLemLength    = mesh%length              )

    do dof=lbound(res,1),ubound(res,1)
      res(dof,1) = -therm_diff*modalCoeffsDiff(dof,1)
    end do
  end subroutine atl_modg_1d_heat_physFlux

  !> Calculate the numerical flux for the Heat equation and MODG scheme
  subroutine atl_modg_1d_heat_numFlux( mesh, equation, facedata, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_1d_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------
    integer :: iDir
    ! --------------------------------------------------------------------------
    ! Numerical flux for faces in all 1 spatial face directions (x dir)

    !>TODO add if statement -> nfaces/thread_num

    do iDir = 1,1
      call modg_1d_heat_oneDim_numFlux(                                        &
        &                       equation = equation ,                          &
        &                       facedata = facedata,                           &
        &                         scheme = scheme ,                            &
        &                          faces = mesh%faces%faces(iDir)%computeFace, &
        &                        faceDir = iDir,                               &
        &                      elem_len  = mesh%length                         )
    end do

  end subroutine atl_modg_1d_heat_numFlux




  !> Numerical flux calculation for the Heat equation across the faces in a
  !! single spatial direction.
  subroutine modg_1d_heat_oneDim_numFlux( equation, facedata, scheme, faces, &
                                        & faceDir, elem_len  )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face state if the equation
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_1d_scheme_type), intent(in) :: scheme
    !> The faces to calculate the fluxes for.
    type(tem_faceIterator_type), intent(in) :: faces
    !> The spatial direction of the faces you calc the fluxes for,
    !! use the following:
    !! 1 --> x direction.
    integer, intent(in) :: faceDir
    !> Length of the element
    real(kind=rk), intent(in) :: elem_len
    ! --------------------------------------------------------------------------!
    ! Loop vars
    ! Modal coefficients for elements left and right of the face:
    ! First dimension is the number of modal coefficients on the face, second
    ! is the number of variables.
    real(kind=rk) :: leftModalCoeff(1,equation%varSys%nScalars ), &
                   & rightModalCoeff(1,equation%varSys%nScalars)
    real(kind=rk) :: leftModalCoeffDiff(1,equation%varSys%nScalars ), &
                   & rightModalCoeffDiff(1,equation%varSys%nScalars)
    real(kind=rk) :: Numflux1(1,equation%varSys%nScalars )
    real(kind=rk) :: Numflux2(1,equation%varSys%nScalars )
    ! Loop var for the faces.
    integer :: iside
    ! Element positions of the left and right element of the face.
    integer :: left_neighbor, right_neighbor
    ! Loop over variables (due to single variable FPTs)
    integer :: iVar
    ! Outer Normals
    real(kind=rk) :: outerNormalLeft, outerNormalRight, Const, C_IP
    real(kind=rk) ::  therm_diff
    ! --------------------------------------------------------------------------
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)
    therm_diff = equation%heat%k
    C_IP = 20.0
    ! for the term \int \delta \jump{u_h} \cdot \jump{v} ds
    Const = C_IP*(scheme%maxPolyDegree**2)/ elem_len

    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, size(faces%leftPos)

      ! Get the fluid neighbors for this face.
      left_neighbor = faces%leftPos(iside)
      right_neighbor = faces%rightPos(iside)

      ! for the left element, we have to access the right face values
      ! and for the right, we have to acess the left face values.
      leftModalCoeffDiff(:,:)  = 0.0_rk
      rightModalCoeffDiff(:,:) = 0.0_rk
      leftModalCoeff(:,:)      = 0.0_rk
      rightModalCoeff(:,:)     = 0.0_rk

      ! get the modal representation of the state
      leftModalCoeff(:,1)      = facedata%faceRep(faceDir)     &
        &                             %dat(left_neighbor,:,1,2)
      rightModalCoeff(:,1)     = facedata%faceRep(faceDir)     &
        &                            %dat(right_neighbor,:,1,1)

      ! get modal representation of the derivative of the state
      leftModalCoeffDiff(:,1)  = facedata%faceRep(faceDir)     &
        &            %dat(left_neighbor,:,1+faceDir*equation%varSys%nScalars,2)

      rightModalCoeffDiff(:,1) = facedata%faceRep(faceDir)     &
       &             %dat(right_neighbor,:,1+faceDir*equation%varSys%nScalars,1)

      !term1!
      ! for the term \int \average{\nabla u_h} \cdot \jump{v} ds
      ! where average and jump are respective operators and v is test func.
      do iVar = 1, equation%varSys%nScalars
        numFlux1(:,iVar) =   LeftmodalCoeffDiff(:,iVar)   +          &
          &                 RightmodalCoeffDiff(:,iVar)
        facedata%faceFlux(faceDir)%dat(left_neighbor,:,iVar, 2) =      &
          &                             -numFlux1(:,iVar)*0.5*therm_diff
        facedata%faceFlux(faceDir)%dat(right_neighbor,:,iVar, 1) =     &
          &                             -numFlux1(:,iVar)*0.5*therm_diff
      end do

      !term 2 & 3!
      !const=1
      do iVar = 1, equation%varSys%nScalars

        numFlux2(:,iVar) =  LeftmodalCoeff(:,iVar)* outerNormalRight   +     &
          &                 RightmodalCoeff(:,iVar)* outerNormalLeft

        facedata%faceFlux(faceDir)%dat(left_neighbor,:,iVar, 2)  =     &
          &                            const*numFlux2(:,iVar)*therm_diff
        facedata%faceFlux(faceDir)%dat(right_neighbor,:,iVar, 1) =     &
          &                            const*numFlux2(:,iVar)*therm_diff


        facedata%faceFlux(faceDir)%dat( left_neighbor,:,               &
          & equation%varSys%nScalars + iVar, 2) =                           &
          &                        -numFlux2(:,iVar)*0.5*therm_diff
        facedata%faceFlux(faceDir)%dat( right_neighbor,:,              &
          & equation%varSys%nScalars + iVar, 1) =                           &
          &                        -numFlux2(:,iVar)*0.5*therm_diff

      end do

    end do FaceLoop

  end subroutine modg_1d_heat_oneDim_numFlux

!! ******************************************************************************!

end module atl_modg_1d_heat_kernel_module
