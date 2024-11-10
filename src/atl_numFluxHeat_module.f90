! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2014, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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
module atl_numFluxHeat_module
  use env_module,               only: rk

  use tem_faceData_module,      only: tem_left,  &
    &                                 tem_right, &
    &                                 tem_faceIterator_type

  use atl_equation_module,      only: atl_equations_type
  use atl_facedata_module,      only: atl_facedata_type, &
    &                                 atl_elemfaceToNormal_prp

  use ply_poly_project_module,  only: assignment(=)


  implicit none

  private

  public :: atl_modg_heat_numFlux_sipg

contains

  !> Numerical flux calculation for Heat equation across the faces in a single
  !! spatial direction.
  subroutine atl_modg_heat_numFlux_sipg( equation, facedata, faces, faceDir, &
    &                                    dofs, elem_len, maxPolyDeg          )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face state if the equation
    type(atl_facedata_type), intent(inout) :: facedata
    !> The faces to calculate the fluxes for.
    type(tem_faceIterator_type), intent(in) :: faces
    !> The spatial direction of the faces you calc the fluxes for, use the following:
    !! 1 --> x direction. \n
    !! 2 --> y direction. \n
    !! 3 --> z direction.
    integer, intent(in) :: faceDir
    !> Parameter for used projection
    integer, intent(in) :: dofs
    !> Length of the element
    real(kind=rk), intent(in) :: elem_len
    !> Max polynomial degree
    integer, intent(in) :: maxPolyDeg
    ! --------------------------------------------------------------------------!
    ! Modal coefficients for elements left and right of the face:
    ! First dimension is the number of modal coefficients on the face, second
    ! is the number of variables.
    real(kind=rk), allocatable :: leftModalCoeffs(:,:), &
                               &  rightModalCoeffs(:,:)
    real(kind=rk), allocatable  :: leftModalCoeffsDiff(:,: ), &
                   & rightModalCoeffsDiff(:,:)
    ! Loop var for the faces.
    integer :: iside, iVar
    ! Element positions of the left and right element of the face.
    integer :: left_neighbor, right_neighbor
    ! Loop over variables (due to single variable FPTs)
    !integer :: nquadpoints, oversamp_dofs
    real(kind=rk) :: outerNormalLeft, outerNormalRight, Const, C_IP
    real(kind=rk) ::  therm_diff
    real(kind=rk) , allocatable :: Numflux1(:,:)
    real(kind=rk) , allocatable :: Numflux2(:,:)
    ! --------------------------------------------------------------------------
    outerNormalLeft = atl_elemfaceToNormal_prp(tem_left)
    outerNormalRight = atl_elemfaceToNormal_prp(tem_right)
    therm_diff = equation%heat%k
    C_IP = 20.0

    allocate( leftModalCoeffs(dofs,equation%varSys%nScalars       ))
    allocate( rightModalCoeffs(dofs,equation%varSys%nScalars      ))
    allocate( leftModalCoeffsDiff(dofs,equation%varSys%nScalars   ))
    allocate( rightModalCoeffsDiff(dofs,equation%varSys%nScalars  ))
    allocate( Numflux1(dofs,equation%varSys%nScalars              ))
    allocate( Numflux2(dofs,equation%varSys%nScalars              ))


    ! Loop over all fluid the faces in x direction
    FaceLoop: do iside = 1, size(faces%leftPos)

      ! Get the fluid neighbors for this face.

      left_neighbor = faces%leftPos(iside)
      right_neighbor = faces%rightPos(iside)

      leftModalCoeffs(:,:)       = 0.0_rk
      rightModalCoeffs(:,:)      = 0.0_rk
      leftModalCoeffsDiff(:,:)   = 0.0_rk
      rightModalCoeffsDiff(:,:)  = 0.0_rk


      leftModalCoeffs(:,1) =         &
        &           facedata%faceRep(faceDir)%dat(left_neighbor,:,1,2)
      RightModalCoeffs(:,1) =        &
        &           facedata%faceRep(faceDir)%dat(right_neighbor,:,1,1)

      leftModalCoeffsDiff(:,1) =               &
         facedata%faceRep(faceDir)%dat(left_neighbor,:,1+faceDir*equation%varSys%nScalars,2)
      RightModalCoeffsDiff(:,1) =              &
         facedata%faceRep(faceDir)%dat(right_neighbor,:,1+faceDir*equation%varSys%nScalars,1)


      do iVar = 1, equation%varSys%nScalars
        numFlux1(:,iVar) =   LeftmodalCoeffsDiff(:,1)   +          &
          &                 RightmodalCoeffsDiff(:,1)

        facedata%faceFlux(faceDir)%dat(left_neighbor,:,iVar, 2) =      &
          &                             -numFlux1(:,iVar)*0.5*therm_diff
        facedata%faceFlux(faceDir)%dat(right_neighbor,:,iVar, 1) =     &
          &                             -numFlux1(:,iVar)*0.5*therm_diff
      end do

      !term 2 & 3!
      ! for the term \int \delta \jump{u_h} \cdot \jump{v} ds
      Const = C_IP*(maxPolyDeg**2)/ elem_len

     ! const = 1.0
      do iVar = 1, equation%varSys%nScalars

        numFlux2(:,iVar) =  LeftmodalCoeffs(:,iVar)* outerNormalRight   +     &
          &                 RightmodalCoeffs(:,iVar)* outerNormalLeft


        facedata%faceFlux(faceDir)%dat(left_neighbor,:,iVar, 2)  =     &
          &                            const*numFlux2(:,iVar)*therm_diff
        facedata%faceFlux(faceDir)%dat(right_neighbor,:,iVar, 1) =     &
          &                            const*numFlux2(:,iVar)*therm_diff


         facedata%faceFlux(faceDir)%dat(left_neighbor,:,  &
          &  iVar+equation%varSys%nScalars, 2) =          &
          &                        -numFlux2(:,iVar)*0.5*therm_diff
         facedata%faceFlux(faceDir)%dat(right_neighbor,:, &
          &  iVar+equation%varSys%nScalars, 1) =          &
          &                        -numFlux2(:,iVar)*0.5*therm_diff

      end do
    end do FaceLoop
  end subroutine atl_modg_heat_numFlux_sipg

end module atl_numFluxHeat_module
