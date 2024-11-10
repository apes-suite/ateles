! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> module that holds all routines to calculate the flux for
!! hyperbolic Maxwell equations (2d, transverse electric mode formulation -TE).
module atl_maxwell_flux_2d_module
  ! Treelm modules
  use env_module, only: rk

  ! Ateles modules
  use ply_poly_project_module,  only: ply_poly_project_type, &
    &                                 assignment(=),         &
    &                                 ply_poly_project_m2n,  &
    &                                 ply_poly_project_n2m
  implicit none

  private

  public :: atl_maxwell_flux_2d

  !> Interface for fluxes of pure Maxwell equations.
  interface atl_maxwell_flux_2d
    module procedure maxwell_flux_cube_2d
    module procedure maxwell_flux_cube_vec_2d
    module procedure maxwell_flux_nonconst_cube_vec_2d
  end interface atl_maxwell_flux_2d


contains

  !> summary: Subroutine to calculate the flux for pure Maxwell equations without
  !! any divergence cleaning on the reference cubic face in 2D.
  !!
  !! This subroutine calculates the flux of the Maxwell equation on the reference
  !! cubic face. The underlying 2D formulation is transverse electric mode formulation
  !! - TE mode.
  subroutine maxwell_flux_cube_2d(left, right, left_mu, left_epsi, right_mu, right_epsi, flux)
    ! --------------------------------------------------------------------------
    !> Left state vector (as conservative variables). The order of this vector
    !! has to be \f$ (D_x, D_y, B_3) \f$ where E and B denoted
    !! electric field vetor and magnetic field (also called magnetic induction) vector.
    real(kind=rk), intent(in)  :: left(7)
    !> Right state vector (as conservative variables). The order of this vector
    !! has to be (D_x, D_y, B_3) where E and B denoted the electric
    !! field vetor and magnetic field (also called magnetic induction) vector.
    real(kind=rk), intent(in)  :: right(7)
    !>  The magnetic permeability of the left element.
    real(kind=rk), intent(in)  :: left_mu
    !>  The electric permitivity of the left element.
    real(kind=rk), intent(in)  :: left_epsi
    !>  The magnetic permeability of the right element.
    real(kind=rk), intent(in)  :: right_mu
    !>  The electric permitivity of the right element.
    real(kind=rk), intent(in)  :: right_epsi
    !> The flux between left and right cell. The order of this vector is the
    !! same as the input arguments.
    real(kind=rk), intent(out) :: flux(7)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: left_speedOfLight, right_speedOfLight
    real(kind=rk) :: inv_denom_mu, inv_denom_epsi
    ! --------------------------------------------------------------------------

    ! The speed of light in the left and right element
    left_speedOfLight = 1.0_rk / sqrt( left_mu * left_epsi )
    right_speedOfLight = 1.0_rk / sqrt( right_mu * right_epsi )

    ! The inverse of the denominators
    inv_denom_mu = 1.0_rk / ((-1.0_rk)*left_speedOfLight*left_mu - right_speedOfLight*right_mu)
    inv_denom_epsi = 1.0_rk / (left_speedOfLight*left_epsi + right_speedOfLight*right_epsi)

    ! D_x
    flux(1) = 0.0_rk

    ! PML
    flux(4:7) = 0.0_rk

    ! the flux for D_y
    flux(2) = ( &
                &    ( (-1.0_rk*left(2) / left_epsi)        &
                &      - (-1.0_rk*right(2) / right_epsi)  ) &
                & -  ( left_speedOfLight * left(3)          &
                &      + right_speedOfLight * right(3) ))
    ! the flux for B_z
    flux(3) = ( &
                &   ( left_speedOfLight * left(2)           &
                &   + right_speedOfLight * right(2) )       &
                & + ( ( left(3) / left_mu )                 &
                &   - ( right(3) / right_mu ) ))

     ! Normalize the calculated fluxes for D_y and B_z
     flux(2) = inv_denom_mu * flux(2)
     flux(3) = inv_denom_epsi * flux(3)

  end subroutine maxwell_flux_cube_2d

  !> summary: calculate flux of pure maxwell equation directly on the face-vector
  !! (formulation is based on TE mode formulation for 2D).
  !!
  !! This subroutine assumes the Maxwell equations with D and B as input
  !! variables. Furthermore, it is able to handle arbitray material parameters
  !! by a pseudo-spectral technique.
  subroutine maxwell_flux_nonconst_cube_vec_2d(nTotalFaces, nSides, nFaceDofs, &
    &                                       faceRep, faceFlux,              &
    &                                       leftPos, rightPos, var,         &
    &                                       material_left, material_right,  &
    &                                       poly_proj, modalCoeffs, pntVal, &
    &                                       nodalNumFlux, numFluxBuffer     )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: nTotalFaces, nFaceDofs, nSides
    !> The modal representation on the faces, left and right trace.
    real(kind=rk), intent(in) :: faceRep(nTotalFaces,nFaceDofs,7,2)
    !> The fluxes for all faces, for left and right elements.
    real(kind=rk), intent(inout) :: faceFlux(nTotalFaces,nFaceDofs,7,2)
    !> Positions for the left and right elements of all faces
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    !> Variable rotation indices
    integer, intent(in) :: var(7)
    !> Material parameters for the left faces.
    real(kind=rk), intent(in) :: material_left(nSides,nFaceDofs,3)
    !> Material parameters for the right faces.
    real(kind=rk), intent(in) :: material_right(nSides,nFaceDofs,3)
    !!> FPT to convert between nodes and modes.
    !type(atl_legFpt_2D_type), intent(inout) :: fpt
    !> Data for projection method
    type(ply_poly_project_type) :: poly_proj
    !> Working array for the left and right modal coefficients
    real(kind=rk), intent(inout) :: modalCoeffs(:,:,:)
    !> Working array for the left and right point values
    real(kind=rk), intent(inout) :: pntVal(:,:,:)
    !> Working array for the nodal flux
    real(kind=rk), intent(inout) :: nodalNumFlux(:,:)
    !> Working array for the modal numerical flux
    real(kind=rk),  intent(inout) :: numFluxBuffer(:,:)
    ! --------------------------------------------------------------------------
    integer :: iSide, left, right, iPoint
    real(kind=rk) :: left_mu, right_mu
    real(kind=rk) :: left_epsi, right_epsi
    real(kind=rk) :: flux(7)
    integer :: nquadpoints, ndofs, maxPolyDegree, oversamp_degree, oversamp_dofs
    integer :: iDegX
    ! --------------------------------------------------------------------------
    nquadpoints = poly_proj%body_1d%nquadpoints
    ndofs = poly_proj%body_2D%ndofs
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs
    oversamp_degree = poly_proj%oversamp_degree
    maxpolyDegree = poly_proj%maxPolyDegree

    ! Loop over all the faces and calculate the non-zero fluxes
    faceLoop: do iSide = 1, nSides
      ! Get indices of left and right adjacent element
      left = leftPos(iSide)
      right = rightPos(iSide)

      ! for the left element, we have to access the right face values
      ! and for the right, we have to acess the left face values.
      ! --> modal space
      modalCoeffs(:,:,1) = 0.0_rk
      modalCoeffs(:,:,2) = 0.0_rk
      do iDegX = 1, poly_proj%min_degree+1
        modalCoeffs(iDegX,:,1) = faceRep(left,iDegX,:,2)
        modalCoeffs(iDegX,:,2) = faceRep(right,iDegX,:,1)
      end do
      ! --> oversmap modal space

      ! transform the 1D modal representation to nodal surface points
      call ply_poly_project_m2n(me = poly_proj, &
       &                       dim = 1 , &
       &                       nVars = 3, &
       &                       nodal_data= pntVal(:,:,1), &
       &                       modal_data= modalCoeffs(:,:,1))
      call ply_poly_project_m2n(me = poly_proj, &
       &                       dim = 1 , &
       &                       nVars = 3, &
       &                       nodal_data= pntVal(:,:,2), &
       &                       modal_data= modalCoeffs(:,:,2))
      ! ... the PML is not involved here
      pntVal(:,4:7,1) = 0.0_rk
      pntVal(:,4:7,2) = 0.0_rk


      ! for each of the surface points calculate the numerical flux
      do iPoint = 1, nquadpoints

        ! Get the material of the left and right element at the current
        ! face.
        left_mu = material_left(iSide,iPoint,1)
        left_epsi = material_left(iSide,iPoint,2)
        right_mu = material_right(iSide,iPoint,1)
        right_epsi = material_right(iSide,iPoint,2)

        ! Calculate the flux for this cube
        call maxwell_flux_cube_2d( left = pntVal(iPoint,var,1), &
                              & left_mu = left_mu, &
                              & left_epsi = left_epsi, &
                              & right = pntVal(iPoint,var,2), &
                              & right_mu = right_mu, &
                              & right_epsi = right_epsi, &
                              & flux = flux )

        ! Rotate back to the origianl face direction
        nodalNumFlux(iPoint,var(1:3)) = flux(1:3)

      end do

      ! transform back to modal space (facial polynomial)
      call ply_poly_project_n2m(me = poly_proj, &
        &                       dim = 1 , &
        &                       nVars = 3, &
        &                       nodal_data= nodalNumFlux(:,:), &
        &                       modal_data= numFluxBuffer(:,:) )
      ! --> oversamp modal space

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      ! check which basisType
      do iDegX = 1, poly_proj%min_degree+1
        ! Fluxes for D and B
        faceFlux(left,iDegX,1:3,2)  = numFluxBuffer(iDegX,1:3)
        faceFlux(right,iDegX,1:3,1) = numFluxBuffer(iDegX,1:3)
        ! No flux for the PML
        faceFlux(left,iDegX,4:7,2)  = 0.0_rk
        faceFlux(right,iDegX,4:7,1) = 0.0_rk
      end do

    end do faceLoop

  end subroutine maxwell_flux_nonconst_cube_vec_2d

  !> summary: calculate flux of pure maxwell equation directly on the face-vector
  !! (formulation is based on TE mode formulation for 2D).
  !!
  !! This subroutine assumes the Maxwell equations with D and B as input
  !! variables. Furthermore, it is able to handle jumping material parameters.
  subroutine maxwell_flux_cube_vec_2d(nTotalFaces, nSides, nFaceDofs, faceRep, &
    &                              faceFlux, leftPos, rightPos, var,        &
    &                              material_left, material_right            )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: nTotalFaces, nFaceDofs, nSides
    real(kind=rk), intent(in) :: faceRep(nTotalFaces,nFaceDofs,7,2)
    real(kind=rk), intent(inout) :: faceFlux(nTotalFaces,nFaceDofs,7,2)
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    integer, intent(in) :: var(7)
    real(kind=rk), intent(in) :: material_left(nSides,1,3)
    real(kind=rk), intent(in) :: material_right(nSides,1,3)
    ! --------------------------------------------------------------------------
    integer :: iSide, left, right, iDof
    real(kind=rk) :: left_mu, right_mu
    real(kind=rk) :: left_epsi, right_epsi
    real(kind=rk) :: left_speedOfLight, right_speedOfLight
    real(kind=rk) :: inv_denom_mu, inv_denom_epsi
    ! --------------------------------------------------------------------------

    ! flux for D_x, B_z, D_y
    do iDof = 1, nFaceDofs

      do iSide = 1, nSides

        ! The position of the left and right element in the state
        ! vector.
        left = leftPos(iSide)
        right = rightPos(iSide)

        ! The material parameters of the left and right element
        left_mu = material_left(iSide,1,1)
        left_epsi = material_left(iSide,1,2)
        left_speedOfLight = 1.0_rk / sqrt( left_mu * left_epsi )
        right_mu = material_right(iSide,1,1)
        right_epsi = material_right(iSide,1,2)
        right_speedOfLight = 1.0_rk / sqrt( right_mu * right_epsi )

        ! The inverse of the denominators
        inv_denom_mu = 1.0_rk / ((-1.0_rk)*left_speedOfLight*left_mu - right_speedOfLight*right_mu)
        inv_denom_epsi = 1.0_rk / (left_speedOfLight*left_epsi + right_speedOfLight*right_epsi)

        ! flux for D_x
        faceFlux(left,iDof,var(1),2) = 0.0_rk

        ! flux for PML
        faceFlux(left,iDof,var(4:7),2) = 0.0_rk

        ! the flux for D_y
        faceFlux(left,iDof,var(2),2) = ( &
          &                  ( ((-1.0_rk)*faceRep(left,iDof,var(2),2) / left_epsi)        &
          &                    - ((-1.0_rk)*faceRep(right,iDof,var(2),1) / right_epsi)  ) &
          &               -  ( left_speedOfLight * faceRep(left,iDof,var(3),2)          &
          &                    + right_speedOfLight * faceRep(right,iDof,var(3),1) ))


        ! the flux for B_z
        faceFlux(left,iDof,var(3),2) = ( &
          &                  ( left_speedOfLight * faceRep(left,iDof,var(2),2)     &
          &                  + right_speedOfLight * faceRep(right,iDof,var(2),1) )&
          &                + ( ( faceRep(left,iDof,var(3),2) / left_mu )     &
          &                  - ( faceRep(right,iDof,var(3),1) / right_mu ) ))

        ! Normalize the calculated fluxes
        faceFlux(left,iDof,var(2),2) = inv_denom_mu * faceFlux(left,iDof,var(2),2)
        faceFlux(left,iDof,var(3),2) = inv_denom_epsi * faceFlux(left,iDof,var(3),2)

        ! Assign the same flux for both adjacent elements
        faceFlux(right,iDof,:,1) = faceFlux(left,iDof,:,2)

      end do

    end do

  end subroutine maxwell_flux_cube_vec_2d

end module atl_maxwell_flux_2d_module
