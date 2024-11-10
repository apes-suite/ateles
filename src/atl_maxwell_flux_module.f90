! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012-2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> author: Jens Zudrop
!! module that holds all routines to calculate the flux for
!! hyperbolic Maxwell equations.
module atl_maxwell_flux_module
  ! Treelm modules
  use env_module,               only: rk

  ! Polynomials modules
  use ply_poly_project_module,  only: ply_poly_project_type, &
    &                                 assignment(=),         &
    &                                 ply_poly_project_m2n,  &
    &                                 ply_poly_project_n2m
  use ply_oversample_module,    only: ply_convert2oversample, &
    &                                 ply_convertFromoversample
  implicit none

  private

  public :: atl_maxwell_flux, atl_maxwell_hc_flux, atl_physFluxMaxwellDivCor

  !> Interface for fluxes of pure Maxwell equations.
  interface atl_maxwell_flux
    module procedure maxwell_flux_cube
    module procedure maxwell_flux_cube_vec
    module procedure maxwell_flux_nonconst_cube_vec
  end interface atl_maxwell_flux

  !> Interface for fluxes of Maxwell equations with
  !! hyperbolic divergence cleaning.
  interface atl_maxwell_hc_flux
    module procedure maxwell_hc_flux_cube
    module procedure maxwell_hc_flux_cube_vec
    module procedure maxwell_hc_flux_nonconst_cube_vec
  end interface atl_maxwell_hc_flux

contains

  ! ************************************************************************ !
  !> Subroutine to calculate the flux for pure Maxwell equations without
  !! any divergence cleaning on the reference cubic face.
  !!
  !! This subroutine calculates the flux of the Maxwell equation on the
  !! reference cubic face. This implementation is based on the this article:
  !! A three-dimensional finite-volume solver for Maxwell equations with
  !! divergence cleaning on unstructured meshes, C.D. Munz, P. Ommes,
  !! R. Schneider, Computer Physiscs communications 130, 83-117, 1999.
  !! Additionally we splitted the correction technique from the Maxwell fluxes
  !! themself to be able to use the fluxes for Maxwell in combination with an
  !! arbitrary divergence correction technique. Please notice that this flux
  !! function assumes constant material parameters in both cells.
  subroutine maxwell_flux_cube( left, right, left_mu, left_epsi, right_mu, &
    &                           right_epsi, flux                           )
    ! -------------------------------------------------------------------- !
    !> Left state vector (as conservative variables). The order of this vector
    !! has to be \f$ (D_x, D_y, D_z, B_1, B_2, B_3) \f$ where E and B denoted
    !! electric field vetor and magnetic field (also called magnetic induction)
    !! vector.
    real(kind=rk), intent(in)  :: left(6)
    !> Right state vector (as conservative variables). The order of this vector
    !! has to be (D_x, D_y, D_z, B_1, B_2, B_3) where E and B denoted the
    !! electric field vetor and magnetic field (also called magnetic induction)
    !! vector.
    real(kind=rk), intent(in)  :: right(6)
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
    real(kind=rk), intent(out) :: flux(6)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: left_speedOfLight, right_speedOfLight
    real(kind=rk) :: inv_denom_mu, inv_denom_epsi
    ! -------------------------------------------------------------------- !

    ! The speed of light in the left and right element
    left_speedOfLight = 1.0_rk / sqrt( left_mu * left_epsi )
    right_speedOfLight = 1.0_rk / sqrt( right_mu * right_epsi )

    ! The inverse of the denominators
    inv_denom_mu = 1.0_rk / ( -left_speedOfLight * left_mu &
      &                       - right_speedOfLight * right_mu         )
    inv_denom_epsi = 1.0_rk / ( left_speedOfLight * left_epsi     &
      &                         + right_speedOfLight * right_epsi )

    ! D_x
    flux(1) = 0.0_rk
    ! B_x
    flux(4) = 0.0_rk

    ! the flux for D_y
    flux(2) = ( ( -left(2) / left_epsi )        &
      &         - ( -right(2) / right_epsi ) )  &
      &       - ( left_speedOfLight * left(6)   &
      &         + right_speedOfLight * right(6) )
    ! the flux for B_z
    flux(6) = ( left_speedOfLight * left(2)       &
      &         + right_speedOfLight * right(2) ) &
      &       + ( ( left(6) / left_mu )           &
      &         - ( right(6) / right_mu ) )
    ! the flux for D_z
    flux(3) = ( ( -left(3) / left_epsi )        &
      &         - ( -right(3) / right_epsi ) )  &
      &       + ( left_speedOfLight * left(5)   &
      &         + right_speedOfLight * right(5) )

    ! the flux for B_y
    flux(5) = ( -left_speedOfLight * left(3)      &
      &         - right_speedOfLight * right(3) ) &
      &       + ( ( left(5) / left_mu )           &
      &         - ( right(5) / right_mu) )

     ! Normalize the calculated fluxes
     flux(2:3) = inv_denom_mu * flux(2:3)
     flux(5:6) = inv_denom_epsi * flux(5:6)

  end subroutine maxwell_flux_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> summary: calculate flux of pure maxwell equation directly on the
  !! face-vector
  !!
  !! This subroutine assumes the Maxwell equations with D and B as input
  !! variables. Furthermore, it is able to handle arbitray material parameters
  !! by a pseudo-spectral technique.
  subroutine maxwell_flux_nonconst_cube_vec( nTotalFaces, nSides, nFaceDofs, &
    &                                        faceRep, faceFlux,              &
    &                                        leftPos, rightPos, var,         &
    &                                        material_left, material_right,  &
    &                                        poly_proj, modalCoeffs, pntVal, &
    &                                        nodalNumFlux, numFluxBuffer     )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: nTotalFaces, nFaceDofs, nSides
    !> The modal representation on the faces, left and right trace.
    real(kind=rk), intent(in) :: faceRep(nTotalFaces,nFaceDofs,6,2)
    !> The fluxes for all faces, for left and right elements.
    real(kind=rk), intent(inout) :: faceFlux(nTotalFaces,nFaceDofs,6,2)
    !> Positions for the left and right elements of all faces
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    !> Variable rotation indices
    integer, intent(in) :: var(6)
    !> Material parameters for the left faces.
    real(kind=rk), intent(in) :: material_left(nSides,nFaceDofs,2)
    !> Material parameters for the right faces.
    real(kind=rk), intent(in) :: material_right(nSides,nFaceDofs,2)
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
    ! -------------------------------------------------------------------- !
    integer :: iSide, left, right, iPoint
    real(kind=rk) :: left_mu, right_mu
    real(kind=rk) :: left_epsi, right_epsi
    real(kind=rk) :: flux(6)
    integer :: nquadpoints, oversamp_dofs
    ! -------------------------------------------------------------------- !
    nquadpoints = poly_proj%body_2d%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    ! Loop over all the faces and calculate the non-zero fluxes
    faceLoop: do iSide = 1, nSides
      ! Get indices of left and right adjacent element
      left = leftPos(iSide)
      right = rightPos(iSide)

      ! --> modal space
      call ply_convert2oversample( state       = faceRep(left,:,:,2), &
        &                          poly_proj   = poly_proj,           &
        &                          nDim        = 2,                   &
        &                          modalCoeffs = modalCoeffs(:,:,1)   )
      call ply_convert2oversample( state       = faceRep(right,:,:,1), &
        &                          poly_proj   = poly_proj,            &
        &                          nDim        = 2,                    &
        &                          modalCoeffs = modalCoeffs(:,:,2)    )
      ! --> oversampling modal space

      ! transform the 2D modal representation to nodal surface points
      call ply_poly_project_m2n( me         = poly_proj,         &
        &                        dim        = 2,                 &
        &                        nVars      = 6,                 &
        &                        nodal_data = pntVal(:,:,1),     &
        &                        modal_data = modalCoeffs(:,:,1) )
      call ply_poly_project_m2n( me         = poly_proj,         &
        &                        dim        = 2,                 &
        &                        nVars      = 6,                 &
        &                        nodal_data = pntVal(:,:,2),     &
        &                        modal_data = modalCoeffs(:,:,2) )


      ! for each of the surface points calculate the numerical flux
      do iPoint = 1, nquadpoints

        ! Get the material of the left and right element at the current
        ! face.
        left_mu = material_left(iSide,iPoint,1)
        left_epsi = material_left(iSide,iPoint,2)
        right_mu = material_right(iSide,iPoint,1)
        right_epsi = material_right(iSide,iPoint,2)

        ! Calculate the flux for this cube
        call maxwell_flux_cube( left       = pntVal(iPoint,var,1), &
          &                     left_mu    = left_mu,              &
          &                     left_epsi  = left_epsi,            &
          &                     right      = pntVal(iPoint,var,2), &
          &                     right_mu   = right_mu,             &
          &                     right_epsi = right_epsi,           &
          &                     flux       = flux                  )

        ! Rotate back to the origianl face direction
        nodalNumFlux(iPoint,var) = flux

      end do

      ! transform back to modal space (facial polynomial)
      call ply_poly_project_n2m( me         = poly_proj,    &
        &                        dim        = 2,            &
        &                        nVars      = 6,            &
        &                        nodal_data = nodalNumFlux, &
        &                        modal_data = numFluxBuffer )
      ! --> oversamp modal space

      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      ! --> oversampled modal space
      call ply_convertFromOversample( modalCoeffs = numFluxBuffer,       &
        &                             poly_proj   = poly_proj,           &
        &                             nDim        = 2,                   &
        &                             state       = faceFlux(left,:,:,2) )

      ! To avoid using the convert routines with the loops again, we just
      ! copy the state also to the right face flux
      faceFlux(right,:,:,1) = faceFlux(left,:,:,2)
      ! --> modal space

    end do faceLoop

  end subroutine maxwell_flux_nonconst_cube_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> summary: calculate flux of pure maxwell equation directly on the
  !! face-vector
  !!
  !! This subroutine assumes the Maxwell equations with D and B as input
  !! variables. Furthermore, it is able to handle jumping material parameters.
  subroutine maxwell_flux_cube_vec( nTotalFaces, nSides, nFaceDofs, faceRep, &
    &                               faceFlux, leftPos, rightPos, var,        &
    &                               material_left, material_right            )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: nTotalFaces, nFaceDofs, nSides
    real(kind=rk), intent(in) :: faceRep(nTotalFaces,nFaceDofs,6,2)
    real(kind=rk), intent(inout) :: faceFlux(nTotalFaces,nFaceDofs,6,2)
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    integer, intent(in) :: var(6)
    real(kind=rk), intent(in) :: material_left(nSides,1,2)
    real(kind=rk), intent(in) :: material_right(nSides,1,2)
    ! -------------------------------------------------------------------- !
    integer :: iSide, left, right, iDof
    real(kind=rk) :: left_mu, right_mu
    real(kind=rk) :: left_epsi, right_epsi
    real(kind=rk) :: left_speedOfLight, right_speedOfLight
    real(kind=rk) :: inv_denom_mu, inv_denom_epsi
    ! -------------------------------------------------------------------- !

    ! flux for D_X, B_Z, D_y and B_z
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
        inv_denom_mu = 1.0_rk                                              &
          & / (-left_speedOfLight * left_mu - right_speedOfLight * right_mu)
        inv_denom_epsi = 1.0_rk                                           &
          & / (left_speedOfLight*left_epsi + right_speedOfLight*right_epsi)

        ! flux for D_x, B_x
        faceFlux(left,iDof,var(1),2) = 0.0_rk
        faceFlux(left,iDof,var(4),2) = 0.0_rk

        ! the flux for D_y
        faceFlux(left,iDof,var(2),2) =                            &
          & ( ( -faceRep(left,iDof,var(2),2) / left_epsi)         &
          &   - (-faceRep(right,iDof,var(2),1) / right_epsi) )    &
          & - ( left_speedOfLight * faceRep(left,iDof,var(6),2)   &
          &   + right_speedOfLight * faceRep(right,iDof,var(6),1) )


        ! the flux for B_z
        faceFlux(left,iDof,var(6),2) =                              &
          & ( left_speedOfLight * faceRep(left,iDof,var(2),2)       &
          &   + right_speedOfLight * faceRep(right,iDof,var(2),1) ) &
          & + ( faceRep(left,iDof,var(6),2) / left_mu               &
          &   - faceRep(right,iDof,var(6),1) / right_mu )

        ! the flux for D_z
        faceFlux(left,iDof,var(3),2) =                            &
          & ( -faceRep(left,iDof,var(3),2) / left_epsi            &
          &   - (-faceRep(right,iDof,var(3),1) / right_epsi) )    &
          & + ( left_speedOfLight * faceRep(left,iDof,var(5),2)   &
          &   + right_speedOfLight * faceRep(right,iDof,var(5),1) )

        ! the flux for B_y
        faceFlux(left,iDof,var(5),2) =                              &
          & ( -left_speedOfLight * faceRep(left,iDof,var(3),2)      &
          &   - right_speedOfLight * faceRep(right,iDof,var(3),1) ) &
          & + ( faceRep(left,iDof,var(5),2) / left_mu               &
          &   - faceRep(right,iDof,var(5),1) / right_mu)

        ! Normalize the calculated fluxes
        faceFlux(left,iDof,var(2),2) =                &
          & inv_denom_mu * faceFlux(left,iDof,var(2),2)
        faceFlux(left,iDof,var(6),2) =                  &
          & inv_denom_epsi * faceFlux(left,iDof,var(6),2)
        faceFlux(left,iDof,var(3),2) =                &
          & inv_denom_mu * faceFlux(left,iDof,var(3),2)
        faceFlux(left,iDof,var(5),2) =                  &
          & inv_denom_epsi * faceFlux(left,iDof,var(5),2)

        ! Assign the same flux for both adjacent elements
        faceFlux(right,iDof,:,1) = faceFlux(left,iDof,:,2)

      end do

    end do

  end subroutine maxwell_flux_cube_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> summary: Subroutine to calculate the flux for pure Maxwell equations with
  !! hyperbolic divergence cleaning on the reference cubic face.
  subroutine maxwell_hc_flux_cube(left, right, mat_left, &
                              & mat_right, flux)
    ! -------------------------------------------------------------------- !
    !> Left state vector (as conservative variables). The order of this vector
    !! has to be \f$ (D_x, D_y, D_z, B_1, B_2, B_3, phi, psi) \f$ where E and B
    !! denoted electric field vetor and magnetic field (also called magnetic
    !! induction) vector.
    real(kind=rk), intent(in)  :: left(8)
    !> Right state vector (as conservative variables). The order of this vector
    !! has to be (D_x, D_y, D_z, B_1, B_2, B_3, phi, psi) where E and B denoted
    !! the electric field vetor and magnetic field (also called magnetic
    !! induction) vector.
    real(kind=rk), intent(in)  :: right(8)
    !> Material for the left face
    real(kind=rk), intent(in)  :: mat_left(4)
    !> Material for the right face
    real(kind=rk), intent(in)  :: mat_right(4)
!!    !>  The magnetic permeability of the left element.
!!    real(kind=rk), intent(in)  :: left_mu
!!    !>  The electric permitivity of the left element.
!!    real(kind=rk), intent(in)  :: left_epsi
!!    !>  Parameter for the magnetic correction on the left element.
!!    real(kind=rk), intent(in)  :: left_gam
!!    !>  Parameter for the electric correction on the left element.
!!    real(kind=rk), intent(in)  :: left_chi
!!    !>  The magnetic permeability of the right element.
!!    real(kind=rk), intent(in)  :: right_mu
!!    !>  The electric permitivity of the right element.
!!    real(kind=rk), intent(in)  :: right_epsi
!!    !>  Parameter for the magnetic correction on the right element.
!!    real(kind=rk), intent(in)  :: right_gam
!!    !>  Parameter for the electric correction on the right element.
!!    real(kind=rk), intent(in)  :: right_chi
    !> The flux between left and right cell. The order of this vector is the
    !! same as the input arguments.
    real(kind=rk), intent(out) :: flux(8)
    ! -------------------------------------------------------------------- !
!! JZ: Old implementation of the flux. There must be an error somewhere. I have to
!!     check my solution for the Riemann problem agian. So, I replaced the flux
!!     by a simple Lax-Friedrich type flux.
!!
!!    real(kind=rk) :: left_speedOfLight, right_speedOfLight
!!    real(kind=rk) :: inv_denom_mu, inv_denom_epsi
!!    ! --------------------------------------------------------------------------
!!
!!    ! The speed of light in the left and right element
!!    left_speedOfLight = 1.0_rk / sqrt( left_mu * left_epsi )
!!    right_speedOfLight = 1.0_rk / sqrt( right_mu * right_epsi )
!!
!!    ! The inverse of the denominators
!!    inv_denom_mu = 1.0_rk / ((-1.0_rk)*left_speedOfLight*left_mu - right_speedOfLight*right_mu)
!!    inv_denom_epsi = 1.0_rk / (left_speedOfLight*left_epsi + right_speedOfLight*right_epsi)
!!
!!    ! D_x
!!    flux(1) = inv_denom_mu * (                                   &
!!              &       (-1.0_rk)*left_chi*left(1)/left_epsi       &
!!              &     - (-1.0_rk)*right_chi*right(1)/right_epsi    &
!!              &              )                                   &
!!              & + (                                              &
!!              &        left_chi*left_chi*left(7)                 &
!!              &      + right_chi*right_chi*right(7)              &
!!              &   ) / ( left_epsi*left_mu + right_epsi*right_mu  )
!!    ! B_x
!!    flux(4) = inv_denom_epsi * (                                 &
!!              &        (-1.0_rk)*left_gam*left(4)/left_mu        &
!!              &      - (-1.0_rk)*right_gam*right(4)/right_mu)    &
!!              & + (                                              &
!!              &       left_gam*left_gam*left(8)                  &
!!              &     + right_gam*right_gam*right(8)               &
!!              &   ) / ( left_epsi*left_mu + right_epsi*right_mu  )
!!
!!    ! the flux for phi (electric correction)
!!    flux(7) = (                                            &
!!              &    left_speedOfLight*left(1)               &
!!              &  + right_speedOfLight*right(1)             &
!!              &  + left_speedOfLight*left_chi*left(7)      &
!!              &  - right_speedOfLight*right_chi*right(7)   &
!!              & ) / (left_speedOfLight + right_speedOfLight)
!!
!!    ! the flux for psi (magnetic correction)
!!    flux(8) = (                                            &
!!              &    left_speedOfLight*left(4)               &
!!              &  + right_speedOfLight*right(4)             &
!!              &  + left_speedOfLight*left_gam*left(8)      &
!!              &  - right_speedOfLight*right_gam*right(8)   &
!!              & ) / (left_speedOfLight + right_speedOfLight)
!!
!!    ! the flux for D_y
!!    flux(2) = ( &
!!                &    ( (-1.0_rk*left(2) / left_epsi)        &
!!                &      - (-1.0_rk*right(2) / right_epsi)  ) &
!!                & -  ( left_speedOfLight * left(6)          &
!!                &      + right_speedOfLight * right(6) ))
!!    ! the flux for B_z
!!    flux(6) = ( &
!!                &   ( left_speedOfLight * left(2)           &
!!                &   + right_speedOfLight * right(2) )       &
!!                & + ( ( left(6) / left_mu )                 &
!!                &   - ( right(6) / right_mu ) ))
!!    ! the flux for D_z
!!    flux(3) = ( &
!!                &   ( ( -1.0_rk * left(3) / left_epsi )     &
!!                &   - ( -1.0_rk * right(3) / right_epsi ) ) &
!!                & + ( left_speedOfLight * left(5)           &
!!                &   + right_speedOfLight * right(5) )       &
!!                &             )
!!
!!    ! the flux for B_y
!!    flux(5) = ( &
!!                &   ( -1.0_rk * left_speedOfLight * left(3) &
!!                &   - right_speedOfLight * right(3) )       &
!!                & + ( ( left(5) / left_mu )                 &
!!                &   - ( right(5) / right_mu) )              &
!!                &             )
!!
!!     ! Normalize the calculated fluxes
!!     flux(2:3) = inv_denom_mu * flux(2:3)
!!     flux(5:6) = inv_denom_epsi * flux(5:6)

    real(kind=rk) :: maxSpeed, maxSpeedLeft, maxSpeedRight

    ! Upper bound on the local wave speed
    maxSpeedLeft = max( max( 1.0 / sqrt(mat_left(1) * mat_left(2)), &
      &   mat_left(3) / sqrt(mat_left(1) * mat_left(2) ) ),         &
      & mat_left(4) / sqrt(mat_left(1) * mat_left(2))               )
    maxSpeedRight = max( max( 1.0 / sqrt(mat_right(1) * mat_right(2)), &
      &   mat_right(3) / sqrt(mat_right(1) * mat_right(2)) ),          &
      & mat_right(4) / sqrt(mat_right(1) * mat_right(2))               )
    maxSpeed = max( maxSpeedLeft, maxSpeedRight )

    ! Lax-Friedrich flux
    flux(:) =                                                    &
    & ( atl_physFluxMaxwellDivCor(left, mat_left)                &
    &   + atl_physFluxMaxwellDivCor(right, mat_right) ) / 2.0_rk &
    & + (maxSpeed / 2.0_rk) * ( left - right )

  end subroutine maxwell_hc_flux_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> summary: calculate flux of maxwell equation with hyperbolic divergence
  !! cleaning directly on the face-vector
  !!
  !! This subroutine assumes the Maxwell equations with D and B as input
  !! variables. Furthermore, it is able to handle jumping material parameters.
  subroutine maxwell_hc_flux_cube_vec( nTotalFaces, nSides, nFaceDofs,faceRep, &
    &                                  faceFlux, leftPos, rightPos, var,       &
    &                                  material_left, material_right           )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: nTotalFaces, nFaceDofs, nSides
    real(kind=rk), intent(in) :: faceRep(nTotalFaces,nFaceDofs,8,2)
    real(kind=rk), intent(inout) :: faceFlux(nTotalFaces,nFaceDofs,8,2)
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    integer, intent(in) :: var(8)
    real(kind=rk), intent(in) :: material_left(nSides,1,4)
    real(kind=rk), intent(in) :: material_right(nSides,1,4)
    ! -------------------------------------------------------------------- !
!! JZ: Old implementation of the flux. There must be an error somewhere. I have to
!!     check my solution for the Riemann problem agian. So, I replaced the flux
!!     by a simple Lax-Friedrich type flux.
!!
!!  integer :: iSide, left, right, iDof
!!  real(kind=rk) :: left_mu, right_mu
!!  real(kind=rk) :: left_epsi, right_epsi
!!  real(kind=rk) :: left_gam, right_gam
!!  real(kind=rk) :: left_chi, right_chi
!!  real(kind=rk) :: left_speedOfLight, right_speedOfLight
!!  real(kind=rk) :: inv_denom_mu, inv_denom_epsi
!!  ! --------------------------------------------------------------------------
!!
!!  ! flux for D_X, B_Z, D_y and B_z
!!  do iDof = 1, nFaceDofs
!!
!!    do iSide = 1, nSides
!!
!!      ! The position of the left and right element in the state
!!      ! vector.
!!      left = leftPos(iSide)
!!      right = rightPos(iSide)
!!
!!      ! The material parameters of the left and right element
!!      left_mu = material_left(iSide,1,1)
!!      left_epsi = material_left(iSide,1,2)
!!      left_gam = material_left(iSide,1,3)
!!      left_chi = material_left(iSide,1,4)
!!      left_speedOfLight = 1.0_rk / sqrt( left_mu * left_epsi )
!!      right_mu = material_right(iSide,1,1)
!!      right_epsi = material_right(iSide,1,2)
!!      right_gam = material_right(iSide,1,3)
!!      right_chi = material_right(iSide,1,4)
!!      right_speedOfLight = 1.0_rk / sqrt( right_mu * right_epsi )
!!
!!      ! The inverse of the denominators
!!      inv_denom_mu = 1.0_rk / ((-1.0_rk)*left_speedOfLight*left_mu - right_speedOfLight*right_mu)
!!      inv_denom_epsi = 1.0_rk / (left_speedOfLight*left_epsi + right_speedOfLight*right_epsi)
!!
!!      ! flux for D_x
!!      faceFlux(left,iDof,var(1),2) = inv_denom_mu * (                                   &
!!                                     &       (-1.0_rk)*left_chi*faceRep(left,iDof,var(1),2)/left_epsi    &
!!                                     &     - (-1.0_rk)*right_chi*faceRep(right,iDof,var(1),1)/right_epsi &
!!                                     &              )                                   &
!!                                     & + (                                              &
!!                                     &        left_chi*left_chi*faceRep(left,iDof,var(7),2) &
!!                                     &      + right_chi*right_chi*faceRep(right,iDof,var(7),1) &
!!                                     &   ) / ( left_epsi*left_mu + right_epsi*right_mu  )
!!
!!      ! flux for B_x
!!      faceFlux(left,iDof,var(4),2) = inv_denom_epsi * (                                 &
!!                                     &        (-1.0_rk)*left_gam*faceRep(left,iDof,var(4),2)/left_mu &
!!                                     &      - (-1.0_rk)*right_gam*faceRep(right,iDof,var(4),1)/right_mu) &
!!                                     & + (                                              &
!!                                     &       left_gam*left_gam*faceRep(left,iDof,var(8),2) &
!!                                     &     + right_gam*right_gam*faceRep(right,iDof,var(8),1) &
!!                                     &   ) / ( left_epsi*left_mu + right_epsi*right_mu  )
!!
!!      ! flux for phi (electric correction)
!!      faceFlux(left,iDof,var(7),2) = (                                            &
!!                                     &    left_speedOfLight*faceRep(left,iDof,var(1),2) &
!!                                     &  + right_speedOfLight*faceRep(right,iDof,var(1),1) &
!!                                     &  + left_speedOfLight*left_chi*faceRep(left,iDof,var(7),2) &
!!                                     &  - right_speedOfLight*right_chi*faceRep(right,iDof,var(7),1)   &
!!                                     & ) / (left_speedOfLight + right_speedOfLight)
!!
!!      ! flux for psi (magnetic correction)
!!      faceFlux(left,iDof,var(8),2) = (                                            &
!!                                     &    left_speedOfLight*faceRep(left,iDof,var(4),2)               &
!!                                     &  + right_speedOfLight*faceRep(right,iDof,var(4),1)             &
!!                                     &  + left_speedOfLight*left_gam*faceRep(left,iDof,var(8),2)      &
!!                                     &  - right_speedOfLight*right_gam*faceRep(right,iDof,var(8),1)   &
!!                                     & ) / (left_speedOfLight + right_speedOfLight)
!!
!!      ! the flux for D_y
!!      faceFlux(left,iDof,var(2),2) = ( &
!!        &                  ( ((-1.0_rk)*faceRep(left,iDof,var(2),2) / left_epsi)        &
!!        &                    - ((-1.0_rk)*faceRep(right,iDof,var(2),1) / right_epsi)  ) &
!!        &               -  ( left_speedOfLight * faceRep(left,iDof,var(6),2)          &
!!        &                    + right_speedOfLight * faceRep(right,iDof,var(6),1) ))
!!
!!
!!      ! the flux for B_z
!!      faceFlux(left,iDof,var(6),2) = ( &
!!        &                  ( left_speedOfLight * faceRep(left,iDof,var(2),2)     &
!!        &                  + right_speedOfLight * faceRep(right,iDof,var(2),1) )&
!!        &                + ( ( faceRep(left,iDof,var(6),2) / left_mu )     &
!!        &                  - ( faceRep(right,iDof,var(6),1) / right_mu ) ))
!!
!!      ! the flux for D_z
!!      faceFlux(left,iDof,var(3),2) = ( &
!!        &                  ( ( (-1.0_rk) * faceRep(left,iDof,var(3),2) / left_epsi )     &
!!        &                  - ( (-1.0_rk) * faceRep(right,iDof,var(3),1) / right_epsi ) ) &
!!        &                + ( left_speedOfLight * faceRep(left,iDof,var(5),2)           &
!!        &                  + right_speedOfLight * faceRep(right,iDof,var(5),1) )       &
!!        &                            )
!!
!!      ! the flux for B_y
!!      faceFlux(left,iDof,var(5),2) = ( &
!!        &                  ( (-1.0_rk) * left_speedOfLight * faceRep(left,iDof,var(3),2) &
!!        &                  - right_speedOfLight * faceRep(right,iDof,var(3),1) )       &
!!        &                + ( ( faceRep(left,iDof,var(5),2) / left_mu )                 &
!!        &                  - ( faceRep(right,iDof,var(5),1) / right_mu) )              &
!!        &                            )
!!
!!      ! Normalize the calculated fluxes
!!      faceFlux(left,iDof,var(2),2) = inv_denom_mu * faceFlux(left,iDof,var(2),2)
!!      faceFlux(left,iDof,var(6),2) = inv_denom_epsi * faceFlux(left,iDof,var(6),2)
!!      faceFlux(left,iDof,var(3),2) = inv_denom_mu * faceFlux(left,iDof,var(3),2)
!!      faceFlux(left,iDof,var(5),2) = inv_denom_epsi * faceFlux(left,iDof,var(5),2)
!!
!!      ! Assign the same flux for both adjacent elements
!!      faceFlux(right,iDof,:,1) = faceFlux(left,iDof,:,2)
!!
!!    end do
!!
!!  end do

    integer :: iSide, left, right, iDof
    real(kind=rk) :: maxSpeedLeft, maxSpeedRight, maxSpeed

    do iDof = 1, nFaceDofs

      do iSide = 1, nSides

        ! The position of the left and right element in the state
        ! vector.
        left = leftPos(iSide)
        right = rightPos(iSide)

        ! Upper bound on the local wave speed
        maxSpeedLeft =                            &
          & max(                                  &
          &   max(                                &
          &     1.0                               &
          &       / sqrt(material_left(iSide,1,1) &
          &       * material_left(iSide,1,2)),    &
          &     material_left(iSide,1,3)          &
          &       / sqrt(material_left(iSide,1,1) &
          &       * material_left(iSide,1,2)) ),  &
          &   material_left(iSide,1,4)            &
          &     / sqrt(material_left(iSide,1,1)   &
          &     * material_left(iSide,1,2))       )
        maxSpeedRight =                              &
          & max(                                     &
          &   max(                                   &
          &     1.0                                  &
          &       / sqrt( material_right(iSide,1,1)  &
          &         * material_right(iSide,1,2) ),   &
          &     material_right(iSide,1,3)            &
          &       / sqrt(material_right(iSide,1,1)   &
          &         * material_right(iSide,1,2) ) ), &
          &   material_right(iSide,1,4)              &
          &     / sqrt(material_right(iSide,1,1)     &
          &       * material_right(iSide,1,2) )      )
        maxSpeed = max( maxSpeedLeft, maxSpeedRight)

        ! Lax-Friedrich flux
        faceFlux(left,iDof,var,2) =                                      &
          & ( atl_physFluxMaxwellDivCor( faceRep(left,iDof,var,2),       &
          &                              material_left(iSide,1,:) )      &
          &   + atl_physFluxMaxwellDivCor( faceRep(right,iDof,var,1),    &
          &                                material_right(iSide,1,:) ) ) &
          & / 2.0_rk                                                     &
          & + (maxSpeed/2.0_rk)                                          &
          & * ( faceRep(left,iDof,var,2) - faceRep(right,iDof,var,1)     )

        !! Assign the same flux for left and right element
        faceFlux(right,iDof,:,1) = faceFlux(left,iDof,:,2)

      end do
    end do

  end subroutine maxwell_hc_flux_cube_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> summary: calculate flux of maxwell equation with hyperbolic divergence
  !! cleaning directly on the face-vector
  !!
  !! This subroutine assumes the Maxwell equations with D and B as input
  !! variables. Furthermore, it is able to handle arbitray material parameters
  !! by a pseudo-spectral technique.
  subroutine maxwell_hc_flux_nonconst_cube_vec( nTotalFaces, nSides,       &
    & nFaceDofs, faceRep, faceFlux, leftPos, rightPos, var, material_left, &
    & material_right, poly_proj, left_modalCoeffs, right_modalCoeffs,      &
    & left_pntVal, right_pntVal, nodalNumFlux, numFluxBuffer               )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: nTotalFaces, nFaceDofs, nSides
    !> The modal representation on the faces, left and right trace.
    real(kind=rk), intent(in) :: faceRep(nTotalFaces,nFaceDofs,8,2)
    !> The fluxes for all faces, for left and right elements.
    real(kind=rk), intent(inout) :: faceFlux(nTotalFaces,nFaceDofs,8,2)
    !> Positions for the left and right elements of all faces
    integer, intent(in) :: leftPos(nSides), rightPos(nsides)
    !> Variable rotation indices
    integer, intent(in) :: var(8)
    !> Material parameters for the left faces.
    real(kind=rk), intent(in) :: material_left(nSides,nFaceDofs,4)
    !> Material parameters for the right faces.
    real(kind=rk), intent(in) :: material_right(nSides,nFaceDofs,4)
   ! !> FPT to convert between nodes and modes.
   ! type(atl_legFpt_2D_type), intent(inout) :: fpt
    !> Data for projection method
    type(ply_poly_project_type) :: poly_proj
!!    !> Working array for the left and right modal coefficients
!!    real(kind=rk), intent(inout) :: left_modalCoeffs((fpt%nQuadPoints)**2,8)
!!    real(kind=rk), intent(inout) :: right_modalCoeffs((fpt%nQuadPoints)**2,8)
!!    !> Working array for the left and right point values
!!    real(kind=rk), intent(inout) :: left_pntVal((fpt%nQuadPoints)**2,8)
!!    real(kind=rk), intent(inout) :: right_pntVal((fpt%nQuadPoints)**2,8)
!!    !> Working array for the nodal flux
!!    real(kind=rk), intent(inout) :: nodalNumFlux((fpt%nQuadPoints)**2,8)
!!    !> Working array for the modal numerical flux
!!    real(kind=rk), intent(inout) :: numFluxBuffer((fpt%nQuadPoints)**2,8)
    !> Working array for the left and right modal coefficients
    real(kind=rk), allocatable, intent(inout) :: left_modalCoeffs(:,:)
    real(kind=rk), allocatable, intent(inout) :: right_modalCoeffs(:,:)
    !> Working array for the left and right point values
    real(kind=rk), allocatable, intent(inout) :: left_pntVal(:,:)
    real(kind=rk), allocatable, intent(inout) :: right_pntVal(:,:)
    !> Working array for the nodal flux
    real(kind=rk), allocatable, intent(inout) :: nodalNumFlux(:,:)
    !> Working array for the modal numerical flux
    real(kind=rk), allocatable, intent(inout) :: numFluxBuffer(:,:)
    ! -------------------------------------------------------------------- !
    integer :: iSide, left, right, iPoint
    !!real(kind=rk) :: left_mu, right_mu
    !!real(kind=rk) :: left_epsi, right_epsi
    !!real(kind=rk) :: left_gam, right_gam
    !!real(kind=rk) :: left_chi, right_chi
    real(kind=rk) :: flux(8)
    integer :: nquadpoints, ndofs, oversamp_dofs
    ! -------------------------------------------------------------------- !

    nquadpoints = poly_proj%body_2D%nquadpoints
    ndofs = poly_proj%body_2D%ndofs
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    allocate(left_modalCoeffs(oversamp_dofs,8))
    allocate(right_modalCoeffs(oversamp_dofs,8))
    allocate(left_pntVal(nquadpoints,8))
    allocate(right_pntVal(nquadpoints,8))
    allocate(nodalNumFlux(nquadpoints,8))
    allocate(numFluxBuffer(ndofs,8))

    ! Loop over all the faces and calculate the non-zero fluxes
    faceLoop: do iSide = 1, nSides
      ! Get indices of left and right adjacent element
      left = leftPos(iSide)
      right = rightPos(iSide)

      ! --> modal space
      call ply_convert2oversample( state       = faceRep(left,:,:,2), &
        &                          poly_proj   = poly_proj,           &
        &                          nDim        = 2,                   &
        &                          modalCoeffs = left_modalCoeffs     )
      call ply_convert2oversample( state       = faceRep(right,:,:,1), &
        &                          poly_proj   = poly_proj,            &
        &                          nDim        = 2,                    &
        &                          modalCoeffs = right_modalCoeffs     )
      ! --> oversampled modal space

      ! transform the 2D modal representation to nodal surface points
      call ply_poly_project_m2n( me         = poly_proj,       &
        &                        dim        = 2,               &
        &                        nVars      = 8,               &
        &                        nodal_data = left_pntVal,     &
        &                        modal_data = left_modalCoeffs )
      call ply_poly_project_m2n( me         = poly_proj,        &
        &                        dim        = 2,                &
        &                        nVars      = 8,                &
        &                        nodal_data = right_pntVal,     &
        &                        modal_data = right_modalCoeffs )
      ! --> oversamp nodal space


      ! for each of the surface points calculate the numerical flux
      do iPoint = 1, nquadpoints

        !!! Get the material of the left and right element at the current
        !!! face.
        !!! ... left material
        !!left_mu = material_left(iSide,iPoint,1)
        !!left_epsi = material_left(iSide,iPoint,2)
        !!left_gam = material_left(iSide,iPoint,3)
        !!left_chi = material_left(iSide,iPoint,4)
        !!! ... right material
        !!right_mu = material_right(iSide,iPoint,1)
        !!right_epsi = material_right(iSide,iPoint,2)
        !!right_gam = material_right(iSide,iPoint,3)
        !!right_chi = material_right(iSide,iPoint,4)

        ! Calculate the flux for this cube
        call maxwell_hc_flux_cube( left      = left_pntVal(iPoint,var),        &
          &                        mat_left  = material_left(iSide,iPoint,:),  &
          &                        mat_right = material_right(iSide,iPoint,:), &
                             !! & left_mu = left_mu, &
                             !! & left_epsi = left_epsi, &
                             !! & left_gam = left_gam, &
                             !! & left_chi = left_chi, &
          &                        right     = right_pntVal(iPoint,var),       &
                             !! & right_mu = right_mu, &
                             !! & right_epsi = right_epsi, &
                             !! & right_gam = right_gam, &
                             !! & right_chi = right_chi, &
          &                        flux      = flux                            )

        ! Rotate back to the origianl face direction
        nodalNumFlux(iPoint,var) = flux

      end do

      ! transform back to modal space (facial polynomial)
      call ply_poly_project_n2m( me         = poly_proj,    &
        &                        dim        = 2,            &
        &                        nVars      = 8,            &
        &                        nodal_data = nodalNumFlux, &
        &                        modal_data = numFluxBuffer )
      ! --> oversamp modal space
      ! Store the modal coefficients of the numerical flux. For the left
      ! element we have calculated the flux on the right face and vice versa.
      call ply_convertFromOversample( modalCoeffs = numFluxBuffer,       &
        &                             poly_proj   = poly_proj,           &
        &                             nDim        = 2,                   &
        &                             state       = faceFlux(left,:,:,2) )
      ! To avoid using the convert routines with the loops again, we just
      ! copy the state also to the right face flux
      faceFlux(right,:,:,1) = faceFlux(left,:,:,2)
      ! --> modal space

    end do faceLoop

  end subroutine maxwell_hc_flux_nonconst_cube_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Function for physical flux of the Maxwell equations in terms of D and B.
  function atl_physFluxMaxwellDivCor(state, material) result(flux)
    ! -------------------------------------------------------------------- !
    !> State to compute the fluxes from (D,B).
    real(kind=rk), intent(in) :: state(8)
    !> Material parameters (mu, epsilon) the flux calculation
    real(kind=rk), intent(in) :: material(4)
    !> The resulting flux in x direction
    real(kind=rk) :: flux(8)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: inv_mu, inv_epsi, chi, gam
    ! -------------------------------------------------------------------- !

    ! Get magnetic permeability
    inv_mu = 1.0_rk / material(1)
    ! Get electric permitivity
    inv_epsi = 1.0_rk / material(2)
    ! gamma, parameter for magnetic correction
    gam = material(3)
    ! chi, parameter for electric correction
    chi = material(4)

    ! 1...3 displacement field D
    flux(1) = state(7) * inv_mu * inv_epsi * chi * chi
    flux(2) = state(6) * inv_mu
    flux(3) = state(5) * (-1.0_rk) * inv_mu
    ! 4...6 magnetic field B
    flux(4) = state(8) * inv_mu * inv_epsi * gam * gam
    flux(5) = state(3) * (-1.0_rk) * inv_epsi
    flux(6) = state(2) * inv_epsi
    ! 7 electric correction
    flux(7) = state(1)
    ! 8 magnetic correction
    flux(8) = state(4)

  end function atl_physFluxMaxwellDivCor
  ! ************************************************************************ !


end module atl_maxwell_flux_module
