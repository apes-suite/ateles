! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!! Nernst-Planck equations.
module atl_nerplanck_flux_module
  use env_module, only: rk

  implicit none

  private

  public :: atl_nerplanck_numflux_solve, atl_nerplanck_numflux_preprocess, &
    & atl_nerplanck_physflux_solve, atl_nerplanck_physflux_preprocess


  interface atl_nerplanck_numflux_solve
    module procedure nerplanck_numflux_concentration
  end interface atl_nerplanck_numflux_solve

  interface atl_nerplanck_numflux_preprocess
    module procedure nerplanck_numflux_diffusiveFlux
  end interface atl_nerplanck_numflux_preprocess

  interface atl_nerplanck_physflux_solve
    module procedure nerplanck_physflux_concentration
  end interface atl_nerplanck_physflux_solve

  interface atl_nerplanck_physflux_preprocess
    module procedure nerplanck_physflux_diffusiveFlux
  end interface atl_nerplanck_physflux_preprocess
contains

  !> summary: Subroutine to calculate the numerical flux for the first equation (u)
  !! of the Nernst-Planck equations.
  !!
  !! This subroutine calculates the flux of the Nernst-Planck equation on the reference
  !! cubic face.
  subroutine nerplanck_numflux_concentration(left, right, diffusivity, flux)
    ! --------------------------------------------------------------------------
    !> Left state vector (as conservative variables). The order of this vector
    !! has to be \( (u, \rho_x, \rho_y, \rho_z) \) where \( u \) and
    !! \( \rho \) denoted concentration (scalar) and diffusive fluxes (vector).
    real(kind=rk), intent(in)  :: left(4)
    !> Right state vector (as conservative variables). The order of this vector
    !! has to be \( (u, \rho_x, \rho_y, \rho_z) \) where \( u \) and
    !! \( \rho \) denoted concentration (scalar) and diffusive fluxes (vector).
    real(kind=rk), intent(in)  :: right(4)
    !> The diffusivity in the left and right cell. Since this is
    !! assumed to be equal in both cells, this flux function cannot be used
    !! to calculate the flux for cells with different material properties
    real(kind=rk), intent(in)  :: diffusivity
    !> The flux between left and right cell. The order of this vector is the
    !! same as the input arguments.
    real(kind=rk), intent(out) :: flux
    ! --------------------------------------------------------------------------
    real(kind=rk) :: diffusivitySqrt
    ! --------------------------------------------------------------------------

    diffusivitySqrt = sqrt(diffusivity)

    ! u
    flux = diffusivitySqrt * ( 1.5_rk * right(2) - 0.5_rk * left(2) )

  end subroutine nerplanck_numflux_concentration


  !> summary: Subroutine to calculate the flux for Nernst-Planck equations.
  !!
  !! This subroutine calculates the flux of the Nernst-Planck equation on the reference
  !! cubic face.
  subroutine nerplanck_numflux_diffusiveFlux(left, right, diffusivity, flux)
    ! --------------------------------------------------------------------------
    !> Left state vector (as conservative variables). The order of this vector
    !! has to be \( (u, \rho_x, \rho_y, \rho_z) \) where \( u \) and
    !! \( \rho \) denoted concentration (scalar) and diffusive fluxes (vector).
    real(kind=rk), intent(in)  :: left(4)
    !> Right state vector (as conservative variables). The order of this vector
    !! has to be \( (u, \rho_x, \rho_y, \rho_z) \) where \( u \) and
    !! \( \rho \) denoted concentration (scalar) and diffusive fluxes (vector).
    real(kind=rk), intent(in)  :: right(4)
    !> The diffusivity in the left and right cell. Since this is
    !! assumed to be equal in both cells, this flux function cannot be used
    !! to calculate the flux for cells with different material properties
    real(kind=rk), intent(in)  :: diffusivity
    !> The flux between left and right cell. The order of this vector is the
    !! same as the input arguments.
    real(kind=rk), intent(out) :: flux(3)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: diffusivitySqrt
    ! --------------------------------------------------------------------------

    diffusivitySqrt = sqrt(diffusivity)

    ! rho_x
    flux(1) = diffusivitySqrt * ( 1.5_rk * left(1) - 0.5_rk * right(1) )

    ! rho_y
    flux(2) = 0.0_rk

    ! rho_z
    flux(3) = 0.0_rk

  end subroutine nerplanck_numflux_diffusiveFlux


  !> summary: Subroutine to calculate the numerical flux for the first equation (u)
  !! of the Nernst-Planck equations.
  !!
  !! This subroutine calculates the flux of the Nernst-Planck equation on the reference
  !! cubic element.
  subroutine nerplanck_physflux_concentration(state, diffusivity, flux)
    ! --------------------------------------------------------------------------
    !> State vector (as conservative variables). The order of this vector
    !! has to be \( (u, \rho_x, \rho_y, \rho_z) \) where \( u \) and
    !! \( \rho \) denoted concentration (scalar) and diffusive fluxes (vector).
    real(kind=rk), intent(in)  :: state(4)
    !> The diffusivity in the cell.
    real(kind=rk), intent(in)  :: diffusivity
    !> The flux inside the cell. The order of this vector is the
    !! same as the input arguments.
    real(kind=rk), intent(out) :: flux(3)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: diffusivitySqrt
    ! --------------------------------------------------------------------------

    diffusivitySqrt = sqrt(diffusivity)

    flux(1) = diffusivitySqrt * state(2)
    flux(2) = diffusivitySqrt * state(3)
    flux(3) = diffusivitySqrt * state(4)

  end subroutine nerplanck_physflux_concentration


  !> summary: Subroutine to calculate the flux for Nernst-Planck equations.
  !!
  !! This subroutine calculates the flux of the Nernst-Planck equation on the reference
  !! cubic element.
  subroutine nerplanck_physflux_diffusiveFlux(state, diffusivity, flux)
    ! --------------------------------------------------------------------------
    !> State vector (as conservative variables). The order of this vector
    !! has to be \( (u, \rho_x, \rho_y, \rho_z) \) where \( u \) and
    !! \( \rho \) denoted concentration (scalar) and diffusive fluxes (vector).
    real(kind=rk), intent(in)  :: state(4)
    !> The diffusivity in the cell.
    real(kind=rk), intent(in)  :: diffusivity
    !> The flux inside the cell. The order of this vector is the
    !! same as the input arguments.
    real(kind=rk), intent(out) :: flux
    ! --------------------------------------------------------------------------
    real(kind=rk) :: diffusivitySqrt
    ! --------------------------------------------------------------------------

    diffusivitySqrt = sqrt(diffusivity)

    flux = diffusivitySqrt * state(1)
  end subroutine nerplanck_physflux_diffusiveFlux

end module atl_nerplanck_flux_module
