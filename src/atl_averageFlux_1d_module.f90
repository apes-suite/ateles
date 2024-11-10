! Copyright (c) 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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
!! Module collects all Lax-Friedrich flux for different types of equations.
module atl_averageFlux_1d_module
  ! Treelm modules
  use env_module,                   only: rk
  ! Ateles modules
  use atl_physFluxEuler_1d_module,  only: atl_physFluxEuler_1d
  use atl_equation_module,          only: atl_equations_type

  implicit none
  private

  public :: atl_averageFluxEuler_1d

contains

  !> A most basic flux function which uses the average of both sides for
  !! 1D Euler.
  !!
  !! This interface has to match the abstract definition compute_numflux in the
  !! atl_equation_module.
  subroutine atl_averageFluxEuler_1d( equation, state_left, state_right,      &
    &                                 material_left, material_right, nPoints, &
    &                                 flux                                    )
    ! ---------------------------------------------------------------------------
    class(atl_equations_type), intent(in) :: equation
    !> The state on the face from its left limit (in conservative variables).
    real(kind=rk), intent(in) :: state_left(:,:)
    !> The state on the face from its right limit (in conservative variables).
    real(kind=rk), intent(in) :: state_right(:,:)
    !> The left value of the characteristic function (stemming from
    !! penalization)
    real(kind=rk), intent(in) :: material_left(:,:)
    !> The right value of the characteristic function (stemming from
    !! penalization)
    real(kind=rk), intent(in) :: material_right(:,:)
    !> Number of points to evaluate the flux at.
    integer, intent(in) :: nPoints
    !> Resulting flux for the left element (in conservative variables).
    real(kind=rk), intent(out) :: flux(:,:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: isen_coeff
    real(kind=rk) :: state_avg(3)
    real(kind=rk) :: material(ubound(material_left,1))
    integer :: iPoint
    integer :: matpoint
    integer :: mm
    ! ---------------------------------------------------------------------------

    mm = ubound(material_left,1)
    isen_coeff = equation%euler%isen_coef
    material = 0.5_rk * (material_left(:,1) + material_right(:,1)) &
      &               * (1.0_rk/equation%euler%porosity - 1.0_rk)

    do iPoint = 1, nPoints
      matpoint = min(iPoint, mm)
      state_avg = 0.5_rk * ( state_left(iPoint, :) + state_right(iPoint, :) )
      flux(iPoint,:) = atl_physFluxEuler_1d(                      &
        & state           = state_avg,                            &
        & isenCoeff       = isen_coeff,                           &
        & U_o             = 0.5_rk * ( material_left(matpoint, 2) &
        &                     + material_right(matpoint, 2) ),    &
        & penalty_scaling = material(matpoint)                    )
    end do

  end subroutine atl_averageFluxEuler_1d

end module atl_averageFlux_1d_module
