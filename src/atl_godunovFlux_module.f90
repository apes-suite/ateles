! Copyright (c) 2013, 2015-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017-2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

module atl_GodunovFlux_module
  use env_module,                     only: rk
  use atl_eqn_euler_module,           only: atl_euler_type
  use atl_physFluxEuler_module,       only: atl_physFluxEuler
  use atl_physFluxEuler_1d_module,    only: atl_physFluxEuler_1d
  use atl_physFluxEuler_2d_module,    only: atl_physFluxEuler_2d
  use atl_exact_riemann_euler_module, only: atl_ere_solState1D_type, &
    &                                       atl_ere_init_consts, &
    &                                       atl_ere_eval_onEdge

  implicit none

  private

  public :: atl_GodunovEuler
  public :: atl_GodunovEuler2D
  public :: atl_GodunovEuler1D
  public :: atl_flux_initGodunov

  type(atl_ere_solState1D_type) :: eqconsts


contains


  !> Godunov flux for the Euler equation.
  !!
  !! @todo Implement correct Riemann problem for the equation with porosity.
  !!       Right now the plain Euler problem is solved and the porosity and
  !!       penalty_char are just averaged centrally.
  subroutine atl_GodunovEuler(euler, state_left, state_right, &
    &                         material_left, material_right, nPoints, flux)
    ! ---------------------------------------------------------------------------
    class(atl_euler_type), intent(in) :: euler
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
    real(kind=rk) :: p_left, p_right
    real(kind=rk) :: v_left(3), v_right(3)
    real(kind=rk) :: estate(5), pstate(5)
    real(kind=rk) :: isen_coeff
    real(kind=rk) :: icm1
    integer :: iPoint
    integer :: matpoint
    integer :: mm

    isen_coeff = euler%isen_coef
    icm1 = isen_coeff-1.0_rk
    mm = ubound(material_left,1)

    do iPoint = 1, nPoints
      matpoint = min(iPoint, mm)
      p_Left  = icm1 * ( state_left(iPoint,5)                   &
        &               - 0.5_rk*sum(state_left(iPoint,2:4)**2) &
        &                       / state_left(iPoint,1)          )

      p_Right = icm1 * ( state_right(iPoint,5)                   &
        &               - 0.5_rk*sum(state_right(iPoint,2:4)**2) &
        &                       / state_right(iPoint,1)          )

      v_left = state_left(iPoint,2:4) / state_left(iPoint,1)
      v_right = state_right(iPoint,2:4) / state_right(iPoint,1)

      call atl_ere_eval_onEdge(me = eqconsts, &
        &                      rho_left = state_left(iPoint,1), &
        &                      vn_left = v_left(1), &
        &                      v1_left = v_left(2), v2_left = v_left(3), &
        &                      p_left = p_left, &
        &                      rho_right = state_right(iPoint,1), &
        &                      vn_right = v_right(1), &
        &                      v1_right = v_right(2), v2_right = v_right(3), &
        &                      p_right = p_right, &
        &                      rho = pstate(1), vn = pstate(2), v1 = pstate(3), &
        &                      v2 = pstate(4), p = pstate(5))
      estate(1)   = pstate(1)
      estate(2:4) = pstate(1)*pstate(2:4)
      estate(5)   = pstate(5)/icm1 + 0.5_rk*sum(pstate(2:4)**2) * pstate(1)

      flux(iPoint,:) = atl_physFluxEuler(                                  &
        &                  state        = estate,                          &
        &                  isenCoeff    = isen_coeff,                      &
        &                  porosity     = euler%porosity,                  &
        &                  penalty_char = 0.5_rk                           &
        &                                * (material_left(matpoint,1)      &
        &                                   + material_right(matpoint,1)), &
        &                  U_o          = 0.5_rk                           &
        &                                * (material_left(matpoint,2)      &
        &                                   + material_right(matpoint,2))  )
    end do

  end subroutine atl_GodunovEuler


  !> Godunov flux for the 2D Euler equation.
  !!
  !! @todo Implement correct Riemann problem for the equation with porosity.
  !!       Right now the plain Euler problem is solved and the porosity and
  !!       penalty_char are just averaged centrally.
  subroutine atl_GodunovEuler2D(euler, state_left, state_right, &
    &                           material_left, material_right, nPoints, flux)
    ! ---------------------------------------------------------------------------
    class(atl_euler_type), intent(in) :: euler
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
    real(kind=rk) :: p_left, p_right
    real(kind=rk) :: v_left(2), v_right(2)
    real(kind=rk) :: estate(4), pstate(4)
    real(kind=rk) :: isen_coeff
    real(kind=rk) :: icm1
    real(kind=rk) :: dummy
    integer :: iPoint
    integer :: matpoint
    integer :: mm

    isen_coeff = euler%isen_coef
    icm1 = isen_coeff-1.0_rk
    mm = ubound(material_left,1)

    do iPoint = 1, nPoints
      matpoint = min(iPoint, mm)
      p_Left  = icm1 * ( state_left(iPoint,4)                   &
        &               - 0.5_rk*sum(state_left(iPoint,2:3)**2) &
        &                       / state_left(iPoint,1)          )

      p_Right = icm1 * ( state_right(iPoint,4)                   &
        &               - 0.5_rk*sum(state_right(iPoint,2:3)**2) &
        &                       / state_right(iPoint,1)          )

      v_left = state_left(iPoint,2:3) / state_left(iPoint,1)
      v_right = state_right(iPoint,2:3) / state_right(iPoint,1)

      call atl_ere_eval_onEdge(me = eqconsts, &
        &                      rho_left = state_left(iPoint,1), &
        &                      vn_left = v_left(1), &
        &                      v1_left = v_left(2), v2_left = 0.0_rk, &
        &                      p_left = p_left, &
        &                      rho_right = state_right(iPoint,1), &
        &                      vn_right = v_right(1), &
        &                      v1_right = v_right(2), v2_right = 0.0_rk, &
        &                      p_right = p_right, &
        &                      rho = pstate(1), vn = pstate(2), v1 = pstate(3), &
        &                      v2 = dummy, p = pstate(4))
      estate(1) = pstate(1)
      estate(2:3) = pstate(1)*pstate(2:3)
      estate(4) = pstate(4)/icm1 + 0.5_rk*sum(pstate(2:3)**2) * pstate(1)

      flux(iPoint,:) = atl_physFluxEuler_2d(                                &
        &                  state        = estate,                           &
        &                  isenCoeff    = isen_coeff,                       &
        &                  porosity     = euler%porosity,                   &
        &                  penalty_char = 0.5_rk                            &
        &                                * (material_left(matpoint,1)       &
        &                                   + material_right(matpoint,1)),  &
        &                  U_o          = 0.5_rk                            &
        &                                * (material_left(matpoint,2)       &
        &                                   + material_right(matpoint,2)) )
    end do

  end subroutine atl_GodunovEuler2D


  !> Godunov flux for the 1D Euler equation.
  !!
  !! @todo Implement correct Riemann problem for the equation with porosity.
  !!       Right now the plain Euler problem is solved and the porosity and
  !!       penalty_char are just averaged centrally.
  subroutine atl_GodunovEuler1D(euler, state_left, state_right, &
    &                           material_left, material_right, nPoints, flux)
    ! ---------------------------------------------------------------------------
    class(atl_euler_type), intent(in) :: euler
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
    real(kind=rk) :: p_left, p_right
    real(kind=rk) :: v_left, v_right
    real(kind=rk) :: estate(3), pstate(3)
    real(kind=rk) :: isen_coeff
    real(kind=rk) :: icm1
    real(kind=rk) :: dummy, dummy_2
    real(kind=rk) :: material(ubound(material_left,1))
    integer :: iPoint
    integer :: matpoint
    integer :: mm

    mm = ubound(material_left,1)

    material = 0.5_rk * (material_left(:,1) + material_right(:,1)) &
      &                        * (1.0_rk/euler%porosity - 1.0_rk)

    isen_coeff = euler%isen_coef
    icm1 = isen_coeff-1.0_rk

    do iPoint = 1, nPoints
      matpoint = min(iPoint, mm)
      p_Left  = icm1 * ( state_left(iPoint,3)              &
        &               - 0.5_rk*(state_left(iPoint,2)**2) &
        &                       / state_left(iPoint,1)     )

      p_Right = icm1 * ( state_right(iPoint,3)              &
        &               - 0.5_rk*(state_right(iPoint,2)**2) &
        &                       / state_right(iPoint,1)     )

      v_left = state_left(iPoint,2) / state_left(iPoint,1)
      v_right = state_right(iPoint,2) / state_right(iPoint,1)

      call atl_ere_eval_onEdge( me        = eqconsts,              &
        &                       rho_left  = state_left(iPoint,1),  &
        &                       vn_left   = v_left,                &
        &                       v1_left   = 0.0_rk,                &
        &                       v2_left   = 0.0_rk,                &
        &                       p_left    = p_left,                &
        &                       rho_right = state_right(iPoint,1), &
        &                       vn_right  = v_right,               &
        &                       v1_right  = 0.0_rk,                &
        &                       v2_right  = 0.0_rk,                &
        &                       p_right   = p_right,               &
        &                       rho       = pstate(1),             &
        &                       vn        = pstate(2),             &
        &                       v1        = dummy,                 &
        &                       v2        = dummy_2,               &
        &                       p         = pstate(3)              )
      estate(1) = pstate(1)
      estate(2) = pstate(1)*pstate(2)
      estate(3) = pstate(3)/icm1 + 0.5_rk*(pstate(2)**2) * pstate(1)

      flux(iPoint,:) = atl_physFluxEuler_1D(                          &
        & state           = estate,                                   &
        & isenCoeff       = isen_coeff,                               &
        & U_o             = 0.5_rk * (material_left(matpoint,2)       &
        &                              + material_right(matpoint,2)), &
        & penalty_scaling = material(matpoint)                        )
    end do

  end subroutine atl_GodunovEuler1D


  subroutine atl_flux_initGodunov(isen_coeff)
    real(kind=rk) :: isen_coeff

    call atl_ere_init_consts(me = eqconsts, gam = isen_coeff)
  end subroutine atl_flux_initGodunov

end module atl_GodunovFlux_module
