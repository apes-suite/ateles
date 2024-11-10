! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> author: Jens Zudrop
!! Collects all functions related to the physical fluxes of the compressible Navier-Stokes equations.
module atl_physFluxFilNvrStk_module
  ! Treelm modules
  use env_module,            only: rk
  use tem_aux_module,        only: tem_abort
  use atl_eqn_nvrstk_module, only: atl_Navier_stokes_rans_type
  use atl_eqn_filNvrStk_var_module, &
    &                        only: atl_get_pointwise_velocity_gradient_2D,  &
    &                              atl_get_pointwise_visc_stress_tensor_2D, &
    &                              atl_get_lower_bound_turb_disscipation

  implicit none
  private


  public :: atl_viscPhysFluxRans
  public :: atl_PhysFluxRans
  public :: atl_PhysFluxRans_2d
  public :: atl_viscPhysFluxRans_2d,                      &
          & atl_mult_nu11_Rans_2d, atl_mult_nu21_Rans_2d, &
          & atl_mult_nu12_Rans_2d, atl_mult_nu22_Rans_2d

contains

  !> Physical flux calculation along x direction for the
  !  filtered Navier Stokes equation.
  function atl_physFluxRans(state, isenCoeff, penalty_char, porosity)       &
    &                       result(physFlux)
    ! ---------------------------------------------------------------------------
    !> The state in nodal space. Dimension is the number of vars, i.e. 5 for Euler
    !  and here =7 for renolds averaged navier stokes
    real(kind=rk), intent(in) :: state(:)
    !> Adiabatice index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff
    !> The value of the characteristic function (stemming from penalization)
    real(kind=rk), intent(in) :: penalty_char
    !> The porosity at the current point
    real(kind=rk), intent(in) :: porosity
    !> The physical flux along the x axis for all variables
    real(kind=rk) :: physFlux(7)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: pressure, velocity(1:3), pressure_eff
    ! ---------------------------------------------------------------------------

    ! calculate pressure
    pressure = (isenCoeff-1.0_rk) * ( &
                 & state(5) - 0.5_rk*(sum(state(2:4)**2))/state(1) &
                           - state(6)           )
    pressure_eff = pressure + 2.0_rk * state(6) / 3.0_rk

    !> @todo JZ: here, we divide by a polynomial, we should be careful! We are leaving
    !! the polynomial space here!
    velocity(1:3) = state(2:4)/state(1)

    ! calculate the nonlinear term for different varibales now.
    ! ... density
    physFlux(1) = (1.0_rk + ((1.0_rk/porosity)-1.0_rk)*penalty_char) * state(2)
    ! ... x-velocity
    physFlux(2) = pressure_eff + state(2)*velocity(1)
    ! ... y-velocity
    physFlux(3) = state(2)*velocity(2)
    ! ... z-velocity
    physFlux(4) = state(2)*velocity(3)
    ! ... total energy
    physFlux(5) = velocity(1) * ( state(5) + pressure_eff )
    ! ... turbulent KE
    physFlux(6) = state(2)*state(6)/ state(1)
    ! ... spec dissipation rate
    physFlux(7) = state(2)*state(7)/ state(1)

  end function atl_physFluxRans


  function atl_physFluxRans_2d(state, isenCoeff, penalty_char, porosity)       &
    &                       result(physFlux)
    ! ---------------------------------------------------------------------------
    !> The state in nodal space. Dimension is the number of vars, i.e. 5 for Euler
    !  and here =7 for renolds averaged navier stokes
    real(kind=rk), intent(in) :: state(:)
    !> Adiabatice index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff
    !> The value of the characteristic function (stemming from penalization)
    real(kind=rk), intent(in) :: penalty_char
    !> The porosity at the current point
    real(kind=rk), intent(in) :: porosity
    !> The physical flux along the x axis for all variables
    real(kind=rk) :: physFlux(6)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: pressure, velocity(1:2), pressure_eff
    ! ---------------------------------------------------------------------------

    ! calculate pressure
    pressure = (isenCoeff-1.0_rk) * ( &
                 & state(4) - 0.5_rk*(sum(state(2:3)**2))/state(1) &
                           - state(5)           )
    pressure_eff = pressure + 2.0_rk * state(5) / 3.0_rk

    !> @todo JZ: here, we divide by a polynomial, we should be careful! We are leaving
    !! the polynomial space here!
    velocity(1:2) = state(2:3)/state(1)

    ! calculate the nonlinear term for different varibales now.
    ! ... density
    physFlux(1) = (1.0_rk + ((1.0_rk/porosity)-1.0_rk)*penalty_char) * state(2)
    ! ... x-velocity
    physFlux(2) = pressure_eff + state(2)*velocity(1)
    ! ... y-velocity
    physFlux(3) = state(2)*velocity(2)
    ! ... total energy
    physFlux(4) = velocity(1) * ( state(4) + pressure_eff )
    ! ... turbulent KE
    physFlux(5) = velocity(1)*state(5)
    ! ... spec dissipation rate
    physFlux(6) = velocity(1)*state(6)

  end function atl_physFluxRans_2d


  !> Physical flux calculation along x direction for Euler equation.
  function atl_viscPhysFluxRans(state, state_gradient, isenCoeff, mu, &
                                   & lambda, thermCond, heatCap ) result(physFlux)
    ! ------------------------------------------------------------------------------------
    !> The state in nodal space. Dimension is the number of vars, i.e. 5 for Navier-Stokes.
    real(kind=rk), intent(in) :: state(:)
    !> The state in nodal space. First dimension is the number of vars, i.e. 5 for Navier-Stokes.
    !! Second dimension is the dimension, e.g. 3 in two dimensions.
    real(kind=rk), intent(in) :: state_gradient(:,:)
    !> Adiabatice index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The physical flux along the x axis for all variables
    real(kind=rk) :: physFlux(7)
    ! -------------------------------------------------------------------------------
    real(kind=rk) :: velocity(1:3)
    ! The turbulent viscosity (dynamic)
    real(kind=rk) :: mu_turb
    ! The turbulent viscosity
    real(kind=rk) :: lam_turb
    ! The effective viscosity (dynamic)
    real(kind=rk) :: mu_eff
    ! The effective turbulent viscosity
    real(kind=rk) :: lam_eff
    real(kind=rk) :: energy_coeff
    real(kind=rk) :: sig_k
    real(kind=rk) :: sig_omg
    real(kind=rk) :: turbulent_prandtl_number
    ! -------------------------------------------------------------------------------

    velocity(1:3) = state(2:4)/state(1)
    !mu_turb = rans_params%c_mu*state(5)/exp(state(6)/state(1))

    ! The code above was deactivated without adding a proper initialization
    ! for the variables. I'm letting the code  abort when it is executed to
    ! inform the user that it is not yet usable.
    call tem_abort( "This routine doesn't use properly initialized variables" )

    lam_turb = 2.0*mu_turb/3.0
    mu_eff = mu + mu_turb
    lam_eff = lambda + lam_turb
    turbulent_prandtl_number = 0.9_rk
    sig_k = 0.5_rk
    sig_omg = 0.5_rk
    energy_coeff = thermCond/(heatCap*state(1))      &
    &                 + mu_turb*isenCoeff/(turbulent_prandtl_number*state(1))

    ! Viscous flux for density
    physFlux(1) = 0.0_rk

    ! Viscous flux for momentum in x
    physFlux(2) =                                                                 &
      ! (nu_{1,1})_{2,1} (\nabla u)_{1,1} : k = 1, i = 1
      & ((-2.0_rk * mu_eff + lam_eff) / state(1))*velocity(1)*state_gradient(1,1) &
      ! (nu_{1,2})_{2,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + (lam_eff / state(1))*velocity(2)*state_gradient(1,2)                    &
      ! (nu_{1,3})_{2,1} (\nabla u)_{1,3} : k = 3, i = 1
      & + (lam_eff / state(1))*velocity(3)*state_gradient(1,3)                    &
      ! (nu_{1,1})_{2,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + ((2.0_rk * mu_eff - lam_eff) / state(1))*state_gradient(2,1)            &
      ! (nu_{1,2})_{2,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + 0.0_rk                                                                  &
      !(nu_{1,3})_{2,2} (\nabla u)_{2,3} : k = 3, i = 2
      & + 0.0_rk                                                                  &
      ! (nu_{1,1})_{2,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + 0.0_rk                                                                  &
      ! (nu_{1,2})_{2,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + (-lam_eff/state(1))*state_gradient(3,2)                                 &
      ! (nu_{1,3})_{2,3} (\nabla u)_{3,3} : k = 2, i = 3
      & + 0.0_rk                                                                  &
      ! (nu_{1,1})_{2,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + 0.0_rk                                                                  &
      ! (nu_{1,2})_{2,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk                                                                  &
      ! (nu_{1,3})_{2,4} (\nabla u)_{4,3} : k = 3, i = 4
      & + (-lam_eff/state(1))*state_gradient(4,3)                                 &
      ! (nu_{1,1})_{2,5} (\nabla u)_{5,1} : k = 1, i = 5
      & + 0.0_rk                                                                  &
      ! (nu_{1,2})_{2,5} (\nabla u)_{5,2} : k = 2, i = 5
      & + 0.0_rk                                                                  &
      ! (nu_{1,3})_{2,5} (\nabla u)_{5,3} : k = 3, i = 5
      & + 0.0_rk                                                                  &
      ! (nu_{1,1})_{2,6} (\nabla u)_{6,1} : k = 1, i = 6
      & + 0.0_rk                                                                  &
      ! (nu_{1,2})_{2,6} (\nabla u)_{6,2} : k = 2, i = 6
      & + 0.0_rk                                                                  &
      ! (nu_{1,3})_{2,6} (\nabla u)_{6,3} : k = 3, i = 6
      & + 0.0_rk                                                                  &
      ! (nu_{1,1})_{2,7} (\nabla u)_{7,1} : k = 1, i = 7
      & + 0.0_rk                                                                  &
      ! (nu_{1,2})_{2,7} (\nabla u)_{7,2} : k = 2, i = 7
      & + 0.0_rk                                                                  &
      ! (nu_{1,3})_{2,7} (\nabla u)_{7,3} : k = 3, i = 7
      & + 0.0_rk

    ! Viscous flux for momentum in y
    physFlux(3) =                                                              &
      ! (nu_{1,1})_{3,1} (\nabla u)_{1,1} : k = 1, i = 1
      &   (-mu_eff / state(1))*velocity(2)*state_gradient(1,1)                 &
      ! (nu_{1,2})_{3,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + (-mu_eff / state(1))*velocity(1)*state_gradient(1,2)                 &
      ! (nu_{1,3})_{3,1} (\nabla u)_{1,3} : k = 3, i = 1
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{3,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{3,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + (mu_eff/state(1))*state_gradient(2,2)                                &
      ! (nu_{1,3})_{3,2} (\nabla u)_{2,3} : k = 3, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{3,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + (mu_eff/state(1))*state_gradient(3,1)                                &
      ! (nu_{1,2})_{3,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{3,3} (\nabla u)_{3,3} : k = 3, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{3,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{3,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{3,4} (\nabla u)_{4,3} : k = 3, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{3,5} (\nabla u)_{5,1} : k = 1, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{3,5} (\nabla u)_{5,2} : k = 2, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{3,5} (\nabla u)_{5,3} : k = 3, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{3,6} (\nabla u)_{6,1} : k = 1, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{3,6} (\nabla u)_{6,2} : k = 2, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{3,6} (\nabla u)_{6,3} : k = 3, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{3,7} (\nabla u)_{7,1} : k = 1, i = 7
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{3,7} (\nabla u)_{7,2} : k = 2, i = 7
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{3,7} (\nabla u)_{7,3} : k = 3, i = 7
      & + 0.0_rk


    ! Viscous flux for momentum in z
    physFlux(4) =                                                              &
      ! (nu_{1,1})_{4,1} (\nabla u)_{1,1} : k = 1, i = 1
      &   (-mu_eff / state(1))*velocity(3)*state_gradient(1,1)                 &
      ! (nu_{1,2})_{4,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{4,1} (\nabla u)_{1,3} : k = 3, i = 1
      & + (-mu_eff / state(1))*velocity(1)*state_gradient(1,3)                 &
      ! (nu_{1,1})_{4,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{4,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{4,2} (\nabla u)_{2,3} : k = 3, i = 2
      & + (mu_eff/state(1))*state_gradient(2,3)                                &
      ! (nu_{1,1})_{4,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{4,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{4,3} (\nabla u)_{3,3} : k = 3, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{4,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + (mu_eff/state(1))*state_gradient(4,1)                                &
      ! (nu_{1,2})_{4,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{4,4} (\nabla u)_{4,3} : k = 3, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{4,5} (\nabla u)_{5,1} : k = 1, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{4,5} (\nabla u)_{5,2} : k = 2, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{4,5} (\nabla u)_{5,3} : k = 3, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{4,6} (\nabla u)_{6,1} : k = 1, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{4,6} (\nabla u)_{6,2} : k = 2, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{4,6} (\nabla u)_{6,3} : k = 3, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{4,7} (\nabla u)_{7,1} : k = 1, i = 7
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{4,7} (\nabla u)_{7,2} : k = 2, i = 7
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{4,7} (\nabla u)_{7,3} : k = 3, i = 7
      & + 0.0_rk

    ! Viscous flux for Energy
    physFlux(5) = &
      ! (nu_{1,1})_{5,1} (\nabla u)_{1,1} : k = 1, i = 1
      & ( ((-2.0_rk * mu_eff + lam_eff) / state(1))*(velocity(1)**2.0_rk)               &
      & +(-mu_eff/state(1))*(velocity(2)**2.0_rk)                                       &
      & +(-mu_eff/state(1))*(velocity(3)**2.0_rk)                                       &
      & -energy_coeff*(state(5)/state(1)-sum(velocity(:)**2.0_rk) - state(6)/state(1))) &
      &  * state_gradient(1,1)                                                          &
      ! (nu_{1,2})_{5,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + ((-mu_eff+lam_eff)/state(1))*velocity(1)*velocity(2)*state_gradient(1,2)      &
      ! (nu_{1,3})_{5,1} (\nabla u)_{1,3} : k = 3, i = 1
      & + ((-mu_eff+lam_eff)/state(1))*velocity(1)*velocity(3)*state_gradient(1,3)      &
      ! (nu_{1,1})_{5,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + ((2.0_rk*mu_eff-lam_eff)/state(1) - energy_coeff)*velocity(1)                 &
      &  * state_gradient(2,1)                                                          &
      ! (nu_{1,2})_{5,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + (mu_eff/state(1))*velocity(2)*state_gradient(2,2)                             &
      ! (nu_{1,3})_{5,2} (\nabla u)_{2,3} : k = 3, i = 2
      & + (mu_eff/state(1))*velocity(3)*state_gradient(2,3)                             &
      ! (nu_{1,1})_{5,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + ((mu_eff/state(1))- energy_coeff)*velocity(2)*state_gradient(3,1)             &
      ! (nu_{1,2})_{5,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + (-lam_eff/state(1))*velocity(1)*state_gradient(3,2)                           &
      ! (nu_{1,3})_{5,3} (\nabla u)_{3,3} : k = 3, i = 3
      & + 0.0_rk                                                                        &
      ! (nu_{1,1})_{5,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + ((mu_eff/state(1))- energy_coeff)*velocity(3)*state_gradient(4,1)             &
      ! (nu_{1,2})_{5,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk                                                                        &
      ! (nu_{1,3})_{5,4} (\nabla u)_{4,3} : k = 3, i = 4
      & + (-lam_eff/state(1))*velocity(1)*state_gradient(4,3)                           &
      ! (nu_{1,1})_{5,5} (\nabla u)_{5,1} : k = 1, i = 5
      & + energy_coeff*state_gradient(5,1)                                              &
      ! (nu_{1,2})_{5,5} (\nabla u)_{5,2} : k = 2, i = 5
      & + 0.0_rk                                                                        &
      ! (nu_{1,3})_{5,5} (\nabla u)_{5,3} : k = 3, i = 5
      & + 0.0_rk                                                                        &
      ! (nu_{1,1})_{5,6} (\nabla u)_{6,1} : k = 1, i = 6
      & + ((mu + sig_k*mu_turb)/state(1) - energy_coeff)*state_gradient(6,1)            &
      ! (nu_{1,2})_{5,6} (\nabla u)_{6,2} : k = 2, i = 6
      & + 0.0_rk                                                                        &
      ! (nu_{1,3})_{5,6} (\nabla u)_{6,3} : k = 3, i = 6
      & + 0.0_rk                                                                        &
      ! (nu_{1,1})_{5,7} (\nabla u)_{7,1} : k = 1, i = 7
      & + 0.0_rk                                                                        &
      ! (nu_{1,2})_{5,7} (\nabla u)_{7,2} : k = 2, i = 7
      & + 0.0_rk                                                                        &
      ! (nu_{1,3})_{5,7} (\nabla u)_{7,3} : k = 3, i = 7
      & + 0.0_rk


    ! Viscous flux for turbulent KE
    physFlux(6) =                                                              &
      ! (nu_{1,1})_{6,1} (\nabla u)_{1,1} : k = 1, i = 1
      & ( -(mu + sig_k*mu_turb)/state(1))*(state(6)/state(1))                  &
      &                                  *state_gradient(1,1)                  &
      ! (nu_{1,2})_{6,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{6,1} (\nabla u)_{1,3} : k = 3, i = 1
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{6,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{6,2} (\nabla u)_{2,3} : k = 3, i = 2
      & + (mu_eff/state(1))*state_gradient(2,3)                                &
      ! (nu_{1,1})_{6,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{6,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{6,3} (\nabla u)_{3,3} : k = 3, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + (mu_eff/state(1))*state_gradient(4,1)                                &
      ! (nu_{1,2})_{6,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{6,4} (\nabla u)_{4,3} : k = 3, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,5} (\nabla u)_{5,1} : k = 1, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{6,5} (\nabla u)_{5,2} : k = 2, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{6,5} (\nabla u)_{5,3} : k = 3, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,6} (\nabla u)_{6,1} : k = 1, i = 6
      & + (mu + sig_k*mu_turb)/state(1)*state_gradient(6,1)                    &
      ! (nu_{1,2})_{6,6} (\nabla u)_{6,2} : k = 2, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{6,6} (\nabla u)_{6,3} : k = 3, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,7} (\nabla u)_{7,1} : k = 1, i = 7
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{6,7} (\nabla u)_{7,2} : k = 2, i = 7
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{6,7} (\nabla u)_{7,3} : k = 3, i = 7
      & + 0.0_rk

    ! Viscous flux for specific Dissipation rate (\Omega)
    physFlux(7) =                                                              &
      ! (nu_{1,1})_{7,1} (\nabla u)_{1,1} : k = 1, i = 1
      & ( -(mu + sig_k*mu_turb)/state(1))*(state(7)/state(1))                  &
      &                                  *state_gradient(1,1)                  &
      ! (nu_{1,2})_{7,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{7,1} (\nabla u)_{1,3} : k = 3, i = 1
      & + 0.0_rk                &
      ! (nu_{1,1})_{7,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{7,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{7,2} (\nabla u)_{2,3} : k = 3, i = 2
      & + (mu_eff/state(1))*state_gradient(2,3)                                &
      ! (nu_{1,1})_{7,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{7,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{7,3} (\nabla u)_{3,3} : k = 3, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{7,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + (mu_eff/state(1))*state_gradient(4,1)                                &
      ! (nu_{1,2})_{7,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{7,4} (\nabla u)_{4,3} : k = 3, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{7,5} (\nabla u)_{5,1} : k = 1, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{7,5} (\nabla u)_{5,2} : k = 2, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{7,5} (\nabla u)_{5,3} : k = 3, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{7,6} (\nabla u)_{6,1} : k = 1, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{7,6} (\nabla u)_{6,2} : k = 2, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{7,6} (\nabla u)_{6,3} : k = 3, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{7,7} (\nabla u)_{7,1} : k = 1, i = 7
      & + ((mu + sig_omg*mu_turb)/state(1))*state_gradient(7,1)                &
      ! (nu_{1,2})_{7,7} (\nabla u)_{7,2} : k = 2, i = 7
      & + 0.0_rk                                                               &
      ! (nu_{1,3})_{7,7} (\nabla u)_{7,3} : k = 3, i = 7
      & + 0.0_rk


  end function atl_viscPhysFluxRans


  function atl_viscPhysFluxRans_2d(state, state_gradient, isenCoeff, mu, &
                                   & lambda, thermCond, rans_params, heatCap ) result(physFlux)
    ! ------------------------------------------------------------------------------------
    !> The state in nodal space. Dimension is the number of vars, i.e. 4 for Navier-Stokes.
    real(kind=rk), intent(in) :: state(:)
    !> The state in nodal space. First dimension is the number of vars, i.e. 4 for Navier-Stokes.
    !! Second dimension is the dimension, e.g. 2 in two dimensions.
    real(kind=rk), intent(in) :: state_gradient(:,:)
    !> Adiabatice index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The physical flux along the x axis for all variables
    real(kind=rk) :: physFlux(6)
    !> The constants for the Rans eqn
    type(atl_Navier_stokes_rans_type), intent(in) :: rans_params
    ! -----------------------------------------------------------------------------
    real(kind=rk) :: velocity(1:2)
    ! The turbulent viscosity (dynamic)
    real(kind=rk) :: mu_turb
    ! The turbulent viscosity
    real(kind=rk) :: lam_turb
    ! The effective viscosity (dynamic)
    real(kind=rk) :: mu_eff
    ! The effective turbulent viscosity
    real(kind=rk) :: lam_eff
    real(kind=rk) :: energy_coeff
    real(kind=rk) :: k_eff
    real(kind=rk) :: turbPrNum, limited_eddy_visc, turb_coeff1, turb_coeff2
    real(kind=rk) :: velGrad(2,2), ViscStressTensor(2,2),  omega_r
    ! -----------------------------------------------------------------------------
    mu_turb = rans_params%c_mu*state(5)/exp(state(6)/state(1))

    lam_turb = 2.0*mu_turb/3.0
    mu_eff = mu + mu_turb
    lam_eff = lambda + lam_turb
    turbPrNum = rans_params%turb_prandtl_num

    k_eff = max(state(6),0.0_rk)
    energy_coeff = thermCond/(heatCap*state(1))               &
    &                 + mu_turb*isenCoeff/(turbPrNum*state(1) )
    velocity(1:2) = state(2:3)/state(1)

    ! Now to calculate the coefficient \omega_r: We need to do the following
    ! Step-1: Get the velocity gradient from the state gradient present
    call atl_get_pointwise_velocity_gradient_2D ( state_gradient, &
     &                           state, velGrad                   )

    ! Step-2: Calculate the viscous stress tensor
    call atl_get_pointwise_visc_stress_tensor_2D(         &
      &          velGrad = velGrad,                       &
      &          S       = ViscStressTensor               )

    ! Step-3: Get the value of \omega_r
    call atl_get_lower_bound_turb_disscipation(           &
      &          S       = ViscStressTensor,              &
      &        c_mu      = rans_params%c_mu,              &
      &        omega     = state(6)/ state(1),            &
      &        omega_r   = omega_r                        )


    limited_eddy_visc = rans_params%alpha*k_eff*exp(-omega_r)
    turb_coeff1 = (mu + rans_params%sig_k*limited_eddy_visc)/ state(1)
    turb_coeff2 = (mu + rans_params%sig_omg*limited_eddy_visc)/ state(1)

    ! Viscous flux for density
    physFlux(1) = 0.0_rk

    ! Viscous flux for momentum in x
    physFlux(2) =                                                               &
      ! (nu_{1,1})_{2,1} (\nabla u)_{1,1} : k = 1, i = 1
      & ((-2.0_rk * mu_eff + lam_eff) / state(1))*velocity(1)*state_gradient(1,1) &
      ! (nu_{1,2})_{2,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + (lam_eff / state(1))*velocity(2)*state_gradient(1,2)                  &
      ! (nu_{1,1})_{2,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + ((2.0_rk * mu_eff - lam_eff) / state(1))*state_gradient(2,1)         &
      ! (nu_{1,2})_{2,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{2,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{2,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + (-lam_eff/state(1))*state_gradient(3,2)                              &
      ! (nu_{1,1})_{2,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{2,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk

    ! Viscous flux for momentum in y
    physFlux(3) =                                                              &
      ! (nu_{1,1})_{3,1} (\nabla u)_{1,1} : k = 1, i = 1
      &   (-mu_eff / state(1))*velocity(2)*state_gradient(1,1)                 &
      ! (nu_{1,2})_{3,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + (-mu_eff / state(1))*velocity(1)*state_gradient(1,2)                 &
      ! (nu_{1,1})_{3,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{3,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + (mu_eff/state(1))*state_gradient(2,2)                                &
      ! (nu_{1,1})_{3,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + (mu_eff/state(1))*state_gradient(3,1)                                &
      ! (nu_{1,2})_{3,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{3,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{3,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk

    ! Viscous flux for Energy
    physFlux(4) =                                                                   &
     ! (nu_{1,1})_{4,1} (\nabla u)_{1,1} : k = 1, i = 1
     & ( ((-2.0_rk * mu_eff + lam_eff) / state(1))*(velocity(1)**2.0_rk)            &
     & + (-mu_eff/state(1))*(velocity(2)**2.0_rk)                                   &
     & -energy_coeff*(state(4)/state(1)-sum(velocity(:)**2.0_rk) - state(5)/state(1)) &
     & - turb_coeff1*(state(5)/state(1)) )                                          &
     &  * state_gradient(1,1)                                                       &
     ! (nu_{1,2})_{4,1} (\nabla u)_{1,2} : k = 2, i = 1
     & + ((-mu_eff+lam_eff)/state(1))*velocity(1)*velocity(2)*state_gradient(1,2)   &
     ! (nu_{1,1})_{4,2} (\nabla u)_{2,1} : k = 1, i = 2
     & + ((2.0_rk*mu_eff-lam_eff)/state(1) - energy_coeff)*velocity(1)              &
     &  * state_gradient(2,1)                                                       &
     ! (nu_{1,2})_{4,2} (\nabla u)_{2,2} : k = 2, i = 2
     & + (mu_eff/state(1))*velocity(2)*state_gradient(2,2)                          &
     ! (nu_{1,1})_{4,3} (\nabla u)_{3,1} : k = 1, i = 3
     & + ((mu_eff/state(1)) - energy_coeff)*velocity(2)*state_gradient(3,1)    &
     ! (nu_{1,2})_{4,3} (\nabla u)_{3,2} : k = 2, i = 3
     & + (-lam_eff/state(1))*velocity(1)*state_gradient(3,2)                  &
     ! (nu_{1,1})_{4,4} (\nabla u)_{4,1} : k = 1, i = 4
     & + energy_coeff*state_gradient(4,1)                                     &
     ! (nu_{1,2})_{4,4} (\nabla u)_{4,2} : k = 2, i = 4
     & + 0.0_rk                                                               &
     ! (nu_{1,1})_{4,5} (\nabla u)_{5,1} : k = 1, i = 5
     & + (-energy_coeff + turb_coeff1)*state_gradient(5,1)                    &
     ! (nu_{1,2})_{4,5} (\nabla u)_{5,2} : k = 2, i = 5
     & + 0.0_rk                                                               &
     ! (nu_{1,1})_{4,6} (\nabla u)_{6,1} : k = 1, i = 6
     & + 0.0_rk                                                               &
     ! (nu_{1,2})_{4,6} (\nabla u)_{6,2} : k = 2, i = 6
     & + 0.0_rk

    physFlux(5) =                                                              &
      ! (nu_{1,1})_{5,1} (\nabla u)_{1,1} : k = 1, i = 1
      & -turb_coeff1*(state(5)/state(1)) *state_gradient(1,1)                  &
      ! (nu_{1,2})_{5,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{5,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{5,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{5,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{5,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{5,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{5,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{5,5} (\nabla u)_{5,1} : k = 1, i = 5
      & + turb_coeff1*state_gradient(5,1)                                      &
      ! (nu_{1,2})_{5,5} (\nabla u)_{5,2} : k = 2, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{5,6} (\nabla u)_{6,1} : k = 1, i = 6
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{5,6} (\nabla u)_{6,2} : k = 2, i = 6
      & + 0.0_rk

    physFlux(6) =                                                              &
      ! (nu_{1,1})_{6,1} (\nabla u)_{1,1} : k = 1, i = 1
      & ( -turb_coeff2*state(6)/state(1) )*state_gradient(1,1)                 &
      ! (nu_{1,2})_{6,1} (\nabla u)_{1,2} : k = 2, i = 1
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,2} (\nabla u)_{2,1} : k = 1, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{6,2} (\nabla u)_{2,2} : k = 2, i = 2
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,3} (\nabla u)_{3,1} : k = 1, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{6,3} (\nabla u)_{3,2} : k = 2, i = 3
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,4} (\nabla u)_{4,1} : k = 1, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{6,4} (\nabla u)_{4,2} : k = 2, i = 4
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,5} (\nabla u)_{5,1} : k = 1, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,2})_{6,5} (\nabla u)_{5,2} : k = 2, i = 5
      & + 0.0_rk                                                               &
      ! (nu_{1,1})_{6,6} (\nabla u)_{6,1} : k = 1, i = 6
      & + turb_coeff2*state_gradient(6,1)                                      &
      ! (nu_{1,2})_{6,6} (\nabla u)_{6,2} : k = 2, i = 6
      & + 0.0_rk

  end function atl_viscPhysFluxRans_2d


  ! Multiplies the viscous flux matrux nu_11 with a given vector
  function atl_mult_nu11_Rans_2d( state, velocity, inVec,                     &
    &                isenCoeff, mu, lambda, thermCond, rans_params, heatCap)  &
    & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The velocity
    real(kind=rk), intent(in) :: velocity(2)
    !> The state array
    real(kind=rk), intent(in) :: state(6)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(6)
    !> Adiabatice index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The constants for the Rans eqn
    type(atl_Navier_stokes_rans_type), intent(in) :: rans_params
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(6)
    ! -------------------------------------------------------------------------!
    !> mu_turb
    real(kind=rk) :: mu_turb, mu_eff, lam_eff, limited_eddy_visc, energy_coeff
    real(kind=rk) :: turb_coeff1, turb_coeff2
    ! -------------------------------------------------------------------------!
    mu_turb = rans_params%c_mu*state(5)/exp(state(6)/state(1))

    mu_eff = mu + mu_turb
    lam_eff = lambda + 2.0*mu_turb/3.0

    ! Stuff neeeded
    energy_coeff = thermCond/(heatCap*state(1))                               &
    &             + mu_turb*isenCoeff/(rans_params%turb_prandtl_num*state(1))

!NA!    ! Now to calculate the coefficient \omega_r: We need to do the following
!NA!    ! Step-1: Get the velocity gradient from the state gradient present
!NA!    call atl_get_pointwise_velocity_gradient_2D ( inVec, &
!NA!     &                           state, velGrad               )
!NA!
!NA!    ! Step-2: Calculate the viscous stress tensor
!NA!    call atl_get_pointwise_visc_stress_tensor_2D(             &
!NA!      &          velGrad = velGrad,                       &
!NA!      &          S       = ViscStressTensor               )
!NA!
!NA!    ! Step-3: Get the value of \omega_r
!NA!    call atl_get_lower_bound_turb_disscipation(                          &
!NA!      &          S       = ViscStressTensor,                         &
!NA!      &        c_mu      = rans_params%c_mu,                         &
!NA!      &        omega     = state(6)/ state(1),                       &
!NA!      &        omega_r   = omega_r                                   )


    ! @todo : NA : Calculate limited_eddy_visc Appropriately
!NA!    limited_eddy_visc = rans_params%alpha*max(state(5),0.0)*exp(-omega_r)
    limited_eddy_visc = rans_params%alpha*max(state(5),0.0_rk)              &
      &                                  *exp(state(6)/state(1))
    turb_coeff1 = (mu + rans_params%sig_k*limited_eddy_visc)/ state(1)
    turb_coeff2 = (mu + rans_params%sig_omg*limited_eddy_visc)/ state(1)

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = ((-2.0_rk*mu_eff + lam_eff)/state(1))*velocity(1) * inVec(1)   &
            & + ((2.0_rk*mu_eff - lam_eff)/state(1)) * inVec(2)

    ! Third row
    outVec(3) = (-mu_eff/state(1)) * velocity(2) * inVec(1)                    &
            & + (mu_eff/state(1)) * inVec(3)

    ! Fourth row
    outVec(4) = (                                                              &
      &         ((-2.0_rk*mu_eff + lam_eff)/state(1))*velocity(1)*velocity(1)  &
      &           + (-mu_eff/state(1)) * velocity(2)*velocity(2)               &
      &           - energy_coeff*(state(4)/state(1) - sum(velocity(:)**2)      &
      &                                      - state(5)/state(1) )             &
      &           - turb_coeff1*(state(5)/state(1))                )* inVec(1) &
      & + (                                                                    &
      &     ((2.0_rk*mu_eff-lam_eff)/state(1) - energy_coeff)*velocity(1)      &
      &   ) * inVec(2)                                                         &
      & + (  (mu_eff/state(1) - energy_coeff)*velocity(2))* inVec(3)           &
      & + (  energy_coeff  ) * inVec(4)                                        &
      & + (  -energy_coeff + turb_coeff1 )* inVec(5)

    ! Fifth Row
    outVec(5) = (                                                              &
      &     -turb_coeff1*(state(5)/state(1)) )* inVec(1)                        &
      &     + turb_coeff1*inVec(5)

    ! Sixth Row
    outVec(6) = ( -turb_coeff2*state(6)/state(1)  )* inVec(1)                   &
      &         +  ( turb_coeff2 )* inVec(6)

  end function atl_mult_nu11_Rans_2d

  ! Multiplies the viscous flux matrux nu_21 with a given vector
  function atl_mult_nu21_Rans_2d( state, velocity, inVec, mu, lambda, &
    &                             rans_params)  &
    & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The velocity
    real(kind=rk), intent(in) :: velocity(2)
    !> The state array
    real(kind=rk), intent(in) :: state(6)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(6)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The constants for the Rans eqn
    type(atl_Navier_stokes_rans_type), intent(in) :: rans_params
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(6)
    ! -------------------------------------------------------------------------!
    !> mu_turb
    real(kind=rk) :: mu_turb, mu_eff, lam_eff
    ! -------------------------------------------------------------------------!
    mu_turb = rans_params%c_mu*state(5)/exp(state(6)/state(1))

    mu_eff = mu + mu_turb
    lam_eff = lambda + 2.0*mu_turb/3.0

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (-mu_eff/state(1)) * velocity(2) * inVec(1) &
            & + (mu_eff/state(1)) * inVec(3)

    ! Third row
    outVec(3) = (lam_eff/state(1)) * velocity(1) * inVec(1) &
            & + (-lam_eff/state(1)) * inVec(2)

    ! Fourth row
    outVec(4) = (                                                                           &
            &      ((-mu_eff+lam_eff)/state(1))*velocity(1)*velocity(2)                     &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &      (-lam_eff/state(1))*velocity(2)                                          &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &      (mu_eff/state(1))*velocity(1)                                            &
            &   ) * inVec(3)

    outVec(5) = 0.0_rk


    outVec(6) = 0.0_rk

  end function atl_mult_nu21_Rans_2d


  ! Multiplies the viscous flux matrux nu_12 with a given vector
  function atl_mult_nu12_Rans_2d( state, velocity, inVec, mu, lambda, &
    &                             rans_params )  &
    & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The velocity
    real(kind=rk), intent(in) :: velocity(2)
    !> The state array
    real(kind=rk), intent(in) :: state(6)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(6)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The constants for the Rans eqn
    type(atl_Navier_stokes_rans_type), intent(in) :: rans_params
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(6)
    ! -------------------------------------------------------------------------!
    !> mu_turb
    real(kind=rk) :: mu_turb, mu_eff, lam_eff
    ! -------------------------------------------------------------------------!
    mu_turb = rans_params%c_mu*state(5)/exp(state(6)/state(1))

    mu_eff = mu + mu_turb
    lam_eff = lambda + 2.0*mu_turb/3.0

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (lam_eff/state(1)) * velocity(2) * inVec(1) &
      &         + (-lam_eff/state(1)) * inVec(3)

    ! Third row
    outVec(3) = (-mu_eff/state(1)) * velocity(1) * inVec(1) &
       &        + (mu_eff/state(1)) * inVec(2)

    ! Fourth row
    outVec(4) = (                                                                           &
            &      ((-mu_eff+lam_eff)/state(1))*velocity(1)*velocity(2)                     &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &      (mu_eff/state(1))*velocity(2)                                            &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &      (-lam_eff/state(1))*velocity(1)                                          &
            &   ) * inVec(3)


     outVec(5) = 0.0_rk
     outVec(6) = 0.0_rk

  end function atl_mult_nu12_Rans_2d

  ! Multiplies the viscous flux matrux nu_22 with a given vector
  function atl_mult_nu22_Rans_2d( state, velocity, inVec,                     &
    &                isenCoeff, mu, lambda, thermCond, rans_params, heatCap)  &
    & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The velocity
    real(kind=rk), intent(in) :: velocity(2)
    !> The state array
    real(kind=rk), intent(in) :: state(6)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(6)
    !> Adiabatice index, also known as isentropic expansion factor.
    real(kind=rk), intent(in) :: isenCoeff
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The constants for the Rans eqn
    type(atl_Navier_stokes_rans_type), intent(in) :: rans_params
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(6)
    ! -------------------------------------------------------------------------!
    !> mu_turb
    real(kind=rk) :: mu_turb, mu_eff, lam_eff, limited_eddy_visc, energy_coeff
    real(kind=rk) :: turb_coeff1, turb_coeff2
    ! -------------------------------------------------------------------------!
    mu_turb = rans_params%c_mu*state(5)/exp(state(6)/state(1))
    mu_eff = mu + mu_turb
    lam_eff = lambda + 2.0*mu_turb/3.0

    ! Stuff neeeded
    energy_coeff = thermCond/(heatCap*state(1))                               &
    &             + mu_turb*isenCoeff/(rans_params%turb_prandtl_num*state(1)  )

    ! @todo : NA : Calculate limited_eddy_visc Appropriately
    limited_eddy_visc = rans_params%alpha*max(state(5),0.0_rk)              &
      &                                  *exp(state(6)/state(1))
    turb_coeff1 = (mu + rans_params%sig_k*limited_eddy_visc)/ state(1)
    turb_coeff2 = (mu + rans_params%sig_omg*limited_eddy_visc)/ state(1)

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (-mu_eff/state(1)) * velocity(1) * inVec(1) &
            & + (mu_eff/state(1)) * inVec(2)

    ! Third row
    outVec(3) = ((-2.0_rk*mu_eff + lam_eff)/state(1))*velocity(2) * inVec(1) &
            & + ((2.0_rk*mu_eff - lam_eff)/state(1)) * inVec(3)

    ! Fourth row
    outVec(4) = (                                                                &
            &    (-mu_eff/state(1)) * velocity(1)*velocity(1)                    &
            &    + ((-2.0_rk*mu_eff + lam_eff)/state(1))*velocity(2)*velocity(2) &
            &    - energy_coeff*(state(4)/state(1) - sum(velocity(:)**2)         &
            &                    -state(5)/state(1)                    )         &
            &    - turb_coeff1*state(5)/state(1)        ) * inVec(1)             &
            & + (                                                                &
            &     (mu_eff/state(1) - energy_coeff)*velocity(1)                   &
            &   ) * inVec(2)                                                     &
            & + (                                                                &
            &     ((2.0_rk*mu_eff-lam_eff)/state(1) - energy_coeff )*velocity(2) &
            &   ) * inVec(3)                                                     &
            & + (                                                                &
            &     energy_coeff                                                   &
            &   ) * inVec(4)                                                     &
            & + (-energy_coeff + turb_coeff1 )* invec(5)

    ! Fifth row
    outVec(5) = ( -turb_coeff1*(state(5)/state(1)) )* inVec(1)                  &
      &           + ( turb_coeff1 )* inVec(5)

    ! Sixt Row
    outVec(6) = ( -turb_coeff2*state(6)/state(1)  ) * inVec(1)                   &
      &           + ( turb_coeff2 )* inVec(6)

  end function atl_mult_nu22_Rans_2d

end module atl_physFluxFilNvrStk_module

