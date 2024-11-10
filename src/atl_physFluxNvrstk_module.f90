! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
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
!! Collects all functions related to the physical fluxes of the compressible
!! Navier-Stokes equations.
module atl_physFluxNvrStk_module
  ! Treelm modules
  use env_module,            only: rk

  implicit none
  private


  public :: atl_viscPhysFluxNavierStokes, &
    &       atl_mult_nu11_NavierStokes, &
    &       atl_mult_nu12_NavierStokes, &
    &       atl_mult_nu13_NavierStokes, &
    &       atl_mult_nu21_NavierStokes, &
    &       atl_mult_nu22_NavierStokes, &
    &       atl_mult_nu23_NavierStokes, &
    &       atl_mult_nu31_NavierStokes, &
    &       atl_mult_nu32_NavierStokes, &
    &       atl_mult_nu33_NavierStokes


contains


  !> Physical flux calculation along x direction for Euler equation.
  function atl_viscPhysFluxNavierStokes( state, state_gradient, mu, lambda, &
    &                                    thermCond, heatCap                 ) &
    &      result(physFlux)
    ! -------------------------------------------------------------------- !
    !> The state in nodal space. Dimension is the number of state variables
    !! i.e. 5 for Navier-Stokes.
    real(kind=rk), intent(in) :: state(:)
    !> The state in nodal space. First dimension is the number of variables
    !! i.e. 5 for Navier-Stokes.
    !! Second dimension is the dimension, e.g. 3 in two dimensions.
    real(kind=rk), intent(in) :: state_gradient(:,:)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The physical flux along the x axis for all variables
    real(kind=rk) :: physFlux(5)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: velocity(1:3)
    ! -------------------------------------------------------------------- !

    velocity = state(2:4)/state(1)

    ! Viscous flux for density
    physFlux(1) = 0.0_rk

    ! Viscous flux for momentum in x
    physFlux(2) =                                                              &
                ! (nu_{1,1})_{2,1} (\nabla u)_{1,1} : k = 1, i = 1
      &         ((-2.0_rk * mu + lambda) / state(1))  &
      &         * velocity(1)*state_gradient(1,1)     &
                ! (nu_{1,2})_{2,1} (\nabla u)_{1,2} : k = 2, i = 1
      &         + (lambda / state(1))                 &
      &           * velocity(2)*state_gradient(1,2)   &
                ! (nu_{1,3})_{2,1} (\nabla u)_{1,3} : k = 3, i = 1
      &         + (lambda / state(1))                 &
      &           * velocity(3)*state_gradient(1,3)   &

                ! (nu_{1,1})_{2,2} (\nabla u)_{2,1} : k = 1, i = 2
      &         + ((2.0_rk * mu - lambda) / state(1)) &
      &           * state_gradient(2,1)               &
                ! (nu_{1,2})_{2,2} (\nabla u)_{2,2} : k = 2, i = 2
      &         + 0.0_rk                              &
                ! (nu_{1,3})_{2,2} (\nabla u)_{2,3} : k = 3, i = 2
      &         + 0.0_rk                              &

                ! (nu_{1,1})_{2,3} (\nabla u)_{3,1} : k = 1, i = 3
      &         + 0.0_rk                              &
                ! (nu_{1,2})_{2,3} (\nabla u)_{3,2} : k = 2, i = 3
      &         + (-lambda/state(1))                  &
      &           * state_gradient(3,2)               &
                ! (nu_{1,3})_{2,3} (\nabla u)_{3,3} : k = 2, i = 3
      &         + 0.0_rk                              &

                ! (nu_{1,1})_{2,4} (\nabla u)_{4,1} : k = 1, i = 4
      &         + 0.0_rk                              &
                ! (nu_{1,2})_{2,4} (\nabla u)_{4,2} : k = 2, i = 4
      &         + 0.0_rk                              &
                ! (nu_{1,3})_{2,4} (\nabla u)_{4,3} : k = 3, i = 4
      &         + (-lambda/state(1))                  &
      &           * state_gradient(4,3)               &
                ! (nu_{1,1})_{2,5} (\nabla u)_{5,1} : k = 1, i = 5
      &         + 0.0_rk                              &
                ! (nu_{1,2})_{2,5} (\nabla u)_{5,2} : k = 2, i = 5
      &         + 0.0_rk                              &
                ! (nu_{1,3})_{2,5} (\nabla u)_{5,3} : k = 3, i = 5
      &         + 0.0_rk

    ! Viscous flux for momentum in y
    physFlux(3) =                                                &
              ! (nu_{1,1})_{3,1} (\nabla u)_{1,1} : k = 1, i = 1
      &       (-mu / state(1))*velocity(2)*state_gradient(1,1)   &
              ! (nu_{1,2})_{3,1} (\nabla u)_{1,2} : k = 2, i = 1
      &       + (-mu / state(1))*velocity(1)*state_gradient(1,2) &
              ! (nu_{1,3})_{3,1} (\nabla u)_{1,3} : k = 3, i = 1
      &       + 0.0_rk                                           &

              ! (nu_{1,1})_{3,2} (\nabla u)_{2,1} : k = 1, i = 2
      &       + 0.0_rk                                           &
              ! (nu_{1,2})_{3,2} (\nabla u)_{2,2} : k = 2, i = 2
      &       + (mu/state(1))*state_gradient(2,2)                &
              ! (nu_{1,3})_{3,2} (\nabla u)_{2,3} : k = 3, i = 2
      &       + 0.0_rk                                           &

              ! (nu_{1,1})_{3,3} (\nabla u)_{3,1} : k = 1, i = 3
      &       + (mu/state(1))*state_gradient(3,1)                &
              ! (nu_{1,2})_{3,3} (\nabla u)_{3,2} : k = 2, i = 3
      &       + 0.0_rk                                           &
              ! (nu_{1,3})_{3,3} (\nabla u)_{3,3} : k = 3, i = 3
      &       + 0.0_rk                                           &

              ! (nu_{1,1})_{3,4} (\nabla u)_{4,1} : k = 1, i = 4
      &       + 0.0_rk                                           &
              ! (nu_{1,2})_{3,4} (\nabla u)_{4,2} : k = 2, i = 4
      &       + 0.0_rk                                           &
              ! (nu_{1,3})_{3,4} (\nabla u)_{4,3} : k = 3, i = 4
      &       + 0.0_rk                                           &

              ! (nu_{1,1})_{3,5} (\nabla u)_{5,1} : k = 1, i = 5
      &       + 0.0_rk                                           &
              ! (nu_{1,2})_{3,5} (\nabla u)_{5,2} : k = 2, i = 5
      &       + 0.0_rk                                           &
              ! (nu_{1,3})_{3,5} (\nabla u)_{5,3} : k = 3, i = 5
      &       + 0.0_rk

    ! Viscous flux for momentum in z
    physFlux(4) =                                                  &
              ! (nu_{1,1})_{4,1} (\nabla u)_{1,1} : k = 1, i = 1
      &       (-mu / state(1))*velocity(3)*state_gradient(1,1)     &
              ! (nu_{1,2})_{4,1} (\nabla u)_{1,2} : k = 2, i = 1
      &       + 0.0_rk                                             &
              ! (nu_{1,3})_{4,1} (\nabla u)_{1,3} : k = 3, i = 1
      &       + (-mu / state(1))*velocity(1)*state_gradient(1,3)   &

              ! (nu_{1,1})_{4,2} (\nabla u)_{2,1} : k = 1, i = 2
      &       + 0.0_rk                                             &
              ! (nu_{1,2})_{4,2} (\nabla u)_{2,2} : k = 2, i = 2
      &       + 0.0_rk                                             &
              ! (nu_{1,3})_{4,2} (\nabla u)_{2,3} : k = 3, i = 2
      &       + (mu/state(1))*state_gradient(2,3)                  &

              ! (nu_{1,1})_{4,3} (\nabla u)_{3,1} : k = 1, i = 3
      &       + 0.0_rk                                             &
              ! (nu_{1,2})_{4,3} (\nabla u)_{3,2} : k = 2, i = 3
      &       + 0.0_rk                                             &
              ! (nu_{1,3})_{4,3} (\nabla u)_{3,3} : k = 3, i = 3
      &       + 0.0_rk                                             &

              ! (nu_{1,1})_{4,4} (\nabla u)_{4,1} : k = 1, i = 4
      &       + (mu/state(1))*state_gradient(4,1)                  &
              ! (nu_{1,2})_{4,4} (\nabla u)_{4,2} : k = 2, i = 4
      &       + 0.0_rk                                             &
              ! (nu_{1,3})_{4,4} (\nabla u)_{4,3} : k = 3, i = 4
      &       + 0.0_rk                                             &

              ! (nu_{1,1})_{4,5} (\nabla u)_{5,1} : k = 1, i = 5
      &       + 0.0_rk                                             &
              ! (nu_{1,2})_{4,5} (\nabla u)_{5,2} : k = 2, i = 5
      &       + 0.0_rk                                             &
              ! (nu_{1,3})_{4,5} (\nabla u)_{5,3} : k = 3, i = 5
      &       + 0.0_rk

    ! Viscous flux for Energy
    physFlux(5) = &
              ! (nu_{1,1})_{5,1} (\nabla u)_{1,1} : k = 1, i = 1
      &       ( ((-2.0_rk * mu + lambda) / state(1))*(velocity(1)**2.0_rk) &
      &         + (-mu/state(1))*(velocity(2)**2.0_rk)                     &
      &         + (-mu/state(1))*(velocity(3)**2.0_rk)                     &
      &         - (thermCond/(heatCap*state(1)))                           &
      &           * ( state(5)/state(1)                                    &
      &               - sum(velocity(:)**2.0_rk) )                         &
      &       ) * state_gradient(1,1)                                      &
              ! (nu_{1,2})_{5,1} (\nabla u)_{1,2} : k = 2, i = 1
      &       + ((-mu+lambda)/state(1))                                    &
      &         * velocity(1)*velocity(2)*state_gradient(1,2)              &
              ! (nu_{1,3})_{5,1} (\nabla u)_{1,3} : k = 3, i = 1
      &       + ((-mu+lambda)/state(1))                                    &
      &         * velocity(1)*velocity(3)*state_gradient(1,3)              &

              ! (nu_{1,1})_{5,2} (\nabla u)_{2,1} : k = 1, i = 2
      &       + ( (2.0_rk*mu-lambda)/state(1)                              &
      &           - thermCond/(heatCap*state(1))                           &
      &         ) * velocity(1)                                            &
      &         * state_gradient(2,1)                                      &
              ! (nu_{1,2})_{5,2} (\nabla u)_{2,2} : k = 2, i = 2
      &       + (mu/state(1))*velocity(2)*state_gradient(2,2)              &
              ! (nu_{1,3})_{5,2} (\nabla u)_{2,3} : k = 3, i = 2
      &       + (mu/state(1))*velocity(3)*state_gradient(2,3)              &

              ! (nu_{1,1})_{5,3} (\nabla u)_{3,1} : k = 1, i = 3
      &       + ((mu/state(1))-thermCond/(heatCap*state(1)))               &
      &         * velocity(2)*state_gradient(3,1)                          &
              ! (nu_{1,2})_{5,3} (\nabla u)_{3,2} : k = 2, i = 3
      &       + (-lambda/state(1))*velocity(1)*state_gradient(3,2)         &
              ! (nu_{1,3})_{5,3} (\nabla u)_{3,3} : k = 3, i = 3
      &       + 0.0_rk                                                     &

              ! (nu_{1,1})_{5,4} (\nabla u)_{4,1} : k = 1, i = 4
      &       + ((mu/state(1))-thermCond/(heatCap*state(1)))               &
      &         * velocity(3)*state_gradient(4,1)                          &
              ! (nu_{1,2})_{5,4} (\nabla u)_{4,2} : k = 2, i = 4
      &       + 0.0_rk                                                     &
              ! (nu_{1,3})_{5,4} (\nabla u)_{4,3} : k = 3, i = 4
      &       + (-lambda/state(1))*velocity(1)*state_gradient(4,3)         &

              ! (nu_{1,1})_{5,5} (\nabla u)_{5,1} : k = 1, i = 5
      &       + (thermCond/(heatCap*state(1)))*state_gradient(5,1)         &
              ! (nu_{1,2})_{5,5} (\nabla u)_{5,2} : k = 2, i = 5
      &       + 0.0_rk                                                     &
              ! (nu_{1,3})_{5,5} (\nabla u)_{5,3} : k = 3, i = 5
      &       + 0.0_rk

  end function atl_viscPhysFluxNavierStokes

  ! Multiplies the viscous flux matrux nu_11 with a given vector
  function atl_mult_nu11_NavierStokes(density, velocity, totEnergy, inVec, &
                                        & mu, lambda, thermCond, heatCap) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> The total energy
    real(kind=rk), intent(in) :: totEnergy
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = ((-2.0_rk*mu + lambda)/density)*velocity(1) * inVec(1) &
            & + ((2.0_rk*mu - lambda)/density) * inVec(2)

    ! Third row
    outVec(3) = (-mu/density) * velocity(2) * inVec(1) &
            & + (mu/density) * inVec(3)

    ! Fourth row
    outVec(4) = (-mu/density) * velocity(3) * inVec(1) &
            & + (mu/density) * inVec(4)

    ! Fifth row
    outVec(5) = (                                                                           &
            &       ((-2.0_rk*mu + lambda)/density)*velocity(1)*velocity(1)                 &
            &     + (-mu/density) * velocity(2)*velocity(2)                                 &
            &     + (-mu/density) * velocity(3)*velocity(3)                                 &
            &     - (thermCond/(heatCap*density))*(totEnergy/density - sum(velocity(:)**2)) &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &     ((2.0_rk*mu-lambda)/density - thermCond/(heatCap*density))*velocity(1)    &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &     (mu/density - thermCond/(heatCap*density))*velocity(2)                    &
            &   ) * inVec(3)                                                                &
            & + (                                                                           &
            &     (mu/density - thermCond/(heatCap*density))*velocity(3)                    &
            &   ) * inVec(4)                                                                &
            & + (                                                                           &
            &     (thermCond/(heatCap*density))                                             &
            &   ) * inVec(5)

  end function atl_mult_nu11_NavierStokes

  ! Multiplies the viscous flux matrux nu_21 with a given vector
  function atl_mult_nu21_NavierStokes(density, velocity, inVec, mu, lambda) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (-mu/density) * velocity(2) * inVec(1) &
            & + (mu/density) * inVec(3)

    ! Third row
    outVec(3) = (lambda/density) * velocity(1) * inVec(1) &
            & + (-lambda/density) * inVec(2)

    ! Fourth row has zeros only
    outVec(4) = 0.0_rk

    ! Fifth row
    outVec(5) = (                                                                           &
            &      ((-mu+lambda)/density)*velocity(1)*velocity(2)                           &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &      (-lambda/density)*velocity(2)                                            &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &      (mu/density)*velocity(1)                                                 &
            &   ) * inVec(3)

  end function atl_mult_nu21_NavierStokes


  ! Multiplies the viscous flux matrux nu_31 with a given vector
  function atl_mult_nu31_NavierStokes(density, velocity, inVec, mu, lambda) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (-mu/density) * velocity(3) * inVec(1) &
            & + (mu/density) * inVec(4)

    ! Third row has zeros only
    outVec(3) = 0.0_rk

    ! Fourth row
    outVec(4) = (lambda/density) * velocity(1) * inVec(1) &
            & + (-lambda/density) * inVec(2)

    ! Fifth row
    outVec(5) = (                                                                           &
            &      ((-mu+lambda)/density)*velocity(1)*velocity(3)                           &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &      (-lambda/density)*velocity(3)                                            &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &      (mu/density)*velocity(1)                                                 &
            &   ) * inVec(4)

  end function atl_mult_nu31_NavierStokes


  ! Multiplies the viscous flux matrux nu_12 with a given vector
  function atl_mult_nu12_NavierStokes(density, velocity, inVec, mu, lambda) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (lambda/density) * velocity(2) * inVec(1) &
            & + (-lambda/density) * inVec(3)

    ! Third row
    outVec(3) = (-mu/density) * velocity(1) * inVec(1) &
            & + (mu/density) * inVec(2)

    ! Fourth row has zeros only
    outVec(4) = 0.0_rk

    ! Fifth row
    outVec(5) = (                                                                           &
            &      ((-mu+lambda)/density)*velocity(1)*velocity(2)                           &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &      (mu/density)*velocity(2)                                                 &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &      (-lambda/density)*velocity(1)                                            &
            &   ) * inVec(3)

  end function atl_mult_nu12_NavierStokes

  ! Multiplies the viscous flux matrux nu_13 with a given vector
  function atl_mult_nu13_NavierStokes(density, velocity, inVec, mu, lambda) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (lambda/density) * velocity(3) * inVec(1) &
            & + (-lambda/density) * inVec(4)

    ! Third row has zeros only
    outVec(3) = 0.0_rk

    ! Fourth row
    outVec(4) = (-mu/density) * velocity(1) * inVec(1) &
            & + (mu/density) * inVec(2)

    ! Fifth row
    outVec(5) = (                                                                           &
            &      ((-mu+lambda)/density)*velocity(1)*velocity(3)                           &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &      (mu/density)*velocity(3)                                                 &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &      (-lambda/density)*velocity(1)                                            &
            &   ) * inVec(4)

  end function atl_mult_nu13_NavierStokes

  ! Multiplies the viscous flux matrux nu_23 with a given vector
  function atl_mult_nu23_NavierStokes(density, velocity, inVec, mu, lambda) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row has zeros only
    outVec(2) = 0.0_rk

    ! Third row
    outVec(3) = (lambda/density) * velocity(3) * inVec(1) &
            & + (-lambda/density) * inVec(4)

    ! Fourth row
    outVec(4) = (-mu/density) * velocity(2) * inVec(1) &
            & + (mu/density) * inVec(3)

    ! Fifth row
    outVec(5) = (                                                                           &
            &      ((-mu+lambda)/density)*velocity(2)*velocity(3)                           &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &      (mu/density)*velocity(3)                                                 &
            &   ) * inVec(3)                                                                &
            & + (                                                                           &
            &      (-lambda/density)*velocity(2)                                            &
            &   ) * inVec(4)

  end function atl_mult_nu23_NavierStokes

  ! Multiplies the viscous flux matrux nu_32 with a given vector
  function atl_mult_nu32_NavierStokes(density, velocity, inVec, mu, lambda) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row has zeros only
    outVec(2) = 0.0_rk

    ! Third row
    outVec(3) = (-mu/density) * velocity(3) * inVec(1) &
            & + (mu/density) * inVec(4)

    ! Fourth row
    outVec(4) = (lambda/density) * velocity(2) * inVec(1) &
            & + (-lambda/density) * inVec(3)

    ! Fifth row
    outVec(5) = (                                                                           &
            &      ((-mu+lambda)/density)*velocity(2)*velocity(3)                           &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &      (-lambda/density)*velocity(3)                                            &
            &   ) * inVec(3)                                                                &
            & + (                                                                           &
            &      (mu/density)*velocity(2)                                                 &
            &   ) * inVec(4)

  end function atl_mult_nu32_NavierStokes

  ! Multiplies the viscous flux matrux nu_22 with a given vector
  function atl_mult_nu22_NavierStokes(density, velocity, totEnergy, inVec, &
                                        & mu, lambda, thermCond, heatCap) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> The total energy
    real(kind=rk), intent(in) :: totEnergy
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (-mu/density) * velocity(1) * inVec(1) &
            & + (mu/density) * inVec(2)

    ! Third row
    outVec(3) = ((-2.0_rk*mu + lambda)/density)*velocity(2) * inVec(1) &
            & + ((2.0_rk*mu - lambda)/density) * inVec(3)

    ! Fourth row
    outVec(4) = (-mu/density) * velocity(3) * inVec(1) &
            & + (mu/density) * inVec(4)

    ! Fifth row
    outVec(5) = (                                                                           &
            &       (-mu/density) * velocity(1)*velocity(1)                                 &
            &     + ((-2.0_rk*mu + lambda)/density)*velocity(2)*velocity(2)                 &
            &     + (-mu/density) * velocity(3)*velocity(3)                                 &
            &     - (thermCond/(heatCap*density))*(totEnergy/density - sum(velocity(:)**2)) &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &     (mu/density - thermCond/(heatCap*density))*velocity(1)                    &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &     ((2.0_rk*mu-lambda)/density - thermCond/(heatCap*density))*velocity(2)    &
            &   ) * inVec(3)                                                                &
            & + (                                                                           &
            &     (mu/density - thermCond/(heatCap*density))*velocity(3)                    &
            &   ) * inVec(4)                                                                &
            & + (                                                                           &
            &     (thermCond/(heatCap*density))                                             &
            &   ) * inVec(5)

  end function atl_mult_nu22_NavierStokes

  ! Multiplies the viscous flux matrux nu_33 with a given vector
  function atl_mult_nu33_NavierStokes(density, velocity, totEnergy, inVec, &
                                        & mu, lambda, thermCond, heatCap) &
                                        & result( outVec )
    ! -------------------------------------------------------------------------!
    !> The density
    real(kind=rk), intent(in) :: density
    !> The velocity
    real(kind=rk), intent(in) :: velocity(3)
    !> The total energy
    real(kind=rk), intent(in) :: totEnergy
    !> Vector to be multiplied with nu11
    real(kind=rk), intent(in) :: inVec(5)
    !> Dynamic Viscosity
    real(kind=rk), intent(in) :: mu
    !> Viscosity
    real(kind=rk), intent(in) :: lambda
    !> The thermal cond
    real(kind=rk), intent(in) :: thermCond
    !> The specific heat capacity (per mass unit mass, at constant volume)
    real(kind=rk), intent(in) :: heatCap
    !> The result of the matrix vector product
    real(kind=rk) :: outVec(5)
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    ! First row has zeros only
    outVec(1) = 0.0_rk

    ! Second row
    outVec(2) = (-mu/density) * velocity(1) * inVec(1) &
            & + (mu/density) * inVec(2)

    ! Third row
    outVec(3) = (-mu/density) * velocity(2) * inVec(1) &
            & + (mu/density) * inVec(3)

    ! Fourth row
    outVec(4) = ((-2.0_rk*mu + lambda)/density)*velocity(3) * inVec(1) &
            & + ((2.0_rk*mu - lambda)/density) * inVec(4)

    ! Fifth row
    outVec(5) = (                                                                           &
            &       (-mu/density) * velocity(1)*velocity(1)                                 &
            &     + (-mu/density) * velocity(2)*velocity(2)                                 &
            &     + ((-2.0_rk*mu + lambda)/density)*velocity(3)*velocity(3)                 &
            &     - (thermCond/(heatCap*density))*(totEnergy/density - sum(velocity(:)**2)) &
            &   ) * inVec(1)                                                                &
            & + (                                                                           &
            &     (mu/density - thermCond/(heatCap*density))*velocity(1)                    &
            &   ) * inVec(2)                                                                &
            & + (                                                                           &
            &     (mu/density - thermCond/(heatCap*density))*velocity(2)                    &
            &   ) * inVec(3)                                                                &
            & + (                                                                           &
            &     ((2.0_rk*mu-lambda)/density - thermCond/(heatCap*density))*velocity(3)    &
            &   ) * inVec(4)                                                                &
            & + (                                                                           &
            &     (thermCond/(heatCap*density))                                             &
            &   ) * inVec(5)

  end function atl_mult_nu33_NavierStokes

end module atl_physFluxNvrStk_module

