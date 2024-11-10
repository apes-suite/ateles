! Copyright (c) 2013-2014, 2016-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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

!> Description of the Euler equation system for inviscid compressible flows.
!!
!! The flow field is described in terms of the conservative variables
!! density, momentum and energy.
!! The Euler equations require the definition of the following parameters:
!!
!! * Isentropic expansion coefficient `isen_coef`
!! * Gas constant `r`
!! * Heat capacity at constant volume `cv`
!! * Numerical flux `numflux`
!! * Linearization limit `linear_limit`
!! * Linearization indicator `linearization_indicator`
!! * Wether to ensure positivity modally `ensure_positivity`
!! * Porosity for modeling geometries as porous materials `porosity`
!! * Viscous permeability for the porous material `viscous_permeability`
!! * Thermal permeability for the porous material `thermal_permeability`
!! * Definition of the material (to model geometries) `material`
!!
!! There are three different implementations for the numerical flux available:
!!
!! * 'hll' Harten/Lax/van-Leer flux
!! * 'godunov' Solve the nonlinear Riemann problem
!! * 'lax_friedrich' Lax/Friedrichs flux, not recommended and may lead to
!!   instabilities, especially at boundaries
!!
!! It is possible to adaptively linearize the flux in elements according to its
!! state. This is activated by setting the `linear_limit` to some value larger
!! than 0. The default is to use no adaptive linearization (`linear_limit=0`).
!! There are three different indicators available that may be used to decide the
!! linearization of elements. `linearization_indicator` may be:
!!
!! * 'error' A modal estimation for the terms neglected by linearization.
!!   This involves all variables and is the default.
!! * 'energy' Use the deviation of the energy to decide the linearization.
!! * 'density' Use the deviation of the density to decide the linearization.
!!
!! The deviations are considered in relation to the cell mean value. Thus,
!! setting `linear_limit = 0.01` and `linearization_indicator = 'density'`
!! would only compute the linearized flux in elements where the maximal
!! deviation of the state is guaranteed to be less or equal to one percent
!! of the mean density in that cell.
!!
!! For the material the following variables need to be defined:
!!
!! * `characteristic` masking function, has to 1 where a geometry is to be found
!!   and 0 elsewhere.
!! * `relax_velocity` the velocity of the object
!! * `relax_temperature` the temperature of the object
!!
!! If the `characteristic` is 0 everywhere, the other settings are irrelevant.
!! Each material component is of [[tem_variable_type]]. However, these can be
!! also simply defined as space-time functions, which can just be constants.
!! See the [[tem_variable_module]] and [[tem_spacetime_fun_module]] for more
!! details.
!!
!! If `ensure_positivity = true` is set, the state will be scaled such that the
!! deviation of the state is guaranteed to be smaller than the mean for density
!! and energy such that these are positive.
!! While this is a computationally cheap measure, you should be aware that it is
!! not sufficient to avoid unphysical states, as the pressure may still be
!! negative.
!!
!! The following equation names are implementing the Euler equations and require
!! the settings described above:
!!
!! * `euler` (3D)
!! * `euler_2d` (2D)
!! * `euler_1d` (1D)
!! * `loclineuler` Locally linearized (3D)
!! * `loclineuler_1d` Locally linearized (1D)
!!
!! The locally linearized variants use fluxes linearized to the cell means
!! within the elements, while between elements the nonlinear flux is computed.
!!
!! An example configuration for Euler equations is shown below:
!!
!!```lua
!!  gamma = 1.4 -- isentropic expansion coefficient
!!  R     = 1.0 -- gas constant
!!  equation = {
!!    name              = 'euler',
!!    isen_coef         = gamma,
!!    r                 = R,
!!    cv                = R/(gamma-1),
!!    numflux           = 'hll',
!!
!!    linear_limit      = 0.001,
!!    linearization_indicator = 'error',
!!
!!    ensure_positivity = false, -- Ensure that density and energy remain
!!                               -- positive by only considering higher modes up
!!                               -- to the point where positivity is guaranteed.
!!
!!    porosity          = 1.0,   -- Porosity to use in material modelling for
!!                               -- wall representation in elements.
!!    viscous_permeability = 1.0e-6, -- Viscous permeability for the porous
!!                                   -- medium to represent wall geometries in
!!                                   -- elements.
!!    thermal_permeability = 1.0e-3, -- Thermal permeability for the material to
!!                                   -- represent walls.
!!    material = {
!!      -- Description of the material distribution to define obstacles inside
!!      -- the domain.
!!      characteristic = 0.0, -- Masking function (may be a variable) that
!!                            -- describes where Material is to be found
!!                            -- (Chi(x,y,z)), should be 1 inside material and 0
!!                            -- everywhere else.
!!      relax_velocity = {0, 0, 0}, -- Velocity of the obstacle.
!!      relax_temperature = 0.0     -- Temperature of the obstacle.
!!    }
!!  }
!!
module atl_eqn_euler_module
  ! Treelm modules
  use env_module,             only: rk

  use tem_aux_module,         only: tem_Abort
  use tem_logging_module,     only: logUnit

  ! Aotus modules
  use aotus_module,           only: flu_State, aot_get_val
  use aot_out_module,         only: aot_out_type, &
    &                               aot_out_val

  implicit none

  private

  public :: atl_euler_type
  public :: atl_load_euler
  public :: atl_save_euler

  !> The Euler equation properties are stored here
  type atl_euler_type
    ! Ideal gas phase parameters
    !> Isentropic coefficient
    real(kind=rk) :: isen_coef
    !> Ideal gas constant
    real(kind=rk) :: r
    !> Heat capacity
    real(kind=rk) :: cv

    !> Limit in energy variation up to which a linearization is to be
    !! allowed in elements.
    !!
    !! Set this to 0 to never allow linearization (default).
    !! Otherwise a linear flux computation will be used within elements,
    !! if the energy is guaranteed to not deviate beyond this limit from
    !! the mean energy level in the element.
    !!
    !! That is the sum of absolute values of all higher modes is to be
    !! smaller then `linear_limit` multiplied with the first mode.
    !!
    !! The energy is used as all characteristic variables contribute to it.
    real(kind=rk) :: linear_limit
    !> Ensure that polynomials for density and energy are everywhere positive
    !! after transfer to oversampled space.
    !!
    !! Modes will only be considered up to the point where the sum of their
    !! absolute values is smaller than the first mode.
    logical :: ensure_positivity

    !> Flag for adaptive timestep calculation
    logical :: adaptive_timestep

    ! Parameter for the Brinkamnn penalization
    !> Porosity (for the continuity equation)
    real(kind=rk) :: porosity
    !> Porosity for the momentum transport
    real(kind=rk) :: viscous_permeability
    !> Porosity for the heat transport
    real(kind=rk) :: thermal_permeability

    !> Procedure to compute the numerical flux for the equation at hand.
    !!
    !! What kind of fluxes are available depends on the equation that is
    !! being solved.
    procedure(euler_numflux), pointer, pass(euler) :: numflux => NULL()

    !> Function to decide, whether linearized flux computation would
    !! be tolerable.
    !!
    !! Set in the configuration by linearization_indicator.
    !! Available are:
    !!
    !! - 'density' use density deviation as indicator
    !! - 'energy' use energy deviation as indicator
    !! - 'error' use an error estimate as indicator
    !! Default is to use the error estimate.
    procedure(linearization_indicator), pointer, pass(euler) :: linear => NULL()

  end type atl_euler_type
  ! ----------------------------------------------------------------------------


  ! ----------------------------------------------------------------------------
  abstract interface
    !> Interface definition for numerical fluxes.
    !! @todo HK: should be vectorized to reduce overheads.
    subroutine euler_numflux(euler, state_left, state_right, &
      &                        material_left, material_right, nPoints, flux)
      import :: atl_euler_type, rk
      !> The global equation specific settings.
      class(atl_euler_type), intent(in) :: euler

      !> State on the left side of the interface.
      !!
      !! First index has to correspond with nPoints, second with number of
      !! variables.
      real(kind=rk), intent(in)  :: state_left(:,:)

      !> State on the right side of the interface.
      !!
      !! First index has to correspond with nPoints, second with number of
      !! variables.
      real(kind=rk), intent(in)  :: state_right(:,:)

      !> Material definition on the left of the interface.
      !!
      !! First index corresponds to nPoints or is 1 for constant material,
      !! second index has to correspond with number of material settings.
      real(kind=rk), intent(in)  :: material_left(:,:)

      !> Material definition on the right of the interface.
      !!
      !! First index corresponds to nPoints or is 1 for constant material,
      !! second index has to correspond with number of material settings.
      real(kind=rk), intent(in)  :: material_right(:,:)

      !> Number of points to evaluate the flux at.
      integer, intent(in) :: nPoints

      !> Resulting flux.
      !!
      !! First index has to correspond with nPoints, second with number of
      !! variables.
      real(kind=rk), intent(out) :: flux(:,:)
    end subroutine euler_numflux


    !> Computation of indicator, whether a linear flux is sufficient.
    pure function linearization_indicator(euler, mean, deviation) &
      &        result(islinear)
      import :: atl_euler_type, rk
      ! ---------------------------------------------------------------- !
      !> Description of the equation
      class(atl_euler_type), intent(in) :: euler

      !> The mean value of each state
      real(kind=rk), intent(in) :: mean(:)

      !> Maximal deviation of each state
      real(kind=rk), intent(in) :: deviation(:)

      logical :: islinear
      ! ---------------------------------------------------------------- !
    end function linearization_indicator

  end interface
  ! ----------------------------------------------------------------------------


contains


  ! ****************************************************************************
  !> Subroutine to initialize an equation of type euler equation
  !! as defined in the configuration file
  subroutine atl_load_euler( euler, conf, eq_table )
    ! --------------------------------------------------------------------------
    !> Resulting description of the Euler equation parameters.
    type(atl_Euler_type), intent(out) :: euler

    !> Handle to the configuration script, to load the parameters from.
    type(flu_State) :: conf

    !> Handle to the table containing the description for the equation
    !! system.
    integer, intent(in) :: eq_table

    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------

    !read the data from the equation table of the lua file
    call aot_get_val( L       = conf,            &
      &               thandle = eq_table,        &
      &               key     = 'isen_coef',     &
      &               val     = euler%isen_coef, &
      &               ErrCode = iError           )
    if(iError.ne.0) then
      call tem_abort( 'ERROR: not able to find isentropic coefficient in ' &
        & // 'equation table, stopping ...'                                )
    end if

    call aot_get_val( L       = conf,     &
      &               thandle = eq_table, &
      &               key     = 'r',      &
      &               val     = euler%r,  &
      &               ErrCode = iError    )
    if(iError.ne.0) then
      call tem_abort( 'ERROR: not able to find ideal gas constant value r in ' &
        & // 'equation table, stopping ...'                                    )
    end if

    call aot_get_val( L       = conf,     &
      &               thandle = eq_table, &
      &               key     = 'cv',     &
      &               val     = euler%cv, &
      &               ErrCode = iError    )
    if(iError.ne.0) then
      call tem_abort( 'ERROR: not able to find volumetric heat capacity in ' &
        & // 'equation table, stopping ...'                                  )
    end if

    call aot_get_val( L       = conf,               &
      &               thandle = eq_table,           &
      &               key     = 'linear_limit',     &
      &               val     = euler%linear_limit, &
      &               default = 0.0_rk,             &
      &               ErrCode = iError              )

    call aot_get_val( L       = conf,                    &
      &               thandle = eq_table,                &
      &               key     = 'ensure_positivity',     &
      &               val     = euler%ensure_positivity, &
      &               default = .false.,                 &
      &               ErrCode = iError                   )

    call aot_get_val( L       = conf,           &
      &               thandle = eq_table,       &
      &               key     = 'porosity',     &
      &               val     = euler%porosity, &
      &               default = 1.0_rk,         &
      &               ErrCode = iError          )
    call aot_get_val( L       = conf,                       &
      &               thandle = eq_table,                   &
      &               key     = 'viscous_permeability',     &
      &               val     = euler%viscous_permeability, &
      &               default = 1.0_rk,                     &
      &               ErrCode = iError                      )
    call aot_get_val( L       = conf,                       &
      &               thandle = eq_table,                   &
      &               key     = 'thermal_permeability',     &
      &               val     = euler%thermal_permeability, &
      &               default = 1.0_rk,                     &
      &               ErrCode = iError                      )

    write(logUnit(1),*) 'Loaded Euler equation'
    write(logUnit(1),*) ' * isen_coef: ', Euler%isen_coef
    write(logUnit(1),*) ' * R: ', Euler%r
    write(logUnit(1),*) ' * cv: ', Euler%cv
    write(logUnit(1),*) ' * porosity: ', Euler%porosity
    write(logUnit(1),*) ' * viscous_permeability: ', Euler%viscous_permeability
    write(logUnit(1),*) ' * thermal_permeability: ', Euler%thermal_permeability
    if (euler%linear_limit > 0.0_rk) then
      write(logunit(1),*) ' * Equation is chosen adaptively'
      write(logunit(1),*) '   In elements, where the indicator does not vary by'
      write(logunit(1),*) '   more than ', euler%linear_limit*100, '%,'
      write(logunit(1),*) '   locally linearized physical flux will be used!'
      write(logUnit(1),*) ''
    end if
    if (euler%ensure_positivity) then
      write(logUnit(1),*) ''
      write(logunit(1),*) ' * Positivity for density and energy is ensured'
      write(logunit(1),*) '   by limiting the order of the polynomials when'
      write(logunit(1),*) '   transferring the modes to the oversampled space.'
      write(logunit(1),*) '   Higher modes will only be considered up to the'
      write(logunit(1),*) '   point, where the sum of their absolute values'
      write(logunit(1),*) '   is smaller than the first mode.'
    end if

  end subroutine atl_load_euler
  ! ****************************************************************************


  ! ****************************************************************************
  ! Save the equation variables into the lua file
  subroutine atl_save_euler(me, eqn_name, conf, withName)
    ! --------------------------------------------------------------------------
    type(atl_euler_type), intent(in) :: me
    character(len=*), intent(in) :: eqn_name
    type(aot_out_type) :: conf
    !> Defines whether or not to print the equation name. Default is true.
    logical, intent(in), optional :: withName
    ! --------------------------------------------------------------------------
    logical :: loc_wName
    ! --------------------------------------------------------------------------

    ! The equation name is printed by default, i.e. when withName is not present
    ! If it is present, it's value is used to decide whether or not to print the
    ! equatio name.
    loc_wName = .true.
    if (present(withName)) loc_wName = withName
    if ( loc_wName ) then
      call aot_out_val(put_conf = conf, vname = 'name', val = trim(eqn_name))
    end if

    ! Dump equation Properties
    call aot_out_val(put_conf = conf, vname = 'isen_coef', val = me%isen_coef)
    call aot_out_val(put_conf = conf, vname = 'r',         val = me%r)
    call aot_out_val(put_conf = conf, vname = 'cv',        val = me%cv)

    ! For the penalizations
    call aot_out_val( put_conf = conf,       &
      &               vname    = 'porosity', &
      &               val      = me%porosity )
    call aot_out_val( put_conf = conf,                   &
      &               vname    = 'viscous_permeability', &
      &               val      = me%viscous_permeability )
    call aot_out_val( put_conf = conf,                   &
      &               vname    = 'thermal_permeability', &
      &               val = me%thermal_permeability      )

  end subroutine atl_save_euler
  ! ****************************************************************************


end module atl_eqn_euler_module
