! Copyright (c) 2013-2016, 2019-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

!> Compressible Navier-Stokes equations for viscous flows.
!!
!! The flow field is represented in terms of the conservative variables
!! density, momentum and energy.
!! The Navier-Stokes equations require the definitions of the inviscous Euler
!! equations, as described in [[atl_eqn_euler_module]] and in addition the
!! following parameters:
!!
!! * Dynamic viscosity `mu`
!! * Thermal conductivity `therm_cond`
!! * Internal penalization (for the diffusive term implementation in DG)
!!
!! Optionally you can also provide a `visc_limit` to adaptively switch off
!! the computation of the viscous terms if their estimated magnitude falls
!! below this limit for an element. If `visc_limit` is 0 (the default) no
!! adaptivity will be considered.
!!
!! The following equation names are implementing the Navier-Stokes equations
!! and require the settings described above.
!!
!! * `navier_stokes` (3D)
!! * `navier_stokes_2d` (2D)
!! * `filtered_navier_stokes` (3D) with turbulence modelling
!! * `filtered_navier_stokes_2d` (2D) with turbulence modelling
!!
!! A complete example for the configuration of the Navier-Stokes equations
!! is given below:
!!
!!```lua
!!  equation = {
!!    name      = 'navier_stokes',
!!    isen_coef = 1.4,
!!    r         = 287,
!!    material = {
!!      characteristic = 0.0,
!!      relax_velocity = {0.0, 0.0, 0.0},
!!      relax_temperature = 0.0
!!    },
!!    -- Viscous parameters
!!    therm_cond = 0.5,
!!    mu         = 2.0,
!!    ip_param   = 4,
!!    visc_limit = 0
!!  }
!!  equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)
!!```
!!
!! The settings for `isen_coef`, `r`, `material` and `cv` are inherited from the
!! Euler equations, see [[atl_eqn_euler_module]]. Please note, that there are more
!! parameters that may be set as described in [[atl_eqn_euler_module]]
!!
module atl_eqn_nvrstk_module
  ! Treelm modules
  use env_module,             only: rk
  use tem_logging_module,     only: logUnit
  use tem_aux_module,         only: tem_abort

  ! Aotus modules
  use aotus_module,           only: flu_State, aot_get_val
  use aot_out_module,         only: aot_out_type, aot_out_val

  ! Ateles modules
  use atl_eqn_euler_module,   only: atl_euler_type, atl_load_euler, atl_save_euler

  implicit none

  private

  public :: atl_navierStokes_type
  public :: atl_filtNavierStokes_type
  public :: atl_load_navierStokes
  public :: atl_save_navierStokes
  public :: atl_load_filtNS
  public :: atl_Navier_stokes_rans_type
  public :: atl_filNvrStk_source_data_type

  !> The Navier-Stokes equation properties  are stored here
  type :: atl_navierStokes_type

    !> shear viscosity
    real(kind=rk) :: mu

    !> bulk viscosity
    real(kind=rk) :: lambda

    !> thermal conductivity
    real(kind=rk) :: therm_cond

    !> The penalty parameter of the Interior Penalty paramter
    real(kind=rk) :: ip_param

    !> Limiter to decide computation of viscous fluxes within
    !! elements.
    !!
    !! The viscous terms will only be considered if the limiter
    !! exceeds this setting.
    real(kind=rk) :: visc_limit = 0.0_rk

    procedure(invisc_indicator), pointer, pass(nvrstk) :: inviscous => NULL()

  end type atl_navierStokes_type

  !> This data-type stores the properties required
  ! for the reynolds averaged navier stokes type
  type atl_navier_stokes_rans_type

    !> Turbulent Prandtl Number
    real(kind=rk) :: turb_prandtl_num
    ! -------Below the coefficients for the high reynolds number
    ! ------ wilcox k-\omega model are defined ---------------- !
    !> Sigma_k
    real(kind=rk) :: sig_k
    !> \beta_k
    real(kind=rk) :: beta_k
    !> Sigma_omg
    real(kind=rk) :: sig_omg
    !> alpha_omg
    real(kind=rk) :: alpha_omg
    !> \beta_omega
    real(kind=rk) :: beta_omg
    !> \c_{\mu}
    real(kind=rk) :: c_mu
    !> \alpha used in limited turbulent eddy viscosity
    real(kind=rk) :: alpha

  end type atl_navier_stokes_rans_type

  ! This Data-type is to be used as the method data in the eval source
  ! interface routines. The data stored here is required in the evaluation
  ! of certain source terms in RANS
  type atl_filNvrStk_source_data_type
    type(atl_NavierStokes_type), pointer :: nvrStk_type => null()
    type(atl_Navier_stokes_rans_type), pointer :: filNvrStkRans_type => null()
  end type atl_filNvrStk_source_data_type

  !> The Smagorinsky model for LES
  !!
  !! parameters and variables for Smagorinsky model
  !type Smagorinsky_type
  !  real(kind=rk) :: Cs          !< Smagorinsky constant
  !  real(kind=rk) :: Ci          !< Yoshizawa model constant
  !  real(kind=rk) :: prandtl_sgs !< Prandtl number for subgrid-scale (SGS)
  !end type Smagorinsky_type


  !> The Germano-Lilly model for LES
  !type GermanoLilly_type
!> @todo JZ: outdated    !< parameters and variables for Germano-Lilly model
!> @todo JZ: outdated    !> \todo this variable was added because the intel compiler complains about
!> @todo JZ: outdated    !! empty datatypes.
!> @todo JZ: outdated    integer :: n
  !end type GermanoLilly_type

  !> The Filtered Navier-Stokes equation properties are stored here
  type :: atl_FiltNavierStokes_type
!> @todo JZ: outdated     character(len=32)       :: model_type = ''  !< choosing the model for LES
!> @todo JZ: outdated     type(Smagorinsky_type)  :: Smagorinsky
!> @todo JZ: outdated     type(GermanoLilly_type) :: Germano
    ! Choosing the models such as LES/Rans
    character(len=32)       :: model_type = ''
    type(atl_navier_stokes_rans_type) :: rans
  end type atl_FiltNavierStokes_type

  abstract interface
    pure function invisc_indicator(nvrstk, mean, deviation, grad) &
      &   result(isinviscous)
      import :: atl_navierstokes_type, rk
      ! -------------------------------------------------------------------- !
      !> Description of the equation
      class(atl_navierstokes_type), intent(in) :: nvrstk

      !> The mean of each state variable.
      real(kind=rk), intent(in) :: mean(:)

      !> Estimation of maximal deviation of each state.
      real(kind=rk), intent(in) :: deviation(:)

      !> Estimation of maximal gradient of each state.
      real(kind=rk), intent(in) :: grad(:)

      !> Resulting indication whether viscous terms can be neglected.
      logical :: isinviscous
      ! -------------------------------------------------------------------- !
    end function invisc_indicator
  end interface

contains

  ! ------------------------------------------------------------------------ !
  !> Subroutine to initialize an equation of type navier stokes equation
  !! as defined in the configuration file
  subroutine atl_load_navierStokes( navierStokes, euler, conf, eq_table )
    ! -------------------------------------------------------------------- !
    !> Output of this subroutine, the data inside this argument will be set to
    !! the values of the configuration file
    type(atl_navierStokes_type), intent(out) :: navierStokes
    type(atl_euler_type), intent(out) :: euler
    !> Variable of type flu_State which is a flu binding to the configuration
    !! file (input)
    type(flu_State), intent(in) :: conf
    !> A handle to the equation table inside the configuration file
    integer, intent(in) :: eq_table
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! first init euler type
    call atl_load_euler( euler, conf, eq_table)

    ! read shear viscosity from lua file
    call aot_get_val( L       = conf,            &
      &               thandle = eq_table,        &
      &               key     = 'mu',            &
      &               val     = navierStokes%mu, &
      &               ErrCode = iError           )
    if(iError.ne.0) then
      call tem_abort( 'ERROR in init_NavierStokes: no entry for mu ' &
        & // 'found in lua configuration file, stopping...'          )
    end if

    ! set the bulk viscosity to the standard value for an ideal gas phase
    NavierStokes%lambda = 2.0_rk * NavierStokes%mu / 3.0_rk

    call aot_get_val( L       = conf,                    &
      &               thandle = eq_table,                &
      &               key     = 'therm_cond',            &
      &               val     = navierStokes%therm_cond, &
      &               ErrCode = iError                   )
    if(iError.ne.0) then
      call tem_abort( 'ERROR in init_NavierStokes: no entry for therm_cond ' &
        & // 'found in lua configuration file, stopping...'                  )
    end if

    call aot_get_val( L       = conf,                  &
      &               thandle = eq_table,              &
      &               key     = 'ip_param',            &
      &               val     = navierStokes%ip_param, &
      &               ErrCode = iError                 )
    if(iError.ne.0) then
      call tem_abort( 'ERROR in init_NavierStokes: no entry for ip_param ' &
        & // 'found in lua configuration file, stopping...'                )
    end if

    call aot_get_val( L       = conf,                    &
      &               thandle = eq_table,                &
      &               key     = 'visc_limit',            &
      &               val     = navierStokes%visc_limit, &
      &               default = 0.0_rk,                  &
      &               ErrCode = iError                   )

    write(logUnit(1),*) 'Loaded Navier-Stokes equation'
    write(logUnit(1),*) ' * mu: ', navierStokes%mu
    write(logUnit(1),*) ' * therm_cond: ', navierStokes%therm_cond
    write(logUnit(1),*) ' * ip_param: ', navierStokes%ip_param

    if (navierStokes%visc_limit > 0.0_rk) then
      write(logUnit(1),*) ' Viscous terms within elements will only be'
      write(logUnit(1),*) ' considered if the limit of ', &
        &                 navierStokes%visc_limit
      write(logUnit(1),*) ' is exceeded by the indicator.'
      write(logUnit(1),*) ' An estimate for the maximal gradient in all'
      write(logUnit(1),*) ' elements will be computed.'
    end if

  end subroutine atl_load_navierStokes
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! dump the equation variables into the lua file
  subroutine atl_save_navierStokes(me_euler, me_nvrstk, eqn_name, conf)
    ! -------------------------------------------------------------------- !
    type(atl_euler_type), intent(in) :: me_euler
    type(atl_navierStokes_type), intent(in) :: me_nvrstk
    character(len=*), intent(in) :: eqn_name
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !

    call aot_out_val( put_conf = conf,          &
      &               vname    = 'name',        &
      &               val      = trim(eqn_name) )

    ! Dump Euler properties
    call atl_save_euler( me       = me_euler, &
      &                  eqn_name = eqn_name, &
      &                  conf     = conf,     &
      &                  withName = .false.   )

    ! Dump Navier-Stokes properties
    call aot_out_val( put_conf = conf,        &
      &               vname    = 'mu',        &
      &               val      = me_nvrstk%mu )
    call aot_out_val( put_conf = conf,                &
      &               vname    = 'therm_cond',        &
      &               val      = me_nvrstk%therm_cond )
    call aot_out_val( put_conf = conf,              &
      &               vname    = 'ip_param',        &
      &               val      = me_nvrstk%ip_param )
    call aot_out_val( put_conf = conf,                &
      &               vname    = 'visc_limit',        &
      &               val      = me_nvrstk%visc_limit )


  end subroutine atl_save_navierStokes
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> subroutine to initialize an equation of type filtered-navier-stokes
  !! equation (turbulence modelling) as defined in the configuration file
  !!
  subroutine atl_load_filtNS( filtNavierStokes, navierStokes, euler, conf, &
      &                       eq_table)
    ! -------------------------------------------------------------------- !
    !> Output of this subroutine, the data inside this argument will be set to
    !! the values of the configuration file
    type(atl_FiltNavierStokes_type), intent(out) :: filtNavierStokes
    type(atl_NavierStokes_type), intent(out) :: navierStokes
    type(atl_Euler_type), intent(out) :: euler
    !> Variable of type flu_State which is a flu binding to the configuration
    !! file
    type(flu_State), intent(in) :: conf
    !> A handle to the equation table inside the configuration file
    integer, intent(in) :: eq_table
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    call atl_load_navierstokes( navierStokes = navierStokes, &
      &                         euler        = euler,        &
      &                         conf         = conf,         &
      &                         eq_table     = eq_table      )

    call aot_get_val( L       = conf,                        &
      &               thandle = eq_table,                    &
      &               key     = 'turbulence_model',          &
      &               val     = filtNavierStokes%model_type, &
      &               ErrCode = iError)

    write(logUnit(1),*) ' * using turbulence model: ', filtNavierStokes%model_type


  end subroutine atl_load_filtNS
  ! ------------------------------------------------------------------------ !

end module atl_eqn_nvrstk_module
