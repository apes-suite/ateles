! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
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
!! Module containing routines and datatypes for the spectral viscosity
!! stabilization.
module atl_spectral_viscosity_module
  use env_module,               only: rk

  use aotus_module,             only: flu_State, aot_get_val

  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit

  implicit none
  private

  integer, parameter :: atl_exp_spectral_visc_prp = 1
  integer, parameter :: atl_poly_spectral_visc_prp = 2

  !> Datatype representing the spectral viscosity stabilization.
  !!
  !! The spectral viscosity stabilzation applies a damping of the
  !! modal coefficients (in spectral space).
  !! The filter applies an expontential damping:
  !!
  !! If the modal expansion is given by:
  !! \( u(x) = \sum_{k=0}^N u_k L_k(x) \)
  !!
  !! The spectral viscosity applies the following damping:
  !! \( u(x) = \sum_{k=0}^N \sigma_k u_k L_k(x) \)
  !!
  !! The damping coefficient \( \sigma_k \) is given by:
  !! \( \sigma_k = exp( -\alpha (|k|_{L2}/|cut_order|_{L2})^{order}  )\)
  type atl_spectral_visc_type
    real(kind=rk) :: alpha = 0.0_rk
    real(kind=rk) :: order = 0.0_rk
    real(kind=rk) :: cut_order = -1.0_rk
    integer :: kind = atl_exp_spectral_visc_prp
    logical :: isAdaptive = .false.
    real(kind=rk) :: recovery_order = 1.0_rk
    real(kind=rk) :: recovery_density = 1e-2_rk
    real(kind=rk) :: recovery_pressure = 1e-2_rk
    real(kind=rk) :: recovery_mach = 2.5_rk
  end type atl_spectral_visc_type

  public :: atl_ini_spectral_visc, atl_spectral_visc_type
  public :: atl_exp_spectral_visc_prp, atl_poly_spectral_visc_prp

contains

  !> Subroutine to load configuration file options for the spectral
  !! viscosity method.
  subroutine atl_ini_spectral_visc(conf, parent_table, filter)
    ! --------------------------------------------------------------------------
    !> flu binding to lua configuration file.
    type(flu_State), intent(in) :: conf
    !> The parent table
    integer, intent(in) :: parent_table
    !> The filter to initialize
    type(atl_spectral_visc_type), intent(out) :: filter
    ! --------------------------------------------------------------------------
    integer :: iError
    character(len=128) :: filterkind
    ! --------------------------------------------------------------------------

    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'alpha', &
      &              val = filter%alpha, &
      &              ErrCode = iError)
    if(iError.ne.0) then
      write(logUnit(1),*) 'ERROR in atl_ini_spectral_visc: not able to read alpha ' // &
        & 'parameter from configuration file, stopping ...'
      call tem_abort()
    end if

    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'order', &
      &              val = filter%order, &
      &              ErrCode = iError)
    if(iError.ne.0) then
      write(logUnit(1),*) 'ERROR in atl_ini_spectral_visc: not able to read order ' // &
        & 'parameter from configuration file, stopping ...'
      call tem_abort()
    end if

    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'isAdaptive', &
      &              val = filter%isAdaptive, &
      &              default = .false., &
      &              ErrCode = iError)
    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'recovery_order', &
      &              val = filter%recovery_order, &
      &              default = 1.0_rk, &
      &              ErrCode = iError)
    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'recovery_density', &
      &              val = filter%recovery_density, &
      &              default = 1e-2_rk, &
      &              ErrCode = iError)
    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'recovery_pressure', &
      &              val = filter%recovery_pressure, &
      &              default = 1e-2_rk, &
      &              ErrCode = iError)
    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'recovery_mach', &
      &              val = filter%recovery_mach, &
      &              default = 2.5_rk, &
      &              ErrCode = iError)

    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'cut_order', &
      &              val = filter%cut_order, &
      &              default = -1.0_rk, &
      &              ErrCode = iError)

    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'kind', &
      &              val = filterkind, &
      &              default = 'exp', &
      &              ErrCode = iError)

    select case(filterkind)
    case('exp')
      filter%kind = atl_exp_spectral_visc_prp
    case('poly')
      filter%kind = atl_poly_spectral_visc_prp
    case default
      write(logUnit(1),*) 'Unknown spectral filter kind: Select exp or poly, stopping now ...'
      call tem_abort()
    end select

    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(1),*) 'Spectral viscosity parameter:'
    write(logUnit(1),*) '* Alpha: ', filter%alpha
    write(logUnit(1),*) '* Order: ', filter%order
    write(logUnit(1),*) '* Cut Order: ', filter%cut_order
    write(logUnit(1),*) '* Kind: ', filter%kind
    write(logUnit(1),*) '* Is adaptive: ', filter%isAdaptive
    write(logUnit(1),*) '* Recovery Order: ', filter%recovery_order
    call tem_horizontalSpacer(fUnit=logUnit(1))

  end subroutine atl_ini_spectral_visc


end module atl_spectral_viscosity_module
