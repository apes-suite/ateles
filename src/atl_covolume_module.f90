! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> author: Jens Zudrop
!! Module containing routines and datatypes for the covolume
!! stabilization.
module atl_covolume_module
  use env_module,                     only: rk, zero_rk

  use aotus_module,                   only: flu_State, aot_get_val

  use tem_tools_module,               only: tem_horizontalSpacer
  use tem_aux_module,                 only: tem_abort
  use tem_logging_module,             only: logUnit

  use atl_spectral_viscosity_module,  only: atl_exp_spectral_visc_prp, &
    &                                       atl_poly_spectral_visc_prp

  implicit none
  private

  !> Datatype representing the covolume filter stabilization
  type atl_covolume_type
    real(kind=rk) :: alpha = 0.0_rk
    real(kind=rk) :: order = 0.0_rk
    real(kind=rk) :: cut_order = -1.0_rk
    real(kind=rk) :: beta = 0.0_rk
    integer :: kind = atl_exp_spectral_visc_prp
    logical :: isAdaptive = .false.
    real(kind=rk) :: recovery_order = 1.0_rk
    real(kind=rk) :: recovery_density = 1e-2_rk
    real(kind=rk) :: recovery_pressure = 1e-2_rk
    real(kind=rk) :: recovery_mach = 2.5_rk
  end type atl_covolume_type

  public :: atl_ini_covolume, atl_covolume_type

contains

  !> Subroutine to load configuration file options for the
  !! covolume filter method.
  subroutine atl_ini_covolume(conf, parent_table, filter)
    ! --------------------------------------------------------------------------
    !> flu binding to lua configuration file.
    type(flu_State), intent(in) :: conf
    !> The parent table
    integer, intent(in) :: parent_table
    !> The filter to initialize
    type(atl_covolume_type), intent(out) :: filter
    ! --------------------------------------------------------------------------
    integer :: iError
    character(len=128) :: filterkind
    ! --------------------------------------------------------------------------

    call aot_get_val(L       = conf,         &
      &              thandle = parent_table, &
      &              key     = 'alpha',      &
      &              val     = filter%alpha, &
      &              default = 0.0_rk,       &
      &              ErrCode = iError        )

    call aot_get_val(L       = conf,         &
      &              thandle = parent_table, &
      &              key     = 'order',      &
      &              val     = filter%order, &
      &              default = 0.0_rk,       &
      &              ErrCode = iError        )

    if (filter%order > zero_rk) then
      write(logUnit(2),*) 'Using spectral viscosity for covolumes.'
      if (filter%alpha <= zero_rk) then
        write(logUnit(1),*) 'NOTE: No parameter alpha provided for spectral' &
          &                 // ' viscosity in covolumes.'
        write(logunit(1),*) '      Using default of 36.'
        filter%alpha = 36._rk
      end if
    else
      if (filter%alpha > zero_rk) then
        write(logunit(1),*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(logUnit(1),*) 'ERROR in atl_ini_covolume: Alpha for spectral' &
          &                 // ' viscosity provided,'
        write(logUnit(1),*) '       but no order given to use for the filter!'
        write(logUnit(1),*) '       Provide a order to activate spectral' &
          &                 // ' viscosity in covolumes,'
        write(logUnit(1),*) '       or remove the definition of alpha to' &
          &                 // ' deactivate it.'
        write(logUnit(1),*) ' ! STOPPING ! '
        write(logunit(1),*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        call tem_abort()
      end if
    end if

    if (filter%alpha > zero_rk) then
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
    end if

    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'beta', &
      &              val = filter%beta, &
      &              default = 0.0_rk, &
      &              ErrCode = iError)
    if (filter%alpha > zero_rk) then
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
    end if

    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(1),*) 'Covolume filter parameter:'
    write(logUnit(1),*) '* Beta: ', filter%beta
    if (filter%alpha > zero_rk) then
      write(logUnit(1),*) '* Spectral filter on covolumes:'
      write(logUnit(1),*) '  + Order: ', filter%order
      write(logUnit(1),*) '  + Alpha: ', filter%alpha
      write(logUnit(1),*) '  + Cut Order: ', filter%cut_order
      write(logUnit(1),*) '  + Kind: ', filter%kind
      write(logUnit(1),*) '  + Is adaptive: ', filter%isAdaptive
      write(logUnit(1),*) '  + Recovery Order: ', filter%recovery_order
    else
      write(logUnit(1),*) '* NO Spectral filter on covolumes'
    end if
    call tem_horizontalSpacer(fUnit=logUnit(1))

  end subroutine atl_ini_covolume

end module atl_covolume_module
