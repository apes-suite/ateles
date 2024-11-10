! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
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
!!
!! Module containing routines and datatypes for the positivity
!! preserving scheme.
!!
!! The limiter is a conservative, high order, positivity
!! preserving limiter described in:
!! Zhang, X., & Shu, C.-W. (2010).
!! On positivity-preserving high order discontinuous
!! Galerkin schemes for compressible Euler equations
!! on rectangular meshes.
!! Journal of Computational Physics, 229(23), 8918–8934.
!! doi:10.1016/j.jcp.2010.08.016
module atl_cons_positivity_preserv_module
  use env_module,               only: rk

  use aotus_module,             only: flu_State, aot_get_val

  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit

  implicit none
  private

  !> Datatype representing the positivity preserving
  !! limiter. Should be applied after all other limiters
  !! are applied and after each step of a SSP timeintegrator.
  type atl_cons_positivity_preserv_type
    !> Smallness parameter. If denisty or pressure are below this
    !! value, the state is considered as unphysical.
    real(kind=rk) :: eps = 1.0e-13
  end type atl_cons_positivity_preserv_type

  public :: atl_ini_cons_positivity_preserv, atl_cons_positivity_preserv_type

contains

  !> Subroutine to load configuration file options for the conservative
  !! positivity preserving limiter.
  subroutine atl_ini_cons_positivity_preserv(conf, parent_table, filter)
    ! --------------------------------------------------------------------------
    !> flu binding to lua configuration file.
    type(flu_State), intent(in) :: conf
    !> The parent table
    integer, intent(in) :: parent_table
    !> The filter to initialize
    type(atl_cons_positivity_preserv_type), intent(out) :: filter
    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------

    call aot_get_val(L = conf, thandle = parent_table, &
      &              key = 'eps', &
      &              val = filter%eps, &
      &              ErrCode = iError)
    if(iError.ne.0) then
      write(logUnit(1),*) 'ERROR in atl_ini_cons_positivity_preserv: not able to read order ' // &
        & 'parameter from configuration file, stopping ...'
      call tem_abort()
    end if

    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(1),*) 'Conservative positivity preserving parameter:'
    write(logUnit(1),*) '* Eps: ', filter%eps
    call tem_horizontalSpacer(fUnit=logUnit(1))

  end subroutine atl_ini_cons_positivity_preserv


end module atl_cons_positivity_preserv_module
