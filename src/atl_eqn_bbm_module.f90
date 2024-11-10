! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Module describing a simple, basic membrane model
module atl_eqn_bbm_module
  use env_module,         only: rk
  use tem_logging_module, only: logUnit
  use aotus_module,       only: flu_State, aot_get_val

  implicit none

  private

  public :: atl_BBMEM_type
  public :: atl_load_BBMEM

  !> Parameters describing the model
  type atl_BBMEM_type
    real(kind=rk) :: tna
    real(kind=rk) :: tcl
    real(kind=rk) :: kappa
  end type atl_BBMEM_type

contains

  !> summary: subroutine to intialize BBM
  !! material parameters are to be filled here.
  subroutine atl_load_BBMEM(BBM, conf, eq_table)
    ! --------------------------------------------------------------------------!
    type(atl_BBMEM_type), intent(out) :: BBM
    type(flu_State)             :: conf
    integer, intent(in)         :: eq_table
    ! --------------------------------------------------------------------------!
    integer :: iError
    ! --------------------------------------------------------------------------!

    write(logUnit(1),*) 'Initializing BBM'
    call aot_get_val(L = conf, thandle = eq_table, key = 'tna', &
      &              val = BBM%tna, &
      &              ErrCode = iError)
    write(logUnit(1),*) '  tna  : ', BBM%tna

    call aot_get_val(L = conf, thandle = eq_table, key = 'tcl', &
      &              val = BBM%tcl, &
      &              ErrCode = iError)
    write(logUnit(1),*) '  tcl  : ', BBM%tcl

    call aot_get_val(L = conf, thandle = eq_table, key = 'kappa', &
      &              val = BBM%kappa, &
      &              ErrCode = iError)
    write(logUnit(1),*) '  kappa: ', BBM%kappa

  end subroutine atl_load_BBMEM

end module atl_eqn_bbm_module
