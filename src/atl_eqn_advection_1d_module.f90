! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> A module describing the adection equation system in 1D.
module atl_eqn_advection_1d_module
  use env_module,     only: rk
  use aotus_module,   only: flu_State, aot_get_val
  use aot_out_module, only: aot_out_type, aot_out_val

  implicit none

  private

  public :: atl_advection_1d_type
  public :: atl_load_advection_1d
  public :: atl_save_advection_1d

  type atl_advection_1d_type
    real(kind=rk) :: velocity  !< advection velocity
  end type atl_advection_1d_type

contains

  !> subroutine to initialize an equation of type advection equation
  !! as defined in the configuration file
  subroutine atl_load_advection_1d(adv, conf, eq_table)
    ! --------------------------------------------------------------------------!
    !> Resulting description of the advection equation parameters.
    type(atl_advection_1d_type), intent(out) :: adv

    !> Handle to the configuration script, to load the parameters from.
    type(flu_State) :: conf

    !> Handle to the table containing the description for the equation
    !! system.
    integer, intent(in) :: eq_table
    ! --------------------------------------------------------------------------!
    integer :: iError
    ! --------------------------------------------------------------------------!

    !read the data from the equation table of the lua file
    call aot_get_val(L = conf, thandle = eq_table, key = 'velocity', &
      &              val = adv%velocity, &
      &              ErrCode = iError)

  end subroutine atl_load_advection_1d

  ! dump the equation variables into the lua file
  subroutine atl_save_advection_1d(me, eqn_name, conf)
    ! --------------------------------------------------------------------------
    type(atl_advection_1d_type), intent(in) :: me
    character(len=*), intent(in) :: eqn_name
    type(aot_out_type) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_val( put_conf = conf, vname = 'name', val = trim(eqn_name) )

    ! Dump equation Properties
    call aot_out_val( put_conf = conf, vname = 'velocity', val = me%velocity  )

  end subroutine atl_save_advection_1d

end module atl_eqn_advection_1d_module
