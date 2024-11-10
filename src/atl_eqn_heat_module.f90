! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> A module describing the heat equation.
!!
!! The heat equation is the simple scalar equation that utilizes the Laplacian
!! operator in space:
!! \[ \frac{\partial u}{\partial t} = k \Delta u \]
!!
!! Where \(u\) represents the temperature.
!! The only parameter that needs to be configured for this equation is the
!! diffusivity parameter \(k\).
!! Thus the equation table for the heat equation takes the following form:
!!
!!```lua
!!  equation = {
!!    name = 'heat',
!!    k = 1
!!  }
!!```
!!
!! The heat equation is implemented for 1, 2 and 3 dimensions, which are
!! configured by choosing the corresponding name:
!!
!! * `heat` (3D)
!! * `heat_2d` (2D)
!! * `heat_1d` (1D)
module atl_eqn_heat_module
  use env_module,             only: rk
  use aotus_module,           only: flu_State, aot_get_val

  use aot_out_module,         only: aot_out_type, aot_out_val

  ! Ateles modules
  use atl_materialFun_module, only: atl_materialFun_type

  implicit none

  private

  public :: atl_heat_type
  public :: atl_load_heat
  public :: atl_save_heat

  type atl_heat_type
    !> Thermal diffusivity of the medium.
    real(kind=rk) :: k

    !> Penalization terms
    type(atl_materialFun_type) :: penalization
  end type atl_heat_type


contains


  ! ------------------------------------------------------------------------ !
  !> subroutine to initialize an equation of type heat equation
  !! as defined in the configuration file
  subroutine atl_load_heat( heat, conf, eq_table )
    ! -------------------------------------------------------------------- !
    !> Resulting description of the heat equation parameters.
    type(atl_heat_type), intent(out) :: heat

    !> Handle to the configuration script, to load the parameters from.
    type(flu_State)               :: conf

    !> Handle to the table containing the description for the equation
    !! system.
    integer, intent(in)           :: eq_table
    ! -------------------------------------------------------------------- !
    integer                       :: iError
    ! -------------------------------------------------------------------- !

    !read the data from the equation table of the lua file
    call aot_get_val(L = conf, thandle = eq_table, key = 'k', &
      &              val = heat%k, &
      &              ErrCode = iError)

  end subroutine atl_load_heat
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! dump the equation variables into the lua file
  subroutine atl_save_heat(me, eqn_name, conf)
    ! -------------------------------------------------------------------- !
    type(atl_heat_type), intent(in) :: me
    character(len=*), intent(in) :: eqn_name
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !

    call aot_out_val(put_conf = conf, vname = 'name', val = trim(eqn_name))

    ! Dump equation Properties
    call aot_out_val(put_conf = conf, vname = 'k', val = me%k)

  end subroutine atl_save_heat
  ! ------------------------------------------------------------------------ !

end module atl_eqn_heat_module
