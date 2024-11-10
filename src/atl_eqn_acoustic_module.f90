! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

!> Module to describe the isothermal acoustic equation system for a medium at
!! rest.
!!
!! This is a simplification of the linearized Euler equations, where the medium
!! is assumed to be at rest and isothermal. These assumptions eliminate one
!! equation and only density and velocity are required to describe the wave
!! transport. Pressure and density are equivalent in this system and differ only
!! by a constant scaling factor.
!!
!! The medium is defined by its pressure and density, which define the speed
!! of sound. They are referred to as background state, and accordingly a
!! `background` table needs to be defined for the acoustic equation table.
!!
!!```lua
!!  equation = {
!!    name = 'acoustic',
!!    background = {
!!      density = 1,
!!      pressure = 1
!!    }
!!  }
!!```
!!
!! The speed of sound \(c\) derives from the provided background state by
!!
!! \[ c = \sqrt{p/\rho} \]
!!
!! Where \(p\) is the background pressure in the medium and \(\rho\) its
!! density.
!!
!! The acoustic equations are available for 3D (`acoustic`) and 2D
!! (`acoustic_2d`).
module atl_eqn_acoustic_module
  ! Treelm modules
  use env_module,             only: rk

  ! Aotus modules
  use aotus_module,           only: flu_State, aot_get_val
  use aot_out_general_module, only: aot_out_type
  use aot_out_module,         only: aot_out_val
  use aot_table_module,       only: aot_table_open, &
    &                               aot_table_close
  ! Treelm modules
  use tem_tools_module,       only: tem_horizontalSpacer
  use tem_logging_module,     only: logUnit

  ! Ateles modules
  use atl_materialFun_module, only: atl_materialFun_type

  implicit none

  private

  public :: atl_acoustic_type
  public :: atl_load_acoustic
  public :: atl_save_acoustic
  public :: atl_dump_acoustic_eqn

  !> Datatype for acoustic equations.
  !!
  type atl_acoustic_type
    !> Number of spatial dimensions to consider.
    !! Implemented are 2 and 3 dimensions.
    integer :: ndims

    !> Density of the wave carrying medium.
    real(kind=rk) :: density_0

    !> Velocity of the medium.
    real(kind=rk),allocatable :: velocity_0(:)

    !> Pressure of the medium.
    real(kind=rk) :: pressure_0

    !> Speed of sound of the medium.
    real(kind=rk) :: speedOfSound

    !> Penalization terms to describe obstacles.
    type(atl_materialFun_type) :: penalization
  end type atl_acoustic_type


contains


  ! ------------------------------------------------------------------------ !
  !> Load the configuration for acoustic equations from the Lua script.
  subroutine atl_load_acoustic( acoustic, conf, eq_table )
    ! -------------------------------------------------------------------- !
    type(atl_acoustic_type), intent(inout) :: acoustic
    type(flu_State)                 :: conf
    integer, intent(in)             :: eq_table
    ! -------------------------------------------------------------------- !
    integer                       :: iError, acoustic_table
    ! -------------------------------------------------------------------- !
    ! allocate the dimension of background velocity array according to dimension
    allocate(acoustic%velocity_0(acoustic%ndims))

    call tem_horizontalSpacer(funit=logUnit(5))
    write(logUnit(5),*) 'Reading background parameters for Acoustic equation . '

    ! Open subtable for backgroudn properties
    call aot_table_open( L       = conf,           &
      &                  parent  = eq_table,       &
      &                  tHandle = acoustic_table, &
      &                  key     = 'background'    )

    if(acoustic_table.eq.0) then
      write(logUnit(1),*) 'ERROR in init_acoustic: no background properties &
        &                  defined, stopping ...'
    end if

    !read the data from the equation table of the lua file
    call aot_get_val( l       = conf,               &
      &               thandle = acoustic_table,     &
      &               key     = 'density',          &
      &               val     = acoustic%density_0, &
      &               errcode = ierror              )

    acoustic%velocity_0 = 0.0

    call aot_get_val( l       = conf,                &
      &               thandle = acoustic_table,      &
      &               key     = 'pressure',          &
      &               val     = acoustic%pressure_0, &
      &               errcode = ierror               )

    ! c = sqrt( p/rho - 1/2 * v^2)
    acoustic%SpeedOfSound = sqrt( acoustic%pressure_0 / acoustic%density_0)


    ! Close the Lua table with the background information
    call aot_table_close( L = conf, tHandle = acoustic_table )

  end subroutine atl_load_acoustic
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> dump the equation variables into the lua file
  subroutine atl_save_acoustic(me, eqn_name, conf)
    ! -------------------------------------------------------------------- !
    type(atl_acoustic_type), intent(in) :: me
    character(len=*), intent(in) :: eqn_name
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    call aot_out_val(put_conf= conf, vname = 'name', val = trim(eqn_name))

    ! Dump equation Properties
    call aot_out_val( put_conf = conf,        &
      &               vname    = 'density',   &
      &               val      = me%density_0 )

    call aot_out_val( put_conf = conf,         &
      &               vname    = 'pressure',   &
      &               val      = me%pressure_0 )

  end subroutine atl_save_acoustic
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! write the equation variables into a file represented by dumpUnit.
  subroutine atl_dump_acoustic_eqn( me, dumpUnit )
    ! -------------------------------------------------------------------- !
    type(atl_acoustic_type), intent(in) :: me
    integer :: dumpUnit
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

     write(dumpUnit,*) 'Background variable for acoustic equation.'
     write(dumpUnit,*) 'Background Density =', me%density_0
     write(dumpUnit,*) 'Background Pressure =', me%pressure_0

  end subroutine atl_dump_acoustic_eqn
  ! ------------------------------------------------------------------------ !

end module atl_eqn_acoustic_module
