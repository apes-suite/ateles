! Copyright (c) 2013-2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Module to describe the Maxwell equation system for electrodynamics.
!!
!! This implements the equation system for electric-magnetic fields of the
!! following form
!! Change of (electric) displacement field \(D = \epsilon E\):
!! \[ \partial_t D = \nabla \times \frac{B}{\mu} - J \]
!!
!! Change of magnetic field \(B\):
!! \[ \partial_t B  = - \nabla \times \frac{D}{\epsilon} \]
!!
!! Nonexistence of magnetic monopoles:
!! \[  \nabla \cdot B = 0 \]
!!
!! Sources of electric field:
!! \[ \nabla \cdot D = \rho \]
!!
!! Where \( \rho \) is the spatial charge density and \(J = \sigma E \) the free
!! The magnetic `permeability` \( \mu \), the electric
!! `permittivity` \( \epsilon \) and the conductivity \( \sigma \) need to be
!! provided as variables in a material table.
!! See [[tem_variable_module]] for details on the definition of variables.
!! In their simplest form they may just be scalars.
!! Thus, a typical definition for the Maxwell equation system takes the
!! following form:
!!
!!```lua
!!  equation = {
!!    name = 'maxwell',
!!    material = {
!!      permeability = 1.0,
!!      permittivity = 1.0,
!!      conductivity = 0.0
!!    }
!!  }
!!```
!!
!! There are various variants of the Maxwell equations available in Ateles:
!!
!! * `maxwell`
!! * `maxwell_2d`
!! * `maxwelldivcorr`
!!
!! For the divergence correction conductivity is not available. Instead the
!! parameters `gam` and `chi` have to be provided for the divergence correction.
!!
module atl_eqn_maxwell_module
  ! Aotus modules
  use aot_out_module, only: aot_out_type, &
    &                       aot_out_val

  implicit none

  private

  public :: atl_maxwell_type
  public :: atl_save_maxwell

  !> Datatype for Maxwell equations.
  !!
  !! This datatype describes the time dependent Maxwell equations of the
  !! following form (where \( \rho \)  is the spatial charge density and
  !! \( J \) is the free current (not including bound current) density ):
  !!
  !! (Electric) discplacement field:
  !! \[ \partial_t D = \nabla \times \frac{B}{\mu}
  !!                  - J \]
  !!
  !! Magnetic field: \( \partial_t B  = - \nabla \times \frac{D}{\epsilon} \)
  !!
  !! Nonexistence of magnetic monopoles: \(  \nabla \cdot B = 0 \)
  !!
  !! Sources of electric field: \( \nabla \cdot D = \rho \)
  !!
  !! The magnetic permeability \( \mu \)  and the electric
  !! \( \epsilon \) permitivity are assumed constant or variable.
  !! D is the (electric) displacement field
  !! vector and B is the magnetic field vector (also called magnetic induction).
  type atl_maxwell_type
  end type atl_maxwell_type


contains


  ! ------------------------------------------------------------------------ !
  ! Dump the equation variables into the lua file
  subroutine atl_save_maxwell(eqn_name, conf)
    ! -------------------------------------------------------------------- !
    character(len=*), intent(in) :: eqn_name
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !

    call aot_out_val( put_conf = conf, vname = 'name', val = trim(eqn_name) )

    ! @todo PV 20150727: Reactive this call
    !call atl_dump_materialFun( conf, me%permeaPermit, 'material' )

  end subroutine atl_save_maxwell
  ! ------------------------------------------------------------------------ !

end module atl_eqn_maxwell_module
