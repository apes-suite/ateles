! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Module to describe the Nernst-Planck equation system.
module atl_eqn_nerplanck_module
  use env_module,     only: rk

  use aotus_module,   only: flu_State, &
    &                       aot_get_val

  implicit none

  private

  public :: atl_nernstPlanck_type
  public :: atl_load_NernstPlanck

  !> Datatype for Nernst-Planck equations.
  !!
  !! This datatype describes the time dependent Nernst-Planck equations of the
  !! following form (where \( D \) is the diffusivity):
  !!
  !! Concentration:
  !! \( \partial_t u = \sqrt{(D)} \nabla \times \rho \)
  !!
  !! Diffusive flux:
  !! \( \rho = \sqrt{(D)} \nabla \times u \)
  !!
  !! The diffusivity is assumed constant. \( u\in\mathbb{R} \) is the
  !! concentration and \( \rho\in\mathbb{R}^3 \) is the diffusive flux
  !! vector field.
  type atl_nernstPlanck_type
    !> diffusivity in SI units. This is assumed to be constant
    !! over the whole domain.
    real(kind=rk)     :: D
  end type atl_nernstPlanck_type


contains

  !> summary: subroutine to intialize Nernst-Planck equation with constant
  !! diffusivity.
  subroutine atl_load_nernstPlanck(nerplanck, conf, eq_table)
    ! --------------------------------------------------------------------------
    type(atl_nernstPlanck_type), intent(out) :: nerplanck
    type(flu_State) :: conf
    integer, intent(in) :: eq_table
    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------

    call aot_get_val( L       = conf,        &
      &               thandle = eq_table,    &
      &               key     = 'D',         &
      &               val     = nerplanck%D, &
      &               ErrCode = iError       )

  end subroutine atl_load_nernstPlanck

end module atl_eqn_nerplanck_module
