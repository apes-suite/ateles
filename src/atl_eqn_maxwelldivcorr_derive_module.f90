! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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

!> Module for derived variables of the Maxwell equations with hyperbolic divergence
!! cleaning.
module atl_eqn_maxwelldivcorr_derive_module

  ! Treelm modules
  use env_module,             only: rk

  ! Ateles modules
  use atl_equation_module,    only: atl_equations_type

  implicit none
  private

  public :: atl_eqn_maxwelldivcorr_cons2prim
  public :: atl_eqn_maxwelldivcorr_prim2cons

contains

  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_maxwelldivcorr_prim2cons( equation, instate, outstate, &
    &                                          material                     )
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional, intent(in) :: material(:,:)
    ! --------------------------------------------------------------------------!

    if( size(material,1).eq.1 ) then
      if(present(outstate)) then
        ! displacement field, i.e. D from
        outstate(:,1) = instate(:,1) * material(1,1)
        outstate(:,2) = instate(:,2) * material(1,1)
        outstate(:,3) = instate(:,3) * material(1,1)

        ! magnetic field, i.e. B from H
        outstate(:,4) = instate(:,4) * material(1,2)
        outstate(:,5) = instate(:,5) * material(1,2)
        outstate(:,6) = instate(:,6) * material(1,2)

        ! corrections are conservative and primitive
        outstate(:,7) = instate(:,7)
        outstate(:,8) = instate(:,8)
      else
        ! displacement field, i.e. D from
        instate(:,1) = instate(:,1) * material(1,1)
        instate(:,2) = instate(:,2) * material(1,1)
        instate(:,3) = instate(:,3) * material(1,1)

        ! magnetic field, i.e. B from H
        instate(:,4) = instate(:,4) * material(1,2)
        instate(:,5) = instate(:,5) * material(1,2)
        instate(:,6) = instate(:,6) * material(1,2)

        ! corrections are conservative and primitive
        instate(:,7) = instate(:,7)
        instate(:,8) = instate(:,8)
      end if
    else
      if(present(outstate)) then
        ! displacement field, i.e. D from
        outstate(:,1) = instate(:,1) * material(:,1)
        outstate(:,2) = instate(:,2) * material(:,1)
        outstate(:,3) = instate(:,3) * material(:,1)

        ! magnetic field, i.e. B from H
        outstate(:,4) = instate(:,4) * material(:,2)
        outstate(:,5) = instate(:,5) * material(:,2)
        outstate(:,6) = instate(:,6) * material(:,2)

        ! corrections are conservative and primitive
        outstate(:,7) = instate(:,7)
        outstate(:,8) = instate(:,8)
      else
        ! displacement field, i.e. D from
        instate(:,1) = instate(:,1) * material(:,1)
        instate(:,2) = instate(:,2) * material(:,1)
        instate(:,3) = instate(:,3) * material(:,1)

        ! magnetic field, i.e. B from H
        instate(:,4) = instate(:,4) * material(:,2)
        instate(:,5) = instate(:,5) * material(:,2)
        instate(:,6) = instate(:,6) * material(:,2)

        ! corrections are conservative and primitive
        instate(:,7) = instate(:,7)
        instate(:,8) = instate(:,8)
      end if
    end if

  end subroutine atl_eqn_maxwelldivcorr_prim2cons


  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_maxwelldivcorr_cons2prim( equation, instate, outstate, &
    &                                          material                     )
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout) :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in) :: material(:,:)
    ! --------------------------------------------------------------------------!

    if(size(material,1).eq.1) then
      if(present(outstate)) then
        ! electric field, i.e. E from D
        outstate(:,1) = instate(:,1) / material(1,1)
        outstate(:,2) = instate(:,2) / material(1,1)
        outstate(:,3) = instate(:,3) / material(1,1)

        ! magnetizing field, i.e. H from B
        outstate(:,4) = instate(:,4) / material(1,2)
        outstate(:,5) = instate(:,5) / material(1,2)
        outstate(:,6) = instate(:,6) / material(1,2)

        ! corrections are conservative and primitive
        outstate(:,7) = instate(:,7)
        outstate(:,8) = instate(:,8)
      else
        ! electric field, i.e. E from D
        instate(:,1) = instate(:,1) / material(1,1)
        instate(:,2) = instate(:,2) / material(1,1)
        instate(:,3) = instate(:,3) / material(1,1)

        ! magnetizing field, i.e. H from B
        instate(:,4) = instate(:,4) / material(1,2)
        instate(:,5) = instate(:,5) / material(1,2)
        instate(:,6) = instate(:,6) / material(1,2)

        ! corrections are conservative and primitive
        instate(:,7) = instate(:,7)
        instate(:,8) = instate(:,8)
      end if
    else
      if(present(outstate)) then
        ! electric field, i.e. E from D
        outstate(:,1) = instate(:,1) / material(:,1)
        outstate(:,2) = instate(:,2) / material(:,1)
        outstate(:,3) = instate(:,3) / material(:,1)

        ! magnetizing field, i.e. H from B
        outstate(:,4) = instate(:,4) / material(:,2)
        outstate(:,5) = instate(:,5) / material(:,2)
        outstate(:,6) = instate(:,6) / material(:,2)

        ! corrections are conservative and primitive
        outstate(:,7) = instate(:,7)
        outstate(:,8) = instate(:,8)
      else
        ! electric field, i.e. E from D
        instate(:,1) = instate(:,1) / material(:,1)
        instate(:,2) = instate(:,2) / material(:,1)
        instate(:,3) = instate(:,3) / material(:,1)

        ! magnetizing field, i.e. H from B
        instate(:,4) = instate(:,4) / material(:,2)
        instate(:,5) = instate(:,5) / material(:,2)
        instate(:,6) = instate(:,6) / material(:,2)

        ! corrections are conservative and primitive
        instate(:,7) = instate(:,7)
        instate(:,8) = instate(:,8)
      end if
    end if
  end subroutine atl_eqn_maxwelldivcorr_cons2prim


end module atl_eqn_maxwelldivcorr_derive_module
