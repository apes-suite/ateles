! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013, 2017-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

!> author: Jens Zudrop
!! author Peter Vitt 2014, 2017, 2018
!! Module collecting all data types and subroutines related to the space basis.
!!
!! The polynomial basis defined here is used for the cubic modal discontinuous
!! Galerkin (CuMoDiG) scheme and the reconstructed polynomials in the Finite
!! Volume schemes within Ateles.
!! It makes use of a module variable for the basis, which is thereby accessible
!! by all modules, which use this module.
module atl_space_basis

  use tem_aux_module,             only: tem_abort

  use atl_scheme_module,          only: atl_scheme_type,        &
    &                                   atl_modg_scheme_prp,    &
    &                                   atl_modg_2d_scheme_prp, &
    &                                   atl_modg_1d_scheme_prp
  use ply_modg_basis_module,      only: ply_init_modg_multilevelCoeffs, &
    &                                   ply_modg_refine_type,           &
    &                                   ply_modg_covolume_type,         &
    &                                   ply_init_modg_covolumeCoeffs


  implicit none
  private


  public :: atl_init_spacebasis

contains

  !> Initialize the space basis, this subroutine has to be called before
  !! using the module variable space_basis#basis.
  !!
  !! It reads the spatial part of the scheme table in the configuration
  !! script and fills the basis accordingly.
  subroutine atl_init_spacebasis( scheme_list, minlevel, maxlevel )
    !> The minumum refinement level of the mesh
    integer, intent(in) :: minlevel
    !> The maximum refinement level of the mesh
    integer, intent(in) :: maxlevel
    !> The list of schemes on the different level
    type(atl_scheme_type), intent(inout) :: scheme_list(minlevel:maxlevel)
    ! --------------------------------------------------------------------------
    integer :: degree
    type(ply_modg_covolume_type) :: covolumeBaseCoeff
    type(ply_modg_refine_type) :: refineBaseCoeff
    ! --------------------------------------------------------------------------

    ! Check if we use the same scheme globally
    if (any(scheme_list(:)%scheme.ne.scheme_list(minlevel)%scheme)) then
      call tem_abort( 'ERROR in atl_init_spacebasis: ' &
        & // ' only global schemes are supported'      )
    end if

    degree = 0 ! default value
    select case(scheme_list(minlevel)%scheme)
    case(atl_modg_scheme_prp)
      degree = maxval(scheme_list(:)%modg%maxPolyDegree)

    case(atl_modg_2d_scheme_prp)
      degree = maxval(scheme_list(:)%modg_2d%maxPolyDegree)

    case(atl_modg_1d_scheme_prp)
      degree = maxval(scheme_list(:)%modg_1d%maxPolyDegree)

    case default
      call tem_abort( 'ERROR in atl_init_spacebasis: unknown spatial scheme' )

    end select

    call ply_init_modg_covolumeCoeffs( nPoints  = degree * 2 + 1,   &
      &                                nFunc    = degree + 1,       &
      &                                integral = covolumeBaseCoeff )
    scheme_list(:)%modg_basis%covolumeBaseCoeff = covolumeBaseCoeff

    call ply_init_modg_multilevelCoeffs( nPoints  = degree * 2 + 1, &
      &                                  nFunc    = degree + 1,     &
      &                                  integral = refineBaseCoeff )
    scheme_list(:)%modg_basis%refineBaseCoeff = refineBaseCoeff

  end subroutine atl_init_spacebasis

end module atl_space_basis
