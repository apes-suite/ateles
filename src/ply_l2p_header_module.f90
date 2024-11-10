! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2013-2014,2017,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013,2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Peter Vitt, Verena Krupp,
! Nikhil Anand, Daniel Petró and Neda Ebrahimi Pour for University of Siegen.
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

!> Parameters for the plain L2 projection method to transform between Legendre
!! modes and nodal representation.
!!
!! This method utilizes the L2 projection from Legendre to Lagrange polynomials
!! or the other way around. A numerical Gauss-Legendre Quadrature is used to
!! compute the integral over the product of both functions.
!! The Lagrange polynomials can be defined on any nodeset, see also
!! [[ply_nodeset_module]].
!!
!! Available options for the nodes to project onto are:
!!
!! * `'gauss-legendre'` these are the Gauss-Legendre integration points that
!!   are also used for the numerical integration (this is the default).
!! * `'chebyshev'` these are the nodes from the Chebyshev integration.
!!
!! The set of nodes to use is configured by the `nodes_kind` option, and if
!! `nodes_kind = 'chebyshev'` it is also possible to make use of Lobatto points
!! to include the interval boundaries in the nodal representation.
!! This is achieved by setting `lobattoPoints = true`, by default this is false.
!!
!! The configuration table for a projection with L2P may, for example, look as
!! follows:
!!
!!```lua
!!  projection = {
!!    kind = 'l2p',
!!    factor = 1.5,
!!    nodes_kind = 'chebyshev',
!!    lobattoPoints = true
!!  }
!!```
!!
!! The example illustrates the three possible settings for the L2P
!! transformation method:
!!
!! * `factor` - Oversampling factor to avoid aliasing.
!! * `nodes_kind` - Selection of set of nodes to use in the nodal
!!   representation.
!! * `lobattoPoints` - Whether to include interval bounds, only
!!   available for Chebyshev nodes.
!!
module ply_l2p_header_module

  use env_module,              only: rk, labelLen
  use aotus_module,            only: flu_State, aot_get_val
  use aot_out_module,          only: aot_out_type, aot_out_val

  use tem_tools_module,        only: upper_to_lower
  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit
  use tem_float_module

  use ply_nodes_header_module

  implicit none

  private

  !> l2p projection header type, consisting of the node header which give
  !! information about the type and number of points for the projection
  type ply_l2p_header_type
    type(ply_nodes_header_type) :: nodes_header
    real(kind=rk) :: factor
  end type ply_l2p_header_type

  interface assignment(=)
    module procedure Copy_l2p_header
  end interface

  interface operator(==)
    module procedure isEqual
  end interface

  interface operator(/=)
    module procedure isUnequal
  end interface

  interface operator(<)
    module procedure isSmaller
  end interface

  interface operator(<=)
    module procedure isSmallerOrEqual
  end interface

  interface operator(>)
    module procedure isGreater
  end interface

  interface operator(>=)
    module procedure isGreaterOrEqual
  end interface

  public :: operator(==), operator(/=), operator(<), operator(<=)
  public :: operator(>), operator(>=)
  public :: assignment(=)
  public :: ply_l2p_header_type
  public :: ply_l2p_header_load,  ply_l2p_header_display
  public :: ply_l2p_header_define
  public :: ply_l2p_header_out


contains


  ! ************************************************************************ !
  pure subroutine Copy_l2p_header( left, right )
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_l2p_header_type), intent(out) :: left
    !> fpt to copy from
    type(ply_l2p_header_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%factor = right%factor
    left%nodes_header = right%nodes_header

  end subroutine Copy_l2p_header
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load settings to describe a projection method from a Lua table.
  subroutine ply_l2p_header_load( me, conf, thandle )
    ! -------------------------------------------------------------------- !
    type(ply_l2p_header_type), intent(out) :: me
    type(flu_State), intent(inout) :: conf
    integer, intent(in) :: thandle
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: nodes_kind
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! for l2p gauss-legendre points are used
    me%nodes_header%nodes_kind = 'gauss-legendre'

    ! fill up l2p header
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'factor',  &
      &               val     = me%factor, &
      &               default = 1.0_rk,    &
      &               ErrCode = iError     )
    call aot_get_val( L       = conf,                          &
      &               thandle = thandle,                       &
      &               key     = 'lobattoPoints',               &
      &               val     = me%nodes_header%lobattoPoints, &
      &               ErrCode = iError,                        &
      &               default = .false.                        )
    call aot_get_val( L       = conf,            &
      &               thandle = thandle,         &
      &               key     = 'nodes_kind',    &
      &               val     = nodes_kind,      &
      &               ErrCode = iError,          &
      &               default = 'gauss-legendre' )

    me%nodes_header%nodes_kind = upper_to_lower(nodes_kind)
    select case (trim(me%nodes_header%nodes_kind))
    case ('gauss-legendre')

      if ( me%nodes_header%lobattoPoints ) then
        write(logUnit(1),*) 'ERROR in loading projection: Legendre nodes' &
          &                 // ' do not support Lobatto points'
        write(logUnit(1),*) 'But you configured to use Lobatto Points!'
        write(logUnit(1),*) 'Probably you want to use Chebyshev nodes instead!'
        call tem_abort()
      end if

    case ('chebyshev')
      ! Nothing to do

    case default
      write(logUnit(1),*) 'ERROR in loading projection: nodes_kind'
      write(logUnit(1),*) trim(me%nodes_header%nodes_kind)
      write(logUnit(1),*) 'The nodes_kind needs to be one of:'
      write(logUnit(1),*) ' * "gauss-legendre"'
      write(logUnit(1),*) ' * "chebyshev"'
      write(logUnit(1),*) ''
      write(logUnit(1),*) 'Aborting...'
      call tem_abort()
    end select

    if (me%factor <= 0.0_rk) then
      write(logUnit(1),*) 'ERROR in loading projection: factor for ' &
        &                 // 'projection has to be larger than 0!'
      write(logUnit(1),*) 'But it is set to ', me%factor
      call tem_abort()
    end if

  end subroutine ply_l2p_header_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_l2p_header_define( me, factor, nodes_kind, lobattoPoints )
    ! -------------------------------------------------------------------- !
    !> L2P header to define.
    type(ply_l2p_header_type), intent(out) :: me

    !> Oversampling factor to use in the projection, defaults to 1.
    integer, optional, intent(in)     :: factor

    !> Set of nodes to use in the nodal representation.
    !!
    !! May be 'gauss-legendre' or 'chebyshev', defaults to 'gauss-legendre'
    character(len=*), optional :: nodes_kind

    !> Wether to use Lobatto points (include interval bounds) when using the
    !! chebyshev nodes, defaults to .false..
    logical, optional, intent(in)     :: lobattoPoints
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ! Set defaults:
    me%nodes_header%nodes_kind = 'gauss-legendre'
    me%nodes_header%lobattoPoints = .false.
    me%factor = 1.0_rk
    if (present(nodes_kind)) me%nodes_header%nodes_kind &
      &                      = upper_to_lower(nodes_kind)
    if (present(lobattoPoints)) me%nodes_header%lobattoPoints = lobattoPoints
    if (present(factor)) me%factor = factor


    select case (trim(me%nodes_header%nodes_kind))
    case ('gauss-legendre')

      if ( me%nodes_header%lobattoPoints ) then
        write(logUnit(1),*) 'ERROR in defining projection: Legendre nodes' &
          &                 // ' do not support Lobatto points'
        write(logUnit(1),*) 'But you configured to use Lobatto Points!'
        write(logUnit(1),*) 'Probably you want to use Chebyshev nodes instead!'
        call tem_abort()
      end if

    case ('chebyshev')
      ! Nothing to do

    case default
      write(logUnit(1),*) 'ERROR in defining projection: nodes_kind'
      write(logUnit(1),*) trim(me%nodes_header%nodes_kind)
      write(logUnit(1),*) 'The nodes_kind needs to be one of:'
      write(logUnit(1),*) ' * "gauss-legendre"'
      write(logUnit(1),*) ' * "chebyshev"'
      write(logUnit(1),*) ''
      write(logUnit(1),*) 'Aborting...'
      call tem_abort()
    end select

    if (me%factor <= 0.0_rk) then
      write(logUnit(1),*) 'ERROR in defining projection: factor for ' &
        &                 // 'projection has to be larger than 0!'
      write(logUnit(1),*) 'But it is set to ', me%factor
      call tem_abort()
    end if

  end subroutine ply_l2p_header_define
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write L2P settings into a Lua table.
  subroutine ply_l2p_header_out( me, conf )
    ! -------------------------------------------------------------------- !
    type(ply_l2p_header_type), intent(in) :: me
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !

    call aot_out_val( put_conf = conf,     &
      &               vname    = 'factor', &
      &               val      = me%factor )

    call aot_out_val( put_conf = conf,                            &
      &               vname    = 'nodes_kind',                    &
      &               val      = trim(me%nodes_header%nodes_kind) )

    call aot_out_val( put_conf = conf,                         &
      &               vname    = 'lobattoPoints',              &
      &               val      = me%nodes_header%lobattoPoints )

  end subroutine ply_l2p_header_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_l2p_header_display( me )
    ! -------------------------------------------------------------------- !
    type(ply_l2p_header_type), intent(in) :: me
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) ' * Kind of projection method = l2p'
    write(logUnit(1),*) ' * Factor to use in projection = ', me%factor
    write(logUnit(1),*) ' * using this set of nodes: ' &
      &                 // trim(me%nodes_header%nodes_kind)
    write(logUnit(1),*) ' * using LobattoPoints =', &
      &                 me%nodes_header%lobattoPoints

    if (me%factor < 2.0_rk) then
      write(logUnit(1),*) ''
      write(logUnit(1),*)                                         &
        & '+-----------------------------------------------------+'
      write(logUnit(1),*)                                         &
        & '| WARNING, the oversampling factor is smaller than 2! |'
      write(logUnit(1),*)                                         &
        & '|        this might lead to bad projections!          |'
      write(logUnit(1),*)                                         &
        & '+-----------------------------------------------------+'
    end if

  end subroutine ply_l2p_header_display
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for equality of two projections.
  !!
  !! Two l2p header are considered to be equal, if their node_header,
  !! and the factor are equal.
  pure function isEqual( left, right ) result( equality )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is equal??
    logical :: equality
    ! -------------------------------------------------------------------- !

    equality = ( left%nodes_header == right%nodes_header ) &
      & .and. ( left%factor .feq. right%factor )

  end function isEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for unequality of two projections.
  !!
  !! Two l2p header are considered to be unequal, if their node_header,
  !! or the factor are not equal.
  pure function isUnequal( left, right ) result( unequality )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    ! -------------------------------------------------------------------- !

    unequality = ( left%nodes_header /= right%nodes_header ) &
      & .or. ( left%factor .fne. right%factor )

  end function isUnequal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a < comparison of two projections.
  !!
  !! Sorting of l2p header is given by node_header and by the factor.
  pure function isSmaller( left, right ) result( small )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !

    small = .false.
    if (left%nodes_header < right%nodes_header) then
      small = .true.
    else
      if (left%nodes_header == right%nodes_header) then
        small = (left%factor < right%factor)
      end if
    end if

  end function isSmaller
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a <= comparison of two projections.
  !!
  !! Sorting of l2p header is given by node_header, l2p_blocksize and
  !! last by factor.
  pure function isSmallerOrEqual( left, right ) result( small )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !

    small = .false.
    if (left%nodes_header < right%nodes_header) then
      small = .true.
    else
      if (left%nodes_header == right%nodes_header) then
        small = (left%factor <= right%factor)
      end if
    end if

  end function isSmallerOrEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a > comparison of two projections.
  !!
  !! Sorting of l2p header is given by node_header, l2p_blocksize and
  !! last by factor.
  pure function isGreater( left, right ) result( great )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !

    great = .false.
    if (left%nodes_header > right%nodes_header) then
      great = .true.
    else
      if (left%nodes_header == right%nodes_header) then
        great = (left%factor > right%factor)
      end if
    end if

  end function isGreater
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a >= comparison of two projections.
  !!
  !! Sorting of l2p header is given by node_header, l2p_blocksize and
  !! last by factor.
  pure function isGreaterOrEqual( left, right ) result( great )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !

    great = .false.
    if (left%nodes_header > right%nodes_header) then
      great = .true.
    else
      if (left%nodes_header == right%nodes_header) then
        great = (left%factor >= right%factor)
      end if
    end if

  end function isGreaterOrEqual
  ! ************************************************************************ !

end module ply_l2p_header_module
