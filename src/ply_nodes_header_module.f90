! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2013,2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Verena Krupp, Peter Vitt, Nikhil Anand
! and Harald Klimach for University of Siegen.
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
module ply_nodes_header_module
  use env_module,        only: labelLen

  implicit none

  private

  type ply_nodes_header_type
    !> which kind of nodes are used. For l2p projection the nodes are
    !! legendre-gauss and chebyshev nodes for using fpt
    character(len=labelLen) :: nodes_kind
    !> Logical to indicate whether Chebyshev-Lobatto points or simple
    !! Chebyshev points are used
    logical :: lobattoPoints = .false.
  end type ply_nodes_header_type

  interface assignment(=)
    module procedure Copy_nodes_header
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
  public :: ply_nodes_header_type


contains

  ! ************************************************************************ !
  pure subroutine Copy_nodes_header( left, right )
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    !type(ply_legFpt_2D_type), intent(out) :: left
    type(ply_nodes_header_type), intent(out) :: left
    !> fpt to copy from
    !type(ply_legFpt_2D_type), intent(in) :: right
    type(ply_nodes_header_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%nodes_kind = right%nodes_kind
    left%lobattoPoints = right%lobattoPoints

  end subroutine Copy_nodes_header
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for equality of two nodes descriptions.
  pure function isEqual( left, right ) result(equality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_nodes_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_nodes_header_type), intent(in) :: right
    !> is equal??
    logical :: equality
    ! -------------------------------------------------------------------- !

    equality = (left%nodes_kind == right%nodes_kind)              &
      &        .and. (left%lobattopoints .eqv. right%lobattopoints)

  end function isEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for unequality of two nodes descriptions.
  pure function isUnequal( left, right ) result(unequality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_nodes_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_nodes_header_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    ! -------------------------------------------------------------------- !

    unequality = ( left%nodes_kind /= right%nodes_kind ) &
      &        .or. (left%lobattopoints .neqv. right%lobattopoints)

  end function isUnequal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a < comparison of two nodes descriptions.
  pure function isSmaller( left, right ) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_nodes_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_nodes_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !
    small = .false.

    if (left%lobattopoints) then
      left_log = 1
    else
      left_log = 0
    end if
    if (right%lobattopoints) then
      right_log = 1
    else
      right_log = 0
    end if

    if (left%nodes_kind < right%nodes_kind) then
      small =.true.
    else
      if (left%nodes_kind == right%nodes_kind) then
        small = (left_log < right_log)
      end if
    end if

  end function isSmaller
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a <= comparison of two nodes descriptions.
  pure function isSmallerOrEqual( left, right ) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_nodes_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_nodes_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !
    small = .false.

    if (left%lobattopoints) then
      left_log = 1
    else
      left_log = 0
    end if
    if (right%lobattopoints) then
      right_log = 1
    else
      right_log = 0
    end if

    if  (left%nodes_kind == right%nodes_kind) then
      small =.true.
    else
      if (left%nodes_kind == right%nodes_kind) then
        small = (left_log <= right_log)
      end if
    end if

  end function isSmallerOrEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a > comparison of nodes descriptions.
  pure function isGreater( left, right ) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_nodes_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_nodes_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !
    great = .false.

    if (left%lobattopoints) then
      left_log = 1
    else
      left_log = 0
    end if
    if (right%lobattopoints) then
      right_log = 1
    else
      right_log = 0
    end if

    if  (left%nodes_kind > right%nodes_kind) then
      great =.true.
    else
      if (left%nodes_kind == right%nodes_kind) then
        great = (left_log > right_log)
      end if
    end if

  end function isGreater
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a >= comparison of two nodes descriptions.
  pure function isGreaterOrEqual( left, right ) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_nodes_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_nodes_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !
    great = .false.

    if (left%lobattopoints) then
      left_log = 1
    else
      left_log = 0
    end if
    if (right%lobattopoints) then
      right_log = 1
    else
      right_log = 0
    end if

    if  (left%nodes_kind == right%nodes_kind) then
      great =.true.
    else
      if (left%nodes_kind == right%nodes_kind) then
        great = (left_log >= right_log)
      end if
    end if

  end function isGreaterOrEqual
  ! ************************************************************************ !

end module ply_nodes_header_module
