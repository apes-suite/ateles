! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt for University of Siegen.
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

program ply_fpt_header_test

  use env_module,               only: rk

  use ply_fpt_header_module,    only: ply_fpt_header_type, &
    &                                 operator(==),        &
    &                                 operator(/=),        &
    &                                 operator(<),         &
    &                                 operator(<=),        &
    &                                 operator(>),         &
    &                                 operator(>=)

  implicit none

  if ( isEqual()                &
    & .and. isUnequal()         &
    & .and. isSmaller()         &
    & .and. isSmallerOrEqual()  &
    & .and. isGreater()         &
    & .and. isGreaterOrEqual () ) then

    write(*,*) 'PASSED'

  end if

contains

  subroutine init_fptHeader(header)
    type(ply_fpt_header_type), intent(inout) :: header
    header%nodes_header%nodes_kind = 'my kind'
    header%nodes_header%lobattoPoints = .true.
    header%factor = 3_rk
    header%blocksize = 4711
    header%approx_terms = 815
    header%striplen = 42
    header%subblockingWidth = 23
    header%adapt_factor_pow2 = .true.
  end subroutine

  logical function isEqual()
    type(ply_fpt_header_type) :: norm
    type(ply_fpt_header_type) :: nodes_header
    type(ply_fpt_header_type) :: factor
    type(ply_fpt_header_type) :: blocksize
    type(ply_fpt_header_type) :: approx_terms
    type(ply_fpt_header_type) :: striplen
    type(ply_fpt_header_type) :: subblockingWidth
    type(ply_fpt_header_type) :: adapt_factor_pow2

    call init_fptHeader( norm )
    call init_fptHeader( nodes_header )
    call init_fptHeader( factor )
    call init_fptHeader( blocksize )
    call init_fptHeader( approx_terms )
    call init_fptHeader( striplen )
    call init_fptHeader( subblockingWidth )
    call init_fptHeader( adapt_factor_pow2 )

    nodes_header%nodes_header%lobattopoints = .false.
    factor%factor = -1_rk
    blocksize%blocksize = -1
    approx_terms%approx_terms = -1
    striplen%striplen = -1
    subblockingWidth%subblockingWidth = -1
    adapt_factor_pow2%adapt_factor_pow2 = .false.

    isEqual = norm == norm &
      & .and..not. norm == nodes_header &
      & .and..not. nodes_header == norm

    isEqual = isEqual .and. norm == norm &
      & .and..not. norm == factor &
      & .and..not. factor == norm

    isEqual = isEqual .and. norm == norm &
      & .and..not. norm == blocksize &
      & .and..not. blocksize == norm

    isEqual = isEqual .and. norm == norm &
      & .and..not. norm == approx_terms &
      & .and..not. approx_terms == norm

    isEqual = isEqual .and. norm == norm &
      & .and..not. norm == striplen &
      & .and..not. striplen == norm

    isEqual = isEqual .and. norm == norm &
      & .and..not. norm == subblockingWidth &
      & .and..not. subblockingWidth == norm

    isEqual = isEqual .and. norm == norm &
      & .and..not. norm == adapt_factor_pow2 &
      & .and..not. adapt_factor_pow2 == norm

    if( .not. isEqual ) write(*,*) 'isEqual failed'

  end function isEqual

  logical function isUnequal()
    type(ply_fpt_header_type) :: norm
    type(ply_fpt_header_type) :: nodes_header
    type(ply_fpt_header_type) :: factor
    type(ply_fpt_header_type) :: blocksize
    type(ply_fpt_header_type) :: approx_terms
    type(ply_fpt_header_type) :: striplen
    type(ply_fpt_header_type) :: subblockingWidth
    type(ply_fpt_header_type) :: adapt_factor_pow2

    call init_fptHeader( norm )
    call init_fptHeader( nodes_header )
    call init_fptHeader( factor )
    call init_fptHeader( blocksize )
    call init_fptHeader( approx_terms )
    call init_fptHeader( striplen )
    call init_fptHeader( subblockingWidth )
    call init_fptHeader( adapt_factor_pow2 )

    nodes_header%nodes_header%lobattopoints = .false.
    factor%factor = -1_rk
    blocksize%blocksize = -1
    approx_terms%approx_terms = -1
    striplen%striplen = -1
    subblockingWidth%subblockingWidth = -1
    adapt_factor_pow2%adapt_factor_pow2 = .false.

    isUnequal = .not. norm /= norm &
      & .and. norm /= nodes_header &
      & .and. nodes_header /= norm

    isUnequal = isUnequal .and..not. norm /= norm &
      & .and. norm /= factor &
      & .and. factor /= norm

    isUnequal = isUnequal .and..not. norm /= norm &
      & .and. norm /= blocksize &
      & .and. blocksize /= norm

    isUnequal = isUnequal .and..not. norm /= norm &
      & .and. norm /= approx_terms &
      & .and. approx_terms /= norm

    isUnequal = isUnequal .and..not. norm /= norm &
      & .and. norm /= striplen &
      & .and. striplen /= norm

    isUnequal = isUnequal .and..not. norm /= norm &
      & .and. norm /= subblockingWidth &
      & .and. subblockingWidth /= norm

    isUnequal = isUnequal .and..not. norm /= norm &
      & .and. norm /= adapt_factor_pow2 &
      & .and. adapt_factor_pow2 /= norm

    if( .not. isUnequal ) write(*,*) 'isUnequal failed'

  end function isUnequal

  logical function isSmaller()
    type(ply_fpt_header_type) :: norm
    type(ply_fpt_header_type) :: nodes_header
    type(ply_fpt_header_type) :: factor
    type(ply_fpt_header_type) :: blocksize
    type(ply_fpt_header_type) :: approx_terms
    type(ply_fpt_header_type) :: striplen
    type(ply_fpt_header_type) :: subblockingWidth
    type(ply_fpt_header_type) :: adapt_factor_pow2

    call init_fptHeader( norm )
    call init_fptHeader( nodes_header )
    call init_fptHeader( factor )
    call init_fptHeader( blocksize )
    call init_fptHeader( approx_terms )
    call init_fptHeader( striplen )
    call init_fptHeader( subblockingWidth )
    call init_fptHeader( adapt_factor_pow2 )

    nodes_header%nodes_header%lobattopoints = .false.
    factor%factor = -1_rk
    blocksize%blocksize = -1
    approx_terms%approx_terms = -1
    striplen%striplen = -1
    subblockingWidth%subblockingWidth = -1
    adapt_factor_pow2%adapt_factor_pow2 = .false.

    isSmaller = .not. norm < norm

    isSmaller = isSmaller .and..not. norm < nodes_header &
      & .and. nodes_header < norm

    isSmaller = isSmaller .and..not. norm < factor .and. factor < norm

    isSmaller = isSmaller .and..not. norm < blocksize .and. blocksize < norm

    isSmaller = isSmaller .and..not. norm < approx_terms &
      & .and. approx_terms < norm

    isSmaller = isSmaller .and..not. norm < striplen .and. striplen < norm

    isSmaller = isSmaller .and..not. norm < subblockingWidth &
      & .and. subblockingWidth < norm

    isSmaller = isSmaller .and..not. norm < adapt_factor_pow2 &
      & .and. adapt_factor_pow2 < norm

    if( .not. isSmaller ) write(*,*) 'isSmaller failed'

  end function isSmaller

  logical function isSmallerOrEqual()
    type(ply_fpt_header_type) :: norm
    type(ply_fpt_header_type) :: nodes_header
    type(ply_fpt_header_type) :: factor
    type(ply_fpt_header_type) :: blocksize
    type(ply_fpt_header_type) :: approx_terms
    type(ply_fpt_header_type) :: striplen
    type(ply_fpt_header_type) :: subblockingWidth
    type(ply_fpt_header_type) :: adapt_factor_pow2

    call init_fptHeader( norm )
    call init_fptHeader( nodes_header )
    call init_fptHeader( factor )
    call init_fptHeader( blocksize )
    call init_fptHeader( approx_terms )
    call init_fptHeader( striplen )
    call init_fptHeader( subblockingWidth )
    call init_fptHeader( adapt_factor_pow2 )

    nodes_header%nodes_header%lobattopoints = .false.
    factor%factor = -1_rk
    blocksize%blocksize = -1
    approx_terms%approx_terms = -1
    striplen%striplen = -1
    subblockingWidth%subblockingWidth = -1
    adapt_factor_pow2%adapt_factor_pow2 = .false.

    isSmallerOrEqual = norm <= norm

    isSmallerOrEqual = isSmallerOrEqual .and..not. norm <= nodes_header &
      & .and. nodes_header <= norm

    isSmallerOrEqual = isSmallerOrEqual .and..not. norm <= factor &
      & .and. factor <= norm

    isSmallerOrEqual = isSmallerOrEqual .and..not. norm <= blocksize &
      & .and. blocksize <= norm

    isSmallerOrEqual = isSmallerOrEqual .and..not. norm <= approx_terms &
      & .and. approx_terms <= norm

    isSmallerOrEqual = isSmallerOrEqual .and..not. norm <= striplen &
      & .and. striplen <= norm

    isSmallerOrEqual = isSmallerOrEqual .and..not. norm <= subblockingWidth &
      & .and. subblockingWidth <= norm

    isSmallerOrEqual = isSmallerOrEqual .and..not. norm <= adapt_factor_pow2 &
      & .and. adapt_factor_pow2 <= norm

    if( .not. isSmallerOrEqual ) write(*,*) 'isSmallerOrEqual failed'

  end function isSmallerOrEqual

  logical function isGreater()
    type(ply_fpt_header_type) :: norm
    type(ply_fpt_header_type) :: nodes_header
    type(ply_fpt_header_type) :: factor
    type(ply_fpt_header_type) :: blocksize
    type(ply_fpt_header_type) :: approx_terms
    type(ply_fpt_header_type) :: striplen
    type(ply_fpt_header_type) :: subblockingWidth
    type(ply_fpt_header_type) :: adapt_factor_pow2

    call init_fptHeader( norm )
    call init_fptHeader( nodes_header )
    call init_fptHeader( factor )
    call init_fptHeader( blocksize )
    call init_fptHeader( approx_terms )
    call init_fptHeader( striplen )
    call init_fptHeader( subblockingWidth )
    call init_fptHeader( adapt_factor_pow2 )

    nodes_header%nodes_header%lobattopoints = .false.
    factor%factor = -1_rk
    blocksize%blocksize = -1
    approx_terms%approx_terms = -1
    striplen%striplen = -1
    subblockingWidth%subblockingWidth = -1
    adapt_factor_pow2%adapt_factor_pow2 = .false.

    isGreater = .not. norm > norm

    isGreater = isGreater .and. norm > nodes_header &
      & .and..not. nodes_header > norm

    isGreater = isGreater .and. norm > factor .and..not. factor > norm

    isGreater = isGreater .and. norm > blocksize .and..not. blocksize > norm

    isGreater = isGreater .and. norm > approx_terms &
      & .and..not. approx_terms > norm

    isGreater = isGreater .and. norm > striplen .and..not. striplen > norm

    isGreater = isGreater .and. norm > subblockingWidth &
      & .and..not. subblockingWidth > norm

    isGreater = isGreater .and. norm > adapt_factor_pow2 &
      & .and..not. adapt_factor_pow2 > norm

    if( .not. isGreater ) write(*,*) 'isGreater failed'

  end function isGreater

  logical function isGreaterOrEqual()
    type(ply_fpt_header_type) :: norm
    type(ply_fpt_header_type) :: nodes_header
    type(ply_fpt_header_type) :: factor
    type(ply_fpt_header_type) :: blocksize
    type(ply_fpt_header_type) :: approx_terms
    type(ply_fpt_header_type) :: striplen
    type(ply_fpt_header_type) :: subblockingWidth
    type(ply_fpt_header_type) :: adapt_factor_pow2

    call init_fptHeader( norm )
    call init_fptHeader( nodes_header )
    call init_fptHeader( factor )
    call init_fptHeader( blocksize )
    call init_fptHeader( approx_terms )
    call init_fptHeader( striplen )
    call init_fptHeader( subblockingWidth )
    call init_fptHeader( adapt_factor_pow2 )

    nodes_header%nodes_header%lobattopoints = .false.
    factor%factor = -1_rk
    blocksize%blocksize = -1
    approx_terms%approx_terms = -1
    striplen%striplen = -1
    subblockingWidth%subblockingWidth = -1
    adapt_factor_pow2%adapt_factor_pow2 = .false.

    isGreaterOrEqual = norm >= norm

    isGreaterOrEqual = isGreaterOrEqual .and. norm >= nodes_header &
      & .and..not. nodes_header >= norm

    isGreaterOrEqual = isGreaterOrEqual .and. norm >= factor &
      & .and..not. factor >= norm

    isGreaterOrEqual = isGreaterOrEqual .and. norm >= blocksize &
      & .and..not. blocksize >= norm

    isGreaterOrEqual = isGreaterOrEqual .and. norm >= approx_terms &
      & .and..not. approx_terms >= norm

    isGreaterOrEqual = isGreaterOrEqual .and. norm >= striplen &
      & .and..not. striplen >= norm

    isGreaterOrEqual = isGreaterOrEqual .and. norm >= subblockingWidth &
      & .and..not. subblockingWidth >= norm

    isGreaterOrEqual = isGreaterOrEqual .and. norm >= adapt_factor_pow2 &
      & .and..not. adapt_factor_pow2 >= norm

    if( .not. isGreaterOrEqual ) write(*,*) 'isGreaterOrEqual failed'

  end function isGreaterOrEqual

end program ply_fpt_header_test

