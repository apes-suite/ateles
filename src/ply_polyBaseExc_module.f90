! Copyright (c) 2012-2014,2016-2018,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Langhammer Kay <kay.langhammer@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop and Harald Klimach
! for German Research School for Simulation Sciences GmbH.
!
! Parts of this file were written by Verena Krupp, Harald Klimach, Peter Vitt,
! Nikhil Anand, Kay Langhammer, Daniel Petró and Neda Ebrahimi Pour for
! University of Siegen.
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

! Copyright (c) 2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Langhammer Kay <kay.langhammer@student.uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt, Harald Klimach and
! Kay Langhammer for University of Siegen.
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









! Make sure unroll is defined to an integer value.
! Usually, this should be done on the command line with -Dunroll=<n>

!> Module to change bases functions of a modal representation.
!!
!! \author{Jens Zudrop}
!!
!! This module provides routines for a fast basis exchange between Legendre
!! and Chebyshev polynomials. The alogorihtm is fast as it completes this task
!! in O(n log(n)) operations (where n is the number of modal coefficients).
!! The alogrithm is based on the following publication:
!! Alpert, B., & Rokhlin, V. (1991).
!! A fast algorithm for the evaluation of Legendre expansions.
!! SIAM Journal on Scientific and Statistical Computing.
!! There also an alternative implementation with O(n) operations is described.
!!
!! Some recommendations to achieve a fast transformation.
!! The number of minimal blocks (n/s) should be a power of two plus 1:
!! (n/s) = 2^k + 1
!! This yields the minimal number of blocks, that need to be computed.
!! If this is not possible, it is probably good to at lease have an odd number
!! of blocks (n/s = 2*k + 1), to reduce the number of smallest blocks.
!! The remainder mod(n, s) + s, should be even, as otherwise there is an
!! additional diagonal that needs to be computed.
!! Similarily also s itself should probably be even.
module ply_polyBaseExc_module
  use, intrinsic :: iso_c_binding
  use env_module,            only: rk
  use tem_float_module
  use tem_param_module,      only: pi
  use tem_gamma_module
  use tem_logging_module,    only: logUnit
  use ply_fpt_header_module, only: ply_fpt_default_subblockingWidth
  use fftw_wrap

  implicit none
  private

  type ply_sub_vec
    real(kind=rk), allocatable :: dat(:)
  end type ply_sub_vec

  !> Expansion coefficients for a certain submatrix.
  type ply_matrixExpCoeff_type
    !> Polynomials expansion coefficients.
    real(kind=rk), allocatable :: coeff(:)
  end type ply_matrixExpCoeff_type

  !> Information for a set of local rows in the current block
  type ply_coldata_type
    !> The Chebyshev expansion coefficients for a set of block local
    !! rows.
    type(ply_matrixExpCoeff_type), allocatable :: rowDat(:)
  end type ply_coldata_type

  !> Sparse data for information of a column in a sub matrix
  type ply_rowdata_type
    !> Column data. Dimension is the different number of column rows in the
    !! sub matrix.
    type(ply_coldata_type), allocatable :: subCol(:)
  end type ply_rowdata_type

  !> Sparse data for information of a row in a sub matrix.
  type ply_submatrix_type
    !> Row data. Dimension is the different number of block rows in the sub matrix.
    type(ply_rowdata_type), allocatable :: subRow(:)
  end type ply_submatrix_type

  !> Expansion coefficients for a certain submatrix.
  type ply_matrixExpCoeffOddEven_type
    !> Polynomials expansion coefficients.
    real(kind=rk), allocatable :: coeff(:,:)
  end type ply_matrixExpCoeffOddEven_type

  !> Sparse data for a subvector
  type ply_subvector_type
    !> Expansion coefficients for a column
    type(ply_matrixExpCoeffOddEven_type), allocatable :: col(:)
  end type ply_subvector_type

  type ply_trafo_params_type
    !> Lagrange polynomials evaluated at the Chebyshev points on [0,+1].
    type(ply_sub_vec), allocatable :: u(:,:)

    !> The array to store the diagonals of the matrix in.
    !!
    !! The first index are the rows, the second index are the columns.
    real(kind=rk), allocatable :: diag(:,:)

    !> The array to store the adapters between diagonal and
    !! blocks in.
    !!
    !! The first index are the rows, the second the diagonals,
    !! and the third the adapter.
    real(kind=rk), allocatable :: adapter(:,:,:)

    !> Data of all sub matrices (separated from the diagonal).
    !! Size is the number of different sub matrix sizes, i.e. h.
    type(ply_submatrix_type), allocatable :: sub(:)

    !> Number of blocks in one direction of the matrix.
    integer :: nBlocks

    !> Length of stripes to use in the matrix operation.
    integer :: striplen

    !> Remaining columns close to the diagonal after subdividing
    !! the matrix into blocks
    integer :: remainder

    !> Number of full diagonals that need to be considered close
    !! the diagonal of the matrix.
    integer :: nDiagonals

    !> Number of diagonals in triangle blocks, that remain between
    !! blocks and full diagonals.
    integer :: nBlockDiagonals

    !> The number of modal coefficients to convert
    integer :: n

    !> The number of Cheb coefficients to approximate M
    integer :: k

    !> The size of the smallest subblock of M
    integer :: s

    !> The number of subblocks (per direction) in M
    integer :: h

    !> The width of the subblocks used during the unrolled base exchange to
    !! ensure a better cache usage.
    integer :: subblockingWidth

    !> The transformation type
    integer :: trafo

    !> Conversion data structure used for fpt.
    type(ply_subvector_type), allocatable :: b(:)
  end type ply_trafo_params_type

  interface assignment(=)
    module procedure Copy_trafo_params
  end interface

  integer, parameter :: ply_legToCheb_param = 1
  integer, parameter :: ply_chebToLeg_param = 2

  public :: ply_fpt_init
  public :: ply_fpt_exec_striped
  public :: ply_fpt_exec
  public :: ply_fpt_single
  public :: ply_trafo_params_type
  public :: ply_legToCheb_param, ply_chebToLeg_param
  public :: ply_lambda
  public :: assignment(=)


contains


  ! ************************************************************************ !
  subroutine Copy_trafo_params( left, right )
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_trafo_params_type), intent(out) :: left
    !> fpt to copy from
    type(ply_trafo_params_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%trafo = right%trafo
    left%h = right%h
    left%s = right%s
    left%k = right%k
    left%n = right%n
    left%subblockingWidth = right%subblockingWidth

    left%nBlocks = right%nBlocks

    left%striplen = right%striplen

    if(allocated(left%sub))then
      deallocate(left%sub)
    end if
    allocate(left%sub(1:right%h))
    left%sub = right%sub
    if(allocated(left%diag))then
      deallocate(left%diag)
    end if
    if(allocated(left%adapter))then
      deallocate(left%adapter)
    end if

    allocate(left%diag(right%n, (right%s+mod(right%s,2)/2)))
    left%diag = right%diag
    allocate( left%adapter(right%s, (right%s+mod(right%s,2)/2), right%nBlocks) )
    left%adapter = right%adapter
    if(allocated(left%u))then
      deallocate(left%u)
    end if
    allocate(left%u(1:right%h,0:right%k-1))
    left%u = right%u

  end subroutine Copy_trafo_params
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_fpt_init( n, params, trafo, blocksize, approx_terms, &
    &                      striplen, subblockingWidth                 )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: n
    type(ply_trafo_params_type), intent(inout) :: params
    integer, intent(in) :: trafo

    !> Smallest block that is approximated by approx_terms coefficients.
    !!
    !! Please note, that this has to be larger than 2*approx_terms to result
    !! in a reduced number of operations. Default is 64.
    integer, optional, intent(in) :: blocksize

    !> Number of approximation terms used to compute off-diagonal products.
    !!
    !! Defaults to 18, which is the suggested accuracy for double precision.
    integer, optional, intent(in) :: approx_terms

    !> Length to use in vectorization, this is the number of independent
    !! matrix multiplications that are to be done simultaneously.
    integer, intent(in) :: striplen

    !> The width of the subblocks used during the unrolled base exchange to
    !! ensure a better cache usage.
    integer, optional, intent(in) :: subblockingWidth
    ! -------------------------------------------------------------------- !
    integer :: r, l, k, s, h, i, j, m, diagonals, blockdiagonals
    integer :: remainder
    integer :: diag_off, block_off
    integer :: nRows
    integer :: ub_row, row_rem
    integer :: rowsize
    real(kind=rk), allocatable :: den(:), t(:)
    real(kind=rk) :: x
    type(ply_submatrix_type), allocatable :: sub(:)
    type(ply_sub_vec), allocatable :: u(:,:)
    ! -------------------------------------------------------------------- !

    write(logUnit(3),*) 'Initializing FPT'

    params%trafo = trafo
    params%striplen = striplen

    if(present(subblockingWidth)) then
      params%subblockingWidth = subblockingWidth
    else
      params%subblockingWidth = ply_fpt_default_subblockingWidth
    endif

    write(logUnit(3),*) 'subblockingWidth = ', params%subblockingWidth

    ! The number of polynomials to approximate the entries in M.
    ! As double precision is required, the paper suggests k = 18.
    if (present(approx_terms)) then
      k = approx_terms
    else
      k = 18
    end if
    write(logUnit(3),*) &
      & 'Number of Cheb. coeffs for matrix function approximation: ', k

    ! The minimum size of a submatrix in M, the paper suggests:
    ! approx. 4*k and power of 2
    if (present(blocksize)) then
      s = blocksize
    else
      s = 64
    end if
    write(logUnit(3),*) 'Smallest subblock size for evaluation of M: ',s

    ! Check if the number of coefficients is smaller than the smallest subblock.
    ! if this is true, we set the smallest subblock to n, to ensure that the
    ! direct multiplication with the entries near the diagonal is working
    ! correctly.
    if (n<s) then
      s = n
      write(logUnit(3),*) &
        & 'Corrected smallest subblock size for evaluation of M: ',s
    end if

    params%nBlocks = n/s

    ! Logarithm of the maximal block size in terms of minimal block size s.
    ! The formula in the Alpert&Rohklin paper is only valid if n is a power of
    ! 2.
    if (params%nBlocks > 1) then
      ! The first block is always a diagonal one, and not available for the
      ! approximation blocks.
      ! Of the remaining ones only half of them can be used in the largest
      ! block, as that block otherwise would not be detached from the diagonal.
      h = int( log( real( params%nBlocks - 1, rk ) ) / log(2.0_rk) ) - 1
    else
      h = -1
    end if

    allocate(t(0:k-1))
    allocate(den(0:k-1))

    ! All diagonals close to the main diagonal are given by the remainder
    ! of the division of the overall matrix size by the block size.
    !
    ! The remainder are the first diagonals close to the main diagonal and have
    ! at least the length of one block.
    remainder = n - s * (params%nBlocks-1)
    params%remainder = remainder

    ! The non-zero diagonals within the remainder (only every second diagonal).
    ! Obviously it is beneficial to have an even remainder, as otherwise there
    ! is an additional diagonal to take into account.
    diagonals = (remainder + mod(remainder,2)) / 2
    params%nDiagonals = diagonals

    !HK: All approximation related initializations are actually only required
    !HK: if params%nBlocks > 2, maybe we should check this here. However, the
    !HK: initialization should also not hurt too much in this case.

    ! Step 1
    ! Comment
    ! [Construct Chebyshev nodes to, tl,..., tk- on the
    ! interval [0, 1]
    do r = 0,k-1
      t(r) = 0.5_rk * ( 1.0_rk - cos( ( r + 0.5_rk ) * pi / k ) )
    end do

    ! Step 2
    ! Comment 1
    ! Evaluate the denominators in the expressions for the Chebyshev
    ! interpolation coefficients u_0, u_1, ... , u_{k-1}
    !NEC$ novector
    do r =0,k-1
      den(r) = 1.0_rk
      do l = 0,k-1
        if(l.ne.r) then
          den(r) = den(r) * ( t(r) - t(l) )
        end if
      end do
    enddo
    ! Comment 2
    ! Evaluate the Chebyshev interpolation coefficients uO, u1, ..., Uk-1
    ! at the uniformly spaced nodes 0,1/s,2/s,...,(s-1)/s
    allocate(u(0:h,0:k-1))
    !NEC$ novector
    do l = 0,h
      rowsize = s * 2**l
      !NEC$ novector
      do r = 0,k-1
        allocate(u(l,r)%dat(0:rowsize-1))
      end do
      !NEC$ novector
      do j = 0, rowsize-1
        x = 1.0_rk
        !NEC$ novector
        do m=0,k-1
          x = x * ( ( real(j,rk)/real(rowsize,rk)) - t(m) )
        end do
        do r = 0,k-1
          u(l,r)%dat(j) = (x/((real(j,rk)/real(rowsize,rk))-t(r))) / den(r)
        end do
      end do
    end do

    ! Step 3
    allocate(sub(0:h))
    if (trafo.eq.ply_legToCheb_param) then
      !NEC$ novector
      do l= 0,h
        rowsize = s * 2**l
        ! The Alpert & Rohklin paper has the wrong formula for the number of
        ! rows, needed on the different levels. It only works properly for
        ! even number of blocks in the matrix.
        ! To be correct, we have to exclude the first minimal block, as this
        ! is always the diagonal and then one block of the current level, as
        ! that would overlap the diagonal, thus we end up with:
        nRows = (params%nBlocks - 1) / (2**l) - 1
        ub_row = 3 - mod(nRows,2)
        row_rem = mod(n-remainder, rowsize) + remainder
        allocate(sub(l)%subRow(0:nRows-1))
        !NEC$ novector
        do i = 0, nRows - 1
          allocate(sub(l)%subRow(i)%subCol(i+2:i+ub_row-mod(i,2)))
          !NEC$ novector
          do j = i+2, i+ub_row-mod(i,2)
            allocate(sub(l)%subRow(i)%subCol(j)%rowDat(0:rowsize-1))
            !NEC$ novector
            do m = 0, rowsize-1
              allocate(sub(l)%subRow(i)%subCol(j)%rowDat(m)%coeff(0:k-1))
              do r = 0, k-1
                sub(l)%subRow(i)%subCol(j)%rowDat(m)%coeff(r)                 &
                  & = ply_m( m + real(i,rk) * real(rowsize,rk),               &
                  &          row_rem + (real(j-1,rk)+t(r)) * real(rowsize,rk) )
              end do ! r
            end do ! m
          end do ! j
        end do ! i
      end do ! l
    else
      !NEC$ novector
      do l= 0,h
        rowsize = s * 2**l
        ! The Alpert & Rokhlin paper has the wrong formula for the number of
        ! rows, needed on the different levels. It only works properly for
        ! even number of blocks in the matrix.
        ! To be correct, we have to exclude the first minimal block, as this
        ! is always the diagonal and then one block of the current level, as
        ! that would overlap the diagonal, thus we end up with:
        nRows = (params%nBlocks - 1) / (2**l) - 1
        ub_row = 3 - mod(nRows,2)
        row_rem = mod(n-remainder, rowsize) + remainder
        allocate(sub(l)%subRow(0:nRows-1))
        !NEC$ novector
        do i = 0, nRows - 1
          allocate(sub(l)%subRow(i)%subCol(i+2:i+ub_row-mod(i,2)))
          !NEC$ novector
          do j = i+2, i+ub_row-mod(i,2)
            allocate(sub(l)%subRow(i)%subCol(j)%rowDat(0:rowsize-1))
            !NEC$ novector
            do m = 0, rowsize-1
              allocate(sub(l)%subRow(i)%subCol(j)%rowDat(m)%coeff(0:k-1))
              do r = 0, k-1
                sub(l)%subRow(i)%subCol(j)%rowDat(m)%coeff(r)                &
                  & = ply_l( m + real(i,rk) * real(rowsize,rk),              &
                  &          row_rem + (real(j-1,rk)+t(r))* real(rowsize,rk) )
              end do ! r
            end do ! m
          end do ! j
        end do ! i
      end do ! l
    end if


    ! Step 4
    ! Allocate the array to store the virtual matrix diagonals without sparse
    ! diagoanls. The length of the array is the number of elements per block
    ! divided by 2, because every second diagonal is 0. Only the diagonals of
    ! the first block are calculated using the striped matrix vector
    ! multiplication so we only need space for this one block's diagonals.

    allocate( params%diag(n, diagonals) )

    ! We need nBlocks-1 adapters, to fill the triangles between the main
    ! diagonals and the stair cases of the approximated blocks.
    ! Each of them follows the same storage scheme as the main diagonals,
    ! however we only have a overall length of one block for each of them.
    ! Example for nBlocks = 3 (n = 3*s):
    !
    ! | block | block || block |
    ! --------------------------
    ! \   d   \  ad   ||       | -
    !  \   i   \  ap  || appr. |
    !   \   a   \  te ||       | b
    !    \   g   \  r ||       | l
    !     \   o   \   || block | o
    !      \   n   \  ||       | c
    !       \   a   \ ||       | k
    !        \   l   \||       |
    !         \   s   \|_______| -
    !          \       \  ad   | -
    !           \       \  ap  |
    !            \       \  te | b
    !             \       \  r | l
    !              \       \   | o
    !               \       \  | c
    !                \       \ | k
    !                 \       \|
    !                  \       | -
    !                   \      | -
    !                    \     |
    !                     \    | b
    !                      \   | l
    !                       \  | o
    !                        \ | c
    !                         \| k
    !                          \ _
    !
    ! If the remainder is not even, the adapters start with a 0 diagonal,
    ! and therefore have one diagonal less to store.
    ! odd remainder (5), odd block (3), 11x11 Matrix:
    ! \0\0\0x0MMM
    !  \0\0\0xMMM
    !   \0\0\0MMM
    !    \0\0\0x0
    !     \0\0\0x
    !      \0\0\0
    !       \0\0\
    !        \0\0
    !         \0\
    !          \0
    !           \
    ! Number of diagonals in each block:
    blockdiagonals = (s+remainder + mod(s+remainder,2)) / 2 - diagonals
    params%nBlockDiagonals = blockdiagonals

    allocate( params%adapter(s, blockdiagonals, params%nBlocks-1) )

    !HK: Not really required, as the 0 entries are supposedly never used
    !    later on.
    params%diag = 0.0_rk

    ! Copy the value of the diagonal into the rectangular representation:
    ! --------------------------------
    !  \0,0\0,1\0,2\...                        |0,0|0,2|...
    !      \1,1\1,2\1,3\...                    |1,1|1,3|...
    !          \2,2\2,3\2,4\...           -->  |2,2|2,4|...
    !              \3,3\3,4\3,5\..             |3,3|3,5|..
    !                  \4,4\4,5\4,6\.          |4,4|4,6|
    !                      \5,5\5,6\5,7|       |5,5|5,7|
    !                          \6,6\6,7|       |6,6| 0 |
    !                              \7,7|       |7,7| 0 |
    !                                 -|
    ! Take care: The virtual source matrix is zero-based in all dimensions, the
    ! target array is one-based, thus we have to transform the indices as
    ! well.
    if (trafo == ply_legToCheb_param) then

      ! Diagonals close to the main diagonal (complete length of the matrix)
      !NEC$ novector
      do i=1,diagonals
        diag_off = (i-1)*2
        ! Loop over all lines. The line end is shorter by 1 with every diagonal.
        ! As the matrix has 0 on every other diagonal, we only read every
        ! second diagonal, what leads to a diagonal shortening by 2.
        do m = 1, n - diag_off
          params%diag(m,i) = ply_m_int(m-1, diag_off + m-1)
        end do
      end do

      ! All adapter blocks between complete diagonals and staircase of
      ! approximated blocks
      !NEC$ novector
      do j=1,params%nBlocks-1
        block_off = (j-1)*s
        !NEC$ novector
        do i=1,blockdiagonals
          diag_off = (i-1)*2 + mod(remainder,2)
          ! Loop over all lines. The line end is shorter by 1 with every
          ! diagonal.
          ! As the matrix has 0 on every other diagonal, we only read every
          ! second diagonal, what leads to a diagonal shortening by 2.
          do m = 1, s - diag_off
            params%adapter(m,i,j) = ply_m_int(           &
              & block_off + m - 1,                       &
              & remainder + block_off + diag_off + m - 1 )
          end do
        end do
      end do

    else

      ! Diagonals close to the main diagonal (complete length of the matrix)
      !NEC$ novector
      do i=1,diagonals
        diag_off = (i-1)*2
        ! Loop over all lines. The line end is shorter by 1 with every diagonal.
        ! As far as the matrix has 0 on every other diagonal, we only read every
        ! second diagonal, what leads to a diagonal shortening by 2.
        do m = 1, n - diag_off
          params%diag(m,i) = ply_l_int(m-1, diag_off + m-1)
        end do
      end do

      ! All adapter blocks between complete diagonals and staircase of
      ! approximated blocks
      !NEC$ novector
      do j=1,params%nBlocks-1
        block_off = (j-1)*s
        !NEC$ novector
        do i=1,blockdiagonals
          diag_off = (i-1)*2 + mod(remainder,2)
          ! Loop over all lines. The line end is shorter by 1 with every
          ! diagonal.
          ! As the matrix has 0 on every other diagonal, we only read every
          ! second diagonal, what leads to a diagonal shortening by 2.
          do m = 1, s - diag_off
            params%adapter(m,i,j) = ply_l_int(block_off + m-1,                 &
              &                               remainder + block_off + diag_off &
              &                                         + m - 1)
          end do
        end do
      end do

    end if

    allocate(params%u(0:h,0:k-1))
    params%u(0:h,0:k-1) = u(0:h,0:k-1)
    allocate(params%sub(0:h))
    params%sub(0:h) = sub(0:h)
    params%n = n
    params%k = k
    params%s = s
    params%h = h

    ! Allocate the coefficients array for conversion
    allocate(params%b(0:params%h))
    !NEC$ novector
    do l=0, h
      nRows = (params%nBlocks - 1) / (2**l) - 1
      allocate(params%b(l)%col(2:nRows+1))
      !NEC$ novector
      do j=2,nRows+1
        allocate(params%b(l)%col(j)%coeff(0:k-1,0:1))
      end do
    end do

  end subroutine ply_fpt_init
  ! ************************************************************************ !


  ! ************************************************************************ !
  elemental function ply_m( iReal, jReal ) result( mVal )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: iReal, jReal
    real(kind=rk) :: mVal
    ! -------------------------------------------------------------------- !

    mVal = (2.0_rk/pi) * ply_lambda( 0.5_rk*(jReal - iReal) ) &
      &                * ply_lambda( 0.5_rk*(jReal + iReal) )

  end function ply_m
  ! ************************************************************************ !


  ! ************************************************************************ !
  elemental function ply_m_int( iReal, jReal ) result( mVal )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: iReal, jReal
    real(kind=rk) :: mVal
    ! -------------------------------------------------------------------- !

    if (mod(iReal+jReal,2).eq.0) then
      if (iReal.eq.0) then
        mVal = (1.0_rk/pi) * ply_lambda(0.5_rk*jReal)**2
      else
        mVal = (2.0_rk/pi) * ply_lambda( 0.5_rk*real(jReal-iReal, rk) ) &
          &                * ply_lambda( 0.5_rk*real(jReal+iReal, rk) )
      end if
    else
    !HK: deactivated to allow declaration as elemental function
    !HK!     write(*,*) iReal, jReal, 'error!'
    !HK!     stop
      mVal = 0.0_rk
    end if

  end function ply_m_int
  ! ************************************************************************ !


  ! ************************************************************************ !
  elemental function ply_l( iReal, jReal ) result( lVal )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: iReal, jReal
    real(kind=rk) :: lVal
    ! -------------------------------------------------------------------- !

    lVal = ( (-1.0_rk)*jReal*(iReal+0.5_rk)          &
      &     / ((jReal+iReal+1.0_rk)*(jReal-iReal)) ) &
      &  * ply_lambda(0.5_rk*(jReal-iReal-2.0_rk))   &
      &  * ply_lambda(0.5_rk*(jReal+iReal-1.0_rk))

  end function ply_l
  ! ************************************************************************ !


  ! ************************************************************************ !
  elemental function ply_l_int( i, j ) result( lVal )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: i, j
    real(kind=rk) :: lVal
    ! -------------------------------------------------------------------- !

    if (i.eq.0 .and. j.eq.0) then
      lVal = 1.0_rk
    elseif (i.eq.j) then
      lVal = sqrt(pi) / ( 2.0_rk * ply_lambda(real(i, rk)) )
    elseif (mod(i+j,2).eq.0) then
      lVal = ply_l(real(i,rk), real(j,rk))
    else
      lVal = 0.0_rk
    end if

  end function ply_l_int
  ! ************************************************************************ !


  ! ************************************************************************ !
  elemental function ply_lambda( val ) result( funcVal )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: val
    real(kind=rk) :: funcVal
    ! -------------------------------------------------------------------- !
    real(kind=rk), parameter :: lb_poly = 1.0_rk/0.067_rk
    real(kind=rk), parameter :: inter1 = 0.058_rk
    real(kind=rk), parameter :: inter2 = 0.04_rk
    real(kind=rk), parameter :: invbound = 0.02_rk
    real(kind=rk) :: invVal, a0, a1, a2, a3, a4, a5
    ! -------------------------------------------------------------------- !

    !>\todo: as we use a relation of gamma, it might be better to use the
    !!       gammln function provided by the numerical recipes, and just
    !!       use the difference in an exponential function.
    if (val <= lb_poly) then

      funcVal = Gamma(val+0.5_rk) / Gamma(val+1.0_rk)

    else

      invVal = 1.0_rk/val

      if (invVal >= inter1) then
        a0 =  0.99999999996378690_rk
        a1 = -0.12499999657282661_rk
        a2 =  0.78123659464717666e-02_rk
        a3 =  0.48855685911602214e-02_rk
        a4 = -0.67176366234107532e-03_rk
        a5 = -0.13533949520771154e-02_rk
      elseif (invVal >= inter2) then
        a0 =  0.99999999999298844_rk
        a1 = -0.12499999914033463_rk
        a2 =  0.78124565447111342e-02_rk
        a3 =  0.48839648427432678e-02_rk
        a4 = -0.65752249058053233e-03_rk
        a5 = -0.14041419931494052e-02_rk
      elseif (invVal >= invbound) then
        a0 =  0.99999999999974725_rk
        a1 = -0.12499999994706490_rk
        a2 =  0.78124954632722315e-02_rk
        a3 =  0.48830152125039076e-02_rk
        a4 = -0.64579205161159155e-03_rk
        a5 = -0.14628616278637035e-02_rk
      else
        a0 =  0.99999999999999999_rk
        a1 = -0.12499999999996888_rk
        a2 =  0.78124999819509772e-02_rk
        a3 =  0.48828163023526451e-02_rk
        a4 = -0.64122689844951054e-03_rk
        a5 = -0.15070098356496836e-02_rk
      end if

      funcVal = (((((a5*invVal + a4)*invVal &
        &                      + a3)*invVal &
        &                      + a2)*invVal &
        &                      + a1)*invVal &
        &                      + a0)*sqrt(invVal)

    end if

  end function ply_lambda
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_calculate_coeff_strip( n, s, gam, matrix, alph,             &
    &                                   nDiagonals, block_offset, remainder, &
    &                                   strip_lb, strip_ub, subblockingWidth )
    ! -------------------------------------------------------------------- !
    ! The overall number of modal coefficients
    integer, intent(in) :: n
    ! Size of the smallest block
    integer, intent(in) :: s
    !> Modal coefficients of the Chebyshev expansion.
    !! Size has to be: (1:indeps*params%n,nVars)
    !!
    !! Note, that the resulting array will have changed layout, and the
    !! transformed direction will run slowest in the array.
    real(kind=rk), intent(inout) :: gam(:)
    !> The arraz that holds the coefficients to calculate.
    real(kind=rk), intent(in) :: matrix(:,:)
    !> Modal coefficients of the Legendre expansion.
    !! Size has to be: (1:params%n*indeps,nVars)
    !!
    !! The direction which is to be transformed has to run fastest in
    !! the array.
    real(kind=rk), intent(in) :: alph(:)
    !> The number of diagonals to calculcate
    integer, intent(in) :: nDiagonals
    !> The offset of the block relative to the origin of the whole matrix.
    integer, intent(in) :: block_offset
    !> The diagonals that are not covered by any block.
    integer, intent(in) :: remainder
    !> The lower bound of the strip to calculate.
    integer, intent(in) :: strip_lb
    !> The upper bound of the strip to calculate.
    integer, intent(in) :: strip_ub
    !> The subblocking width defines the size of the blocking of the diagonal
    !! strides, i.e. in y direction. This subblocking is used to get a better
    !! cache locality.
    integer, intent(in) :: subblockingWidth
    ! -------------------------------------------------------------------- !
    ! The index for the line loop
    integer :: m
    ! The next index to calculate
    integer :: m_next
    integer :: iFun, iVal, indep
    integer :: row_off, col_off
    ! The index for the diagonal loop
    integer :: iDiag
    ! The next diagonal to calculate
    integer :: iDiag_next
    integer :: diag_off
    !> The start of the current block of m
    integer :: m_blocking
    ! -------------------------------------------------------------------- !
    iDiag_next = 1
    m = 1 ! default value for the index of the line loop

    ! Loop over all diagonals with a stride for unrolling. Stop when no more
    ! full stride can be calculated.
    !NEC$ novector
    do iDiag=1, nDiagonals-7, 8
      diag_off = (iDiag-1)*2 + mod(remainder,2)

      ! We are calculating two diagonals at once. Thus, we also have to stop two
      ! rows before the end. Every diagonal is one row shorter than the
      ! previous diagonal, and in matrix only every second diagonal is
      ! contained.
      ! We start with the first line and get back (in m_next) the first line
      ! that wasn't calculatable with the current loop unrolling length.
      !nec$ novector
      do m_blocking = 1, s-diag_off-14, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-14)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                  &
              & + matrix(m, idiag)     * alph(ifun + col_off)      &
              & + matrix(m, idiag + 1) * alph(ifun + col_off +  2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off +  4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off +  6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off +  8) &
              & + matrix(m, idiag + 5) * alph(ifun + col_off + 10) &
              & + matrix(m, idiag + 6) * alph(ifun + col_off + 12) &
              & + matrix(m, idiag + 7) * alph(ifun + col_off + 14)
          end do
        end do
      end do
      m_next = m
      ! This line in m_next is now the first line for the next loop with a
      ! smaller loop unrolling length.
      !nec$ novector
      do m_blocking = m_next, s-diag_off-12, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-12)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                  &
              & + matrix(m, idiag)     * alph(ifun + col_off)      &
              & + matrix(m, idiag + 1) * alph(ifun + col_off +  2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off +  4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off +  6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off +  8) &
              & + matrix(m, idiag + 5) * alph(ifun + col_off + 10) &
              & + matrix(m, idiag + 6) * alph(ifun + col_off + 12)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-10, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-10)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                  &
              & + matrix(m, idiag)     * alph(ifun + col_off)      &
              & + matrix(m, idiag + 1) * alph(ifun + col_off +  2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off +  4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off +  6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off +  8) &
              & + matrix(m, idiag + 5) * alph(ifun + col_off + 10)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-8, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-8)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off + 8)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-6, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-6)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-4, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-4)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-2, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-2)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival) + matrix(m, idiag) * alph(ifun + col_off)
          end do
        end do
      end do
    end do

    ! The next diagonal to calculate is the same as iDiag from the previous
    ! loop because the previous loop stops when it's stride is greater than the
    ! remaining diagonals. But as in the last loop all diagonals up to the
    ! current iDiag have been calculated, iDiag is the first diagonal that
    ! still has to be done.
    iDiag_next = iDiag

    !NEC$ novector
    do iDiag=iDiag_next, nDiagonals-6, 7
      diag_off = (iDiag-1)*2 + mod(remainder,2)
      !nec$ novector
      do m_blocking = 1, s-diag_off-12, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-12)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                  &
              & + matrix(m, idiag)     * alph(ifun + col_off)      &
              & + matrix(m, idiag + 1) * alph(ifun + col_off +  2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off +  4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off +  6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off +  8) &
              & + matrix(m, idiag + 5) * alph(ifun + col_off + 10) &
              & + matrix(m, idiag + 6) * alph(ifun + col_off + 12)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-10, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-10)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                  &
              & + matrix(m, idiag)     * alph(ifun + col_off)      &
              & + matrix(m, idiag + 1) * alph(ifun + col_off +  2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off +  4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off +  6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off +  8) &
              & + matrix(m, idiag + 5) * alph(ifun + col_off + 10)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-8, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-8)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off + 8)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-6, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-6)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-4, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-4)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-2, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-2)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival) + matrix(m, idiag) * alph(ifun + col_off)
          end do
        end do
      end do
    end do

    iDiag_next = iDiag

    !NEC$ novector
    do iDiag=iDiag_next, nDiagonals-5, 6
      diag_off = (iDiag-1)*2 + mod(remainder,2)
      !nec$ novector
      do m_blocking = 1, s-diag_off-10, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-10)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                  &
              & + matrix(m, idiag)     * alph(ifun + col_off)      &
              & + matrix(m, idiag + 1) * alph(ifun + col_off +  2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off +  4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off +  6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off +  8) &
              & + matrix(m, idiag + 5) * alph(ifun + col_off + 10)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-8, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-8)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off + 8)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-6, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-6)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-4, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-4)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-2, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-2)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival) + matrix(m, idiag) * alph(ifun + col_off)
          end do
        end do
      end do
    end do

    iDiag_next = iDiag

    !NEC$ novector
    do iDiag=iDiag_next, nDiagonals-4, 5
      diag_off = (iDiag-1)*2 + mod(remainder,2)
      !nec$ novector
      do m_blocking = 1, s-diag_off-8, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-8)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6) &
              & + matrix(m, idiag + 4) * alph(ifun + col_off + 8)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-6, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-6)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-4, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-4)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-2, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-2)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival) + matrix(m, idiag) * alph(ifun + col_off)
          end do
        end do
      end do
    end do

    iDiag_next = iDiag

    !NEC$ novector
    do iDiag=iDiag_next, nDiagonals-3, 4
      diag_off = (iDiag-1)*2 + mod(remainder,2)
      !nec$ novector
      do m_blocking = 1, s-diag_off-6, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-6)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4) &
              & + matrix(m, idiag + 3) * alph(ifun + col_off + 6)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-4, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-4)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-2, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-2)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival) + matrix(m, idiag) * alph(ifun + col_off)
          end do
        end do
      end do
    end do

    iDiag_next = iDiag

    !NEC$ novector
    do iDiag=iDiag_next, nDiagonals-2, 3
      diag_off = (iDiag-1)*2 + mod(remainder,2)
      !nec$ novector
      do m_blocking = 1, s-diag_off-4, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-4)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2) &
              & + matrix(m, idiag + 2) * alph(ifun + col_off + 4)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off-2, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-2)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival) + matrix(m, idiag) * alph(ifun + col_off)
          end do
        end do
      end do
    end do

    iDiag_next = iDiag
    do iDiag=iDiag_next, nDiagonals-1, 2
      !NEC$ novector
      diag_off = (iDiag-1)*2 + mod(remainder,2)
      !nec$ novector
      do m_blocking = 1, s-diag_off-2, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off-2)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival)                                 &
              & + matrix(m, idiag)     * alph(ifun + col_off)     &
              & + matrix(m, idiag + 1) * alph(ifun + col_off + 2)
          end do
        end do
      end do
      m_next = m
      !nec$ novector
      do m_blocking = m_next, s-diag_off, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival) + matrix(m, idiag) * alph(ifun + col_off)
          end do
        end do
      end do
    end do

    iDiag_next = iDiag

    !NEC$ novector
    do iDiag=iDiag_next, nDiagonals
      diag_off = (iDiag-1)*2 + mod(remainder,2)
      !nec$ novector
      do m_blocking = 1, s-diag_off, subblockingwidth
        !nec$ novector
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingwidth - 1, s-diag_off)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !nec$ ivdep
          do indep = strip_lb + 1, strip_ub
            ifun = (indep - 1) * n
            ival = ifun + row_off
            gam(ival) = gam(ival) + matrix(m, idiag) * alph(ifun + col_off)
          end do
        end do
      end do
    end do

  end subroutine ply_calculate_coeff_strip
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Convert strip of coefficients of a modal representation in terms of
  !! Legendre polynomials to modal coefficients in terms of Chebyshev
  !! polynomials.
  subroutine ply_fpt_exec( alph, gam, params, nIndeps )
    ! -------------------------------------------------------------------- !
    !> Number of values that can be computed independently.
    integer, intent(in) :: nIndeps

    !> Modal coefficients of the Legendre expansion.
    !! Size has to be: (1:params%n*nIndeps)
    !!
    !! The direction which is to be transformed has to run fastest in
    !! the array.
    real(kind=rk), intent(inout) :: alph(:)

    !> Modal coefficients of the Chebyshev expansion.
    !! Size has to be: (1:params%n*nIndeps)
    real(kind=rk), intent(out) :: gam(:)

    !> The parameters of the fast polynomial transformation.
    type(ply_trafo_params_type), intent(inout) :: params
    ! -------------------------------------------------------------------- !
    integer :: j, r, i, l, k, h, n, s, m
    integer :: iFun, indep
    integer :: iVal
    integer :: odd
    integer :: remainder
    integer :: nRows
    integer :: ub_row, row_rem
    integer :: rowsize
    integer :: block_off
    integer :: iBlock
    ! -------------------------------------------------------------------- !

    n = params%n
    k = params%k
    s = params%s
    h = params%h

    remainder = params%remainder

    !NEC$ novector
    indeploop: do indep = 1, nIndeps

      iFun = (indep-1)*n

      ! Set the output to zero
      gam(iFun+1:iFun+n) = 0.0_rk

      ! Calculate bs for all columns
      !NEC$ novector
      blockSizeLoop: do l = 0,h
        rowsize = s * 2**l
        nRows = (params%nBlocks - 1) / (2**l) - 1
        ub_row = 3 - mod(nRows,2)
        row_rem = mod(n-remainder, rowsize) + remainder + iFun
        !NEC$ novector
        blockColLoop: do j = 2, nRows+1, 1+mod(nRows,2)
          !NEC$ novector
          do r = 0, k-1
            params%b(l)%col(j)%coeff(r,0) = 0.0_rk
            params%b(l)%col(j)%coeff(r,1) = 0.0_rk
            !NEC$ novector
            do m = 0, rowsize-1
              odd = mod(row_rem + m + (j-1)*rowsize,2)
              params%b(l)%col(j)%coeff(r,odd) &
                &  = params%b(l)%col(j)%coeff(r,odd) &
                &    + params%u(l,r)%dat(m) &
                &      * alph(row_rem + m + (j-1)*rowsize + 1)
            end do
          end do
        end do blockColLoop

        ! Multiply with the blocks that are separated from the diagonal
        !NEC$ novector
        do i = 0, nRows - 1
          block_off = i*rowsize
          !NEC$ novector
          do j = i+2, i+ub_row - mod(i,2)
            !NEC$ novector
            do m = 0, rowsize - 1
              odd = mod(m+block_off,2)
              iVal = iFun + m + block_off+1
              !NEC$ novector
              do r = 0, k-1
                gam(iVal) = gam(iVal)      &
                  &       + params%sub(l)%subRow(i)%subCol(j)%rowDat(m)&
                  &               %coeff(r) &
                  &         * params%b(l)%col(j)%coeff(r,odd)
              end do ! r
            end do ! m
          end do ! j
        end do ! i

      end do blockSizeLoop

      if (params%trafo == ply_legToCheb_param) then
        ! Divide the first row in gam by 2, if we transform from Legendre
        ! to Chebyshev
        gam((indep-1)*n+1) = 0.5_rk*gam((indep-1)*n+1)
      end if

    end do indeploop

    ! Multiply with the entries near the diagonal
    call ply_calculate_coeff_strip(                &
      & n                = params%n,               &
      & s                = params%n,               &
      & gam              = gam,                    &
      & matrix           = params%diag,            &
      & alph             = alph,                   &
      & nDiagonals       = params%nDiagonals,      &
      & block_offset     = 0,                      &
      & remainder        = 0,                      &
      & strip_lb         = 0,                      &
      & strip_ub         = nIndeps,                &
      & subblockingWidth = params%subblockingWidth )

    ! Multiply with entries in the adapters
    !NEC$ novector
    do iBlock=1,params%nBlocks-1

      block_off = (iBlock-1)*params%s

      call ply_calculate_coeff_strip(                    &
        & n                = params%n,                   &
        & s                = params%s,                   &
        & gam              = gam,                        &
        & matrix           = params%adapter(:,:,iBlock), &
        & alph             = alph,                       &
        & nDiagonals       = params%nBlockDiagonals,     &
        & block_offset     = block_off,                  &
        & remainder        = remainder,                  &
        & strip_lb         = 0,                          &
        & strip_ub         = nIndeps,                    &
        & subblockingWidth = params%subblockingWidth     )

    end do

  end subroutine ply_fpt_exec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Convert strip of coefficients of a modal representation in terms of
  !! Legendre polynomials to modal coefficients in terms of Chebyshev
  !! polynomials.
  subroutine ply_fpt_single( alph, gam, params )
    ! -------------------------------------------------------------------- !
    !> The parameters of the fast polynomial transformation.
    type(ply_trafo_params_type), intent(inout) :: params

    !> Modal coefficients of the Legendre expansion.
    !! Size has to be: params%n
    !!
    !! The direction which is to be transformed has to run fastest in
    !! the array.
    real(kind=rk), intent(inout) :: alph(params%n)

    !> Modal coefficients of the Chebyshev expansion.
    !! Size has to be: params%n
    real(kind=rk), intent(out) :: gam(params%n)
    ! -------------------------------------------------------------------- !
    integer :: j, r, i, l, k, h, n, s, m
    integer :: iVal
    integer :: odd
    integer :: remainder
    integer :: nRows
    integer :: ub_row, row_rem
    integer :: rowsize
    integer :: block_off
    integer :: iBlock
    ! -------------------------------------------------------------------- !

    n = params%n
    k = params%k
    s = params%s
    h = params%h

    remainder = params%remainder

    ! Set the output to zero
    gam = 0.0_rk

    ! Calculate bs for all columns
    !NEC$ novector
    blockSizeLoop: do l = 0,h
      rowsize = s * 2**l
      nRows = (params%nBlocks - 1) / (2**l) - 1
      ub_row = 3 - mod(nRows,2)
      row_rem = mod(n-remainder, rowsize) + remainder
      !NEC$ novector
      blockColLoop: do j = 2, nRows+1, 1+mod(nRows,2)
        do r = 0, k-1
          params%b(l)%col(j)%coeff(r,0) = 0.0_rk
          params%b(l)%col(j)%coeff(r,1) = 0.0_rk
        end do
        !NEC$ novector
        do m = 0, rowsize-1
          odd = mod(row_rem + m + (j-1)*rowsize,2)
          do r = 0, k-1
            params%b(l)%col(j)%coeff(r,odd) &
              &  = params%b(l)%col(j)%coeff(r,odd) &
              &    + params%u(l,r)%dat(m) &
              &      * alph(row_rem + m + (j-1)*rowsize + 1)
          end do
        end do
      end do blockColLoop

      ! Multiply with the blocks that are separated from the diagonal
      !NEC$ novector
      do i = 0, nRows - 1
        block_off = i*rowsize
        !NEC$ novector
        do j = i+2, i+ub_row - mod(i,2)
          !NEC$ novector
          do m = 0, rowsize - 1
            odd = mod(m+block_off,2)
            iVal = m + block_off+1
            do r = 0, k-1
              gam(iVal) = gam(iVal)      &
                &       + params%sub(l)%subRow(i)%subCol(j)%rowDat(m)&
                &               %coeff(r) &
                &         * params%b(l)%col(j)%coeff(r,odd)
            end do ! r
          end do ! m
        end do ! j
      end do ! i

    end do blockSizeLoop

    if (params%trafo == ply_legToCheb_param) then
      ! Divide the first row in gam by 2, if we transform from Legendre
      ! to Chebyshev
      gam(1) = 0.5_rk*gam(1)
    end if

    ! Multiply with the entries near the diagonal
    call ply_calculate_coeff_strip(                &
      & n                = params%n,               &
      & s                = params%n,               &
      & gam              = gam,                    &
      & matrix           = params%diag,            &
      & alph             = alph,                   &
      & nDiagonals       = params%nDiagonals,      &
      & block_offset     = 0,                      &
      & remainder        = 0,                      &
      & strip_lb         = 0,                      &
      & strip_ub         = 1,                      &
      & subblockingWidth = params%subblockingWidth )

    ! Multiply with entries in the adapters
    !NEC$ novector
    do iBlock=1,params%nBlocks-1

      block_off = (iBlock-1)*params%s

      call ply_calculate_coeff_strip(                    &
        & n                = params%n,                   &
        & s                = params%s,                   &
        & gam              = gam,                        &
        & matrix           = params%adapter(:,:,iBlock), &
        & alph             = alph,                       &
        & nDiagonals       = params%nBlockDiagonals,     &
        & block_offset     = block_off,                  &
        & remainder        = remainder,                  &
        & strip_lb         = 0,                          &
        & strip_ub         = 1,                          &
        & subblockingWidth = params%subblockingWidth     )

    end do

  end subroutine ply_fpt_single
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Convert coefficients of a modal representation in terms of Legendre
  !! polynomials to modal coefficients in terms of Chebyshev polynomials.
  subroutine ply_fpt_exec_striped( nIndeps, alph, gam, params )
    ! -------------------------------------------------------------------- !
    !> Number of values that can be computed independently.
    integer, intent(in) :: nIndeps

    !> The parameters of the fast polynomial transformation.
    type(ply_trafo_params_type), intent(inout) :: params

    !> Modal coefficients of the Legendre expansion.
    !! Size has to be: (1:params%n*indeps,nVars)
    !!
    !! The direction which is to be transformed has to run fastest in
    !! the array.
    real(kind=rk), intent(in) :: alph(:)

    !> Modal coefficients of the Chebyshev expansion.
    !! Size has to be: (1:indeps*params%n,nVars)
    !!
    !! Note, that the resulting array will have changed layout, and the
    !! transformed direction will run slowest in the array.
    real(kind=rk), intent(out) :: gam(:)
    ! -------------------------------------------------------------------- !
    integer :: j, r, i, l, k, h, n, s, m
    integer :: iStrip, iFun, indep
    integer :: odd
    integer :: strip_ub
    integer :: striplen
    integer :: remainder
    integer :: nRows
    integer :: ub_row, row_rem(params%striplen)
    integer :: stripend
    integer :: rowsize
    integer :: block_off
    integer :: stripoff
    integer :: iBlock
    ! -------------------------------------------------------------------- !

    n = params%n
    k = params%k
    s = params%s
    h = params%h
    striplen = params%striplen

    remainder = params%remainder

    ! Set the output to zero
    gam = 0.0_rk

    ! Loop over all strips
    !NEC$ novector
    do iStrip = 0,nIndeps-1,striplen
      ! Calculate the upper bound of the current strip
      strip_ub = min(iStrip + striplen, nIndeps)
      stripend = strip_ub - iStrip

      ! Calculate bs for all columns
      !NEC$ novector
      blockSizeLoop: do l =0,h
        rowsize = s * 2**l
        nRows = (params%nBlocks - 1) / (2**l) - 1
        ub_row = 3 - mod(nRows,2)
        do indep = 1, stripend
          iFun = (indep+iStrip-1)*params%n
          row_rem(indep) = mod(n-remainder, rowsize) + remainder + iFun
        end do

        !NEC$ novector
        blockColLoop: do j = 2, nRows+1, 1+mod(nRows,2)
          !NEC$ novector
          do r = 0, k-1
            params%b(l)%col(j)%coeff(r,0) = 0.0_rk
            params%b(l)%col(j)%coeff(r,1) = 0.0_rk
            !NEC$ novector
            do m = 0, rowsize-1
              do indep = 1, stripend
                odd = mod(row_rem(indep) + m + (j-1)*rowsize,2)
                params%b(l)%col(j)%coeff(r,odd) &
                  &  = params%b(l)%col(j)%coeff(r,odd) &
                  &    + params%u(l,r)%dat(m) &
                  &      * alph(row_rem(indep) + m + (j-1)*rowsize + 1)
              end do
            end do
          end do
        end do blockColLoop

        ! Multiply with the blocks that are separated from the diagonal
        !NEC$ novector
        do i = 0, nRows - 1
          block_off = i*rowsize
          !NEC$ novector
          do j = i+2, i+ub_row - mod(i,2)
            !NEC$ novector
            do m = 0, rowsize - 1
              odd = mod(m+block_off,2)
              stripoff = (m+block_off)*nIndeps + iStrip
              !NEC$ novector
              do r = 0, k-1
                do indep = 1, stripend
                  gam(indep+stripoff) = gam(indep+stripoff)                           &
                    &                   + params%sub(l)%subRow(i)%subCol(j)%rowDat(m) &
                    &                           %coeff(r)                             &
                    &                     * params%b(l)%col(j)%coeff(r,odd)
                end do
              end do ! r
            end do ! m
          end do ! j
        end do ! i

      end do blockSizeLoop

      if (params%trafo == ply_legToCheb_param) then
        ! Divide the first row in gam by 2, if we transform from legendre
        ! to chebyshev
        do indep = iStrip+1, strip_ub
          gam(indep) = 0.5_rk*gam(indep)
        end do
      end if

      ! Multiply with the entries near the diagonal
      call ply_calculate_coeff_strip(                &
        & n                = params%n,               &
        & s                = params%n,               &
        & gam              = gam,                    &
        & matrix           = params%diag,            &
        & alph             = alph,                   &
        & nDiagonals       = params%nDiagonals,      &
        & block_offset     = 0,                      &
        & remainder        = 0,                      &
        & strip_lb         = iStrip,                 &
        & strip_ub         = strip_ub,               &
        & subblockingWidth = params%subblockingWidth )


      ! Multiply with entries in the adapters
      do iBlock=1,params%nBlocks-1

        block_off = (iBlock-1)*params%s

        call ply_calculate_coeff_strip(                    &
          & n                = params%n,                   &
          & s                = params%s,                   &
          & gam              = gam,                        &
          & matrix           = params%adapter(:,:,iBlock), &
          & alph             = alph,                       &
          & nDiagonals       = params%nBlockDiagonals,     &
          & block_offset     = block_off,                  &
          & remainder        = remainder,                  &
          & strip_lb         = iStrip,                     &
          & strip_ub         = strip_ub,                   &
          & subblockingWidth = params%subblockingWidth     )

      end do
    end do ! iStrip

  end subroutine ply_fpt_exec_striped
  ! ************************************************************************ !


end module ply_polyBaseExc_module
