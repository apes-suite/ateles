! Copyright (c) 2017,2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2018-2019 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach and Peter Vitt for
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

!> This module provides the functionality to split Legendre polynomials into
!! a left and right subinterval with transformed coordinates.
!!
!! The original polynomial is defined on the interval [-1,1] and the two
!! new polynomial representations are computed in the intervals [-1,0] and
!! [0,1] but in the changed coordinate system, with the interval [-1,1] for
!! each.
!! Thus, if we refer to the coordinates in the original (coarse) element as
!! x and to the respective coordinates in the two halfed elements as xi_left
!! and xi_right, we compute the modal representation of the original Legendre
!! polynomial series under these transformations:
!! \[ x = 0.5 \cdot \xi_{left} - 0.5 \]
!! \[ x = 0.5 \cdot \xi_{right} + 0.5 \]
!!
!! This is needed when refining elements.
module ply_split_legendre_module
  use env_module, only: rk

  implicit none

  private

  public :: ply_split_legendre_matrix
  public :: ply_split_legendre_test


contains


  ! ------------------------------------------------------------------------ !
  !> Compute the transformation matrix for a projection to the left and right
  !! half-interval of Legendre polynomials for the given maximal number of
  !! modes.
  !!
  !! Note: The transformation matrices to each subinterval are triangular, and
  !!       the diagonal entries are the same. To save memory both matrices are
  !!       stored in a single 2 dimensional array of size
  !!       (nModes, nModes).
  !!
  !! This matrix only needs to be computed once for a sufficiently high order,
  !! as submatices out of it can by used to perform the transformation for
  !! any lower polynomial degree.
  !!
  !! The upper triangular matrix is created for the right subinterval,
  !! while the lower triangular matrix is used to store the rotated version
  !! for the left subinterval.
  !! For the right interval we interpret the first index as row index
  !! and the second as column. For the left interval this is reverted and
  !! we interpret the first index as columns of the matrix.
  !!
  !>@note Why is this a function? The reasoning for making this a function
  !!      is that we need to return exactly one thing (the split matrix).
  !!      It is then quite natural to refer to this by
  !!      ply_split_legendre_matrix. A subroutine on the other hand usually
  !!      describes something that should be done. Thus the name for a
  !!      subroutine would then be ply_split_legendre_compute_matrix
  !!      (describing the action performed by the subroutine).
  !!      When using OpenMP it sometimes is better to use subroutines, even
  !!      though it would be more natural to use a function. However, here
  !!      we do not expect this to be the case, as this is expected to be
  !!      called only once.
  !!@endnote
  pure function ply_split_legendre_matrix(nModes) result(split_matrix)
    ! -------------------------------------------------------------------- !
    !> The maximal number of modes to compute the transformation for.
    !!
    !! The resulting matrix v will be max_modes x max_modes large and can
    !! be used for the transformation of all polynomials with up to this
    !! many modes.
    integer, intent(in) :: nModes

    real(kind=rk) :: split_matrix(nModes, nModes)
    ! -------------------------------------------------------------------- !
    integer :: m
    integer :: orig
    integer :: sign_factor
    ! We only consider the split into the two halves. Thus, the according
    ! coordinate transformation factors in x = scaling*xi + shifting, are
    ! constant.
    real(kind=rk), parameter :: shifting = 0.5_rk
    real(kind=rk), parameter :: scaling = 0.5_rk
    ! -------------------------------------------------------------------- !

    ! The split matrix looks like this:
    ! indicing: split_matrix(row, column) for the right half.
    !
    ! (1,:)   [1.0  --  --     shift=0.5   ]
    ! (2,:)   | |  0.5  --                 |
    ! (3,:)   | |   |  0.25             ...|
    ! (4,:)   |             0.125          |
    ! (i,:)   | shift=-0.5        0.0625   |
    !   :     [    :                    ...]
    !
    ! To compute the matrix we use the Clenshaw algorithm to iteratively
    ! compute the modal representation step by step for the right halfinterval
    ! and we store each step as a column in the upper triangular matrix.
    ! The left halfinterval can then easily obtained from that due to the
    ! symmetry of this problem, as only the sign in the shifting is changed.

    ! The first step (column) is always just 1 in the first mode (row).
    if (nModes > 0) split_matrix(1,1) = 1.0_rk

    atleast2: if (nModes > 1) then

      ! In the next step (second column) we multiply the previous column
      ! by scaling and shift all modes one up. (There would also be a
      ! down-shifting but we now that all higher modes of the previous step are
      ! 0 anyway.)
      ! Then we add the shifting in the first mode. For the second step
      ! we can easily write this down explicitly:
      split_matrix(1,2) = shifting
      split_matrix(2,2) = scaling

      atleast3: if (nModes > 2) then

        ! For all higher modes we now actually need to compute something.
        do orig=3,nModes
          m = 1
          ! Skip terms from below 1 (we shift 0 in).
          split_matrix(m, orig) =            beta(orig-1)                &
            &                                  * split_matrix(m, orig-2) &
            &                   + shifting * alpha(orig-1)               &
            &                                  * split_matrix(m, orig-1) &
            &                   - scaling  * alpha_beta(m+1, orig-1)     &
            &                                 * split_matrix(m+1, orig-1)
          do m=2,orig-2
            split_matrix(m, orig) =            beta(orig-1)                  &
              &                                  * split_matrix(m, orig-2)   &
              &                   + shifting * alpha(orig-1)                 &
              &                                  * split_matrix(m, orig-1)   &
              &                   - scaling  * alpha_beta(m+1, orig-1)       &
              &                                  * split_matrix(m+1, orig-1) &
              &                   + scaling  * alpha_frac(m-1, orig-1)       &
              &                              * split_matrix(m-1, orig-1)
          end do
          ! Skip all terms from beyond the diagonal (they are 0).
          m = orig-1
          split_matrix(m, orig) = shifting * alpha(orig-1)                 &
            &                                  * split_matrix(m, orig-1)   &
            &                   + scaling  * alpha_frac(m-1, orig-1)       &
            &                              * split_matrix(m-1, orig-1)
          m = orig
          split_matrix(m, orig) = scaling  * alpha_frac(m-1, orig-1) &
            &                              * split_matrix(m-1,orig-1)
        end do

      end if atleast3

      ! Due to the symmetry of the problem (the left subinterval has just
      ! the shifting with a changed sign), we can fill the other half of
      ! the matrix by copying the already computed values accordingly with
      ! a change in sign, as needed (alternatingly).
      do orig=1,nModes
        sign_factor = mod(orig,2)*2 - 1
        !NEC$ ivdep
        do m=1, orig-1, 2
          split_matrix(orig, m) = sign_factor * split_matrix(m, orig)
        end do
        !NEC$ ivdep
        do m=2, orig-1, 2
          split_matrix(orig, m) = -sign_factor * split_matrix(m, orig)
        end do
      end do

    end if atleast2

  end function ply_split_legendre_matrix
  ! ======================================================================== !


  ! !!!!!!! !
  ! private !
  ! !!!!!!! !

  ! ------------------------------------------------------------------------ !
  !> Coefficient alpha from the recursive formulation of Legendre polynomials,
  !! for the Legendre mode 'mode'.
  !!
  !! \[L_n(x) = \alpha \cdot x \cdot L_{n-1}(x) + \beta \cdot L_{n-2}(x)\]
  !!
  !! For Legendre polynomials we have:
  !!
  !! \[\alpha = \frac{2 n - 1}{n}\]
  elemental function alpha(mode)
    ! -------------------------------------------------------------------- !
    !> The Legendre mode to compute \(\alpha\) for.
    integer, intent(in) :: mode

    real(kind=rk) :: alpha
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    if (mode > 0) then
      alpha = real((2 * mode - 1), kind=rk) / real(mode, kind=rk)
    else
      alpha = 0.0_rk
    end if

  end function alpha
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Coefficient beta from the recursive formulation of Legendre polynomials,
  !! for the Legendre mode 'mode'.
  !!
  !! \[L_n(x) = \alpha \cdot x \cdot L_{n-1}(x) + \beta \cdot L_{n-2}(x)\]
  !!
  !! For Legendre polynomials we have:
  !!
  !! \[\beta = \frac{1 - n}{n}\]
  !!
  !>@note This is negative for all modes > 1.
  elemental function beta(mode)
    ! -------------------------------------------------------------------- !
    !> The Legendre mode to compute \(\beta\) for.
    integer, intent(in) :: mode

    real(kind=rk) :: beta
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    if (mode > 1) then
      beta = real((1 - mode), kind=rk) / real(mode, kind=rk)
    else
      beta = 0.0_rk
    end if

  end function beta
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Quotient of two alpha values.
  !!
  !! This function computes alpha(numerator)/alpha(denominator).
  !!
  !>@note This is intended to keep as many integer operations together as
  !!      possible.
  elemental function alpha_frac(denominator, numerator)
    ! -------------------------------------------------------------------- !
    !> Legendre mode of the \(\alpha\) to use in the denominator.
    integer, intent(in) :: denominator

    !> Legendre mode of the \(\alpha\) to use in the numeratorr.
    integer, intent(in) :: numerator

    real(kind=rk) :: alpha_frac
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    if ( denominator > 0 .and. numerator > 0 ) then
      if ( denominator == numerator ) then
        alpha_frac = 1.0_rk
      else
        alpha_frac = real((2*numerator-1)*denominator, kind=rk)  &
          &          / real((2*denominator-1)*numerator, kind=rk)
      end if
    else
      alpha_frac = 0.0_rk
    end if
  end function alpha_frac
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Prodcut of alpha(numerator) * beta(denominator) / alpha(denominator) as
  !! needed by the Clenshaw algorithm in ply_split_legendre_matrix.
  elemental function alpha_beta(denominator, numerator)
    ! -------------------------------------------------------------------- !
    !> Legendre mode for the \(\alpha\) in the numerator.
    integer, intent(in) :: numerator

    !> Legendre mode for the \(\alpha\) in the denominator and the \(\beta\).
    integer, intent(in) :: denominator

    real(kind=rk) :: alpha_beta
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    if (numerator > 0 .and. denominator > 0) then
      if ( denominator == numerator ) then
        alpha_beta = real((1-denominator), kind=rk) &
          &          / real(denominator, kind=rk)
      else
        alpha_beta = real((2*numerator-1)*(1-denominator), kind=rk) &
          &          / real(numerator*(2*denominator-1), kind=rk)
      end if
    else
      alpha_beta = 0.0_rk
    end if

  end function alpha_beta
  ! ======================================================================== !


  ! !!!!!!! !
  ! testing !
  ! !!!!!!! !

  ! ------------------------------------------------------------------------ !
  !> A small testing routine to check the functions of this module.
  !!
  subroutine ply_split_legendre_test(success)
    ! -------------------------------------------------------------------- !
    !> Indication whether the tests were completed successfully.
    logical, intent(out) :: success
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: large_mat(:,:)
    real(kind=rk), allocatable :: small_mat(:,:)
    real(kind=rk) :: res
    integer :: iiter
    integer :: jiter
    ! -------------------------------------------------------------------- !

    success = .true.

    ! Expected values for alpha:
    res = alpha(1)
    if (abs(res - 1.0_rk)  > epsilon(res)) success = .false.
    res = alpha(2)
    if (abs(res - 1.5_rk)  > epsilon(res)) success = .false.
    res = alpha(4)
    if (abs(res - 1.75_rk) > epsilon(res)) success = .false.
    res = alpha(5)
    if (abs(res - 1.8_rk)  > epsilon(res)) success = .false.
    res = alpha(10)
    if (abs(res - 1.9_rk)  > epsilon(res)) success = .false.
    res = alpha(100)
    if (abs(res - 1.99_rk) > epsilon(res)) success = .false.
    if (.not. success) then
      write(*,*) 'Alpha check failed'
      RETURN
    end if

    ! Expected values for beta:
    res = beta(1)
    if (abs(res) > tiny(res)) success = .false.
    res = beta(2)
    if (abs(res + 0.5_rk)  > epsilon(res)) success = .false.
    res = beta(4)
    if (abs(res + 0.75_rk) > epsilon(res)) success = .false.
    res = beta(5)
    if (abs(res + 0.8_rk)  > epsilon(res)) success = .false.
    res = beta(10)
    if (abs(res + 0.9_rk)  > epsilon(res)) success = .false.
    res = beta(100)
    if (abs(res + 0.99_rk) > epsilon(res)) success = .false.
    if (.not. success) then
      write(*,*) 'Beta check failed'
      RETURN
    end if

    ! Expected properties for alpha_frac:
    ! * It should be the same as dividing the corresponding values of alpha
    ! * If both modes are identic, it should be one
    ! * If the modes are exchanged the value should be inversed
    do iiter=1,100
      res = alpha_frac(iiter, iiter)
      if (abs(res - 1.0_rk)  > epsilon(res)) success = .false.
      do jiter=iiter+1,100
        res = alpha_frac(iiter, jiter)
        if (abs(res - (alpha(jiter)/alpha(iiter))) > epsilon(res)) &
          & success = .false.
        if (abs(res - (1._rk/alpha_frac(jiter, iiter))) > epsilon(res)) &
          & success = .false.
      end do
    end do
    ! Some expected values for alpha_frac:
    res = alpha_frac(1,2)
    if (abs(res - 1.5_rk)  > epsilon(res)) success = .false.
    res = alpha_frac(1,4)
    if (abs(res - 1.75_rk) > epsilon(res)) success = .false.
    res = alpha_frac(1,5)
    if (abs(res - 1.8_rk)  > epsilon(res)) success = .false.
    res = alpha_frac(2,4)
    if (abs(res - (1.75_rk/1.5_rk)) > epsilon(res)) success = .false.
    res = alpha_frac(2,5)
    if (abs(res - (1.8_rk/1.5_rk))  > epsilon(res)) success = .false.
    res = alpha_frac(4,5)
    if (abs(res - (1.8_rk/1.75_rk)) > epsilon(res)) success = .false.
    if (.not. success) then
      write(*,*) 'alpha_frac check failed'
      RETURN
    end if

    ! Expected properties for alpha_beta:
    ! * It should be the same as when using alpha and beta accordingly
    ! * If both modes are identic, it should be the same as beta(mode)
    do iiter=1,100
      res = alpha_beta(iiter, iiter)
      if (abs(res - beta(iiter)) > epsilon(res)) success = .false.
      do jiter=iiter+1,100
        res = alpha_beta(iiter, jiter)
        if ( abs(res - ((alpha(jiter)*beta(iiter))/alpha(iiter))) &
          &  > epsilon(res) ) success = .false.
      end do
    end do
    ! Some expected values for alpha_beta:
    res = alpha_beta(1,2)
    if (abs(res) > tiny(res)) success = .false.
    res = alpha_beta(1,4)
    if (abs(res) > tiny(res)) success = .false.
    res = alpha_beta(1,5)
    if (abs(res) > tiny(res)) success = .false.
    res = alpha_beta(2,4)
    if (abs(res + (7.0_rk/12.0_rk))  > epsilon(res)) success = .false.
    res = alpha_beta(2,5)
    if (abs(res + 0.6_rk)            > epsilon(res)) success = .false.
    res = alpha_beta(4,5)
    if (abs(res + (27.0_rk/35.0_rk)) > epsilon(res)) success = .false.
    if (.not. success) then
      write(*,*) 'alpha_beta check failed'
      RETURN
    end if


    ! Some expected properties of the ply_split_legendre_matrix function
    allocate(large_mat(100,100))
    large_mat = ply_split_legendre_matrix(100)

    ! * The diagonal should be 0.5**(row-1)
    do iiter=1,100
      if ((large_mat(iiter, iiter) - 0.5_rk**(iiter-1)) > epsilon(res)) &
        & success = .false.
    end do

    ! * For smaller modes, we should just get submatrices of the larger one.
    do iiter=1,99
      allocate(small_mat(iiter, iiter))
      small_mat =  ply_split_legendre_matrix(iiter)
      if (maxval(abs(large_mat(:iiter, :iiter) - small_mat)) > epsilon(res)) &
        & success = .false.
      deallocate(small_mat)
    end do
    if (.not. success) then
      write(*,*) 'ply_split_legendre_matrix check failed'
      RETURN
    end if

  end subroutine ply_split_legendre_test
  ! ======================================================================== !

end module ply_split_legendre_module
