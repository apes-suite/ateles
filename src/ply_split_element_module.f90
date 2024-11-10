! Copyright (c) 2017, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach for University of Siegen.
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
!> This module provides the methods to project the polynomial representation in
!! elements onto the representations in their halves in each dimension.
!!
!! To perform the projection for Legendre polynomials we will use the computed
!! coefficients for the Clenshaw algorithm from [[ply_split_legendre_module]].
!! With those the transformation is just a simple triangular matrix
!! multiplication, but we need to take care of the orthogonal degrees of freedom
!! as we want to handle all of them at the same time.
!! Further we want to allow the transformation to be performed for multiple
!! elements at once.
!!
!! In each dimension we need to perform the following coordinate transformation:
!!
!! \[ x = 0.5 \cdot \xi_{left} - 0.5 \]
!! \[ x = 0.5 \cdot \xi_{right} + 0.5 \]
!!
!! Where \(x\) refers to the coordinate in the original (coarse) element, and
!! \(\xi\) to the coordinates in the two (left and right) halves of the element.
module ply_split_element_module
  use env_module,                 only: rk
  use ply_split_legendre_module,  only: ply_split_legendre_matrix
  use ply_modg_basis_module,      only: ply_legendre_1d

  implicit none

  private

  public :: ply_split_element_singleD
  public :: ply_split_element
  public :: ply_split_element_1D
  public :: ply_split_element_2D
  public :: ply_split_element_3D
  public :: ply_split_element_init

  public :: ply_split_element_test

  abstract interface
    !> Split elements of degree parent_degree into elements with polynomials of
    !! degree child_degree.
    subroutine ply_split_element( parent_degree, child_degree, parent_data, &
      &                           child_data, ignore_highmodes )
      ! -------------------------------------------------------------------- !
      import :: rk
      !> Polynomial degree in the parent element.
      integer, intent(in) :: parent_degree

      !> Polynomial degree in the child elements.
      integer, intent(in) :: child_degree

      !> Polynomial data in the parent element. The first index describes the
      !! degrees of freedom. The second index refers to the elements to split.
      real(kind=rk), intent(in) :: parent_data(:,:)

      !> Polynomial data in the child elements. The first index describes the
      !! degrees of freedom. The second index refers to the elements, there
      !! needs to be four times as many elements than in the parent_data.
      !!
      !! Elements follow the ordering of the Z space filling curve.
      real(kind=rk), intent(out) :: child_data(:,:)

      !> Whether to ignore high modes from the parent element.
      !!
      !! This can be used as a simple lowpass filter by ignoring all higher
      !! modes from the parent element, that exceed the target polynomial
      !! degree. Thus, the polynomials are filtered before projection,
      !! instead of cutting them only of after refinement.
      !! Defaults to false (no filtering).
      logical, optional, intent(in) :: ignore_highmodes
      ! -------------------------------------------------------------------- !
    end subroutine ply_split_element
  end interface

  !> Precomputed matrix to hold the transformation operation to project
  !! Legendre polynomials to its two half intervals.
  !!
  !! This is computed by [[ply_split_legendre_matrix]], see there for details.
  !! There are two triangular matrices stored in this array, one for the
  !! projection to the left half (-1,0) , and one for the projection to the
  !! right half (0,1).
  !!
  !! This is a module variable, as it is only needed to be computed once with
  !! sufficient size. All lower orders are just subarrays out of the larger one.
  real(kind=rk), allocatable :: split_legendre(:,:)


contains


  ! ------------------------------------------------------------------------ !
  !> Initialization of the module.
  !! This needs to be performed before any call of the actual transformation
  !! [[ply_split_element_1D]].
  !!
  !! The initialization will compute the transformation matrix for Legendre
  !! polynomials with at least nMaxModes. If the initialization was already
  !! called before with the same or larger nMaxModes, the matrix will not be
  !! changed. Thus, calling this routine will only increase the size of the
  !! module variable split_legendre, never decrease it.
  subroutine ply_split_element_init(nMaxModes)
    ! -------------------------------------------------------------------- !
    !> Maximal number of expected modes to perform the splitting for.
    integer, intent(in) :: nMaxModes
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    if (allocated(split_legendre)) then
      if (size(split_legendre, 1) < nMaxModes) deallocate(split_legendre)
    end if

    if (.not. allocated(split_legendre)) then
      allocate(split_legendre(nMaxModes, nMaxModes))
      split_legendre = ply_split_legendre_matrix(nMaxModes)
    end if

  end subroutine ply_split_element_init
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Project a polynomial representation in elements in one dimension to its
  !! two halves in that direction.
  !!
  !! For each parent element the projection on the two respective child elements
  !! (half intervals) are computed for one dimension.
  !!
  !>@note Preliminary data layout and interface planning.
  !! It might be that we should rather split the index into the direction
  !! in which we perform the operation and all the other directions normal
  !! to that. For a dense matrix this may allow the compiler to detect
  !! the matrix multiply. However, here for the triangular matrix it is not
  !! so sure, whether this would be possible.
  !!@endnote
  !!
  !>@note After discussions with Stephan Walter, it looks like the separate
  !! indices would most likely be better.
  !! Maybe, using explicit shaped arrays and therby allowing more dimensions
  !! in the input, while keeping the interface to two dimensions for all
  !! cases (the normal direction and all independent degrees of freedom).
  !! For vectorization on x86 it also is necessary to have a stride-1 access
  !! only in reading and writing.
  !! The rotation of data might not be the best option because of this.
  !! Instead, it may be that we need to have different routines for each
  !! direction.
  !! Or, maybe, we need to use the elements as first index and vectorize
  !! over those.
  !!@endnote
  !!
  !! As we need to perform this operation in all dimensions, it would be good
  !! to shift the indices around. When doing this, we can stick to the same
  !! implementation for all directions, without the need to put any logic in
  !! here to decide on the current direction.
  !! In 3D we would end up with this chain:
  !! (x,y,z) -> split_element for Z -> (z,x,y)
  !!         -> split_element for Y -> (y,z,x)
  !!         -> split_element for X -> (x,y,z)
  !! Thus, the logic is that we perform the split on the last dimension, and
  !! cycle the indices in the output.
  !!
  !! We can generalize this to arbitrary dimensions.
  !! In 2D it would look like this:
  !! (x,y) -> split_element for Y -> (y,x)
  !!       -> split_element for X -> (x,y)
  !! And in 1D, we just need to perform one transformation:
  !! (x) -> split_element for X -> (x)
  !!
  !! As we allow for a changed number of polynomial degrees in the input and
  !! output, we need to take care of different lengths for each direction.
  !! Thus, we need: the dimensionality, two 1D arrays with the length of this
  !! dimensionality to provide the number of degrees of freedom for each
  !! direction (once for the input, and once for the output).
  !!
  !! We need: nDofs in the direction where the transformation is to be done
  !!          and the nDofs for all normal directions.
  subroutine ply_split_element_singleD( nDims, inLen, outLen, parent_data, &
    &                                   child_data, ignore                 )
    ! -------------------------------------------------------------------- !
    !> Number of dimensions of the polynomial data.
    integer, intent(in) :: nDims

    !> Number degrees of freedom for each direction in parent_Data.
    !!
    !! The first index of parent_data needs to have a length equal to the
    !! product of all inLen components.
    !! The splitting operation will be done in the last dimension.
    integer, intent(in) :: inLen(nDims)

    !> Number degrees of freedom for each direction in child_Data.
    !!
    !! The first index of child_data needs to have a length equal to the
    !! product of all outLen components.
    !! The data will be cyclicly exchanged. Thus, the last dimension in
    !! parent_data corresponds to the first in one in child_data and all
    !! other components are shifted once to the right.
    integer, intent(in) :: outLen(nDims)

    !> Polynomial representation in the parent elements.
    !!
    !! The first index are the degrees of freedom in elements, the second index
    !! are the elements.
    !! In the first index the shape of data has to be in the form
    !! (inLen(1), inLen(2), ... , inLen(nDims)).
    !! The splitting operation is performed on the last dimension in that
    !! data.
    real(kind=rk), intent(in) :: parent_data(:,:)

    !> Whether to ignore high modes that exceed the target maximal polynomial
    !! degree.
    !!
    !! This can be used as a simple lowpass filter that cuts off the highest
    !! modes in the parent elements prior to mapping to child elements.
    logical, intent(in) :: ignore

    !> Computed projection of the polynomial representation in the child
    !! elements.
    !!
    !! Again, the first index refers to the degrees of freedom, while the
    !! second index are the elements. There need to be twice as many elements
    !! as in the parent_data.
    !! Left childs are stored in iChild = (iParent*2 - 1), and the right
    !! childs in iParent*2.
    !!
    !! In the first index the shape of the data has to be in the form
    !! (outLen(1), outLen(2), ... , outLen(nDims)), the data is rotated
    !! in comparison to parent_data and the splitted direction has to be
    !! the first one in child_data (while it was the last in parent_data),
    !! and all other dimensions are shifted by one to the right.
    real(kind=rk), intent(out) :: child_data(:,:)
    ! -------------------------------------------------------------------- !
    integer :: iDir
    integer :: iParent, Lchild, Rchild
    integer :: parentMode, childMode
    integer :: maxrow
    integer :: maxcol
    integer :: indep
    integer :: nIndeps
    integer :: nParents
    integer :: parentpos, childpos
    ! -------------------------------------------------------------------- !

    nParents = size(parent_data,2)

    ! Use split_legendre to compute the two child_data elements for each
    ! parent_data element.
    ! We store the left childs in iChild = (iParent*2 - 1), and the right
    ! childs in iParent*2.

    child_data = 0.0_rk

    ! The number of independent modes (in normal directions) is given
    ! by the product of the length in all directions, except the last one.
    nIndeps = 1
    do iDir=1,nDims-1
      nIndeps = nIndeps*inLen(iDir)
    end do

    if (ignore) then
      maxcol = min(outLen(1), inLen(nDims))
    else
      maxcol = inLen(nDims)
    end if
    oldmodes: do parentMode=1,maxcol
      ! Maximal number modes to compute, as this is a triangular matrix
      ! it is limited by the diagonal (parentMode). However, it may be
      ! that the target polynomial space in the output is smaller, in this
      ! case we cap the computations and no more than outLen(1) entries
      ! are to be computed.
      maxrow = min(parentMode, outLen(1))

      elemloop: do iParent=1,nParents
        Rchild = iParent*2
        Lchild = Rchild - 1
        newmodes: do childMode=1,maxrow

          do indep=1,nIndeps
            parentpos = indep + nIndeps*(parentMode-1)
            childpos = childmode + (indep-1)*outLen(1)
            child_data(childpos, Lchild) = child_data(childpos, Lchild) &
              &                          + split_legendre( parentmode,  &
              &                                            childmode )  &
              &                            * parent_data(parentpos, iParent)
            child_data(childpos, Rchild) = child_data(childpos, Rchild) &
              &                          + split_legendre( childmode,   &
              &                                            parentmode ) &
              &                            * parent_data(parentpos, iParent)
          end do

        end do newmodes
      end do elemloop

    end do oldmodes

  end subroutine ply_split_element_singleD
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Split one-dimensional elements of degree parent_degree into two elements
  !! with polynomials of degree child_degree.
  subroutine ply_split_element_1D( parent_degree, child_degree, parent_data, &
    &                              child_data, ignore_highmodes              )
    ! -------------------------------------------------------------------- !
    !> Polynomial degree in the parent element.
    integer, intent(in) :: parent_degree

    !> Polynomial degree in the child elements.
    integer, intent(in) :: child_degree

    !> Polynomial data in the parent element. The first index describes the
    !! degrees of freedom. The second index refers to the elements to split.
    real(kind=rk), intent(in) :: parent_data(:,:)

    !> Polynomial data in the child elements. The first index describes the
    !! degrees of freedom. The second index refers to the elements, there
    !! needs to be four times as many elements than in the parent_data.
    !!
    !! Elements follow the ordering of the Z space filling curve.
    real(kind=rk), intent(out) :: child_data(:,:)

    !> Whether to ignore high modes from the parent element.
    !!
    !! This can be used as a simple lowpass filter by ignoring all higher
    !! modes from the parent element, that exceed the target polynomial
    !! degree. Thus, the polynomials are filtered before projection,
    !! instead of cutting them only of after refinement.
    !! Defaults to false (no filtering).
    logical, optional, intent(in) :: ignore_highmodes
    ! -------------------------------------------------------------------- !
    logical :: ignore
    integer :: pardofs
    integer :: childdofs
    ! -------------------------------------------------------------------- !

    ignore = .false.
    if (present(ignore_highmodes)) then
      ignore = ignore_highmodes
    end if

    pardofs = parent_degree + 1
    childdofs = child_degree + 1

    call ply_split_element_singleD( nDims       = 1,            &
      &                             inLen       = [pardofs],    &
      &                             outLen      = [childdofs],  &
      &                             ignore      = ignore,       &
      &                             parent_data = parent_data,  &
      &                             child_data  = child_data    )

  end subroutine ply_split_element_1D
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Split two-dimensional elements of degree parent_degree into four elements
  !! with polynomials of degree child_degree.
  subroutine ply_split_element_2D( parent_degree, child_degree, parent_data, &
    &                              child_data, ignore_highmodes              )
    ! -------------------------------------------------------------------- !
    !> Polynomial degree in the parent element.
    integer, intent(in) :: parent_degree

    !> Polynomial degree in the child elements.
    integer, intent(in) :: child_degree

    !> Polynomial data in the parent element. The first index describes the
    !! degrees of freedom. The second index refers to the elements to split.
    real(kind=rk), intent(in) :: parent_data(:,:)

    !> Polynomial data in the child elements. The first index describes the
    !! degrees of freedom. The second index refers to the elements, there
    !! needs to be four times as many elements than in the parent_data.
    !!
    !! Elements follow the ordering of the Z space filling curve.
    real(kind=rk), intent(out) :: child_data(:,:)

    !> Whether to ignore high modes from the parent element.
    !!
    !! This can be used as a simple lowpass filter by ignoring all higher
    !! modes from the parent element, that exceed the target polynomial
    !! degree. Thus, the polynomials are filtered before projection,
    !! instead of cutting them only of after refinement.
    !! Defaults to false (no filtering).
    logical, optional, intent(in) :: ignore_highmodes
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: ysplit(:,:)
    logical :: ignore
    integer :: pardofs
    integer :: childdofs
    ! -------------------------------------------------------------------- !

    ignore = .false.
    if (present(ignore_highmodes)) then
      ignore = ignore_highmodes
    end if

    pardofs = parent_degree + 1
    childdofs = child_degree + 1

    allocate(ysplit(childdofs*pardofs, 2))

    call ply_split_element_singleD( nDims       = 2,                     &
      &                             inLen       = [pardofs, pardofs],    &
      &                             outLen      = [childdofs, pardofs],  &
      &                             ignore      = ignore,                &
      &                             parent_data = parent_data,           &
      &                             child_data  = ysplit                 )

    call ply_split_element_singleD( nDims       = 2,                       &
      &                             inLen       = [childdofs, pardofs],    &
      &                             outLen      = [childdofs, childdofs],  &
      &                             ignore      = ignore,                &
      &                             parent_data = ysplit,                  &
      &                             child_data  = child_data               )

    deallocate(ysplit)

  end subroutine ply_split_element_2D
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Split three-dimensional elements of degree parent_degree into eight
  !! elements with polynomials of degree child_degree.
  subroutine ply_split_element_3D( parent_degree, child_degree, parent_data, &
    &                              child_data, ignore_highmodes              )
    ! -------------------------------------------------------------------- !
    !> Polynomial degree in the parent element.
    integer, intent(in) :: parent_degree

    !> Polynomial degree in the child elements.
    integer, intent(in) :: child_degree

    !> Polynomial data in the parent element. The first index describes the
    !! degrees of freedom. The second index refers to the elements to split.
    real(kind=rk), intent(in) :: parent_data(:,:)

    !> Polynomial data in the child elements. The first index describes the
    !! degrees of freedom. The second index refers to the elements, there
    !! needs to be four times as many elements than in the parent_data.
    !!
    !! Elements follow the ordering of the Z space filling curve.
    real(kind=rk), intent(out) :: child_data(:,:)

    !> Whether to ignore high modes from the parent element.
    !!
    !! This can be used as a simple lowpass filter by ignoring all higher
    !! modes from the parent element, that exceed the target polynomial
    !! degree. Thus, the polynomials are filtered before projection,
    !! instead of cutting them only of after refinement.
    !! Defaults to false (no filtering).
    logical, optional, intent(in) :: ignore_highmodes
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: ysplit(:,:)
    real(kind=rk), allocatable :: zsplit(:,:)
    logical :: ignore
    integer :: pardofs
    integer :: childdofs
    ! -------------------------------------------------------------------- !

    pardofs = parent_degree + 1
    childdofs = child_degree + 1

    ignore = .false.
    if (present(ignore_highmodes)) then
      ignore = ignore_highmodes
    end if


    allocate(zsplit(childdofs * pardofs**2, 2))
    allocate(ysplit(childdofs**2 * pardofs, 4))

    call ply_split_element_singleD( nDims       = 3,                     &
      &                             inLen       = [ pardofs, pardofs,    &
      &                                             pardofs           ], &
      &                             outLen      = [ childdofs, pardofs,  &
      &                                             pardofs           ], &
      &                             ignore      = ignore,                &
      &                             parent_data = parent_data,           &
      &                             child_data  = zsplit                 )

    call ply_split_element_singleD( nDims       = 3,                      &
      &                             inLen       = [ childdofs, pardofs,   &
      &                                             pardofs           ],  &
      &                             outLen      = [ childdofs, childdofs, &
      &                                             pardofs           ],  &
      &                             ignore      = ignore,                 &
      &                             parent_data = zsplit,                 &
      &                             child_data  = ysplit                  )

    call ply_split_element_singleD( nDims       = 3,                       &
      &                             inLen       = [ childdofs, childdofs,  &
      &                                             pardofs             ], &
      &                             outLen      = [ childdofs, childdofs,  &
      &                                             childdofs           ], &
      &                             ignore      = ignore,                  &
      &                             parent_data = ysplit,                  &
      &                             child_data  = child_data               )

    deallocate(ysplit)
    deallocate(zsplit)

  end subroutine ply_split_element_3D
  ! ======================================================================== !


  ! !!!!!!! !
  ! testing !
  ! !!!!!!! !

  ! To test the transformation, we check various mode combinations.
  ! For those, with the same number of modes in the children as in the
  ! parents, the resulting polynomials in the children should coincide.
  ! We check this by creating random parent polynomials, and then testing
  ! a random set of points after the split operation.
  !
  ! When modes are cutted off, we only check the first mode, to see, whether
  ! the resulting integral mean, is the same as in the parent polynomial.
  !
  ! To identify the children we use the following terminology to refer to
  ! the directions:
  ! west   = -x, east  = +x
  ! south  = -y, north = +y
  ! bottom = -z, top   = +z
  !
  ! The children are expected to have the following layout:
  !
  !                      |   west  |   east  |
  !
  !                      +---------+---------+
  !                      |         |         |
  !                north |    7    |    8    |
  !                      |         |         |
  !    top layer:        +---------+---------+
  !                      |         |         |
  !                south |    5    |    6    |
  !                      |         |         |
  !                      +---------+---------+
  !
  !
  !                      +---------+---------+
  !                      |         |         |
  !                north |    3    |    4    |
  !                      |         |         |
  ! bottom layer:        +---------+---------+
  !                      |         |         |
  !                south |    1    |    2    |
  !                      |         |         |
  !                      +---------+---------+

  ! ------------------------------------------------------------------------ !
  !> Testing the 1D splitting.
  !!
  !! We test all combinations, even though a higher number of modes in the
  !! children is probably not that relevant, it should still be possible.
  subroutine ply_split_element_1D_test(nModes, success)
    ! -------------------------------------------------------------------- !
    !> Number of modes in the (1D) polynomials to use in the check.
    integer, intent(in) :: nModes

    !> Indication whether the tests were completed successfully.
    logical, intent(out) :: success
    ! -------------------------------------------------------------------- !
    integer :: parentModes, childmodes
    integer :: iPoint
    integer :: iElem
    real(kind=rk) :: xi(nModes)
    real(kind=rk) :: x_left(nModes)
    real(kind=rk) :: x_right(nModes)
    real(kind=rk) :: legchild(nModes, nModes)
    real(kind=rk) :: legleft(nModes, nModes)
    real(kind=rk) :: legright(nModes, nModes)
    real(kind=rk) :: rootval(nModes, 2)
    real(kind=rk) :: childval
    real(kind=rk), allocatable :: rootelem(:,:)
    real(kind=rk), allocatable :: childelem(:,:)
    real(kind=rk) :: tolerance
    ! -------------------------------------------------------------------- !

    call ply_split_element_init(nModes)

    tolerance = 8*epsilon(1.0_rk)*nmodes**2
    success = .true.

    ! Some random points to check the resulting child polynomials.
    call random_number(xi)

    legchild = ply_legendre_1d(xi, nModes-1)

    ! The corresponding positions in the left and right half of the root
    ! element.
    x_right = 0.5_rk*xi + 0.5_rk
    x_left  = 0.5_rk*xi - 0.5_rk

    legleft = ply_legendre_1d(x_left, nModes-1)
    legright = ply_legendre_1d(x_right, nModes-1)

    do parentmodes=1,nModes
      allocate(rootelem(parentModes,1))
      call random_number(rootelem)
      do iPoint=1,nModes
        rootval(iPoint,1) = sum( rootelem(:,1)                  &
          &                    * legleft(:parentModes,iPoint) )
        rootval(iPoint,2) = sum( rootelem(:,1)                   &
          &                     * legright(:parentModes,iPoint) )
      end do
      do childmodes=1,parentModes-1
        allocate(childelem(childmodes,2))
        call ply_split_element_singleD( nDims       = 1,             &
          &                             inLen       = [parentModes], &
          &                             outLen      = [childModes],  &
          &                             ignore      = .false.,       &
          &                             parent_data = rootelem,      &
          &                             child_data  = childelem      )
        success = success                                          &
          &       .and. ( 0.5_rk*(childelem(1,1) + childelem(1,2)) &
          &               - rootelem(1,1) < tolerance              )
        deallocate(childelem)
      end do
      do childmodes=parentmodes,nModes
        allocate(childelem(childmodes,2))
        call ply_split_element_singleD( nDims       = 1,             &
          &                             inLen       = [parentModes], &
          &                             outLen      = [childModes],  &
          &                             ignore      = .false.,       &
          &                             parent_data = rootelem,      &
          &                             child_data  = childelem      )
        do iElem=1,2
          do iPoint=1,nModes
            childval = sum( childelem(:,iElem)             &
                            * legchild(:childmodes,iPoint) )
            success = success &
              &       .and. ( abs(rootval(iPoint,iElem) - childval) &
              &               < tolerance                           )
          end do
        end do
        deallocate(childelem)
      end do
      deallocate(rootelem)
    end do

  end subroutine ply_split_element_1D_test
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Testing the 2D splitting.
  !!
  !! In two dimensions we only check the downsized polynomial splitting
  !! (child_degree <= parent_degree), upsized splitting is checked for 1D
  !! operations already.
  !! For child_degree == parent_degree the resulting polynomials are probed
  !! at a set of random points to ensure the polynomials coincide with the
  !! parent polynomial.
  !! For those, where modes are cut off, we check the integral mean to be
  !! maintained.
  subroutine ply_split_element_2D_test(nModes, success)
    ! -------------------------------------------------------------------- !
    !> Number of modes in the (1D) polynomials to use in the check.
    integer, intent(in) :: nModes

    !> Indication whether the tests were completed successfully.
    logical, intent(out) :: success
    ! -------------------------------------------------------------------- !
    integer :: parentModes, childmodes
    integer :: iPoint
    integer :: iElem
    integer :: iMode
    real(kind=rk) :: xi(nModes,2)
    real(kind=rk) :: x_southwest(nModes,2)
    real(kind=rk) :: x_southeast(nModes,2)
    real(kind=rk) :: x_northwest(nModes,2)
    real(kind=rk) :: x_northeast(nModes,2)
    real(kind=rk) :: legchild(nModes, nModes, 2)
    real(kind=rk) :: legsouthwest(nModes, nModes, 2)
    real(kind=rk) :: legsoutheast(nModes, nModes, 2)
    real(kind=rk) :: legnorthwest(nModes, nModes, 2)
    real(kind=rk) :: legnortheast(nModes, nModes, 2)
    real(kind=rk) :: rootvaly(nModes, 4)
    real(kind=rk) :: rootval(nModes, 4)
    real(kind=rk) :: childval
    real(kind=rk) :: childvaly(nModes)
    real(kind=rk), allocatable :: rootelem(:,:)
    real(kind=rk), allocatable :: childelem(:,:)
    real(kind=rk) :: tolerance
    ! -------------------------------------------------------------------- !

    call ply_split_element_init(nModes)

    tolerance = 8*epsilon(1.0_rk)*nmodes**2
    success = .true.

    ! Some random points to check the resulting child polynomials.
    call random_number(xi)

    legchild(:,:,1) = ply_legendre_1d(xi(:,1), nModes-1)
    legchild(:,:,2) = ply_legendre_1d(xi(:,2), nModes-1)

    ! The corresponding positions in the left and right half of the root
    ! element.
    x_southwest = 0.5_rk*xi - 0.5_rk
    x_northeast = 0.5_rk*xi + 0.5_rk

    x_southeast(:,1) = 0.5_rk*xi(:,1) + 0.5_rk
    x_southeast(:,2) = 0.5_rk*xi(:,2) - 0.5_rk

    x_northwest(:,1) = 0.5_rk*xi(:,1) - 0.5_rk
    x_northwest(:,2) = 0.5_rk*xi(:,2) + 0.5_rk

    legsouthwest(:,:,1) = ply_legendre_1d(x_southwest(:,1), nModes-1)
    legsouthwest(:,:,2) = ply_legendre_1d(x_southwest(:,2), nModes-1)

    legsoutheast(:,:,1) = ply_legendre_1d(x_southeast(:,1), nModes-1)
    legsoutheast(:,:,2) = ply_legendre_1d(x_southeast(:,2), nModes-1)

    legnorthwest(:,:,1) = ply_legendre_1d(x_northwest(:,1), nModes-1)
    legnorthwest(:,:,2) = ply_legendre_1d(x_northwest(:,2), nModes-1)

    legnortheast(:,:,1) = ply_legendre_1d(x_northeast(:,1), nModes-1)
    legnortheast(:,:,2) = ply_legendre_1d(x_northeast(:,2), nModes-1)

    do parentmodes=1,nModes
      allocate(rootelem(parentModes**2,1))
      call random_number(rootelem)
      do iPoint=1,nModes
        do iMode=1,parentmodes
          rootvaly(iMode,1) = sum( rootelem((iMode-1)*parentmodes+1    &
            &                               :iMode*parentmodes,1)      &
            &                      * legsouthwest(:parentModes,iPoint, &
            &                                     1)                   )
          rootvaly(iMode,2) = sum( rootelem((iMode-1)*parentmodes+1    &
            &                               :iMode*parentmodes,1)      &
            &                      * legsoutheast(:parentModes,iPoint, &
            &                                     1)                   )
          rootvaly(iMode,3) = sum( rootelem((iMode-1)*parentmodes+1    &
            &                               :iMode*parentmodes,1)      &
            &                      * legnorthwest(:parentModes,iPoint, &
            &                                     1)                   )
          rootvaly(iMode,4) = sum( rootelem((iMode-1)*parentmodes+1    &
            &                               :iMode*parentmodes,1)      &
            &                      * legnortheast(:parentModes,iPoint, &
            &                                     1)                   )
        end do
        rootval(iPoint,1) = sum( rootvaly(:parentmodes,1)              &
          &                      * legsouthwest(:parentModes,iPoint,2) )
        rootval(iPoint,2) = sum( rootvaly(:parentmodes,2)              &
          &                      * legsoutheast(:parentModes,iPoint,2) )
        rootval(iPoint,3) = sum( rootvaly(:parentmodes,3)              &
          &                      * legnorthwest(:parentModes,iPoint,2) )
        rootval(iPoint,4) = sum( rootvaly(:parentmodes,4)              &
          &                      * legnortheast(:parentModes,iPoint,2) )
      end do
      do childmodes=1,parentModes-1
        allocate(childelem(childmodes**2,4))
        call ply_split_element_2D( parent_degree = parentModes-1, &
          &                        child_degree  = childModes-1,  &
          &                        parent_data   = rootelem,      &
          &                        child_data    = childelem      )
        success = success                             &
          &       .and. ( 0.25_rk*sum(childelem(1,:)) &
          &               - rootelem(1,1) < tolerance )
        deallocate(childelem)
      end do
      childmodes = parentmodes
      allocate(childelem(childmodes**2,4))
      call ply_split_element_2D( parent_degree = parentModes-1, &
        &                        child_degree  = childModes-1,  &
        &                        parent_data   = rootelem,      &
        &                        child_data    = childelem      )
      do iElem=1,4
        do iPoint=1,nModes
          do iMode=1,childmodes
            childvaly(iMode) = sum( childelem((iMode-1)*childmodes+1   &
              &                               :iMode*childmodes,iElem) &
              &                     * legchild(:childmodes,iPoint,1)   )
          end do
          childval = sum( childvaly(:childmodes)           &
            &             * legchild(:childmodes,iPoint,2) )
          success = success &
            &       .and. ( abs(rootval(iPoint,iElem) - childval) &
            &               < tolerance                           )
        end do
      end do
      deallocate(childelem)
      deallocate(rootelem)
    end do

  end subroutine ply_split_element_2D_test
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Testing the 3D splitting.
  !!
  !! In three dimensions we only check the splitting to polynomials of same
  !! polynomial degree.
  !! Downsizing is checked in 2D, and upsizing is only checked for 1D.
  !! The resulting polynomials are probed at a set of random points to ensure
  !  the polynomials coincide with the parent polynomial.
  subroutine ply_split_element_3D_test(nModes, success)
    ! -------------------------------------------------------------------- !
    !> Number of modes in the (1D) polynomials to use in the check.
    integer, intent(in) :: nModes

    !> Indication whether the tests were completed successfully.
    logical, intent(out) :: success
    ! -------------------------------------------------------------------- !
    integer :: parentModes, childmodes
    integer :: iPoint
    integer :: iElem
    integer :: iMode, jMode
    integer :: ij
    integer :: iX, iY, iZ
    real(kind=rk) :: xi(nModes,3)
    real(kind=rk) :: x(nModes,3,8)
    real(kind=rk) :: legchild(nModes, nModes, 3)
    real(kind=rk) :: legparent(nModes, nModes, 3, 8)
    real(kind=rk) :: rootvalz(nModes, 8)
    real(kind=rk) :: rootvaly(nModes, 8)
    real(kind=rk) :: rootval(nModes, 8)
    real(kind=rk) :: childval
    real(kind=rk) :: childvaly(nModes)
    real(kind=rk) :: childvalz(nModes)
    real(kind=rk), allocatable :: rootelem(:,:)
    real(kind=rk), allocatable :: childelem(:,:)
    real(kind=rk) :: tolerance
    ! -------------------------------------------------------------------- !

    call ply_split_element_init(nModes)

    tolerance = 8*epsilon(1.0_rk)*nmodes**3
    success = .true.

    ! Some random points to check the resulting child polynomials.
    call random_number(xi)

    legchild(:,:,1) = ply_legendre_1d(xi(:,1), nModes-1)
    legchild(:,:,2) = ply_legendre_1d(xi(:,2), nModes-1)
    legchild(:,:,3) = ply_legendre_1d(xi(:,3), nModes-1)

    ! The corresponding positions in the left and right half of the root
    ! element.
    do iZ=0,1
      do iY=0,1
        do iX=0,1
          iElem = 1 + iX + (iY + iZ*2)*2
          x(:, 1, iElem) = 0.5_rk*xi(:,1) + 0.5_rk*(iX*2 - 1)
          x(:, 2, iElem) = 0.5_rk*xi(:,2) + 0.5_rk*(iY*2 - 1)
          x(:, 3, iElem) = 0.5_rk*xi(:,3) + 0.5_rk*(iZ*2 - 1)
        end do
      end do
    end do

    do iElem=1,8
      legparent(:,:,1,iElem) = ply_legendre_1d(x(:,1,iElem), nModes-1)
      legparent(:,:,2,iElem) = ply_legendre_1d(x(:,2,iElem), nModes-1)
      legparent(:,:,3,iElem) = ply_legendre_1d(x(:,3,iElem), nModes-1)
    end do

    do parentmodes=1,nModes
      allocate(rootelem(parentModes**3,1))
      call random_number(rootelem)
      do iPoint=1,nModes
        ! Evaluation in Y direction
        do jMode=1,parentmodes
          do iElem=1,8
            ! For each Y mode evaluate the polynomial in X (in each of the
            ! children elements)
            do iMode=1,parentmodes
              ij = (jMode-1)*parentmodes + iMode
              rootvaly(iMode,iElem) = sum( rootelem((ij-1)*parentmodes+1      &
                &                                   :ij*parentmodes,1)        &
                &                          * legparent(:parentModes,iPoint,1, &
                &                                      iElem)                 )
            end do
          end do
          ! Evaluate the current Y-Mode to get the 1D polynomial in Z at the
          ! xy coordinates of iPoint.
          do iElem=1,8
            rootvalz(jMode,iElem) = sum( rootvaly(:parentmodes,iElem)       &
              &                          * legparent(:parentModes,iPoint,2, &
              &                                      iElem)                 )
          end do
        end do
        ! Finally evaluate the Z polynomial at the position of iPoint.
        do iElem=1,8
          rootval(iPoint,iElem) = sum( rootvalz(:parentmodes,iElem)       &
            &                          * legparent(:parentModes,iPoint,3, &
            &                                      iElem)                 )
        end do
      end do
      childmodes = parentmodes
      allocate(childelem(childmodes**3,8))
      call ply_split_element_3D( parent_degree = parentModes-1, &
        &                        child_degree  = childModes-1,  &
        &                        parent_data   = rootelem,      &
        &                        child_data    = childelem      )
      do iElem=1,8
        do iPoint=1,nModes
          do jMode=1,childmodes
            do iMode=1,childmodes
              ij = (jMode-1)*childmodes + iMode
              childvaly(iMode) = sum( childelem((ij-1)*childmodes+1   &
                &                               :ij*childmodes,iElem) &
                &                     * legchild(:childmodes,iPoint,1)   )
            end do
            childvalz(jMode) = sum( childvaly(:childmodes)           &
              &                     * legchild(:childmodes,iPoint,2) )
          end do
          childval = sum( childvalz(:childmodes)           &
            &             * legchild(:childmodes,iPoint,3) )
          success = success &
            &       .and. ( abs(rootval(iPoint,iElem) - childval) &
            &               < tolerance                           )
        end do
      end do
      deallocate(childelem)
      deallocate(rootelem)
    end do

  end subroutine ply_split_element_3D_test
  ! ======================================================================== !



  ! ------------------------------------------------------------------------ !
  !> Testing routine for the functions of this module.
  subroutine ply_split_element_test(success)
    ! -------------------------------------------------------------------- !
    !> Indication whether the tests were completed successfully.
    logical, intent(out) :: success
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    call ply_split_element_init(80)

    ! The split_legendre matrix generation is already checked by the
    ! ply_split_legendre_test routine.

    call ply_split_element_1D_test(nModes = 30, success = success)

    if (.not. success) then
      write(*,*) 'Check for 1D splitting FAILED!'
      RETURN
    end if

    call ply_split_element_2D_test(nModes = 20, success = success)

    if (.not. success) then
      write(*,*) 'Check for 2D splitting FAILED!'
      RETURN
    end if

    call ply_split_element_3D_test(nModes = 10, success = success)

    if (.not. success) then
      write(*,*) 'Check for 3D splitting FAILED!'
      RETURN
    end if

  end subroutine ply_split_element_test
  ! ======================================================================== !

end module ply_split_element_module
