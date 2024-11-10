! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> This module provides methods to filter polynomial representation in
!! elements based on their shape.
!!
!! The main goal of this filtering is to smooth out Gibbs oscillations while
!! maintaining the strong gradients at discontinuities.
module ply_filter_element_module
  use env_module, only: rk, labelLen

  use aotus_module, only: flu_State, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_err_module, only: aoterr_fatal

  use tem_aux_module, only: tem_abort
  use tem_logging_module, only: logunit
  use tem_tools_module, only: upper_to_lower

  implicit none

  private

  public :: ply_filter_element
  public :: ply_filter_element_oddfract
  public :: ply_filter_element_type
  public :: ply_filter_element_load

  integer, parameter :: filter_strat_none = 0
  integer, parameter :: filter_strat_oddfract = 1

  !> Paramaters describing the filtering to apply to elemental polynomial data.
  type ply_filter_element_type
    !> Filter strategy to use.
    integer :: strategy = filter_strat_none

    !> Maximal order for exponential spectral filtering to use
    !! where little filtering is to be done.
    integer :: max_order

    !> Minimal order for exponential spectral filtering to use.
    integer :: min_order

    !> Exponent to use for the fraction.
    integer :: fract_exponent

    !> Function pointer for 1D filtering
    procedure(ply_filter_element), pointer :: filter1D => NULL()

    !> Function pointer for 2D filtering
    procedure(ply_filter_element), pointer :: filter2D => NULL()

    !> Function pointer for 3D filtering
    procedure(ply_filter_element), pointer :: filter3D => NULL()
  end type ply_filter_element_type


  abstract interface
    !> Filter the polynomial data in a given element.
    subroutine ply_filter_element( me, element_degree, element_data )
      ! -------------------------------------------------------------------- !
      import :: rk, ply_filter_element_type
      !> Parameters of the filter.
      class(ply_filter_element_type), intent(in) :: me

      !> Polynomial degree in the parent element.
      integer, intent(in) :: element_degree

      !> Polynomial data in element. The first index describes the
      !! degrees of freedom. The second index refers to the elements to filter.
      real(kind=rk), intent(inout) :: element_data(:,:)
      ! -------------------------------------------------------------------- !
    end subroutine ply_filter_element
  end interface


contains


  ! ------------------------------------------------------------------------ !
  !> Loading parameters for the filtering from the configuration script.
  !! This needs to be performed before any call of the actual transformation
  !! [[ply_split_element_1D]].
  !!
  !! The initialization will compute the transformation matrix for Legendre
  !! polynomials with at least nMaxModes. If the initialization was already
  !! called before with the same or larger nMaxModes, the matrix will not be
  !! changed. Thus, calling this routine will only increase the size of the
  !! module variable split_legendre, never decrease it.
  subroutine ply_filter_element_load(me, conf, parent)
    ! -------------------------------------------------------------------- !
    !> Data structure that holds the filter parameters.
    type(ply_filter_element_type), intent(out) :: me

    !> Lua script to get the filter parameters from.
    type(flu_state) :: conf

    !> Table handle to a possible parent, that contains the filter table
    !! to load.
    integer, optional, intent(in) :: parent
    ! -------------------------------------------------------------------- !
    integer :: thandle
    character(len=labelLen) :: stratname
    integer :: iError
    ! -------------------------------------------------------------------- !

    nullify(me%filter1D)
    nullify(me%filter2D)
    nullify(me%filter3D)

    call aot_table_open( L       = conf,            &
      &                  parent  = parent,          &
      &                  thandle = thandle,         &
      &                  key     = 'filter_element' )

    call aot_get_val( L       = conf,       &
      &               thandle = thandle,    &
      &               key     = 'strategy', &
      &               val     = stratname,  &
      &               ErrCode = iError,     &
      &               default = 'none'      )

    me%strategy = filter_strat_none
    stratname = upper_to_lower(trim(stratname))
    select case(trim(stratname))
    case ('oddfract')
      me%strategy = filter_strat_oddfract

    case ('none')
      me%strategy = filter_strat_none

    case default
      write(logunit(1),*) 'ERROR: Strategy ', trim(stratname), ' for filter' &
        &                 // ' element not known!'
      write(logunit(1),*) 'Available options are:'
      write(logunit(1),*) '* oddfract: filter based on the fraction of odd'
      write(logunit(1),*) '            modes in the spectral energy'
      write(logunit(1),*) '* none: deactivate element filtering'
      call tem_abort()
    end select

    select case(me%strategy)
    case(filter_strat_oddfract)
      ! Get parameters for the odd mode fraction filtering strategy.
      call aot_get_val( L       = conf,         &
        &               thandle = thandle,      &
        &               key     = 'max_order',  &
        &               val     = me%max_order, &
        &               ErrCode = iError,       &
        &               default = 10            )
      if ( btest(iError, aoterr_Fatal) ) then
        write(logunit(1),*) 'ERROR: max_order for filter_element needs to be' &
          &                 // ' an integer!'
        write(logunit(1),*) 'You provided max_order but not in a form that' &
          &                 // ' could be interpreted as a number.'
        write(logunit(1),*) 'Aborting the execution, please check your config!'
        call tem_abort()
      end if

      call aot_get_val( L       = conf,         &
        &               thandle = thandle,      &
        &               key     = 'min_order',  &
        &               val     = me%min_order, &
        &               ErrCode = iError,       &
        &               default = 2             )
      if ( btest(iError, aoterr_Fatal) ) then
        write(logunit(1),*) 'ERROR: min_order for filter_element needs to be' &
          &                 // ' an integer!'
        write(logunit(1),*) 'You provided min_order but not in a form that' &
          &                 // ' could be interpreted as a number.'
        write(logunit(1),*) 'Aborting the execution, please check your config!'
        call tem_abort()
      end if

      call aot_get_val( L       = conf,              &
        &               thandle = thandle,           &
        &               key     = 'fract_exponent',  &
        &               val     = me%fract_exponent, &
        &               ErrCode = iError,            &
        &               default = 2                  )
      if ( btest(iError, aoterr_Fatal) ) then
        write(logunit(1),*) 'ERROR: fract_exponent for filter_element needs' &
          &                 // ' to be an integer!'
        write(logunit(1),*) 'You provided fract_exponent but not in a form' &
          &                 // ' that could be interpreted as a number.'
        write(logunit(1),*) 'Aborting the execution, please check your config!'
        call tem_abort()
      end if

      write(logunit(3),*) 'Using element filtering with the odd fraction'
      write(logunit(3),*) 'strategy of data before each refinement.'
      write(logunit(3),*) 'Following parameters are used:'
      write(logunit(3),*) '* min_order=', me%min_order
      write(logunit(3),*) '* max_order=', me%max_order
      write(logunit(3),*) '* fract_exponent=', me%fract_exponent

      me%filter1D => ply_filter_oddfract_1D
      me%filter2D => ply_filter_oddfract_2D
      me%filter3D => ply_filter_oddfract_3D
    end select

    call aot_table_close( L       = conf,   &
      &                   thandle = thandle )

  end subroutine ply_filter_element_load
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Filter a polynomial representation in elements in one dimension according
  !! to its odd mode fraction.
  !!
  !! Odd and even modes are weighed with their polynomial degree (so it's a
  !! little like actually using the derivative), squared and summed
  !! respectively.
  !! We then compute a spectral damping order with fraction of the odd mode
  !! energy in the total modal energy. The larger the fraction, the larger
  !! damping order and the weaker the filtering.
  !! `uscale` is given by `(me%max_order-me%min_order+1)`.
  !! `odd_fraction` is given by `odd / (odd+even)` and the damping order is
  !! then computed by `me%min_order + floor(uscale * odd_fraction**me%fract_ex)`
  !! The `fract_exponent` provides a mechanism to take larger odd fractions
  !! stronger into account than smaller ones.
  !!
  !! As we need to perform this operation in all dimensions, it would be good
  !! to shift the indices around. When doing this, we can stick to the same
  !! implementation for all directions, without the need to put any logic in
  !! here to decide on the current direction.
  !! In 3D we would end up with this chain:
  !! (x,y,z) -> filter_element for Z -> (z,x,y)
  !!         -> filter_element for Y -> (y,z,x)
  !!         -> filter_element for X -> (x,y,z)
  !! Thus, the logic is that we perform the filter on the last dimension, and
  !! cycle the indices in the output.
  !!
  !! We can generalize this to arbitrary dimensions.
  !! In 2D it would look like this:
  !! (x,y) -> filter_element for Y -> (y,x)
  !!       -> filter_element for X -> (x,y)
  !! And in 1D, we just need to perform one transformation:
  !! (x) -> filter_element for X -> (x)
  !!
  !! We need: nDofs in the direction where the transformation is to be done
  !!          and the nDofs for all normal directions.
  subroutine ply_filter_element_oddfract( me, nDims, inLen, element_data, &
    &                                     filtered_data                   )
    ! -------------------------------------------------------------------- !
    !> Filter parameters.
    class(ply_filter_element_type), intent(in) :: me

    !> Number of dimensions of the polynomial data.
    integer, intent(in) :: nDims

    !> Number degrees of freedom for each direction in element_data.
    !!
    !! The first index of element_data needs to have a length equal to the
    !! product of all inLen components.
    !! The splitting operation will be done in the last dimension.
    integer, intent(in) :: inLen(nDims)

    !> Polynomial representation in the elements.
    !!
    !! The first index are the degrees of freedom in elements, the second index
    !! are the elements.
    !! In the first index the shape of data has to be in the form
    !! (inLen(1), inLen(2), ... , inLen(nDims)).
    !! The filtering operation is performed on the last dimension in that
    !! data.
    real(kind=rk), intent(in) :: element_data(:,:)

    !> The filtered polynomial modes.
    !!
    !! The ordering is rotated, such, that the filtered dimension becomes the
    !! first one, and all others are shifted by one to the right.
    !! Thus, the new data has the layout (inLen(nDims), inLen(1), inLen(2), ...)
    real(kind=rk), intent(out) :: filtered_data(:,:)
    ! -------------------------------------------------------------------- !
    integer :: iDir
    integer :: iParent
    integer :: parentMode
    integer :: maxcol
    integer :: indep
    integer :: nIndeps
    integer :: nParents
    integer :: parentpos, filterpos
    integer :: upper_scale
    real(kind=rk) :: dof_fract
    real(kind=rk) :: damping
    real(kind=rk), allocatable :: even(:,:)
    real(kind=rk), allocatable :: odd(:,:)
    real(kind=rk) :: spenergy
    integer, allocatable :: damp_ord(:,:)
    ! -------------------------------------------------------------------- !

    nParents = size(element_data,2)

    ! The number of independent modes (in normal directions) is given
    ! by the product of the length in all directions, except the last one.
    nIndeps = 1
    do iDir=1,nDims-1
      nIndeps = nIndeps*inLen(iDir)
    end do

    maxcol = inLen(nDims)

    allocate(damp_ord(nIndeps, nParents))
    allocate(even(nIndeps, nParents))
    allocate(odd(nIndeps, nParents))
    even = 0.0_rk
    odd = 0.0_rk

    parmodes: do parentMode=2,maxcol-1,2
      do iParent=1,nParents
        do indep=1,nIndeps
          parentpos = indep + nIndeps*(parentMode-1)
          odd(indep, iParent) = odd(indep, iParent)                  &
            &                 + ( (parentMode-1)                     &
            &                     * element_data(parentpos, iParent) &
            &                   )**2
          even(indep, iParent) = even(indep, iParent)                         &
            &                  + ( (parentMode-1)                             &
            &                      * element_data(parentpos+nIndeps, iParent) &
            &                    )**2
        end do
      end do
    end do parmodes
    if ((mod(maxcol,2) == 0)) then
      do iParent=1,nParents
        do indep=1,nIndeps
          parentpos = indep + nIndeps*(maxcol-1)
          even(indep, iParent) = even(indep, iParent)                &
            &                 + ( (parentMode-1)                     &
            &                     * element_data(parentpos, iParent) &
            &                   )**2
        end do
      end do
    end if

    upper_scale = me%max_order - me%min_order + 1
    do iParent=1,nParents
      do indep=1,nIndeps
        spenergy = even(indep, iParent)+odd(indep,iParent)
        if (spenergy > epsilon(damping)**2) then
          damp_ord(indep,iParent) = me%min_order                         &
            &                     +floor( upper_scale                    &
            &                             * (odd(indep,iParent)          &
            &                             / spenergy)**me%fract_exponent )
        else
          damp_ord(indep,iParent) = me%max_order
        end if
      end do
    end do

    oldmodes: do parentMode=1,maxcol
      dof_fract = real(parentMode-1,kind=rk)/real(max(maxcol-1,1),kind=rk)

      elemloop: do iParent=1,nParents
        do indep=1,nIndeps
          parentpos = indep + nIndeps*(parentMode-1)
          filterpos = parentMode + maxcol*(indep-1)
          spenergy = even(indep, iParent)+odd(indep,iParent)
          if (spenergy > epsilon(damping)**2) then
            damp_ord = me%min_order                         &
              &      +floor( (me%max_order+1)               &
              &              * (odd(indep,iParent)          &
              &              / spenergy)**me%fract_exponent )
          else
            damp_ord = me%max_order
          end if
          damping = exp(-36 * (dof_fract**damp_ord(indep,iParent)))
          filtered_data(filterpos, iParent) = damping &
            &                               * element_data(parentpos, iParent)
        end do

      end do elemloop

    end do oldmodes

    deallocate(damp_ord)
    deallocate(even)
    deallocate(odd)

  end subroutine ply_filter_element_oddfract
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Filter one-dimensional elements of degree element_degree.
  subroutine ply_filter_oddfract_1D( me, element_degree, element_data )
    ! -------------------------------------------------------------------- !
    !> Filter parameters.
    class(ply_filter_element_type), intent(in) :: me

    !> Polynomial degree in the parent element.
    integer, intent(in) :: element_degree

    !> Polynomial data in the parent element. The first index describes the
    !! degrees of freedom. The second index refers to the elements to split.
    real(kind=rk), intent(inout) :: element_data(:,:)
    ! -------------------------------------------------------------------- !
    integer :: pardofs
    real(kind=rk), allocatable :: filtered_data(:,:)
    ! -------------------------------------------------------------------- !

    pardofs = element_degree + 1
    allocate( filtered_data(pardofs, size(element_data,2)) )
    call ply_filter_element_oddfract( me            = me,           &
      &                               nDims         = 1,            &
      &                               inLen         = [pardofs],    &
      &                               element_data  = element_data, &
      &                               filtered_data = filtered_data )
    element_data = filtered_data
    deallocate(filtered_data)

  end subroutine ply_filter_oddfract_1D
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Filter two-dimensional elements of degree element_degree.
  subroutine ply_filter_oddfract_2D( me, element_degree, element_data )
    ! -------------------------------------------------------------------- !
    !> Filter parameters.
    class(ply_filter_element_type), intent(in) :: me

    !> Polynomial degree in the parent element.
    integer, intent(in) :: element_degree

    !> Polynomial data in the parent element. The first index describes the
    !! degrees of freedom. The second index refers to the elements to split.
    real(kind=rk), intent(inout) :: element_data(:,:)
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: ysplit(:,:)
    integer :: pardofs
    ! -------------------------------------------------------------------- !

    pardofs = element_degree + 1

    allocate(ysplit(pardofs*pardofs, size(element_data,2)))

    call ply_filter_element_oddfract( me            = me,                 &
      &                               nDims         = 2,                  &
      &                               inLen         = [pardofs, pardofs], &
      &                               element_data  = element_data,       &
      &                               filtered_data = ysplit              )
    call ply_filter_element_oddfract( me            = me,                 &
      &                               nDims         = 2,                  &
      &                               inLen         = [pardofs, pardofs], &
      &                               element_data  = ysplit,             &
      &                               filtered_data = element_data        )

    deallocate(ysplit)

  end subroutine ply_filter_oddfract_2D
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Filter three-dimensional elements of degree element_degree.
  subroutine ply_filter_oddfract_3D( me, element_degree, element_data )
    ! -------------------------------------------------------------------- !
    !> Filter parameters.
    class(ply_filter_element_type), intent(in) :: me

    !> Polynomial degree in the parent element.
    integer, intent(in) :: element_degree

    !> Polynomial data in the parent element. The first index describes the
    !! degrees of freedom. The second index refers to the elements to split.
    real(kind=rk), intent(inout) :: element_data(:,:)
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: ysplit(:,:)
    integer :: pardofs
    ! -------------------------------------------------------------------- !

    pardofs = element_degree + 1

    allocate(ysplit(pardofs*pardofs, size(element_data,2)))

    call ply_filter_element_oddfract( me            = me,                 &
      &                               nDims         = 3,                  &
      &                               inLen         = [ pardofs, pardofs, &
      &                                                 pardofs ],        &
      &                               element_data  = element_data,       &
      &                               filtered_data = ysplit              )

    call ply_filter_element_oddfract( me            = me,                 &
      &                               nDims         = 3,                  &
      &                               inLen         = [ pardofs, pardofs, &
      &                                                 pardofs ],        &
      &                               element_data  = ysplit,             &
      &                               filtered_data = element_data        )

    call ply_filter_element_oddfract( me            = me,                 &
      &                               nDims         = 3,                  &
      &                               inLen         = [ pardofs, pardofs, &
      &                                                 pardofs ],        &
      &                               element_data  = element_data,       &
      &                               filtered_data = ysplit              )
    element_data = ysplit


    deallocate(ysplit)

  end subroutine ply_filter_oddfract_3D
  ! ======================================================================== !


end module ply_filter_element_module
