! Copyright (c) 2014 Verena Krupp
! Copyright (c) 2014, 2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
!
! Parts of this file were written by Verena Krupp, Harald Klimach, Peter Vitt
! and Neda Ebrahimi Pour for University of Siegen.
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
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012, 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Daniel Harlacher <d.harlacher@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015-2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Simon Zimny and Manuel
! Hasert for German Research School for Simulation Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Kannan Masilamani,
! Daniel Harlacher, Kartik Jain, Verena Krupp, Jiaxing Qi, Peter Vitt,
! Daniel Fleischer, Tobias Girresser and Daniel Petró for University Siegen.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! This file contains the source code for growing and dynamic arrays.
! This is used for arrays of primitives (int, long_int, real, ...) as well as
! for arrays of derived datatypes (tem_variable_type,...).
!
! To use these macros include the following to your source file.
!
! Smart growing array (GA) for ?tstring?
! Growing Arrays:
!
! declaration
!
!
! implementation
!

! -----------------------------------------------------------------
! 2d Array, which can grow in second dimension only (GA2d)
! tname ... indicates type of dynamic array (long, int, real, ...)


!
!------------------------------------------------------------------------------
!
! dynamic Arrays:
!
! declaration
!
!
! implementation
!
!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
module ply_dynarray_project_module

  use env_module,             only: minLength, &
    &                               zeroLength
  use aotus_module,           only: flu_State
  use tem_logging_module,     only: logUnit

  use ply_prj_header_module


  implicit none

  private

  !> Projection definition.
  type ply_prj_init_type
    !> Polynomial basis type.
    !!
    !! 3D Monomials have the form x^i * y^j * z^k
    !! - Q_space: quadratic polynomial space (i,j,k) <= maxPolyDegree
    !! - P_space: polynomial space i+j+k <= maxPolyDegree
    integer :: basisType
    !> The maximal polynomial degree per spatial direction.
    integer :: maxPolyDegree
    !> projection header consits of general information like which kind
    !! of projection is used
    type(ply_prj_header_type) :: header
  end type ply_prj_init_type

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

  interface assignment(=)
    module procedure Copy_ply_prj_init
  end interface


! \brief smart dynamic array (da) for type(ply_prj_init_type)
!
! this datatype implements a dynamic array, which is capable of
! growing and adding of unique elements. it is available for
! various types, here we deal with $tstring$.
! sorted array contains the pointers to val array, instead of
! the actual values in val array. for example:
! val:     8, 6, 7, 5
! sorted:  4, 2, 3, 1
  !> dynamic array (da) type for type(ply_prj_init_type)
  type dyn_projectionarray_type
    integer :: nvals = 0
    integer :: containersize = 0
    type(ply_prj_init_type), allocatable :: val(:)
    integer, allocatable :: sorted(:) !< pointers, not values
  end type

  !> initialize the dynamic array
  interface init
    module procedure init_da_projection
  end interface

  !> destroy the dynamic array
  interface destroy
    module procedure destroy_da_projection
  end interface

  !> append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_da_projection
    module procedure append_da_vecprojection
  end interface

  !> truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_da_projection
  end interface

  !> empty the array, reset nvals to be 0
  interface empty
    module procedure empty_da_projection
  end interface


  !> fix the dynamic array, meaning:
  !! store the array in the sorted order and cut
  !! off the trailing empty entries
  interface sorttruncate
    module procedure sorttruncate_da_projection
  end interface

  !> increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_da_projection
  end interface

  !> return the position of a given value
  !! in the array val, which is what you usually want to know.
  !! it is the index of a given value
  interface positionofval
    module procedure posofval_projection
  end interface

  !> return the position of a given value
  !! in the list 'sorted'. this is mainly for internal usage.
  !! the sorted list is only a pointer list to the actual values
  !! thus, in order to get the index of a given value, you
  !! need to look up the entry in the sorted list.
  !! this is done by the positionofval routine
  interface sortedposofval
    module procedure sortposofval_projection
  end interface

  public :: init, append, dyn_projectionArray_type, truncate, empty
  public :: sortTruncate, positionOfVal
  public :: operator(==), operator(/=), operator(<), operator(<=)
  public :: operator(>), operator(>=)
  public :: ply_prj_init_define, ply_prj_init_type
  public :: ply_fill_dynProjectArray


contains

! ***************************************************************************** !
  !> initialization of a dynamic array
  !!
  !! before a dynamic array can be used, it has to be initialized
  !! with this routine. the initial length provided here, can
  !! avoid reallocations and memory copying, if approximated
  !! correctly enough. if none is specified, the provided container
  !! initially will be of size 0.
  subroutine init_da_projection(me, length)
    !-----------------------------------------------------------------
    type(dyn_projectionarray_type), intent(out) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container
    !-----------------------------------------------------------------

    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if

    ! deallocate ...
    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
    ! ... and reallocate
    allocate(me%val(me%containersize))
    allocate(me%sorted(me%containersize))
    me%nvals = 0

  end subroutine init_da_projection

  !> destruction of a dynamic array
  !!
  !! this subroutine takes care of a proper destruction of a
  !! dynamic array, it frees the allocated memory and resets
  !! the internal counts to 0.
  subroutine destroy_da_projection(me)
    type(dyn_projectionarray_type), intent(inout) :: me !< dynamic array to init

    me%containersize = 0
    me%nvals         = 0

    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
  end subroutine destroy_da_projection
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a value to the dynamic array
  !!
  !! with this subroutine, a given value can be added to the
  !! dynamic array. the actual position of this value in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_projection(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type) :: me   !< array to append the value to
    type(ply_prj_init_type), intent(in)     :: val  !< value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, if the value is found
    integer, intent(out), optional :: pos
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded
    !------------------------------------------------------------------------
    integer :: foundpos
    integer :: i
    !------------------------------------------------------------------------

    ! do a binary search on existing entries (returns closest entry next to
    ! it if not found).
    foundpos = sortedposofval(me, val, .true.)
    if( present( wasadded ) ) wasadded = .false.

    ! if it found the value, the position is smaller than nvals
    if (foundpos <= me%nvals) then

      ! the returned position might actually be the right entry already or
      ! not, check for this here.
      if ( me%val(me%sorted(foundpos)) == val ) then

        ! found the value in a list of unique values,
        ! nothing to do, just return its position.
        if( present( pos ) ) pos = me%sorted(foundpos)

      else

        ! need to append a new value!

        if (me%nvals == huge(me%nvals)) then
           write(*,*) "reached end of integer range for dynamic array!"
           write(*,*) "aborting!!"
           stop
        end if

        if( present( wasadded ) ) wasadded = .true.
        if (me%nvals == me%containersize) then

          ! container is full, need to expand it
          call expand(me = me, length = length)
        end if
        me%nvals = me%nvals + 1

        ! put the new value into the last position in the
        ! array.
        me%val(me%nvals) = val
        do while( foundpos < me%nvals )
          if(me%val(me%sorted(foundpos)) /= val) then
            exit
          end if
          ! in case of multiple entries with the same value
          ! move on to the first differing entry.
          foundpos = foundpos + 1
        end do
        ! shift the sorted list of indices, to create a
        ! whole for the value to be inserted, at position
        ! foundpos.
        do i=me%nvals-1,foundpos,-1
          me%sorted(i+1) = me%sorted(i)
        end do
        ! put the index of the new value into the
        ! sorted list at the now freed position.
        me%sorted(foundpos) = me%nvals
        if( present( pos ) ) pos = me%nvals

      end if

    else

      ! value to append is larger than all existing ones,
      ! just put it to the end of the list, this captures
      ! also the case of empty lists.
      ! in this case foundpos = me%nvals + 1 holds.
      if( present( wasadded ) ) wasadded = .true.
      if (foundpos > me%containersize) then
        ! expand the array, if its boundary is reached
        call expand(me = me, length = length)
      end if
      me%nvals = foundpos
      me%val(foundpos) = val
      me%sorted(foundpos) = foundpos
      if( present( pos ) ) pos = foundpos

    end if

  end subroutine append_da_projection
! ***************************************************************************** !


! ***************************************************************************** !
  !> appending a sorted list of values to the dynamic array
  !!
  !! with this subroutine, a given list of sorted values can be added to the
  !! dynamic array. the actual positions of these values in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. with the wasadded flag, it is indicated,\n
  !! wasadded = true,  if this entry had to be added,\n
  !! wasadded = false, if this was already found in the array.
  subroutine append_da_vecprojection(me, val, length, pos, wasadded )
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type) :: me   !< array to append the value to
    type(ply_prj_init_type), intent(in)     :: val(:)  !< values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> position in the array, the values are found at.
    integer, intent(out), optional :: pos(:)
    !> flag to indicate, if val was newly added
    logical, intent(out), optional :: wasadded(:)
    !------------------------------------------------------------------------
    type(ply_prj_init_type) :: lastval
    logical :: addedval(size(val))
    integer :: i
    integer :: veclen
    integer :: maxlen
    integer :: nappend
    integer :: rem_app
    integer :: curval, ival, iold, iadd
    integer, allocatable :: newsorted(:)
    !------------------------------------------------------------------------

    if (size(val) == 0) return

    veclen = size(val)
    maxlen = veclen + me%nvals

    allocate(newsorted(maxlen))

    addedval = .false.

    iold = 1
    iadd = 1

    nappend = 0
    curval = 0

    ! select the first entry before the loop unconditionally without checks
    ! for uniqueness (nothing to check against yet).
    if ( me%val(me%sorted(iold)) <= val(iadd) ) then
      curval = curval + 1
      newsorted(curval) = me%sorted(iold)
      lastval = me%val(me%sorted(iold))
      iold = iold + 1
    else
      curval = curval + 1
      nappend = nappend + 1
      newsorted(curval) = me%nvals + nappend
      lastval = val(iadd)
      if (present(pos))  pos(iadd) = newsorted(curval)
      addedval(iadd) = .true.
      iadd = iadd + 1
    end if

    do ival=2,maxlen

      if ( (iadd <= veclen) .and. (iold <= me%nvals) ) then

        if ( me%val(me%sorted(iold)) <= val(iadd) ) then

          ! the original list's values are appended to newsorted before
          ! the additional list is appended.
          curval = curval + 1
          newsorted(curval) = me%sorted(iold)
          lastval = me%val(me%sorted(iold))
          iold = iold + 1

        else

          ! only append the value to unique lists, if it is not yet in the list.
          ! (if it is already in the list, it has to be the previous (curval-1)
          !  entry.)
          if ( lastval < val(iadd) ) then
            nappend = nappend + 1
            curval = curval + 1
            newsorted(curval) = me%nvals + nappend
            lastval = val(iadd)
            addedval(iadd) = .true.
          end if
          if (present(pos)) pos(iadd) = newsorted(curval)
          iadd = iadd + 1

        end if

      else

        ! reached the end of one or both of the sorted lists.
        exit

      end if

    end do

    if (iold <= me%nvals) then
      ! still some values from the original list left.
      newsorted(curval+1:me%nvals+nappend) = me%sorted(iold:me%nvals)
    end if

    if (iadd <= veclen) then
      ! still some values from the list to append left.
      rem_app = iadd
      do i = rem_app,veclen
        if ( lastval < val(iadd) ) then
          nappend = nappend + 1
          curval = curval + 1
          newsorted(curval) = me%nvals + nappend
          lastval = val(iadd)
          addedval(iadd) = .true.
        end if
        if (present(pos)) pos(iadd) = newsorted(curval)
        iadd = iadd + 1
      end do
    end if

    if (me%nvals > huge(me%nvals)-nappend) then
       write(*,*) "reached end of integer range for dynamic array!"
       write(*,*) "aborting!!"
       stop
    end if

    if (me%nvals + nappend > me%containersize) then
      call expand( me        = me,      &
        &          increment = nappend, &
        &          length    = length   )
    end if
    me%sorted(:me%nvals+nappend) = newsorted(:me%nvals+nappend)
    curval = me%nvals
    do iadd=1,veclen
      if (addedval(iadd)) then
        curval = curval + 1
        me%val(curval) = val(iadd)
      end if
    end do
    me%nvals = me%nvals + nappend

    if( present( wasadded ) ) wasadded = addedval

  end subroutine append_da_vecprojection
! ***************************************************************************** !


! ***************************************************************************** !
  !> truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !!
  subroutine truncate_da_projection(me)
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    type(ply_prj_init_type), allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------

    if (me%nvals < me%containersize) then
      allocate(swpval(me%nvals))
      allocate(swpsort(me%nvals))

      swpval = me%val(:me%nvals)
      swpsort = me%sorted(:me%nvals)

      call move_alloc(swpval, me%val)
      call move_alloc(swpsort, me%sorted)

      me%containersize = me%nvals
    end if

  end subroutine truncate_da_projection
! ***************************************************************************** !


! ***************************************************************************** !
  !> empty all contents of the array without changing the size or status of any
  !! array
  !!
  subroutine empty_da_projection(me)
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    ! reset the number of entries
    me%nvals = 0
  end subroutine empty_da_projection
! ***************************************************************************** !


! ***************************************************************************** !
  !> fixing the dynamic array
  !!
  !! truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !! store the array in the sorted order according to the sorted( ) array
  !!
  subroutine sorttruncate_da_projection(me)
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type) :: me !< array to sorttruncate
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type) :: tarray !< temporary array
    integer :: ival
    integer :: dpos
    !------------------------------------------------------------------------
    ! allocate the temporary array
    call init( me = tarray, length = me%nvals )
    ! copy the entries in a sorted fashion into the temporary array
    do ival = 1, me%nvals
      call append( me = tarray, val = me%val( me%sorted( ival )), &
           &       pos = dpos)
    enddo
    call destroy( me = me )

    me = tarray
    call destroy( me = tarray )

  end subroutine sorttruncate_da_projection
! ***************************************************************************** !


! ***************************************************************************** !
  !> expanding the dynamic array
  !!
  !! this is a helping subroutine, which doubles the container
  !! of the given dynamic array. as the container might be
  !! initially 0-sized, a module variable minlength has been introduced, which
  !! is used here, to at least create a container of this size.
  subroutine expand_da_projection(me, increment, length)
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type) :: me !< array to resize
    integer, optional :: increment !< used for vector append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !------------------------------------------------------------------------
    type(ply_prj_init_type), allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------
    integer :: addvals, explen
    !------------------------------------------------------------------------

    addvals = 1
    if (present(increment)) addvals = increment

    if (addvals > 0) then

      ! if length is present, use that, otherwise double the size
      if( present( length ) ) then
        explen = length
      else
        ! set the global minimum length, if doubling would be smaller than that
        explen = max(me%containersize, minlength)
      end if

      ! check whether all elements will fit
      if( addvals > explen ) then
        explen = addvals
      end if

      ! check whether the new size will exceed the max container size.
      if( (huge(me%containersize) - explen) <= me%containersize ) then
        ! if so, expand to the maximum size
        me%containersize = huge(me%containersize)
      else
        ! if not, expand to the calculated size
        me%containersize = me%containersize + explen
      end if

      ! only need to do something, if there are actually values to append.
      if (me%nvals > 0) then

        allocate(swpval(me%containersize))
        swpval(1:me%nvals) = me%val(1:me%nvals)
        call move_alloc( swpval, me%val )

        allocate(swpsort(me%containersize))
        swpsort(1:me%nvals) = me%sorted(1:me%nvals)
        call move_alloc( swpsort, me%sorted )

      else ! me%nvals == 0

        if( allocated(me%val) ) &
          deallocate(me%val)
        allocate(me%val(me%containersize))
        if( allocated(me%sorted) ) &
          deallocate(me%sorted)
        allocate(me%sorted(me%containersize))

      end if

    end if

  end subroutine expand_da_projection
! ***************************************************************************** !


! ***************************************************************************** !
  !> return the sorted position of a value in the given dynamic array
  !!
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  function sortposofval_projection(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type), intent(in) :: me !< dynamic array
    type(ply_prj_init_type), intent(in) :: val !< value to look for
    !> flag to indicate, if the next entry in the list should be returned,
    !! if the searched one is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position of val in the sorted list, 0 if not found
    !------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    type(ply_prj_init_type) :: lb_val, ub_val
    type(ply_prj_init_type) :: mid_val
    !------------------------------------------------------------------------

    retnext = .false.
    if (present(nextifnotfound)) retnext = nextifnotfound

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    !> binary search on sorted list
    do while(ub >= lb)
      lb_val = me%val(me%sorted(lb))

      if (val < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me%val(me%sorted(ub))

      if (val > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      ! safe guard against integer limit overflow
      mid = lb + (ub-lb) / 2
      mid_val = me%val(me%sorted(mid))
      if (val == mid_val) then
        pos = mid
        exit
      end if
      if (val > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function sortposofval_projection
! ***************************************************************************** !


! ***************************************************************************** !
  !> the actual position of a given value in the dynamic array
  !!
  !! most likely this is what you need in codes, using this
  !! data structure, it first does the binary search on the sorted
  !! values with sortposofval_projection and then returns the looked
  !! up position in the original unsorted array, which corresponds
  !! to the position returned by the append routine.
  function posofval_projection(me, val, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_projectionarray_type), intent(in) :: me !< dynamic array
    type(ply_prj_init_type), intent(in) :: val !< value to search for
    !> flag to indicate, if the position of the next entry in the sorted
    !! list should be returned instead, if val is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position in the array of the searche value, 0 if not found
    !------------------------------------------------------------------------
    integer :: sortpos
    integer :: lb, ub
    !------------------------------------------------------------------------

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    sortpos = sortedposofval(me, val, nextifnotfound, lb, ub)

    ! if result (sorted pos)
    if ((sortpos <= me%nvals) .and. (sortpos > 0)) then
      pos = me%sorted(sortpos)
    else
      pos = sortpos
    end if
  end function posofval_projection
! ***************************************************************************** !


  ! ************************************************************************ !
  pure subroutine Copy_ply_prj_init(left,right)
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_prj_init_type), intent(out) :: left
    !> fpt to copy from
    type(ply_prj_init_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%header = right%header

    left%maxPolyDegree = right%maxPolyDegree
    left%basisType = right%basisType

  end subroutine copy_ply_prj_init
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Define a projection, without filling its body.
  pure subroutine ply_prj_init_define(me, header, maxPolyDegree, basisType)
    ! -------------------------------------------------------------------- !
    type(ply_prj_init_type), intent(inout) :: me
    type(ply_prj_header_type), intent(in) :: header
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: basisType
    ! -------------------------------------------------------------------- !

    me%header = header
    me%maxPolyDegree = maxPolyDegree
    me%basisType = basisType

  end subroutine ply_prj_init_define
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load settings to describe a projection method from a Lua table.
  subroutine ply_fill_dynProjectArray( proj_pos, dyn_projectionArray,         &
    &                                  basisType, maxPolyDegree, conf, parent )
    ! -------------------------------------------------------------------- !
    integer, intent(inout) :: proj_pos
    type(dyn_ProjectionArray_type), intent(inout) :: dyn_projectionArray
    type(flu_State), intent(in) :: conf
    integer, intent(in) :: basisType
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: parent
    ! -------------------------------------------------------------------- !
    type(ply_prj_header_type) :: header
    type(ply_prj_init_type) :: proj_init
    ! -------------------------------------------------------------------- !

    ! check if it is general projection (no parent is present) or individuall
    ! projection
    call ply_prj_header_load( me   = header,  &
      &                       conf = conf,    &
      &                       parent = parent )
    ! --> projection table exist, we build up the init_type and check if
    ! it is already in the DA
    ! define the init_type
    call  ply_prj_init_define( me            = proj_init,     &
      &                        header        = header,        &
      &                        maxPolyDegree = maxPolyDegree, &
      &                        basisType     = basisType      )

    ! Call to the unqiue list -> dynamic array for the projection method
    ! chack if the list needs to be appended...store the position
    write(logUnit(5),*) 'The projection_initialization type is'
    write(logUnit(5),*) ' kind=', proj_init%header%kind
    write(logUnit(5),*) ' degree=', proj_init%maxPolydegree
    write(logUnit(5),*) ' basisType=', proj_init%basisType
    write(logUnit(5),*) 'Now append the dynamic list for poly projection'
    call append( me = dyn_projectionArray,  &
      &          val = proj_init,           &
      &          pos = proj_pos             )

  end subroutine ply_fill_dynProjectArray
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for equality of two projections.
  !!
  !! Two projections are considered to be equal, if their kind, nodes_kind,
  !! maxPolyDegree and oversampling are equal.
  pure function isEqual(left, right) result(equality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is equal??
    logical :: equality
    ! -------------------------------------------------------------------- !

    equality = ( left%header == right%header )              &
      & .and. ( left%maxPolyDegree == right%maxPolyDegree ) &
      & .and. ( left%basisType == right%basisType )

  end function isEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for unequality of two projections.
  !!
  !! Two projections are considered to be unequal, if their kind, nodes_kind,
  !! maxpolydegree or factor are not equal.
  pure function isUnequal(left, right) result(unequality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    ! -------------------------------------------------------------------- !

    unequality = ( left%header /= right%header )           &
      & .or. ( left%maxPolyDegree /= right%maxPolyDegree ) &
      & .or. ( left%basisType /= right%basisType )

  end function isUnequal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a < comparison of two projections.
  !!
  !! Sorting of projections is given by maxPolyDegree, kind, nodes_kind and
  !! last by factor.
  pure function isSmaller(left, right) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !

    small = .false.
    if (left%header < right%header) then
      small = .true.
    else
      if (left%header == right%header) then
        if (left%maxPolyDegree < right%maxPolyDegree) then
          small = .true.
        else
          if (left%maxPolyDegree == right%maxPolyDegree) then
            small = (left%basisType < right%basisType)
          end if
        end if
      end if
    end if

  end function isSmaller
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a <= comparison of two projections.
  !!
  !! Sorting of projections is given by maxPolyDegree, kind, nodes_kind and
  !! last by factor.
  pure function isSmallerOrEqual(left, right) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !

    small = .false.
    if (left%header < right%header) then
      small = .true.
    else
      if (left%header == right%header) then
        if (left%maxPolyDegree < right%maxPolyDegree) then
          small = .true.
        else
          if (left%maxPolyDegree == right%maxPolyDegree) then
            small = (left%basisType <= right%basisType)
          end if
        end if
      end if
    end if

  end function isSmallerOrEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a > comparison of two projections.
  !!
  !! Sorting of projections is given by maxPolyDegree, kind, nodes_kind and
  !! last by factor.
  pure function isGreater(left, right) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !

    great = .false.
    if (left%header > right%header) then
      great = .true.
    else
      if (left%header == right%header) then
        if (left%maxPolyDegree > right%maxPolyDegree) then
          great = .true.
        else
          if (left%maxPolyDegree == right%maxPolyDegree) then
            great = (left%basisType > right%basisType)
          end if
        end if
      end if
    end if

  end function isGreater
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a >= comparison of two projections.
  !!
  !! Sorting of projections is given by maxPolyDegree, kind, nodes_kind and
  !! last by factor.
  pure function isGreaterOrEqual(left, right) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !

    great = .false.
    if (left%header > right%header) then
      great = .true.
    else
      if (left%header == right%header) then
        if (left%maxPolyDegree > right%maxPolyDegree) then
          great = .true.
        else
          if (left%maxPolyDegree == right%maxPolyDegree) then
            great = (left%basisType >= right%basisType)
          end if
        end if
      end if
    end if

  end function isGreaterOrEqual
  ! ************************************************************************ !

end module ply_dynarray_project_module
