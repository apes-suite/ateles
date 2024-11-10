! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2013-2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
!
! Parts of this file were written by Verena Krupp, Peter Vitt, Harald Klimach,
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
module ply_prj_header_module
  use env_module,               only: labelLen

  use fftw_wrap,                only: fftw_available

  use aotus_module,             only: flu_State, aot_get_val
  use aot_out_module,           only: aot_out_type, aot_out_val

  use tem_aux_module,           only: tem_abort
  use tem_tools_module,         only: upper_to_lower
  use tem_logging_module,       only: logUnit

  use ply_fpt_header_module
  use ply_l2p_header_module
  use ply_fxt_header_module

  implicit none

  private

  !> Configurable projection settings.
  type ply_prj_header_type
    !> Kind of projection. Currently available:
    !! - 'l2p', L2-Projection
    !! - 'fpt', Fast Polynomial Transformation. Requires the FFTW.
    character(len=labelLen) :: kind
    type(ply_fpt_header_type) :: fpt_header
    type(ply_l2p_header_type) :: l2p_header
    type(ply_fxt_header_type) :: fxt_header
  end type ply_prj_header_type

  interface assignment(=)
    module procedure copy_poly_proj_header
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

  public :: ply_prj_header_type
  public :: ply_prj_header_load
  public :: ply_prj_header_out
  public :: assignment(=)
  public :: operator(==), operator(/=), operator(<), operator(<=)
  public :: operator(>), operator(>=)


contains

  ! ************************************************************************ !
  pure subroutine Copy_poly_proj_header( left, right )
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_prj_header_type), intent(out) :: left
    !> fpt to copy from
    type(ply_prj_header_type), intent(in) :: right
    ! -------------------------------------------------------------------- !
    left%kind = right%kind
    left%fpt_header = right%fpt_header
    left%l2p_header = right%l2p_header
    left%fxt_header = right%fxt_header
  end subroutine copy_poly_proj_header
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load settings to describe a projection method from a Lua table.
  subroutine ply_prj_header_load( me, conf, parent )
    ! -------------------------------------------------------------------- !
    type(ply_prj_header_type), intent(out) :: me
    type(flu_State) :: conf
    !> A parent Lua table, in which the boundary conditions are to be found.
    integer, intent(in) :: parent
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    if (parent /= 0) then

      call aot_get_val(L       = conf,    &
        &              thandle = parent,  &
        &              key     = 'kind',  &
        &              val     = me%kind, &
        &              default = 'l2p',   &
        &              ErrCode = iError   )
      me%kind = upper_to_lower(me%kind)

      select case(trim(me%kind))
      case('l2p')
        call ply_l2p_header_load( me      = me%l2p_header, &
          &                       conf    = conf,          &
          &                       thandle = parent         )
        call ply_l2p_header_display(me = me%l2p_header)

      case('fxt')
        call ply_fxt_header_load( me      = me%fxt_header, &
          &                       conf    = conf,          &
          &                       thandle = parent         )
        call ply_fxt_header_display(me = me%fxt_header     )

      case('fpt')
        if (fftw_available) then
          call ply_fpt_header_load( me      = me%fpt_header, &
            &                       conf    = conf,          &
            &                       thandle = parent         )
          call ply_fpt_header_display (me = me%fpt_header)
        else
          write(logUnit(1),*) ''
          write(logUnit(1),*)                                       &
            & '+===================================================+'
          write(logUnit(1),*)                                       &
            & '!! FFTW NOT available but necessary for FPT!       !!'
          write(logUnit(1),*)                                       &
            & '!! WARNING: Deactivating fast polynomial transform !!'
          write(logUnit(1),*)                                       &
            & '!!          for this projection!                   !!'
          write(logUnit(1),*)                                       &
            & '!!                                                 !!'
          write(logUnit(1),*)                                       &
            & '!! Falling back to L2 Projection.                  !!'
          write(logUnit(1),*)                                       &
            & '+===================================================+'
          write(logUnit(1),*) ''
          me%kind='l2p'
          call ply_l2p_header_load( me      = me%l2p_header, &
             &                      conf    = conf,          &
             &                      thandle = parent         )
          call ply_l2p_header_display(me = me%l2p_header)
        end if

      case default
        write(logUnit(1),*) 'ERROR while loading projection:'
        write(logUnit(1),*) '      Unknown projection method ' &
          & // trim(me%kind) // '!'
        write(logUnit(1),*) '      Available methods are:'
        write(logUnit(1),*) '      * l2p - L2 Projection'
        write(logUnit(1),*) '      * fxt - FXTPACK: Fast Multipole Method'
        if (fftw_available) then
          write(logUnit(1),*) '      * fpt - Fast Polynomial Transformation'
        end if
        call tem_abort()
      end select
    else
      write(logUnit(1),*) 'No projection provided, using defaults.'
      me%kind = 'l2p'
      call ply_l2p_header_load( me      = me%l2p_header, &
        &                       conf    = conf,          &
        &                       thandle = parent         )
      call ply_l2p_header_display( me = me%l2p_header )
    end if

  end subroutine ply_prj_header_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load settings to describe a projection method from a Lua table.
  subroutine ply_prj_header_out( me, conf )
    ! -------------------------------------------------------------------- !
    type(ply_prj_header_type), intent(in) :: me
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !

    call aot_out_val( put_conf = conf,   &
      &               vname    = 'kind', &
      &               val      = me%kind )

    select case(trim(me%kind))
    case('l2p')
      call ply_l2p_header_out( me   = me%l2p_header, &
        &                      conf = conf           )

    case('fpt')
      call ply_fpt_header_out( me   = me%fpt_header, &
        &                      conf = conf           )

    case('fxt')
      call ply_fxt_header_out( me   = me%fxt_header, &
        &                      conf = conf           )

    end select

  end subroutine ply_prj_header_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for equality of the header for two
  !! projections.
  !!
  !! The headers are considered to be equal, if their kind and the corresponding
  !! headers are equal. For unknown kinds, the headers are not taken into
  !! consideration.
  pure function isEqual( left, right ) result(equality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_header_type), intent(in) :: right
    !> is equal??
    logical :: equality
    ! -------------------------------------------------------------------- !

    equality = ( left%kind == right%kind )
    select case(left%kind)
      case ('fpt')
        equality = equality .and. ( left%fpt_header == right%fpt_header )
      case ('l2p')
        equality = equality .and. ( left%l2p_header == right%l2p_header )
      case ('fxt')
        equality = equality .and. ( left%fxt_header == right%fxt_header )
    end select

  end function isEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for unequality of the header of two
  !! projections.
  !!
  !! Two projections are considered to be unequal, if their kind, their
  !! fpt-header or l2p_header are not equal.
  pure function isUnequal( left, right ) result(unequality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_header_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    ! -------------------------------------------------------------------- !

    unequality = ( left%kind /= right%kind )

    if (.not. unequality) then
      select case(left%kind)
        case ('fpt')
          unequality = ( left%fpt_header /= right%fpt_header )
        case ('l2p')
          unequality = ( left%l2p_header /= right%l2p_header )
        case ('fxt')
          unequality = ( left%fxt_header /= right%fxt_header )
      end select
    end if

  end function isUnequal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a < comparison of the header of two projections.
  !!
  !! Sorting of projections is given by the kind, fpt_header and
  !! last by l2p_header.
  pure function isSmaller(left, right) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !

    small = .false.
    if (left%kind < right%kind) then
      small = .true.
    else
      if (left%kind == right%kind) then
        select case(left%kind)
          case ('fpt')
            small = ( left%fpt_header < right%fpt_header )
          case ('l2p')
            small = ( left%l2p_header < right%l2p_header )
          case ('fxt')
            small = ( left%fxt_header < right%fxt_header )
        end select
      end if
    end if

  end function isSmaller
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a <= comparison of the header of two projections.
  !!
  !! Sorting of projections is given by kind, fpt_header and
  !! last by the l2p header.
  pure function isSmallerOrEqual( left, right ) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !

    small = .false.
    if (left%kind < right%kind) then
      small = .true.
    else
      if (left%kind == right%kind) then
        select case(left%kind)
          case ('fpt')
            small = ( left%fpt_header <= right%fpt_header )
          case ('l2p')
            small = ( left%l2p_header <= right%l2p_header )
          case ('fxt')
            small = ( left%fxt_header <= right%fxt_header )
        end select
      end if
    end if

  end function isSmallerOrEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a > comparison of the header of two projections.
  !!
  !! Sorting of projections is given by kind, fpt_header and
  !! last by l2p_header.
  pure function isGreater( left, right ) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !

    great = .false.
    if (left%kind > right%kind) then
      great = .true.
    else
      if (left%kind == right%kind) then
        select case(left%kind)
          case ('fpt')
            great = ( left%fpt_header > right%fpt_header )
          case ('l2p')
            great = ( left%l2p_header > right%l2p_header )
          case ('fxt')
            great = ( left%fxt_header > right%fxt_header )
        end select
      end if
    end if

  end function isGreater
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a >= comparison of the header of two projections.
  !!
  !! Sorting of projections is given by kind, fpt_header and
  !! last by l2p_header.
  pure function isGreaterOrEqual( left, right ) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_prj_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !

    great = .false.
    if (left%kind > right%kind) then
      great = .true.
    else
      if (left%kind == right%kind) then
        select case(left%kind)
          case ('fpt')
            great = ( left%fpt_header >= right%fpt_header )
          case ('l2p')
            great = ( left%l2p_header >= right%l2p_header )
          case ('fxt')
            great = ( left%fxt_header >= right%fxt_header )
        end select
      end if
    end if

  end function isGreaterOrEqual
  ! ************************************************************************ !

end module ply_prj_header_module
