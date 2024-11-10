! Copyright (c) 2015, 2019 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2015 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach and Kay Langhammer
! for University of Siegen.
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

!> This module provides the ISO_C_Binding interfaces to the fxtpack routines.
!!
!! falt: Fast Associated Legendre Transform (spherical harmonics)
!! flpt: Fast Legendre Polynomial Transform
!! vecld: Array of doubles
module fxt_fif
  use, intrinsic :: iso_c_binding
  implicit none

  interface

    ! FALTLD routines:
    ! /*** load fast spherical harmonic transform ***/
    ! fxt_faltld* fxt_faltld_load(char *fname)
!    function fxt_faltld_load(fname) result(falt) bind(c)
!      use, intrinsic :: iso_c_binding
!      character(c_char) :: fname
!      type(c_ptr) :: falt
!    end function fxt_faltld_load

    subroutine fxt_faltld_preproc(p, n, mv, prec, fname) bind(c)
      use, intrinsic :: iso_c_binding
      integer(kind=c_long), value :: p
      integer(kind=c_long), value :: n
      type(c_ptr) :: mv
      real(kind=c_double), value :: prec
      character(kind=c_char) :: fname
    end subroutine fxt_faltld_preproc

! TODO: write the routine in C
!    type(c_ptr) function fxt_faltld_init(p, n, prec) bind(c)
!      use, intrinsic :: iso_c_binding
!      integer(c_long), value :: p
!      integer(c_long), value :: n
!      real(c_double), value :: prec
!    end function fxt_faltld_init

    ! deallocate fast spherical harmonic transform
    ! void fxt_faltld_del(fxt_faltld *falt);
    subroutine fxt_faltld_del(falt) bind (c)
     use, intrinsic :: iso_c_binding
     type(c_ptr), value :: falt           !fxt_faltld
    end subroutine fxt_faltld_del

    ! size of working array
    ! long fxt_faltld_wsize(fxt_faltld *falt, long m);
    function fxt_faltld_wsize(falt, m) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: falt            !fxt_faltld
      integer(kind=c_long), value :: m
      integer(kind=c_long) :: fxt_faltld_wsize
    end function fxt_faltld_wsize

    !  /*** maximum size of working array ***/
    ! long fxt_faltld_wsizemax(fxt_faltld *falt);
    function fxt_faltld_wsizemax(falt) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: falt          !fxt_faltld
      integer(kind=c_long) :: fxt_faltld_wsizemax
    end function fxt_faltld_wsizemax

    ! /*** evaluate fast spherical harmonic transform ***/
    ! void fxt_faltld_evl(fxt_vecld *v, fxt_faltld *falt, long m,
    !                     fxt_vecld *u, fxt_vecld *w);
    subroutine fxt_faltld_evl(v, falt, m, u, w) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: u, v, w                 !fxt_vecld
       type(c_ptr) :: falt                   !fxt_faltld
       integer(kind=c_long), value :: m
    end subroutine fxt_faltld_evl

    !  /*** expand fast spherical harmonic transform ***/
    ! void fxt_faltld_exp(fxt_vecld *u, fxt_faltld *falt, long m,
    !                     fxt_vecld *v, fxt_vecld *w);
    subroutine fxt_faltld_exp(u, falt, m, v, w) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: u, v, w            !fxt_vecld
       type(c_ptr) :: falt              !fxt_faltld
       integer(kind=c_long), value :: m
    end subroutine fxt_faltld_exp
    ! ------------------------------------------------------------------------ !


    ! ........................................................................ !
    ! FLPTLD routines:
    subroutine fxt_flptld_preproc(p, n, prec, fname) bind(c)
      use, intrinsic :: iso_c_binding
      integer(kind=c_long), value :: p
      integer(c_long), value :: n
      real(c_double), value :: prec
      character(c_char) :: fname
    end subroutine fxt_flptld_preproc

    function fxt_flptld_init(p, n, prec) bind(c)
      use, intrinsic :: iso_c_binding
      integer(kind=c_long), value :: p
      integer(kind=c_long), value :: n
      real(c_double), value :: prec
      type(c_ptr) :: fxt_flptld_init
    end function fxt_flptld_init

    ! deallocate fast Legendre polynomial transform
    ! void fxt_flptld_del(fxt_flptld *flpt);
    subroutine fxt_flptld_del(flpt) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: flpt             !fxt_flptld
    end subroutine fxt_flptld_del

    ! size of working array
    ! long fxt_flptld_wsize(fxt_flptld *flpt);
    function fxt_flptld_wsize(flpt) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: flpt              !fxt_flptld
      integer(kind=c_long) :: fxt_flptld_wsize
    end function fxt_flptld_wsize

    ! evaluate fast Legendre Polynomial transform
    ! void fxt_flptld_evl(fxt_vecld *v, fxt_flptld *flpt,
    !                     fxt_vecld *u, fxt_vecld *w);
    subroutine fxt_flptld_evl(v, flpt, u, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: v, u, w              ! fxt_vecld
      type(c_ptr), value :: flpt               ! fxt_flptld
    end subroutine fxt_flptld_evl

    ! expand fast Legendre Polynomial transform
    ! void fxt_flptld_exp(fxt_vecld *u, fxt_flptld *flpt,
    !                     fxt_vecld *v, fxt_vecld *w);
    subroutine fxt_flptld_exp(u, flpt, v, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value ::  u, v, w        ! fxt_vecld
      type(c_ptr), value :: flpt            ! fxt_flptld
    end subroutine fxt_flptld_exp
    ! ------------------------------------------------------------------------ !


    ! ........................................................................ !
    ! VECL / VECLD
    function fxt_vecl_new(size) bind(c)
      use, intrinsic :: iso_c_binding
      integer(kind=c_long) :: size
      type(c_ptr) :: fxt_vecl_new
    end function fxt_vecl_new

    function fxt_vecld_new(size) bind(c)
      use, intrinsic :: iso_c_binding
      integer(kind=c_long), value :: size
      type(c_ptr) :: fxt_vecld_new
    end function fxt_vecld_new
    ! ------------------------------------------------------------------------ !

  end interface

end module fxt_fif
