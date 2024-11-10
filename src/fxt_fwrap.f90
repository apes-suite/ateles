! Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Nikhil Anand
! and Kay Langhammer for University of Siegen.
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

module fxt_fwrap
  use, intrinsic :: iso_c_binding
  use env_module, only: rk
  use fxt_fif, only: fxt_flptld_init, fxt_flptld_wsize, fxt_vecld_new

  implicit none

  private

  !> This datatype provides a handle to the information that FXTPACK needs
  !! to have about the transformation.
  type fxtf_flptld_type
    !> Handle for the fast Legendre polynomial transformation data in FXTPACK.
    type(c_ptr) :: handle

    !> Pointer to the working array, that is required by the transformations.
    type(c_ptr) :: work
  end type fxtf_flptld_type


  !> Interface declarations to the fxtf_wrapper.c routines.
  !!
  !! Those routines enable the passing of Fortran arrays to the FXTPACK and
  !! take care of putting the data into the fxt_vecld data structures, which are
  !! then passed on to the actual fxt_* routines.
  interface
    subroutine fxtf_flptld_evl(v, vn, flpt, u, un, w) bind(c)
      use, intrinsic :: iso_c_binding
      real(c_double), dimension(*) :: v
      integer(c_int), value :: vn
      type(c_ptr), value :: flpt
      real(c_double), dimension(*) :: u
      integer(c_int), value :: un
      type(c_ptr), value :: w
    end subroutine fxtf_flptld_evl

    subroutine fxtf_flptld_exp(u, un, flpt, v, vn, w) bind(c)
      use, intrinsic :: iso_c_binding
      real(c_double), dimension(*)  ::  u
      integer(c_int), value  :: un
      type(c_ptr), value  :: flpt
      real(c_double), dimension(*)  :: v
      integer(c_int), value :: vn
      type(c_ptr), value  :: w
    end subroutine fxtf_flptld_exp

    subroutine fxtf_faltld_evl(v, vn, falt, m, u, un, w) bind(c)
      use, intrinsic :: iso_c_binding
      real(c_double), dimension(*) :: v
      integer(c_int), value :: vn
      type(c_ptr), value :: falt
      integer(kind=c_long), value :: m
      real(c_double), dimension(*) :: u
      integer(c_int), value :: un
      type(c_ptr), value :: w
    end subroutine fxtf_faltld_evl

    subroutine fxtf_faltld_exp(u, un, falt, m, v, vn, w) bind(c)
      use, intrinsic :: iso_c_binding
      real(c_double), dimension(*) :: u
      integer(c_int), value :: un
      type(c_ptr), value :: falt
      integer(kind=c_long), value :: m
      real(c_double), dimension(*) :: v
      integer(c_int), value :: vn
      type(c_ptr), value :: w
    end subroutine fxtf_faltld_exp

  end interface

  public :: fxtf_flptld_type
  public :: fxtf_flptld_init
  public :: fxtf_flptld_evl
  public :: fxtf_flptld_exp
  public :: fxtf_flptld_m2n, fxtf_flptld_n2m


contains


  !> Convert modal data to nodal data using flpt.
  !!
  !! This encapsualtes the pure C-Interface, with extraction of the array
  !! sizes and dealing with the flpt data.
  !!
  !! Note: The modal and nodal data array sizes need to match the flpt
  !! definitions, provided in the fxtf_flptld_init call.
  subroutine fxtf_flptld_m2n(flpt, modal_data, nodal_data)
    !> Description of the Fast Legendre Polynomial Transform
    type(fxtf_flptld_type), intent(in) :: flpt
    !> Modal data
    real(kind=rk), target :: modal_data(:)
    !> Nodal data
    real(kind=rk), target :: nodal_data(:)

    integer(kind=c_int) :: nNodes, nModes

    nNodes = size(nodal_data)
    nModes = size(modal_data)

    call fxtf_flptld_evl( nodal_data, nNodes, flpt%handle, &
      &                   modal_data, nModes, flpt%work    )

  end subroutine fxtf_flptld_m2n


  !> Convert nodal data to modal data using flpt.
  !!
  !! This encapsualtes the pure C-Interface, with extraction of the array
  !! sizes and dealing with the flpt data.
  !!
  !! Note: The modal and nodal data array sizes need to match the flpt
  !! definitions, provided in the fxtf_flptld_init call.
  subroutine fxtf_flptld_n2m(flpt, nodal_data, modal_data)
    !> Description of the Fast Legendre Polynomial Transform
    type(fxtf_flptld_type) :: flpt
    !> Nodal data
    real(kind=rk), target :: nodal_data(:)
    !> Modal data
    real(kind=rk), target :: modal_data(:)

    integer(kind=c_int) :: nNodes, nModes

    nNodes = size(nodal_data)
    nModes = size(modal_data)

    call fxtf_flptld_exp( modal_data, nModes, flpt%handle, &
      &                   nodal_data, nNodes, flpt%work    )

  end subroutine fxtf_flptld_n2m


  !> Initialize the flpt data structure for fast legendre polynomial
  !! transformation via the fxtpack.
  subroutine fxtf_flptld_init(flpt, degree, nPoints, prec)
    !> Handle to the resulting fast polynomial table.
    type(fxtf_flptld_type), intent(out) :: flpt

    !> Polynomial degree.
    integer, intent(in) :: degree

    !> Number of points.
    !!
    !! Optional, defaults to degree+1.
    integer, intent(in), optional :: nPoints

    !> Required precision for the transformation.
    !!
    !! Optional, defaults to 8 times the precision of c_double.
    real(kind=c_double), intent(in), optional :: prec

    integer(c_long) :: wsize
    integer(c_long) :: p
    integer(c_long) :: n

    real(kind=c_double) :: lprec

    n = degree
    if (present(nPoints)) then
      p = nPoints
    else
      p = degree + 1
    end if

    if (present(prec)) then
      lprec = prec
    else
      lprec = 8*epsilon(lprec)
    end if

    flpt%handle = fxt_flptld_init(p, n, lprec)
    wsize = fxt_flptld_wsize(flpt%handle)
    flpt%work = fxt_vecld_new(wsize)

  end subroutine fxtf_flptld_init

end module fxt_fwrap
