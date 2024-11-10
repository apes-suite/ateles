! Copyright (c) 2015 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
!
! Parts of this file were written by Kay Langhammer, Nikhil Anand, Harald
! Klimach and Robin Weihe for University of Siegen.
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

program fxtp_test
  use iso_c_binding
  use env_module, only: rk
  use fxt_fwrap, only: fxtf_flptld_type, fxtf_flptld_init, &
    &                  fxtf_flptld_exp, fxtf_flptld_evl

  implicit none

  type(fxtf_flptld_type), target :: flpt

  real(kind=rk) :: nodal_orig(10)
  real(c_double) :: modal(10)
  real(c_double) :: nodal(10)
  integer(c_int):: modalLen, nodalLen

  integer :: nPoints        ! number of points
  integer :: maxDegree      ! maximum degree
  real(kind=rk) :: prec

  nPoints = 10
  maxDegree = 9
  prec = 8*epsilon(prec)
  modalLen = nPoints
  nodalLen = nPoints

  call fxtf_flptld_init(flpt    = flpt,      &
    &                   degree  = maxDegree, &
    &                   nPoints = nPoints,   &
    &                   prec    = prec       )

  call random_number(nodal_orig)
  write(*,*) 'orig :', nodal_orig
  nodal = nodal_orig

  write(*,*) "passing nodal values = ", nodal
  write(*,*) "modalLen, nodalLen = ", modalLen, nodalLen
  ! ther e....
  ! transform from physical to wave space
  call fxtf_flptld_exp( modal, modalLen, flpt%handle, &
    &                   nodal, nodalLen, flpt%work    )

  write(*,*) "modal val", modal
  nodal = 0.0

  ! ...and back again
  ! transform from wave to physical space
  call fxtf_flptld_evl( nodal, nodalLen, flpt%handle, &
    &                   modal, modalLen, flpt%work    )

  write(*,*) 'trafo:', nodal
  write(*,*) 'Should be the same as orig.'

  if (all(abs(nodal - nodal_orig) < 8.*epsilon(1.0_c_double))) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'Data does not match after conversion:'
    write(*,*) 'max error:', maxval(abs(nodal - nodal_orig))
    write(*,*) 'FAILED'
  end if

end program fxtp_test
