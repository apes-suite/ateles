! Copyright (c) 2015 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
! Copyright (c) 2015,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Kay Langhammer, Harald Klimach, Nikhil
! Anand, Tobias Girresser and Peter Vitt for University of Siegen.
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
!> Fast polynomial transformation using the FXTPACK implementation of a
!! fast multipole method.
module ply_fxt_module
  use env_module, only: rk
  use fxt_fwrap, only: fxtf_flptld_type, &
    &                  fxtf_flptld_n2m,  &
    &                  fxtf_flptld_m2n,  &
    &                  fxtf_flptld_init
  use ply_fxt_header_module, only: ply_fxt_header_type

  implicit none

  private

  type ply_fxt_type
    type(fxtf_flptld_type) :: flpt
    real(kind=rk) :: prec
    integer :: ndims
  end type ply_fxt_type


  public :: ply_fxt_type
  public :: ply_init_fxt
  public :: ply_fxt_m2n_1D, ply_fxt_m2n_2D,ply_fxt_m2n_3D
  public :: ply_fxt_n2m_1D, ply_fxt_n2m_2D,ply_fxt_n2m_3D


contains


  ! ************************************************************************ !
  !> Initialize the flpt data structure for fast legendre polynomial
  !! transformation via the fxtpack.
  subroutine ply_init_fxt( fxt, header, degree )
    ! -------------------------------------------------------------------- !
    !> Handle to the resulting fast polynomial table.
    type(ply_fxt_type), intent(out) :: fxt
    type(ply_fxt_header_type), intent(in) :: header
     !> Polynomial degree.
    integer, intent(in) :: degree
    ! -------------------------------------------------------------------- !
    integer :: nPoints
    ! -------------------------------------------------------------------- !

    nPoints = degree + 1

    call fxtf_flptld_init( flpt    = fxt%flpt,   &
      &                    degree  = degree,     &
      &                    nPoints = nPoints,    &
      &                    prec    = header%prec )

  end subroutine ply_init_fxt
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Convert modal data to nodal data in 1D using flpt.
  !!
  !! This encapsualtes the pure C-Interface, with extraction of the array
  !! sizes and dealing with the flpt data.
  !!
  !! Note: The modal and nodal data array sizes need to match the flpt
  !! definitions, provided in the fxtf_flptld_init call.
  subroutine ply_fxt_m2n_1D( fxt, modal_data, nodal_data )
    ! -------------------------------------------------------------------- !
    !> Description of the Fast Legendre Polynomial Transform
    type(ply_fxt_type) :: fxt
    !> Nodal data
    real(kind=rk), target, intent(inout) :: nodal_data(:)
    !> Modal data
    real(kind=rk), target, intent(inout) :: modal_data(:)
    ! -------------------------------------------------------------------- !

    call fxtf_flptld_m2n( flpt       = fxt%flpt,   &
      &                   modal_data = modal_data, &
      &                   nodal_data = nodal_data  )

  end subroutine ply_fxt_m2n_1D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Convert modal data to nodal data in 2D using flpt.
  subroutine ply_fxt_m2n_2D( fxt, modal_data, nodal_data, oversamp_degree )
    ! -------------------------------------------------------------------- !
    !> Description of the Fast Legendre Polynomial Transform
    type(ply_fxt_type) :: fxt
    !> Nodal data
    real(kind=rk), target, intent(inout) :: nodal_data(:)
    !> Modal data
    real(kind=rk), target, intent(inout) :: modal_data(:)
    integer, intent(in) :: oversamp_degree
    ! -------------------------------------------------------------------- !
    integer :: ub, lb, iLine, iColumn, nModesPerDim, msq
    ! -------------------------------------------------------------------- !

    nModesPerDim = (oversamp_degree+1)
    msq = nModesPerDim*nModesPerDim

    do iLine = 1, oversamp_degree+1
      lb = (iLine-1) * (oversamp_degree+1) + 1
      ub = lb + oversamp_degree
      call fxtf_flptld_m2n( flpt       = fxt%flpt,          &
        &                   modal_data = modal_data(lb:ub), &
        &                   nodal_data = nodal_data(lb:ub)  )
    end do

    do iColumn = 1, oversamp_degree+1
      lb = iColumn
      call fxtf_flptld_m2n( flpt       = fxt%flpt,                             &
        &                   modal_data = nodal_data(lb:msq:oversamp_degree+1), &
        &                   nodal_data = modal_data(lb:msq:oversamp_degree+1)  )
    end do
    nodal_data = modal_data

  end subroutine ply_fxt_m2n_2D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Convert modal data to nodal data in 3D using flpt.
  subroutine ply_fxt_m2n_3D( fxt, modal_data, nodal_data, oversamp_degree )
    ! -------------------------------------------------------------------- !
    !> Description of the Fast Legendre Polynomial Transform
    type(ply_fxt_type) :: fxt
    !> Nodal data
    real(kind=rk), target, intent(inout) :: nodal_data(:)
    !> Modal data
    real(kind=rk), target, intent(inout) :: modal_data(:)
    integer, intent(in) :: oversamp_degree
    ! -------------------------------------------------------------------- !
    integer :: ub, lb, iLine, iColumn, nModesPerDim, msq, ntotalDofs
    real(kind=rk), pointer :: tmp_in(:), tmp_out(:)
    ! -------------------------------------------------------------------- !

    nModesPerDim = (oversamp_degree+1)
    msq = nModesPerDim*nModesPerDim
    nTotalDofs =  (oversamp_degree+1)**3
    allocate(tmp_in(nModesPerDim))
    allocate(tmp_out(nModesPerDim))
    tmp_in = -42
    tmp_out = -42

    ! The loop for msq stripes for independent x Dir evaluations
    do iLine = 1, msq
      lb = (iLine-1) * (oversamp_degree+1) + 1
      ub = lb + oversamp_degree
      tmp_in = modal_data(lb:ub)
      call fxtf_flptld_m2n( flpt       = fxt%flpt, &
        &                   modal_data = tmp_in,   &
        &                   nodal_data = tmp_out   )
      nodal_data(lb:ub) = tmp_out
    end do

    ! The loop for msq stripes for independent y Dir evaluations
    do iColumn = 1, msq
      lb = int( (iColumn-1 ) / nModesPerDim ) * msq &
        & + mod( iColumn-1, nModesPerDim )          &
        & + 1
      ub = lb + msq - 1
      call fxtf_flptld_m2n( flpt       = fxt%flpt,                       &
        &                   modal_data = nodal_data(lb:ub:nModesPerDim), &
        &                   nodal_data = modal_data(lb:ub:nModesPerDim)  )
    end do

    ! The loop for msq stripes for independent z Dir evaluations
    ub = nTotalDofs
    do iColumn = 1, msq
      lb = iColumn
      call fxtf_flptld_m2n( flpt       = fxt%flpt,              &
        &                   modal_data = modal_data(lb:ub:msq), &
        &                   nodal_data = nodal_data(lb:ub:msq)  )
    end do

  end subroutine ply_fxt_m2n_3D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Convert nodal data to modal data using flpt.
  !!
  !! This encapsualtes the pure C-Interface, with extraction of the array
  !! sizes and dealing with the flpt data.
  !!
  !! Note: The modal and nodal data array sizes need to match the flpt
  !! definitions, provided in the fxtf_flptld_init call.
  subroutine ply_fxt_n2m_1D( fxt, nodal_data, modal_data )
    ! -------------------------------------------------------------------- !
    !> Description of the Fast Legendre Polynomial Transform
    type(ply_fxt_type) :: fxt
    !> Nodal data
    real(kind=rk), target, intent(inout) :: nodal_data(:)
    !> Modal data
    real(kind=rk), target, intent(inout) :: modal_data(:)
    ! -------------------------------------------------------------------- !

    call fxtf_flptld_n2m( flpt       = fxt%flpt,   &
      &                   nodal_data = nodal_data, &
      &                   modal_data = modal_data  )

  end subroutine ply_fxt_n2m_1D
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_fxt_n2m_2D( fxt, nodal_data, modal_data, oversamp_degree )
    ! -------------------------------------------------------------------- !
    !> Description of the Fast Legendre Polynomial Transform
    type(ply_fxt_type) :: fxt
    !> Nodal data
    real(kind=rk), target :: nodal_data(:)
    !> Modal data
    real(kind=rk), target :: modal_data(:)
    integer, intent(in) :: oversamp_degree
    ! -------------------------------------------------------------------- !
    integer :: ub, lb, iLine, iColumn, nModesPerDim, msq
    ! -------------------------------------------------------------------- !

    nModesPerDim = (oversamp_degree+1)
    msq = nModesPerDim*nModesPerDim

    do iLine = 1, oversamp_degree+1
      lb = (iLine-1) * (oversamp_degree+1) + 1
      ub = lb + oversamp_degree
      call fxtf_flptld_n2m( flpt       = fxt%flpt,          &
        &                   nodal_data = nodal_data(lb:ub), &
        &                   modal_data = modal_data(lb:ub)  )
    end do

    do iColumn = 1, oversamp_degree+1
      lb = iColumn
      call fxtf_flptld_n2m( flpt       = fxt%flpt,                             &
        &                   nodal_data = modal_data(lb:msq:oversamp_degree+1), &
        &                   modal_data = nodal_data(lb:msq:oversamp_degree+1)  )
    end do
    modal_data = nodal_data
  end subroutine ply_fxt_n2m_2D
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_fxt_n2m_3D( fxt, nodal_data, modal_data, oversamp_degree )
    ! -------------------------------------------------------------------- !
    !> Description of the Fast Legendre Polynomial Transform
    type(ply_fxt_type) :: fxt
    !> Nodal data
    real(kind=rk), target, intent(inout) :: nodal_data(:)
    !> Modal data
    real(kind=rk), target, intent(inout) :: modal_data(:)
    integer, intent(in) :: oversamp_degree
    ! -------------------------------------------------------------------- !
    integer :: ub, lb, iLine, iColumn, nModesPerDim, msq, ntotalDofs
    ! -------------------------------------------------------------------- !

    nModesPerDim = (oversamp_degree+1)
    msq = nModesPerDim*nModesPerDim
    nTotalDofs =  (oversamp_degree+1)**3

    ! The loop for msq stripes for independent x Dir evaluations
    do iLine = 1, msq
      lb = (iLine-1) * (oversamp_degree+1) + 1
      ub = lb + oversamp_degree
      call fxtf_flptld_n2m( flpt       = fxt%flpt,          &
        &                   nodal_data = nodal_data(lb:ub), &
        &                   modal_data = modal_data(lb:ub)  )
    end do

    ! The loop for msq stripes for independent y Dir evaluations
    do iColumn = 1, msq
      lb = int( (iColumn-1) / nModesPerDim ) * msq &
        & + mod( iColumn-1, nModesPerDim )         &
        & + 1
      ub = lb + msq - 1
      call fxtf_flptld_n2m( flpt       = fxt%flpt,                       &
        &                   nodal_data = modal_data(lb:ub:nModesPerDim), &
        &                   modal_data = nodal_data(lb:ub:nModesPerDim)  )
    end do

    ! The loop for msq stripes for independent z Dir evaluations
    ub = nTotalDofs
    do iColumn = 1, msq
      lb = iColumn
      call fxtf_flptld_n2m( flpt       = fxt%flpt,              &
        &                   nodal_data = nodal_data(lb:ub:msq), &
        &                   modal_data = modal_data(lb:ub:msq)  )
    end do

  end subroutine ply_fxt_n2m_3D
  ! ************************************************************************ !


end module ply_fxt_module
