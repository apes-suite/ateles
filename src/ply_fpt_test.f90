! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2016, 2019-2020 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop for German Research School
! for Simulation Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Peter Vitt, Verena Krupp,
! and Nikhil Anand for University of Siegen.
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

!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program ply_fpt_test
  use env_module,            only: rk, fin_env
  use tem_param_module,      only: PI
  use tem_logging_module,    only: logUnit
  use ply_fpt_header_module, only: ply_fpt_header_type, ply_fpt_header_define
  use ply_legFpt_module,     only: ply_init_legFpt, &
    &                              ply_legFpt_type
  use ply_modg_basis_module, only: ply_legendre_1D
  use tem_general_module,    only: tem_general_type, tem_start

  !mpi!nprocs = 1

  implicit none

  integer :: iDegree, iBSize
  real(kind=rk) :: res, newRes
  type(tem_general_type) :: general

  integer, parameter :: nDegrees = 14
  integer, parameter :: degree(nDegrees) = [  1,  2,  3,   4,   5,   9,  31, &
    &                                        63, 71, 72, 107, 127, 179, 255  ]
  integer, parameter :: nBSizes = 3
  integer, parameter :: blocksize(nBSizes) = [36, 37, 64]

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'ply_fpt_test', &
    &            version  = 'utest',        &
    &            general  = general         )


  res = 0.0_rk
  do iBSize=1,nBSizes
    do iDegree=1,nDegrees
      call ply_check_legToPnt( blocksize     = blocksize(iBsize), &
        &                      maxpolydegree = degree(iDegree),   &
        &                      res           = newRes             )
      if (newRes > res) then
        res = newRes
      end if
  end do
  end do

  if (res < 1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if

  call fin_env()

contains

  subroutine ply_check_legToPnt(blocksize, maxpolydegree, res)
    integer, intent(in) :: blocksize
    integer, intent(in) :: maxPolyDegree
    real(kind=rk), intent(out) :: res
    integer :: iPoint, iPoly
    real(kind=rk), allocatable :: legCoeffs(:)
    real(kind=rk), allocatable :: pntVal(:), refVal(:)
    real(kind=rk), allocatable :: chebPnt(:)
    real(kind=rk), allocatable :: legValChebPnt(:,:)
    type(ply_fpt_header_type) :: header
    type(ply_legFpt_type) :: fpt

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    write(logUnit(10),*) '------- Number of Legendre coefficients: ', maxPolyDegree+1

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs(maxPolyDegree+1))
    legCoeffs(:) = 1.0_rk

    ! Create the Chebyshev nodes on the interval [-1,+1]
    allocate(chebPnt(maxPolyDegree+1))
    do iPoint = 1, maxPolyDegree+1
      chebPnt(iPoint) = (-1.0_rk) * cos(PI/(maxPolyDegree+1)*((iPoint-1.0_rk)+1.0_rk/2.0_rk))
      !write(*,*) 'Cehbyshev point', iPoint, ' is at: ', chebPnt(iPoint)
    end do

    ! define the reference results for the point values (Chebyshev nodes)
    allocate( legValChebPnt(maxPolyDegree+1,maxPolyDegree+1) )
    legValChebPnt(:,:) = ply_legendre_1D(chebPnt, maxPolyDegree)
    allocate(refVal(maxPolyDegree+1))
    refVal(:) = 0.0_rk
    write(logUnit(10),*) 'Calculating reference results ...'
    do iPoly = 1, maxPolyDegree+1
      refVal(:) = refVal(:) + legValChebPnt(iPoly,:) * legCoeffs(iPoly)
    end do
    write(logUnit(10),*) 'Finished'

    ! Init the FPT
    call ply_fpt_header_define( me        = header,   &
      &                         blocksize = blocksize )
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree, &
      &                   nIndeps       = 1,             &
      &                   fpt           = fpt,           &
      &                   header        = header         )

    ! now transform to the Chebyshev nodes
    allocate(pntVal(1:maxPolyDegree+1))
    write(logUnit(10),*) 'Calculating FPT ...'
    call fpt%legToPnt( legCoeffs = legCoeffs,       &
      &                pntVal = pntVal, nIndeps = 1 )
    write(logUnit(10),*) 'Finished'

    !!do iPoint = 1, maxPolyDegree+1
    !!  write(*,*) 'Point: ', chebPnt(iPoint), &
    !!           & ' FPT: ', pntVal(iPoint), &
    !!           & ' Ref.: ', refVal(iPoint), &
    !!           & ' error: ', pntVal(iPoint)-refVal(iPoint)
    !!end do

    ! Write out the point with the largest absolute error
    write(*,*) 'Cheb-Point ',maxloc(abs(pntVal(:) - refVal(:))),  ' has largest error of: ' ,maxval(abs(pntVal(:) - refVal(:)))

    res = maxval(abs(pntVal(:) - refVal(:)))

  end subroutine

end program ply_fpt_test
