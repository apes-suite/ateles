! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
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

!> Unit test to check functionallity of fast polynomial transformations
!! to nodal values at Lobatto-Chebyshev-Nodes.
!! \author{Jens Zudrop}
program ply_ifpt_lobattoNodes_test
  use env_module,            only: rk, fin_env
  use tem_param_module,      only: PI
  use tem_logging_module,    only: logUnit
  use tem_aux_module,        only: tem_abort
  use tem_general_module,    only: tem_general_type, tem_start
  use ply_fpt_header_module, only: ply_fpt_header_type, ply_fpt_header_define
  use ply_legFpt_module,     only: ply_init_legFpt, ply_legFpt_type
  use ply_modg_basis_module, only: ply_legendre_1D

  !mpi!nprocs = 1

  implicit none

  integer :: iPower
  real(kind=rk) :: res, newRes
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = general             )

  res = 0.0_rk
  do iPower = 1,4
    call ply_check_pntToLeg(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if
  call fin_env()

contains

  subroutine ply_check_pntToLeg(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    integer :: maxPolyDegree, iPoint, iPoly
    real(kind=rk), allocatable :: legCoeffs(:)
    real(kind=rk), allocatable :: pntVal(:), legVal(:)
    real(kind=rk), allocatable :: chebPnt(:)
    real(kind=rk), allocatable :: legValChebPnt(:,:)
    type(ply_fpt_header_type) :: header
    type(ply_legFpt_type) :: fpt

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree =  2**power-1  ! maxPolyDegree+1 has to be a power of 2
    write(logUnit(10),*) '------- Number of Legendre coefficients: ', maxPolyDegree+1

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs(1:maxPolyDegree+1))
    legCoeffs(:) = 1.0_rk

    ! Create the Chebyshev nodes on the interval [-1,+1]
    allocate(chebPnt(maxPolyDegree+1))
    do iPoint = 1, maxPolyDegree+1
      chebPnt(iPoint) = cos((iPoint-1.0_rk)*PI/maxPolyDegree);
      !write(*,*) 'Cehbyshev point', iPoint, ' is at: ', chebPnt(iPoint)
    end do

    ! define the point values (Lobatto-Chebyshev-nodes)
    allocate( legValChebPnt(maxPolyDegree+1,maxPolyDegree+1) )
    legValChebPnt(:,:) = ply_legendre_1D(chebPnt, maxPolyDegree)
    allocate(pntVal(maxPolyDegree+1))
    pntVal(:) = 0.0_rk
    write(logUnit(10),*) 'Calculating point values (input) ...'
    do iPoly = 1, maxPolyDegree+1
      pntVal(:) = pntVal(:) + legValChebPnt(iPoly,:) * legCoeffs(iPoly)
    end do
    write(logUnit(10),*) 'Finished'

    ! Init the FPT
    call ply_fpt_header_define( me = header,           &
      &                         lobattoPoints = .true. )
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree, &
      &                   nIndeps       = 1,             &
      &                   fpt           = fpt,           &
      &                   header        = header         )

    ! now transform to the Legendre coefficients
    allocate(legVal(1:maxPolyDegree+1))
    write(logUnit(10),*) 'Calculating inverse FPT ...'
    call fpt%pntToLeg( pntVal = pntVal, legCoeffs = legVal, &
      &                nIndeps = 1                          )
    write(logUnit(10),*) 'Finished'

    !!do iPoly = 1, maxPolyDegree+1
    !!  write(*,*) 'Poly degree: ', iPoly, &
    !!           & ' iFPT: ', legVal(iPoly), &
    !!           & ' Ref.: ', legCoeffs(iPoly), &
    !!           & ' error: ', legVal(iPoly)-legCoeffs(iPoly)
    !!end do

    ! Write out the polynomial coefficient with the largest absolute error
    write(*,*) 'Leg-Coefficient ',maxloc(abs(legVal(:) - legCoeffs(:))),  &
              & ' has largest error of: ' ,maxval(abs(legVal(:) - legCoeffs(:)))

    res = maxval(abs(legVal(:) - legCoeffs(:)))

  end subroutine

end program ply_ifpt_lobattoNodes_test
