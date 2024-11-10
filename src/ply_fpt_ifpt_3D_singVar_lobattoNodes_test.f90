! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016, 2018-2020 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2013-2014 Peter Vitt <peter.vitt2@uni-siegen.de>
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
program ply_fpt_ifpt_3D_singVar_lobattoNodes_test
  use env_module,               only: rk, fin_env
  use tem_general_module,       only: tem_general_type, tem_start
  use tem_logging_module,       only: logUnit
  use ply_fpt_header_module,    only: ply_fpt_header_type, ply_fpt_header_define
  use ply_legFpt_module,        only: ply_legFpt_type, ply_init_legFPT
  use ply_legFpt_3D_module,     only: ply_legToPnt_3D, &
    &                                 ply_pntToLeg_3D

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
  do iPower = 1,3
    call ply_check_legToPnt_3D(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if

  call fin_env()


contains


  subroutine ply_check_legToPnt_3D(power,res)
    integer, intent(in) :: power
    real(kind=rk) :: res
    integer :: maxPolyDegree, iVar, nVars, iDof
    real(kind=rk), allocatable :: legCoeffs(:,:), legCoeffsIn(:,:)
    real(kind=rk), allocatable :: pntVal(:,:), legVal(:,:)
    type(ply_fpt_header_type) :: header
    type(ply_legFpt_type) :: fpt

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree =  2**power-1  ! maxPolyDegree+1 has to be a power of 2
    nVars = 3
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (per dim): ', maxPolyDegree+1
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (total): ',(maxPolyDegree+1)**3

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs((maxPolyDegree+1)**3,nVars))
    allocate(legCoeffsIn((maxPolyDegree+1)**3,nVars))
    do iVar = 1, nVars
      legCoeffs(:,iVar) = real(iVar, rk)
    end do

    ! Init the FPT
    call ply_fpt_header_define( me = header,           &
      &                         lobattoPoints = .true. )
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree,        &
      &                   nIndeps       = (maxPolyDegree+1)**2, &
      &                   fpt           = fpt,                  &
      &                   header        = header                )

    ! now transform to the Chebyshev nodes
    allocate(pntVal( (maxPolyDegree+1)**3, nVars ))
    legCoeffsIn = legCoeffs
    write(logUnit(10),*) 'Calculating FPT ...'
    do iVar=1,nVars
      call ply_legToPnt_3D( fpt = fpt, legCoeffs = legCoeffsIn(:,iVar), &
        &                   pntVal = pntVal(:,iVar) )
    end do
    write(logUnit(10),*) 'Finished'

    ! now transform back to Legendre coefficients
    allocate(legVal( (maxPolyDegree+1)**3,nVars ))
    write(logUnit(10),*) 'Calculating inverse FPT ...'
    do iVar=1,nVars
      call ply_pntToLeg_3D( fpt = fpt, pntVal = pntVal(:,iVar), &
        &                   legCoeffs = legVal(:,iVar)          )
    end do
    write(logUnit(10),*) 'Finished'


    ! Write out the coefficient with the largest absolute error
    do iVar = 1, nVars
      do iDof = 1, size(legVal(:,iVar))
        write(logUnit(10),*) legVal(iDof,iVar), legCoeffs(iDof,iVar)
      end do
      write(logUnit(10),*) 'For var ', iVar, &
               & ' Leg-Coeff ',maxloc(abs(legVal(:,iVar) - legCoeffs(:,iVar))),  &
               & ' has largest error of: ' ,maxval(abs(legVal(:,iVar) - legCoeffs(:,iVar)))
    end do

    res = maxval(abs(legVal(:,:) - legCoeffs(:,:)))

  end subroutine

end program ply_fpt_ifpt_3D_singVar_lobattoNodes_test
