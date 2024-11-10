! Copyright (c) 2013-2016,2019 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2013-2014, 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Peter Vitt, Verena Krupp
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

! This is the unit test for the projection module.
program ply_project_fpt_test
  use env_module,              only: rk, fin_env
  use tem_logging_module,      only: logUnit
  use ply_dof_module,          only: Q_space
  use ply_dynArray_project_module, only: ply_prj_init_define, &
    &                                    ply_prj_init_type
  use ply_poly_project_module,     only: ply_poly_project_fillbody, &
    &                                    ply_poly_project_m2n,      &
    &                                    ply_poly_project_n2m,      &
    &                                    ply_poly_project_type
  use ply_prj_header_module,       only: ply_prj_header_type
  use ply_fpt_header_module,       only: ply_fpt_header_define
  use tem_general_module,          only: tem_general_type, tem_start

  !mpi!nprocs = 1

  implicit none


  real(kind=rk) :: res, newRes
  integer :: power
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = general             )

  res = 0.0_rk

  !>\todo Check reading of the projection configuration.
  !!      The input has to be configured on the fly!

  !>\todo Put various projections into a projection descriptor.

  !>\todo Check those projections by doing m2n and n2m, where
  !!      input and output should be the same within certain
  !!      bounds.

  ! check poly project module with fpt
  do power = 1,7
    write(logUnit(10),*) '---------------------------   CHECKING CHEB->LEG FPT TRAFO FOR ', 2**power
    call check_fpt(power, newRes)
    if (newRes.gt.res) then
      res = newRes
    end if
    write(logUnit(10),*) '--------------------------- DONE'
  end do

  !>\todo If everything worked fine, write PASSED on the very last line of output, to
  !!      indicate a successful run of the unit test:
  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if
  call fin_env()

contains

  subroutine check_fpt(power, res)
    !---------------------------------------
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    !---------------------------------------
    type(ply_poly_project_type) :: me
    type(ply_prj_init_type) :: prj_init
    type(ply_prj_header_type) :: header
    real(kind=rk),allocatable     :: nodal_data(:,:)
    real(kind=rk), allocatable    :: modal_data(:,:)
    real(kind=rk), allocatable    :: ref_modes(:)
    !-----------for init
    integer :: basisType, maxdegree, i
    !-----------for oversamp
    real(kind=rk), allocatable    :: oversamp_modal(:,:)
    integer :: iDegX, iDegY, iDegZ,dof, dofOverSamp

    basisType = Q_space
    maxdegree = power
    header%kind = 'fpt'
    call ply_fpt_header_define( me = header%fpt_header, &
      &                         approx_terms = 18,      &
      &                         striplen     = 256,     &
      &                         lobattoPoints = .false. )

    ! define poly projection init type
    call ply_prj_init_define( me            = prj_init,  &
      &                       header        = header,    &
      &                       maxPolyDegree = maxdegree, &
      &                       basisType     = basistype  )

    ! fill the projection body according to the header
    call ply_poly_project_fillbody(me = me, proj_init=prj_init, scheme_dim=3)

    allocate(ref_modes(1:me%body_3d%ndofs))
    allocate(modal_data(1:me%body_3d%ndofs,1))
    allocate(oversamp_modal(1:me%body_3d%oversamp_dofs,1))
    allocate(nodal_data(1:me%body_3d%nQuadPoints,1))

    ref_modes = 0.0_rk
    do i=1, me%body_3d%ndofs
       ref_modes(i)=1.0/real(i, kind=rk)
    end do

    modal_data(:,:) = 0.0_rk
    oversamp_modal(:,:) = 0.0_rk
    ! oversampling of modes
    do iDegX = 1, maxdegree+1
      do iDegY = 1, maxdegree+1
        do iDegZ = 1, maxdegree+1
          dof = 1 + (iDegX-1) + (iDegY-1)*(maxdegree+1) + &
           &    (iDegZ-1) * (maxdegree+1) * (maxdegree+1)
          dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(me%oversamp_degree+1)+ &
           &            (iDegZ-1) * (me%oversamp_degree+1)*&
           &            (me%oversamp_degree+1)
          oversamp_modal(dofOverSamp,1) = ref_modes(dof)
        end do
      end do
    end do

    call ply_poly_project_m2n(me = me ,                &
      &                       dim = 3 ,                &
      &                       nVars = 1,               &
      &                       nodal_data=nodal_data,   &
      &                       modal_data=oversamp_modal)
    call ply_poly_project_n2m(me = me,                  &
      &                       dim = 3 ,                 &
      &                       nVars = 1,                &
      &                       nodal_data=nodal_data,    &
      &                       modal_data= oversamp_modal)

    do iDegX = 1, maxdegree+1
      do iDegY = 1, maxdegree+1
        do iDegZ = 1, maxdegree+1
          dof = 1 + (iDegX-1) + (iDegY-1)*(maxdegree+1)+ &
             &    (iDegZ-1) * (maxdegree+1) * (maxdegree+1)
          dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(me%oversamp_degree+1)+ &
           &            (iDegZ-1) * (me%oversamp_degree+1)*&
           &            (me%oversamp_degree+1)
          modal_data(dof,:) = oversamp_modal(dofOverSamp,:)
        end do
      end do
    end do

    res= maxval( abs(ref_modes-modal_data(:,1) ) )

    write(*,*) 'power=', power
    write(*,*) 'res=', res
  end subroutine check_fpt

end program ply_project_fpt_test
