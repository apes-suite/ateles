! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2014-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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

! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt and Harald Klimach for
! University of Siegen.
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
!
! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the number of degrees of freedom for Q polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * A variable to store the number of degrees of freedom for a P tensor product
!   polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.

! The x ansatz degree is turned into the degree of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x ansatz degree is turned into the degree of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
!> author: Jens Zudrop
!! Module with all covolume projections.
module atl_covolume_projection_module
  use env_module,                    only: rk, zero_rk

  ! Treelm modules
  use tem_aux_module,                only: tem_abort
  use tem_logging_module,            only: logUnit

  ! Ateles modules
  use atl_scheme_module,             only: atl_scheme_type
  use atl_covolume_module,           only: atl_covolume_type
  use atl_spectral_viscosity_module, only: atl_poly_spectral_visc_prp

  implicit none
  private

  public :: atl_primal_to_covolume_projection
  public :: atl_covolume_to_primal_projection
  public :: atl_primal_to_covolume_projection_2d
  public :: atl_covolume_to_primal_projection_2d
  public :: atl_primal_to_covolume_projection_1d
  public :: atl_covolume_to_primal_projection_1d

contains

  !> summary: Project two elements onto single co-volume element.
  !!
  !! This routine projects to elements (left and right) onto its co-volume
  !! element. The geometrical setup is as follows:
  !!
  !!          left                 right
  !! |---------------------||---------------------|
  !!                 \               /
  !!                  \             /
  !!                  _\|         |/_
  !!             |---------------------|
  !!                    co-volume
  !!
  !! The transformation to the co-volume is carried out by a simple (but efficient)
  !! L2-projection.
  subroutine atl_primal_to_covolume_projection( left, right, dir, filter,     &
    &                                           scheme, maxPolyDeg, covolume, &
    &                                           order                         )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: left(:,:)
    real(kind=rk), intent(in) :: right(:,:)
    integer, intent(in) :: dir
    type(atl_covolume_type), intent(in) :: filter
    !> The numerical schemes for the current level to get the modg basis
    type(atl_scheme_type), intent(in) :: scheme
    integer, intent(in) :: maxPolyDeg
    real(kind=rk), intent(out) :: covolume(:,:)
    real(kind=rk), intent(in) :: order
    ! --------------------------------------------------------------------------
    integer :: iDegX, iDegY, iDegZ, dofX, dofY, dofZ, iHelpVar, iter
    integer :: dofAbsX, dofAbsY, dofAbsZ
    real(kind=rk) :: dofAbs, damping, cut
    integer :: mpd1, pos, pos_primal
    logical :: use_damping
    ! --------------------------------------------------------------------------

    use_damping = (filter%alpha > zero_rk)

    if (filter%cut_order <= 0.0_rk) then
      cut = real(maxPolyDeg,rk)
    else
      cut = filter%cut_order
    end if

    mpd1 = maxPolyDeg+1


    do iter=lbound(covolume,2),ubound(covolume,2)
      covolume(:,iter) = 0.0_rk
    end do

    select case(dir)
    ! co-volume projection in x-direction
    case(1)
      ! for each y-z-mode we perform x-covol. proj.
      do pos=1,mpd1**3
        iDegZ = (pos-1)/(mpd1)**2 + 1
        iHelpVar = pos - (iDegZ-1)*(mpd1)**2
        iDegY = (iHelpVar-1)/(mpd1)+1
        dofX = mod(pos-1, mpd1) +1

        do iDegX=1,mpd1
  pos_primal = idegx                                      &
    &      + ( ( idegy-1)                             &
    &      + (idegz-1)*(maxpolydeg+1))*(maxpolydeg+1)
          covolume(pos,:) = covolume(pos,:)                                  &
          & + left(pos_primal,:)                                             &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dofX,iDegX,1) &
          & + right(pos_primal,:)                                            &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dofX,iDegX,2)
        end do

        if (use_damping) then
          if (filter%kind == atl_poly_spectral_visc_prp) then

            ! Modal polynomial damping factor
            dofAbsX = int(min((real(dofX,rk)  - 1.0_rk) / cut, 1.0_rk))
            dofAbsY = int(min((real(iDegY,rk) - 1.0_rk) / cut, 1.0_rk))
            dofAbsZ = int(min((real(iDegZ,rk) - 1.0_rk) / cut, 1.0_rk))
            damping =   (1.0_rk - (dofAbsX**order)) &
                    & * (1.0_rk - (dofAbsY**order)) &
                    & * (1.0_rk - (dofAbsZ**order))

          else
            ! Modal exponential damping factor
            dofAbs = sqrt( &
                           & ( &
                           &      ((real(dofX,rk)-1.0_rk)**2) &
                           &    + ((real(iDegY,rk)-1.0_rk)**2) &
                           &    + ((real(iDegZ,rk)-1.0_rk)**2) &
                           & ) &
                           &  / (cut**2) &
                           &  )
            damping = exp( -filter%alpha * (real(dofAbs,rk)**order) )

          end if

          ! Apply modal damping factor
          covolume(pos, :) = damping * covolume(pos, :)

        end if

      end do !pos

    ! co-volume projection in y-direction
    case(2)
      ! for each x-z-mode we perform y-covol. proj.
      do pos=1,mpd1**3
        iDegZ = (pos-1)/(mpd1)**2 + 1
        iHelpVar = pos - (iDegZ-1)*(mpd1)**2
        dofY = (iHelpVar-1)/(mpd1)+1
        iDegX = mod(pos-1, mpd1) +1

        do iDegY=1,mpd1
  pos_primal = idegx                                      &
    &      + ( ( idegy-1)                             &
    &      + (idegz-1)*(maxpolydeg+1))*(maxpolydeg+1)
          covolume(pos, :) = covolume(pos, :)                                &
          & + left(pos_primal,:)                                             &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dofY,iDegY,1) &
          & + right(pos_primal,:)                                            &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dofY,iDegY,2)
        end do

        if (use_damping) then
          if (filter%kind .eq. atl_poly_spectral_visc_prp) then

            ! Modal polynomial damping factor
            dofAbsX = int(min((real(iDegX,rk) - 1.0_rk) / cut, 1.0_rk))
            dofAbsY = int(min((real(dofY,rk)  - 1.0_rk) / cut, 1.0_rk))
            dofAbsZ = int(min((real(iDegZ,rk) - 1.0_rk) / cut, 1.0_rk))
            damping =   (1.0_rk - (dofAbsX**order)) &
                    & * (1.0_rk - (dofAbsY**order)) &
                    & * (1.0_rk - (dofAbsZ**order))

          else
            ! Modal exponential damping factor
            dofAbs = sqrt( &
                           & ( &
                           &      ((real(iDegX,rk)-1.0_rk)**2) &
                           &    + ((real(dofY,rk)-1.0_rk)**2) &
                           &    + ((real(iDegZ,rk)-1.0_rk)**2) &
                           & ) &
                           &  / (cut**2) &
                           &  )
            damping = exp( -filter%alpha * (real(dofAbs,rk)**order) )
          end if

          ! Apply modal damping factor
          covolume(pos, :) = damping * covolume(pos, :)

        end if

     end do ! pos

    ! co-volume projection in z-direction
    case(3)
      ! for each x-y-mode we perform z-covol. proj.
      do pos=1,mpd1**3
        dofZ = (pos-1)/(mpd1)**2 + 1
        iHelpVar = pos - (dofZ-1)*(mpd1)**2
        iDegY = (iHelpVar-1)/(mpd1)+1
        iDegX = mod(pos-1, mpd1) +1

        do iDegZ=1,mpd1
  pos_primal = idegx                                      &
    &      + ( ( idegy-1)                             &
    &      + (idegz-1)*(maxpolydeg+1))*(maxpolydeg+1)
          covolume(pos, :) = covolume(pos, :)                                &
          & + left(pos_primal,:)                                             &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dofZ,iDegZ,1) &
          & + right(pos_primal,:)                                            &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dofZ,iDegZ,2)
        end do

        if (use_damping) then
          if(filter%kind .eq. atl_poly_spectral_visc_prp) then

            ! Modal polynomial damping factor
            dofAbsX = int(min((real(iDegX,rk) - 1.0_rk) / cut, 1.0_rk))
            dofAbsY = int(min((real(iDegY,rk) - 1.0_rk) / cut, 1.0_rk))
            dofAbsZ = int(min((real(dofZ,rk)  - 1.0_rk) / cut, 1.0_rk))
            damping =   (1.0_rk - (dofAbsX**order)) &
                    & * (1.0_rk - (dofAbsY**order)) &
                    & * (1.0_rk - (dofAbsZ**order))

          else

            ! Modal exponential damping factor
            dofAbs = sqrt( &
                           & ( &
                           &      ((real(iDegX,rk)-1.0_rk)**2) &
                           &    + ((real(iDegY,rk)-1.0_rk)**2) &
                           &    + ((real(dofZ,rk)-1.0_rk)**2) &
                           & ) &
                           &  / (cut**2) &
                           &  )
            damping = exp( -filter%alpha * (real(dofAbs,rk)**order) )

          end if

          ! Apply modal damping factor
          covolume(pos, :) = damping * covolume(pos, :)

        end if

     end do ! iter

    case default
      write(logUnit(1),*) 'ERROR in atl_primal_to_covolume_projection:'
      write(logUnit(1),*) 'Unknown spatial co-volume direction, stopping ...'
      call tem_abort()
    end select

  end subroutine atl_primal_to_covolume_projection


  !> summary: Project two co-volume elements onto single a single element.
  !!
  !! This routine projects two co-volume elements (left and right) onto its
  !! primal element. The geometrical setup is as follows:
  !!
  !!          left (co-vol.)       right (co-vol.)
  !! |---------------------||---------------------|
  !!                 \               /
  !!                  \             /
  !!                  _\|         |/_
  !!             |---------------------|
  !!                    primal element
  !!
  !! The transformation to the primal element is carried out by a simple (but efficient)
  !! L2-projection.
  function atl_covolume_to_primal_projection( left, right, dir, filter,     &
    &                                         scheme, maxPolyDeg, nScalars, &
    &                                         state )                       &
    &                                         result( primal )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: left(:,:)
    real(kind=rk), intent(in) :: right(:,:)
    integer, intent(in) :: dir
    type(atl_covolume_type), intent(in) :: filter
    !> The numerical schemes for the current level to get the modg basis
    type(atl_scheme_type), intent(in) :: scheme
    integer, intent(in) :: maxPolyDeg
    integer, intent(in) :: nScalars
    real(kind=rk), intent(in) :: state(:,:)
    real(kind=rk) :: primal((maxPolyDeg+1)**3,nScalars)
    ! --------------------------------------------------------------------------
    integer :: iDegX, iDegY, iDegZ, iHelpVar, dof_covolume
    integer :: mpd1
    integer :: pos_primal, pos_covol
    ! --------------------------------------------------------------------------

    mpd1 = maxPolyDeg+1


    ! Weighted original state.
    primal =  (1.0_rk - filter%beta)*state(:,:)

    select case(dir)

    ! Co-volume projection in x direction
    case(1)

      ! Project back to primary grid
      do pos_primal =1,mpd1**3
        iDegZ = (pos_primal-1)/(mpd1)**2 + 1
        iHelpVar = pos_primal - (iDegZ-1)*(mpd1)**2
        iDegY = (iHelpVar-1)/(mpd1)+1
        iDegX = mod(pos_primal-1, mpd1) +1

        covolXLoop: do dof_covolume = 1, mpd1
  pos_covol = dof_covolume                                      &
    &      + ( ( idegy-1)                             &
    &      + (idegz-1)*(maxpolydeg+1))*(maxpolydeg+1)
          primal(pos_primal,:) = primal( pos_primal,:)             &
          & + left(pos_covol,:)                                    &
          & * filter%beta                                          &
          & * scheme%modg_basis%covolumeBaseCoeff                  &
          &                    %anz_anzShift(iDegX,dof_covolume,1) &
          & + right(pos_covol,:)                                   &
          & * filter%beta                                          &
          & * scheme%modg_basis%covolumeBaseCoeff                  &
          &                    %anz_anzShift(iDegX,dof_covolume,2)
        end do covolXLoop

      end do !iter

    ! Co-volume projection in y direction
    case(2)

      ! Project back to primary grid
      do pos_primal =1,mpd1**3
        iDegZ = (pos_primal-1)/(mpd1)**2 + 1
        iHelpVar = pos_primal - (iDegZ-1)*(mpd1)**2
        iDegY = (iHelpVar-1)/(mpd1)+1
        iDegX = mod(pos_primal-1, mpd1) +1

        covolYLoop: do dof_covolume = 1, mpd1
  pos_covol = idegx                                      &
    &      + ( ( dof_covolume-1)                             &
    &      + (idegz-1)*(maxpolydeg+1))*(maxpolydeg+1)
          primal(pos_primal,:) = primal( pos_primal,:)             &
          & + left(pos_covol,:)                                    &
          & * filter%beta                                          &
          & * scheme%modg_basis%covolumeBaseCoeff                  &
          &                    %anz_anzShift(iDegY,dof_covolume,1) &
          & + right(pos_covol,:)                                   &
          & * filter%beta                                          &
          & * scheme%modg_basis%covolumeBaseCoeff                  &
          &                    %anz_anzShift(iDegY,dof_covolume,2)
        end do covolYLoop

      end do !iter

    ! Co-volume projection in z direction
    case(3)

      ! Project back to primary grid
      do pos_primal =1,mpd1**3
        iDegZ = (pos_primal-1)/(mpd1)**2 + 1
        iHelpVar = pos_primal - (iDegZ-1)*(mpd1)**2
        iDegY = (iHelpVar-1)/(mpd1)+1
        iDegX = mod(pos_primal-1, mpd1) +1

        covolZloop: do dof_covolume = 1, mpd1
  pos_covol = idegx                                      &
    &      + ( ( idegy-1)                             &
    &      + (dof_covolume-1)*(maxpolydeg+1))*(maxpolydeg+1)
          primal(pos_primal,:) = primal( pos_primal,:)             &
          & + left(pos_covol,:)                                    &
          & * filter%beta                                          &
          & * scheme%modg_basis%covolumeBaseCoeff                  &
          &                    %anz_anzShift(iDegZ,dof_covolume,1) &
          & + right(pos_covol,:)                                   &
          & * filter%beta                                          &
          & * scheme%modg_basis%covolumeBaseCoeff                  &
          &                    %anz_anzShift(iDegZ,dof_covolume,2)
        end do covolZloop

      end do !iter

    case default
      write(logUnit(1),*) 'ERROR in atl_covolume_to_primal_projection:'
      write(logUnit(1),*) 'Unknown spatial co-volume direction, stopping ...'
      call tem_abort()
    end select


  end function atl_covolume_to_primal_projection


  !> summary: Project two elements onto single co-volume element.
  !!
  !! This routine projects to elements (left and right) onto its co-volume
  !! element. The geometrical setup is as follows:
  !!
  !!          left                 right
  !! |---------------------||---------------------|
  !!                 \               /
  !!                  \             /
  !!                  _\|         |/_
  !!             |---------------------|
  !!                    co-volume
  !!
  !! The transformation to the co-volume is carried out by a simple (but efficient)
  !! L2-projection.
  subroutine atl_primal_to_covolume_projection_2d( left, right, dir, filter, &
    &                                              scheme, maxPolyDeg,       &
    &                                              covolume, order           )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: left(:,:)
    real(kind=rk), intent(in) :: right(:,:)
    integer, intent(in) :: dir
    type(atl_covolume_type), intent(in) :: filter
    !> The numerical schemes for the current level to get the modg basis
    type(atl_scheme_type), intent(in) :: scheme
    integer, intent(in) :: maxPolyDeg
    real(kind=rk), intent(out) :: covolume(:,:)
    real(kind=rk), intent(in) :: order
    ! --------------------------------------------------------------------------
    integer :: iDegX, iDegY, dof, iter
    real(kind=rk) :: dofAbs, damping, cut, dofAbsX, dofAbsY
    integer :: mpd1, pos, pos_primal
    logical :: use_damping
    ! --------------------------------------------------------------------------

    use_damping = (filter%alpha > zero_rk)

    if (filter%cut_order <= 0.0_rk) then
      cut = real(maxPolyDeg,rk)
    else
      cut = filter%cut_order
    end if

    mpd1 = maxPolyDeg+1


    do iter=lbound(covolume,2),ubound(covolume,2)
    covolume(:,iter) = 0.0_rk
      end do

    select case(dir)

    ! co-volume projection in x-direction
    case(1)

    ! PrimalYLoop and covolXLoop collapsed
    ! for each y-mode we perform x-covol. proj.
    do iter=1, mpd1**2
      idegY= ((iter-1)/mpd1)+1
      dof=mod(iter-1,mpd1)+1
  pos = dof                                      &
    &      + ( ( idegy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
      do iDegX=1,mpd1
  pos_primal = idegx                                      &
    &      + ( ( idegy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
        covolume(pos, :) = covolume(pos, :)                                 &
          & + left(pos_primal,:)                                            &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dof,iDegX,1) &
          & + right(pos_primal,:)                                           &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dof,iDegX,2)
      end do

      if (use_damping) then
        if (filter%kind == atl_poly_spectral_visc_prp) then
          ! Modal polynomial damping factor
          dofAbsX = min((real(dof,rk)   - 1.0_rk) / cut, 1.0_rk)
          dofAbsY = min((real(iDegY,rk) - 1.0_rk) / cut, 1.0_rk)
          damping =  (1.0_rk -  (dofAbsX**order)) * (1.0_rk - (dofAbsY**order))
        else
          ! Modal exponential damping factor
          dofAbs = sqrt((((real(dof,rk) - 1.0_rk)**2) &
            & + ((real(iDegY,rk)-1.0_rk)**2)) / (cut**2) )
          damping = exp( (-filter%alpha * (dofAbs**order) ) )
        end if
        ! Apply modal damping factor
        covolume(pos, :) = damping * covolume(pos, :)
      end if
    end do

    ! co-volume projection in y-direction
    case(2)

    ! PrimalXLoop and covolYLoop collapsed
    ! for each x-mode we perform y-covol. proj.
    do iter=1, mpd1**2
      idegX= ((iter-1)/mpd1)+1
      dof=mod(iter-1,mpd1)+1
  pos = idegx                                      &
    &      + ( ( dof-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
      do iDegY=1,mpd1
  pos_primal = idegx                                      &
    &      + ( ( idegy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
        covolume(pos, :) = covolume(pos, :)                                 &
          & + left(pos_primal,:)                                            &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dof,iDegY,1) &
          & + right(pos_primal,:)                                           &
          & * scheme%modg_basis%covolumeBaseCoeff%anz_anzShift(dof,iDegY,2)
      end do

      if (use_damping) then
        if (filter%kind == atl_poly_spectral_visc_prp) then
          ! Modal polynomial damping factor
          dofAbsX = min((real(iDegX,rk) - 1.0_rk) / cut, 1.0_rk)
          dofAbsY = min((real(dof,rk)   - 1.0_rk) / cut, 1.0_rk)
          damping = (1.0_rk - (dofAbsX**order)) * (1.0_rk - (dofAbsY**order))
        else
          ! Modal exponential damping factor
          dofAbs = sqrt((((real(iDegX,rk) - 1.0_rk)**2) &
            & + ((real(dof,rk)-1.0_rk)**2)) / (cut**2) )
          damping = exp( (-filter%alpha * (dofAbs**order) ) )
        end if

        ! Apply modal damping factor
        covolume(pos, :) = damping * covolume(pos, :)
      end if
    end do

    case default
      write(logUnit(1),*) 'ERROR in atl_primal_to_covolume_projection_2d:'
      write(logUnit(1),*) 'Unknown spatial co-volume direction, stopping ...'
      call tem_abort()
    end select


  end subroutine atl_primal_to_covolume_projection_2d



  !> summary: Project two co-volume elements onto single a single element.
  !!
  !! This routine projects two co-volume elements (left and right) onto its
  !! primal element. The geometrical setup is as follows:
  !!
  !!          left (co-vol.)       right (co-vol.)
  !! |---------------------||---------------------|
  !!                 \               /
  !!                  \             /
  !!                  _\|         |/_
  !!             |---------------------|
  !!                    primal element
  !!
  !! The transformation to the primal element is carried out by a simple (but efficient)
  !! L2-projection.
  function atl_covolume_to_primal_projection_2d( left, right, dir, filter,     &
    &                                            scheme, maxPolyDeg, nScalars, &
    &                                            state )                       &
    &                                            result( primal )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: left(:,:)
    real(kind=rk), intent(in) :: right(:,:)
    integer, intent(in) :: dir
    type(atl_covolume_type), intent(in) :: filter
    !> The numerical schemes for the current level to get the modg basis
    type(atl_scheme_type), intent(in) :: scheme
    integer, intent(in) :: maxPolyDeg
    integer, intent(in) :: nScalars
    real(kind=rk), intent(in) :: state(:,:)
    real(kind=rk) :: primal((maxPolyDeg+1)**2,nScalars)
    real(kind=rk) :: prim_loc(nScalars)
    ! --------------------------------------------------------------------------
    integer :: iDegX, iDegY, dof_covolume, iter, iHelpVar
    integer :: mpd1
    integer :: pos_primal, pos_covol
    ! --------------------------------------------------------------------------

    mpd1 = maxPolyDeg+1


    select case(dir)

    ! Co-volume projection in x direction
    case(1)

      do iter=1,mpd1**2
        iDegY = (iter-1)/mpd1 + 1
        iHelpVar = iter - (iDegY-1)*mpd1**2
        iDegX = mod(iter-1,mpd1) + 1
  pos_primal = idegx                                      &
    &      + ( ( idegy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)

        ! Weighted original state.
        prim_loc =  (1.0_rk - filter%beta)*state(pos_primal,:)
        do dof_covolume=1,mpd1
          ! Project back to primary grid
  pos_covol = dof_covolume                                      &
    &      + ( ( idegy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
          prim_loc = prim_loc                                              &
            &      + left(pos_covol,:) * filter%beta                       &
            &                          * scheme%modg_basis                 &
            &                                  %covolumeBaseCoeff          &
            &                                  %anz_anzShift(iDegX,        &
            &                                                dof_covolume, &
            &                                                1)            &
            &      + right(pos_covol,:)* filter%beta                       &
            &                          * scheme%modg_basis                 &
            &                                  %covolumeBaseCoeff          &
            &                                  %anz_anzShift(iDegX,        &
            &                                                dof_covolume, &
            &                                                2)
        end do
        primal(pos_primal,:) = prim_loc

      end do

    ! Co-volume projection in y direction
    case(2)

      do iter=1,mpd1**2
        iDegY = (iter-1)/mpd1 + 1
        iHelpVar = iter - (iDegY-1)*mpd1**2
        iDegX = mod(iter-1,mpd1) + 1
  pos_primal = idegx                                      &
    &      + ( ( idegy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)

        ! Weighted original state.
        prim_loc =  (1.0_rk - filter%beta)*state(pos_primal,:)
        do dof_covolume=1,mpd1
          ! Project back to primary grid
  pos_covol = idegx                                      &
    &      + ( ( dof_covolume-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
          prim_loc = prim_loc                                              &
            &      + left(pos_covol,:) * filter%beta                       &
            &                          * scheme%modg_basis                 &
            &                                  %covolumeBaseCoeff          &
            &                                  %anz_anzShift(iDegY,        &
            &                                                dof_covolume, &
            &                                                1)            &
            &      + right(pos_covol,:)* filter%beta                       &
            &                          * scheme%modg_basis                 &
            &                                  %covolumeBaseCoeff          &
            &                                  %anz_anzShift(iDegY,        &
            &                                                dof_covolume, &
            &                                                2)
        end do
        primal(pos_primal,:) = prim_loc

      end do

    case default
      write(logUnit(1),*) 'ERROR in atl_covolume_to_primal_projection_2d:'
      write(logUnit(1),*) 'Unknown spatial co-volume direction, stopping ...'
      call tem_abort()
    end select


  end function atl_covolume_to_primal_projection_2d



  !> summary: Project two elements onto single co-volume element.
  !!
  !! This routine projects to elements (left and right) onto its co-volume
  !! element. The geometrical setup is as follows:
  !!
  !!          left                 right
  !! |---------------------||---------------------|
  !!                 \               /
  !!                  \             /
  !!                  _\|         |/_
  !!             |---------------------|
  !!                    co-volume
  !!
  !! The transformation to the co-volume is carried out by a simple (but efficient)
  !! L2-projection.
  subroutine atl_primal_to_covolume_projection_1d( left, right, filter, &
    &                                              scheme, maxPolyDeg,  &
    &                                              covolume             )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: left(:,:)
    real(kind=rk), intent(in) :: right(:,:)
    type(atl_covolume_type), intent(in) :: filter
    !> The numerical schemes for the current level to get the modg basis
    type(atl_scheme_type), intent(in) :: scheme
    integer, intent(in) :: maxPolyDeg
    real(kind=rk), intent(out) :: covolume(:,:)
    ! --------------------------------------------------------------------------
    integer :: dof, dof_primary
    real(kind=rk) :: dofAbs, damping, cut
    integer :: mpd1
    logical :: use_damping
    ! --------------------------------------------------------------------------

    use_damping = (filter%alpha > zero_rk)

    if (filter%cut_order <= 0.0_rk) then
      cut = real(maxPolyDeg,rk)
    else
      cut = filter%cut_order
    end if

    mpd1 = maxPolyDeg+1


    do dof=lbound(covolume,2),ubound(covolume,2)
      covolume(:,dof) = 0.0_rk
    end do

    do dof=1,mpd1
      do dof_primary = 1, mpd1
        covolume(dof, :) = covolume(dof, :)                   &
        & + left(dof_primary,:)                               &
        & * scheme%modg_basis%covolumeBaseCoeff               &
        &                    %anz_anzShift(dof,dof_primary,1) &
        & + right(dof_primary,:)                              &
        & * scheme%modg_basis%covolumeBaseCoeff               &
        &                    %anz_anzShift(dof,dof_primary,2)
      end do

      if (use_damping) then
        if (filter%kind == atl_poly_spectral_visc_prp) then
          ! Modal polynomial damping factor
          dofAbs = min((real(dof,rk) - 1.0_rk) / cut, 1.0_rk)
          damping = 1.0_rk - dofAbs**filter%order
        else
          ! Modal exponential damping factor
          dofAbs = sqrt( ((real(dof,rk)-1.0_rk)**2) / (cut**2) )
          damping = exp( (-filter%alpha * (dofAbs**filter%order) ) )
        end if
        ! Apply modal damping factor
        covolume(dof, :) = damping * covolume(dof, :)
      end if

    end do


  end subroutine atl_primal_to_covolume_projection_1d

  !> summary: Project two co-volume elements onto single a single element.
  !!
  !! This routine projects two co-volume elements (left and right) onto its
  !! primal element. The geometrical setup is as follows:
  !!
  !!          left (co-vol.)       right (co-vol.)
  !! |---------------------||---------------------|
  !!                 \               /
  !!                  \             /
  !!                  _\|         |/_
  !!             |---------------------|
  !!                    primal element
  !!
  !! The transformation to the primal element is carried out by a simple (but efficient)
  !! L2-projection.
  function atl_covolume_to_primal_projection_1d( left, right, filter, scheme,  &
    &                                            maxPolyDeg, nScalars, state ) &
    &                                            result( primal )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: left(:,:)
    real(kind=rk), intent(in) :: right(:,:)
    type(atl_covolume_type), intent(in) :: filter
    !> The numerical schemes for the current level to get the modg basis
    type(atl_scheme_type), intent(in) :: scheme
    integer, intent(in) :: maxPolyDeg
    integer, intent(in) :: nScalars
    real(kind=rk), intent(in) :: state(:,:)
    real(kind=rk) :: primal(maxPolyDeg+1,nScalars)
    ! --------------------------------------------------------------------------
    integer :: dof, dof_covolume
    integer :: mpd1
    ! --------------------------------------------------------------------------

    mpd1 = maxPolyDeg+1


    ! Project back to primary grid
    primal = (1.0_rk - filter%beta)*state(:,:)

    do dof = 1, mpd1
      do dof_covolume = 1, mpd1
        primal(dof,:) = primal( dof,:)                           &
          & + left(dof_covolume,:)                               &
          & * filter%beta                                        &
          & * scheme%modg_basis%covolumeBaseCoeff                &
          &                    %anz_anzShift(dof,dof_covolume,1) &
          & + right(dof_covolume,:)                              &
          & * filter%beta                                        &
          & * scheme%modg_basis%covolumeBaseCoeff                &
          &                    %anz_anzShift(dof,dof_covolume,2)
      end do
    end do


  end function atl_covolume_to_primal_projection_1d

end module atl_covolume_projection_module
