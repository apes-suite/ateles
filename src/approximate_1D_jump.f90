! Copyright (c) 2015 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2016, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach and Peter Vitt
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

!> This is a small utility to approximate a discontinuity at a given
!! location in the interval [-1,1].
!!
!! The jump is located at a given point jota, and the target function is
!! 0 for x <= jota and 1 for x > jota.
!! The interval is subdivided by level bisections towards jota.
!! Intervals for which the lower bound is smaller than jota are assumed to
!! be 0, all others are 1.
!! This interval values are used to determine the values at Chebyshev nodes
!! and a FPT is done to find the Legendre modes.
!!
!! Required settings are therefore:
!! - jota (location of the jump)
!! - level (number of bisections towards jota)
!! - target polynomial degree
!! - oversampling factor
!!
!! These have to be provided as command line parameters in this order.
!!
!! Resulting output is:
!! - Interval subdivision
!! - Approximated jump location
!! - Chebyshev nodes
!! - Legendre modes
!! - L2-Error
program approximate_1D_jump
  use env_module,                   only: rk

  use tem_dyn_array_module,         only: dyn_realArray_type, init, append

  use ply_dof_module,               only: q_space
  use ply_dynArray_project_module,  only: ply_prj_init_type
  use ply_modg_basis_module,        only: ply_modg_refine_type,           &
    &                                     ply_init_modg_multilevelCoeffs, &
    &                                     ply_integrateLeg,               &
    &                                     ply_legendre_1D,                &
    &                                     ply_scalprodleg
  use ply_poly_project_module,      only: ply_poly_project_type, &
    &                                     ply_poly_project_fillbody, &
    &                                     ply_poly_project_n2m

  implicit none

  character(len=32) :: argstring
  real(kind=rk) :: jota
  integer :: level
  integer :: polydegree

  integer :: ibis
  integer :: ival
  integer :: ipoint
  integer :: imode
  integer :: iCoarse, iFine
  integer :: intmode
  integer :: bis_lb
  integer :: current, last
  integer :: lastside
  integer :: bpos
  real(kind=rk) :: interval_jump
  real(kind=rk) :: ofact
  real(kind=rk) :: bis_half
  real(kind=rk) :: x_point
  real(kind=rk) :: cur_lb
  real(kind=rk) :: ivldistance
  real(kind=rk) :: l2err
  real(kind=rk) :: optierr
  real(kind=rk) :: voldiff
  real(kind=rk), allocatable :: bisect(:)
  real(kind=rk), allocatable :: nodal_data(:,:)
  real(kind=rk), allocatable :: modal_data(:,:)
  real(kind=rk), allocatable :: legmodes(:)
  real(kind=rk), allocatable :: exact(:)
  real(kind=rk), allocatable :: const1_left(:)
  real(kind=rk), allocatable :: const1_right(:)
  real(kind=rk), allocatable :: integral(:)
  real(kind=rk), allocatable :: intatjota(:,:)
  integer, allocatable :: ivcolor(:)
  type(ply_prj_init_type) :: project
  type(ply_poly_project_type) :: polypro
  type(ply_modg_refine_type) :: refine
  type(dyn_realArray_type) :: ivlen
  character(len=10) :: method

  call get_command_argument(1,method)

  call get_command_argument(2,argstring)
  read(argstring,*) jota
  call get_command_argument(3,argstring)
  read(argstring,*) level
  call get_command_argument(4,argstring)
  read(argstring,*) polydegree

  if (trim(method) == 'fast') then
    call get_command_argument(5,argstring)
    read(argstring,*) ofact
  else
    ofact = 1
  end if

  if ( (jota <= -1.0_rk) .or. (jota >= 1.0_rk) ) then
    write(*,*) 'First parameter (jota) needs to be a real value in'
    write(*,*) 'the interval (-1,1), but is:', jota
    write(*,*) 'Stopping!'
    STOP
  end if

  if (level < 0) then
    write(*,*) 'Second parameter (level) needs to be non-negative!'
    write(*,*) 'But it is:', level
    write(*,*) 'Stopping!'
    STOP
  end if

  if (polydegree < 0) then
    write(*,*) 'Third parameter (polydegree) needs to be non-negative!'
    write(*,*) 'But it is:', polydegree
    write(*,*) 'Stopping!'
    STOP
  end if

  if (trim(method) == 'fast') then
    if (ofact < 0.0_rk) then
      write(*,*) 'Fourth parameter (oversampling) needs to be non-negative!'
      write(*,*) 'But it is:', ofact
      write(*,*) 'Stopping!'
      STOP
    end if
  end if

  ! Initialize current to a value outside it's normal range to check whether it
  ! was set during runtime.
  current = -1

  ! Voxelization (in 1D just a bisection algorithm...)
  allocate(bisect(0:level+1))

  bisect = 1.0_rk
  bisect(0) = -1.0_rk
  bis_lb = 0
  do ibis=1,level

    bis_half = 0.5_rk * ( bisect(bis_lb+1) - bisect(bis_lb) )
    do ival=level,bis_lb+1,-1
      bisect(ival+1) = bisect(ival)
    end do
    bisect(bis_lb+1) = bisect(bis_lb) + bis_half
    if (bisect(bis_lb+1) < jota) then
      bis_lb = bis_lb+1
    end if

  end do

  call init( me     = ivlen,   &
    &        length = level+1  )
  allocate(ivcolor(level+1))
  do ibis=1,level+1
    ivldistance = bisect(ibis)-bisect(ibis-1)
    call append( me  = ivlen,       &
      &          val = ivldistance, &
      &          pos = bpos         )
    if (bisect(bpos-1) < jota) then
      ivcolor(bpos) = 0
    else
      ivcolor(bpos) = 1
    end if
  end do

  write(*,*) 'Bisection points:'
  do ibis=0,level+1
    write(*,*) bisect(ibis)
  end do
  write(*,*) '-----------------'

  ival = 0
  do
    if (bisect(iVal) > jota) EXIT
    iVal = iVal + 1
  end do
  interval_jump = bisect(ival)

  write(*,*) 'Approximated jump location:', interval_jump

  if (trim(method) == 'fast') then
    project%basisType = q_space
    project%maxPolyDegree = polydegree
    project%header%kind = 'fpt'
    project%header%fpt_header%factor = ofact
    project%header%fpt_header%approx_terms = 18
    project%header%fpt_header%striplen = 128
    project%header%fpt_header%adapt_factor_pow2 = .false.
    project%header%fpt_header%nodes_header%lobattopoints = .false.
    project%header%fpt_header%nodes_header%nodes_kind = 'chebyshev'
    call ply_poly_project_fillbody( me         = polypro, &
      &                             proj_init  = project, &
      &                             scheme_dim = 1        )
  end if


  if (trim(method) == 'iterative') then

    call ply_init_modg_multilevelCoeffs( nPoints  = 128,          &
      &                                  nFunc    = polydegree+1, &
      &                                  integral = refine        )

    ! Projection of constant 1 to the half of a coarser element.
    ! (Precomputation)
    allocate(const1_left(polydegree+1))
    allocate(const1_right(polydegree+1))
    do iCoarse=1,min(polydegree+1,2)
      !! sqnorm = 2.0_rk / (2.0_rk*iCoarse - 1.0_rk)
      !! jacobidetfinetocoarse = 0.5
      !! const = anz_anzShift * jacobidetfinetocoarse / sqnorm
      const1_left(iCoarse) = 0.25_rk * (2.0_rk*iCoarse-1) &
        &                            * refine%anz_anzShift(1,iCoarse,1)
      const1_right(iCoarse) = 0.25_rk * (2.0_rk*iCoarse-1) &
        &                             * refine%anz_anzShift(1,iCoarse,2)
    end do

    allocate(modal_data(polydegree+1,0:1))

    modal_data = 0.0_rk

    ! set the current lower bound to the lower bound of the smallest interval.
    cur_lb = bisect(ivlen%sorted(1)-1)

    ! Fill the coarser element according to the status of the two small
    ! intervals.
    if (ivcolor(ivlen%sorted(1)) == 1) then
      modal_data(:,0) = const1_left
    end if
    if (ivcolor(ivlen%sorted(2)) == 1) then
      modal_data(:,0) = modal_data(:,0) + const1_right
    end if

    do ibis=3,level+1
      current = mod(ibis,2)
      last = mod(ibis-1,2)
      if ( bisect(ivlen%sorted(ibis)-1) < cur_lb) then
        ! New constant part is left.
        lastside = 2
        if (ivcolor(ivlen%sorted(ibis)) == 1) then
          modal_data(:,current) = const1_left
        else
          modal_data(:,current) = 0.0_rk
        end if
        cur_lb = bisect(ivlen%sorted(ibis)-1)
      else
        ! New constant part is right.
        lastside = 1
        if (ivcolor(ivlen%sorted(ibis)) == 1) then
          modal_data(:,current) = const1_right
        else
          modal_data(:,current) = 0.0_rk
        end if
      end if
      do iCoarse=1,polydegree+1
        !! sqnorm = 2.0_rk / (2.0_rk*iCoarse - 1.0_rk)
        !! jacobidetfinetocoarse = 0.5
        !! const = anz_anzShift * jacobidetfinetocoarse / sqnorm
        do iFine=1,polydegree+1
          modal_data(iCoarse,current) = modal_data(iCoarse,current)      &
            &                  + 0.25_rk * (2.0_rk*iCoarse-1)            &
            &                            * refine%anz_anzShift(iFine,    &
            &                                                  iCoarse,  &
            &                                                  lastside) &
            &                            * modal_data(iFine,last)
        end do
      end do
    end do

  end if


  if (trim(method) == 'fast') then
    current = 1
    allocate(nodal_data(polypro%nQuadPointsPerDir,1))
    allocate(modal_data(polypro%nQuadPointsPerDir,1))

    write(*,*)
    write(*,*) 'Chebyshev points:'
    do ipoint=1,polypro%nQuadPointsPerDir
      x_point = polypro%body_1D%nodes(ipoint,1)
      write(*,*) x_point
      if (x_point >= interval_jump) then
        nodal_data(ipoint,1) = 1.0_rk
      else
        nodal_data(ipoint,1) = 0.0_rk
      end if
    end do
    write(*,*) '-----------------'

    call ply_poly_project_n2m( me         = polypro,    &
      &                        dim        = 1,          &
      &                        nVars      = 1,          &
      &                        nodal_data = nodal_data, &
      &                        modal_data = modal_data  )

  end if

  allocate(legmodes(polydegree+1))
  allocate(exact(polydegree+1))
  allocate(integral(polydegree+2))
  allocate(intatjota(polydegree+2,1))

  do iMode=1,polydegree+1
    legmodes = 0.0_rk
    legmodes(iMode) = 1.0_rk
    integral(:iMode+1) = ply_integrateleg( integrand = legmodes(:iMode), &
      &                                    maxdegree = iMode             )
    intatjota = ply_legendre_1D( points = [jota], degree = polydegree+1 )
    exact(iMode) = 0.0_rk
    do intmode=1,iMode+1
      exact(iMode) = exact(iMode) &
        &          + integral(intmode)*(1.0_rk - intatjota(intmode,1))
    end do
    exact(iMode) = exact(iMode) / ply_scalprodleg(iMode)
  end do

  if (current == -1) then
    write(*,*) "Variable current is used unintialized. Results may be invalid."
  endif

  l2err = 0.0_rk
  optierr = 0.0_rk
  write(*,*) 'Legendre Modes:'
  do iMode=1,polydegree+1
    write(*,*) modal_data(iMode,current), exact(iMode)
    optierr = optierr + exact(iMode)**2 * ply_scalprodleg(iMode)
    l2err = l2err + ply_scalprodleg(iMode)*(modal_data(iMode,current) &
      &                                 - exact(iMode))**2
  end do
  write(*,*) '-----------------'

  voldiff = abs(2*modal_data(1,current) - (1.0_rk-jota))

  write(*,*) 'N=', 2**level
  write(*,*) 'L2 Error numeric vs optimal: ', sqrt(l2err)
  write(*,*) 'Projection Error:', sqrt(1-jota - optierr)
  write(*,*) 'Vol Error: ', voldiff


end program approximate_1D_jump
