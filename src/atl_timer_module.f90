! Copyright (c) 2011-2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2014-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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

! *******************************************************************************
!> This module provides a convenient way to setup timers in the code.
!!
!! To define a new one, declare an integer module variable here, and
!! add a tem_timer_module::TEM_addTimer statement for it in the ATL_addTimers
!! routine with a label for identification. Then use this module in the relevant
!! code part and call tem_timer_module::TEM_startTimer and
!! tem_timer_module::TEM_stopTimer by passing in this integer.
!!
!! Upon finishing execution, Ateles will then print the measured time, for
!! this code block into timing.res.
module atl_timer_module
  use env_module,              only: PathLen, rk, long_k, &
    &                                my_status_int, newunit

  use tem_logging_module,      only: logunit
  use tem_timer_module,        only: tem_addTimer, tem_getTimerVal, &
    &                                tem_resetTimer, &
    &                                tem_getMaxTimerVal, tem_getTimerName
  use soi_revision_module,     only: soi_solver_revision
  use tem_general_module,      only: tem_general_type
  use treelmesh_module,        only: treelmesh_type
  use tem_timer_module,        only: tem_timer_type, tem_appendTimers


  implicit none

  private

  public :: atl_timerHandles
  public :: atl_elemTimers
  public :: atl_cpl_elemTimers
  public :: atl_addTimers
  public :: atl_resetTimers
  public :: atl_addTimers_perElem
  public :: atl_dumpTimers
  public :: atl_timer_handle_type
  public :: atl_get_timerHandles
  public :: atl_set_timerHandles
  public :: atl_get_elemTimers
  public :: atl_set_elemTimers
  public :: atl_get_cpl_elemTimers
  public :: atl_set_cpl_elemTimers

  !> Handles for timer objects to measure the time for some code parts
  type atl_timer_handle_type
    integer :: init
    integer :: preciceInit
    integer :: simLoop
    integer :: commState
    integer :: wRestart
    integer :: preprocessKernel
    integer :: projectToFace
    integer :: setBnd
    integer :: invMassMatrix
    integer :: numFlux
    integer :: physFlux
    integer :: pF_initStateDer
    integer :: pF_projConv
    integer :: pF_pen
    integer :: pF_eval
    integer :: pF_projectTestFunc
    integer :: projectTestFunc
    integer :: stabalize
    integer :: updateBackground
    integer :: get_timestep
    integer :: TimeStepInfo
    !integer :: syncUpdate
    integer :: convergeCheck
    integer :: checkVal
    integer :: preciceAdv
    integer :: readBC
    integer :: varElem
    integer :: constElem
    integer :: gradient
    integer :: preciceWrite
    !> Set first timer handle in ateles
    integer :: first = 0
    !> Set last timer handle in ateles
    integer :: last = -1
!VK!    integer :: precice
  end type atl_timer_handle_type

  type(atl_timer_handle_type), save :: atl_timerHandles
  ! element timer for multi-level
  type(tem_timer_type), save :: atl_elemTimers
  ! element timer for coupling
  type(tem_timer_type), save :: atl_cpl_elemTimers


contains


  ! ***************************************************************************
  !> Setup timers to assess the runtime of various parts of Ateles.
  subroutine atl_addTimers()
    !--------------------------------------------------------------------------

    ! Create some timer objects
    call tem_addTimer(timerName   = 'initialize',         &
      &               timerHandle = atl_timerHandles%init )

    call tem_addTimer(timerName   = 'preciceInit',         &
      &               timerHandle = atl_timerHandles%preciceinit )

    call tem_addTimer(timerName   = 'simLoop',               &
      &               timerHandle = atl_timerHandles%simLoop )

    call tem_addTimer(timerName   = 'commState',               &
      &               timerHandle = atl_timerHandles%commState )

    call tem_addTimer(timerName   = 'Output',                 &
      &               timerHandle = atl_timerHandles%wRestart )

    call tem_addTimer(timerName   = 'preprocKern',                    &
      &               timerHandle = atl_timerHandles%preprocessKernel )

    call tem_addTimer(timerName   = 'projToFace',                  &
      &               timerHandle = atl_timerHandles%projectToFace )

    call tem_addTimer(timerName   = 'setBnd',               &
      &               timerHandle = atl_timerHandles%setBnd )

    call tem_addTimer(timerName   = 'invMassMat',                  &
      &               timerHandle = atl_timerHandles%invMassMatrix )

    call tem_addTimer(timerName   = 'numFlux',               &
      &               timerHandle = atl_timerHandles%numFlux )

    call tem_addTimer(timerName   = 'physFlux',               &
      &               timerHandle = atl_timerHandles%physFlux )

    call tem_addTimer(timerName   = 'pF_initState',                  &
      &               timerHandle = atl_timerHandles%pF_initStateDer )

    call tem_addTimer(timerName   = 'pF_projConv',               &
      &               timerHandle = atl_timerHandles%pF_projConv )

    call tem_addTimer(timerName   = 'pF_pen',               &
      &               timerHandle = atl_timerHandles%pF_pen )

    call tem_addTimer(timerName   = 'pF_eval',               &
      &               timerHandle = atl_timerHandles%pF_eval )

    call tem_addTimer(timerName   = 'pF_projTestFunc',                  &
      &               timerHandle = atl_timerHandles%pF_projectTestFunc )

    call tem_addTimer(timerName   = 'projTestFun',                   &
      &               timerHandle = atl_timerHandles%projectTestFunc )

    call tem_addTimer(timerName   = 'stabalize',                      &
      &               timerHandle = atl_timerHandles%stabalize )

!VK not used at the moment
!VK    call tem_addTimer(timerName   = 'invMassT',                              &
!VK           &               timerHandle = atl_timerHandles%transposedInvMassMatrix )
!VK
!VK    call tem_addTimer(timerName   = 'localProj',                     &
!VK      &               timerHandle = atl_timerHandles%localProjection )

    call tem_addTimer(timerName   = 'updateBG',                       &
      &               timerHandle = atl_timerHandles%updateBackground )

    call tem_addTimer(timerName   = 'getTimestep',                &
      &               timerHandle = atl_timerHandles%get_timestep )

    call tem_addTimer(timerName   = 'TimeStepInfo',               &
      &               timerHandle = atl_timerHandles%TimeStepInfo )

   !! call tem_addTimer(timerName   = 'syncUpdate',               &
   !!   &               timerHandle = atl_timerHandles%syncUpdate )

    call tem_addTimer(timerName   = 'convCheck',                   &
      &               timerHandle = atl_timerHandles%convergeCheck )

    call tem_addTimer(timerName   = 'checkVal',               &
      &               timerHandle = atl_timerHandles%checkVal )

    call tem_addTimer(timerName   = 'preciceAdv',               &
      &               timerHandle = atl_timerHandles%preciceAdv )

    call tem_addTimer(timerName   = 'readBC',               &
      &               timerHandle = atl_timerHandles%readBC )

    call tem_addTimer(timerName   = 'constElem',               &
      &               timerHandle = atl_timerHandles%constElem )

    call tem_addTimer(timerName   = 'varElem',               &
      &               timerHandle = atl_timerHandles%varElem )

    call tem_addTimer(timerName   = 'gradient',               &
      &               timerHandle = atl_timerHandles%gradient )

    call tem_addTimer(timerName   = 'preciceWrite',               &
      &               timerHandle = atl_timerHandles%preciceWrite )

    ! Set first and last timer handle to get timer val in dump timing
    atl_timerHandles%first = atl_timerHandles%init
    atl_timerHandles%last = atl_timerHandles%preciceWrite
  end subroutine atl_addTimers
  ! ***************************************************************************


  ! ***************************************************************************
  subroutine atl_resetTimers()
    !--------------------------------------------------------------------------
    ! This routine helps for debugging and finding routins, that might 
    ! more or less time during computation. With that we can reset the timer 
    ! e.g. after each iteration
    ! Create some timer objects
    call tem_resetTimer(timerHandle = atl_timerHandles%init )

    call tem_resetTimer(timerHandle = atl_timerHandles%preciceinit )

    call tem_resetTimer(timerHandle = atl_timerHandles%simLoop )

    call tem_resetTimer(timerHandle = atl_timerHandles%commState )

    call tem_resetTimer(timerHandle = atl_timerHandles%wRestart )

    call tem_resetTimer(timerHandle = atl_timerHandles%preprocessKernel )

    call tem_resetTimer(timerHandle = atl_timerHandles%projectToFace )

    call tem_resetTimer(timerHandle = atl_timerHandles%setBnd )

    call tem_resetTimer(timerHandle = atl_timerHandles%invMassMatrix )

    call tem_resetTimer(timerHandle = atl_timerHandles%numFlux )

    call tem_resetTimer(timerHandle = atl_timerHandles%physFlux )

    call tem_resetTimer(timerHandle = atl_timerHandles%pF_initStateDer )

    call tem_resetTimer(timerHandle = atl_timerHandles%pF_projConv )

    call tem_resetTimer(timerHandle = atl_timerHandles%pF_pen )

    call tem_resetTimer(timerHandle = atl_timerHandles%pF_eval )

    call tem_resetTimer(timerHandle = atl_timerHandles%pF_projectTestFunc )

    call tem_resetTimer(timerHandle = atl_timerHandles%projectTestFunc )

    call tem_resetTimer(timerHandle = atl_timerHandles%stabalize )

    call tem_resetTimer(timerHandle = atl_timerHandles%updateBackground )

    call tem_resetTimer(timerHandle = atl_timerHandles%get_timestep )

    call tem_resetTimer(timerHandle = atl_timerHandles%TimeStepInfo )

    !call tem_resetTimer(timerHandle = atl_timerHandles%syncUpdate )

    call tem_resetTimer(timerHandle = atl_timerHandles%convergeCheck )

    call tem_resetTimer(timerHandle = atl_timerHandles%checkVal )

    call tem_resetTimer(timerHandle = atl_timerHandles%preciceAdv )

    call tem_resetTimer(timerHandle = atl_timerHandles%readBC )

    call tem_resetTimer(timerHandle = atl_timerHandles%preciceWrite )

  end subroutine atl_resetTimers
  ! ***************************************************************************

  ! ***************************************************************************
  !> Setup timers to assess the runtime of various parts of Ateles.
  subroutine atl_addTimers_perElem(tree)
    !--------------------------------------------------------------------------
    !> Mesh data in treelmesh format.
    type( treelmesh_type ), intent(in) :: tree
    !--------------------------------------------------------------------------
    ! add timer for get point routine for LB - elementwise
    write(logUnit(5),*) 'Add the elementwise timer for writing LB weights'

    ! at the moment we use only one element timer
    ! if adding more timer here you need to be carefull
    ! since the remaining code rely one the elemPos to acccess these timers
    call tem_appendTimers( me    = atl_elemTimers , &
        &                  nVals = tree%nElems      )

    call tem_appendTimers( me    = atl_cpl_elemTimers, &
        &                  nVals = tree%nElems         )


  end subroutine atl_addTimers_perElem
  ! ***************************************************************************


  ! ***************************************************************************!
  !> This function returns local modular variable atl_timerHandles to apesmate
  function atl_get_timerHandles() result(timerHandles)
    !---------------------------------------------------------------------------
    type(atl_timer_handle_type) :: timerHandles
    !---------------------------------------------------------------------------
    timerHandles = atl_timerHandles
  end function atl_get_timerHandles
  !****************************************************************************!


  ! ***************************************************************************!
  !> This function returns local modular variable atl_elemTimers to apesmate
  function atl_get_cpl_elemTimers() result(elemTimers)
    !---------------------------------------------------------------------------
    type(tem_timer_type) :: elemTimers
    !---------------------------------------------------------------------------
    elemTimers = atl_elemTimers
  end function atl_get_cpl_elemTimers
  !****************************************************************************!


  ! ***************************************************************************!
  !> This routine sets elementTimers passed by apesmate
  subroutine atl_set_cpl_elemTimers(elemTimers)
    !---------------------------------------------------------------------------
    type(tem_timer_type), intent(in) :: elemTimers
    !---------------------------------------------------------------------------
    atl_elemTimers = elemTimers
  end subroutine atl_set_cpl_elemTimers
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> This function returns local modular variable atl_elemTimers to apesmate
  function atl_get_elemTimers() result(elemTimers)
    !---------------------------------------------------------------------------
    type(tem_timer_type) :: elemTimers
    !---------------------------------------------------------------------------
    elemTimers = atl_elemTimers
  end function atl_get_elemTimers
  !****************************************************************************!


  ! ***************************************************************************!
  !> This routine sets atl_timerHandles passed by apesmate
  subroutine atl_set_timerHandles(timerHandles)
    !---------------------------------------------------------------------------
    type(atl_timer_handle_type), intent(in) :: timerHandles
    !---------------------------------------------------------------------------
    atl_timerHandles = timerHandles
  end subroutine atl_set_timerHandles
  ! ***************************************************************************!

  ! ***************************************************************************!
  !> This routine sets elementTimers passed by apesmate
  subroutine atl_set_elemTimers(elemTimers)
    !---------------------------------------------------------------------------
    type(tem_timer_type), intent(in) :: elemTimers
    !---------------------------------------------------------------------------
    atl_elemTimers = elemTimers
  end subroutine atl_set_elemTimers
  ! ***************************************************************************!

  ! ***************************************************************************
  !> Performance results are written to a file for statistical review
  !! The file-format is simple can be evaluated with gnuplot
  subroutine atl_dumptimers(general, nElems, nDofs, nVars)
    ! ---------------------------------------------------------------------------
    !> Parameters of the current simulation
    type(tem_general_type), intent(in) :: general
    integer(kind=long_k), intent(in) :: nElems
    !> The number of degrees of freedom per element per scalar variable
    integer, intent(in) :: nDofs
    !> The number of scalar variables
    integer, intent(in) :: nVars
    ! ---------------------------------------------------------------------------
    !< Memory usage (resident set size, high water mark)
    integer :: memRss, memHwm
    logical              :: file_exists
    character(len=pathLen)        :: filename
    integer              :: fileunit, iTimer
    real(kind=rk),allocatable        :: timerVal(:)
    character(len=40),allocatable    :: timerLabel(:)
    character(len=PathLen)    :: header
    character(len=PathLen)    :: output
    integer              :: iterations
    integer(kind=long_k) :: totaldofs
    real(kind=rk)        :: kweight
    integer :: nTimers, counter
    real(kind=rk) :: tAteles
    real(kind=rk) :: tSyncUpdate
    ! ---------------------------------------------------------------------------
    ! number of timers measured internally in Ateles
    nTimers = atl_timerHandles%last - atl_timerHandles%first + 1
    allocate( timerVal( nTimers ) )
    ! first and last handle convers all ateles handles which are added
    ! contigously
    counter = 0
    ! Get MaxTimer outside if isRoot since tem_getMaxTimerVal uses
    ! mpi_allreduce
    do iTimer = atl_timerHandles%first, atl_timerHandles%last
      counter = counter+1
      timerVal(counter)   = tem_getMaxTimerVal(timerHandle = iTimer,           &
        &                                      comm        = general%proc%comm )
    end do

    if (general%proc%isRoot) then
      !>@todo HK: Make mem-stuff configurable.
      !!          Maybe reduce values from all processes
      !!          to find global maximum.
      memRss = my_status_int('VmRSS:')
      memHwm = my_status_int('VmHWM:')

      write(header,'(a1,1x,a12,1x,a20,1x,a8,1x,a8,1x,a12,1x,a15,6(1x,a12), &
        &            a1,1x,a12,a1)')&
        & '#', 'Revision', &
        & 'Casename', &
        & 'nProcs',&     ! The number of proc (i.e. MPI ranks)
        & 'threads', &   ! The number of OMP threads per proc
        & 'DomSize', &   ! The number of elements
        & 'Dofs', &      ! Total dofs = DofPEPV * elements * nVars
        & 'DofPE', &     ! Dofs per element = DofPEPV * nVars = Dof / DomSize
        & 'DofPEPV', &   ! Dofs per element per scalar variable
        & 'nVars', &     ! The number of scalar variables in the equ. system
        & 'KEUPS',&      ! Thousand element updates per second
        & 'KDUPS', &     ! Thousand dof updates per second
        & 'maxIter', '|', &
        & 'timeAteles', '|' ! Total time taken for Ateles

      allocate( timerLabel( nTimers ) )

      ! first and last handle convers all ateles handles which are added
      ! contigously
      counter = 0
      do iTimer = atl_timerHandles%first, atl_timerHandles%last
        counter = counter+1
        timerLabel(counter) = trim(tem_getTimerName(timerHandle = iTimer))
        write(header,'(a,a13,a1)') trim(header), trim(timerLabel(counter)), '|'
      enddo
        write(header,'(a,a13,a1)') trim(header), 'tem_syncUpdate', '|'

      write(header,'(a,2(1x,a12))') trim(header), &
        & 'MemRSS', &    ! memory usage in sim loop
        & 'MemHWM'       ! memory usage max

      !>@todo HK: ensure, that timing is actually now, and it is valid to use
      !!          the iter component of it as the overall number of iterations
      !!          (might be different after restart?)
      iterations = general%simControl%now%iter &
        &        - general%simControl%timeControl%min%iter
      totaldofs = nDofs*nVars*nElems
      kweight = 1.0_rk/(tem_getTimerVal(timerHandle = atl_timerHandles%simLoop) &
        &               *1000.0_rk)
      ! total time taken for Ateles
      tAteles = tem_getTimerVal( timerHandle = general%solver%timerHandle )
      tSyncUpdate= tem_getTimerVal( timerHandle = general%simControl%syncUpdate_timer )
      write(output, '(1x,a13,1x,a20,1x,i8,1x,i8,1x,i12,1x,' &
        &           // 'i15,1x,i12,1x,i12,1x,i12,1x,en12.3,1x,en12.3,' &
        &           // '1x,i12,1x,en12.3)') &
        &   trim(soi_solver_revision), &
        &   trim(general%solver%simName), &
        &   general%proc%comm_size, &
        &   general%proc%nThreads, &
        &   nElems, &
        &   nDofs*nVars*nElems, & ! Total dofs = DofPEPV * elements * nVars
        &   nDofs*nVars, & ! Dofs per element = DofPEPV * nVars = Dof / DomSize
        &   nDofs, & ! Dofs per element per scalar variable
        &   nVars, & ! The number of scalar variables in the equ. system
        &   (nElems*iterations)*kweight, & ! Thousand element updates per second
        &   (totaldofs*iterations)*kweight, & ! Thousand dof updates per second
        &   iterations, &
        &   tAteles ! Time spend on Ateles

      do iTimer = 1, nTimers
        write(output,'(a,1x,en13.3)') trim(output), timerVal( iTimer )
      enddo
      write(output,'(a,1x,en13.3)') trim(output), tSyncUpdate

      write(output,'(a,i12,i12)') trim(output), memRss, memHwm

      filename = trim(general%timingFile)
      write(logunit(2),*) 'Writing timing information to ', trim(filename), '.'
      inquire(file=filename, exist=file_exists)
      fileunit = newunit()
      open(unit=fileunit, file=trim(filename), position='append')

      if (.not. file_exists ) then
         write(fileunit,'(a)') trim(header)
      end if
      write(fileunit,'(a)') trim(output)
      close(fileunit)
    end if


  end subroutine atl_dumptimers
! ****************************************************************************** !

end module atl_timer_module
! ******************************************************************************
