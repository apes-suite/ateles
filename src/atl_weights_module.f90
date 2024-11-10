! Copyright (c) 2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
! Copyright (c) 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> This module include the routine required for element wie dumping weight
!! for better load balancing
!! This is intialize when add 'write_weights = 'weight' ' to the config file
!! and are based on element wise time measurements
module atl_weights_module
  use mpi
  use env_module,              only: PathLen, rk, tem_create_EndianSuffix
  use aotus_module,            only: flu_State, aot_get_val

  use tem_logging_module,      only: logunit
  use treelmesh_module,        only: treelmesh_type, tem_dump_weights
  use tem_timer_module,        only: tem_getTimerVal

  use atl_timer_module,        only: atl_timerHandles, atl_cpl_elemTimers, &
    &                                atl_elemTimers, atl_addTimers_perElem

  implicit none

  private

  public :: atl_initWeights
  public :: atl_dumpWeights


contains


  ! ************************************************************************** !
  subroutine atl_initWeights( tree )
    ! ---------------------------------------------------------------------- !
    !> Handle to the configuration script, to load the parameters from.
    !> Mesh data in treelmesh format.
    type( treelmesh_type ), intent(in) :: tree
    ! ---------------------------------------------------------------------- !

    ! init the weight timer per element
    ! @TODO VK:
    ! this need to be done despite weights are written or not
    ! since the timer will be started and stoped in getStateForPoint anyway
    call atl_addTimers_perElem(tree)

  end subroutine atl_initWeights
  ! ************************************************************************** !


  ! ************************************************************************** !
  subroutine atl_getWeights( tree, weights, comp_weights, multiLevel_weights, &
      &                      cpl_weights                                      )
    ! ---------------------------------------------------------------------- !
    type(treelmesh_type), intent(in) :: tree
    real(kind=rk), intent(inout), allocatable :: weights(:)
    real(kind=rk), intent(inout), allocatable :: multiLevel_weights(:)
    real(kind=rk), intent(inout), allocatable :: cpl_weights(:)
    real(kind=rk), intent(inout), allocatable :: comp_weights(:)
    ! ---------------------------------------------------------------------- !
    real(kind=rk) :: numFlux, physFlux, projTestFun
    real(kind=rk) :: projtoFace, invMassMat, compute_perElem
    real(kind=rk) :: preprocKern, checkVal, updateBG, stabalize
    integer :: iElem
    ! ---------------------------------------------------------------------- !

    ! get weights of the computation
    preprocKern = tem_getTimerVal( timerHandle = atl_timerHandles &
      &                                          %preprocessKernel)
    numFlux     = tem_getTimerVal( timerHandle = atl_timerHandles%numFlux )
    physFlux    = tem_getTimerVal( timerHandle = atl_timerHandles%physFlux )
    projTestFun = tem_getTimerVal(timerHandle = atl_timerHandles%projectTestFunc)
    projtoFace = tem_getTimerVal( timerHandle = atl_timerHandles%projectToFace)
    invMassMat = tem_getTimerVal( timerHandle = atl_timerHandles%invMassMatrix)
    checkVal   = tem_getTimerVal( timerHandle = atl_timerHandles%checkVal)
    updateBG   = tem_getTimerVal( timerHandle =                     &
      &                           atl_timerHandles%updateBackground )
    stabalize = tem_getTimerVal( timerHandle = atl_timerHandles%stabalize)
    !> boundary timer should only be measured for boundary elements
!!VK    setBnd = tem_getTimerVal( timerHandle = atl_timerHandles%setBnd)
!!VK    readBC = tem_getTimerVal( timerHandle = atl_timerHandles%readBC)

    compute_perElem = (                                                 &
      &                 numFlux + projTestFun   &
      &                 + projtoFace + invMassMat + checkVal + updateBG &
      &                 + stabalize                                     &
      &               )                                                 &
      &               / tree%nElems

    ! get the weight due to coupling and calling getPnt routine
    ! and add the compute part per element
    allocate ( weights(tree%nElems) )
    allocate ( multiLevel_weights(tree%nElems) )
    allocate ( cpl_weights(tree%nElems) )
    allocate ( comp_weights(tree%nElems) )
    weights = 0.0_rk
    do iElem = 1, tree%nElems
      comp_weights(iElem) = compute_perElem
      multiLevel_weights(iElem) = tem_getTimerVal(atl_ElemTimers, iElem)
      cpl_weights(iElem) = tem_getTimerVal(atl_cpl_ElemTimers, iElem)
      weights(iElem) = compute_perElem &
        &              + cpl_weights(iElem) + multilevel_weights(iElem)
    end do

  end subroutine atl_getWeights
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Dump weights to a file.
  subroutine atl_dumpWeights(tree)
    ! ---------------------------------------------------------------------- !
    type(treelmesh_type), intent(in) :: tree
    ! ---------------------------------------------------------------------- !
    character(len=pathLen) :: weights_file
    real(kind=rk),allocatable :: weights(:)
    ! weights due to pure computation
    real(kind=rk),allocatable :: comp_weights(:)
    ! weights due to multi level
    real(kind=rk),allocatable :: multiLevel_weights(:)
    ! weights due to coupling
    real(kind=rk),allocatable :: cpl_weights(:)
    character(len=PathLen) :: filename
    character(len=4) :: EndianSuffix
    ! ---------------------------------------------------------------------- !

    EndianSuffix = tem_create_EndianSuffix()
    call atl_getWeights(tree, weights, comp_weights, multiLevel_weights, &
      &                 cpl_weights)

    ! write weights file
    weights_file = trim(tree%write_weights)
    filename = trim(weights_file)//EndianSuffix
    call tem_dump_weights( filename = filename, &
      &                    me       = tree,     &
      &                    weights  = weights   )
    deallocate(weights)

    ! write comp weights file
    filename = trim(weights_file)//'_compute'//EndianSuffix
    call tem_dump_weights( filename = filename,    &
      &                     me       = tree,        &
      &                     weights  = comp_weights )
    deallocate(comp_weights)

    ! write muliLevel weights file
    EndianSuffix = tem_create_EndianSuffix()
    filename = trim(weights_file)//'_multiLevel'//EndianSuffix
    call tem_dump_weights( filename = filename,          &
      &                     me       = tree,              &
      &                     weights  = multilevel_weights )
    deallocate(multiLevel_weights)

    ! write coupling weights file
    EndianSuffix = tem_create_EndianSuffix()
    filename = trim(weights_file)//'_coupling'//EndianSuffix
    call tem_dump_weights( filename = filename,   &
      &                     me       = tree,       &
      &                     weights  = cpl_weights )
    deallocate(cpl_weights)

  end subroutine atl_dumpWeights
  ! ************************************************************************** !


end module atl_weights_module
! ****************************************************************************** !
