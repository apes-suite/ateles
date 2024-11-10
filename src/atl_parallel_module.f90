! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2013, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2015, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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

!> author: Jens Zudrop
!! module to keep all routines and data types related to parallel
!! execution.
module atl_parallel_module
  use tem_comm_module,        only: tem_commPattern_type, &
                                  & tem_communication_type
  use tem_construction_module, only: tem_levelDesc_type
  use tem_faceData_module,    only: tem_face_type

  use atl_scheme_module,      only: atl_scheme_type
  use atl_cube_elem_module,   only: atl_cube_elem_type
  use atl_boundary_module,    only: atl_level_boundary_type

  implicit none

  private

  public :: atl_init_parallel_module
  public :: atl_init_cellstatebuffer, atl_init_faceStateBuffer


  type intarray_type
    integer, allocatable :: val(:)
  end type

contains

  !> initialize the parallel module to make it usable in ATELES
  subroutine atl_init_parallel_module( scheme, nValsElem, nValsStateFace,   &
    &                                  nValsFluxFace, cube, boundary,       &
    &                                  createCellBuffer, createFaceBuffer,  &
    &                                  createStabFaceBuffer,                &
    &                                  createStabElemBuffer, nBndStabElems, &
    &                                  minLevel, maxLevel, commPattern      )
    ! --------------------------------------------------------------------------
    !> The minimal refinement level of your mesh.
    integer, intent(in)                 :: minLevel
    !> The maximum refinement level of your mesh.
    integer, intent(in)                 :: maxLevel
    !> levelwise list of schemes.
    type(atl_scheme_type),intent(in)    :: scheme(minLevel:maxLevel)
    !> list of cubic meshes you want to build the buffers for.
    type(atl_cube_elem_type), intent(inout) :: cube(minLevel:maxLevel)
    !> The boundary description for the faces on the current level.
    type(atl_level_boundary_type), intent(in) :: boundary(minLevel:maxLevel)
    !> the number of scalar values for each element.
    integer, intent(in)                 :: nValsElem
    !> the number of scalar values on the face for the state
    integer, intent(in) :: nValsStateFace
    !> the number of scalar values on the face for the flux
    integer, intent(in) :: nValsFluxFace
    !> Boolean to indicate if cell state buffers a required or not.
    logical, intent(in)                 :: createCellBuffer
    !> Boolean to indicate if face buffers a required or not.
    logical, intent(in)                 :: createFaceBuffer
    !> Boolean to indicate if face buffers for stabilization are required or not.
    logical, intent(in)                 :: createStabFaceBuffer
    !> Boolean to indicate if elem buffers for stabilization are required or not.
    logical, intent(in)                 :: createStabElemBuffer
    !> The number of boundary elements for the stabilization element buffer
    integer, intent(in)                 :: nBndStabElems(minLevel:maxLevel,1:3)
    !> mpi communication pattern type
    type(tem_commPattern_type), intent(in) :: commPattern
    ! --------------------------------------------------------------------------
    integer :: level, iDir
    ! --------------------------------------------------------------------------

    ! init the buffers for each level
    do level = minLevel , maxLevel

      ! we init the buffers to transfer the cell states inside the stencils
      if( createCellBuffer ) then
        call atl_init_cellStateBuffer( scheme      = scheme(Level),          &
          &                            nVars       = nValsElem,              &
          &                            levelDesc   = cube(level)%descriptor, &
          &                            nBndElems   = 0,                      &
          &                            commPattern = commPattern             )
      end if

      if( createStabElemBuffer ) then
        do iDir = 1, 3
          call atl_init_cellStateBuffer(                               &
            & scheme      = scheme(Level),                             &
            & nVars       = nValsElem,                                 &
            & levelDesc   = cube(level)%faces_stab%dimByDimDesc(iDir), &
            & nBndElems   = nBndStabElems(level,iDir),                 &
            & commPattern = commPattern                                )
        end do
      end if

      if ( createStabFaceBuffer ) then
        call atl_init_faceStateBuffer( nFaceDofs   = scheme(level)%nDofs,    &
          &                            faces       = cube(level)%faces_stab, &
          &                            nValsState  = nValsElem,              &
          &                            nValsFlux   = nValsElem,              &
          &                            boundary    = boundary(level),        &
          &                            commPattern = commPattern             )
      end if

      ! we init the buffers to transfer the face values for the flux calculation
      if( createFaceBuffer ) then
        call atl_init_faceStateBuffer( nFaceDofs   = scheme(level)%nFaceDofs, &
          &                            faces       = cube(level)%faces,       &
          &                            nValsState  = nValsStateFace,          &
          &                            nValsFlux   = nValsFluxFace,           &
          &                            boundary    = boundary(level),         &
          &                            commPattern = commPattern              )
      end if
    end do

  end subroutine atl_init_parallel_module

  !> summary: initializes the face buffers for communication.
  subroutine atl_init_faceStateBuffer( nFaceDofs, faces, nValsState,     &
    &                                  nValsFlux, boundary, commPattern, &
    &                                  initRealBuf                       )
    ! --------------------------------------------------------------------------
    !> The number of degrees of freedoms per scalar variable per face.
    integer,intent(in)        :: nFaceDofs
    !> the number of scalar values on the face for the state
    integer, intent(in) :: nValsState
    !> the number of scalar values on the face for the flux
    integer, intent(in) :: nValsFlux
    !> list of cubic meshes you want to build the buffers for.
    type(tem_face_type), intent(inout) :: faces
    !> The boundary description for the faces on the current level.
    type(atl_level_boundary_type), intent(in) :: boundary
    !> mpi communication pattern type
    type(tem_commPattern_type), intent(in) :: commPattern
    !> Init real buffer (default), if not the integer buffer is initialized.
    logical, optional :: initRealBuf
    ! -----------------------------------------------------------------------
    integer :: iDir, iFace, maxBufSizeState, maxBufSizeFlux, nTotalFaces
    integer :: iProc, elemPosSize, iElem, iVar, iDof, nDofs, posIndex
    integer :: maxvalrecvState, maxvalsendState
    integer :: maxvalrecvFlux, maxvalsendFlux
    integer, allocatable :: bufEIDState(:), bufEIDFlux(:)
    logical :: realBuf
    ! -----------------------------------------------------------------------

    realBuf = .true.
    if(present(initRealBuf)) then
      realBuf = initRealBuf
    end if

    ! Loop over all the directions and left and right faces and init all the buffers
    ! for these parameters.
    directionLoop: do iDir = 1,3

      ! The number of faces in the current direction
      nTotalFaces = faces%dimByDimDesc(iDir)%nElems &
                   & + sum( boundary%bnd(:)%faces(iDir,1)%facePos%nVals ) &
                   & + sum( boundary%bnd(:)%faces(iDir,2)%facePos%nVals )

      faceLoop: do iFace = 1, 2

        maxValRecvState = 0
        maxValSendState = 0
        maxValRecvFlux = 0
        maxValSendFlux = 0

        ! init the face buffer itself
        call allocate_comm_buffer( &
          &        comm_type = faces%faces(iDir)%sendbuffer_state(iFace), &
          &        maxSize = maxValSendState, &
          &        isReal = realBuf )
        call allocate_comm_buffer( &
          &        comm_type = faces%faces(iDir)%sendbuffer_flux(iFace), &
          &        maxSize = maxValSendFlux, &
          &        isReal = realBuf )
        call allocate_comm_buffer( &
          &        comm_type = faces%faces(iDir)%recvbuffer_state(iFace), &
          &        maxSize = maxValRecvState, &
          &        isReal = realBuf )
        call allocate_comm_buffer( &
          &        comm_type = faces%faces(iDir)%recvbuffer_flux(iFace), &
          &        maxSize = maxValRecvFlux, &
          &        isReal = realBuf )


        ! Determine the size of the buffers
        nDofs = nFaceDofs
        ! ... for the state
        maxBufSizeState = max(maxvalRecvState, maxvalSendState) * nValsState * nDofs
        if(allocated(bufEIDState)) then
          deallocate(bufEIDState)
        end if
        allocate(bufEIDState(maxBufSizeState))
        ! ... for the flux
        maxBufSizeFlux = max(maxvalRecvFlux, maxvalSendFlux) * nValsFlux * nDofs
        if(allocated(bufEIDFlux)) then
          deallocate(bufEIDFlux)
        end if
        allocate(bufEIDFlux(maxBufSizeFlux))

        ! for sending the state
        do iProc = 1, faces%faces(iDir)%sendbuffer_state(iFace)%nProcs
          !>... cell states
          if( allocated( faces%faces(iDir)             &
            &                 %sendbuffer_state(iFace) &
            &                 %elemPos(iProc)%val )    ) then

            elemPosSize = faces%faces(iDir)             &
              &                %sendbuffer_state(iFace) &
              &                %elemPos(iProc)          &
              &                %nvals

          else
            !>@todo HK: can we not cancel out 0-sized messages?
            !! (is this part ever executed?)
            elemPosSize = 0
          end if
          posIndex = 0
          do iElem = 1, elemPosSize
            do iVar = 1, nValsState
              do iDof = 1, nDofs
                ! ... positions of cell state
                posIndex = posIndex + 1
                bufEIDState(posIndex) = 1                         &
                  & + (faces%faces(iDir)                          &
                  &         %sendbuffer_state(iFace)              &
                  &         %elempos(iProc)                       &
                  &         %val(ielem)                           &
                  &   - 1 )                                       &
                  & + (iDof-1) * nTotalFaces                      &
                  & + (iVar-1) * nTotalFaces * nDofs              &
                  & + (iFace-1) * nTotalFaces * nDofs * nValsState
              end do
            end do
          end do
          if(realBuf) then
            call commPattern%initBuf_real(            &
              & me    = faces%faces(iDir)             &
              &              %sendbuffer_state(iFace) &
              &              %buf_real(iProc),        &
              & pos   = bufEIDState,                  &
              & nVals = posIndex                      )
          else
            call commPattern%initBuf_int(             &
              & me    = faces%faces(iDir)             &
              &              %sendbuffer_state(iFace) &
              &              %buf_int(iProc),         &
              & pos   = bufEIDState,                  &
              & nVals = posIndex                      )
          end if
        end do
        ! for sending the flux
        do iProc = 1, faces%faces(iDir)%sendbuffer_flux(iFace)%nProcs
          if( allocated( faces%faces(iDir)            &
            &                 %sendbuffer_flux(iFace) &
            &                 %elemPos(iProc)         &
            &                 %val)                   ) then

            elemPosSize = faces%faces(iDir)            &
              &                %sendbuffer_flux(iFace) &
              &                %elemPos(iProc)         &
              &                %nvals

          else
            !>@todo HK: can we not cancel out 0-sized messages?
            !! (is this part ever executed?)
            elemPosSize = 0
          end if
          posIndex = 0
          do iElem = 1, elemPosSize
            do iVar = 1, nValsFlux
              do iDof = 1, nDofs
                ! ... positions of cell state
                posIndex = posIndex + 1
                bufEIDFlux(posIndex) = 1                          &
                  & + ( faces%faces(iDir)                         &
                  &          %sendbuffer_flux(iFace)              &
                  &          %elempos(iProc)                      &
                  &          %val(ielem)                          &
                  &   - 1 )                                       &
                  & + (iDof-1) * nTotalFaces                      &
                  & + (iVar-1) * nTotalFaces * nDofs              &
                  & + (iFace-1) * nTotalFaces * nDofs * nValsFlux
              end do
            end do
          end do
          if(realBuf) then
            call commPattern%initBuf_real(           &
              & me    = faces%faces(iDir)            &
              &              %sendbuffer_flux(iFace) &
              &              %buf_real(iProc),       &
              & pos   = bufEIDFlux,                  &
              & nVals = posIndex                     )
          else
            call commPattern%initBuf_int(            &
              & me    = faces%faces(iDir)            &
              &              %sendbuffer_flux(iFace) &
              &              %buf_int(iProc),        &
              & pos   = bufEIDFlux,                  &
              & nVals = posIndex                     )
          end if
        end do

        ! for receiving the state
        do iProc = 1, faces%faces(iDir)%recvbuffer_state(iFace)%nProcs
          if( allocated( faces%faces(iDir)             &
            &                 %recvbuffer_state(iFace) &
            &                 %elemPos(iProc)%val )    ) then

            elemPosSize = faces%faces(iDir)             &
              &                %recvbuffer_state(iFace) &
              &                %elemPos(iProc)          &
              &                %nvals

          else
            !>@todo HK: can we not cancel out 0-sized messages?
            !! (is this part ever executed?)
            elemPosSize = 0
          end if
          posIndex = 0
          do iElem = 1, elemPosSize
            do iVar = 1, nValsState
              do iDof = 1, nDofs
                posIndex = posIndex + 1
                bufEIDState(posIndex) = 1                         &
                  & + ( faces%faces(iDir)                         &
                  &          %recvbuffer_state(iFace)             &
                  &          %elempos(iProc)                      &
                  &          %val(ielem)                          &
                  &   - 1 )                                       &
                  & + (iDof-1) * nTotalFaces                      &
                  & + (iVar-1) * nTotalFaces * nDofs              &
                  & + (iFace-1) * nTotalFaces * nDofs * nValsState
              end do
            end do
          end do
          if(realBuf) then
            call commPattern%initBuf_real(            &
              & me    = faces%faces(iDir)             &
              &              %recvbuffer_state(iFace) &
              &              %buf_real(iProc),        &
              & pos   = bufEIDState,                  &
              & nVals = posIndex                      )
          else
            call commPattern%initBuf_int(             &
              & me    = faces%faces(iDir)             &
              &              %recvbuffer_state(iFace) &
              &              %buf_int(iProc),         &
              & pos   = bufEIDState,                  &
              & nVals = posIndex                      )
          end if
        end do
        ! for receiving the flux
        do iProc = 1, faces%faces(iDir)%recvbuffer_flux(iFace)%nProcs
          if( allocated( faces%faces(iDir)            &
            &                 %recvbuffer_flux(iFace) &
            &                 %elemPos(iProc)         &
            &                 %val )                  ) then
            elemPosSize = faces%faces(iDir)            &
              &                %recvbuffer_flux(iFace) &
              &                %elemPos(iProc)         &
              &                %nvals
          else
            !>@todo HK: can we not cancel out 0-sized messages?
            !! (is this part ever executed?)
            elemPosSize = 0
          end if
          posIndex = 0
          do iElem = 1, elemPosSize
            do iVar = 1, nValsFlux
              do iDof = 1, nDofs
                posIndex = posIndex + 1
                bufEIDFlux(posIndex) = 1                          &
                  & + ( faces%faces(iDir)                         &
                  &          %recvbuffer_flux(iFace)              &
                  &          %elempos(iProc)                      &
                  &          %val(ielem)                          &
                  &   - 1 )                                       &
                  & + (iDof-1) * nTotalFaces                      &
                  & + (iVar-1) * nTotalFaces * nDofs              &
                  & + (iFace-1) * nTotalFaces * nDofs * nValsFlux
              end do
            end do
          end do
          if(realBuf) then
            call commPattern%initBuf_real(           &
              & me    = faces%faces(iDir)            &
              &              %recvbuffer_flux(iFace) &
              &              %buf_real(iProc),       &
              & pos   = bufEIDFlux,                  &
              & nVals = posIndex                     )
          else
            call commPattern%initBuf_int(            &
              & me    = faces%faces(iDir)            &
              &              %recvbuffer_flux(iFace) &
              &              %buf_int(iProc),        &
              & pos   = bufEIDFlux,                  &
              & nVals = posIndex                     )
          end if
        end do

      end do faceLoop
    end do directionLoop

  contains
    subroutine allocate_comm_buffer( comm_type, isReal, maxSize )
      type( tem_communication_type ), intent(inout) :: comm_type
      logical, intent(in) :: isReal
      integer, intent(out) :: maxSize

      ! ... allocate memory for the receive buffer of the flux
      !@todo: instead of allocating explicitly here,
      !       one should call tem_comm_init to allocate buf_real or buf_int.
      if(comm_type%nProcs > 0) then
        maxSize = maxval( comm_type%elemPos(:)%nVals)
        if ( isReal ) then
          allocate( comm_type%Buf_real( comm_type%nProcs) )
        else
          allocate( comm_type%Buf_int(  comm_type%nProcs) )
        end if
        if (.not. allocated( comm_type%rqHandle ) ) then
          allocate( comm_type%rqHandle(comm_type%nProcs) )
        end if
      else
        maxSize = 0
      end if

    end subroutine allocate_comm_buffer

  end subroutine atl_init_faceStateBuffer

  !> Initialize the parallel module to make it usable in ATELES.
  subroutine atl_init_cellStateBuffer( scheme, nVars, levelDesc, nBndElems, &
    &                                  commPattern                          )
    ! --------------------------------------------------------------------------
    !> the the data of the kernel on the current level.
    type(atl_scheme_type),intent(in)      :: scheme
    !> the number of variables in our equation.
    integer, intent(in)                 :: nVars
    !> the buffer for the cell state transfer you want to be initialized.
    type(tem_levelDesc_type), intent(inout) :: levelDesc
    !> The number of boundary elements.
    integer, intent(in) :: nBndElems
    !> mpi communication pattern type
    type(tem_commPattern_type), intent(in) :: commPattern
    ! --------------------------------------------------------------------------
    integer :: iProc, elemPosSize, iElem, iVar, iDof, nDofs, posIndex
    integer, allocatable :: bufEID(:)
    integer :: maxBufSize
    ! the number of total elements (including fluid, ghost, halo and bnd cells).
    integer :: nTotalElems
    ! --------------------------------------------------------------------------

    nTotalElems = size(levelDesc%total) + nBndElems

    nDofs = scheme%nDofs
    maxBufSize = 0

    maxBufSize = max(maxval(levelDesc%sendbuffer%elemPos(:)%nVals), &
      &              maxval(levelDesc%recvbuffer%elemPos(:)%nVals), &
      &              maxval(levelDesc%sendbufferFromFiner%elemPos(:)%nVals), &
      &              maxval(levelDesc%recvbufferFromFiner%elemPos(:)%nVals), &
      &              maxval(levelDesc%sendbufferFromCoarser%elemPos(:)%nVals), &
      &              maxval(levelDesc%recvbufferFromCoarser%elemPos(:)%nVals)) &
      &         * nVars*nDofs
    allocate(bufEID(maxBufSize))

    ! for sending
    do iProc = 1, levelDesc%sendbuffer%nProcs
      !>... cell states
      if(allocated(levelDesc%sendbuffer%elemPos(iProc)%val)) then
        elemPosSize = levelDesc%sendbuffer%elemPos(iProc)%nvals
      else
        elemPosSize = 0
      end if
      posIndex = 0
      do iElem = 1, elemPosSize
        do iVar = 1, nVars
          do iDof = 1, nDofs
            ! ... positions of cell state
            posIndex = posIndex + 1
            bufEID(posIndex) = nTotalElems * nDofs * (iVar-1) &
              &              + nTotalElems * (iDof-1) &
              &              + levelDesc%sendbuffer%elempos(iProc)%val(ielem)
          end do
        end do
      end do
      call commPattern%initBuf_real( &
        &    me = levelDesc%sendbuffer%buf_real(iProc), &
        &    pos = bufEID, &
        &    nVals = posIndex)
    end do
    ! for sending (from finer)
    do iProc = 1, levelDesc%sendbufferFromFiner%nProcs
      !>... cell states
      if(allocated(levelDesc%sendbufferFromFiner%elemPos(iProc)%val)) then
        elemPosSize = levelDesc%sendbufferFromFiner%elemPos(iProc)%nvals
      else
        elemPosSize = 0
      end if
      posIndex = 0
      do iElem = 1, elemPosSize
        do iVar = 1, nVars
          do iDof = 1, nDofs
            ! ... positions of cell state
            posIndex = posIndex + 1
            bufEID(posIndex) = nTotalElems * nDofs * (iVar-1) &
              & + nTotalElems * (iDof-1) &
              & + levelDesc%sendbufferFromFiner%elempos(iProc)%val(ielem)
          end do
        end do
      end do
      call commPattern%initBuf_real( &
        &    me = levelDesc%sendbufferFromFiner%buf_real(iProc), &
        &    pos = bufEID, &
        &    nVals = posIndex)
    end do
    ! for sending (from coarser)
    do iProc = 1, levelDesc%sendbufferFromCoarser%nProcs
      !>... cell states
      if(allocated(levelDesc%sendbufferFromCoarser%elemPos(iProc)%val)) then
        elemPosSize = levelDesc%sendbufferFromCoarser%elemPos(iProc)%nvals
      else
        elemPosSize = 0
      end if
      posIndex = 0
      do iElem = 1, elemPosSize
        do iVar = 1, nVars
          do iDof = 1, nDofs
            ! ... positions of cell state
            posIndex = posIndex + 1
            bufEID(posIndex) = nTotalElems * nDofs * (iVar-1) &
              & + nTotalElems * (iDof-1) &
              & + levelDesc%sendbufferFromCoarser%elempos(iProc)%val(ielem)
          end do
        end do
      end do
      call commPattern%initBuf_real( &
        &    me = levelDesc%sendbufferFromCoarser%buf_real(iProc), &
        &    pos = bufEID, &
        &    nVals = posIndex)
    end do

    ! for receiving
    do iProc = 1, levelDesc%recvbuffer%nProcs
      if(allocated(levelDesc%recvbuffer%elemPos(iProc)%val)) then
        elemPosSize = levelDesc%recvbuffer%elemPos(iProc)%nvals
      else
        elemPosSize = 0
      end if
      posIndex = 0
      do iElem = 1, elemPosSize
        do iVar = 1, nVars
          do iDof = 1, nDofs
            posIndex = posIndex + 1
            bufEID(posIndex) = nTotalElems * nDofs * (iVar-1) &
              &              + nTotalElems * (iDof-1) &
              &              + levelDesc%recvbuffer%elempos(iProc)%val(ielem)
          end do
        end do
      end do
      call commPattern%initBuf_real( &
        &    me = levelDesc%recvbuffer%buf_real(iProc), &
        &    pos = bufEID, &
        &    nVals = posIndex)
    end do
    ! for receiving (from finer)
    do iProc = 1, levelDesc%recvbufferFromFiner%nProcs
      if(allocated(levelDesc%recvbufferFromFiner%elemPos(iProc)%val)) then
        elemPosSize = levelDesc%recvbufferFromFiner%elemPos(iProc)%nvals
      else
        elemPosSize = 0
      end if
      posIndex = 0
      do iElem = 1, elemPosSize
        do iVar = 1, nVars
          do iDof = 1, nDofs
            posIndex = posIndex + 1
            bufEID(posIndex) = nTotalElems * nDofs * (iVar-1) &
              & + nTotalElems * (iDof-1) &
              & + levelDesc%recvbufferFromFiner%elempos(iProc)%val(ielem)
          end do
        end do
      end do
      call commPattern%initBuf_real( &
        &    me = levelDesc%recvbufferFromFiner%buf_real(iProc), &
        &    pos = bufEID, &
        &    nVals = posIndex)
    end do
    ! for receiving (from Coarser)
    do iProc = 1, levelDesc%recvbufferFromCoarser%nProcs
      if(allocated(levelDesc%recvbufferFromCoarser%elemPos(iProc)%val)) then
        elemPosSize = levelDesc%recvbufferFromCoarser%elemPos(iProc)%nvals
      else
        elemPosSize = 0
      end if
      posIndex = 0
      do iElem = 1, elemPosSize
        do iVar = 1, nVars
          do iDof = 1, nDofs
            posIndex = posIndex + 1
            bufEID(posIndex) = nTotalElems * nDofs * (iVar-1) &
              & + nTotalElems * (iDof-1) &
              & + levelDesc%recvbufferFromCoarser%elempos(iProc)%val(ielem)
          end do
        end do
      end do
      call commPattern%initBuf_real( &
        &    me = levelDesc%recvbufferFromCoarser%buf_real(iProc), &
        &    pos = bufEID, &
        &    nVals = posIndex)
    end do

    deallocate(bufEID)

  end subroutine atl_init_cellStateBuffer

end module atl_parallel_module
