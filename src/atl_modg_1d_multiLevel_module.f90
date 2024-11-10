! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
!! Collections of operations and datatypes related to multilevel simulations for
!! the MODG scheme.
module atl_modg_1d_multilevel_module
  use env_module,                only: rk

  use tem_topology_module,       only: tem_childNumber
  use tem_element_module,        only: eT_fluid,            &
    &                                  eT_ghostFromCoarser, &
    &                                  eT_ghostFromFiner
  use tem_faceData_module,       only: tem_invFace_map

  use atl_scheme_module,         only: atl_scheme_type
  use atl_modg_1d_scheme_module, only: atl_modg_1d_scheme_type
  use atl_cube_elem_module,      only: atl_cube_elem_type
  use atl_facedata_module,       only: atl_facedata_type
  use atl_kerneldata_module,     only: atl_statedata_type

  use ply_modg_basis_module,     only: ply_modg_basis_type

  implicit none
  private

  public :: atl_modg_1d_coarseToFineFace, atl_modg_1d_fineToCoarseFace
  public :: atl_modg_coarseToFineElem_1d, atl_modg_fineToCoarseElem_1d

contains



  !> \brief Interpolate modal face representation from coarse to next finer faces (level
  !! difference between coarser and finer faces has to be 1).
  !!
  !! Interpolates functions defined on faces of the current level to faces
  !! of the next finer level. \n
  !! \n
  !!      faces on fine                            face on current           \n
  !!          level                                     level                \n
  !! ------------------------                 ------------------------       \n
  !! |          |           |                 |                      |       \n
  !! |    3     |     4     |                 |                      |       \n
  !! |          |           |                 |                      |       \n
  !! ------------------------   <<----------  |          5           |       \n
  !! |          |           |                 |                      |       \n
  !! |    1     |     2     |                 |                      |       \n
  !! |          |           |                 |                      |       \n
  !! ------------------------                 ------------------------       \n
  !! \n
  !! This is accomplished with lower complexity (with respect to the polynomial
  !! degree) by a dimension by dimension approach: \n
  !! \n
  !!      faces on fine                            face on current           \n
  !!          level                                     level                \n
  !! ------------------------                 ------------------------       \n
  !! |          |           |                 |                      |       \n
  !! |    3     |     4     |                 |                      |       \n
  !! |          |           |                 |                      |       \n
  !! ------------------------   <<----------  |          5           |       \n
  !! |          |           |                 |                      |       \n
  !! |    1     |     2     |                 |                      |       \n
  !! |          |           |                 |                      |       \n
  !! ------------------------                 ------------------------       \n
  !!                 \                               /                       \n
  !!                  \                             /                        \n
  !!                   \                           /                         \n
  !!                    \                         /                          \n
  !!                     ------------------------                            \n
  !!                     |                      |                            \n
  !!                     |          b           |                            \n
  !!                     |                      |                            \n
  !!                     ------------------------                            \n
  !!                     |                      |                            \n
  !!                     |          a           |                            \n
  !!                     |                      |                            \n
  !!                     ------------------------                            \n
  !!
  subroutine atl_modg_1d_coarseToFineFace( minLevel, maxLevel, currentLevel, &
                                  & mesh, facedata, nScalars)
    ! --------------------------------------------------------------------------
    !> The minumum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The current level (i.e. the coarse level).
    integer, intent(in) :: currentLevel
    !> The mesh representation.
    type(atl_cube_elem_type), intent(in) :: mesh(minLevel:maxLevel)
    !> The face representations (finer faces are interpolated from coarser ones).
    type(atl_facedata_type), intent(inout) :: facedata(minLevel:maxLevel)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    ! --------------------------------------------------------------------------
    integer :: iDir, iAlign, iFace, invAlign
    integer :: nDofsCoarse, nDofsFine, nDofsInter
    integer :: elemPos, childPos(4)
    real(kind=rk), allocatable :: faceDat(:,:)
    real(kind=rk), allocatable :: aFace(:,:), bFace(:,:)
    real(kind=rk), allocatable :: firstFace(:,:), secondFace(:,:)
    real(kind=rk), allocatable :: thirdFace(:,:), fourthFace(:,:)
    ! --------------------------------------------------------------------------

    ! The number of degrees of freedom on a face for the coarse, fine and intermediate level
    nDofsCoarse = 1
    nDofsFine = 1
    nDofsInter = 1

    ! Create the intermediate and resulting arrays
    allocate( faceDat(nDofsCoarse,nScalars), &
            & aFace(nDofsInter,nScalars), bFace(nDofsInter,nScalars), &
            & firstFace(nDofsFine,nScalars), secondFace(nDofsFine,nScalars), &
            & thirdFace(nDofsFine,nScalars), fourthFace(nDofsFine,nScalars) )


    ! Iterate over all the from finer faces of the current level and project
    ! the polynomials downwards to the next finer level.
    do iDir = 1,1
      do iAlign = 1,2
        do iFace = 1, size(mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%elemPos)

          ! The face we have to interpolate. If the face is refined from its left element
          ! (i.e. the right face of the element element is refined, -> iAlign == 2), we have
          ! to work on the left face of the right element.
          invAlign = tem_invFace_map(iAlign)

          ! Get the face representation (i.e. face 5)
          elemPos = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%elemPosOp(iFace)
          childPos(1) = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%childPosOp(1,iFace)
          childPos(2) = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%childPosOp(2,iFace)
          childPos(3) = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%childPosOp(3,iFace)
          childPos(4) = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%childPosOp(4,iFace)
          faceDat = facedata(currentLevel)%faceRep(iDir)%dat(elemPos,:,:,invAlign)
          facedata(currentLevel+1)%faceRep(iDir)%dat(childPos(1),:,:,invAlign) = faceDat
          facedata(currentLevel+1)%faceRep(iDir)%dat(childPos(2),:,:,invAlign) = faceDat
          facedata(currentLevel+1)%faceRep(iDir)%dat(childPos(3),:,:,invAlign) = faceDat
          facedata(currentLevel+1)%faceRep(iDir)%dat(childPos(4),:,:,invAlign) = faceDat

        end do
      end do
    end do


  end subroutine atl_modg_1d_coarseToFineFace




  !> summary: Interpolate modal face representation from next finer faces to coarse level (level
  !! difference between coarser and finer faces has to be 1).
  !!
  !! Interpolates functions defined on finer faces to faces of the current level.\n
  !! \n
  !!      faces on fine                            face on current
  !!          level                                     level
  !! ------------------------                 ------------------------
  !! |          |           |                 |                      |
  !! |    3     |     4     |                 |                      |
  !! |          |           |                 |                      |
  !! ------------------------   ---------->>  |          5           |
  !! |          |           |                 |                      |
  !! |    1     |     2     |                 |                      |
  !! |          |           |                 |                      |
  !! ------------------------                 ------------------------
  !! \n
  !! This is accomplished with lower complexity (with respect to the polynomial
  !! degree) by a dimension by dimension approach: \n
  !! \n
  !!      faces on fine                            face on current           \n
  !!          level                                     level                \n
  !! ------------------------                 ------------------------       \n
  !! |          |           |                 |                      |       \n
  !! |    3     |     4     |                 |                      |       \n
  !! |          |           |                 |                      |       \n
  !! ------------------------   ---------->>  |          5           |       \n
  !! |          |           |                 |                      |       \n
  !! |    1     |     2     |                 |                      |       \n
  !! |          |           |                 |                      |       \n
  !! ------------------------                 ------------------------       \n
  !!                 \                               /                       \n
  !!                  \                             /                        \n
  !!                   \                           /                         \n
  !!                    \                         /                          \n
  !!                     ------------------------                            \n
  !!                     |                      |                            \n
  !!                     |          b           |                            \n
  !!                     |                      |                            \n
  !!                     ------------------------                            \n
  !!                     |                      |                            \n
  !!                     |          a           |                            \n
  !!                     |                      |                            \n
  !!                     ------------------------                            \n
  !!
  subroutine atl_modg_1d_fineToCoarseFace( minLevel, maxLevel, currentLevel, &
                                  & mesh, facedata, nScalars)
    ! --------------------------------------------------------------------------
    !> The minumum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The current level (i.e. the coarse level).
    integer, intent(in) :: currentLevel
    !> The mesh representation.
    type(atl_cube_elem_type), intent(in) :: mesh(minLevel:maxLevel)
    !> The face representations (finer faces are interpolated from coarser ones).
    type(atl_facedata_type), intent(inout) :: facedata(minLevel:maxLevel)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    ! --------------------------------------------------------------------------
    integer :: iDir, iAlign, iFace, invAlign
    integer :: nDofsCoarse, nDofsFine, nDofsInter
    integer :: elemPos, childPos(4)
    real(kind=rk), allocatable :: faceDat(:,:), aFace(:,:), bFace(:,:)
    real(kind=rk), allocatable :: firstDat(:,:), secondDat(:,:)
    real(kind=rk), allocatable :: thirdDat(:,:), fourthDat(:,:)
    ! --------------------------------------------------------------------------

    ! The number of degrees of freedom on a face
    nDofsCoarse = 1
    nDofsFine = 1
    nDofsInter = 1

    ! Create the intermediate and resulting arrays
    allocate( faceDat(nDofsCoarse,nScalars), &
            & aFace(nDofsInter,nScalars), bFace(nDofsInter,nScalars), &
            & firstDat(nDofsFine,nScalars), secondDat(nDofsFine,nScalars), &
            & thirdDat(nDofsFine,nScalars), fourthDat(nDofsFine,nScalars) )


    ! Iterate over all the faces and project from the faces of the finer level
    ! to the current level.
    do iDir = 1,1
      do iAlign = 1,2
        do iFace = 1, size(mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%elemPos)

          ! The face we have to interpolate. If the face is refined from its left element
          ! (i.e. the right face of the element element is refined, -> iAlign == 2), we have
          ! to work on the left face of the right element.
          invAlign = tem_invFace_map(iAlign)

          ! Get the data of the faces 1,2,3 and 4
          ! ... face 1
          elemPos = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%elemPosOp(iFace)
          childPos(1) = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%childPosOp(1,iFace)
          childPos(2) = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%childPosOp(2,iFace)
          childPos(3) = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%childPosOp(3,iFace)
          childPos(4) = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%childPosOp(4,iFace)
          firstDat = facedata(currentLevel+1)%faceFlux(iDir)%dat(childPos(1),:,:,invAlign)
          secondDat = facedata(currentLevel+1)%faceFlux(iDir)%dat(childPos(2),:,:,invAlign)
          thirdDat = facedata(currentLevel+1)%faceFlux(iDir)%dat(childPos(3),:,:,invAlign)
          fourthDat = facedata(currentLevel+1)%faceFlux(iDir)%dat(childPos(4),:,:,invAlign)
          facedata(currentLevel)%faceFlux(iDir)%dat(elemPos,:,:,invAlign) = &
              & ( firstDat + secondDat + thirdDat + fourthDat ) / 4.0_rk

        end do
      end do
    end do


  end subroutine atl_modg_1d_fineToCoarseFace

  !> Project data from 8 smaller elements to its parent element in terms
  !! of L2 projections.
  subroutine atl_modg_fineToCoarseElem_1d( minLevel, maxLevel, currentLevel, &
    &                                      iDir, mesh, state_stab, scheme,   &
    &                                      nScalars                          )
    ! --------------------------------------------------------------------------
    !> The minumum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The current level (i.e. the coarse level).
    integer, intent(in) :: currentLevel
    !> The direction to interpolate.
    integer, intent(in) :: iDir
    !> The mesh representation.
    type(atl_cube_elem_type), intent(in) :: mesh(minLevel:maxLevel)
    !> The face representations (finer faces are interpolated from coarser ones).
    type(atl_statedata_type), intent(inout) :: state_stab(minLevel:maxLevel,1:3)
    !> The schemes on the different levels.
    type(atl_scheme_type), intent(in) :: scheme(minLevel:maxLevel)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    ! --------------------------------------------------------------------------
    integer :: iElem, iChild
    integer :: iRefineX, iRefineY
    integer :: nDofsCoarse, nDofsFine
    integer :: elemPos, childPos
    integer :: iDof, iVar, i
    real(kind=rk), allocatable :: faceDat(:,:)
    real(kind=rk), allocatable :: childFace(:,:,:)
    real(kind=rk), allocatable :: firstCoarse(:,:,:), secondCoarse(:,:,:)
    ! --------------------------------------------------------------------------

    ! The number of degrees of freedom on a face
    nDofsCoarse = scheme(currentLevel)%nDofs
    nDofsFine = scheme(currentLevel+1)%nDofs

    ! Create the intermediate and resulting arrays
    allocate( childFace(nDofsFine,nScalars,8) )
    allocate( firstCoarse(nDofsFine,nScalars,4) )
    allocate( secondCoarse(nDofsFine,nScalars,2) )
    allocate( faceDat(nDofsCoarse,nScalars) )


    ! Iterate over all the from finer elements and project from the elements of
    ! the finer level to the current level.
    do iElem = 1, mesh(currentLevel)%faces_stab              &
      &                             %dimByDimDesc(iDir)      &
      &                             %elem                    &
      &                             %nElems(eT_ghostFromFiner)

      do iChild = 1,8
        childPos = mesh(currentLevel)%faces_stab         &
          &                          %dimByDimDesc(iDir) &
          &                          %depFromFiner(iElem)%elem%val(iChild)
        do i=1,nDofsFine*nScalars
          iDof = mod(i-1,nDofsFine) + 1
          iVar = (i-1)/nDofsFine + 1
          childFace(iDof,iVar,iChild) = state_stab(currentLevel+1,iDir) &
            &                           %state(childPos,iDof,iVar)
        end do
      end do

      ! First coarsening step in z direction (just average in this direction)
      do iRefineY = 1,2
        do iRefineX = 1,2
          do i=1,nDofsFine*nScalars
            iDof = mod(i-1,nDofsFine) + 1
            iVar = (i-1)/nDofsFine + 1
            firstCoarse(iDof,iVar,(iRefineY-1)*2+iRefineX)                    &
              &  = 0.5_rk * ( childFace(iDof,iVar,(iRefineY-1)*2+iRefineX)    &
              &              + childFace(iDof,iVar,4+(iRefineY-1)*2+iRefineX) )
          end do
        end do
      end do

      ! Second coarsening step in y direction
      do iRefineX = 1,2
        do i=1,nDofsFine*nScalars
          iDof = mod(i-1,nDofsFine) + 1
          iVar = (i-1)/nDofsFine + 1
          secondCoarse(iDof,iVar,iRefineX)                      &
            &  = 0.5_rk * ( firstCoarse(iDof,iVar,iRefineX)     &
            &               + firstCoarse(iDof,iVar,2+iRefineX) )
        end do
      end do

      ! Third coarsening step in x direction
      do i=1,nDofsCoarse*nScalars
        iDof = mod(i-1,nDofsCoarse) + 1
        iVar = (i-1)/nDofsCoarse + 1
        faceDat(iDof,iVar) = 0.0_rk
      end do

      do iRefineX = 1,2
        call modg_semiCoarseElem_1d(                         &
          & modalRepFace  = secondCoarse(:,:,iRefineX),      &
          & modg_basis    = scheme(currentLevel)%modg_basis, &
          & schemeCoarse  = scheme(currentLevel)%modg_1d,    &
          & schemeFine    = scheme(currentLevel+1)%modg_1d,  &
          & fineElemShift = iRefineX,                        &
          & modalCoarsed  = faceDat                          )
      end do

      ! Assign the data (same for both sides of the face)
      elemPos = iElem + mesh(currentLevel)%faces_stab                &
        &                                 %dimByDimDesc(iDir)        &
        &                                 %elem                      &
        &                                 %nElems(eT_fluid)          &
        &             + mesh(currentLevel)%faces_stab                &
        &                                 %dimByDimDesc(iDir)        &
        &                                 %elem                      &
        &                                 %nElems(eT_ghostFromCoarser)

      do i=1,nDofsCoarse*nScalars
        iDof = mod(i-1,nDofsCoarse) + 1
        iVar = (i-1)/nDofsCoarse + 1
        state_stab(currentLevel,iDir)%state(elemPos,iDof,iVar) &
          &  = faceDat(iDof,iVar)
      end do

    end do


  end subroutine atl_modg_fineToCoarseElem_1d


  !> Subroutine to semi-coarsen an element with modal polynomial representation
  !! to its semi-parent.
  subroutine modg_semiCoarseElem_1d( modalRepFace, modg_basis, schemeCoarse, &
    &                                schemeFine, fineElemShift, modalCoarsed )
    ! --------------------------------------------------------------------------
    !> Modal representation of a function on one of refined element.
    !! Which fine element is determined by fineFaceShift
    !! Dimensions are: (modg%maxPolyDegree+1)^3 for the first dimension
    !! and nScalars for the second dimension.
    real(kind=rk), intent(in) :: modalRepFace(:,:)
    !> The polynomial basis for the  current level of the modg scheme
    type(ply_modg_basis_type), intent(in) :: modg_basis
    !> The parameters of your MODG scheme on the coarse level.
    type(atl_modg_1d_scheme_type), intent(in) :: schemeCoarse
    !> The parameters of your MODG scheme on the fint level.
    type(atl_modg_1d_scheme_type), intent(in) :: schemeFine
    !> The semi-refined element you want to obtain.
    integer, intent(in) :: fineElemShift
    !> The modal representation of modalRepFace on the coarser element, restricted
    !! to the given fine element.
    real(kind=rk), intent(inout) :: modalCoarsed(:,:)
    ! --------------------------------------------------------------------------
    integer :: iFunc, iDegX
    integer :: funcPos, dof
    real(kind=rk) :: coarseSqNorm, jacobiDetFineToCoarse
    integer :: mpd1
    ! --------------------------------------------------------------------------

    ! Ratio of determinants of the Jacobians of the mapping from fine to reference and
    ! coarse to reference (one-dimensional)
    jacobiDetFineToCoarse = 0.5_rk

    mpd1 = schemeCoarse%maxPolyDegree+1

    do dof = 1, mpd1
      iDegX = dof
      coarseSqNorm = 2.0_rk / ( 2.0_rk * iDegX - 1.0_rk)
      ! ... loop over the ansatz function of the semi-refined element.
      do iFunc = 1, schemeFine%maxPolyDegree+1
  funcpos = ifunc                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(schemefine%maxpolydegree+1))*(schemefine%maxpolydegree+1)
        modalCoarsed(dof,:) = modalCoarsed(dof,:)                  &
          & + modg_basis%refineBaseCoeff                           &
          &             %anz_anzShift(iFunc, iDegX, fineElemShift) &
          & * modalRepFace(funcPos,:)                              &
          & * jacobiDetFineToCoarse / ( coarseSqNorm )
      end do
    end do


  end subroutine modg_semiCoarseElem_1d

  !> Project coarse parent element to its 8 finer child elements
  !! by a simple L2 projection.
  subroutine atl_modg_coarseToFineElem_1d( minLevel, maxLevel, currentLevel, &
    &                                      iDir, mesh, state_stab, scheme,   &
    &                                      nScalars                          )
    ! --------------------------------------------------------------------------
    !> The minumum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The current level (i.e. the coarse level).
    integer, intent(in) :: currentLevel
    !> The direction to project
    integer, intent(in) :: iDir
    !> The mesh representation.
    type(atl_cube_elem_type), intent(in) :: mesh(minLevel:maxLevel)
    !> The face representations (finer faces are interpolated from coarser ones).
    type(atl_statedata_type), intent(inout) :: state_stab(minLevel:maxLevel,1:3)
    !> The schemes on the different levels.
    type(atl_scheme_type), intent(in) :: scheme(minLevel:maxLevel)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    ! --------------------------------------------------------------------------
    integer :: iElem, childPos, childNum, parentPos
    integer :: iRefineX, iRefineY, iRefineZ
    integer :: nDofsCoarse
    integer :: nDofsInter_firstRefine
    integer :: i, iVar, iDof
    real(kind=rk), allocatable :: faceDat(:,:)
    real(kind=rk), allocatable :: firstRefine(:,:)
    ! --------------------------------------------------------------------------

    ! The number of degrees of freedom on an element for the coarse, fine and intermediate level
    nDofsCoarse = scheme(currentLevel-1)%nDofs
    nDofsInter_firstRefine = (scheme(currentLevel)%modg_1d%maxPolyDegree+1)

    ! Create the intermediate and resulting arrays
    allocate( faceDat(nDofsCoarse,nScalars) )
    allocate( firstRefine(nDofsInter_firstRefine,nScalars) )


    do iElem = 1, mesh(currentLevel)%faces_stab                &
      &                             %dimByDimDesc(iDir)        &
      &                             %elem                      &
      &                             %nElems(eT_ghostFromCoarser)

      ! Get position and number of the fine (child) element
      childPos = iElem + mesh(currentLevel)%faces_stab         &
        &                                  %dimByDimDesc(iDir) &
        &                                  %elem               &
        &                                  %nElems(eT_fluid)
      childNum = tem_childNumber(mesh(currentLevel)%faces_stab         &
        &                                          %dimByDimDesc(iDir) &
        &                                          %total(childPos))

      ! Get position of the parent element
      parentPos = mesh(currentLevel)%faces_stab            &
        &                           %dimByDimDesc(iDir)    &
        &                           %depFromCoarser(iElem) &
        &                           %elem%val(1)

      ! Get state of the coarser element
      do i=1,nDofsCoarse*nScalars
        iDof = mod(i-1,nDofsCoarse) + 1
        iVar = (i-1)/nDofsCoarse + 1
        faceDat(iDof,iVar) = state_stab(currentLevel-1,iDir)%state(parentPos,iDof,iVar)
      end do

      ! Get the fine element shift from the child number
      iRefineZ = (childNum-1)/4 + 1
      iRefineY = (childNum-1-(iRefineZ-1)*4)/2+1
      iRefineX = mod(childNum-1,2)+1

      ! Do the first refinement step in x direction
      do i=1,nDofsInter_firstRefine*nScalars
        iDof = mod(i-1,nDofsInter_firstRefine) + 1
        iVar = (i-1)/nDofsInter_firstRefine + 1
        firstRefine(iDof,iVar) = 0.0_rk
      end do

      call modg_semiRefineElem_1d(                         &
        & modalRepFace  = faceDat,                         &
        & modg_basis    = scheme(currentLevel)%modg_basis, &
        & schemeCoarse  = scheme(currentLevel-1)%modg_1d,  &
        & schemeFine    = scheme(currentLevel)%modg_1d,    &
        & fineElemShift = iRefineX,                        &
        & modalRefined  = firstRefine                      )

      ! Assign element state for the fine (child) element
      do i=1,nDofsInter_firstRefine*nScalars
        iDof = mod(i-1,nDofsInter_firstRefine) + 1
        iVar = (i-1)/nDofsInter_firstRefine + 1
        state_stab(currentLevel,iDir)%state(childPos,iDof,iVar) = firstRefine(iDof,iVar)
      end do

    end do

  end subroutine atl_modg_coarseToFineElem_1d

  !> Subroutine to semi-refine an element with modal polynomial representation
  !! into its semi-children.
  subroutine modg_semiRefineElem_1d( modalRepFace, modg_basis, schemeCoarse, &
    &                                schemeFine, fineElemShift, modalRefined )
    ! --------------------------------------------------------------------------
    !> Modal representation of a function on the non-refined face.
    !! Dimensions are: (modg%maxPolyDegree+1)^2 for the first dimension
    !! and nScalars for the second dimension.
    real(kind=rk), intent(in) :: modalRepFace(:,:)
    !> The polynomial basis for the  current level of the modg scheme
    type(ply_modg_basis_type), intent(in) :: modg_basis
    !> The parameters of your MODG scheme on the coarse level.
    type(atl_modg_1d_scheme_type), intent(in) :: schemeCoarse
    !> The parameters of your MODG scheme on the fine level.
    type(atl_modg_1d_scheme_type), intent(in) :: schemeFine
    !> The semi-refined element you want to obtain.
    integer, intent(in) :: fineElemShift
    !> The modal representation of modalRepFace restricted to the semi-refined
    !! element.
    real(kind=rk), intent(inout) :: modalRefined(:,:)
    ! --------------------------------------------------------------------------
    integer :: iDegX, iCoarseFunc
    integer :: dof, coarsePos
    real(kind=rk) :: fineSqNorm
    integer :: mpd1
    ! --------------------------------------------------------------------------


    mpd1 = schemeFine%maxPolyDegree+1

    ! loop over all degrees of freedoms on the (semi-)refined element
    do dof = 1, mpd1
      iDegX = dof
      fineSqNorm = 2.0_rk /(2.0_rk * iDegX - 1.0_rk)
      ! ... loop over the ansatz functions of the current (non-refined) element
      do iCoarseFunc = 1, schemeCoarse%maxPolyDegree+1
  coarsepos = icoarsefunc                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(schemecoarse%maxpolydegree+1))*(schemecoarse%maxpolydegree+1)
        ! Weight and add it to the semi refined face representation
        modalRefined(dof,:) = modalRefined(dof,:)                        &
          & + modg_basis%refineBaseCoeff                                 &
          &             %anz_anzShift(iDegX, iCoarseFunc, fineElemShift) &
          & * modalRepFace(coarsePos,:) / ( fineSqNorm )
      end do
    end do

  end subroutine modg_semiRefineElem_1d

end module atl_modg_1d_multilevel_module
