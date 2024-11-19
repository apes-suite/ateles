! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013, 2015-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
!!
!! Due to a compiler problem in Intel 15, we had to replace array syntax by do
!! loops below, and use OMP DO instead of WORKSHARE constructs.
module atl_modg_2d_multilevel_module
  use mpi
  use env_module,                only: rk

  use tem_topology_module,       only: tem_childNumber
  use tem_faceData_module,       only: tem_invFace_map
  use tem_element_module,        only: eT_fluid,            &
    &                                  eT_ghostFromCoarser, &
    &                                  eT_ghostFromFiner
  use tem_timer_module,         only: tem_startTimer, &
    &                                 tem_stopTimer

  use atl_scheme_module,         only: atl_scheme_type
  use atl_modg_2d_scheme_module, only: atl_modg_2d_scheme_type
  use atl_cube_elem_module,      only: atl_cube_elem_type
  use atl_facedata_module,       only: atl_facedata_type
  use atl_kerneldata_module,     only: atl_statedata_type
  use atl_timer_module,          only: atl_elemTimers
  use ply_modg_basis_module,     only: ply_modg_basis_type

  implicit none
  private

  public :: atl_modg_2d_coarseToFineFace, atl_modg_2d_fineToCoarseFace
  public :: atl_modg_coarseToFineElem_2d, atl_modg_fineToCoarseElem_2d


contains


  !> summary: Interpolate modal face representation from coarse to next finer faces (level
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
  subroutine atl_modg_2d_coarseToFineFace( minLevel, maxLevel, currentLevel, &
                                  & mesh, facedata, scheme , nScalars)
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
    !> The schemes on the different levels.
    type(atl_scheme_type), intent(in) :: scheme(minLevel:maxLevel)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    ! --------------------------------------------------------------------------
    integer :: iDir, iAlign, iFace, invAlign
    integer :: iChild
    integer :: i, iDof, iVar
    integer :: nDofsCoarse, nDofsFine, nDofsInter
    integer :: nFluids
    integer :: elemPos, childPos(4)
    integer :: totpos
    integer :: coff, yshift
    integer :: nFaces
    real(kind=rk), allocatable :: faceDat(:,:)
    real(kind=rk), allocatable :: aFace(:,:)
    real(kind=rk), allocatable :: firstFace(:,:), secondFace(:,:)
    ! --------------------------------------------------------------------------

    ! The number of degrees of freedom on a face for the coarse, fine and
    ! intermediate level
    nDofsCoarse = scheme(currentLevel)%modg_2d%maxPolyDegree+1
    nDofsFine = scheme(currentLevel+1)%modg_2d%maxPolyDegree+1
    nDofsInter = (scheme(currentLevel)%modg_2d%maxPolyDegree+1)*1
    nfluids = mesh(currentLevel)%descriptor%elem%nElems(eT_fluid)

    ! Create the intermediate and resulting arrays
    allocate( faceDat(nDofsCoarse,nScalars), &
            & aFace(nDofsInter,nScalars), &
            & firstFace(nDofsFine,nScalars), secondFace(nDofsFine,nScalars) )

    ! Iterate over all the from finer faces of the current level and project
    ! the polynomials downwards to the next finer level.
    do iDir = 1,2
      do iAlign = 1,2
        nFaces = size( mesh(currentLevel)%faces%faces(iDir)           &
          &                                    %fromFinerFace(iAlign) &
          &                                    %elemPos               )
        do iFace = 1,nFaces

          ! The face we have to interpolate. If the face is refined from its
          ! left element (i.e. the right face of the element element is refined,
          !               -> iAlign == 2), we have to work on the left face of
          ! the right element.
          invAlign = tem_invFace_map(iAlign)

          ! Get the face representation (i.e. face 5)
          elemPos = mesh(currentLevel)%faces%faces(iDir)           &
            &                               %fromFinerFace(iAlign) &
            &                               %elemPosOp(iFace)
          if (elempos <= nFluids) then
            totpos = mesh(currentLevel)%descriptor%pntTID(elempos)
            call tem_startTimer( me          = atl_elemTimers, &
              &                  timerHandle = totPos          )
          end if
          do iChild=1,4
            childPos(iChild) = mesh(currentLevel)%faces%faces(iDir)      &
              &                                  %fromFinerFace(iAlign)  &
              &                                  %childPosOp(iChild,iFace)
          end do

          do i=1,nDofsCoarse*nScalars
            iDof = mod(i-1,nDofsCoarse)+1
            iVar = (i-1)/nDofsCoarse + 1
            faceDat(iDof,iVar) = facedata(currentLevel)%faceRep(iDir) &
              &                  %dat(elemPos,iDof,iVar,invAlign)
          end do

          ! 2 Slices of refinement in y
          do yshift=1,2
            do i=1,nDofsInter*nScalars
              iDof = mod(i-1,nDofsInter)+1
              iVar = (i-1)/nDofsInter + 1
              aFace(iDof,iVar) = 0.0_rk
            end do

            do i=1,nDofsFine*nScalars
              iDof = mod(i-1,nDofsFine)+1
              iVar = (i-1)/nDofsFine + 1
              firstFace(iDof,iVar)  = 0.0_rk
              secondFace(iDof,iVar) = 0.0_rk
            end do

            ! Project face 5 to intermediate face:
            call modg_2d_semiRefineFace(                         &
              & modalRepFace  = faceDat,                         &
              & modg_basis    = scheme(currentLevel)%modg_basis, &
              & schemeCoarse  = scheme(currentLevel)%modg_2d,    &
              & schemeFine    = scheme(currentLevel+1)%modg_2d,  &
              & refineDir     = 2,                               &
              & fineFaceShift = yshift,                          &
              & modalRefined  = aFace                            )

            ! Project intermediate face a to fine faces:

            ! ... project a to 1
            call modg_2d_semiRefineFace(                         &
              & modalRepFace  = aFace,                           &
              & modg_basis    = scheme(currentLevel)%modg_basis, &
              & schemeCoarse  = scheme(currentLevel)%modg_2d,    &
              & schemeFine    = scheme(currentLevel+1)%modg_2d,  &
              & refineDir     = 1,                               &
              & fineFaceShift = 1,                               &
              & modalRefined  = firstFace                        )

            ! ... project a to 2
            call modg_2d_semiRefineFace(                         &
              & modalRepFace  = aFace,                           &
              & modg_basis    = scheme(currentLevel)%modg_basis, &
              & schemeCoarse  = scheme(currentLevel)%modg_2d,    &
              & schemeFine    = scheme(currentLevel+1)%modg_2d,  &
              & refineDir     = 1,                               &
              & fineFaceShift = 2,                               &
              & modalRefined  = secondFace                       )

            ! Assign the interpolated data to the child faces.
            coff = (yshift-1)*2
            do i=1,nDofsFine*nScalars
              iDof = mod(i-1,nDofsFine)+1
              iVar = (i-1)/nDofsFine + 1
              facedata(currentLevel+1)%faceRep(iDir)           &
                &                     %dat(childPos(1+coff),   &
                &                          iDof,iVar,invAlign) &
                &  = firstFace(iDof,iVar)
              facedata(currentLevel+1)%faceRep(iDir)           &
                &                     %dat(childPos(2+coff),   &
                &                          iDof,iVar,invAlign) &
                &  = secondFace(iDof,iVar)
            end do
          end do

          if (elemPos <= nFluids) then
            call tem_stopTimer( me          = atl_elemTimers, &
              &                 timerHandle = totPos          )
          end if

        end do
      end do
    end do

  end subroutine atl_modg_2d_coarseToFineFace



  !> summary: Project modal representation on a face to a semi-refined face (i.e. on a face
  !! that is refined in one of the spatial directions.
  !!
  !! Project modal representation on a face to a semi-refined face (i.e. on a face
  !! that is refined in one of the spatial directions. The result is a modal representation
  !! on the semi-refined element. \n
  !! \n
  !! The function is executing one of the following projections: \n
  !!
  !!       face on current                      semi-refined                            \n
  !!            level                               face                                \n
  !!
  !!  ------------------------             ------------------------                     \n
  !!  |                      |             |          |           |                     \n
  !!  |                      |  refineDir  |          |           |                     \n
  !!  |                      |      == 1   | fineFace | fineFace  |                     \n
  !!  |                      | -------->>  |  Shift   |  Shift    |                     \n
  !!  |                      |             |   == 1   |   == 2    |                     \n
  !!  |                      |             |          |           |                     \n
  !!  |                      |             |          |           |                     \n
  !!  ------------------------             ------------------------                     \n
  !!
  !!  or:
  !!
  !!  ------------------------             ------------------------                     \n
  !!  |                      |             |                      |                     \n
  !!  |                      |  refineDir  |     fineFaceShift    |                     \n
  !!  |                      |      == 2   |        == 2          |                     \n
  !!  |                      | -------->>  ------------------------                     \n
  !!  |                      |             |                      |                     \n
  !!  |                      |             |     fineFaceShift    |                     \n
  !!  |                      |             |        == 1          |                     \n
  !!  ------------------------             ------------------------                     \n
  !!
  subroutine modg_2d_semiRefineFace( modalRepFace, modg_basis, schemeCoarse, &
    &                                schemeFine, refineDir, fineFaceShift,   &
    &                                modalRefined                            )
    ! --------------------------------------------------------------------------
    !> Modal representation of a function on the non-refined face.
    !! Dimensions are: (modg%maxPolyDegree+1)^2 for the first dimension
    !! and nScalars for the second dimension.
    real(kind=rk), intent(in) :: modalRepFace(:,:)
    !> Informations about the polynomial basis of the MODG scheme.
    type(ply_modg_basis_type) :: modg_basis
    !> The parameters of your MODG scheme on the coarse level.
    type(atl_modg_2d_scheme_type), intent(in) :: schemeCoarse
    !> The parameters of your MODG scheme on the fine level.
    type(atl_modg_2d_scheme_type), intent(in) :: schemeFine
    !> The direction of the semi-refinement. Either 1 or 2. Have a look at the
    !! function description.
    integer, intent(in) :: refineDir
    !> The semi-refined element you want to obtain.
    integer, intent(in) :: fineFaceShift
    !> The modal representation of modalRepFace restricted to the semi-refined
    !! element.
    real(kind=rk), intent(inout) :: modalRefined(:,:)
    ! --------------------------------------------------------------------------
    integer :: iFunc, iCoarseFunc, comFunc
    integer :: funcPos, coarsePos
    real(kind=rk) :: fineSqNorm
    ! --------------------------------------------------------------------------

    ! Loop over all the ansatz functions of the semi-refined element and calculate
    ! all the modal coefficients for it.

    ! We treat the refinement directions separately to have the if condition
    ! outside of the loop.
    if( refineDir .eq. 1 ) then

      ! ... loop over the ansatz function that are common for the fine and current element
      comFunc = 1
      ! ... loop over the ansatz function of the semi-refined element.
      do iFunc = 1, schemeFine%maxPolyDegree+1
  funcpos = ifunc                                      &
    &      + ( ( comfunc-1)                             &
    &      + (1-1)*(schemefine%maxpolydegree+1))*(schemefine%maxpolydegree+1)
        fineSqNorm = 2.0_rk /(2.0_rk * iFunc - 1.0_rk)
        ! ... loop over the ansatz functions of the current (non-refined) element
        do iCoarseFunc = 1, schemeCoarse%maxPolyDegree+1
          coarsePos = 1                                      &
            & + (iCoarseFunc-1)                              &
            & + (comFunc-1) * (schemeCoarse%maxPolyDegree + 1)
          ! Weight and add it to the semi refined face representation
          modalRefined(funcPos,:) = modalRefined(funcPos,:)                &
            & + modg_basis%refineBaseCoeff                                 &
            &             %anz_anzShift(iFunc, iCoarseFunc, fineFaceShift) &
            & * modalRepFace(coarsePos,:) / ( fineSqNorm )
        end do
      end do

    else

      ! ... loop over the ansatz function that are common for the fine and current element
      do comFunc = 1, schemeCoarse%maxPolyDegree+1
        ! ... loop over the ansatz function of the semi-refined element.
        iFunc = 1
        funcPos = 1 + (comFunc-1) + (iFunc-1)*(schemeCoarse%maxPolyDegree+1)
        fineSqNorm = 2.0_rk /(2.0_rk * iFunc - 1.0_rk)
        ! ... loop over the ansatz functions of the current (non-refined) element
        iCoarseFunc = 1
  coarsepos = comfunc                                      &
    &      + ( ( icoarsefunc-1)                             &
    &      + (1-1)*(schemecoarse%maxpolydegree+1))*(schemecoarse%maxpolydegree+1)
        modalRefined(funcPos,:) = modalRefined(funcPos,:)                &
          & + modg_basis%refineBaseCoeff                                 &
          &             %anz_anzShift(iFunc, iCoarseFunc, fineFaceShift) &
          & * modalRepFace(coarsePos,:) / ( fineSqNorm )
      end do

    end if


  end subroutine modg_2d_semiRefineFace




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
  subroutine atl_modg_2d_fineToCoarseFace( minLevel, maxLevel, currentLevel, &
                                  & mesh, facedata, scheme, nScalars)
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
    !> The schemes on the different levels.
    type(atl_scheme_type), intent(in) :: scheme(minLevel:maxLevel)
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    ! --------------------------------------------------------------------------
    integer :: iDir, iAlign, iFace, invAlign
    integer :: i, iVar, iDof, yshift
    integer :: cp1, cp2
    integer :: nDofsCoarse, nDofsFine, nDofsInter
    integer :: nFluids
    integer :: elemPos, childPos(4)
    integer :: totpos
    real(kind=rk), allocatable :: faceDat(:,:), aFace(:,:)
    real(kind=rk), allocatable :: fineDat(:,:)
    ! --------------------------------------------------------------------------

    nfluids = mesh(currentLevel)%descriptor%elem%nElems(eT_fluid)

    ! The number of degrees of freedom on a face
    nDofsCoarse = scheme(currentLevel)%modg_2d%maxPolyDegree+1
    nDofsFine = scheme(currentLevel+1)%modg_2d%maxPolyDegree+1
    nDofsInter = (scheme(currentLevel)%modg_2d%maxPolyDegree+1)*1

    ! Create the intermediate and resulting arrays
    allocate( faceDat(nDofsCoarse,nScalars), &
            & aFace(nDofsInter,nScalars),    &
            & fineDat(nDofsFine,nScalars)    )

    ! Iterate over all the faces and project from the faces of the finer level
    ! to the current level.
    do iDir = 1,2
      do iAlign = 1,2
        do iFace = 1, size(mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%elemPos)
          ! get elempos
          elemPos = mesh(currentLevel)%faces%faces(iDir)%fromFinerFace(iAlign)%elemPosOp(iFace)

          if (elempos <= nFluids) then
            totpos = mesh(currentLevel)%descriptor%pntTID(elempos)
            ! start element wise timer for LB weights
            call tem_startTimer( me          = atl_elemTimers, &
              &                  timerHandle = totPos          )
          end if


          ! The face we have to interpolate. If the face is refined from its left element
          ! (i.e. the right face of the element element is refined, -> iAlign == 2), we have
          ! to work on the left face of the right element.
          invAlign = tem_invFace_map(iAlign)

          ! Get the data of the faces 1,2,3 and 4
          ! ... face 1
          do i=1,4
            childPos(i) = mesh(currentLevel)%faces%faces(iDir)           &
              &                                   %fromFinerFace(iAlign) &
              &                                   %childPosOp(i,iFace)
          end do

          do i=1,nDofsCoarse*nScalars
            iDof = mod(i-1,nDofsCoarse)+1
            iVar = (i-1)/nDofsCoarse + 1
            faceDat(iDof,iVar) = 0.0_rk
          end do

          do yshift=0,1

            do i=1,nDofsInter*nScalars
              iDof = mod(i-1,nDofsInter)+1
              iVar = (i-1)/nDofsInter + 1
              aFace(iDof,iVar) = 0.0_rk
            end do

            cp1 = childpos(1 + yshift*2)
            cp2 = childpos(2 + yshift*2)
            do i=1,nDofsFine*nScalars
              iDof = mod(i-1,nDofsFine)+1
              iVar = (i-1)/nDofsFine + 1
              fineDat(iDof,iVar) = facedata(currentLevel+1)%faceFlux(iDir) &
                &                   %dat(cp1,iDof,iVar,invAlign)
            end do

            ! ... project 1 to a
            call modg_2d_semiCoarseFace(                         &
              & modalRepFace  = fineDat,                         &
              & modg_basis    = scheme(currentLevel)%modg_basis, &
              & schemeCoarse  = scheme(currentLevel)%modg_2d,    &
              & schemeFine    = scheme(currentLevel+1)%modg_2d,  &
              & coarseDir     = 1,                               &
              & fineFaceShift = 1,                               &
              & modalCoarsed  = aFace                            )

            do i=1,nDofsFine*nScalars
              iDof = mod(i-1,nDofsFine)+1
              iVar = (i-1)/nDofsFine + 1
              fineDat(iDof,iVar) = facedata(currentLevel+1)%faceFlux(iDir) &
                &                   %dat(cp2,iDof,iVar,invAlign)
            end do

            ! ... project 2 to a
            call modg_2d_semiCoarseFace(                         &
              & modalRepFace  = fineDat,                         &
              & modg_basis    = scheme(currentLevel)%modg_basis, &
              & schemeCoarse  = scheme(currentLevel)%modg_2d,    &
              & schemeFine    = scheme(currentLevel+1)%modg_2d,  &
              & coarseDir     = 1,                               &
              & fineFaceShift = 2,                               &
              & modalCoarsed  = aFace                            )

            ! ... project a to coarse
            call modg_2d_semiCoarseFace(                         &
              & modalRepFace  = aFace,                           &
              & modg_basis    = scheme(currentLevel)%modg_basis, &
              & schemeCoarse  = scheme(currentLevel)%modg_2d,    &
              & schemeFine    = scheme(currentLevel+1)%modg_2d,  &
              & coarseDir     = 2,                               &
              & fineFaceShift = 1+yshift,                        &
              & modalCoarsed  = faceDat                          )
          end do

          do i=1,nDofsCoarse*nScalars
            iDof = mod(i-1,nDofsCoarse)+1
            iVar = (i-1)/nDofsCoarse + 1
            facedata(currentLevel)%faceFlux(iDir)%dat(elemPos,iDof,iVar,invAlign) = faceDat(iDof,iVar)
          end do

          if (elempos <= nFluids) then
            ! start element wise timer for LB weights
            call tem_stopTimer( me          = atl_elemTimers, &
              &                 timerHandle = totPos          )
          end if

        end do
      end do
    end do

  end subroutine atl_modg_2d_fineToCoarseFace



  !> summary: Project modal representation of a semi-refined face (i.e. on a face
  !! that is refined in one of the spatial directions) to a coarse element.
  !!
  !! Project modal representation of a semi-refined face (i.e. on a face
  !! that is refined in one of the spatial directions) to a coarser element.
  !! The result is a modal representation on the coarse element that is representing
  !! the polynomial function of the refined element on the coarse element when restricted
  !! to the refined element again. However the modal representation is done in terms
  !! of the ansatz functions of the coarse element. \n
  !! \n
  !! The function is executing one of the following projections: \n
  !!
  !!       face on current                      semi-refined                            \n
  !!            level                               face                                \n
  !!                                                                                    \n
  !!  ------------------------             ------------------------                     \n
  !!  |                      |             |          |           |                     \n
  !!  |                      |  coarseDir  |          |           |                     \n
  !!  |                      |      == 1   | fineFace | fineFace  |                     \n
  !!  |                      | <<--------  |  Shift   |  Shift    |                     \n
  !!  |                      |             |   == 1   |   == 2    |                     \n
  !!  |                      |             |          |           |                     \n
  !!  |                      |             |          |           |                     \n
  !!  ------------------------             ------------------------                     \n
  !!
  !!  or:
  !!
  !!  ------------------------             ------------------------                     \n
  !!  |                      |             |                      |                     \n
  !!  |                      |  coarseDir  |     fineFaceShift    |                     \n
  !!  |                      |      == 2   |        == 2          |                     \n
  !!  |                      | <<--------  ------------------------                     \n
  !!  |                      |             |                      |                     \n
  !!  |                      |             |     fineFaceShift    |                     \n
  !!  |                      |             |        == 1          |                     \n
  !!  ------------------------             ------------------------                     \n
  !!
  subroutine modg_2d_semiCoarseFace( modalRepFace, modg_basis, schemeCoarse, &
    &                                schemeFine, coarseDir, fineFaceShift,   &
    &                                modalCoarsed                            )
    ! --------------------------------------------------------------------------
    !> Modal representation of a function on one of refined face. Which fine face
    !! is determined by fineFaceShift
    !! Dimensions are: (modg%maxPolyDegree+1)^2 for the first dimension
    !! and nScalars for the second dimension.
    real(kind=rk), intent(in) :: modalRepFace(:,:)
    !> Informations about the polynomial basis of the MODG scheme.
    type(ply_modg_basis_type) :: modg_basis
    !> The parameters of your MODG scheme on the coarse level.
    type(atl_modg_2d_scheme_type), intent(in) :: schemeCoarse
    !> The parameters of your MODG scheme on the fint level.
    type(atl_modg_2d_scheme_type), intent(in) :: schemeFine
    !> The direction of the semi-coarsening. Either 1 or 2. Have a look at the
    !! function description.
    integer, intent(in) :: coarseDir
    !> The semi-refined element you want to obtain.
    integer, intent(in) :: fineFaceShift
    !> The modal representation of modalRepFace on the coarser element, restricted
    !! to the given fine element.
    real(kind=rk), intent(inout) :: modalCoarsed(:,:)
    ! --------------------------------------------------------------------------
    integer :: comFunc, iFunc, iCoarseFunc
    integer :: funcPos, coarsePos
    real(kind=rk) :: coarseSqNorm, jacobiDetFineToCoarse
    ! --------------------------------------------------------------------------


    ! Ratio of determinants of the Jacobians of the mapping from fine to reference and
    ! coarse to reference (one-dimensional)
    jacobiDetFineToCoarse = 0.5_rk

    ! Check for the spatial direction of the coarsening.
    if(coarseDir.eq.1)  then

      ! ... loop over the ansatz function that are common for the fine and current element
      comFunc = 1
      ! ... loop over the ansatz functions of the current (non-refined) element
      do iCoarseFunc = 1, schemeCoarse%maxPolyDegree+1
        coarsePos = 1 &
          & + (iCoarseFunc-1) &
          & + (comFunc-1) * (schemeCoarse%maxPolyDegree + 1)
        coarseSqNorm = 2.0_rk / ( 2.0_rk * iCoarseFunc - 1.0_rk)
        ! ... loop over the ansatz function of the semi-refined element.
        do iFunc = 1, schemeFine%maxPolyDegree+1
  funcpos = ifunc                                      &
    &      + ( ( comfunc-1)                             &
    &      + (1-1)*(schemefine%maxpolydegree+1))*(schemefine%maxpolydegree+1)
          modalCoarsed(coarsePos,:) = modalCoarsed(coarsePos,:)            &
            & + modg_basis%refineBaseCoeff                                 &
            &             %anz_anzShift(iFunc, iCoarseFunc, fineFaceShift) &
            & * modalRepFace(funcPos,:)                                    &
            & * jacobiDetFineToCoarse / ( coarseSqNorm )
        end do
      end do

    else

      ! ... loop over the ansatz function that are common for the fine and current element
      do comFunc = 1, schemeCoarse%maxPolyDegree+1
        ! ... loop over the ansatz functions of the current (non-refined) element
        iCoarseFunc = 1
  coarsepos = comfunc                                      &
    &      + ( ( icoarsefunc-1)                             &
    &      + (1-1)*(schemecoarse%maxpolydegree+1))*(schemecoarse%maxpolydegree+1)
        coarseSqNorm = 2.0_rk / ( 2.0_rk * iCoarseFunc - 1.0_rk)
        ! ... loop over the ansatz function of the semi-refined element.
        iFunc = 1
        funcPos = 1 + (comFunc-1) + (iFunc-1)*(schemeCoarse%maxPolyDegree+1)
        modalCoarsed(coarsePos,:) = modalCoarsed(coarsePos,:)                  &
          & + modg_basis%refineBaseCoeff                                       &
          &             %anz_anzShift(iFunc,iCoarseFunc,fineFaceShift)         &
          & * modalRepFace(funcPos,:) * jacobiDetFineToCoarse / ( coarseSqNorm )
      end do

    end if

  end subroutine modg_2d_semiCoarseFace

  !> Project data from 8 smaller elements to its parent element in terms
  !! of L2 projections.
  subroutine atl_modg_fineToCoarseElem_2d( minLevel, maxLevel, currentLevel, &
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
    integer :: nFluids
    integer :: nDofsCoarse, nDofsFine, nDofsInter_firstCoarse, &
      &        nDofsInter_secondCoarse
    integer :: elemPos, childPos
    integer :: totpos
    real(kind=rk) :: starttime, stoptime
    real(kind=rk), allocatable :: faceDat(:,:)
    real(kind=rk), allocatable :: childFace(:,:,:)
    real(kind=rk), allocatable :: firstCoarse(:,:,:), secondCoarse(:,:,:)
    ! --------------------------------------------------------------------------

    nfluids = mesh(currentLevel)%descriptor%elem%nElems(eT_fluid)

    ! The number of degrees of freedom on a face
    nDofsCoarse = scheme(currentLevel)%nDofs
    nDofsFine = scheme(currentLevel+1)%nDofs
    nDofsInter_firstCoarse =  (scheme(currentLevel+1)%modg_2d%maxPolyDegree+1) &
                           & *(scheme(currentLevel+1)%modg_2d%maxPolyDegree+1)
    nDofsInter_secondCoarse =  (scheme(currentLevel+1)%modg_2d%maxPolyDegree+1) &
                            & *(scheme(currentLevel)%modg_2d%maxPolyDegree+1)

    ! Create the intermediate and resulting arrays
    allocate( childFace(nDofsFine,nScalars,8) )
    allocate( firstCoarse(nDofsInter_firstCoarse,nScalars,4) )
    allocate( secondCoarse(nDofsInter_firstCoarse,nScalars,2) )
    allocate( faceDat(nDofsCoarse,nScalars) )

    ! Iterate over all the from finer elements and project from the elements of the finer level
    ! to the current level.
    do iElem = 1, mesh(currentLevel)%faces_stab              &
      &                             %dimByDimDesc(iDir)      &
      &                             %elem                    &
      &                             %nElems(eT_ghostFromFiner)
      elemPos = iElem + mesh(currentLevel)%faces_stab                &
        &                                 %dimByDimDesc(iDir)        &
        &                                 %elem                      &
        &                                 %nElems(eT_fluid)          &
        &             + mesh(currentLevel)%faces_stab                &
        &                                 %dimByDimDesc(iDir)        &
        &                                 %elem                      &
        &                                 %nElems(eT_ghostFromCoarser)

      starttime = MPI_Wtime()

      do iChild = 1,8
        childPos = mesh(currentLevel)%faces_stab%dimByDimDesc(iDir)%depFromFiner(iElem)%elem%val(iChild)
        childFace(:,:,iChild) = state_stab(currentLevel+1,iDir)%state(childPos,:,:)
      end do

      ! First coarsening step in z direction (just average in this direction)
      firstCoarse(:,:,:) = 0.0_rk
      do iRefineY = 1,2
        do iRefineX = 1,2
          firstCoarse(:,:,(iRefineY-1)*2+iRefineX) = &
               & ( childFace(:,:,(iRefineY-1)*2+iRefineX) &
               & + childFace(:,:,4+(iRefineY-1)*2+iRefineX) )/2.0_rk
        end do
      end do

      ! Second coarsening step in y direction
      secondCoarse(:,:,:) = 0.0_rk
      do iRefineY = 1,2
        do iRefineX = 1,2
          call modg_semiCoarseElem_2d(                                  &
            & modalRepFace  = firstCoarse(:,:,(iRefineY-1)*2+iRefineX), &
            & modg_basis    = scheme(currentLevel)%modg_basis,          &
            & schemeCoarse  = scheme(currentLevel)%modg_2d,             &
            & schemeFine    = scheme(currentLevel+1)%modg_2d,           &
            & coarseDir     = 2,                                        &
            & fineElemShift = iRefineY,                                 &
            & modalCoarsed  = secondCoarse(:,:,iRefineX)                )
        end do
      end do

      ! Third coarsening step in x direction
      do iChild=lbound(faceDat,2),ubound(faceDat,2)
        faceDat(:,iChild) = 0.0_rk
      end do
      do iRefineX = 1,2
        call modg_semiCoarseElem_2d(                         &
          & modalRepFace  = secondCoarse(:,:,iRefineX),      &
          & modg_basis    = scheme(currentLevel)%modg_basis, &
          & schemeCoarse  = scheme(currentLevel)%modg_2d,    &
          & schemeFine    = scheme(currentLevel+1)%modg_2d,  &
          & coarseDir     = 1,                               &
          & fineElemShift = iRefineX,                        &
          & modalCoarsed  = faceDat(:,:)                     )
      end do

      ! Assign the data (same for both sides of the face)
      state_stab(currentLevel,iDir)%state(elemPos,:,:) = faceDat(:,:)
      stoptime = MPI_Wtime()

      do iChild = 1,4
        childPos = mesh(currentLevel)%faces_stab%dimByDimDesc(iDir)  &
          &                                     %depFromFiner(iElem) &
          &                                     %elem%val(iChild)

        if ( childpos <= nFluids ) then
          totpos = mesh(currentLevel)%descriptor%pntTID(childpos)
          atl_elemTimers%duration%val(totPos)        &
            &  = atl_elemTimers%duration%val(totPos) &
            &    + 0.25_rk*(stoptime-starttime)
        end if
      end do

    end do

  end subroutine atl_modg_fineToCoarseElem_2d


  !> Subroutine to semi-coarsen an element with modal polynomial representation
  !! to its semi-parent.
  subroutine modg_semiCoarseElem_2d( modalRepFace, modg_basis, schemeCoarse, &
    &                                schemeFine, coarseDir, fineElemShift,   &
    &                                modalCoarsed                            )
    ! --------------------------------------------------------------------------
    !> Modal representation of a function on one of refined element.
    !! Which fine element is determined by fineFaceShift
    !! Dimensions are: (modg%maxPolyDegree+1)^3 for the first dimension
    !! and nScalars for the second dimension.
    real(kind=rk), intent(in) :: modalRepFace(:,:)
    !> Informations about the polynomial basis of the MODG scheme.
    type(ply_modg_basis_type) :: modg_basis
    !> The parameters of your MODG scheme on the coarse level.
    type(atl_modg_2d_scheme_type), intent(in) :: schemeCoarse
    !> The parameters of your MODG scheme on the fint level.
    type(atl_modg_2d_scheme_type), intent(in) :: schemeFine
    !> The direction of the semi-coarsening. Either 1 or 2 or 3. Have a look at the
    !! function description.
    integer, intent(in) :: coarseDir
    !> The semi-refined element you want to obtain.
    integer, intent(in) :: fineElemShift
    !> The modal representation of modalRepFace on the coarser element, restricted
    !! to the given fine element.
    real(kind=rk), intent(inout) :: modalCoarsed(:,:)
    ! --------------------------------------------------------------------------
    integer :: iFunc, iDegX, iDegY
    integer :: funcPos, dof
    real(kind=rk) :: coarseSqNorm, jacobiDetFineToCoarse
    integer :: mpd1, mpd1_square
    ! --------------------------------------------------------------------------

    ! Ratio of determinants of the Jacobians of the mapping from fine to reference and
    ! coarse to reference (one-dimensional)
    jacobiDetFineToCoarse = 0.5_rk

    ! Check for the spatial direction of the coarsening.
    if(coarseDir.eq.1)  then

      mpd1 = schemeCoarse%maxPolyDegree+1
      mpd1_square = (schemeFine%maxPolyDegree+1)*mpd1

      do dof = 1, mpd1_square
        iDegX = mod(dof-1,mpd1)+1
        iDegY = (dof-1)/mpd1+1
        coarseSqNorm = 2.0_rk / ( 2.0_rk * iDegX - 1.0_rk)
        ! ... loop over the ansatz function of the semi-refined element.
        do iFunc = 1, schemeFine%maxPolyDegree+1
  funcpos = ifunc                                      &
    &      + ( ( idegy-1)                             &
    &      + (1-1)*(schemefine%maxpolydegree+1))*(schemefine%maxpolydegree+1)
          modalCoarsed(dof,:) = modalCoarsed(dof,:)                  &
            & + modg_basis%refineBaseCoeff                           &
            &             %anz_anzShift(iFunc, iDegX, fineElemShift) &
            & * modalRepFace(funcPos,:)                              &
            & * jacobiDetFineToCoarse / ( coarseSqNorm )
        end do
      end do

    else

      mpd1 = schemeCoarse%maxPolyDegree+1
      mpd1_square = (schemeCoarse%maxPolyDegree+1)*mpd1

      do dof = 1, mpd1_square
        iDegX = mod(dof-1,mpd1)+1
        iDegY = (dof-1)/mpd1+1
        coarseSqNorm = 2.0_rk / ( 2.0_rk * iDegY - 1.0_rk)
        ! ... loop over the ansatz function of the semi-refined element.
        do iFunc = 1, schemeFine%maxPolyDegree+1
  funcpos = idegx                                      &
    &      + ( ( ifunc-1)                             &
    &      + (1-1)*(schemefine%maxpolydegree+1))*(schemefine%maxpolydegree+1)
          modalCoarsed(dof,:) = modalCoarsed(dof,:)                  &
            & + modg_basis%refineBaseCoeff                           &
            &             %anz_anzShift(iFunc, iDegY, fineElemShift) &
            & * modalRepFace(funcPos,:)                              &
            & * jacobiDetFineToCoarse / ( coarseSqNorm )
        end do
      end do

    end if

  end subroutine modg_semiCoarseElem_2d


  !> Project coarse parent element to its 8 finer child elements
  !! by a simple L2 projection.
  subroutine atl_modg_coarseToFineElem_2d( minLevel, maxLevel, currentLevel, &
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
    integer :: iElem, childPos, childNum, parentPos, ii
    integer :: totpos
    integer :: iRefineX, iRefineY, iRefineZ
    integer :: nDofsCoarse, nDofsFine
    integer :: nDofsInter_firstRefine, nDofsInter_secondRefine
    integer :: nFluids
    real(kind=rk), allocatable :: faceDat(:,:)
    real(kind=rk), allocatable :: firstRefine(:,:), secondRefine(:,:)
    ! --------------------------------------------------------------------------

    ! The number of degrees of freedom on an element for the coarse, fine and intermediate level
    nDofsCoarse = scheme(currentLevel-1)%nDofs
    nDofsFine = scheme(currentLevel)%nDofs
    nDofsInter_firstRefine = &
                         &   (scheme(currentLevel)%modg_2d%maxPolyDegree+1) &
                         & * (scheme(currentLevel-1)%modg_2d%maxPolyDegree+1)
    nDofsInter_secondRefine = &
                         &   (scheme(currentLevel)%modg_2d%maxPolyDegree+1) &
                         & * (scheme(currentLevel)%modg_2d%maxPolyDegree+1)

    nfluids = mesh(currentLevel)%descriptor%elem%nElems(eT_fluid)

    ! Create the intermediate and resulting arrays
    allocate( faceDat(nDofsCoarse,nScalars) )
    allocate( firstRefine(nDofsInter_firstRefine,nScalars) )
    allocate( secondRefine(nDofsInter_secondRefine,nScalars) )

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
        &                                          %total(childPos)    )

      ! For 2D the if childpos would return childnum > 4 then it should
      ! just project the corresponding child from 1-4
      childNum = mod(childNum-1, 4) + 1

      ! Get position of the parent element
      parentPos = mesh(currentLevel)%faces_stab%dimByDimDesc(iDir)%depFromCoarser(iElem)%elem%val(1)

      if (parentPos <= nfluids) then
        totpos = mesh(currentLevel)%descriptor%pntTID(parentpos)
        ! start timer for LB of parent element
        call tem_startTimer( me          = atl_elemTimers, &
          &                  timerHandle = totPos          )
      end if

      ! Get state of the coarser element

      faceDat(:,:) = state_stab(currentLevel-1,iDir)%state(parentPos,:,:)

      ! Get the fine element shift from the child number
      iRefineZ = (childNum-1)/4 + 1
      iRefineY = (childNum-1-(iRefineZ-1)*4)/2+1
      iRefineX = mod(childNum-1,2)+1

      ! Do the first refinement step in x direction
      do ii = lbound(firstRefine,2),ubound(firstRefine,2)
        firstRefine(:, ii) = 0.0_rk
      end do
      call modg_semiRefineElem_2d(                         &
        & modalRepFace  = faceDat,                         &
        & modg_basis    = scheme(currentLevel)%modg_basis, &
        & schemeCoarse  = scheme(currentLevel-1)%modg_2d,  &
        & schemeFine    = scheme(currentLevel)%modg_2d,    &
        & refineDir     = 1,                               &
        & fineElemShift = iRefineX,                        &
        & modalRefined  = firstRefine(:,:)                 )

      ! Do the second refinement step in y direction
      do ii=lbound(secondRefine,2),ubound(secondRefine,2)
        secondRefine(:,ii) = 0.0_rk
      end do
      call modg_semiRefineElem_2d(                         &
        & modalRepFace  = firstRefine(:,:),                &
        & modg_basis    = scheme(currentLevel)%modg_basis, &
        & schemeCoarse  = scheme(currentLevel-1)%modg_2d,  &
        & schemeFine    = scheme(currentLevel)%modg_2d,    &
        & refineDir     = 2,                               &
        & fineElemShift = iRefineY,                        &
        & modalRefined  = secondRefine(:,:)                )

      ! Assign element state for the fine (child) element
      state_stab(currentLevel,iDir)%state(childPos,:,:) = secondRefine(:,:)

      if (parentPos <= nFluids) then
        call tem_stopTimer( me          = atl_elemTimers, &
          &                 timerHandle = totPos          )
      end if

    end do

  end subroutine atl_modg_coarseToFineElem_2d


  !> Subroutine to semi-refine an element with modal polynomial representation
  !! into its semi-children.
  subroutine modg_semiRefineElem_2d( modalRepFace, modg_basis, schemeCoarse, &
    &                                schemeFine, refineDir, fineElemShift,   &
    &                                modalRefined                            )
    ! --------------------------------------------------------------------------
    !> Modal representation of a function on the non-refined face.
    !! Dimensions are: (modg%maxPolyDegree+1)^2 for the first dimension
    !! and nScalars for the second dimension.
    real(kind=rk), intent(in) :: modalRepFace(:,:)
    !> Informations about the polynomial basis of the MODG scheme.
    type(ply_modg_basis_type) :: modg_basis
    !> The parameters of your MODG scheme on the coarse level.
    type(atl_modg_2d_scheme_type), intent(in) :: schemeCoarse
    !> The parameters of your MODG scheme on the fine level.
    type(atl_modg_2d_scheme_type), intent(in) :: schemeFine
    !> The direction of the semi-refinement. Either 1 or 2. Have a look at the
    !! function description.
    integer, intent(in) :: refineDir
    !> The semi-refined element you want to obtain.
    integer, intent(in) :: fineElemShift
    !> The modal representation of modalRepFace restricted to the semi-refined
    !! element.
    real(kind=rk), intent(inout) :: modalRefined(:,:)
    ! --------------------------------------------------------------------------
    integer :: iDegX, iDegY, iCoarseFunc
    integer :: dof, coarsePos
    real(kind=rk) :: fineSqNorm
    integer :: mpd1, mpd1_square
    ! --------------------------------------------------------------------------

    ! We treat the refinement directions separately to have the if condition
    ! outside of the loop.
    if( refineDir .eq. 1 ) then

      mpd1 = schemeCoarse%maxPolyDegree+1
      mpd1_square = mpd1 * (schemeFine%maxPolyDegree+1)

      ! loop over all degrees of freedoms on the (semi-)refined element
      do dof = 1, mpd1_square
        iDegX = mod(dof-1,mpd1)+1
        iDegY = (dof-1)/mpd1+1
        fineSqNorm = 2.0_rk /(2.0_rk * iDegX - 1.0_rk)
        ! ... loop over the ansatz functions of the current (non-refined) element
        do iCoarseFunc = 1, schemeCoarse%maxPolyDegree+1
  coarsepos = icoarsefunc                                      &
    &      + ( ( idegy-1)                             &
    &      + (1-1)*(schemecoarse%maxpolydegree+1))*(schemecoarse%maxpolydegree+1)
          ! Weight and add it to the semi refined face representation
          modalRefined(dof,:) = modalRefined(dof,:)                        &
            & + modg_basis%refineBaseCoeff                                 &
            &             %anz_anzShift(iDegX, iCoarseFunc, fineElemShift) &
            & * modalRepFace(coarsePos,:) / ( fineSqNorm )
        end do
      end do

    elseif( refineDir .eq. 2 ) then

      mpd1 = schemeFine%maxPolyDegree+1
      mpd1_square = mpd1 * (schemeFine%maxPolyDegree+1)

      ! loop over all degrees of freedoms on the (semi-)refined element
      do dof = 1, mpd1_square
        iDegX = mod(dof-1,mpd1)+1
        iDegY = (dof-1)/mpd1+1
        fineSqNorm = 2.0_rk /(2.0_rk * iDegY - 1.0_rk)
        ! ... loop over the ansatz functions of the current (non-refined) element
        do iCoarseFunc = 1, schemeCoarse%maxPolyDegree+1
  coarsepos = idegx                                      &
    &      + ( ( icoarsefunc-1)                             &
    &      + (1-1)*(schemecoarse%maxpolydegree+1))*(schemecoarse%maxpolydegree+1)
          ! Weight and add it to the semi refined face representation
          modalRefined(dof,:) = modalRefined(dof,:)                        &
            & + modg_basis%refineBaseCoeff                                 &
            &             %anz_anzShift(iDegY, iCoarseFunc, fineElemShift) &
            & * modalRepFace(coarsePos,:) / ( fineSqNorm )
        end do
      end do

    end if

  end subroutine modg_semiRefineElem_2d

end module atl_modg_2d_multilevel_module