! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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
!! Collects routines and datatypes related to face information.
module atl_facedata_module
  use env_module,                  only: rk
  use tem_faceData_module,         only: tem_face_type
  use tem_construction_module,     only: tem_levelDesc_type

  use atl_boundary_module,         only: atl_level_boundary_type

  implicit none
  private

  !> summary: Parameter to map left or right face of an element to its corresponding
  !! outer unit surface normal coefficient.
  !!
  !! Maps left or right face to its outer unit surface normal. Because of our
  !! cubical elements we do not need full 3-dimensional vectors here and we
  !! store only the non-zero components.
  !! Use \link tem_face_module:tem_leftface_prp \endlink and
  !! \link tem_face_module:tem_rightface_prp \endlink to access this mapping.
  real(kind=rk), parameter, dimension(2) :: atl_elemfaceToNormal_prp = &
                               & reshape((/ -1.0_rk, 1.0_rk /), (/ 2 /))

  !> Representation of a solution on a set of faces (all of them in
  !! a fixed spatial direction).
  type atl_faceRep_type
    !> The number of faces with a face representation.
    integer :: nFaces = 0
    !> Data for the face representation. The size of this array
    !! is: \n
    !! (nFaces, nFaceDofs, nScalars,2) -> 2 because of left and right face.
    real(kind=rk), allocatable :: dat(:,:,:,:)

  end type atl_faceRep_type

  !> Datatype to represent data defined on the faces of a cubical element.
  type atl_facedata_type
    !> Representation of the state on the face. We handle it for
    !! each spatial direction (x,y,z) separately.
    type(atl_faceRep_type) :: faceRep(3)
    !> Representation of the flux on the face. We handle it for
    !! each spatial direction (x,y,z) separately.
    type(atl_faceRep_type) :: faceFlux(3)
  end type atl_facedata_type

  interface atl_init_facedata
    module procedure atl_init_facedata_sym
    module procedure atl_init_facedata_asym
  end interface atl_init_facedata

  public :: atl_init_facedata, atl_facedata_type, atl_faceRep_type, &
          & atl_elemfaceToNormal_prp

contains

  !> summary: Initializes the face data by a given set of faces.
  !!
  !! Creates the datastructures for modal representations on the faces.
  !! This routine creates symmetric representations, i.e. state and fluxes have
  !! the same number of variables on the face.
  subroutine atl_init_facedata_sym( faces, facedata, minLevel, maxLevel, nDim, &
    &                               nScalars, nFaceDofs, boundary              )
    ! ---------------------------------------------------------------------------
    !> The minimum refinement level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum refinement level of the mesh.
    integer, intent(in) :: maxLevel
    !> The number of dimensions of the used scheme.
    integer, intent(in) :: nDim
    !> Description of the faces (levelwise).
    type(tem_face_type),intent(in)  :: faces(minLevel:maxLevel)
    !> The face data to be initialized (levelwise, running from minLevel to maxLevel).
    type(atl_facedata_type), allocatable, intent(out) :: facedata(:)
    !> The number of scalar variables to store on the faces.
    integer, intent(in) :: nScalars
    !> The number of degrees of freedoms on the face (one value for each level).
    integer, intent(in) :: nFaceDofs(minLevel:maxLevel)
    !> The boundary description for the faces on all levels of the mesh.
    type(atl_level_boundary_type), intent(in) :: boundary(minLevel:maxLevel)
    ! ---------------------------------------------------------------------------
    integer :: iLevel, iDir
    ! ---------------------------------------------------------------------------

    allocate(facedata(minLevel:maxLevel))

    ! We loop over all the level and all three spatial directions to initialize
    ! the face descriptors.
    levelLoop: do iLevel = minLevel, maxLevel

      directionLoop: do iDir = 1, nDim

        ! Intialize each of face representations.
        call atl_init_facerep(                               &
          & dimLevelDesc = faces(iLevel)%dimByDimDesc(iDir), &
          & nScalars     = nScalars,                         &
          & nFaceDofs    = nFaceDofs(iLevel),                &
          & boundary     = boundary(iLevel),                 &
          & dir          = iDir,                             &
          & faceRep      = facedata(iLevel)%faceRep(iDir)    )

        ! Intialize each of face fluxes.
        call atl_init_facerep(                               &
          & dimLevelDesc = faces(iLevel)%dimByDimDesc(iDir), &
          & nScalars     = nScalars,                         &
          & nFaceDofs    = nFaceDofs(iLevel),                &
          & boundary     = boundary(iLevel),                 &
          & dir          = iDir,                             &
          & faceRep      = facedata(iLevel)%faceFlux(iDir)   )

      end do directionLoop

    end do levelLoop

  end subroutine atl_init_facedata_sym



  !> summary: Initializes the face data by a given set of faces (asymmetrically).
  !!
  !! Creates the datastructures for modal representations on the faces.
  !! This routine creates asymmetric representations, i.e. state and fluxes have
  !! different number of variables on the face.
  subroutine atl_init_facedata_asym( faces, facedata, minLevel, maxLevel, &
    &                                nDim, nScalarsState, nScalarsFlux,   &
    &                                nFaceDofs, boundary                  )
    ! ---------------------------------------------------------------------------
    !> The minimum refinement level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum refinement level of the mesh.
    integer, intent(in) :: maxLevel
    !> The number of dimensions of the used scheme.
    integer, intent(in) :: nDim
    !> Description of the faces (levelwise).
    type(tem_face_type),intent(in)  :: faces(minLevel:maxLevel)
    !> The face data to be initialized (levelwise, running from minLevel to
    !! maxLevel).
    type(atl_facedata_type), allocatable, intent(out) :: facedata(:)
    !> The number of scalar variables to store on the faces for the state.
    integer, intent(in) :: nScalarsState
    !> The number of scalar variables to store on the faces for the flux.
    integer, intent(in) :: nScalarsFlux
    !> The number of degrees of freedoms on the face (one value for each level).
    integer, intent(in) :: nFaceDofs(minLevel:maxLevel)
    !> The boundary description for the faces on all levels of the mesh.
    type(atl_level_boundary_type), intent(in) :: boundary(minLevel:maxLevel)
    ! ---------------------------------------------------------------------------
    integer :: iLevel, iDir
    ! ---------------------------------------------------------------------------

    allocate(facedata(minLevel:maxLevel))

    ! We loop over all the level and all three spatial directions to initialize
    ! the face descriptors.
    levelLoop: do iLevel = minLevel, maxLevel

      directionLoop: do iDir = 1, nDim

        ! Intialize each of face representations.
        call atl_init_facerep(                               &
          & dimLevelDesc = faces(iLevel)%dimByDimDesc(iDir), &
          & nScalars     = nScalarsState,                    &
          & nFaceDofs    = nFaceDofs(iLevel),                &
          & boundary     = boundary(iLevel),                 &
          & dir          = iDir,                             &
          & faceRep      = facedata(iLevel)%faceRep(iDir)    )

        ! Intialize each of face fluxes.
        call atl_init_facerep(                               &
          & dimLevelDesc = faces(iLevel)%dimByDimDesc(iDir), &
          & nScalars     = nScalarsFlux,                     &
          & nFaceDofs    = nFaceDofs(iLevel),                &
          & boundary     = boundary(iLevel),                 &
          & dir          = iDir,                             &
          & faceRep      = facedata(iLevel)%faceFlux(iDir)   )

      end do directionLoop

    end do levelLoop

  end subroutine atl_init_facedata_asym

  !> summary: Initializes the face data by a given set of faces.
  subroutine atl_init_facerep( nScalars, nFaceDofs, faceRep, dimLevelDesc, &
    &                          dir, boundary )
    ! ---------------------------------------------------------------------------
    !> The number of scalar variables to store on the faces.
    integer, intent(in) :: nScalars
    !> The number of degrees of freedoms on the face.
    integer, intent(in) :: nFaceDofs
    !> Dimension-by-dimension level description of the mesh.
    type(tem_levelDesc_type), intent(in) :: dimLevelDesc
    !> The face representation to be initialized.
    type(atl_facerep_type), intent(out) :: faceRep
    !> The spatial direciton:\n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction
    integer, intent(in) :: dir
    !> The boundary description for the faces on the current level.
    type(atl_level_boundary_type), intent(in) :: boundary
    ! ---------------------------------------------------------------------------

    ! For each element in the dim-by-dim level descriptor we store face information.
    ! Additionally, we have to create some faces for the boundaries, so we add
    ! them here.
    ! ATTENTION:
    ! If we change the dimension here, we have
    ! to keep in mind that we have to adapt the
    ! parallel module, because this requires to
    ! serialize the dat arrays!
    faceRep%nFaces = dimLevelDesc%nElems &
                   & + sum( boundary%bnd(:)%faces(dir,1)%facePos%nVals ) &
                   & + sum( boundary%bnd(:)%faces(dir,2)%facePos%nVals )
    allocate( faceRep%dat(faceRep%nFaces, nFaceDofs, nScalars, 2) )
    faceRep%dat = 0.0_rk

  end subroutine atl_init_facerep

end module atl_facedata_module
