! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Verena Krupp, Peter Vitt,
! Tobias Girresser, Nikhil Anand and Neda Ebrahimi Pour for University of
! Siegen.
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
!> Description of point sets.
module ply_nodes_module
  use env_module, only: rk, labelLen
  use fftw_wrap, only: fftw_available
  use aotus_module, only: flu_State, aot_get_val

  use tem_aux_module, only: tem_abort

  use ply_nodeset_module, only: ply_nodeset_coords,    &
    &                           ply_nodeset_legendre,  &
    &                           ply_nodeset_chebyshev, &
    &                           ply_nodeset_chebyloba

  use ply_nodes_header_module, only: ply_nodes_header_type, assignment(=)

  implicit none

  private

  !> Datatype to represent facewise nodes
  type ply_faceNodes_type
    !> The number of face nodes
    integer :: nquadPoints
    !> The face nodes.
    !! First index goes from 1 to nPoints and second index
    !! from 1 to 3 for the 3 spatial coordinates.
    real(kind=rk), allocatable :: points(:,:)
  end type ply_faceNodes_type

  public :: ply_nodes_create
  public :: ply_faceNodes_type


contains


  ! ------------------------------------------------------------------------ !
  !> Initialize points with the Chebyshev quadrature points, 3D
  subroutine ply_nodes_create(me, nodes, faces, nQuadPointsPerDir, ndims )
    ! -------------------------------------------------------------------- !
    type(ply_nodes_header_type), intent(in) :: me
    real(kind=rk), allocatable, intent(out)  :: nodes(:,:)
    type(ply_faceNodes_type), allocatable, intent(out)  :: faces(:,:)
    integer, intent(in) :: nQuadPointsPerDir
    integer, intent(in) :: ndims
    ! -------------------------------------------------------------------- !
    procedure(ply_nodeset_coords), pointer :: nodeset => NULL()
    integer :: iDir
    ! -------------------------------------------------------------------- !

    select case(trim(me%nodes_kind))
    case('gauss-legendre')
      nodeset => ply_nodeset_legendre
    case('chebyshev')
      if (me%lobattoPoints) then
        nodeset => ply_nodeset_chebyloba
      else
        nodeset => ply_nodeset_chebyshev
      end if
    end select

    ! Build the Chebyshev nodes on the reference element (volume)
    call ply_nodes_volume_coords(                      &
      &    num_intp_per_direction = nQuadPointsPerDir, &
      &    nDims                  = nDims,             &
      &    nodeset                = nodeset,           &
      &    points                 = nodes              )

    ! Build the Chebyshev nodes on the reference faces
    allocate( faces(ndims,2) )
    faces(:,:)%nquadpoints = nQuadPointsPerDir**(nDims-1)

    do iDir=1,ndims
      call ply_nodes_surface_coords(                         &
        &    num_intp_per_direction = nQuadPointsPerDir,     &
        &    nDims                  = ndims,                 &
        &    nodeset                = nodeset,               &
        &    left                   = faces(iDir,1)%points,  &
        &    right                  = faces(iDir,2)%points,  &
        &    dir                    = iDir                   )
    end do

  end subroutine ply_nodes_create
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create multidimensional points from given 1D set of nodes in the cubic
  !! reference element.
  !!
  !! The points will be created by a tensor product of the provided 1d nodeset
  !! for the given number of dimensions.
  subroutine ply_nodes_volume_coords( num_intp_per_direction, &
    &                                 nDims, nodeset, points  )
    ! -------------------------------------------------------------------- !
    !> Number auf integration points in each direction.
    integer, intent(in) :: num_intp_per_direction

    !> Number of dimensions to create the points for.
    integer, intent(in) :: nDims

    !> Set of node coordinates to use in the element.
    procedure(ply_nodeset_coords) :: nodeset

    !> Resulting list of points. First index runs over all points, second
    !! indicates the coordinate dimension (x=1,y=2,z=3).
    !!
    !! For ndims smaller than 3, the higher dimensions will be set to 0.
    real(kind=rk), allocatable, intent(out) :: points(:,:)
    ! -------------------------------------------------------------------- !
    integer :: n1d
    integer :: numQuadPoints
    ! -------------------------------------------------------------------- !

    n1d = num_intp_per_direction
    numQuadPoints = n1d**nDims

    allocate(points(numQuadPoints,3))

    points = 0.0_rk
    call ply_point_tensor( nPoints1D = n1d,             &
      &                    nDims     = nDims,           &
      &                    nodeset   = nodeset,         &
      &                    points    = points(:,:nDims) )

  end subroutine ply_nodes_volume_coords
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create the integration points on the surface of (cubical) elements.
  subroutine ply_nodes_surface_coords(                 &
    &          num_intp_per_direction, ndims, nodeset, &
    &          left, right, dir                        )
    ! -------------------------------------------------------------------- !
    !> Number of integration points in each direction
    integer, intent(in) :: num_intp_per_direction

    !> Number of dimensions in the element.
    integer, intent(in) :: ndims

    !> Set of node coordinates to use in the element for which the surface
    !! points are to be defined.
    procedure(ply_nodeset_coords) :: nodeset

    !> The points on the left surface.
    real(kind=rk), allocatable, intent(out) :: left(:,:)

    !> The points on the right surface.
    real(kind=rk), allocatable, intent(out) :: right(:,:)

    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction
    integer :: dir
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: surface(:,:)
    integer :: nquadPoints
    integer :: n1d
    integer :: tangents(2,3)
    ! -------------------------------------------------------------------- !

    tangents(:,1) = [2,3]
    tangents(:,2) = [1,3]
    tangents(:,3) = [1,2]

    n1d = num_intp_per_direction
    nQuadPoints = n1d**(ndims-1)
    allocate(left(nQuadPoints, 3))
    allocate(right(nQuadPoints, 3))
    ! Ensure that the higher dimensions are set to 0,
    ! if not all three dimensions are used.
    left = 0.0_rk
    right = 0.0_rk

    left(:,dir)  = -1.0_rk
    right(:,dir) =  1.0_rk

    if (nDims >= 2) then
      allocate(surface(nQuadPoints, nDims-1))

      call ply_point_tensor( nPoints1D = n1d,     &
        &                    nDims     = nDims-1, &
        &                    nodeset   = nodeset, &
        &                    points    = surface  )
      left(:, tangents(:nDims-1, dir))  = surface
      right(:, tangents(:nDims-1, dir)) = surface
    end if

  end subroutine ply_nodes_surface_coords
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Compute a multi-dimensional tensor for the given set of nodes.
  subroutine ply_point_tensor( nPoints1D, nDims, nodeset, points  )
    ! -------------------------------------------------------------------- !
    !> Number auf integration points in each direction.
    integer, intent(in) :: nPoints1D

    !> Number of dimensions to create the points for.
    integer, intent(in) :: nDims

    !> Set of node coordinates to use in the element.
    procedure(ply_nodeset_coords) :: nodeset

    !> Resulting list of points. First index runs over all points, second
    !! indicates the coordinate dimension (x=1,y=2,z=3).
    real(kind=rk), intent(out) :: points(nPoints1D**nDims, nDims)
    ! -------------------------------------------------------------------- !
    integer :: j, k
    integer :: jk
    integer :: n1d
    integer :: nPlane
    real(kind=rk), allocatable :: gaussp1D(:)
    integer :: numQuadPoints
    ! -------------------------------------------------------------------- !

    n1d = nPoints1D
    numQuadPoints = n1d**nDims
    nPlane = n1d**(nDims-1)

    if (numQuadPoints > 1) then
      ! Only if there are multiple integration points to compute look into
      ! the nodeset.
      allocate(gaussp1D(n1d))

      gaussp1D = nodeset(n1d)

      ! X coordinates
      do jk=1,nPlane
        points((jk-1)*n1d+1:jk*n1d, 1) = gaussp1D
      end do

      if (nDims >= 2) then

        ! Y coordinates
        do jk=1,nPlane
          j = mod(jk-1, n1D) + 1
          points((jk-1)*n1d+1:jk*n1d, 2) = gaussp1D(j)
        end do

        if (nDims >= 3) then

          ! Y coordinates
          do k=1,n1D
            points((k-1)*nPlane+1:k*nPlane, 3) = gaussp1D(k)
          end do

        end if

      end if

    else
      ! Just a single integration point, return the center.
      points = 0.0_rk
    end if

  end subroutine ply_point_tensor
  ! ------------------------------------------------------------------------ !

end module ply_nodes_module
