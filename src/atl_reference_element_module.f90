! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
!! Module for operations relating reference and physical elements.
!!
!! This module collects all the functions, subroutines related to operations
!! on the reference element or mappings between the reference and the physical
!! elements.
module atl_reference_element_module
  ! Treelm modules
  use env_module, only: rk

  use treelmesh_module, only: treelmesh_type
  use tem_topology_module, only: tem_coordOfID
  use tem_geometry_module, only: tem_elemSizeLevel

  implicit none
  private

  public :: atl_refToPhysCoord, atl_ref_in_elempos


contains


  !> Subroutine to move points defined on the reference element [-1,+1] to
  !! the physical element coordinates.
  subroutine atl_refToPhysCoord( refPoints, nPoints, baryCoord, elemLength, &
    &                            physPoints                                 )
    !> The points of the reference elements you want to move to their physical
    !! location.
    !! Size is: nPoints for first dimension and 3 for the second dimension (x,y,z
    !! coordinates).
    real(kind=rk), intent(in) :: refPoints(:,:)
    !> The number of points to move.
    integer, intent(in) :: nPoints
    !> The barycentric coordinate of the element you want to map to (x,y,z coord)
    real(kind=rk), intent(in) :: baryCoord(3)
    !> The length of the cubic element to map to.
    real(kind=rk), intent(in) :: elemLength
    !> The physical points. Dimensions are the same as refPoints.
    real(kind=rk), intent(inout) :: physPoints(nPoints, 3)
    ! ---------------------------------------------------------------------------
    integer :: iDir
    ! ---------------------------------------------------------------------------


    ! Now, we iterate over the points and move each of them to the right position.
    ! To do so, we calculate the affin mapping from reference to physical
    ! element in advance.
    do iDir = 1,3
      physPoints(:,iDir) = baryCoord(iDir) + (elemLength/2.0_rk) * refPoints(:,iDir)
    end do


  end subroutine atl_refToPhysCoord


  ! ------------------------------------------------------------------------ !
  !> Transform reference points to physical points in the element
  !! of the tree identified by the provided elempos.
  subroutine atl_ref_in_elempos(refPoints, tree, elempos, physPoints)
    ! -------------------------------------------------------------------- !
    !> Reference points to transform.
    !!
    !! They have the form (npoints,3).
    real(kind=rk), intent(in) :: refPoints(:,:)

    !> Tree in which the element is found.
    type(treelmesh_type), intent(in) :: tree

    !> Position of the element in the list of treeIDs of the tree.
    integer, intent(in) :: elempos

    !> Transformed points.
    !!
    !! Has the same form as refPoints.
    real(kind=rk), intent(out) :: physPoints(:,:)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: bary(3)
    real(kind=rk) :: elemlen
    integer :: coord(4)
    integer :: iDir
    ! -------------------------------------------------------------------- !

    coord = tem_coordOfId(tree%treeID(elempos))
    elemlen = tem_elemSizeLevel(tree, coord(4))

    bary = tree%global%origin + (real(coord(:3), kind=rk) + 0.5_rk)*elemlen

    do iDir = 1,3
      physPoints(:,iDir) = bary(iDir) + (0.5_rk*elemLen) * refPoints(:,iDir)
    end do

  end subroutine atl_ref_in_elempos
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module atl_reference_element_module
