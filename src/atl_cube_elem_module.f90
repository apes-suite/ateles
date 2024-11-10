! Copyright (c) 2011-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2015, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!! Container collecting the subroutines and datatypes which are specific
!! for cubic elements.
module atl_cube_elem_module
  use env_module,              only: rk, long_k

  use tem_geometry_module,     only: tem_baryOfId
  use tem_construction_module, only: tem_levelDesc_type
  use treelmesh_module,        only: treelmesh_type
  use tem_faceData_module,     only: tem_face_type
  use tem_logging_module,      only: logUnit

  implicit none
  private

  !> Container type describing cubic elements on a single refinement level.
  !!
  !! It contains the number of elements as well as neighbor relations and
  !! specific informations of cubic elements on this level like jacobi matrix,
  !! properties of jacobi matrix, length, coordinates of the cubic cells,
  !! area of the sides, etc.
  type atl_cube_elem_type

    !> Description of the faces on this level.
    !!
    !! This description includes list of compute faces, interpolation faces
    !! and communicated faces. Furthermore it holds the relevant buffers
    !! for the interprocess communication with respect to the faces.
    type(tem_face_type) :: faces

    !> Description of the faces (required for stabilization) on this level.
    !!
    !! This description includes list of compute faces, interpolation faces
    !! and communicated faces. Furthermore it holds the relevant buffers
    !! for the interprocess communication with respect to the faces.
    !! The datatype is initialized properly in case a stabilization
    !! requires dimension-by-dimension neighbor information.
    type(tem_face_type) :: faces_stab

    !> Descriptor describing the grid on this refinement level.
    !!
    !! This variable contains all information about neighbor relations, as well
    !! as ghost and halo cells. It also
    !! contains informations about the number of cells on this level and so on.
    type(tem_levelDesc_type) :: descriptor

    !> Length of the cubical element.
    real(kind=rk) :: length

    !> Surface area of the individual sides.
    real(kind=rk) :: side_area

    !> Volume of each element.
    real(kind=rk) :: volume

    !> Bary center coordinate of each element.
    !!
    !! The dimension of the array is: (numFluid+numGhosts+numHalos+numBnds,3).
    !! It is used in stecil mapping for example.
    real(kind=rk), allocatable :: bary_coord(:,:)

    !> Inverse jacobian \f$\partial\xi_i/\partial x_j\f$ (the jacobian is the
    !! jacobi matrix of the mapping from the reference to the physical cell).
    real(kind=rk) :: inv_jacobit

    !> Determinant of the jacobian.
    !!
    !! The jacobian is the jacobi matrix of the mapping from the reference to
    !! the physical cell).
    !! Since we consider only cubic elements the jacobian is a constant
    !! and therefore its determinant is also constant.
    real(kind=rk)  :: jacobit_det
  end type atl_cube_elem_type

  public :: atl_cube_elem_type, atl_init_cube_elem
  public :: atl_get_numberOfElemsPerLevel

contains

  !> Initialize the cubic elements.
  !!
  !! The output element represents all cubic elements given by their tree ids
  !! in the tree.
  subroutine atl_init_cube_elem( element, descriptor, level, tree )
    ! --------------------------------------------------------------------------
    !> This is the output and represenets the cubic elements given
    !! as a subset of tree ids in the complete tree.
    type(atl_cube_elem_type), intent(out) :: element

    !> The descriptor of the element list, describing the connectivity
    !! of the mesh explicitly.
    !!
    !! All the descriptors are derived previously from the complete
    !! treelmesh, and need to be passed in here.
    type(tem_levelDesc_type), intent(in) :: descriptor

    !> The tree representation of your mesh.
    type(treelmesh_type), intent(in) :: tree

    !> Treelm level of the cubes to be initialized.
    integer, intent(in) :: level
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    write(logUnit(1),*) "Initializing Cubic elements ..."

    ! Copy of descriptor information, HK: really necessary?
    element%descriptor = descriptor

    ! Simple geometrical info
    element%length = tree%global%BoundingCubeLength / 2**level
    element%volume = element%length**3
    element%side_area = element%length**2

    ! Provide some information to the user
    write(logUnit(1),*) "--"
    write(logUnit(1),*) "Element information:"
    write(logUnit(1),*) "* Level: ", level
    write(logUnit(1),*) "* Length: ", element%length
    write(logUnit(1),*) "* Surface: ", element%side_area
    write(logUnit(1),*) "* Volume: ", element%volume
    write(logUnit(1),*) "--"

    ! Allocate bary_coord for all elements
    allocate(element%bary_coord(size(descriptor%total, 1), 3))

    ! Calculate physical barycentric coordinates for all cells
    ! we store data in the state vector, i.e. fluid, ghost, halos.
    call calc_barycoord(treeids    = descriptor%total,          &
      &                 nElems     = size(descriptor%total, 1), &
      &                 bary_coord = element%bary_coord ,       &
      &                 tree       = tree                       )


    ! constant jacobians
    element%inv_jacobit = 1.0_rk / element%length
    element%jacobit_det = element%volume

  end subroutine atl_init_cube_elem

  !> Calculate barycentric coordinates of the tree ids given in treeids.
  subroutine calc_barycoord( treeids, nElems, bary_coord, tree )
    ! --------------------------------------------------------------------------
    !> Tree ids you want to build the barycentric
    !! coordinates for.
    integer(kind=long_k), intent(in) :: treeids(:)
    !> the number of tree ids in in treeids.
    integer, intent(in)              :: nElems
    !> Array of barycenteric coordinates for the cells given
    !! in treeidsubset. This array has to be allocated before you can use it
    !! as an input argument. The dimensions are: First dimension is subsetsize
    !! and the second dimension is the space dimension, i.e. 3.
    real(kind=rk), intent(out)       :: bary_coord(:,:)
    !> Tree representation of your mesh.
    type(treelmesh_type),intent(in) :: tree
    ! --------------------------------------------------------------------------
    integer                          ::  ielem
    ! --------------------------------------------------------------------------
    do ielem = 1, nElems
      bary_coord(iElem, :)  = tem_baryOfId(tree, treeids(iElem))
    enddo
  end subroutine calc_barycoord

  !> Subroutine to count the number of elements per level.
  !!
  !! This count includes fluid, ghost and halo elements on each level.
  subroutine atl_get_numberOfElemsPerLevel( descriptor, nCells, tree )
    ! --------------------------------------------------------------------------
    !> The tree representation of the mesh
    type(treelmesh_type), intent(in) :: tree

    !> The number of cells for each levelDescriptor (including fluid, ghost
    !! and halo cells).
    integer, allocatable, intent(out) :: nCells(:)

    !> Array of descriptors for each level in the mesh.
    type(tem_levelDesc_type),intent(in) :: descriptor(tree%global%minLevel &
      &                                                  :tree%global%maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iList
    ! --------------------------------------------------------------------------

    allocate(nCells(tree%global%minLevel:tree%global%maxLevel))

    do iList = tree%global%minLevel, tree%global%maxLevel
      nCells(iList) =   descriptor(iList)%nElems
    end do
  end subroutine atl_get_numberOfElemsPerLevel

end module atl_cube_elem_module
