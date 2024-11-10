! Copyright (c) 2014-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014-2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> author: Jens Zudrop
!! This module provides central Ateles data type, containing the various data
!! of the simulation.
!!
!! The idea is to support different element types in the mesh, and collect
!! them in the container_module#element_container_type.
module atl_cube_container_module

  use flu_binding,                  only: flu_state

  use treelmesh_module,             only: treelmesh_type
  use tem_tools_module,             only: tem_horizontalSpacer
  use tem_logging_module,           only: logUnit
  use tem_bc_prop_module,           only: tem_bc_prop_type
  use tem_construction_module,      only: tem_find_allElements,             &
    &                                     tem_init_elemLevels,              &
    &                                     tem_build_VerticalDependencies,   &
    &                                     tem_build_HorizontalDependencies, &
    &                                     tem_cleanupDependencyArrays,      &
    &                                     tem_levelDesc_type
  use tem_stencil_module,           only: tem_stencilHeader_type, &
    &                                     tem_stencil_dump
  use tem_faceData_module,          only: tem_face_type
  use tem_face_module,              only: tem_build_face_info
  use tem_comm_module,              only: tem_commPattern_type
  use tem_comm_env_module,          only: tem_comm_env_type
  use tem_bc_header_module,         only: tem_bc_header_type

  use atl_boundary_module,          only: atl_init_bndList,  &
    &                                     atl_init_elem_bnd, &
    &                                     atl_level_boundary_type
  use atl_source_types_module,      only: atl_source_type
  use atl_cube_elem_module,         only: atl_cube_elem_type, &
    &                                     atl_init_cube_elem
  use atl_equation_module,          only: atl_equations_type
  use atl_facedata_module,          only: atl_facedata_type
  use atl_kerneldata_module,        only: atl_kerneldata_type, &
    &                                     atl_statedata_type
  use atl_scheme_module,            only: atl_scheme_type, &
    &                                     atl_define_SchemeStencil
  use atl_bc_header_module,         only: atl_boundary_type
  use atl_materialPrp_module,       only: atl_material_type
  use atl_penalization_module,      only: atl_penalizationData_type

  implicit none

  private

  !> Container type collecting all the data of the simulation domain which
  !! have cubic cells.
  !!
  !! It stores informations related to the mesh as well as the states of the
  !! kernels like state variables, scheme type on each part, etc. Please
  !! notice that you can have different schemes in different parts of the cubic
  !! mesh.
  type atl_cube_container_type
    integer :: minLevel !< Minimal level in tree
    integer :: maxLevel !< Maximal level in tree
    !> the number of list we have in this container. This is the number of
    !! levels with nonzero number of elements (i.e. maxLevel - minLevel)
    integer :: nlists = 0

    !> List of the solver state. This includes information about the timepoints
    !! and the description of the equation state (e.g. in model representation).
    !! Indices are running from minLevel to maxLevel.
    type(atl_statedata_type), allocatable :: statedata_list(:)

    !> List of solver state on the faces. Indices are running from minLevel to
    !! maxLevel.
    type(atl_facedata_type), allocatable :: facedata_list(:)

    !> List of boundaries for the current part of the cubic mesh its
    !! stabilization. The boundary informations are stored for each refinement
    !! level seperately. Therefore the dimension of this array is nlists.
    !! Indices are running from minLevel to maxLevel.
    type(atl_level_boundary_type), allocatable :: boundary_stab_list(:)
    !> Statedata of the stabilization.
    !! Indices are running from minLevel to maxLevel.
    type(atl_statedata_type), allocatable :: statedata_stab_list(:,:)

    !> List of the states of the kernels. This includes data that is used in the
    !! kernel itself (e.g. for reconstruction, etc.). It might be necessary
    !! to transfer the data during each computation of the right hand side.
    !! Indices are running from minLevel to maxLevel.
    type(atl_kerneldata_type), allocatable :: kerneldata_list(:)

    !> List of penalization data. This includes all data that is used for
    !! penalizations (e.g. Brinkmann penalizations of the Navier-Stokes
    !! equations). Indices are running from minLevel to maxLevel.
    type(atl_penalizationData_type), allocatable :: penalizationdata_list(:)

    !> List of meshes of the different kernels. We can have multiple kernels
    !! e.g. for each refinement level, polynomial degree, etc. . The states
    !! of the kernels of mesh_list(index) are located inside
    !! kerneldata_list(index). This separation was done to enable a generic
    !! kernel interface. Indices are running from minLevel to maxLevel.
    type(atl_cube_elem_type), allocatable :: mesh_list(:)

    !> List of boundaries for the current part of the cubic mesh. The boundary
    !! informations are stored for each refinement level seperately. Therefore
    !! the dimension of this array is nlists. Indices are running from minLevel
    !! to maxLevel.
    type(atl_level_boundary_type), allocatable :: boundary_list(:)
    !> Global description of the boundary conditions. Size is the number of
    !! boundary conditions.
    type(atl_boundary_type), allocatable :: bc(:)

    type(atl_source_type) :: source

    !> List of material parameter information for the mesh. One entry for
    !! level, running from minlevel to maxlevel.
    type(atl_material_type), allocatable :: material_list(:)
    type(atl_material_type), allocatable :: penalization_list(:)

    !> Pointer from treeIDlist entry to level-wise fluid part of total list.
    !! The length of this vector is the total number of cubic elements.
    integer, allocatable :: levelPointer(:)

    !> All informations about the scheme for each level, e.g. for PNPM we store
    !! in this variable n and m and many others.... For a detailed
    !! description you should have a look at
    !! the documentation of the atl_scheme_type. Indices are runnging from
    !! minLevel to maxLevel.
    type(atl_scheme_type), allocatable :: scheme_list(:)

    !> Gives the position of the levelwise projection method in the unique list
    !! of projections.
    integer, allocatable :: poly_proj_pos(:)

    !> Indication whether to compute the maximal deviation in the
    !! polynomials of each element.
    logical :: need_element_deviations = .false.

  end type atl_cube_container_type

  public :: atl_cube_container_type, atl_init_cube_container

contains


  ! ****************************************************************************
  !> Method to initialize a cube mesh by tree and boundary
  !! definitions obtained by treelm.
  subroutine atl_init_cube_container( tree, boundary, cube_container, conf, &
    &                                 equation, proc, commPattern,          &
    &                                 need_element_deviations               )
    ! --------------------------------------------------------------------------
    !> The container with all cubic elements, generated with this routine.
    type(atl_cube_container_type),intent(inout) :: cube_container

    !> The tree representation of your mesh.
    type(treelmesh_type), intent(inout) :: tree

    !> The boundaries of your simulation domain
    type(tem_bc_prop_type), intent(in) :: boundary

    !> Handle for the Lua config file
    type(flu_State), intent(in) :: conf

    !> The equation you are simulating.
    type(atl_equations_type),intent(inout) :: equation

    !> mpi communication environment with mpi communicator
    type(tem_comm_env_type) :: proc

    !> mpi communication pattern type
    type(tem_commPattern_type) :: commPattern

    logical, intent(in) :: need_element_deviations

    ! --------------------------------------------------------------------------
    integer                               :: iLevel
    type(tem_levelDesc_type), allocatable :: levelDesc(:)
    type(tem_face_type), allocatable      :: faces(:), faces_stab(:)
    type(tem_bc_header_type)              :: bc_header
    type(tem_stencilHeader_type)          :: stencil(1)
    ! --------------------------------------------------------------------------

    cube_container%need_element_deviations = need_element_deviations

    ! Now we can build the list of neighbors.
    ! Then we use the generic tem_construction_module
    ! of treelm to get the right connections and ghost and halos.

    !> now we build the face information of our mesh.
    call tem_horizontalSpacer(fUnit=logUnit(2))
    write(logUnit(2),*) 'Creating face description ... '
    call tem_build_face_info( tree              = tree,        &
      &                       boundary          = boundary,    &
      &                       commPattern       = commPattern, &
      &                       proc              = proc,        &
      &                       faces             = faces,       &
      &                       nEligibleChildren = 4            )
    write(logUnit(2),*) 'Creating face description (for stabilization)... '
    call tem_build_face_info( tree              = tree,        &
      &                       boundary          = boundary,    &
      &                       commPattern       = commPattern, &
      &                       proc              = proc,        &
      &                       faces             = faces_stab,  &
      &                       nEligibleChildren = 8            )

    ! Additional elements that might be required by the scheme.
    !!VK stencil is not required for modg, but for tracking
    call atl_define_schemeStencil( nDims = equation%nDimensions, &
      &                            me    = stencil(1)            )
    call tem_stencil_dump(stencil(1))

    ! Find the treeIDs for all offsets (stencil elements) for all
    ! elements. Direct neighbors, always required.
    write(logUnit(3),*) 'tem_init_elemLevels'
    call tem_init_elemLevels(  me       = levelDesc, &
      &                        boundary = boundary,  &
      &                        tree     = tree,      &
      &                        stencils = stencil    )

    call tem_find_allElements( tree           = tree,                        &
      &                        levelDesc      = levelDesc,                   &
      &                        levelPointer   = cube_container%levelPointer, &
      &                        computeStencil = stencil,                     &
      &                        commpattern    = commpattern,                 &
      &                        proc           = proc                         )

    ! For each of the neighbor lists create the horizontal relations
    write(logUnit(3),*) 'tem_build_HorizontalDependencies'
    call tem_build_HorizontalDependencies( iStencil       = 1,         &
      &                                    levelDesc      = levelDesc, &
      &                                    tree           = tree,      &
      &                                    computeStencil = stencil(1) )
    write(logUnit(3),*) 'tem_build_VerticalDependencies'
    call tem_build_VerticalDependencies( levelDesc = levelDesc,            &
      &                                  minLevel  = tree%global%minLevel, &
      &                                  maxLevel  = tree%global%maxLevel  )

    ! The dependencay arrays are used in buildFaceInformation subroutine
    ! so, cleaning up the arrays must be done after buildFaceInformation!
    call tem_cleanupDependencyArrays( levelDesc = levelDesc )

    ! Allocate memory for cubical elements.
    cube_container%minLevel = tree%global%minLevel
    cube_container%maxLevel = tree%global%maxLevel
    cube_container%nLists = tree%global%maxLevel - tree%global%minLevel + 1

    allocate(cube_container%mesh_list( &
      & tree%global%minLevel:tree%global%maxLevel))
    allocate(cube_container%boundary_list( &
      & tree%global%minLevel:tree%global%maxLevel))
    allocate(cube_container%boundary_stab_list( &
      & tree%global%minLevel:tree%global%maxLevel))
    allocate(cube_container%kerneldata_list( &
      & tree%global%minLevel:tree%global%maxLevel))
    allocate(cube_container%statedata_list( &
      & tree%global%minLevel:tree%global%maxLevel))
    allocate(cube_container%statedata_stab_list( &
      & tree%global%minLevel:tree%global%maxLevel,3))
    allocate(cube_container%material_list( &
      & tree%global%minLevel:tree%global%maxLevel))
    allocate(cube_container%penalizationdata_list( &
      & tree%global%minLevel:tree%global%maxLevel))

    ! Allocate the source data type depending on levels
    allocate(cube_container%source%poly_proj_pos( &
      & tree%global%minLevel:tree%global%maxLevel))

    ! Iterate over all levels and initialize the container for it.
    do iLevel = tree%global%minLevel, tree%global%maxLevel
      call atl_init_cube_elem( element    = cube_container%mesh_list(iLevel), &
        &                      descriptor = levelDesc(iLevel),                &
        &                      level      = iLevel,                           &
        &                      tree       = tree                              )
    end do

    ! Create the levelwise list of boundaries.
    call atl_init_bndList(                            &
      & conf          = conf,                         &
      & tree          = tree,                         &
      & equation      = equation,                     &
      & bc_prop       = boundary,                     &
      & face_list     = faces,                        &
      & boundary_list = cube_container%boundary_list, &
      & bc            = cube_container%bc,            &
      & bc_header     = bc_header,                    &
      & scheme_list   = cube_container%scheme_list    )

    ! create the boundary informations for the dim-by-dim description
    ! for co-volume
    call atl_init_elem_bnd( minLevel      = tree%global%minLevel,              &
      &                     maxLevel      = tree%global%maxLevel,              &
      &                     bc_header     = bc_header,                         &
      &                     face_list     = faces_stab,                        &
      &                     boundary_list = cube_container%boundary_stab_list, &
      &                     scheme_list   = cube_container%scheme_list         )

    do iLevel = tree%global%minLevel, tree%global%maxLevel
      ! This component needs to be explicitely copied as a workaround for
      ! the cray compiler
      cube_container%mesh_list(ilevel)%faces = faces(ilevel)
      cube_container%mesh_list(ilevel)%faces_stab = faces_stab(ilevel)
    end do

  end subroutine atl_init_cube_container
! ******************************************************************************!


end module atl_cube_container_module
