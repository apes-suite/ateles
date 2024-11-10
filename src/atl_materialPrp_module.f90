! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
!! Module for a generic description of materials in ATELES.
module atl_materialPrp_module

  ! Treelm modules
  use env_module,                     only: rk
  use tem_faceData_module,            only: tem_faceIterator_type
  use tem_dyn_array_module,           only: dyn_intArray_type
  use tem_spacetime_fun_module,       only: tem_spacetime_fun_type
  use tem_varMap_module,              only: tem_possible_variable_type
  use tem_stringKeyValuePair_module,  only: grw_stringKeyValuePairArray_type

  ! Ateles modules
  use atl_boundary_module,            only: atl_level_boundary_type

  implicit none

  private

  !> Parameter for a face with purely constant material parameters.
  integer, parameter :: atl_pureConstMat_prp = 1

  !> Parameter for a face with mixed materials (constant-variable or
  !! variable-constant or variable-variable face).
  integer, parameter :: atl_mixedMat_prp = 2

  !> The index used by the atl_material_property_type%computeElems to address
  !! the array with elements that have constant material properties.
  integer, parameter :: atl_ConstMatIdx = 1
  !> The index used by the atl_material_property_type%computeElems to address
  !! the array with elements that have variable material properties.
  integer, parameter :: atl_VarMatIdx = 2

  !> Encapsulates a pointer to an tem_spacetime_fun_type-instace. Used to create
  !! arrays of pointers.
  type atl_spacetime_fun_pointer_type
    !> The position of the material variable in the global varSys.
    integer :: matVarPos
    !> The position of the space-time-function stFunPtr is pointing to in the
    !! material variable's stFunList.
    integer :: stFunPos
    !> Pointer to an instance of tem_spacetime_fun_type.
    type(tem_spacetime_fun_type), pointer :: stFunPtr
  end type atl_spacetime_fun_pointer_type

  !> Description of the material properties adjacent to a faces.
  type atl_face_material_type
    !> The material positions for all the faces of the face desciption for
    !! left and right element. \n
    !! Therefore the dimensions are:
    !! \li First dimension is the number of faces
    !! \li Second dimension is the 2 for left and right element.
    !! \sa tem_left
    !! \sa tem_right
    !! The third index is do determine between the different material
    !! paramters.
    type(atl_spacetime_fun_pointer_type), allocatable :: mat(:,:,:)
  end type atl_face_material_type

  !> List of elements that are relevant for the computation.
  type atl_computeElem_type
    !> The number of elements in this container
    integer :: nElems = 0
    !> Fluid element indices in the total list of level descriptor
    integer, allocatable :: totElemIndices(:)
  end type atl_computeElem_type

  !> Boundary (face) material description.
  type atl_boundaryMaterial_type
    !> Boundary face description.
    type(atl_level_boundary_type) :: boundary
  end type atl_boundaryMaterial_type

  !> Levewise description of the material properties.
  type atl_material_property_type
    !> The material property position for all fluid elements.
    !! The first index is the element ordering, which is the same as in the
    !! total list of the level descriptor.\n
    !! The second index is the material parameter, as they are defined
    !! independently of each other.
    type(atl_spacetime_fun_pointer_type), allocatable :: material_elems(:,:)
    !> Information for the faces and the materials meeting at these interfaces,
    !! for each spatial direction (i.e. x y and z one entry).
    type(atl_face_material_type) :: material_face(3)
    !> The compute list for the material-face-combinations.
    !! First dimension is is number of spatial directions (i.e. x, y, z).
    !! Second dimension is 2 for the two different material face combinations
    !! \sq atl_pureConstMat_prp
    !! \sa atl_mixedMat_prp
    type(tem_faceIterator_type) :: computeFace(3,2)
    !> Create a list with element indices with constant material properties and
    !! those with variable ones. First entry covers all the elements with
    !! constant material paramters, while the second entry covers all the
    !! elements with non-constant material parameters.
    !! \sa atl_ConstMatIdx
    !! \sa atl_VarMatIdx
    type(atl_computeElem_type) :: computeElems(2)
    !> The boundary face information.
    !! First entry on first index covers all boundary faces with constant
    !! material parameters. The second entry on first index covers all the
    !! boundary faces with non-constant material parameters.
    !! \sa atl_ConstMatIdx
    !! \sa atl_VarMatIdx
    type(atl_boundaryMaterial_type) :: bnd_faces(2)
  end type atl_material_property_type

  !> Datatype provides information each
  type atl_elemMaterialData_type
    !> The number of nodal points with material information per spatial
    !! direction.
    integer :: nPointsPerDir
    !> The nodal material properties. \n
    !! First dimension is the number of elements with this material info. \n
    !! Second dimension is the number of nPointsPerDir^3 \n
    !! Third dimension is the number of material components. Attention: We do
    !! not use the number of material parameters, but the number of their
    !! components here, because a material parameter is not limited to consist
    !! of only one component. The ordering of the material components is also
    !! not fixed, but depends on the order of the material parameter definition
    !! in the configuration file.
    real(kind=rk), allocatable :: materialDat(:,:,:)
  end type atl_elemMaterialData_type

  !> Datatype provides information each face
  type atl_faceMaterialData_type
    !> The number of nodal points with material information per spatial
    !! direction.
    integer :: nPointsPerDir
    !> The nodal material properties on the face for the left element. \n
    !! First dimension is the number of faces with this material info. \n
    !! Second dimension is the number of nPointsPerDir^3 \n
    !! Third dimension is the number of material components. Attention: We do
    !! not use the number of material parameters, but the number of their
    !! components, because a material parameter is not limited to consist
    !! of only one component. \n
    real(kind=rk), allocatable :: leftElemMaterialDat(:,:,:)
    !> The nodal material properties on the face for the right element. \n
    !! First dimension is the number of faces with this material info. \n
    !! Second dimension is the number of nPointsPerDir^3 \n
    !! Third dimension is the number of material components. Attention: We do
    !! not use the number of material parameters, but the number of their
    !! components, because a material parameter is not limited to consist
    !! of only one component. \n
    real(kind=rk), allocatable :: rightElemMaterialDat(:,:,:)
  end type atl_faceMaterialData_type

  !> Datatype to provide material parameter information
  type atl_materialData_type
    !> Material parameter information for each fluid element. \n
    !! First entry covers all the constant material parameters, second
    !! entry covers all non-constant material parmaters.
    !! \sa atl_ConstMatIdx
    !! \sa atl_VarMatIdx
    type(atl_elemMaterialData_type) :: elemMaterialData(2)
    !> Material parameter information for the compute faces.\n
    !! First dimension is the number of face directions, i.e. x y and z.
    !! Second dimension covers the material combinations, i.e.:
    !! First entry covers all the constant material parameters, second
    !! entry covers all non-constant material parmaters.
    type(atl_faceMaterialData_type) :: faceMaterialData(3,2)
    !> Indicates whether this element fulfills all prerequisites to be computed
    !! with a lower resolution.
    logical, allocatable :: mode_reducable(:)
  end type atl_materialData_type

  ! -----------------------------------------------------------------------------
  !> This type is used to set the material up. It contains information from
  !! each equation systems that are used to process the information from the
  !! configuration file. Once the materials are added to the variable system,
  !! the information contained in here are not needed anymore, thus are not
  !! stored in another, longer persisting type.
  type atl_init_material_type
    !> The possible material property.
    type(tem_possible_variable_type) :: poss_materialVars

    !> This dictionary is used to map variables from the variable system to
    !! source variables to be used as data source.
    type(grw_stringKeyValuePairArray_type) :: materialDict
  end type atl_init_material_type

  ! -----------------------------------------------------------------------------
  !> Levelwise description of the material parameters in the mesh.
  type atl_material_type
    !> The maximum information propagation speed of all materials (used
    !! for calculation of the timesteps), i.e. the maximum speed of light)
    real(kind=rk) :: maxPropSpeed
    !> Description of all materials in the mesh.
    type(atl_material_property_type) :: material_desc
    !> Material parameter data for all relevant elements and faces in
    !! the mesh.
    type(atl_materialData_type) :: material_dat
    !> Postition of individual projection method in the projection list
    integer :: poly_proj_pos
    !> Postition of individual projection method in the projection list
    integer :: poly_proj_pos_state2Mat
  end type atl_material_type

  public ::                             &
    & atl_material_type,                &
    & atl_materialData_type,            &
    !!& atl_face_material_type,           &
    & atl_faceMaterialData_type,        &
    & atl_material_property_type,       &
    & atl_init_material_type,           &
    & atl_spacetime_fun_pointer_type,   &
    & atl_mixedMat_prp,                 &
    & atl_pureConstMat_prp,             &
    & atl_ConstMatIdx, atl_VarMatIdx

end module atl_materialPrp_module
