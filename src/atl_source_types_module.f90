! Copyright (c) 2014-2015, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> The only purpose of this module is to keep the types contained
!! separated from the routines in atl_source_module to avoid circular
!! references between atl_source_module and atl_equation_module. In a
!! future rework of these structures, this module should be merged back
!! into the atl_source_module.
module atl_source_types_module
  use env_module,               only: rk
  use tem_varSys_module,        only: tem_varSys_type
  use tem_time_module,          only: tem_time_type
  use tem_stringkeyvaluepair_module, only: grw_stringKeyValuePairArray_type
  use tem_varMap_module,        only: tem_possible_variable_type
  use tem_grow_array_module,    only: grw_intArray_type

  use ply_poly_project_module,  only: ply_poly_project_type

  use atl_cube_elem_module,     only: atl_cube_elem_type

  implicit none

  private

  public :: atl_source_type
  public :: atl_init_source_type
  public :: atl_source_op_type
  public :: atl_source_elems_type
  public :: atl_update_source
  public :: atl_eqn_sourceMap_type

  ! -----------------------------------------------------------------------------
  !> Contains information about the elements affected by a source term. The
  !! information contained are for one particular level.
  type atl_source_elems_type
    !> The number of elements affected by the source on this level
    integer :: nElems = 0

    !> Position of elements in state array to apply source terms.
    !! Position in state array is same as position in total list
    !! Size: nElems
    type(grw_intArray_type) :: posInTotal

    !> Index to access point data type to retrieve values from variable
    !! refered for source variable
    type(grw_intArray_type) :: idx

  end type atl_source_elems_type
  ! -----------------------------------------------------------------------------

  !> Description contains list of elements on which source is active and
  !! function pointer to update source
  type atl_source_op_type
    !> Position of this source term variable in the variable system
    integer :: srcTerm_varPos

    !> Position of data variable usally st-fun provided in config file
    !! in the varSys
    integer :: data_varPos

    !> Contains levelwise information about the elements affected by the source
    type(atl_source_elems_type), allocatable :: elems(:)

    !> True for permanent sources, active on global domain so no need to
    !! initialize source elems array.
    logical :: isPermanent = .false.

    !> SourceData to be stored for a particular level
    !> buffer for evaluated source terms (contribution to the right hand side of
    !! the PDE)
    real(kind=rk), allocatable :: val(:,:,:)

    !> Function to update state with source term
    procedure(atl_update_source), pointer :: updateSrc => null()

  end type atl_source_op_type
  ! -----------------------------------------------------------------------------


  ! -----------------------------------------------------------------------------
  !> Description of the new source type in Ateles
  type atl_source_type
    !> Contains source elements position in tree%treeID and
    !! function pointer to update source.
    !! size: varDict%nVals
    type(atl_source_op_type), allocatable :: method(:)

    !> Postition of individual projection method in the projection list.
    !! This list is levelwise, size: (minLevel:maxLevel)
    integer, allocatable :: poly_proj_pos(:)

    !> This dictionary is used to map variables from the variable system to
    !! source variables to be used as data source.
    type(grw_stringKeyValuePairArray_type) :: varDict
  end type atl_source_type
  ! -----------------------------------------------------------------------------

  ! -----------------------------------------------------------------------------
  !> datatype containing mapping of source variables to function pointers. The
  !! function this pointer is pointing to is used to evaluate the space time
  !! function defining the values for the source variable.
  type atl_eqn_sourceMap_type
    procedure(atl_update_source), pointer, nopass :: do => NULL()
  end type atl_eqn_sourceMap_type
  ! -----------------------------------------------------------------------------

  ! -----------------------------------------------------------------------------
  !> This type is used to set the source terms up. It contains information from
  !! the equation systems that are used to process the information from the
  !! configuration file. Once the source terms are added to the variable system,
  !! the information contained in here are not needed anymore, thus are not
  !! stored in another, longer persisting type.
  type atl_init_source_type
    !> The possible source variables
    type(tem_possible_variable_type) :: poss_srcVars

    !> for each possible source (a variable in sources%poss_varSys) a pointer to the
    !! function that adds the necessary source term to the right hand side
    !! vector of the equation system
    type(atl_eqn_sourceMap_type), allocatable :: eval_source(:)
  end type atl_init_source_type
  ! -----------------------------------------------------------------------------


  ! -----------------------------------------------------------------------------
  !> Abstract interface to update state with source terms
  abstract interface
    subroutine atl_update_source( fun, varSys, time, mesh, poly_proj,       &
      &                           currentLevel, state, material, sourcedata )
      ! -------------------------------------------------------------------------
      import :: atl_source_op_type,    &
        &       tem_varSys_type,       &
        &       tem_time_type,         &
        &       atl_cube_elem_type,    &
        &       ply_poly_project_type, &
        &       rk
      ! -------------------------------------------------------------------------

      !> Description of method to update source
      class(atl_source_op_type), intent(in) :: fun

      !> variable system
      type(tem_varSys_type), intent(in) :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)  :: time

      !> Current level mesh information
      type(atl_cube_elem_type), intent(in) :: mesh

      !> Parameters for projection
      type(ply_poly_project_type), intent(inout) :: poly_proj

      !> current level
      integer, intent(in) :: currentLevel

      !> The state in modal space.
      !! This is needed for several source terms that have to be applied to the
      !! current state
      real(kind=rk), intent(in) :: state(:,:,:)

      !> The background material to use for evaluating specific source terms.
      !!
      !! At the moment this background material is used for all elements. At a
      !! stage each element should get it's own material parameters.
      real(kind=rk), intent(in) :: material(:)

      !> The source data to update. When all source terms are added to this
      !! buffer, it is applied to the state.
      real(kind=rk), intent(inout) :: sourcedata(:,:,:)
    end subroutine atl_update_source
  end interface
  ! -----------------------------------------------------------------------------

end module atl_source_types_module

