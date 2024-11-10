! Copyright (c) 2014-2015 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015-2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
!! Module for datatypes and routines related to
!! penalization procedures.
module atl_penalization_module
  use, intrinsic :: iso_c_binding, only: c_f_pointer, &
    &                                    c_loc
  use env_module,                  only: rk

  ! Aotus modules
!PV!  use aotus_module,                only: flu_State

  ! Treelm modules
  use tem_element_module,          only: eT_fluid
  use treelmesh_module,            only: treelmesh_type
  use tem_element_module,          only: eT_fluid
  use tem_dyn_array_module,        only: positionOfVal
  use tem_spacetime_fun_module,    only: tem_st_fun_listElem_type
  use tem_dyn_array_module,        only: PositionOfVal

  ! ATELES modules
  use atl_scheme_module,           only: atl_scheme_type
  use atl_equation_module,         only: atl_equations_type
  use atl_cube_elem_module,        only: atl_cube_elem_type
!PV!  use atl_materialPrp_module,      only: atl_material_type

  implicit none
  private


  !> Datatype collects all volumetric penalization information.
  type atl_penalizationData_type

    !> Logical indicating if this datatype is filled with
    !! valueable data.
    logical :: isActive = .false.

    real(kind=rk), allocatable :: penalization_data(:,:,:)

  end type atl_penalizationData_type


  public :: atl_penalizationData_type, atl_init_penalization

contains

  !> Routine to init container for penalization data.
  subroutine atl_init_penalization( tree, penalizationdata_list, scheme_list, &
    &                               equation, mesh_list )
!PV!    , penalization_list,   &
!PV!    &                               conf, levelPointer                        )
    ! ---------------------------------------------------------------------------
    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(in) :: tree
    !> The penalization data to be initialized.
    type(atl_penalizationData_type), intent(inout) :: &
      & penalizationdata_list( tree%global%minlevel:tree%global%maxlevel )
    !> The list of schemes. One for each level.
    type(atl_scheme_type), intent(in) :: &
      & scheme_list(tree%global%minlevel:tree%global%maxlevel)
    !> The equation to be solved.
    type(atl_equations_type), intent(inout) :: equation
    !> The mesh description for each level
    type(atl_cube_elem_type), intent(in) :: &
      & mesh_list(tree%global%minlevel:tree%global%maxlevel)
!PV!    !> The description of the material properties on the element basis.
!PV!    type(atl_material_type), intent(inout) :: &
!PV!      & penalization_list(tree%global%minlevel:tree%global%maxlevel)
!PV!    !> Flu binding to configuration file.
!PV!    type(flu_State), intent(inout) :: conf
!PV!    !> Pointer from treeIDlist entry to level-wise fluid part of total list.
!PV!    !! The length of this vector is the total number of cubic elements.
!PV!    integer, intent(in)  :: levelPointer(:)
    ! ---------------------------------------------------------------------------
    integer :: iLevel, matVarPos
    type(tem_st_fun_listElem_type), pointer :: stFunList
    ! ---------------------------------------------------------------------------

    ! Check if the equation can be penalized ...
    select case(equation%eq_kind)
    case('euler', 'navier_stokes', 'filtered_navier_stokes', 'euler_2d', &
      &   'navier_stokes_2d', 'euler_1d', 'filtered_navier_stokes_2d'    )

      ! Check if a variable penalization is present. The presence of this
      ! variable activates the penalization.
      matVarPos = positionOfVal( equation%varSys%varname, 'characteristic' )

      if( matVarPos > 0 ) then
        call c_f_pointer(equation%varSys%method%val(matVarPos)%method_data, stFunList)
        if( any(stFunList%val(:)%fun_kind /= 'const') ) then
          do iLevel = tree%global%minlevel, tree%global%maxlevel
            penalizationdata_list(iLevel)%isActive = .true.
            ! Allocate the penalization arraz with the same number of components
            ! the state has to later on use array operations to apply penalization
            ! to the state.
            allocate(penalizationdata_list(iLevel)%penalization_data( &
              & mesh_list(iLevel)%descriptor%elem%nElems(eT_fluid),   &
              & scheme_list(iLevel)%nDofs,                            &
              & equation%varSys%nScalars                            ) )
            penalizationdata_list(iLevel)%penalization_data = 0.0_rk
          end do
        endif
      else
        penalizationdata_list(:)%isActive = .false.
      end if

    !NA! Is it implemented for maxwell ?
    !case('maxwell', 'maxwell_2d')
    case( 'maxwell_2d')

      ! Check if a variable conductivity is present. The presence of this
      ! variable activates the penalization.
      if( positionOfVal( equation%varSys%varname, 'conductivity' ) > 0 ) then
        do iLevel = tree%global%minlevel, tree%global%maxlevel
          penalizationdata_list(iLevel)%isActive = .true.
          ! Allocate the penalization arraz with the same number of components
          ! the state has to later on use array operations to apply penalization
          ! to the state.
          allocate(penalizationdata_list(iLevel)%penalization_data( &
            & mesh_list(iLevel)%descriptor%elem%nElems(eT_fluid),   &
            & scheme_list(iLevel)%nDofs,                            &
            & equation%varSys%nScalars                            ) )
          penalizationdata_list(iLevel)%penalization_data = 0.0_rk
        end do
      else
        penalizationdata_list(:)%isActive = .false.
      end if

    case default
      penalizationdata_list(:)%isActive = .false.
    end select

  end subroutine atl_init_penalization


end module atl_penalization_module

