! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

module atl_physFlux_module

  use env_module,                     only: rk
  use atl_penalization_module,        only: atl_penalizationData_type
  use ply_poly_project_module,        only: ply_poly_project_type
  use atl_materialPrp_module,         only: atl_material_type
  use atl_scheme_module,              only: atl_scheme_type
  use atl_equation_module,            only: atl_equations_type


  implicit none

  private

  public :: atl_physflux_pointer_type
  public :: atl_penalization_pointer_type
  public :: physFlux_interface

  type atl_physflux_pointer_type
    procedure(physFlux_interface), pointer, nopass :: physFlux => null()
  end type atl_physFlux_pointer_type

  type atl_penalization_pointer_type
    procedure(penalization_interface), pointer, nopass :: Pen => null()
  end type atl_penalization_pointer_type


  !> This is the interface of the physical flux pointers. State in nodal or modal
  ! form can be passed and the result can be pbtained in both nodal or modal form
  abstract interface
    subroutine physFlux_interface(  equation, res, state, iElem, iDir,        &
                                  &  penalizationData, poly_proj, material,   &
                                  &  nodal_data, nodal_GradData, nodal_res,   &
                                  &  elemLength, scheme_min, scheme_current   )
      ! --------------------------------------------------------------------------

      import :: atl_equations_type, atl_penalizationData_type,                &
                ply_poly_project_type, atl_material_type, atl_scheme_type, rk
      ! -------------------------------------------------------------------------
      !> The equation system we are working with
      type(atl_equations_type), intent(in) :: equation
      !> The result in the modal form
      real(kind=rk), intent(inout)     :: res(:,:)
      !> The state in the modal form
      real(kind=rk), intent(in), optional :: state(:,:)
      !> The current element index
      integer, intent(in) :: iElem
      !> The current direction
      integer, intent(in) :: iDir
      !> The Penalization data
      type(atl_penalizationData_type), intent(inout) :: penalizationData
      !> The projection datatype for the projection information
      type(ply_poly_project_type), intent(inout) :: poly_proj
      !> The material information
      type(atl_material_type), intent(inout) :: material
      !> The state data in the nodal form
      real(kind=rk), intent(in), optional :: nodal_data(:,:)
      !> The state data in the nodal form
      real(kind=rk), intent(in), optional :: nodal_GradData(:,:,:)
      !> The result in the nodal form
      real(kind=rk), intent(inout)     :: nodal_res(:,:)
      !> The length of the current element
      real(kind=rk), intent(in) :: ElemLength
      !> The scheme information of the min level (This is needed for the temp
      ! buffer array for evaluating the physical fluxes )
      type(atl_scheme_type), intent(inout) :: scheme_min
      !> Information about the current level
      type(atl_scheme_type), intent(inout) :: scheme_current
    end subroutine physFlux_interface
  end interface

  !> To calculate the penalization term for density, momentum and energy
  abstract interface
    subroutine penalization_interface(  equation,  poly_proj, nodal_data,      &
                         &  scheme_min, penalizationData, iElem, material  )
      ! --------------------------------------------------------------------------

      import :: atl_equations_type, atl_penalizationData_type,                &
                ply_poly_project_type, atl_material_type, atl_scheme_type, rk
      ! -------------------------------------------------------------------------
      !> The equation system we are working with
      type(atl_equations_type), intent(in) :: equation
      !> The current element index
      integer, intent(in) :: iElem
      !> The Penalization data
      type(atl_penalizationData_type), intent(inout) :: penalizationData
      !> The projection datatype for the projection information
      type(ply_poly_project_type), intent(inout) :: poly_proj
      !> The material information
      type(atl_material_type), intent(inout) :: material
      !> The state data in the nodal form
      real(kind=rk), intent(in), optional :: nodal_data(:,:)
      !> The scheme information of the min level (This is needed for the temp
      ! buffer array for evaluating the physical fluxes )
      type(atl_scheme_type), intent(inout) :: scheme_min
    end subroutine penalization_interface
  end interface


end module atl_physFlux_module
