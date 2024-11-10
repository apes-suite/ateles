! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014, 2016-2017, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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

!> summary: computational loops to calculate element-local RHS, needed for
!> local predictors (predictor-corrector approach for local timestepping)
module atl_compute_local_module

  ! Treelm modules
  use treelmesh_module,        only: treelmesh_type
  use tem_aux_module,          only: tem_abort
  use tem_general_module,      only: tem_general_type
  use tem_timer_module,        only: tem_startTimer, tem_stopTimer

  ! Ateles modules
  use atl_timer_module,        only: atl_timerHandles
  use atl_cube_elem_module,    only: atl_cube_elem_type
  use atl_kerneldata_module,   only: atl_statedata_type
  use atl_source_types_module, only: atl_source_type
  use atl_scheme_module,       only: atl_scheme_type,    &
    &                                atl_modg_scheme_prp
  use atl_equation_module,     only: atl_equations_type
  use atl_boundary_module,     only: atl_level_boundary_type
  use atl_modg_kernel_module,  only: atl_preprocess_modg_kernel

  use atl_materialPrp_module,  only: atl_material_type
  use ply_poly_project_module, only: ply_poly_project_type, assignment(=)

  implicit none

  private

  public :: atl_preprocess_local_rhs

  interface atl_preprocess_local_rhs
    module procedure preprocess_local_rhs_cubes
  end interface


contains


  ! ************************************************************************ !
  !> summary: should only evaluate sources
  subroutine preprocess_local_rhs_cubes( mesh, statedata, scheme, general, &
    &                                    boundary, poly_proj_list, equation, &
    &                                    material )
    ! -------------------------------------------------------------------- !
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type),  intent(inout) :: mesh
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme
    !> General treelm structures.
    type(tem_general_type), intent(inout) :: general
    !> Data for projection, for each level.
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> The equation you are operating with.
    type(atl_equations_type),intent(inout) :: equation
    !> The material description
    type(atl_material_type), intent(inout) :: material
    type(atl_level_boundary_type), intent(in)  :: boundary
    ! -------------------------------------------------------------------- !

    select case(scheme%scheme)
    case(atl_modg_scheme_prp) ! MODG kernel
      select case(trim(equation%eq_kind))
      case( 'maxwell', 'maxwelldivcorrection', 'euler', 'acoustic', &
        & 'lineareuler', 'loclineuler'                              )

        call tem_startTimer( timerHandle = atl_timerHandles%preprocessKernel )
        ! Update the source terms.
        call atl_preprocess_modg_kernel(                                &
          & equation           = equation,                              &
          & statedata          = statedata,                             &
          & mesh               = mesh,                                  &
          & boundary           = boundary,                              &
          & scheme             = scheme,                                &
          & material           = material,                              &
          & poly_proj_material = poly_proj_list(                        &
          &                        material%poly_proj_pos),             &
          & commPattern        = general%commPattern,                   &
          & proc               = general%proc                           )

        call tem_stopTimer( timerHandle = atl_timerHandles%preprocessKernel )

      case default
        call tem_abort( 'ERROR in preprocess_local_rhs_cubes: unknown' &
          & // ' equation for MODG'                                    )
      end select

    case default
      call tem_abort( 'ERROR in preprocess_local_rhs_cubes: unknown scheme' )
    end select

  end subroutine preprocess_local_rhs_cubes
  ! ************************************************************************ !


end module atl_compute_local_module
