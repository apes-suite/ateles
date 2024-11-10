! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014, 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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
!! Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme for the advection equation. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
module atl_modg_1d_advection_kernel_module
  use env_module,               only: rk
  use treelmesh_module,         only: treelmesh_type
  use tem_faceData_module,      only: tem_faceIterator_type
  use treelmesh_module,         only: treelmesh_type

  use atl_equation_module,      only: atl_equations_type
  use atl_cube_elem_module,     only: atl_cube_elem_type
  use atl_facedata_module,      only: atl_facedata_type

  implicit none

  private

  public :: atl_modg_1d_advection_numflux, atl_modg_1d_advection_physFlux

contains

  !> Calculate the physical flux for the MODG scheme and
  !! advection equation.
  subroutine atl_modg_1d_advection_physFlux( equation,  state, res )
    !--------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type)          :: equation
    real(kind=rk), intent(in) :: state(:,:)
    real(kind=rk), intent(inout) :: res(:,:)
    !--------------------------------------------------------------------------!


    res(:,1) = state(:,1) * equation%advection%velocity



  end subroutine atl_modg_1d_advection_physFlux

  !> Calculate the numerical flux for the advection equation and MODG scheme
  subroutine atl_modg_1d_advection_numFlux( mesh, equation, facedata )
    !--------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !--------------------------------------------------------------------------
    integer :: iDir
    !--------------------------------------------------------------------------

    ! Numerical flux for faces in all 1 spatial face directions (x dir)
    do iDir = 1,1
      call modg_1d_advection_oneDim_numFlux( &
        &    equation = equation , facedata = facedata, &
        &    faces = mesh%faces%faces(iDir)%computeFace, faceDir = iDir )
    end do

  end subroutine atl_modg_1d_advection_numFlux

  !> Numerical flux calculation for the advection equation across the faces in a single
  !! spatial direction.
  subroutine modg_1d_advection_oneDim_numFlux( equation, facedata, faces, &
    &                                          faceDir )
    !--------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face state if the equation
    type(atl_facedata_type), intent(inout) :: facedata
    !> The faces to calculate the fluxes for.
    type(tem_faceIterator_type), intent(in) :: faces
    !> The spatial direction of the faces you calc the fluxes for, use the following:
    !! 1 --> x direction.
    integer, intent(in) :: faceDir
    !--------------------------------------------------------------------------!
    real(kind=rk) :: flux
    ! Loop var for the faces.
    integer :: iside
    ! Element positions of the left and right element of the face.
    integer :: left_neighbor, right_neighbor
    !--------------------------------------------------------------------------

    ! Density is transported to the positive x direction, for upwinding take
    ! the left value
    if(equation%advection%velocity .gt. 0.0_rk ) then
     ! Loop over all fluid the faces in x direction
     do iside = 1, size(faces%leftPos)
       ! Get the fluid neighbors for this face.
       left_neighbor = faces%leftPos(iside)
       right_neighbor = faces%rightPos(iside)

       flux = facedata%faceRep(faceDir)%dat(left_neighbor,1,1,2) * equation%advection%velocity

       facedata%faceFlux(faceDir)%dat(left_neighbor,1,1, 2) = flux
       facedata%faceFlux(faceDir)%dat(right_neighbor,1,1, 1) = flux
      end do
    ! Density is transported to the negative x direction, for upwinding take
    ! the right value
    else
     do iside = 1, size(faces%leftPos)
       ! Get the fluid neighbors for this face.
       left_neighbor = faces%leftPos(iside)
       right_neighbor = faces%rightPos(iside)

       flux = facedata%faceRep(faceDir)%dat(right_neighbor,1,1,1) * equation%advection%velocity

       facedata%faceFlux(faceDir)%dat(left_neighbor,1,1, 2) = flux
       facedata%faceFlux(faceDir)%dat(right_neighbor,1,1, 1) = flux
      end do
    end if


  end subroutine modg_1d_advection_oneDim_numFlux


end module atl_modg_1d_advection_kernel_module
