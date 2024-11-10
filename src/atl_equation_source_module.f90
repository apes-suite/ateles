! Copyright (c) 2015-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015, 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!--------------------------------------------
!    A O S - Array of structures layout new
!-------------------------------------------
! Access to get_point value output
! Access to get_element value output

!!This module contains the routine for evaluation of the source term for
!!common for all the equation system. Here the common stuff is done and
!!equation specific stuff is called as a function pointer
!!This is designed to work for all kinds of source terms for all equations
!!and for all dimentions
module atl_equation_source_module

  use env_module,                   only: rk, long_k
  use tem_varSys_module,            only: tem_varSys_type
  use tem_time_module,              only: tem_time_type

  use ply_oversample_module,        only: ply_convert2oversample,   &
    &                                     ply_convertFromoversample
  use ply_poly_project_module,      only: ply_poly_project_type, &
    &                                     assignment(=),         &
    &                                     ply_poly_project_m2n,  &
    &                                     ply_poly_project_n2m,  &
    &                                     ply_prj_body_type

  use atl_source_types_module,      only: atl_source_op_type

  implicit none

  private

  public :: atl_equation_evaluate_source_modal
  public :: atl_equation_evaluate_source_nodal
  public :: atl_compute_source_interface

  abstract interface

    subroutine atl_compute_source_interface( rhs, source, state, constants )

      ! -----------------------------------------------------------------------!
      use env_module,                   only: rk
      ! -----------------------------------------------------------------------!
      !> The Right Hand side to be updated
      real(kind=rk), intent(inout) :: rhs(:,:)
      !> The source data to be used
      real(kind=rk), intent(in) :: source(:,:)
      !> The state in the modal form
      real(kind=rk), intent(in) :: state(:,:)
      !> Some constants that might be needed for the source
      ! term evaluation
      real(kind=rk), intent(in) :: constants(:)

    end subroutine atl_compute_source_interface

  end interface


contains

! *******************************************************************************
  ! This subroutine evaluates the source term and the evaluation of RHS involves
  ! non-linear operation. i.e nodal values of state and source is to be used
  ! to calculate the RHS
  ! This framework can be used for all the equation system and work for
  ! equations in all Dimensions
  subroutine atl_equation_evaluate_source_nodal( fun, varsys, currentLevel, &
    & nDim, time, eval_rhs, state, poly_proj, polyProjBody, sourceData,     &
    & consts                                                                )
    ! -------------------------------------------------------------------------!

    !> Description of method to update source
    class(atl_source_op_type), intent(in) :: fun

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> The current Level information
    integer, intent(in) :: currentLevel

    !> The dimension information
    integer, intent(in) :: nDim

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> The pointer to the equation specific RHS evaluation
    procedure(atl_compute_source_interface), pointer :: eval_rhs

    ! State information
    real(kind=rk), intent(in) :: state(:,:,:)

    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> the data needed for the projection method
    type(ply_prj_body_type) :: polyProjBody

    ! SourceData information to be modified after the evaluation
    real(kind=rk), intent(inout) :: sourceData(:,:,:)

    ! Array of constants that may be required for the evaluation of source terms
    real(kind=rk), intent(in), optional :: consts(:)
    ! --------------------------------------------------------------------------!
    integer :: iElem, nComps, nDofs, nElems, nQuadPoints, nOverSampDofs
    integer :: posInTotal, idx_start, idx_end
    integer :: iPoint, iComp
    real(kind=rk), allocatable :: constants(:)
    real(kind=rk), allocatable :: source(:), sourceNodal(:,:)
    real(kind=rk), allocatable :: rhsModal(:,:), rhs(:,:), stateModal(:,:)
    real(kind=rk), allocatable :: rhsNodal(:,:), stateNodal(:,:)
    ! --------------------------------------------------------------------------!
    if (present(consts)) then
     allocate(constants(size(consts)))
     constants = consts
    else
     allocate(constants(0))
    endif

    nComps = varSys%method%val(fun%srcTerm_varPos)%nComponents
    nDofs = polyProjBody%nDofs
    nQuadPoints = polyProjBody%nQuadPoints
    nOverSampDofs = polyProjBody%OverSamp_dofs
    nElems = fun%elems(currentLevel)%nElems

    allocate( source(nQuadPoints * nComps ))
    allocate( sourceNodal(nQuadPoints, nComps ))
    allocate( rhsNodal( nQuadPoints, varSys%nScalars ) )
    allocate( rhsModal( nQuadPoints, varSys%nScalars ) )
    allocate( rhs( nDofs, varSys%nScalars ) )
    allocate( stateModal( noversampdofs, varSys%nScalars ) )
    allocate( stateNodal( nQuadPoints, varSys%nScalars ) )


    ! Loop over elements
    do iElem = 1, nElems

      ! 3. Call get_valOfIndex to get the modal values of data variable
      ! defined in config file
      idx_start = (iElem-1)*nQuadPoints+1
      idx_end = iElem*nQuadPoints
      call varSys%method%val(fun%data_varPos)%get_valOfIndex(            &
        & varSys  = varSys,                                              &
        & time    = time,                                                &
        & iLevel  = currentLevel,                                        &
        & idx     = fun%elems(currentLevel)%idx%val(idx_start: idx_end), &
        & nVals   = nQuadPoints,                                         &
        & res     = source                                               )

      ! Transfer the serialized result from get_element into usable array for
      ! the compute_source interface
      do iComp = 1, nComps
        do iPoint = 1, nQuadpoints
          sourceNodal(iPoint, iComp) =                             &
            & source(( 1-1)* ncomps* nquadpoints+( ipoint-1)* ncomps+icomp)
        end do
      end do

      posInTotal = fun%elems(currentLevel)%posInTotal%val(iElem)
      !  Convert the state to  nodal values
      ! --> modal space
      call ply_convert2oversample( state       = state(posInTotal,:,:), &
        &                          poly_proj   = poly_proj,             &
        &                          nDim        = nDim,                  &
        &                          modalCoeffs = stateModal(:,:)        )
      ! --> oversampled modal space

      call ply_poly_project_m2n( me         = poly_proj,       &
        &                        dim        = nDim,            &
        &                        nVars      = varSys%nScalars, &
        &                        nodal_data = stateNodal(:,:), &
        &                        modal_data = stateModal(:,:)  )
      ! --> oversamp nodal space

      ! Call compute_source_pointer_type to compute the rhs (equation specific)
      call eval_rhs(rhsNodal, sourceNodal, stateNodal, constants)

      ! Convert everything back to modal space
      call ply_poly_project_n2m( me         = poly_proj,       &
        &                        dim        = nDim,            &
        &                        nVars      = varSys%nScalars, &
        &                        nodal_data = rhsNodal(:,:),   &
        &                        modal_data = rhsModal(:,:)    )

      ! --> oversamp modal space
      call ply_convertFromoversample( modalCoeffs = rhsModal(:,:), &
        &                             poly_proj   = poly_proj,     &
        &                             nDim        = nDim,          &
        &                             state       = rhs(:,:)       )
      ! --> oversampled modal space

      ! Add rhsmodal to the source representation
      sourcedata(iElem,:,:) = sourcedata(iElem,:,:) + rhs(:,:)

    end do

    deallocate( source )
    deallocate( sourceNodal )
    deallocate( rhsNodal )
    deallocate( rhsModal )
    deallocate( rhs )
    deallocate( stateModal )
    deallocate( stateNodal )

  end subroutine atl_equation_evaluate_source_nodal
! *******************************************************************************


! *******************************************************************************
  ! This subroutine evaluates the source term and the evaluation of RHS involves
  ! linear operation. i.e the modal values of state and source is to be used
  ! to calculate the RHS
  subroutine atl_equation_evaluate_source_modal( fun, varsys, currentLevel, &
    & nDim, time, eval_rhs, state, poly_proj, polyProjBody, sourceData,     &
    & consts                                                                )
    ! -------------------------------------------------------------------------!

    !> Description of method to update source
    class(atl_source_op_type), intent(in) :: fun

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> The current Level information
    integer, intent(in) :: currentLevel

    !> The dimension information
    integer, intent(in) :: nDim

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> The pointer to the equation specific RHS evaluation
    procedure(atl_compute_source_interface), pointer :: eval_rhs

    ! State information
    real(kind=rk), intent(in) :: state(:,:,:)

    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> the data needed for the projection method
    type(ply_prj_body_type) :: polyProjBody

    ! SourceData information to be modified after the evaluation
    real(kind=rk), intent(inout) :: sourceData(:,:,:)

    ! Array of constants that may be required for the evaluation of source terms
    real(kind=rk), intent(in), optional :: consts(:)
    ! --------------------------------------------------------------------------!
    integer :: iElem, nComps, nDofs, nElems, nQuadPoints, nOverSampDofs
    integer :: posInTotal, idx_start, idx_end
    integer :: iPoint, iComp
    real(kind=rk), allocatable :: constants(:)
    real(kind=rk), allocatable :: source(:), sourceNodal(:,:)
    real(kind=rk), allocatable :: src(:,:), sourceModal(:,:)
    real(kind=rk), allocatable :: rhs(:,:)
    ! --------------------------------------------------------------------------!
    if (present(consts)) then
     allocate(constants(size(consts)))
     constants = consts
    else
     allocate(constants(0))
    endif

    nComps = varSys%method%val(fun%srcTerm_varPos)%nComponents
    nElems = fun%elems(currentLevel)%nElems
    nDofs = polyProjBody%nDofs
    nQuadPoints = polyProjBody%nQuadPoints
    nOverSampDofs = polyProjBody%OverSamp_dofs

    allocate( source(nQuadPoints*nComps ))
    allocate( sourceNodal(nQuadPoints, nComps ))
    allocate( src( nOverSampDofs, nComps ) )
    allocate( sourceModal( nDofs, nComps ) )
    allocate( rhs( nDofs, varSys%nScalars ) )

    ! Loop over elements
    do iElem = 1, nElems

      ! Call get_valOfIndex to get the modal values of data variable
      ! defined in config file
      idx_start = (iElem-1)*nQuadPoints+1
      idx_end = iElem*nQuadPoints
      call varSys%method%val(fun%data_varPos)%get_valOfIndex(            &
        & varSys  = varSys,                                              &
        & time    = time,                                                &
        & iLevel  = currentLevel,                                        &
        & idx     = fun%elems(currentLevel)%idx%val(idx_start: idx_end), &
        & nVals   = nQuadPoints,                                         &
        & res     = source                                               )

      ! Transfer the serialized result from get_element into a 3d array to have
      ! it fulfill the interface requirements from ply_convert2oversample
      do iComp = 1, nComps
        do iPoint = 1, nQuadpoints
          sourceNodal(iPoint, iComp) =                             &
            & source(( 1-1)* ncomps* nquadpoints+( ipoint-1)* ncomps+icomp)
        end do
      end do

      posInTotal = fun%elems(currentLevel)%posInTotal%val(iElem)

      ! Convert the source back to modal space
      call ply_poly_project_n2m( me         = poly_proj,        &
        &                        dim        = nDim,             &
        &                        nVars      = nComps,           &
        &                        nodal_data = sourceNodal(:,:), &
        &                        modal_data = src(:,:)          )

      ! --> oversamp modal space
      call ply_convertFromoversample( modalCoeffs = src(:,:),        &
        &                             poly_proj   = poly_proj,       &
        &                             nDim        = nDim,            &
        &                             state       = SourceModal(:,:) )
      ! --> oversampled modal space

      ! Call compute_source_pointer_type to compute the rhs (equation specific)
      call eval_rhs(rhs, sourceModal , state(posInTotal,:,:), constants)

      ! Add rhsmodal to the source representation
      sourcedata(iElem,:,:) = sourcedata(iElem,:,:) + rhs(:,:)

    end do


  end subroutine atl_equation_evaluate_source_modal
! *******************************************************************************

end module atl_equation_source_module
