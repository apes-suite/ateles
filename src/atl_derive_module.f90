! Copyright (c) 2014, 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
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

module atl_derive_module
  use env_module,                   only: rk

  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use tem_varSys_module,            only: tem_varSys_type,   &
    &                                     tem_varSys_op_type
  use tem_topology_module,          only: tem_levelOf

  use atl_varSys_module,            only: atl_varSys_solverData_type
  use ply_oversample_module,        only: ply_convert2oversample,   &
    &                                     ply_convertFromOversample
  use ply_poly_project_module,      only: ply_poly_project_m2n, &
     &                                    ply_poly_project_n2m, &
     &                                    ply_poly_project_type
  use ply_dof_module,               only: ply_dof_2degree

  implicit none

  private

  public :: atl_derive_inputVar_type,         &
    &       atl_derive_fromModalData,         &
    &       atl_generic_fromModal_getElement

  !> This type stores the state data for a given variable and a given element.
  !! This type is used to allocate arrays for arbitrary variable counts, as
  !! jagged arrays are not allowed in Fortran.
  type atl_derive_inputVar_type
    real(kind=rk), allocatable :: data(:,:)
  end type

  abstract interface
    !> This interface describes the arguments to be used for routines that do
    !! the derivation of variables from the variable system. These derive
    !! routines are passed to atl_derive_getElement via a procedure pointer.
    subroutine atl_derive_fromModalData(fun, varsys, tree, iElem, elemPos,   &
      &                                               nodalInput, nodalRes   )
      import :: tem_varSys_op_type, &
        &       tem_varSys_type, &
        &       rk, &
        &       atl_derive_inputVar_type, &
        &       treelmesh_type

      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in) :: varSys

      !> global treelm mesh info
      type(treelmesh_type), intent(in) :: tree

      !> The Current element index
      integer, intent(in) :: iElem

      !> TreeID of the element to get the variable for.
      integer, intent(in) :: elempos(:)

      !> The input data. nodalInput contains one entry for each input variable.
      !! This entry itself contains the nodal data for the dofs and components of
      !! the input variable. These nodal data has to be gained by oversampling
      !! and projecting the modal state into nodal space.
      type(atl_derive_inputVar_type) :: nodalInput(:)
      !> The result in nodal space
      real(kind=rk), allocatable :: nodalRes(:,:)
    end subroutine atl_derive_fromModalData
  end interface


contains


  !> This routine prepares the data for variable derivation or operators. It
  !! gathers all input variables from the variable system, oversamples and
  !! projects them into nodal space, calls the function with the actual
  !! calculation and projects the results back into modal space.
  !! As these projections are common to all elementwise variable accesses, this
  !! generic routine does all necessary operations in a generic way.
  recursive subroutine atl_generic_fromModal_getElement( fun, varsys, elempos, &
      &                                                  time,tree, nElems,    &
      &                                                  nDofs, fnCalcPtr,     &
      &                                                  solverData, res )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    !! Careful: This argument determines how many dofs is returned
    !! in the output of this routine
    integer, intent(in) :: nDofs

    !> Function pointer to perform specific operation.
    !! Sets in routine which calls this routine
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr

    !> Solver data container
    type(atl_varSys_solverData_type), intent(inout) :: solverData

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_derive_inputVar_type), allocatable :: inputState(:,:)
    type(atl_derive_inputVar_type), allocatable :: nodalInput(:)
    real(kind=rk), allocatable :: input(:)
    real(kind=rk), allocatable :: modalInput(:,:)
    real(kind=rk), allocatable :: modalRes(:,:)
    real(kind=rk), allocatable :: nodalRes(:,:), res2d(:,:)
    integer, allocatable :: nComps(:)
    integer :: iElem, iDof, iComp, i, level, pos
    ! This variable would have the actual number of degrees of freedom
    ! of the simulation
    integer :: numDofs
    integer, allocatable :: levelDofs(:)
    integer, allocatable :: overSamp_dofs(:)
    type(ply_poly_project_type) :: projection
    ! --------------------------------------------------------------------------

    allocate(inputState(fun%nInputs,nElems))
    allocate(nodalInput(fun%nInputs))
    allocate(nComps(fun%nInputs))
    allocate(leveldofs(size(solverData%polyProjectPtr)))
    allocate(overSamp_dofs(size(solverData%polyProjectPtr)))
    select case(solverData%equationPtr%nDimensions)
    case(1)
      overSamp_dofs = solverData%polyProjectPtr(:)%body_1d%oversamp_dofs
      leveldofs = solverData%polyProjectPtr(:)%body_1d%nDofs
    case(2)
      overSamp_dofs = solverData%polyProjectPtr(:)%body_2d%oversamp_dofs
      leveldofs = solverData%polyProjectPtr(:)%body_2d%nDofs
    case(3)
      overSamp_dofs = solverData%polyProjectPtr(:)%body_3d%oversamp_dofs
      leveldofs = solverData%polyProjectPtr(:)%body_3d%nDofs
    end select

    allocate(res2d(ndofs,fun%nComponents))

    do i = 1, fun%nInputs
      nComps(i) = varSys%method%val(fun%input_varPos(i))%nComponents
      do iElem = 1, nElems
        level = tem_levelOf(tree%treeID(elemPos(iElem)))
        pos = solverData%poly_proj_posPtr(level)
        numdofs = leveldofs(pos)
        allocate(inputState(i,iElem)%data(numdofs,nComps(i)))

        allocate(input(nComps(i)*numDofs))

        call varSys%method%val(fun%input_varPos(i))%get_element( &
          & varSys  = varSys,                                    &
          & elemPos = elemPos(iElem:iElem),                      &
          & time    = time,                                      &
          & tree    = tree,                                      &
          & nElems  = 1,                                         &
          & nDofs   = numDofs,                                   &
          & res     = input                                      )

        do iComp = 1, nComps(i)
          do iDof = 1, numDofs
            inputState(i,iElem)%data(iDof, iComp) =                &
              & input((( 1-1)* ncomps(i)* numdofs+( idof-1)* ncomps(i)+icomp))
          end do
        end do
        deallocate(input)
      end do
    end do

    do iElem = 1, nElems
      level = tem_levelOf(tree%treeID(elemPos(iElem)))
      pos = solverData%poly_proj_posPtr(level)
      allocate(modalRes(overSamp_dofs(pos), fun%nComponents))
      allocate(nodalRes(overSamp_dofs(pos), fun%nComponents))
      do i = 1, fun%nInputs
        allocate(nodalInput(i)%data(overSamp_dofs(pos), nComps(i)))
        allocate(modalInput(overSamp_dofs(pos), nComps(i)))
        call ply_convert2oversample(                          &
          & state       = inputState(i,iElem)%data,           &
          & poly_proj   = solverData%polyProjectPtr(pos),     &
          & nDim        = solverData%equationPtr%nDimensions, &
          & modalCoeffs = modalInput                          )
        deallocate(inputState(i,iElem)%data)

        nodalInput(i)%data = 0.0_rk
        call ply_poly_project_m2n(                                           &
          & me         = solverData%polyProjectPtr(pos),                     &
          & dim        = solverData%equationPtr%nDimensions,                 &
          & nVars      = varSys%method%val(fun%input_varPos(i))%nComponents, &
          & nodal_data = nodalInput(i)%data,                                 &
          & modal_data = modalInput                                          )
        deallocate(modalInput)
      end do

      !! Call the procedure that does the calculation
      call fnCalcPtr(fun, varSys, tree, iElem, elemPos, nodalInput, nodalRes)

      modalRes = 0.0_rk
      call ply_poly_project_n2m(                            &
        & me         = solverData%polyProjectPtr(pos),      &
        & dim        = solverData%equationPtr%nDimensions,  &
        & nVars      = fun%nComponents,                     &
        & nodal_data = nodalRes,                            &
        & modal_data = modalRes                             )

      res2d = 0.0_rk
      projection = solverData%polyProjectPtr(pos)
      projection%min_degree = ply_dof_2degree(ndofs = ndofs,                 &
        &                                     space = projection%basisType,  &
        &                                     ndims = solverData%equationPtr &
        &                                                       %nDimensions )

      projection%min_degree = nint( ndofs**( 1._rk                     &
        &                                    / solverData%equationPtr  &
        &                                                %nDimensions) ) - 1
      call ply_convertFromOversample(                       &
        & poly_proj   = projection,                         &
        & modalCoeffs = modalRes,                           &
        & nDim        = solverData%equationPtr%nDimensions, &
        & state       = res2d                               )

      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res(( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp) = &
            & res2d(iDof, iComp)
        end do
      end do

      do i = 1, fun%nInputs
        deallocate(nodalInput(i)%data)
      end do

      deallocate(modalRes)
      deallocate(nodalRes)
    end do ! Elems


    deallocate(nodalInput)
    deallocate(nComps)
    deallocate(inputState)
    deallocate(leveldofs)
    deallocate(overSamp_dofs)
    deallocate(res2d)
  end subroutine atl_generic_fromModal_getElement
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module atl_derive_module
