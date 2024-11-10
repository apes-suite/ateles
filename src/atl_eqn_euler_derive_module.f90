! Copyright (c) 2013-2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014, 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> Routines to derive quantities from the state in the Euler equation system.
module atl_eqn_euler_derive_module
  use, intrinsic :: iso_c_binding,  only: c_f_pointer
  use env_module,                   only: rk, long_k

  use tem_aux_module,               only: tem_abort
  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use tem_topology_module,          only: tem_coordOfId, tem_IDofCoord, &
    &                                     tem_levelOf
  use tem_geometry_module,          only: tem_CoordOfReal, &
    &                                     tem_PosofId
  use tem_logging_module,           only: logUnit
  use tem_varSys_module,            only: tem_varSys_type, &
    &                                     tem_varSys_op_type

  use atl_aux_module,               only: atl_bubbleSortArray
  use atl_equation_module,          only: atl_equations_type
  use atl_varSys_module,            only: atl_varSys_data_type
  use atl_derive_module,            only: atl_derive_inputVar_type,      &
    &                                     atl_derive_fromModalData,      &
    &                                     atl_generic_fromModal_getElement

  implicit none

  private

  public :: atl_speedOfSound_getPoint, atl_speedOfSound_getElement
  public :: atl_pressure_getPoint, atl_pressure_getElement
  public :: atl_pressure_getIndex
  public :: atl_temperature_getPoint, atl_temperature_getElement
  public :: atl_machNumber_getPoint, atl_machNumber_getElement
  public :: atl_KineticEnergy_getPoint, atl_kineticEnergy_getElement
  public :: atl_vorticity_getPoint, atl_vorticity_getElement
  public :: atl_QCriterion_getPoint, atl_qCriterion_getElement
  public :: atl_lambda2_getPoint, atl_lambda2_getElement
  public :: atl_linindicator_getPoint, atl_linindicator_getElement
  public :: atl_eqn_euler_cons2prim
  public :: atl_eqn_euler_prim2cons
  public :: atl_eqn_euler_cons2prim_grad
  public :: atl_eqn_euler_prim2cons_grad
  public :: atl_eqn_euler_cons2prim_elems
  public :: atl_eqn_euler_prim2cons_elems
  public :: atl_eqn_euler_cons2primTemp
  public :: atl_eqn_euler_primTemp2cons
  public :: atl_eqn_euler_cons2primVel
  public :: atl_eqn_euler_primVel2cons


contains


  subroutine atl_speedOfSound_getPoint( fun, varsys, point, time,tree, nPnts, &
      &                                 res                                   )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: pressure(nPnts), density(nPnts)
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = pressure                                 )

    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = density                                  )

    res = sqrt(pressure/density*fPtr%solverData%equationPtr%Euler%isen_coef)

  end subroutine atl_speedOfSound_getPoint

  subroutine atl_deriveSpeedOfSound(fun, varsys, tree, iElem, elemPos, &
    &                                        nodalInput, nodalRes      )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer, parameter :: pressure = 1, density = 2
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    nodalRes = sqrt(nodalInput(pressure)%data &
      &   / nodalInput(density)%data          &
      &   * fPtr%solverData%equationPtr%Euler%isen_coef )

  end subroutine atl_deriveSpeedOfSound

  subroutine atl_speedOfSound_getElement(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    fnCalcPtr => atl_deriveSpeedOfSound

    call atl_generic_fromModal_getElement( &
      & fun        = fun,                  &
      & varsys     = varsys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr%solverData,      &
      & res        = res                   )

  end subroutine atl_speedOfSound_getElement

  subroutine atl_pressure_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: density(nPnts), momentum(3*nPnts), energy(nPnts)
    ! --------------------------------------------------------------------------
    ! Initialize the momentum to zero
    momentum = 0.0_rk

    call C_F_POINTER( fun%method_Data, fPtr )

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = density                                  )
    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = momentum                                 )
    call varSys%method%val(fun%input_varPos(3))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = energy                                   )

    res = (fPtr%solverData%equationPtr%Euler%isen_coef - 1.0_rk)             &
      &   * ( energy                                                         &
      &      - ( 0.5_rk / density                                            &
      &        * (momentum(1::3)**2 + momentum(2::3)**2 + momentum(3::3)**2) &
      &        )                                                             &
      &     )

  end subroutine atl_pressure_getPoint


  subroutine atl_pressure_getIndex( fun, varSys, time, iLevel, &
    &                               idx, idxLen, nVals,  res   )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables,
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: n
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: density(nVals), momentum(3*nVals), energy(nVals)
    ! --------------------------------------------------------------------------
    write(logUnit(4),*) 'Get the values of indices for derived variable ',  &
      &                  trim(varSys%varname%val(fun%myPos))

    ! Initialize the momentum to zero
    density = 0.0_rk
    momentum = 0.0_rk
    energy = 0.0_rk

    call C_F_POINTER( fun%method_Data, fPtr )

    call varSys%method%val(fun%input_varPos(1))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = fPtr%opData%input_pntIndex(1)                 &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = density                                       )
    call varSys%method%val(fun%input_varPos(2))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = fPtr%opData%input_pntIndex(2)                 &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = momentum                                      )
    call varSys%method%val(fun%input_varPos(3))%get_valOfIndex( &
      & varSys  = varSys,                                       &
      & time    = time,                                         &
      & iLevel  = iLevel,                                       &
      & idx     = fPtr%opData%input_pntIndex(3)                 &
      &           %indexLvl(iLevel)%val( idx(:) ),              &
      & nVals   = nVals,                                        &
      & res     = energy                                        )
    !!call varSys%method%val(fun%inpuit_varPos(1))%get_point( &
    !!  & varSys  = varSys,                                  &
    !!  & point   = point,                                   &
    !!  & time    = time,                                    &
    !!  & tree    = tree,                                    &
    !!  & nPnts   = nPnts,                                   &
    !!  & res     = density                                  )
    !!call varSys%method%val(fun%input_varPos(2))%get_point( &
    !!  & varSys  = varSys,                                  &
    !!  & point   = point,                                   &
    !!  & time    = time,                                    &
    !!  & tree    = tree,                                    &
    !!  & nPnts   = nPnts,                                   &
    !!  & res     = momentum                                 )
    !!call varSys%method%val(fun%input_varPos(3))%get_point( &
    !!  & varSys  = varSys,                                  &
    !!  & point   = point,                                   &
    !!  & time    = time,                                    &
    !!  & tree    = tree,                                    &
    !!  & nPnts   = nPnts,                                   &
    !!  & res     = energy                                   )

    res = (fPtr%solverData%equationPtr%Euler%isen_coef - 1.0_rk)             &
      &   * ( energy                                                         &
      &      - ( 0.5_rk / density                                            &
      &        * (momentum(1::3)**2 + momentum(2::3)**2 + momentum(3::3)**2) &
      &        )                                                             &
      &     )

  end subroutine atl_pressure_getIndex

  subroutine atl_derivePressure(fun, varsys, tree, iElem, elemPos, &
    &                                        nodalInput, nodalRes  )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer, parameter :: density = 1, momentum = 2, energy = 3
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    nodalRes(:,1) = (fPtr%solverData%equationPtr%Euler%isen_coef - 1.0_rk) &
      &   * (nodalInput(energy)%data(:,1)                                  &
      &   - 0.5_rk                                                         &
      &   / nodalInput(density)%data(:,1)                                  &
      &   * sum(array=nodalInput(momentum)%data**2,dim=2)                  )

  end subroutine atl_derivePressure

  subroutine atl_pressure_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                    nDofs, res                            )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    fnCalcPtr => atl_derivePressure

    call atl_generic_fromModal_getELement( &
      & fun        = fun,                  &
      & varsys     = varsys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr%solverData,      &
      & res        = res                   )

  end subroutine atl_pressure_getElement

  subroutine atl_temperature_getPoint(fun, varsys, point, time,tree, nPnts, res)
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: R_q
    real(kind=rk) :: pressure(nPnts), density(nPnts)
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = pressure                                 )
    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = density                                  )

    R_q = 1.0_rk / fPtr%solverData%equationPtr%Euler%r

    res = R_q * pressure / density

  end subroutine atl_temperature_getPoint

  subroutine atl_deriveTemperature(fun, varsys, tree, iElem, elemPos, &
    &                              nodalInput, nodalRes               )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: R_q
    integer, parameter :: pressure = 1, density = 2
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    R_q = 1.0_rk / fPtr%solverData%equationPtr%Euler%r

    nodalRes = R_q * nodalInput(pressure)%data / nodalInput(density)%data

  end subroutine atl_deriveTemperature

  subroutine atl_temperature_getElement(fun, varsys, elempos, time, tree, &
    &                                   nElems, nDofs, res                )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    fnCalcPtr => atl_deriveTemperature

    call atl_generic_fromModal_getELement( &
      & fun        = fun,                  &
      & varsys     = varsys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr%solverData,      &
      & res        = res                   )

  end subroutine atl_temperature_getElement

  subroutine atl_machNumber_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: density(nPnts), momentum(3*nPnts), speedOfSound(nPnts)
    integer :: iPoint
    ! --------------------------------------------------------------------------
    ! Initialize the momentum to zero
    momentum = 0.0_rk

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = density                                  )
    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = momentum                                 )
    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = speedOfSound                             )

    do iPoint = 1, nPnts
      res(iPoint) = sqrt(sum(momentum(iPoint*3-2:iPoint*3)) &
        & / (density(iPoint))**2)                           &
        & / speedOfSound(iPoint)
    end do

  end subroutine atl_machNumber_getPoint

  subroutine atl_deriveMachNumber(fun, varsys, tree, iElem, elemPos,      &
    &                                              nodalInput, nodalRes   )
    ! --------------------------------------------------------------------------
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
    ! --------------------------------------------------------------------------
    integer, parameter :: density = 1, momentum = 2, speedOfSound = 3
    type(atl_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )
    select case(fPtr%solverData%equationPtr%nDimensions)
    case(1)
      nodalRes(:,1) = sqrt(                       &
        &   ( nodalInput(momentum)%data(:,1)      &
        & / nodalInput(density)%data(:,1)  )**2 ) &
        & / nodalInput(speedOfSound)%data(:,1)


    case(2)
      nodalRes(:,1) = sqrt(                       &
        &   (nodalInput(momentum)%data(:,1)       &
        & / nodalInput(density)%data(:,1))**2     &
        & + (nodalInput(momentum)%data(:,2)       &
        & / nodalInput(density)%data(:,1))**2   ) &
        & / nodalInput(speedOfSound)%data(:,1)


    case(3)
      nodalRes(:,1) = sqrt(                     &
        &   (nodalInput(momentum)%data(:,1)     &
        & / nodalInput(density)%data(:,1))**2   &
        & + (nodalInput(momentum)%data(:,2)     &
        & / nodalInput(density)%data(:,1))**2   &
        & + (nodalInput(momentum)%data(:,3)     &
        & / nodalInput(density)%data(:,1))**2 ) &
        & / nodalInput(speedOfSound)%data(:,1)

    end select
  end subroutine atl_deriveMachNumber

  subroutine atl_machNumber_getElement(fun, varsys, elempos, time, tree, &
    &                                  nElems, nDofs, res                )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    fnCalcPtr => atl_deriveMachNumber

    call atl_generic_fromModal_getELement( &
      & fun        = fun,                  &
      & varsys     = varsys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr%solverData,      &
      & res        = res                   )

  end subroutine atl_machNumber_getElement

  subroutine atl_KineticEnergy_getPoint( fun, varsys, point, time,tree, nPnts, &
      &                                  res                                   )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: density(nPnts), momentum(3*nPnts)
    integer :: iPoint
    ! --------------------------------------------------------------------------
    ! Initialize the momentum to zero
    momentum = 0.0_rk

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = density                                  )
    call varSys%method%val(fun%input_varPos(2))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = momentum                                 )

    do iPoint = 1, nPnts
      res = 0.5* sum(momentum(iPoint*3-2:iPoint*3)**2)     &
        &                        / density(iPoint)
    end do

  end subroutine atl_kineticEnergy_getPoint

  subroutine atl_deriveKineticEnergy(fun, varsys, tree, iElem, elemPos,    &
    &                                          nodalInput, nodalRes        )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> The current element index
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
    ! --------------------------------------------------------------------------
    integer, parameter :: density = 1, momentum = 2
    ! --------------------------------------------------------------------------

    nodalRes(:,1) = 0.5                         &
      & * ( nodalInput(momentum)%data(:,1)**2   &
      &   + nodalInput(momentum)%data(:,2)**2   &
      &   + nodalInput(momentum)%data(:,3)**2 ) &
      & / nodalInput(density)%data(:,1)

  end subroutine atl_deriveKineticEnergy

  subroutine atl_kineticEnergy_getElement(fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res               )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    fnCalcPtr => atl_deriveKineticEnergy

    call atl_generic_fromModal_getELement( &
      & fun        = fun,                  &
      & varsys     = varsys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr%solverData,      &
      & res        = res                   )

  end subroutine atl_kineticEnergy_getElement

  subroutine atl_vorticity_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk), allocatable :: GradVelocity(:)
    integer :: nGradVelComp
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    nGradVelComp = varSys%method%val(fun%input_varPos(1))%nComponents
    allocate(GradVelocity(nGradVelComp*nPnts))

    ! Calculate the Gradient of velocity
    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = GradVelocity                             )

    select case(fPtr%solverData%equationPtr%nDimensions)
    case(1)
      write(logunit(1),*) "Vorticity can't be evaluated for 1D. stopping..."
      call tem_abort()
    case(2)
      res(1) = GradVelocity(3) - GradVelocity (2)
    case(3)
      ! Calculate Vorticity
      res(1) = GradVelocity(8) - GradVelocity(6)
      res(2) = GradVelocity(3) - GradVelocity(7)
      res(3) = GradVelocity(4) - GradVelocity(2)
    end select

    deallocate(GradVelocity)

  end subroutine atl_vorticity_getPoint

  subroutine atl_vorticity_getElement(fun, varsys, elempos, time, tree, nElems,&
    &                                    nDofs, res                            )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk), allocatable  :: Grad_velocity(:)
    real(kind=rk), allocatable  :: ModalGradV(:,:,:), vort(:,:)
    integer :: nGradVelComp, iComp, iElem, iDof
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    nGradVelComp = varSys%method%val(fun%input_varPos(1))%nComponents
    allocate(ModalGradV(nElems,nDofs,nGradVelComp))
    allocate(Grad_velocity(nElems*nDofs*nGradVelComp))
    allocate(Vort(nDofs,fun%nComponents))

    ! Derive the modal velocity gradient
    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elempos = elempos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = Grad_velocity                              )

    do iComp = 1, nGradVelComp
      do iDof = 1, nDofs
        do iElem = 1, nElems
          modalGradV(iElem,iDof, iComp) =                               &
            & Grad_velocity((( ielem-1)* ngradvelcomp* ndofs+( idof-1)* ngradvelcomp+icomp))
        end do
      end do
    end do

    ! The components of modalGradV for 3D are organized as
    ! modalGradV(:,:,1)   - \frac{\partial v_x}{\partial x}
    ! modalGradV(:,:,2)   - \frac{\partial v_x}{\partial y}
    ! modalGradV(:,:,3)   - \frac{\partial v_x}{\partial z}

    ! modalGradV(:,:,4)   - \frac{\partial v_y}{\partial x}
    ! modalGradV(:,:,5)   - \frac{\partial v_y}{\partial y}
    ! modalGradV(:,:,6)   - \frac{\partial v_y}{\partial z}

    ! modalGradV(:,:,7)   - \frac{\partial v_z}{\partial x}
    ! modalGradV(:,:,8)   - \frac{\partial v_z}{\partial y}
    ! modalGradV(:,:,9)   - \frac{\partial v_z}{\partial z}


    ! The components of modalGradV for 2D are organized as
    ! modalGradV(:,:,1)   - \frac{\partial v_x}{\partial x}
    ! modalGradV(:,:,2)   - \frac{\partial v_x}{\partial y}
    ! modalGradV(:,:,3)   - \frac{\partial v_y}{\partial x}
    ! modalGradV(:,:,4)   - \frac{\partial v_y}{\partial y}


    do iElem = 1,nElems

    select case(fPtr%solverData%equationPtr%nDimensions)
    case(1)
      write(logunit(1),*) "Vorticity can't be evaluated for 1D. stopping..."
      call tem_abort()
    case(2)
      Vort(:,1) = ModalGradV(iElem,:,3) - ModalGradV (iElem,:,2)
    case(3)
      ! Calculate Vorticity
      Vort(:,1) = ModalGradV(iElem,:,8) - ModalGradV (iElem,:,6)
      Vort(:,2) = ModalGradV(iElem,:,3) - ModalGradV (iElem,:,7)
      Vort(:,3) = ModalGradV(iElem,:,4) - ModalGradV (iElem,:,2)
    end select


      do iDof = 1, nDofs
        do iComp = 1, fun%nComponents
          res(( ielem-1)* fun%ncomponents* ndofs+( idof-1)* fun%ncomponents+icomp) = &
            & vort(iDof, iComp)
        end do
      end do

    end do

    deallocate(ModalGradV)
    deallocate(Grad_velocity)
    deallocate(Vort)

  end subroutine atl_Vorticity_getElement


  subroutine atl_qCriterion_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: GradVelocity(9*nPnts)
    integer :: iPoint, nGradVelocityComponent
    real(kind=rk) :: eig(3)
    real(kind=rk) :: temp_q(3,3)
    ! --------------------------------------------------------------------------

    nGradVelocityComponent = 9
    ! Calculate the Gradient of velocity
    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = GradVelocity                             )

    do iPoint = 1, nPnts
      ! This routine returns the sum of squared strain tensor and the vorticity
      ! tensor. So, the output matrix temp_q
      ! temp_q(3,3) = S^2 + \Omega^2 for each ndof
      ! where S = 0.5*(\nabla v + \nabla v^{T}) - strain tensor
      ! and  \Omega = 0.5*(\nabla v - \nabla v^{T}) - vorticity tensor
      call calc_matrix_QCriterion(                                        &
        &       gradV = GradVelocity((iPoint-1)*nGradVelocityComponent+1: &
        &                             iPoint*nGradVelocityComponent ),    &
        &                           Q  = temp_q                           )

      ! Calculate the eigenValues of the symmetric matrix temp_q
      call calc_eigenValues_3by3_matrix(temp_q,eig)

      ! Evaluate the q_criterion
      Res(1) = -0.5*(eig(1)+eig(2)+eig(3))

    end do
  end subroutine atl_qCriterion_getPoint

  !> This routine evaluates the q_criterion. The input is the nodal value of
  ! the gradient of velocity. The q_criterion is evaluated from that and the
  ! nodal value is passed back to the routine atl_generic_fromModal_getElement
  ! where it is transferred back to the modal space
  subroutine atl_deriveQcriterion(fun, varsys, tree, iElem, elemPos,    &
    &                                              nodalInput, nodalRes )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> The current element index
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer, parameter :: gradV = 1
    integer :: nPoints, level, pos, iPoint
    real(kind=rk) :: eig(3)
    real(kind=rk) :: temp_q(3,3)
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    level = tem_levelOf(tree%treeID(elemPos(iElem)))
    pos = fptr%solverData%poly_proj_posPtr(level)
    npoints = fPtr%solverData%polyProjectPtr(pos)%body_3d%nQuadPoints

    do iPoint = 1, nPoints

      ! This routine returns the sum of squared strain tensor and the vorticity
      ! tensor. So, the output matrix temp_q
      ! temp_q(3,3) = S^2 + \Omega^2 for each ndof
      ! where S = 0.5*(\nabla v + \nabla v^{T}) - strain tensor
      ! and  \Omega = 0.5*(\nabla v - \nabla v^{T}) - vorticity tensor
      call calc_matrix_QCriterion( gradV = NodalInput(GradV)%data(iPoint,:), &
        &                          Q  = temp_q                               )

      ! Calculate the eigenValues of the symmetric matrix temp_q
      call calc_eigenValues_3by3_matrix(temp_q,eig)

      ! Evaluate the q_criterion
      NodalRes(iPoint,1) = -0.5*(eig(1)+eig(2)+eig(3))

    end do

  end subroutine atl_deriveqCriterion

  subroutine atl_qCriterion_getElement( fun, varsys, elempos, time, tree, &
    &                                   nElems, nDofs, res                )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    fnCalcPtr => atl_deriveQcriterion

    call atl_generic_fromModal_getElement( &
      & fun        = fun,                  &
      & varsys     = varsys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr%solverData,      &
      & res        = res                   )

  end subroutine atl_qCriterion_getElement


  subroutine atl_lambda2_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: GradVelocity(9*nPnts)
    integer :: iPoint, nGradVelocityComponent
    real(kind=rk) :: eig(3)
    real(kind=rk) :: temp_q(3,3)
    ! --------------------------------------------------------------------------

    nGradVelocityComponent = 9
    ! Calculate the Gradient of velocity
    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = GradVelocity                             )

    do iPoint = 1, nPnts
      ! This routine returns the sum of squared strain tensor and the vorticity
      ! tensor. So, the output matrix temp_q
      ! temp_q(3,3) = S^2 + \Omega^2 for each ndof
      ! where S = 0.5*(\nabla v + \nabla v^{T}) - strain tensor
      ! and  \Omega = 0.5*(\nabla v - \nabla v^{T}) - vorticity tensor
      call calc_matrix_QCriterion(                                        &
        &       gradV = GradVelocity((iPoint-1)*nGradVelocityComponent+1: &
        &                             iPoint*nGradVelocityComponent ),    &
        &                           Q  = temp_q                           )

      ! Calculate the eigenValues of the symmetric matrix temp_q
      call calc_eigenValues_3by3_matrix(temp_q,eig)

      ! Evaluate the q_criterion
      Res(1) = eig(2)

    end do

  end subroutine atl_lambda2_getPoint

  !> This routine evaluates the lambda2 criterion. The input is the nodal value
  ! of the gradient of velocity. The lambda2 is evaluated from that and the
  ! nodal value is passed back to the routine atl_generic_fromModal_getElement
  ! where it is transferred back to the modal space
  subroutine atl_deriveLambda2(fun, varsys, tree, iElem, elemPos,       &
    &                                              nodalInput, nodalRes )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> The current element index
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
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer, parameter :: gradV = 1
    integer :: nPoints, level, pos, iPoint
    real(kind=rk) :: eig(3)
    real(kind=rk) :: lam2matrix(3,3)
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    level = tem_levelOf(tree%treeID(elemPos(iElem)))
    pos = fptr%solverData%poly_proj_posPtr(level)
    npoints = fPtr%solverData%polyProjectPtr(pos)%body_3d%nQuadPoints

    do iPoint = 1, nPoints

      ! This routine returns the sum of squared strain tensor and the vorticity
      ! tensor. So, the output matrix lam2matrix
      ! lam2matrix(3,3) = S^2 + \Omega^2 for each ndof
      ! where S = 0.5*(\nabla v + \nabla v^{T}) - strain tensor
      ! and  \Omega = 0.5*(\nabla v - \nabla v^{T}) - vorticity tensor
      call calc_matrix_QCriterion(  gradV = NodalInput(GradV)%data(iPoint,:), &
        &                           Q  = lam2matrix                           )

      ! Calculate the eigenValues of the symmetric matrix lam2matrix
      call calc_eigenValues_3by3_matrix(lam2matrix,eig)

      ! Evaluate the q_criterion
      NodalRes(iPoint,1) = eig(2)

    end do

  end subroutine atl_deriveLambda2

  subroutine atl_lambda2_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                    nDofs, res                           )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    fnCalcPtr => atl_deriveLambda2

    call atl_generic_fromModal_getELement( &
      & fun        = fun,                  &
      & varsys     = varsys,               &
      & elempos    = elempos,              &
      & time       = time,                 &
      & tree       = tree,                 &
      & nElems     = nElems,               &
      & nDofs      = nDofs,                &
      & fnCalcPtr  = fnCalcPtr,            &
      & solverData = fPtr%solverData,      &
      & res        = res                   )

  end subroutine atl_lambda2_getElement


  subroutine atl_linindicator_getPoint( fun, varsys, point, time,tree, nPnts, &
      &                                 res                                   )
    ! --------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer(kind=long_k) :: treeId
    integer :: iPnt
    integer :: coord(4)
    integer :: elempos
    integer :: pos
    integer :: level
    logical :: islinear
    ! --------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    do iPnt=1,nPnts
      coord =  tem_CoordOfReal(tree, point(iPnt,:), tree%global%maxLevel)
      treeId = tem_IdOfCoord(coord)
      ! get the position of treeid or position of the parent treeid
      elemPos = abs(tem_PosofId(treeId, tree%treeID))
      level = tem_levelOf( tree%treeID( elemPos ) )
      Pos = fPtr%solverData%levelPointer(elemPos)
      islinear = fPtr%solverdata%equationPtr%euler%linear(    &
        &          mean = fPtr%solverdata                     &
        &                     %statedata_listPtr(level)       &
        &                     %state(pos, 1, :),              &
        &          deviation = fPtr%solverdata                &
        &                          %kerneldata_listPtr(level) &
        &                          %deviation(pos,:)          )
      if (islinear) then
        res(iPnt) = 0.0_rk
      else
        res(iPnt) = 1.0_rk
      end if
    end do

  end subroutine atl_linindicator_getPoint


  subroutine atl_linindicator_getElement(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
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
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    logical :: islinear
    integer :: iElem
    integer :: level
    integer :: firstdof
    integer :: pos
    ! --------------------------------------------------------------------------
    call C_F_POINTER(fun%method_data, fPtr)

    res = 0.0_rk

    firstdof = 1
    do iElem=1,nElems
      level = tem_levelOf(tree%treeID(elempos(iElem)))
      pos = fptr%solverData%levelPointer(elemPos(iElem))

      islinear = fPtr%solverdata%equationPtr%euler%linear(    &
        &          mean = fPtr%solverdata                     &
        &                     %statedata_listPtr(level)       &
        &                     %state(pos, 1, :),              &
        &          deviation = fPtr%solverdata                &
        &                          %kerneldata_listPtr(level) &
        &                          %deviation(pos,:)          )
      if (islinear) then
        res(firstdof) = 0.0_rk
      else
        res(firstdof) = 1.0_rk
      end if

      firstdof = firstdof + nDofs
    end do

  end subroutine atl_linindicator_getElement



  ! This routine returns the sum of squared strain tensor and the vorticity
  ! tensor. So, the output matrix Q_crit
  ! Q_crit(ndof,3,3) = S^2 + \Omega^2 for each ndof
  ! where S = 0.5*(\nabla v + \nabla v^{T}) - strain tensor
  ! and  \Omega = 0.5*(\nabla v - \nabla v^{T}) - vorticity tensor
  subroutine calc_matrix_QCriterion(gradV, Q)
    real(kind=rk),intent(in) :: gradV(:)
    real(kind=rk),intent(inout) :: Q(:,:)
    ! ------------------------------------------------------------------------

    real(kind=rk) :: s(3,3)
    real(kind=rk) :: omg(3,3)
    integer       :: i,j

    ! ------------------------------------------------------------------------
    ! The 9 components of GradV are organized as
    ! GradV(1)   - \frac{\partial v_x}{\partial x}
    ! GradV(2)   - \frac{\partial v_x}{\partial y}
    ! GradV(3)   - \frac{\partial v_x}{\partial z}

    ! GradV(4)   - \frac{\partial v_y}{\partial x}
    ! GradV(5)   - \frac{\partial v_y}{\partial y}
    ! GradV(6)   - \frac{\partial v_y}{\partial z}

    ! GradV(7)   - \frac{\partial v_z}{\partial x}
    ! GradV(8)   - \frac{\partial v_z}{\partial y}
    ! GradV(9)   - \frac{\partial v_z}{\partial z}


    s(1,1) = gradV(1)
    s(1,2) = 0.5*(gradV(2)+ gradV(4))
    s(1,3) = 0.5*(gradV(3)+ gradV(7))
    s(2,1) = 0.5*(gradV(4)+ gradV(2))
    S(2,2) = gradV(5)
    s(2,3) = 0.5*(gradV(6)+ gradV(8))
    s(3,1) = 0.5*(gradV(7)+ gradV(3))
    s(3,2) = 0.5*(gradV(8)+ gradV(6))
    s(3,3) = gradV(9)

    omg(1,1) = 0.0
    omg(1,2) = 0.5*(gradV(2) - gradV(4))
    omg(1,3) = 0.5*(gradV(3) - gradV(7))
    omg(2,1) = 0.5*(gradV(4) - gradV(2))
    omg(2,2) = 0.0
    omg(2,3) = 0.5*(gradV(6) - gradV(8))
    omg(3,1) = 0.5*(gradV(7) - gradV(3))
    omg(3,2) = 0.5*(gradV(8) - gradV(6))
    omg(3,3) = 0.0

    ! Squares the matrix and sums up
    do i = 1,3
      do j = 1,3
         Q(i,j) = s(i,1)*s(1,j)+s(i,2)*s(2,j)+s(i,3)*s(3,j) + &
           & omg(i,1)*omg(1,j)+omg(i,2)*omg(2,j)+omg(i,3)*omg(3,j)
      enddo
    enddo
  end subroutine calc_matrix_QCriterion

  ! This routine calculates the eigenValues of a 3 x 3 symmetric matrix
  ! and returns eigenValues arranged in ascending order
  subroutine calc_eigenValues_3by3_matrix(mat, eig)
    real(kind=rk) :: mat(3,3)
    real(kind=rk) :: eig(3)
    ! ------------------------------------------------------------------------
    real(kind = rk) :: a00, a01, a02, a11, a12, a22
    real(kind = rk) :: c0,c1,c2
    real(kind = rk) :: c2Div3, aDiv3, mbDiv2, q, angle, cs, sn, magnitude
    real(kind = rk) :: inv3, root3
    ! ------------------------------------------------------------------------
    inv3 = 1.0/3.0
    root3 = sqrt(3.0)

    a00 = mat(1,1)
    a01 = mat(1,2)
    a02 = mat(1,3)
    a11 = mat(2,2)
    a12 = mat(2,3)
    a22 = mat(3,3)

    !Coefficients of the cubic equation%
    !  lam^3 - c2*lam^2 + c1*lam -c0 = -det(A-lam*I)
    c0 = a00*a11*a22 + 2*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01
    c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12
    c2 = a00 + a11 + a22

    c2Div3 = c2/3
    aDiv3 = (c1-c2*c2Div3)*inv3

    if (aDiv3 > 0.0) then
       aDiv3 = 0.0
    endif

    mbDiv2 = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1))
    q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3

    if (q> 0.0) then
        q = 0.0
    endif

    magnitude = sqrt(-aDiv3)
    angle = atan2(sqrt(-q),mbDiv2)*inv3
    cs = cos(angle)
    sn = sin(angle)

    eig(1) = c2Div3 + 2*magnitude*cs
    eig(2) = c2Div3 - magnitude*(cs + root3*sn)
    eig(3) = c2Div3 - magnitude*(cs - root3*sn)

   ! Now arrange these eigenValues in Ascending order
   call atl_bubbleSortArray( eig, 3)

  end subroutine calc_eigenValues_3by3_matrix


  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_prim2cons(equation, instate, outstate, material)
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert. If no outstate is provided, the
    !! conversion is done in place.
    real(kind=rk), intent(inout)         :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)        :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in) :: material(:,:)
    ! -------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,1)*instate(:,2)
      outstate(:,3) = instate(:,1)*instate(:,3)
      outstate(:,4) = instate(:,1)*instate(:,4)
      outstate(:,5) = instate(:,5) &
        &                      / (equation%euler%isen_coef-1.0_rk) &
        &                    + 0.5_rk*instate(:,1) &
        &                            * ( instate(:,2)**2 &
        &                               + instate(:,3)**2 &
        &                               + instate(:,4)**2 )
    else
      instate(:,5) = instate(:,5) &
        &                      / (equation%euler%isen_coef-1.0_rk) &
        &                    + 0.5_rk*instate(:,1) &
        &                            * ( instate(:,2)**2 &
        &                               + instate(:,3)**2 &
        &                               + instate(:,4)**2 )
      instate(:,2) = instate(:,1)*instate(:,2)
      instate(:,3) = instate(:,1)*instate(:,3)
      instate(:,4) = instate(:,1)*instate(:,4)
    end if
  end subroutine atl_eqn_euler_prim2cons


  !> Convert primitive varibales to conservative variables including their
  !! gradients.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_prim2cons_grad(equation, instate, outstate, &
    &                                     material                     )
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional, intent(in)    :: material(:,:)
    ! -------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,1) * instate(:,2)
      outstate(:,3) = instate(:,1) * instate(:,3)
      outstate(:,4) = instate(:,1) * instate(:,4)
      outstate(:,5) = instate(:,5) / (equation%euler%isen_coef-1.0_rk)  &
        & + 0.5_rk * instate(:,1) &
        & * ( instate(:,2)**2 + instate(:,3)**2 + instate(:,4)**2)
      outstate(:,6) = instate(:,6)
      outstate(:,7) = instate(:,2) * instate(:,6) &
        &  + instate(:,1) * instate(:,7)
      outstate(:,8) = instate(:,3) * instate(:,6) &
        &  + instate(:,1) * instate(:,8)
      outstate(:,9) = instate(:,4) * instate(:,6) &
        &  + instate(:,1) * instate(:,9)
      outstate(:,10) = 0.5 * (  instate(:,2)**2 + instate(:,3)**2    &
        &                     + instate(:,4)**2                    ) &
        &  * instate(:,6) + instate(:,1)                             &
        &  * ( instate(:,2)*instate(:,7)+instate(:,3)*instate(:,8) + &
        &      instate(:,4) * instate(:,9) )                &
        &  + 1/(equation%euler%isen_coef-1.0_rk) * instate(:,10)
    else
      instate(:,10) = 0.5 * ( instate(:,2)**2 + instate(:,3)**2      &
        &                    + instate(:,4)**2                     ) &
        &  * instate(:,6) + instate(:,1)                             &
        &  * ( instate(:,2)*instate(:,7)+instate(:,3)*instate(:,8) + &
        &      instate(:,4) * instate(:,9) )                &
        &  + 1/(equation%euler%isen_coef-1.0_rk) * instate(:,10)
      instate(:,7) = instate(:,2) * instate(:,6) &
        &  + instate(:,1) * instate(:,7)
      instate(:,8) = instate(:,3) * instate(:,6) &
        &  + instate(:,1) * instate(:,8)
      instate(:,9) = instate(:,4) * instate(:,6) &
        &  + instate(:,1) * instate(:,9)
      instate(:,5) = instate(:,5) / (equation%euler%isen_coef-1.0_rk) &
        & + 0.5_rk*instate(:,1) &
        & * ( instate(:,2)**2 + instate(:,3)**2 + instate(:,4)**2 )
      instate(:,4) = instate(:,1)*instate(:,4)
      instate(:,3) = instate(:,1)*instate(:,3)
      instate(:,2) = instate(:,1)*instate(:,2)
    end if

  end subroutine atl_eqn_euler_prim2cons_grad

  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_cons2prim(equation, instate, outstate, material)
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert. If no outstate is provided, the
    !! conversion is done in place.
    real(kind=rk), intent(inout)         :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)        :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in) :: material(:,:)
    ! -------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,2)/instate(:,1)
      outstate(:,3) = instate(:,3)/instate(:,1)
      outstate(:,4) = instate(:,4)/instate(:,1)
      outstate(:,5) = (equation%euler%isen_coef - 1.0_rk) &
        &           * ( instate(:,5) - 0.5_rk * instate(:,1) &
        &                                     * (outstate(:,2)**2 &
        &                                        + outstate(:,3)**2 &
        &                                        + outstate(:,4)**2) &
        &             )
    else
      instate(:,2) = instate(:,2)/instate(:,1)
      instate(:,3) = instate(:,3)/instate(:,1)
      instate(:,4) = instate(:,4)/instate(:,1)
      instate(:,5) = (equation%euler%isen_coef - 1.0_rk) &
        &           * ( instate(:,5) - 0.5_rk * instate(:,1) &
        &                                     * (instate(:,2)**2 &
        &                                        + instate(:,3)**2 &
        &                                        + instate(:,4)**2) &
        &             )
    end if
  end subroutine atl_eqn_euler_cons2prim

  !> Convert conservative to primitive variables including the gradients.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_cons2prim_grad( equation, instate, outstate, &
    &                                      material                     )
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional, intent(in)    :: material(:,:)
    ! -------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,2)/instate(:,1)
      outstate(:,3) = instate(:,3)/instate(:,1)
      outstate(:,4) = instate(:,4)/instate(:,1)
      outstate(:,5) = (equation%euler%isen_coef - 1.0_rk)            &
        & * ( instate(:,5)                                           &
        &   - 0.5_rk * instate(:,1)                                  &
        &  *(outstate(:,2)**2 + outstate(:,3)**2 + outstate(:,4)**2) )
      ! the gradients
      outstate(:,6) = instate(:,6)
      outstate(:,7) = - instate(:,2)**2 / instate(:,1)**2 &
        &  * instate(:,5)                                 &
        &  + 1/instate(:,1) * instate(:,7)
      outstate(:,8) = - instate(:,3)**2 / instate(:,1)**2 &
        &  * instate(:,5)                                 &
        &  + 1/instate(:,1) * instate(:,8)
      outstate(:,9) = - instate(:,4)**2 / instate(:,1)**2 &
        &  * instate(:,5)                                 &
        &  + 1/instate(:,1) * instate(:,9)
      outstate(:,10) = (equation%euler%isen_coef-1.0_rk) * 0.5       &
        & * (instate(:,2)**2 + instate(:,3)**2 + instate(:,4)**2)    &
        & / instate(:,1)**2 * instate(:,6)                           &
        & + (1.0_rk-equation%euler%isen_coef)/instate(:,1)           &
        & * ( instate(:,2)*instate(:,7) + instate(:,3)*instate(:,8)  &
        &     + instate(:,4)*instate(:,9)                          ) &
        & + (equation%euler%isen_coef-1.0_rk) * instate(:,10)
    else
      ! the gradients
      instate(:,10) = (equation%euler%isen_coef-1.0_rk) * 0.5        &
        & * (instate(:,2)**2 + instate(:,3)**2 + instate(:,4)**2)    &
        & / instate(:,1)**2 * instate(:,6)                           &
        & + (1.0_rk-equation%euler%isen_coef)/instate(:,1)           &
        & * ( instate(:,2)*instate(:,7) + instate(:,3)*instate(:,8)  &
        &     + instate(:,4)*instate(:,9)                          ) &
        & + (equation%euler%isen_coef-1.0_rk) * instate(:,10)
      outstate(:,9) = - instate(:,4)**2 / instate(:,1)**2 &
        &  * instate(:,5)                                 &
        &  + 1/instate(:,1) * instate(:,9)
      outstate(:,8) = - instate(:,3)**2 / instate(:,1)**2 &
        &  * instate(:,5)                                 &
        &  + 1/instate(:,1) * instate(:,8)
      outstate(:,7) = - instate(:,2)**2 / instate(:,1)**2 &
        &  * instate(:,5)                                 &
        &  + 1/instate(:,1) * instate(:,7)
      ! the state
      instate(:,2) = instate(:,2)/instate(:,1)
      instate(:,3) = instate(:,3)/instate(:,1)
      instate(:,4) = instate(:,4)/instate(:,1)
      instate(:,5) = (equation%euler%isen_coef - 1.0_rk) &
        & * ( instate(:,5)                               &
        &   - 0.5_rk * instate(:,1) * ( instate(:,2)**2  &
        &   + instate(:,3)**2 + instate(:,4)**2 )        )
    end if
  end subroutine atl_eqn_euler_cons2prim_grad

  !> Convert conservative to primitive variables (including temperature
  !! instead of pressure).
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_cons2primTemp(equation, instate, outstate, material)
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)        :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in) :: material(:,:)
    ! -------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,2)/instate(:,1)
      outstate(:,3) = instate(:,3)/instate(:,1)
      outstate(:,4) = instate(:,4)/instate(:,1)
      ! ... pressure
      outstate(:,5) = (equation%euler%isen_coef - 1.0_rk) &
        &           * ( instate(:,5) - 0.5_rk * instate(:,1) &
        &                                     * (outstate(:,2)**2 &
        &                                        + outstate(:,3)**2  &
        &                                        + outstate(:,4)**2 ) &
        &             )
       ! ... temperature
       outstate(:,5) = outstate(:,5) / outstate(:,1) / equation%euler%r
     else
       instate(:,2) = instate(:,2)/instate(:,1)
       instate(:,3) = instate(:,3)/instate(:,1)
       instate(:,4) = instate(:,4)/instate(:,1)
       ! ... pressure
       instate(:,5) = (equation%euler%isen_coef - 1.0_rk) &
         &           * ( instate(:,5) - 0.5_rk * instate(:,1) &
         &                                     * ( instate(:,2)**2 &
         &                                        + instate(:,3)**2  &
         &                                        + instate(:,4)**2 ) &
         &             )
       ! ... temperature
       instate(:,5) = instate(:,5) / instate(:,1) / equation%euler%r
     end if
  end subroutine atl_eqn_euler_cons2primTemp

  !> Convert primitive varibales (including temperature instead
  !! of pressure) to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_primTemp2cons(equation, instate, outstate, material)
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)          :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)  :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in)  :: material(:,:)
    ! -------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,1)*instate(:,2)
      outstate(:,3) = instate(:,1)*instate(:,3)
      outstate(:,4) = instate(:,1)*instate(:,4)
      ! ... convert from temperature to pressure (p = rho * R * T)
      outstate(:,5) = instate(:,1) * instate(:,5) * equation%euler%r
      ! ... total energy
      outstate(:,5) = (outstate(:,5) &
        &                      / (equation%euler%isen_coef-1.0_rk)) &
        &                    + 0.5_rk*instate(:,1) &
        &                            * ( instate(:,2)**2 &
        &                               + instate(:,3)**2 &
        &                               + instate(:,4)**2 )
    else
      ! ... convert from temperature to pressure (p = rho * R * T)
      instate(:,5) = instate(:,1) * instate(:,5) * equation%euler%r
      ! ... total energy
      instate(:,5) = instate(:,5) &
        &                      / (equation%euler%isen_coef-1.0_rk) &
        &                    + 0.5_rk*instate(:,1) &
        &                            * ( instate(:,2)**2 &
        &                               + instate(:,3)**2 &
        &                               + instate(:,4)**2 )
      instate(:,2) = instate(:,1)*instate(:,2)
      instate(:,3) = instate(:,1)*instate(:,3)
      instate(:,4) = instate(:,1)*instate(:,4)
    end if
  end subroutine atl_eqn_euler_primTemp2cons

  !> Convert conservative to conservative variables (including velocity
  !! instead of momentum).
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_cons2primVel(equation, instate, outstate, material)
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)        :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in) :: material(:,:)
    ! -------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,2)/instate(:,1)
      outstate(:,3) = instate(:,3)/instate(:,1)
      outstate(:,4) = instate(:,4)/instate(:,1)
      outstate(:,5) = instate(:,5)
    else
      instate(:,2) = instate(:,2)/instate(:,1)
      instate(:,3) = instate(:,3)/instate(:,1)
      instate(:,4) = instate(:,4)/instate(:,1)
    end if

  end subroutine atl_eqn_euler_cons2primVel

  !> Convert conservative varibales (including velocity instead
  !! of temperature) to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_primVel2cons(equation, instate, outstate, material)
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)        :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in) :: material(:,:)
    ! -------------------------------------------------------------------------!

    if(present(outstate)) then
      outstate(:,1) = instate(:,1)
      outstate(:,2) = instate(:,1)*instate(:,2)
      outstate(:,3) = instate(:,1)*instate(:,3)
      outstate(:,4) = instate(:,1)*instate(:,4)
      outstate(:,5) = instate(:,5)
    else
      instate(:,2) = instate(:,1)*instate(:,2)
      instate(:,3) = instate(:,1)*instate(:,3)
      instate(:,4) = instate(:,1)*instate(:,4)
    end if

  end subroutine atl_eqn_euler_primVel2cons

  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_prim2cons_elems(equation, instate, outstate, nElems)
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:,:)

    !> Number of elements to act on (first index in the state arrays).
    integer, intent(in)               :: nElems
    ! -------------------------------------------------------------------------!

    if (present(outstate)) then
      outstate(1:nElems,:,1) = instate(1:nElems,:,1)
      outstate(1:nElems,:,2) = instate(1:nElems,:,1)*instate(1:nElems,:,2)
      outstate(1:nElems,:,3) = instate(1:nElems,:,1)*instate(1:nElems,:,3)
      outstate(1:nElems,:,4) = instate(1:nElems,:,1)*instate(1:nElems,:,4)
      outstate(1:nElems,:,5) = instate(1:nElems,:,5) &
        &                      / (equation%euler%isen_coef-1.0_rk) &
        &                    + 0.5_rk*instate(1:nElems,:,1) &
        &                            * ( instate(1:nElems,:,2)**2 &
        &                               + instate(1:nElems,:,3)**2 &
        &                               + instate(1:nElems,:,4)**2 )
    else
      instate(1:nElems,:,5) = instate(1:nElems,:,5) &
        &                      / (equation%euler%isen_coef-1.0_rk) &
        &                    + 0.5_rk*instate(1:nElems,:,1) &
        &                            * ( instate(1:nElems,:,2)**2 &
        &                               + instate(1:nElems,:,3)**2 &
        &                               + instate(1:nElems,:,4)**2 )
      instate(1:nElems,:,2) = instate(1:nElems,:,1)*instate(1:nElems,:,2)
      instate(1:nElems,:,3) = instate(1:nElems,:,1)*instate(1:nElems,:,3)
      instate(1:nElems,:,4) = instate(1:nElems,:,1)*instate(1:nElems,:,4)
    end if
  end subroutine atl_eqn_euler_prim2cons_elems


  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_euler_cons2prim_elems(equation, instate, outstate, nElems)
    ! -------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in)  :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:,:)

    !> Number of elements to act on (first index in the state arrays).
    integer,                 intent(in)    :: nElems
    ! -------------------------------------------------------------------------!

    if (present(outstate)) then
      outstate(1:nElems,:,1) = instate(1:nElems,:,1)
      outstate(1:nElems,:,2) = instate(1:nElems,:,2) / instate(1:nElems,:,1)
      outstate(1:nElems,:,3) = instate(1:nElems,:,3) / instate(1:nElems,:,1)
      outstate(1:nElems,:,4) = instate(1:nElems,:,4) / instate(1:nElems,:,1)
      outstate(1:nElems,:,5) = (equation%euler%isen_coef - 1.0_rk) &
        & * ( instate(1:nElems,:,5)                                &
        &   - 0.5_rk * instate(1:nElems,:,1)                       &
        &     * (outstate(1:nElems,:,2)**2                         &
        &       + outstate(1:nElems,:,3)**2                        &
        &       + outstate(1:nElems,:,4)**2))
    else
      instate(1:nElems,:,2) = instate(1:nElems,:,2) / instate(1:nElems,:,1)
      instate(1:nElems,:,3) = instate(1:nElems,:,3) / instate(1:nElems,:,1)
      instate(1:nElems,:,4) = instate(1:nElems,:,4) / instate(1:nElems,:,1)
      instate(1:nElems,:,5) = (equation%euler%isen_coef - 1.0_rk)        &
        & * ( instate(1:nElems,:,5)                                      &
        &   - 0.5_rk * instate(1:nElems,:,1) * (instate(1:nElems,:,2)**2 &
        &   + instate(1:nElems,:,3)**2                                   &
        &   + instate(1:nElems,:,4)**2) )
    end if
  end subroutine atl_eqn_euler_cons2prim_elems


end module atl_eqn_euler_derive_module
