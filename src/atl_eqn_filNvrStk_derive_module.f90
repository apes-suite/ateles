! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
module atl_eqn_filNvrStk_derive_module
  use, intrinsic :: iso_c_binding,  only: c_f_pointer
  use env_module,                   only: rk
  use tem_varSys_module,            only: tem_varSys_type,   &
    &                                     tem_varSys_op_type

  use atl_derive_module,            only: atl_derive_inputVar_type,        &
    &                                     atl_derive_fromModalData,        &
    &                                     atl_generic_fromModal_getElement
  use atl_varSys_module,            only: atl_varSys_data_type
  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type
  use atl_equation_module,          only: atl_equations_type

  implicit none

  private

  public :: atl_Rans_pressure_getPoint
  public :: atl_Rans_pressure_getElement
  public :: atl_eqn_rans_cons2prim_elems
  public :: atl_eqn_rans_prim2cons_elems
  public :: atl_eqn_rans_2d_cons2prim_elems
  public :: atl_eqn_rans_2d_prim2cons_elems

contains



  subroutine atl_Rans_pressure_getPoint( fun, varsys, point, time,tree, nPnts, &
    &                                    res )
    ! ---------------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: density(nPnts), momentum(3*nPnts), energy(nPnts)
    real(kind=rk) :: turbulent_KE(nPnts)
    ! ---------------------------------------------------------------------------

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

    call varSys%method%val(fun%input_varPos(4))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = turbulent_KE                             )

    res = (fPtr%solverData%equationPtr%Euler%isen_coef - 1.0_rk) &
      &   * (energy                                              &
      & - 0.5_rk                                                 &
      &   / density                                              &
      &   * sum(momentum**2)                                     &
      &   - turbulent_KE                                         )

  end subroutine atl_Rans_pressure_getPoint

  subroutine atl_derivePressure(fun, varsys, tree, iElem, elemPos,    &
    &                                        nodalInput, nodalRes     )
    ! ---------------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------------
    type(atl_varSys_data_type), pointer :: fPtr
    integer, parameter :: density=1, momentum=2, energy=3, turbulent_KE=4
    ! ---------------------------------------------------------------------------

    call C_F_POINTER( fun%method_Data, fPtr )

    nodalRes(:,1) = (fPtr%solverData%equationPtr%Euler%isen_coef - 1.0_rk) &
      &   * (nodalInput(energy)%data(:,1)                                  &
      & - 0.5_rk                                                           &
      &   / nodalInput(density)%data(:,1)                                  &
      &   * sum(array=nodalInput(momentum)%data**2,dim=2)                  &
      &   - nodalInput(turbulent_KE)%data(:,1)                             )

  end subroutine atl_derivePressure

  subroutine atl_Rans_pressure_getElement(fun, varsys, elempos, time, tree,  &
    &                                     nElems, nDofs, res                 )
    ! ---------------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------------
    procedure(atl_derive_fromModalData), pointer :: fnCalcPtr
    type(atl_varSys_data_type), pointer :: fPtr
    ! ---------------------------------------------------------------------------
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

  end subroutine atl_Rans_pressure_getElement

  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_rans_cons2prim_elems(equation, instate, outstate, nElems)
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type),   intent(in)    :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:,:)

    !> Number of elements to act on (first index in the state arrays).
    integer,                 intent(in)    :: nElems
    ! --------------------------------------------------------------------------!

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
        &       + outstate(1:nElems,:,4)**2) - instate(1:nElems,:,6))
      outstate(1:nElems,:,6) = instate(1:nElems,:,6)
      outstate(1:nElems,:,7) = instate(1:nElems,:,7)
    else
      instate(1:nElems,:,2) = instate(1:nElems,:,2) / instate(1:nElems,:,1)
      instate(1:nElems,:,3) = instate(1:nElems,:,3) / instate(1:nElems,:,1)
      instate(1:nElems,:,4) = instate(1:nElems,:,4) / instate(1:nElems,:,1)
      instate(1:nElems,:,5) = (equation%euler%isen_coef - 1.0_rk)        &
        & * ( instate(1:nElems,:,5)                                      &
        &   - 0.5_rk * instate(1:nElems,:,1) * (instate(1:nElems,:,2)**2 &
        &   + instate(1:nElems,:,3)**2                                   &
        &   + instate(1:nElems,:,4)**2)  - instate(1:nElems,:,6))
    end if
  end subroutine atl_eqn_rans_cons2prim_elems


! ******************************************************************************!

  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_rans_prim2cons_elems(equation, instate, outstate, nElems)
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:,:)

    !> Number of elements to act on (first index in the state arrays).
    integer, intent(in)               :: nElems
    ! --------------------------------------------------------------------------!

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
        &                               + instate(1:nElems,:,4)**2 )&
        &                    + instate(1:nElems,:,6)
      outstate(1:nElems,:,6) = instate(1:nElems,:,6)
      outstate(1:nElems,:,7) = instate(1:nElems,:,7)
    else
      instate(1:nElems,:,5) = instate(1:nElems,:,5)                  &
        &                      / (equation%euler%isen_coef-1.0_rk)   &
        &                    + 0.5_rk*instate(1:nElems,:,1)          &
        &                            * ( instate(1:nElems,:,2)**2    &
        &                               + instate(1:nElems,:,3)**2   &
        &                               + instate(1:nElems,:,4)**2 ) &
        &                    + instate(1:nElems,:,6)
      instate(1:nElems,:,2) = instate(1:nElems,:,1)*instate(1:nElems,:,2)
      instate(1:nElems,:,3) = instate(1:nElems,:,1)*instate(1:nElems,:,3)
      instate(1:nElems,:,4) = instate(1:nElems,:,1)*instate(1:nElems,:,4)
    end if
  end subroutine atl_eqn_rans_prim2cons_elems


  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_rans_2d_cons2prim_elems(equation, instate, outstate, nElems)
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert. If outstate is not provided, conversion
    !! will take in place of instate.
    real(kind=rk),           intent(inout) :: instate(:,:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out)   :: outstate(:,:,:)

    !> Number of elements to act on (first index in the state arrays).
    integer,                 intent(in)    :: nElems
    ! --------------------------------------------------------------------------!

    if (present(outstate)) then
      outstate(1:nElems,:,1) = instate(1:nElems,:,1)
      outstate(1:nElems,:,2) = instate(1:nElems,:,2) / instate(1:nElems,:,1)
      outstate(1:nElems,:,3) = instate(1:nElems,:,3) / instate(1:nElems,:,1)
      outstate(1:nElems,:,4) = (equation%euler%isen_coef - 1.0_rk) &
        & * ( instate(1:nElems,:,4)                                &
        &   - 0.5_rk * instate(1:nElems,:,1)                       &
        &     * (outstate(1:nElems,:,2)**2                         &
        &       + outstate(1:nElems,:,3)**2) - instate(1:nElems,:,5))
      outstate(1:nElems,:,5) = instate(1:nElems,:,5)
      outstate(1:nElems,:,6) = instate(1:nElems,:,6)
    else
      instate(1:nElems,:,2) = instate(1:nElems,:,2) / instate(1:nElems,:,1)
      instate(1:nElems,:,3) = instate(1:nElems,:,3) / instate(1:nElems,:,1)
      instate(1:nElems,:,4) = (equation%euler%isen_coef - 1.0_rk)        &
        & * ( instate(1:nElems,:,4)                                      &
        &   - 0.5_rk * instate(1:nElems,:,1) * (instate(1:nElems,:,2)**2 &
        &   + instate(1:nElems,:,3)**2) - instate(1:nElems,:,5) )
    end if
  end subroutine atl_eqn_rans_2d_cons2prim_elems


! ******************************************************************************!

  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_rans_2d_prim2cons_elems(equation, instate, outstate, nElems)
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:,:)

    !> Number of elements to act on (first index in the state arrays).
    integer, intent(in)               :: nElems
    ! --------------------------------------------------------------------------!

    if (present(outstate)) then
      outstate(1:nElems,:,1) = instate(1:nElems,:,1)
      outstate(1:nElems,:,2) = instate(1:nElems,:,1)*instate(1:nElems,:,2)
      outstate(1:nElems,:,3) = instate(1:nElems,:,1)*instate(1:nElems,:,3)
      ! e = p/(\gamma-1) + \rho * v*v*0.5 + \rho * k !
      outstate(1:nElems,:,4) = instate(1:nElems,:,4) &
        &                      / (equation%euler%isen_coef-1.0_rk) &
        &                    + 0.5_rk*instate(1:nElems,:,1) &
        &                            * ( instate(1:nElems,:,2)**2 &
        &                               + instate(1:nElems,:,3)**2 ) &
        &                    + instate(1:nElems,:,5)
      outstate(1:nElems,:,5) = instate(1:nElems,:,5)
      outstate(1:nElems,:,6) = instate(1:nElems,:,6)
    else
      instate(1:nElems,:,4) = instate(1:nElems,:,4) &
        &                      / (equation%euler%isen_coef-1.0_rk) &
        &                    + 0.5_rk*instate(1:nElems,:,1) &
        &                            * ( instate(1:nElems,:,2)**2 &
        &                               + instate(1:nElems,:,3)**2 ) &
        &                    + instate(1:nElems,:,5)
      instate(1:nElems,:,2) = instate(1:nElems,:,1)*instate(1:nElems,:,2)
      instate(1:nElems,:,3) = instate(1:nElems,:,1)*instate(1:nElems,:,3)
    end if
  end subroutine atl_eqn_rans_2d_prim2cons_elems


end module atl_eqn_filNvrStk_derive_module
