! Copyright (c) 2012-2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013-2014,2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Jan Hueckelheim, Melven
! Zoellner and Harald Klimach for German Research School for Simulation
! Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Verena Krupp, Peter Vitt,
! Daniel Fleischer and Tobias Girresser for University of Siegen.
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

! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt and Harald Klimach for
! University of Siegen.
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
!
! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the number of degrees of freedom for Q polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * A variable to store the number of degrees of freedom for a P tensor product
!   polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.

! The x ansatz degree is turned into the degree of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x ansatz degree is turned into the degree of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
!> Module provides subroutines, functions and datatypes regarding
!! cell local degrees of freedoms.
module ply_dof_module
  use env_module, only: rk

  implicit none

  private

  public :: ply_dof_2degree, ply_degree_2dof, ply_change_poly_space

  !> Parameter to identify Q polynomials
  integer, public, parameter :: Q_space = 1
  !> Parameter to identify P polynomials
  integer, public, parameter :: P_space = 2


contains


  elemental function ply_dof_2degree(ndofs, space, ndims) result(deg)
    integer, intent(in) :: ndofs
    integer, intent(in) :: space
    integer, intent(in) :: ndims
    integer :: deg

    integer :: estimate

    select case(space)
    case (Q_space)
      deg = nint(ndofs**(1._rk/real(ndims,kind=rk))) - 1
    case (P_space)
      deg = 0
      do
        estimate = ply_degree_2dof(deg, space, nDims)
        if (estimate >= ndofs) then
          EXIT
        end if
        deg = deg + 1
      end do
    end select

  end function ply_dof_2degree

  elemental function ply_degree_2dof(deg, space, nDims) result(nDofs)
    integer, intent(in) :: deg
    integer, intent(in) :: space
    integer, intent(in) :: nDims
    integer :: nDofs

    ndofs = -1

    select case(space)
    case (Q_space)
      nDofs = (deg+1)**nDims
    case (P_space)
      select case (nDims)
      case(3)
        nDofs = ((deg + 1)  &
          &   *  (deg + 2)  &
          &   *  (deg + 3)) &
          &   / 6
      case(2)
        nDofs = ((deg + 1)  &
          &   *  (deg + 2)) &
          &   / 2
      case(1)
        nDofs = (deg + 1)
      end select
    end select

  end function ply_degree_2dof

  !> Subroutine to change the polynomial space (Q or P) of an
  !! atl_statedata_type from Q-space to P-space and vice versa.
  !  The space of the instate (inspace) defines the space of the
  !  outstate.
  subroutine ply_change_poly_space( inspace, instate, outstate,      &
    &                               maxPolyDeg, nElems, nVars, nDims )
    ! -------------------------------------------------------------------- !
    !> Polynomial space of the input state (P_sapce or Q_space)
    integer, intent(in) :: inspace

    !> States of the variables of the input in polynomial space as
    !! prescribed in inspace.
    real(kind=rk), intent(in) :: instate(:,:,:)

    !> States of the variables of the output.
    real(kind=rk), intent(inout) :: outstate(:,:,:)

    integer, intent(in) :: maxPolyDeg

    integer, intent(in) :: nElems

    integer, intent(in) :: nVars

    integer, intent(in) :: nDims
    ! -------------------------------------------------------------------- !
    integer :: iElem, iVar, iAnsX, iAnsY, iAnsZ
    integer :: P_pos, Q_pos
    ! -------------------------------------------------------------------- !

    select case(nDims)
    case(3)
      select case(inspace)
      ! Instate space is P_space so outstate space is Q_space
      ! Copy the dofs in the right order and fill up the higher modes with zeros
      case(P_space)
        outstate = 0.0_rk
        do iElem = 1, nElems
          do iVar = 1, nVars
            do iAnsZ = 1, maxPolyDeg+1
              do iAnsY = 1, maxPolyDeg+1 - (iAnsZ-1)
                do iAnsX = 1, maxPolyDeg+1 - (iAnsZ-1) - (iAnsY-1)
  ! integer divisions are no mistake here.
  p_pos = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
  q_pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
                  outstate( iElem, Q_pos, iVar ) = instate(iElem, P_pos, iVar)
                end do
              end do
            end do
          end do
        end do

      ! Instate space is Q_space so outstate space is P_space
      ! Copy the dofs in the right order and cut off the higher modes
      case(Q_space)
        do iElem = 1, nElems
          do iVar = 1, nVars
            do iAnsZ = 1, maxPolyDeg+1
              do iAnsY = 1, maxPolyDeg+1 - (iAnsZ-1)
                do iAnsX = 1, maxPolyDeg+1 - (iAnsZ-1) - (iAnsY-1)
  ! integer divisions are no mistake here.
  p_pos = (((iansx + iansy + iansz - 3) &
    &     * (iansx + iansy + iansz - 2) &
    &     * (iansx + iansy + iansz - 1)) &
    &   / 6 + 1)             &
    & + ((iansz-1) * (iansx + iansy + iansz -2) &
    &   - ((iansz-2) * (iansz-1)) / 2) &
    & + (iansy-1)
  q_pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
                  outstate( iElem, P_pos, iVar ) = instate(iElem, Q_pos, iVar)
                end do
              end do
            end do
          end do
        end do
      end select

    case(2)
      select case(inspace)
      ! Instate space is P_space so outstate space is Q_space
      ! Copy the dofs in the right order and fill up the higher modes with zeros
      case(P_space)
        outstate = 0.0_rk
        do iElem = 1, nElems
          do iVar = 1, nVars
            do iAnsY = 1, maxPolyDeg+1
              do iAnsX = 1, maxPolyDeg+1 - (iAnsY-1)
  ! integer divisions are no mistake here.
  p_pos = ((((iansx - 1) + (iansy - 1))            &
    &   * (((iansx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
  q_pos = iansx + (iansy-1)*(maxpolydeg+1)
                outstate( iElem, Q_pos, iVar ) = instate(iElem, P_pos, iVar)
              end do
            end do
          end do
        end do

      ! Instate space is Q_space so outstate space is P_space
      ! Copy the dofs in the right order and cut off the higher modes
      case(Q_space)
        do iElem = 1, nElems
          do iVar = 1, nVars
            do iAnsY = 1, maxPolyDeg+1
              do iAnsX = 1, maxPolyDeg+1 - (iAnsY-1)
  ! integer divisions are no mistake here.
  p_pos = ((((iansx - 1) + (iansy - 1))            &
    &   * (((iansx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
  q_pos = iansx + (iansy-1)*(maxpolydeg+1)
                outstate( iElem, P_pos, iVar ) = instate(iElem, Q_pos, iVar)
              end do
            end do
          end do
        end do
      end select
    end select

  end subroutine ply_change_poly_space

end module ply_dof_module
