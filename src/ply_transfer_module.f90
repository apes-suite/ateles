! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Peter Vitt,
! Tobias Girresser and Daniel Petró for University of Siegen.
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
!> This module provides routines to transfer degrees of freedom from one
!! polynomial representation to another.
!!
!! These routines allow the copying of data from one order and polynomial space
!! into another.
module ply_transfer_module
  use env_module,     only: rk

  use ply_dof_module, only: P_Space, &
    &                       Q_space

  implicit none

  private

  public :: ply_transfer_dofs_1D
  public :: ply_transfer_dofs_2D
  public :: ply_transfer_dofs_3D
  public :: ply_transfer_dofs
  public :: ply_transfer_P_dim


contains


  ! ************************************************************************ !
  !> Transfer of degrees of freedom from one polynomial to another in 1D.
  !!
  !! If the indat is larger than outdat, the higher modes are truncated.
  !! If outdat is larger, higher modes are padded with zeros.
  subroutine ply_transfer_dofs_1D( indat, indegree, outdat, outdegree )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Degree of the input polynomial. There are indegree+1 modes expected
    !! in indat.
    integer, intent(in) :: indegree

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Degree of the output polynomial. There are outdegree+1 modes expected
    !! in outdat.
    integer, intent(in) :: outdegree
    ! -------------------------------------------------------------------- !
    integer :: minOrd
    ! -------------------------------------------------------------------- !

    minord = min(outdegree+1, indegree+1)

    outdat(:minord) = indat(:minord)
    outdat(minord+1:) = 0.0_rk

  end subroutine ply_transfer_dofs_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transfer of degrees of freedom from one polynomial to another in 2D.
  !!
  !! If the indat is larger than outdat, the higher modes are truncated.
  !! If outdat is larger, higher modes are padded with zeros.
  !!
  !! When different multidimensional polynomial layouts are used, the modes
  !! are copied to to the corresponding locations.
  subroutine ply_transfer_dofs_2D( indat, inspace, indegree,   &
    &                              outdat, outspace, outdegree )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Multi-dimensional polynomial layout of the input data.
    !!
    !! Has to be either [[ply_dof_module:Q_Space]]
    !! or [[ply_dof_module:P_Space]].
    integer, intent(in) :: inspace

    !> Maximal polynomial degree in the input data.
    integer, intent(in) :: indegree

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Multi-dimensional polynomial layout of the output data.
    !!
    !! Has to be either [[ply_dof_module:Q_Space]]
    !! or [[ply_dof_module:P_Space]].
    integer, intent(in) :: outspace

    !> Maximal polynomial degree in the output data.
    integer, intent(in) :: outdegree
    ! -------------------------------------------------------------------- !
    integer :: outdofs
    integer :: indofs
    integer :: min_dofs
    integer :: iStep, iDof
    integer :: out_X, out_Y
    integer :: in_X, in_Y
    integer :: out_pos, in_pos
    integer :: out_off
    integer :: in_off
    integer :: minOrd
    ! -------------------------------------------------------------------- !

    minord = min(outdegree+1, indegree+1)

    select case(inspace)
    case (Q_Space)
      indofs = (indegree+1)**2

    case (P_Space)
      indofs = ((indegree+1)*(indegree+2)) / 2

    end select

    outdat = 0.0_rk

    ospace: select case(outspace)
    case (Q_Space) ospace
      outdofs = (outdegree+1)**2

      ispace_oq: if (inspace == Q_Space) then

        ! Both, output and input are Q Polynomials
        do out_Y=0,minord-1
          out_off = out_Y*(outdegree+1)
          in_off = out_Y*(indegree+1)
          outdat(out_off+1:out_off+minord) = indat(in_off+1:in_off+minord)
        end do

      else ispace_oq

        ! Output is Q, but input is P
        in_X = 1
        in_Y = 1
        do iDof=1,indofs
  ! integer divisions are no mistake here.
  in_pos = ((((in_x - 1) + (in_y - 1))            &
    &   * (((in_x - 1) + (in_y - 1)) + 1)) / 2 + 1) &
    & + (in_y - 1)
          out_pos = in_X + (outdegree+1)*(in_Y-1)
          outdat(out_pos) = indat(in_pos)
          ! Ensure, that next iteration is in the minord range
          do istep=iDof,indofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (in_x .ne. 1) then
    ! next item
    in_x = in_x - 1
    in_y = in_y + 1
  else
    ! next layer
    in_x = in_y + 1
    in_y = 1
  end if
            if ((in_X <= minord) .and. (in_Y <= minord)) EXIT
          end do
          if ( (in_X > minord) .or. (in_Y > minord) &
             & .or. (in_X+in_Y-2 > indegree)        ) EXIT
        end do

      end if ispace_oq


    case (P_Space) ospace
      outdofs = ((outdegree+1)*(outdegree+2)) / 2

      ispace_op: if (inspace == Q_Space) then

        ! Output is P, input is Q
        out_X = 1
        out_Y = 1
        do iDof=1,outdofs
  ! integer divisions are no mistake here.
  out_pos = ((((out_x - 1) + (out_y - 1))            &
    &   * (((out_x - 1) + (out_y - 1)) + 1)) / 2 + 1) &
    & + (out_y - 1)
          in_pos = out_X + (indegree+1)*(out_Y-1)
          outdat(out_pos) = indat(in_pos)
          ! Ensure, that next iteration is in the target range
          do istep=iDof,outdofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (out_x .ne. 1) then
    ! next item
    out_x = out_x - 1
    out_y = out_y + 1
  else
    ! next layer
    out_x = out_y + 1
    out_y = 1
  end if
            if ((out_X <= minord) .and. (out_Y <= minord)) EXIT
          end do
          if ( (out_X > minord) .or. (out_Y > minord) &
            &  .or. (out_X+out_Y-2 > outdegree)       ) EXIT
        end do

      else ispace_op

        ! Both input and output are P polynomials
        min_dofs = (minord*(minord+1))/2
        outdat(:min_dofs) = indat(:min_dofs)

      end if ispace_op

    end select ospace

  end subroutine ply_transfer_dofs_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transfer of degrees of freedom from one polynomial to another in 3D.
  !!
  !! If the indat is larger than outdat, the higher modes are truncated.
  !! If outdat is larger, higher modes are padded with zeros.
  !!
  !! When different multidimensional polynomial layouts are used, the modes
  !! are copied to to the corresponding locations.
  subroutine ply_transfer_dofs_3D( indat, inspace, indegree,   &
    &                              outdat, outspace, outdegree )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Multi-dimensional polynomial layout of the input data.
    !!
    !! Has to be either [[ply_dof_module:Q_Space]]
    !! or [[ply_dof_module:P_Space]].
    integer, intent(in) :: inspace

    !> Maximal polynomial degree in the input data.
    integer, intent(in) :: indegree

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Multi-dimensional polynomial layout of the output data.
    !!
    !! Has to be either [[ply_dof_module:Q_Space]]
    !! or [[ply_dof_module:P_Space]].
    integer, intent(in) :: outspace

    !> Maximal polynomial degree in the output data.
    integer, intent(in) :: outdegree
    ! -------------------------------------------------------------------- !
    integer :: outdofs
    integer :: indofs
    integer :: min_dofs
    integer :: iStep, iDof
    integer :: out_X, out_Y, out_Z
    integer :: in_X, in_Y, in_Z
    integer :: out_pos, in_pos
    integer :: out_off, out_zoff
    integer :: in_off, in_zoff
    integer :: minOrd
    ! -------------------------------------------------------------------- !

    minord = min(outdegree+1, indegree+1)

    select case(inspace)
    case (Q_Space)
      indofs = (indegree+1)**3

    case (P_Space)
      indofs = ((indegree+1)*(indegree+2)*(indegree+3))/6

    end select

    outdat = 0.0_rk

    ospace: select case(outspace)
    case (Q_Space) ospace
      outdofs = (outdegree+1)**3
      ispace_oq: if (inspace == Q_Space) then

        ! Both, output and input are Q Polynomials
        do out_Z=0,minord-1
          out_zoff = out_Z*(outdegree+1)**2
          in_zoff = out_Z*(indegree+1)**2
          do out_Y=0,minord-1
            out_off = out_Y*(outdegree+1) + out_zoff
            in_off = out_Y*(indegree+1) + in_zoff
            outdat(out_off+1:out_off+minord) &
              & = indat(in_off+1:in_off+minord)
          end do
        end do

      else ispace_oq

        ! Output is Q, but input is P
        in_X = 1
        in_Y = 1
        in_Z = 1
        do iDof=1,indofs
  ! integer divisions are no mistake here.
  in_pos = (((in_x + in_y + in_z - 3) &
    &     * (in_x + in_y + in_z - 2) &
    &     * (in_x + in_y + in_z - 1)) &
    &   / 6 + 1)             &
    & + ((in_z-1) * (in_x + in_y + in_z -2) &
    &   - ((in_z-2) * (in_z-1)) / 2) &
    & + (in_y-1)
          out_pos = in_X + (outdegree+1)*( (in_Y-1) &
            &                             + (outdegree+1)*(in_Z-1))
          outdat(out_pos) = indat(in_pos)
          ! Ensure, that next iteration is in the target range
          do istep=iDof,indofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (in_x .ne. 1) then
    ! next item
    in_x = in_x - 1
    in_y = in_y + 1
  elseif (in_y .ne. 1) then
    ! next block
    in_x = in_y - 1
    in_y = 1
    in_z = in_z + 1
  else
    ! next layer
    in_x = in_z + 1
    in_y = 1
    in_z = 1
  end if
            if ( (in_X <= minord) .and. (in_Y <= minord) &
              &                   .and. (in_Z <= minord) ) EXIT
          end do
          if ( (in_X > minord) .or. (in_Y > minord) &
            &                  .or. (in_Z > minord) &
            &  .or. (in_X+in_Y+in_Z-3 > indegree)   ) EXIT
        end do

      end if ispace_oq

    case (P_Space) ospace
      outdofs = ((outdegree+1)*(outdegree+2)*(outdegree+3))/6

      ispace_op: if (inspace == Q_Space) then
        ! Output is P, input is Q
        out_X = 1
        out_Y = 1
        out_Z = 1
        do iDof=1,outdofs
  ! integer divisions are no mistake here.
  out_pos = (((out_x + out_y + out_z - 3) &
    &     * (out_x + out_y + out_z - 2) &
    &     * (out_x + out_y + out_z - 1)) &
    &   / 6 + 1)             &
    & + ((out_z-1) * (out_x + out_y + out_z -2) &
    &   - ((out_z-2) * (out_z-1)) / 2) &
    & + (out_y-1)
          in_pos = out_X + (indegree+1)*( (out_Y-1) &
            &                                + (indegree+1)*(out_Z-1))
          outdat(out_pos) = indat(in_pos)
          ! Ensure, that next iteration is in the target range
          do istep=iDof,outdofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (out_x .ne. 1) then
    ! next item
    out_x = out_x - 1
    out_y = out_y + 1
  elseif (out_y .ne. 1) then
    ! next block
    out_x = out_y - 1
    out_y = 1
    out_z = out_z + 1
  else
    ! next layer
    out_x = out_z + 1
    out_y = 1
    out_z = 1
  end if
            if ((out_X <= minord) .and. (out_Y <= minord) &
              &                   .and. (out_Z <= minord) ) EXIT
          end do
          if ( (out_X > minord) .or. (out_Y > minord) &
            &                   .or. (out_Z > minord) &
            &  .or. (out_X+out_Y+out_Z-3 > outdegree) ) EXIT
        end do

      else ispace_op

        ! Both input and output are P polynomials
        min_dofs = ( (minord+2)*((minord*(minord+1))/2) ) / 3
        outdat(:min_dofs) = indat(:min_dofs)

      end if ispace_op

    end select ospace

  end subroutine ply_transfer_dofs_3d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Small helping routine to wrap transfers in all allowed dimensions.
  subroutine ply_transfer_dofs( indat, inspace, indegree,          &
    &                           outdat, outspace, outdegree, ndims )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Multi-dimensional polynomial layout of the input data.
    !!
    !! Has to be either [[ply_dof_module:Q_Space]]
    !! or [[ply_dof_module:P_Space]].
    integer, intent(in) :: inspace

    !> Maximal polynomial degree in the input data.
    integer, intent(in) :: indegree

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Multi-dimensional polynomial layout of the output data.
    !!
    !! Has to be either [[ply_dof_module:Q_Space]]
    !! or [[ply_dof_module:P_Space]].
    integer, intent(in) :: outspace

    !> Maximal polynomial degree in the output data.
    integer, intent(in) :: outdegree

    !> Number of dimensions in the polynomials to transfer
    integer, intent(in) :: ndims
    ! -------------------------------------------------------------------- !

    select case(ndims)
    case(1)
      call ply_transfer_dofs_1D( indat     = indat,    &
        &                        indegree  = indegree, &
        &                        outdat    = outdat,   &
        &                        outdegree = outdegree )

    case(2)
      call ply_transfer_dofs_2D( indat     = indat,    &
        &                        inspace   = inspace,  &
        &                        indegree  = indegree, &
        &                        outdat    = outdat,   &
        &                        outspace  = outspace, &
        &                        outdegree = outdegree )

    case(3)
      call ply_transfer_dofs_3D( indat     = indat,    &
        &                        inspace   = inspace,  &
        &                        indegree  = indegree, &
        &                        outdat    = outdat,   &
        &                        outspace  = outspace, &
        &                        outdegree = outdegree )

    end select

  end subroutine ply_transfer_dofs
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transfer the polynomial in P representation from on dimension to
  !! another one.
  !!
  !! Only needed for P polynomials, for Q this trivially copying of
  !! contiguous memory.
  subroutine ply_transfer_P_dim( indat, indim, outdat, outdim, degree )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Dimension of the input polynomial.
    integer, intent(in) :: indim

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Dimension of the output polynomial.
    integer, intent(in) :: outdim

    !> Maximal polynomial degree in the output data.
    integer, intent(in) :: degree
    ! -------------------------------------------------------------------- !
    integer :: nInDofs, nOutDofs
    integer :: iMode
    integer :: iPos
    integer :: iX, iY, iZ
    ! -------------------------------------------------------------------- !

    in_d: select case(indim)
    case(1) in_d

      i1_out_d: select case(outdim)
      case(1) i1_out_d
        outdat = indat

      case(2) i1_out_d
        outdat = 0.0_rk
        do iMode=1,degree+1
  ! integer divisions are no mistake here.
  ipos = ((((imode - 1) + (1 - 1))            &
    &   * (((imode - 1) + (1 - 1)) + 1)) / 2 + 1) &
    & + (1 - 1)
          outdat(iPos) = indat(iMode)
        end do

      case(3) i1_out_d
        outdat = 0.0_rk
        do iMode=1,degree+1
  ! integer divisions are no mistake here.
  ipos = (((imode + 1 + 1 - 3) &
    &     * (imode + 1 + 1 - 2) &
    &     * (imode + 1 + 1 - 1)) &
    &   / 6 + 1)             &
    & + ((1-1) * (imode + 1 + 1 -2) &
    &   - ((1-2) * (1-1)) / 2) &
    & + (1-1)
          outdat(iPos) = indat(iMode)
        end do

      end select i1_out_d


    case(2) in_d

      i2_out_d: select case(outdim)
      case(1) i2_out_d
        do iMode=1,degree+1
  ! integer divisions are no mistake here.
  ipos = ((((imode - 1) + (1 - 1))            &
    &   * (((imode - 1) + (1 - 1)) + 1)) / 2 + 1) &
    & + (1 - 1)
          outdat(iMode) = indat(iPos)
        end do

      case(2) i2_out_d
        outdat = indat

      case(3) i2_out_d
        outdat = 0.0_rk
        nInDofs = ((degree+1)*(degree+2))/2
        iX = 1
        iY = 1
        iZ = 1
        do iMode=1,nInDofs
  ! integer divisions are no mistake here.
  ipos = (((ix + iy + iz - 3) &
    &     * (ix + iy + iz - 2) &
    &     * (ix + iy + iz - 1)) &
    &   / 6 + 1)             &
    & + ((iz-1) * (ix + iy + iz -2) &
    &   - ((iz-2) * (iz-1)) / 2) &
    & + (iy-1)
          outdat(iPos) = indat(iMode)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (ix .ne. 1) then
    ! next item
    ix = ix - 1
    iy = iy + 1
  else
    ! next layer
    ix = iy + 1
    iy = 1
  end if
        end do

      end select i2_out_d


    case(3) in_d

      i3_out_d: select case(outdim)
      case(1) i3_out_d
        do iMode=1,degree+1
  ! integer divisions are no mistake here.
  ipos = (((imode + 1 + 1 - 3) &
    &     * (imode + 1 + 1 - 2) &
    &     * (imode + 1 + 1 - 1)) &
    &   / 6 + 1)             &
    & + ((1-1) * (imode + 1 + 1 -2) &
    &   - ((1-2) * (1-1)) / 2) &
    & + (1-1)
          outdat(iMode) = indat(iPos)
        end do

      case(2) i3_out_d
        nOutDofs = ((degree+1)*(degree+2))/2
        iX = 1
        iY = 1
        iZ = 1
        do iMode=1,nOutDofs
  ! integer divisions are no mistake here.
  ipos = (((ix + iy + iz - 3) &
    &     * (ix + iy + iz - 2) &
    &     * (ix + iy + iz - 1)) &
    &   / 6 + 1)             &
    & + ((iz-1) * (ix + iy + iz -2) &
    &   - ((iz-2) * (iz-1)) / 2) &
    & + (iy-1)
          outdat(iMode) = indat(iPos)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (ix .ne. 1) then
    ! next item
    ix = ix - 1
    iy = iy + 1
  else
    ! next layer
    ix = iy + 1
    iy = 1
  end if
        end do

      case(3) i3_out_d
        outdat = indat

      end select i3_out_d


    end select in_d

  end subroutine ply_transfer_P_dim
  ! ************************************************************************ !


end module ply_transfer_module
