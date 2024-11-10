! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014, 2017-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Nikhil Anand, Harald Klimach,
! Verena Krupp, Peter Vitt, and Tobias Girresser for University of Siegen.
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
!> This module provides functions to transfer polynomials from and to the
!! oversampled representation for nodal treatments.
module ply_oversample_module

  use env_module,               only: rk
  use ply_dof_module,           only: Q_space, P_space
  use ply_poly_project_module,  only: ply_poly_project_type

  implicit none

  private

  public :: ply_convert2oversample
  public :: ply_convertFromOversample

contains


  ! ************************************************************************ !
  !> Copy a single element state into a larger array and pad it with zeroes.
  !!
  !! The additional modes might be necessary for proper dealing with
  !! nonlinear operations.
  !! Please note, that the oversampled representation is always a Q-Space
  !! polynomial.
  subroutine ply_convert2oversample( state, poly_proj, nDim, modalCoeffs, &
    &                                nScalars, ensure_positivity          )
    ! -------------------------------------------------------------------- !
    !> Description of the projection method.
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> State in a single element, to be oversampled.
    real(kind=rk), intent(in) :: state(:,:)

    !> The number of dimensions to determine the correct oversampling routine.
    integer, intent(in) :: nDim

    !> Oversampled array for modal coefficients from state.
    !!
    !! These are always Q-Polynomial representations, as projections only work
    !! with those.
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars

    !> Only use modes up to the point, where we are sure that the resulting
    !! polynomial will be positive everywhere.
    !! This is an array of logicals for each variable, if not given, the
    !! default is false (no positivity ensured).
    logical, intent(in), optional :: ensure_positivity(:)
    ! -------------------------------------------------------------------- !
    select case(nDim)
    case(1)
      call ply_convert2oversample_1d(state, poly_proj, modalCoeffs, nScalars, &
        &                            ensure_positivity                        )
    case(2)
      call ply_convert2oversample_2d(state, poly_proj, modalCoeffs, nScalars, &
        &                            ensure_positivity                        )
    case(3)
      call ply_convert2oversample_3d(state, poly_proj, modalCoeffs, &
        &                            ensure_positivity              )
    end select
  end subroutine ply_convert2oversample
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Truncating an oversampled polynomial representation back to the original
  !! representation.
  !!
  !! Note, that the oversampled array, which is given as an input here is
  !! is always a Q-Polynomial representation, while the target truncated state
  !! might also be a P-Polynomial.
  subroutine ply_convertFromOversample( modalCoeffs, poly_proj, nDim, state, &
    &                                   nScalars                             )
    ! -------------------------------------------------------------------- !
    !> Data of the projection method
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> Oversampled modal array for one element
    real(kind=rk), intent(in) :: modalCoeffs(:,:)

    !> The number of dimensions to determine the correct oversampling routine.
    integer, intent(in) :: nDim

    !> Truncated state for one element obtained from the modalCoeffs
    real(kind=rk), intent(out) :: state(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    ! -------------------------------------------------------------------- !
    select case(nDim)
    case(1)
      call ply_convertFromOversample_1d(modalCoeffs, poly_proj, state, nScalars)
    case(2)
      call ply_convertFromOversample_2d(modalCoeffs, poly_proj, state, nScalars)
    case(3)
      call ply_convertFromOversample_3d(modalCoeffs, poly_proj, state)
    end select
  end subroutine ply_convertFromOversample
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Copy a single element state into a larger array and pad it with zeroes.
  !!
  !! The additional modes might be necessary for proper dealing with
  !! nonlinear operations.
  !! Please note, that the oversampled representation is always a Q-Space
  !! polynomial.
  subroutine ply_convert2oversample_3d( state, poly_proj, modalCoeffs, &
    &                                   ensure_positivity              )
    ! -------------------------------------------------------------------- !
    !> Description of the projection method.
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> State in a single element, to be oversampled.
    real(kind=rk), intent(in) :: state(:,:)

    !> Oversampled array for modal coefficients from state.
    !!
    !! These are always Q-Polynomial representations, as projections only work
    !! with those.
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)

    !> Only use modes up to the point, where we are sure that the resulting
    !! polynomial will be positive everywhere.
    !!
    !! This is an array of logicals for each variable, if not given, the
    !! default is false (no positivity ensured).
    logical, intent(in), optional :: ensure_positivity(:)
    ! -------------------------------------------------------------------- !
    integer :: oversamp_degree
    integer :: nScalars
    integer :: iVar
    integer :: mpd1, mpd1_square, mpd1_cube
    integer :: iDegX, iDegY, iDegZ, idof, dof, dofOverSamp
    real(kind=rk) :: ordersum(3*(poly_proj%min_degree + 1)-2)
    real(kind=rk) :: varsum
    integer :: iOrd
    integer :: maxorders
    integer :: ord_lim
    ! -------------------------------------------------------------------- !
    ! Information for the oversampling loop
    oversamp_degree = poly_proj%oversamp_degree
    mpd1 = poly_proj%min_degree + 1
    mpd1_square = mpd1**2
    mpd1_cube = mpd1**3
    nScalars = size(modalCoeffs,2)
    maxorders = 3*(poly_proj%min_degree + 1)-2

    if (poly_proj%basisType == Q_Space) then
      posQ: if (present(ensure_positivity)) then
        ord_lim = maxorders
        modalCoeffs = 0.0_rk
        varQ: do iVar=1,nScalars
          if (ensure_positivity(iVar)) then
            ordersum = 0.0_rk
            do dof = 1, mpd1_cube
              iDegZ = (dof-1)/mpd1_square + 1
              iDegY = (dof-1-(iDegZ-1)*mpd1_square)/mpd1+1
              iDegX = mod(dof-1,mpd1)+1
              iOrd = iDegX+iDegY+iDegZ-2
              ordersum(iOrd) = ordersum(iOrd) + abs(state(dof,iVar))
            end do
            varsum = 0.0_rk
            do iOrd=2,ord_lim
              varsum = varsum + ordersum(iOrd)
              if (varsum >= ordersum(1)) then
                 ord_lim = iOrd-1
                 EXIT
              end if
            end do
          end if
        end do varQ
        do iVar=1,nScalars
          do dof = 1, mpd1_cube
            iDegZ = (dof-1)/mpd1_square + 1
            iDegY = (dof-1-(iDegZ-1)*mpd1_square)/mpd1+1
            iDegX = mod(dof-1,mpd1)+1
            iOrd = iDegX+iDegY+iDegZ-2
            if (iOrd <= ord_lim) then
              dofOverSamp = iDegX + ( iDegY-1  &
                &                     + (iDegZ-1)*(oversamp_degree+1) &
                &                   ) * (oversamp_degree+1)
              modalCoeffs(dofOverSamp,iVar) = state(dof,iVar)
            end if
          end do
        end do
      else posQ
        if (oversamp_degree == poly_proj%min_degree) then
          modalCoeffs = state
        else
          modalCoeffs = 0.0_rk
          do dof = 1, mpd1_cube
            iDegZ = (dof-1)/mpd1_square + 1
            iDegY = (dof-1-(iDegZ-1)*mpd1_square)/mpd1+1
            iDegX = mod(dof-1,mpd1)+1
            dofOverSamp = iDegX + ( iDegY-1  &
              &                     + (iDegZ-1)*(oversamp_degree+1) &
              &                   ) * (oversamp_degree+1)
            do iVar=1,nScalars
              modalCoeffs(dofOverSamp,iVar) = state(dof,iVar)
            end do
          end do
        end if
      end if posQ

    else !P_Space
      modalCoeffs = 0.0_rk
      posP: if (present(ensure_positivity)) then
        ord_lim = maxorders
        varP: do iVar=1,nScalars
          if (ensure_positivity(iVar)) then
            iDegX = 1
            iDegY = 1
            iDegZ = 1
            ordersum = 0.0_rk
            do idof = 1, poly_proj%body_3d%min_dofs
  ! integer divisions are no mistake here.
  dof = (((idegx + idegy + idegz - 3) &
    &     * (idegx + idegy + idegz - 2) &
    &     * (idegx + idegy + idegz - 1)) &
    &   / 6 + 1)             &
    & + ((idegz-1) * (idegx + idegy + idegz -2) &
    &   - ((idegz-2) * (idegz-1)) / 2) &
    & + (idegy-1)
              iOrd = iDegX+iDegY+iDegZ-2
              ordersum(iOrd) = ordersum(iOrd) + abs(state(dof,iVar))
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (idegx .ne. 1) then
    ! next item
    idegx = idegx - 1
    idegy = idegy + 1
  elseif (idegy .ne. 1) then
    ! next block
    idegx = idegy - 1
    idegy = 1
    idegz = idegz + 1
  else
    ! next layer
    idegx = idegz + 1
    idegy = 1
    idegz = 1
  end if
            end do
            varsum = 0.0_rk
            do iOrd=2,ord_lim
              varsum = varsum + ordersum(iOrd)
              if (varsum >= ordersum(1)) then
                 ord_lim = iOrd-1
                 EXIT
              end if
            end do
          end if
        end do varP
        do iVar=1,nScalars
          iDegX = 1
          iDegY = 1
          iDegZ = 1
          do idof = 1, poly_proj%body_3d%min_dofs
            iOrd = iDegX+iDegY+iDegZ-2
            if (iOrd > ord_lim) EXIT
  ! integer divisions are no mistake here.
  dof = (((idegx + idegy + idegz - 3) &
    &     * (idegx + idegy + idegz - 2) &
    &     * (idegx + idegy + idegz - 1)) &
    &   / 6 + 1)             &
    & + ((idegz-1) * (idegx + idegy + idegz -2) &
    &   - ((idegz-2) * (idegz-1)) / 2) &
    & + (idegy-1)
            dofOverSamp = iDegX + ( iDegY-1  &
              &                     + (iDegZ-1)*(oversamp_degree+1) &
              &                   ) * (oversamp_degree+1)
            modalCoeffs(dofOverSamp,iVar) = state(dof,iVar)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (idegx .ne. 1) then
    ! next item
    idegx = idegx - 1
    idegy = idegy + 1
  elseif (idegy .ne. 1) then
    ! next block
    idegx = idegy - 1
    idegy = 1
    idegz = idegz + 1
  else
    ! next layer
    idegx = idegz + 1
    idegy = 1
    idegz = 1
  end if
          end do
        end do
      else posP
        iDegX = 1
        iDegY = 1
        iDegZ = 1
        do idof = 1, poly_proj%body_3d%min_dofs
  ! integer divisions are no mistake here.
  dof = (((idegx + idegy + idegz - 3) &
    &     * (idegx + idegy + idegz - 2) &
    &     * (idegx + idegy + idegz - 1)) &
    &   / 6 + 1)             &
    & + ((idegz-1) * (idegx + idegy + idegz -2) &
    &   - ((idegz-2) * (idegz-1)) / 2) &
    & + (idegy-1)
          dofOverSamp = iDegX + ( iDegY-1  &
            &                     + (iDegZ-1)*(oversamp_degree+1) &
            &                   ) * (oversamp_degree+1)
          modalCoeffs(dofOverSamp,:) = state(dof,:)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (idegx .ne. 1) then
    ! next item
    idegx = idegx - 1
    idegy = idegy + 1
  elseif (idegy .ne. 1) then
    ! next block
    idegx = idegy - 1
    idegy = 1
    idegz = idegz + 1
  else
    ! next layer
    idegx = idegz + 1
    idegy = 1
    idegz = 1
  end if
        end do
      end if posP

    end if

  end subroutine ply_convert2oversample_3d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Truncating an oversampled polynomial representation back to the original
  !! representation.
  !!
  !! Note, that the oversampled array, which is given as an input here is
  !! is always a Q-Polynomial representation, while the target truncated state
  !! might also be a P-Polynomial.
  subroutine ply_convertFromOversample_3d( modalCoeffs, poly_proj, state )
    ! -------------------------------------------------------------------- !
    !> Data of the projection method
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> Oversampled modal array for one element
    real(kind=rk), intent(in) :: modalCoeffs(:,:)

    !> Truncated state for one element obtained from the modalCoeffs
    real(kind=rk), intent(out) :: state(:,:)
    ! -------------------------------------------------------------------- !
    integer :: oversamp_degree
    integer :: nScalars
    integer :: iVar
    integer :: mpd1, mpd1_square, mpd1_cube
    integer :: iDegX, iDegY, iDegZ, idof, dof, dofOverSamp
    ! -------------------------------------------------------------------- !

    ! Information for the oversampling loop
    oversamp_degree = poly_proj%oversamp_degree
    mpd1 = poly_proj%min_degree + 1
    mpd1_square = mpd1**2
    mpd1_cube = mpd1**3
    nScalars = size(modalCoeffs,2)

    if (poly_proj%basisType == Q_Space) then
      if (oversamp_degree == poly_proj%min_degree) then
        state = modalCoeffs
      else
        do iVar=1,nScalars
          do dof = 1, mpd1_cube
            iDegZ = (dof-1)/mpd1_square + 1
            iDegY = (dof-1-(iDegZ-1)*mpd1_square)/mpd1+1
            iDegX = mod(dof-1,mpd1)+1
            dofOverSamp = iDegX + ( iDegY-1  &
              &                     + (iDegZ-1)*(oversamp_degree+1) &
              &                   ) * (oversamp_degree+1)
            state(dof,iVar) = modalCoeffs(dofOverSamp,iVar)
          end do
        end do
      end if

    else !P_Space

      iDegX = 1
      iDegY = 1
      iDegZ = 1
      do idof = 1, poly_proj%body_3d%min_dofs
  ! integer divisions are no mistake here.
  dof = (((idegx + idegy + idegz - 3) &
    &     * (idegx + idegy + idegz - 2) &
    &     * (idegx + idegy + idegz - 1)) &
    &   / 6 + 1)             &
    & + ((idegz-1) * (idegx + idegy + idegz -2) &
    &   - ((idegz-2) * (idegz-1)) / 2) &
    & + (idegy-1)
        dofOverSamp = iDegX + ( iDegY-1  &
          &                     + (iDegZ-1)*(oversamp_degree+1) &
          &                   ) * (oversamp_degree+1)
        state(dof,:) = modalCoeffs(dofOverSamp,:)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (idegx .ne. 1) then
    ! next item
    idegx = idegx - 1
    idegy = idegy + 1
  elseif (idegy .ne. 1) then
    ! next block
    idegx = idegy - 1
    idegy = 1
    idegz = idegz + 1
  else
    ! next layer
    idegx = idegz + 1
    idegy = 1
    idegz = 1
  end if
      end do
    end if

  end subroutine ply_convertFromoversample_3d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Copy a single 2D element state into a larger array and pad it with zeroes.
  !!
  !! The additional modes might be necessary for proper dealing with
  !! nonlinear operations.
  !! Please note, that the oversampled representation is always a Q-Space
  !! polynomial.
  subroutine ply_convert2oversample_2d( state, poly_proj, modalCoeffs, &
    &                                   nScalars, ensure_positivity    )
    ! -------------------------------------------------------------------- !
    !> Description of the projection method.
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> State in a single element, to be oversampled.
    real(kind=rk), intent(in) :: state(:,:)

    !> Oversampled array for modal coefficients from state.
    !!
    !! These are always Q-Polynomial representations, as projections only work
    !! with those.
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars

    !> Only use modes up to the point, where we are sure that the resulting
    !! polynomial will be positive everywhere.
    !! This is an array of logicals for each variable, if not given, the
    !! default is false (no positivity ensured).
    logical, intent(in), optional :: ensure_positivity(:)
    ! -------------------------------------------------------------------- !
    integer :: oversamp_degree
    integer :: mpd1, mpd1_square
    integer :: iDegX, iDegY, idof, dof, dofOverSamp, nPVars
    integer :: iVar
    real(kind=rk) :: ordersum(2*(poly_proj%min_degree + 1)-1)
    real(kind=rk) :: varsum
    integer :: iOrd
    integer :: maxorders
    integer :: ord_lim
    ! -------------------------------------------------------------------- !
    ! Information for the oversampling loop
    if(present(nScalars)) then
      nPVars = nScalars
    else
      nPVars = size(state,2)
    end if
    maxorders = 2*(poly_proj%min_degree + 1)-1

    ! Information for the oversampling loop
    oversamp_degree = poly_proj%oversamp_degree
    mpd1 = poly_proj%min_degree + 1
    mpd1_square = mpd1**2

    ! Initialize oversampled space correct to 0
    modalCoeffs(:,:) = 0.0_rk

    if (poly_proj%basisType == Q_Space) then
      posQ: if (present(ensure_positivity)) then
        ord_lim = maxorders
        varQ: do iVar=1,nPVars
          if (ensure_positivity(iVar)) then
            ordersum = 0.0_rk
            do dof = 1, mpd1_square
              iDegX = mod(dof-1,mpd1)+1
              iDegY = (dof-1)/mpd1+1
              iOrd = iDegX+iDegY-1
              ordersum(iOrd) = ordersum(iOrd) + abs(state(dof,iVar))
            end do
            varsum = 0.0_rk
            do iOrd=2,ord_lim
              varsum = varsum + ordersum(iOrd)
              if (varsum >= ordersum(1)) then
                 ord_lim = iOrd-1
                 EXIT
              end if
            end do
          end if
        end do varQ
        do iVar=1,nPVars
          do dof = 1, mpd1_square
            iDegX = mod(dof-1,mpd1)+1
            iDegY = (dof-1)/mpd1+1
            iOrd = iDegX+iDegY-1
            if (iOrd <= ord_lim) then
              dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(oversamp_degree+1)
              modalCoeffs(dofOverSamp,iVar) = state(dof,iVar)
            end if
          end do
        end do
      else posQ
        do dof = 1, mpd1_square
          iDegX = mod(dof-1,mpd1)+1
          iDegY = (dof-1)/mpd1+1
          dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(oversamp_degree+1)
          modalCoeffs(dofOverSamp,1:nPVars) = state(dof,1:nPVars)
        end do
      end if posQ

    else !P_Space

      posP: if (present(ensure_positivity)) then
        ord_lim = maxorders
        varP: do iVar=1,nPVars
          if (ensure_positivity(iVar)) then
            iDegX = 1
            iDegY = 1
            ordersum = 0.0_rk
            do idof = 1, poly_proj%body_2d%min_dofs
  ! integer divisions are no mistake here.
  dof = ((((idegx - 1) + (idegy - 1))            &
    &   * (((idegx - 1) + (idegy - 1)) + 1)) / 2 + 1) &
    & + (idegy - 1)
              iOrd = iDegX+iDegY-1
              ordersum(iOrd) = ordersum(iOrd) + abs(state(dof,iVar))
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (idegx .ne. 1) then
    ! next item
    idegx = idegx - 1
    idegy = idegy + 1
  else
    ! next layer
    idegx = idegy + 1
    idegy = 1
  end if
            end do
            varsum = 0.0_rk
            do iOrd=2,ord_lim
              varsum = varsum + ordersum(iOrd)
              if (varsum >= ordersum(1)) then
                 ord_lim = iOrd-1
                 EXIT
              end if
            end do
          end if
        end do varP
        do iVar=1,nPVars
          iDegX = 1
          iDegY = 1
          do idof = 1, poly_proj%body_2d%min_dofs
            iOrd = iDegX+iDegY-1
            if (iOrd > ord_lim) EXIT
  ! integer divisions are no mistake here.
  dof = ((((idegx - 1) + (idegy - 1))            &
    &   * (((idegx - 1) + (idegy - 1)) + 1)) / 2 + 1) &
    & + (idegy - 1)
            dofOverSamp = iDegX + (iDegY-1)*(oversamp_degree+1)
            modalCoeffs(dofOverSamp,1:nPVars) = state(dof,1:nPVars)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (idegx .ne. 1) then
    ! next item
    idegx = idegx - 1
    idegy = idegy + 1
  else
    ! next layer
    idegx = idegy + 1
    idegy = 1
  end if
          end do
        end do
      else posP
        iDegX = 1
        iDegY = 1
        do idof = 1, poly_proj%body_2d%min_dofs
  ! integer divisions are no mistake here.
  dof = ((((idegx - 1) + (idegy - 1))            &
    &   * (((idegx - 1) + (idegy - 1)) + 1)) / 2 + 1) &
    & + (idegy - 1)
          dofOverSamp = iDegX + (iDegY-1)*(oversamp_degree+1)
          modalCoeffs(dofOverSamp,1:nPVars) = state(dof,1:nPVars)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (idegx .ne. 1) then
    ! next item
    idegx = idegx - 1
    idegy = idegy + 1
  else
    ! next layer
    idegx = idegy + 1
    idegy = 1
  end if
        end do
      end if posP

    end if

  end subroutine ply_convert2oversample_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Truncating an oversampled 2D polynomial representation back to the original
  !! representation.
  !!
  !! Note, that the oversampled array, which is given as an input here is
  !! is always a Q-Polynomial representation, while the target truncated state
  !! might also be a P-Polynomial.
  subroutine ply_convertFromOversample_2d( modalCoeffs, poly_proj, state, &
    &                                      nScalars                       )
    ! -------------------------------------------------------------------- !
    !> Data of the projection method
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> Oversampled modal array for one element
    real(kind=rk), intent(in) :: modalCoeffs(:,:)

    !> Truncated state for one element obtained from the modalCoeffs
    real(kind=rk), intent(out) :: state(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    ! -------------------------------------------------------------------- !
    integer :: oversamp_degree
    integer :: mpd1, mpd1_square
    integer :: iDegX, iDegY, iDegZ, idof, dof, dofOverSamp, nPVars
    ! -------------------------------------------------------------------- !
    ! Information for the oversampling loop
    if(present(nScalars)) then
      nPVars = nScalars
    else
      nPVars = size(state,2)
    end if

    ! Information for the oversampling loop
    oversamp_degree = poly_proj%oversamp_degree
    mpd1 = poly_proj%min_degree + 1
    mpd1_square = mpd1**2

    if (poly_proj%basisType == Q_Space) then

      do dof = 1, mpd1_square
        iDegX = mod(dof-1,mpd1)+1
        iDegY = (dof-1)/mpd1+1
        dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(oversamp_degree+1)
        state(dof,1:nPVars) = modalCoeffs(dofOverSamp,1:nPVars)
      end do

    else !P_Space

      iDegX = 1
      iDegY = 1
      iDegZ = 0 ! not used in posOfModgCoeffPTens_2D, nextModgCoeffPTens
      do idof = 1, poly_proj%body_2d%min_dofs
  ! integer divisions are no mistake here.
  dof = ((((idegx - 1) + (idegy - 1))            &
    &   * (((idegx - 1) + (idegy - 1)) + 1)) / 2 + 1) &
    & + (idegy - 1)
        dofOverSamp = iDegX + (iDegY-1)*(oversamp_degree+1)
        state(dof,1:nPVars) = modalCoeffs(dofOverSamp,1:nPVars)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (idegx .ne. 1) then
    ! next item
    idegx = idegx - 1
    idegy = idegy + 1
  else
    ! next layer
    idegx = idegy + 1
    idegy = 1
  end if
      end do

    end if

  end subroutine ply_convertFromoversample_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Copy a single 1D element state into a larger array and pad it with zeroes.
  !!
  !! The additional modes might be necessary for proper dealing with
  !! nonlinear operations.
  subroutine ply_convert2oversample_1d( state, poly_proj, modalCoeffs, &
    &                                   nScalars,  ensure_positivity   )
    ! -------------------------------------------------------------------- !
    !> Description of the projection method.
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> State in a single element, to be oversampled.
    real(kind=rk), intent(in) :: state(:,:)

    !> Oversampled array for modal coefficients from state.
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars

    !> Only use modes up to the point, where we are sure that the resulting
    !! polynomial will be positive everywhere.
    !! This is an array of logicals for each variable, if not given, the
    !! default is false (no positivity ensured).
    logical, intent(in), optional :: ensure_positivity(:)
    ! -------------------------------------------------------------------- !
    integer :: iVar, iPoint, iVP, nPVars
    integer :: nVars
    integer :: ord_lim
    real(kind=rk) :: varsum
    ! -------------------------------------------------------------------- !
    ! Information for the oversampling loop
    if(present(nScalars)) then
      nVars = nScalars
    else
      nVars = size(state,2)
    end if
    nPVars = (poly_proj%maxPolyDegree+1)*nVars

    ! Initialize oversampled space correct to 0
    ModalCoeffs(:,:) = 0.0_rk

    if (present(ensure_positivity)) then
      ord_lim = poly_proj%min_degree+1
      do iVar = 1,nVars
        if (ensure_positivity(iVar)) then
          varSum = 0.0_rk
          do iPoint=2,ord_lim
            varSum = varSum + abs(state(iPoint,iVar))
            if (varSum >= state(1,iVar)) then
              ord_lim = iPoint - 1
              EXIT
            end if
          end do
        end if
      end do
      do iVar=1,nVars
        do iPoint=1,ord_lim
          ModalCoeffs(iPoint,iVar) = state(iPoint,iVar)
        end do
      end do
    else
      do iVP = 1,nPVars
        iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
        iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)
        ModalCoeffs(iPoint,iVar) = state(iPoint,iVar)
      end do
    end if

  end subroutine ply_convert2oversample_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Truncating an oversampled 1D polynomial representation back to the original
  !! representation.
  subroutine ply_convertFromOversample_1d( modalCoeffs, poly_proj, state, &
    &                                      nScalars                       )
    ! -------------------------------------------------------------------- !
    !> Data of the projection method
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> Oversampled modal array for one element
    real(kind=rk), intent(in) :: modalCoeffs(:,:)

    !> Truncated state for one element obtained from the modalCoeffs
    real(kind=rk), intent(out) :: state(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    ! -------------------------------------------------------------------- !
    integer :: iVar, iPoint, iVP, nPVars
    ! -------------------------------------------------------------------- !
    ! Information for the oversampling loop
    if(present(nScalars)) then
      nPVars = (poly_proj%maxPolyDegree+1)*nScalars
    else
      nPVars = (poly_proj%maxPolyDegree+1)*size(state,2)
    end if

    do iVP = 1,nPVars
      iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
      iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)
      state(iPoint,iVar) = modalCoeffs(iPoint,iVar)
    end do

  end subroutine ply_convertFromoversample_1d
  ! ************************************************************************ !

end module ply_oversample_module
