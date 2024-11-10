! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014,2016,2018,2020,2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2013-2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Simon Zimny and Harald Klimach
! for German Research School for Simulation Sciences GmbH.
!
! Parts of this file were written by Verena Krupp, Harald Klimach, Peter Vitt
! and Kay Langhammer for University of Siegen.
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
!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_2D_module
  use, intrinsic :: iso_c_binding


  use fftw_wrap

  use env_module,         only: rk
  use ply_legFpt_module,  only: ply_legFpt_type, &
    &                           assignment(=)

  implicit none

  private

  interface ply_legToPnt_2D
    module procedure ply_legToPnt_2D_singVar
    module procedure ply_legToPnt_2D_multVar
  end interface ply_legToPnt_2D

  interface ply_pntToLeg_2D
    module procedure ply_pntToLeg_2D_singVar
    module procedure ply_pntToLeg_2D_multVar
  end interface ply_pntToLeg_2D

  public :: ply_legToPnt_2D, ply_pntToLeg_2D


contains


  ! ------------------------------------------------------------------------ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt_2D_singVar( fpt, legCoeffs, pntVal )
    ! -------------------------------------------------------------------- !
    !> The FPT parameters.
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The Legendre coefficients to convert to point values (Chebyshev nodes).
    !! \attention Although this array serves as input only, it is modified
    !! inside of this routine by the underlying FPT algorithm. So, when
    !! this routine returns from its call the original values of pntVal will
    !! be modified.
    real(kind=rk), intent(inout) :: legCoeffs(:)
    !> The resulting point values (Chebyshev nodes).
    real(kind=rk), intent(inout) :: pntVal(:)
    ! -------------------------------------------------------------------- !
    integer :: n, iAlph, nIndeps
    real(kind=rk), dimension(:), allocatable :: alph
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n
    nIndeps = n

    allocate(alph(n**2))

    ! original layout (n = 3):
    !  1  2  3
    !  4  5  6
    !  7  8  9

    ! layout after y-trafo:
    !  1  4  7
    !  2  5  8
    !  3  6  9
    ! >>>>> Y-Direction >>>>> !
    ! iAlph is the index of the first element in a line for the transformation
    ! in y-direction.
    do iAlph = 1, n
      alph((iAlph-1)*n+1:iAlph*n) = legCoeffs(iAlph::n)
    end do

    call fpt%legToPnt( nIndeps   = nIndeps,  &
      &                legCoeffs = alph,     &
      &                pntVal    = legCoeffs )
    ! <<<<< Y-Direction <<<<< !

    ! >>>>> X-Direction >>>>> !
    do iAlph = 1, n
      alph((iAlph-1)*n+1:iAlph*n) = legCoeffs(iAlph::n)
    end do

    call fpt%legToPnt( nIndeps   = nIndeps, &
      &                legCoeffs = alph,    &
      &                pntVal    = pntVal   )
    ! <<<<< X-Direction <<<<< !

  end subroutine ply_legToPnt_2D_singVar
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt_2D_multVar( fpt, legCoeffs, pntVal, nVars )
    ! -------------------------------------------------------------------- !
    !> The FPT parameters.
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The Legendre coefficients to convert to point values (Chebyshev nodes).
    !! \attention Although this array serves as input only, it is modified
    !! inside of this routine by the underlying FPT algorithm. So, when
    !! this routine returns from its call the original values of pntVal will
    !! be modified.
    real(kind=rk), intent(inout) :: legCoeffs(:,:)
    !> The resulting point values (Chebyshev nodes).
    real(kind=rk), intent(inout) :: pntVal(:,:)
    !> The number of scalar variables to transform.
    integer, intent(in) :: nVars
    ! -------------------------------------------------------------------- !
    integer :: iVar
    ! -------------------------------------------------------------------- !

    do iVar = 1, nVars
      call ply_legToPnt_2D(fpt, legCoeffs(:,iVar), pntVal(:,iVar))
    end do

  end subroutine ply_legToPnt_2D_multVar
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_2D_singVar( fpt, pntVal, legCoeffs )
    ! -------------------------------------------------------------------- !
    !> Parameters of the Fast Polynomial transformation.
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The point values to transform to 2D modal Legendre expansion.
    !! \attention Although this array serves as input only, it is modified
    !! inside of this routine by the underlying DCT algorithm. So, when
    !! this routine returns from its call the original values of pntVal will
    !! be modified.
    real(kind=rk), intent(inout) :: pntVal(:)
    !> The Legendre coefficients.
    real(kind=rk), intent(inout) :: legCoeffs(:)
    ! -------------------------------------------------------------------- !
    integer :: nIndeps, iAlph, n, n_squared
    real(kind=rk), dimension(:), allocatable :: alph
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n
    nIndeps = n
    n_squared = n**2

    allocate(alph(n_squared))

    ! >>>>> Y-Direction >>>>> !
    do iAlph = 1, n
      alph((iAlph-1)*n+1:iAlph*n) = pntVal(iAlph::n)
    end do

    call fpt%pntToLeg( nIndeps   = nIndeps, &
      &                legCoeffs = pntVal,  &
      &                pntVal    = alph     )
    ! <<<<< Y-Direction <<<<< !

    ! x-direction
    ! >>>>> X-Direction >>>>> !
    do iAlph = 1, n
      !ztrafo
      alph((iAlph-1)*n+1:iAlph*n) = pntVal(iAlph::n)
    end do

    call fpt%pntToLeg( nIndeps   = nIndeps,   &
      &                legCoeffs = legCoeffs, &
      &                pntVal    = alph       )
    ! <<<<< X-Direction <<<<< !

  end subroutine ply_pntToLeg_2D_singVar
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_2D_multVar( fpt, pntVal, legCoeffs, nVars )
    ! -------------------------------------------------------------------- !
    !> Parameters of the Fast Polynomial transformation.
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The point values to transform to 2D modal Legendre expansion.
    !! \attention Although this array serves as input only, it is modified
    !! inside of this routine by the underlying DCT algorithm. So, when
    !! this routine returns from its call the original values of pntVal will
    !! be modified.
    real(kind=rk), intent(inout) :: pntVal(:,:)
    !> The Legendre coefficients.
    real(kind=rk), intent(inout) :: legCoeffs(:,:)
    !> The number of scalar variables to transform.
    integer, intent(in) :: nVars
    ! -------------------------------------------------------------------- !
    integer :: iVar
    ! -------------------------------------------------------------------- !

    do iVar = 1, nVars
      call ply_pntToLeg_2D(fpt, pntVal(:,iVar), legCoeffs(:,iVar))
    end do

  end subroutine ply_pntToLeg_2D_multVar
  ! ------------------------------------------------------------------------ !

end module ply_legFpt_2D_module
