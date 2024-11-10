! Copyright (c) 2012-2014, 2016, 2018 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014,2020 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2016 Langhammer Kay <kay.langhammer@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop and Harald Klimach
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

!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_module
  use, intrinsic :: iso_c_binding
  use env_module,             only: rk
  use tem_compileconf_module, only: vlen
  use ply_polyBaseExc_module, only: ply_trafo_params_type, &
    &                               ply_fpt_init,          &
    &                               ply_fpt_exec_striped,  &
    &                               ply_fpt_single,        &
    &                               ply_legToCheb_param,   &
    &                               ply_chebToLeg_param,   &
    &                               assignment(=)
  use ply_fpt_header_module, only: ply_fpt_header_type, &
    &                              ply_fpt_scalar,      &
    &                              ply_fpt_vector
  use fftw_wrap

  implicit none

  private

  !> Datatype for parameters of the FPT used for 1d, 2d and 3d.
  !!
  !! Stores of the parameters for a fast conversion of a modal
  !! Legendre expansion to point values (located at Chebyshev nodes)
  !! and vice versa. \n
  !! The FPT covers: \n
  !! - Transformation from Legendre expansion to point values
  !!   at Chebyshev nodes \n
  !! - Transformation from point values (Chebyshev nodes) to
  !!   modal Legendre expansion \n
  type ply_legFpt_type
    !> FPT params for the fast base exchange from Legendre to
    !! Chebyshev expansion.
    type(ply_trafo_params_type) :: legToChebParams

    !> FPT params for the fast base exchange from Chebyshev to
    !! Legendre expansion.
    type(ply_trafo_params_type) :: chebToLegParams

    !> FFTW plan for DCT from Chebyshev coefficients to point values.
    type(C_PTR) :: planChebToPnt

    !> FFTW plan for DCT from point values to Chebyshev coefficients.
    type(C_PTR) :: planPntToCheb

    !> Flag whether to use Lobatto points (include boundary points)
    logical :: use_lobatto_points

    procedure(ply_fptm2n), pointer :: legtopnt => NULL()
    procedure(ply_fptn2m), pointer :: pnttoleg => NULL()

  end type ply_legFpt_type

  interface assignment(=)
    module procedure Copy_fpt
  end interface

  public :: ply_legFpt_type, ply_init_legFpt
  public :: assignment(=)


  interface
    subroutine ply_fptm2n( fpt, legCoeffs, pntVal, nIndeps )
      import :: ply_legFpt_type, rk
      real(kind=rk), intent(inout) :: legCoeffs(:)
      class(ply_legFpt_type), intent(inout) :: fpt
      real(kind=rk), intent(inout) :: pntVal(:)
      integer, intent(in) :: nIndeps
    end subroutine
    subroutine ply_fptn2m( fpt, pntVal, legCoeffs, nIndeps )
      import :: ply_legFpt_type, rk
      real(kind=rk), intent(inout) :: pntVal(:)
      class(ply_legFpt_type), intent(inout) :: fpt
      real(kind=rk), intent(inout) :: legCoeffs(:)
      integer, intent(in) :: nIndeps
    end subroutine
  end interface


contains


  ! ------------------------------------------------------------------------ !
  subroutine Copy_fpt( left, right )
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_legFpt_type), intent(out) :: left
    !> fpt to copy from
    type(ply_legFpt_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%legToChebParams = right%legToChebParams
    left%chebToLegParams = right%chebToLegParams

    left%planChebToPnt = right%planChebToPnt
    left%planPntToCheb = right%planPntToCheb

    left%legtopnt => right%legtopnt
    left%pnttoleg => right%pnttoleg

  end subroutine Copy_fpt
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to initialize the fast polynomial transformation
  !! for Legendre expansion.
  subroutine ply_init_legFpt( maxPolyDegree, nIndeps, fpt, header, fft_flags)
    ! -------------------------------------------------------------------- !
    !> Maximal polynomial degree for the transformation.
    integer, intent(in) :: maxPolyDegree

    !> Number of independent values that can be computed simultaneously.
    integer, intent(in) :: nIndeps

    !> The Fast Polynomial Transformation setting to initialize.
    type(ply_legFpt_type), intent(inout) :: fpt

    !> Configuration settings for the projection.
    type(ply_fpt_header_type), intent(in) :: header

    !> Planning flags for the FFT.
    !!
    !! Configuration to how much time to spend on finding an optimal FFT
    !! implementation in the FFTW.
    !! See: http://www.fftw.org/doc/Planner-Flags.html#Planner-Flags
    !!
    !! Defaults to FFTW_MEASURE.
    integer, optional, intent(in) :: fft_flags
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: tmpOut(:), tmpIn(:)
    logical :: lob
    integer :: n
    integer :: maxstriplen
    integer :: planning_flags
    ! -------------------------------------------------------------------- !

    if (present(fft_flags)) then
      planning_flags = fft_flags
    else
      planning_flags = FFTW_MEASURE
    end if

    maxstriplen = min(header%striplen, nIndeps)

    lob = header%nodes_header%lobattoPoints

    fpt%use_lobatto_points = lob

    ! Init the fast Legendre to Chebyshev transformation.
    call ply_fpt_init( n                = maxPolyDegree+1,        &
      &                params           = fpt%legToChebParams,    &
      &                trafo            = ply_legToCheb_param,    &
      &                blocksize        = header%blocksize,       &
      &                approx_terms     = header%approx_terms,    &
      &                striplen         = maxstriplen,            &
      &                subblockingWidth = header%subblockingWidth )

    ! Init the fast Chebyshev to Legendre transformation.
    call ply_fpt_init( n                = maxPolyDegree+1,        &
      &                params           = fpt%chebToLegParams,    &
      &                trafo            = ply_chebToLeg_param,    &
      &                blocksize        = header%blocksize,       &
      &                approx_terms     = header%approx_terms,    &
      &                striplen         = maxstriplen,            &
      &                subblockingWidth = header%subblockingWidth )

    ! Create the buffers for the intermediate arrays
    n = fpt%legToChebParams%n

    ! Temporary arrays to initialize FFTW real->real transformations
    select case(header%implementation)
    case(ply_fpt_scalar)
      allocate( tmpIn(n) )
      allocate( tmpOut(n) )
    case(ply_fpt_vector)
      allocate( tmpIn(n*nIndeps) )
      allocate( tmpOut(n*nIndeps) )
    end select

    if (.not. lob) then

      select case(header%implementation)
      case(ply_fpt_scalar)
        fpt%legtopnt => ply_legtopnt_single
        fpt%pnttoleg => ply_pnttoleg_single
        ! Init the DCT III ( Leg -> Point values )
        !NEC: The NEC FFTW interface use n0 as parameter name instead of n
        !NEC: (nlc 2.0.0).
        !NEC: Omitting keywords to be compatible.
        !!fpt%planPntToCheb = fftw_plan_r2r_1d( n     = n,             &
        !!  &                                   in    = tmpIn,         &
        !!  &                                   out   = tmpOut,        &
        !!  &                                   kind  = FFTW_REDFT10,  &
        !!  &                                   flags = planning_flags )
        fpt%planChebToPnt = fftw_plan_r2r_1d( n,             &
          &                                   tmpIn,         &
          &                                   tmpOut,        &
          &                                   FFTW_REDFT01,  &
          &                                   planning_flags )
        ! Init the DCT II ( Point values -> Leg )
        !NEC: The NEC FFTW interface use n0 as parameter name instead of n
        !NEC: (nlc 2.0.0).
        !NEC: Omitting keywords to be compatible.
        !!fpt%planPntToCheb = fftw_plan_r2r_1d( n     = n,             &
        !!  &                                   in    = tmpIn,         &
        !!  &                                   out   = tmpOut,        &
        !!  &                                   kind  = FFTW_REDFT10,  &
        !!  &                                   flags = planning_flags )
        fpt%planPntToCheb = fftw_plan_r2r_1d( n,             &
          &                                   tmpIn,         &
          &                                   tmpOut,        &
          &                                   FFTW_REDFT10,  &
          &                                   planning_flags )

      case(ply_fpt_vector)
        fpt%legtopnt => ply_legtopnt_vec
        fpt%pnttoleg => ply_pnttoleg_vec
        fpt%planChebToPnt = fftw_plan_many_r2r(         &
          &                   rank    = 1,              &
          &                   n       = [n],            &
          &                   howmany = nIndeps,        &
          &                   in      = tmpIn,          &
          &                   inembed = [n],            &
          &                   istride = 1,              &
          &                   idist   = n,              &
          &                   out     = tmpOut,         &
          &                   onembed = [n],            &
          &                   ostride = 1,              &
          &                   odist   = n,              &
          &                   kind    = [FFTW_REDFT01], &
          &                   flags   = planning_flags  )
        fpt%planPntToCheb = fftw_plan_many_r2r(         &
          &                   rank    = 1,              &
          &                   n       = [n],            &
          &                   howmany = nIndeps,        &
          &                   in      = tmpIn,          &
          &                   inembed = [n],            &
          &                   istride = 1,              &
          &                   idist   = n,              &
          &                   out     = tmpOut,         &
          &                   onembed = [n],            &
          &                   ostride = 1,              &
          &                   odist   = n,              &
          &                   kind    = [FFTW_REDFT10], &
          &                   flags   = planning_flags  )

      end select

    else

      select case(header%implementation)
      case(ply_fpt_scalar)
        fpt%legtopnt => ply_legtopnt_lobatto_single
        fpt%pnttoleg => ply_pnttoleg_lobatto_single
        ! Init the DCT I  (Leg -> nodal):
        !   To be used with a normalization factor for trafo ...
        !NEC: The NEC FFTW interface use n0 as parameter name instead of n
        !NEC: (nlc 2.0.0).
        !NEC: Omitting keywords to be compatible.
        !!fpt%planChebToPnt = fftw_plan_r2r_1d( n     = n,             &
        !!  &                                   in    = tmpIn,         &
        !!  &                                   out   = tmpOut,        &
        !!  &                                   kind  = FFTW_REDFT00,  &
        !!  &                                   flags = planning_flags )
        fpt%planChebToPnt = fftw_plan_r2r_1d( n,             &
          &                                   tmpIn,         &
          &                                   tmpOut,        &
          &                                   FFTW_REDFT00,  &
          &                                   planning_flags )

      case(ply_fpt_vector)
        fpt%legtopnt => ply_legtopnt_lobatto_vec
        fpt%pnttoleg => ply_pnttoleg_lobatto_vec
        fpt%planChebToPnt = fftw_plan_many_r2r(         &
          &                   rank    = 1,              &
          &                   n       = [n],            &
          &                   howmany = nIndeps,        &
          &                   in      = tmpIn,          &
          &                   inembed = [n],            &
          &                   istride = 1,              &
          &                   idist   = n,              &
          &                   out     = tmpOut,         &
          &                   onembed = [n],            &
          &                   ostride = 1,              &
          &                   odist   = n,              &
          &                   kind    = [FFTW_REDFT00], &
          &                   flags   = planning_flags  )

      end select

      ! Init the DCT I  (nodal -> Leg):
      !   To be used with a normalization factor for trafo ...
      fpt%planPntToCheb = fpt%planChebToPnt

    end if

    deallocate( tmpIn, tmpOut )

  end subroutine ply_init_legFpt
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt_single( fpt, legCoeffs, pntVal, nIndeps )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(inout) :: legCoeffs(:)
    class(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n)
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    do iDof = 1, nIndeps*n, n
      call ply_fpt_single( alph   = legCoeffs(iDof:iDof+n-1), &
        &                  gam    = cheb,                     &
        &                  params = fpt%legToChebParams       )

      ! Normalize the coefficients of the Chebyshev polynomials due
      ! to the unnormalized version of DCT in the FFTW.
      cheb(2:n:2) = -0.5_rk * cheb(2:n:2)
      cheb(3:n:2) =  0.5_rk * cheb(3:n:2)

      call fftw_execute_r2r( fpt%planChebToPnt,    &
        &                    cheb,                 &
        &                    pntVal(iDof:iDof+n-1) )
    end do

  end subroutine ply_legToPnt_single
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Vectorizing subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt_vec( fpt, legCoeffs, pntVal, nIndeps )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(inout) :: legCoeffs(:)
    class(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n*nIndeps)
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    call ply_fpt_exec_striped( nIndeps = nIndeps,            &
      &                        alph    = legCoeffs,          &
      &                        gam     = cheb,               &
      &                        params  = fpt%legToChebParams )

    do iDof = 1, nIndeps*n, n
      ! Normalize the coefficients of the Chebyshev polynomials due
      ! to the unnormalized version of DCT in the FFTW.
      cheb(iDof+1:n+iDof-1:2) = -0.5_rk * cheb(iDof+1:n+iDof-1:2)
      cheb(iDof+2:n+iDof-1:2) =  0.5_rk * cheb(iDof+2:n+iDof-1:2)
    end do

    call fftw_execute_r2r( fpt%planChebToPnt, &
      &                    cheb,              &
      &                    pntVal             )

  end subroutine ply_legToPnt_vec
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev-Lobatto nodes.
  subroutine ply_legToPnt_lobatto_single( fpt, legCoeffs, pntVal, nIndeps )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(inout) :: legCoeffs(:)
    class(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n)
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    do iDof = 1, nIndeps*n, n
      call ply_fpt_single( alph   = legCoeffs(iDof:iDof+n-1), &
        &                  gam    = cheb,                     &
        &                  params = fpt%legToChebParams       )

      ! Normalize the coefficients of the Chebyshev polynomials due
      ! to the unnormalized version of DCT in the FFTW.
      cheb(2:n-1) = 0.5_rk * cheb(2:n-1)

      call fftw_execute_r2r( fpt%planChebToPnt,    &
        &                    cheb,                 &
        &                    pntVal(iDof:iDof+n-1) )
    end do

  end subroutine ply_legToPnt_lobatto_single
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Vectorizing subroutine to transform Legendre expansion to point values
  !! at Chebyshev-Lobatto nodes.
  subroutine ply_legToPnt_lobatto_vec( fpt, legCoeffs, pntVal, nIndeps )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(inout) :: legCoeffs(:)
    class(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n*nIndeps)
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    call ply_fpt_exec_striped( nIndeps = nIndeps,            &
      &                        alph    = legCoeffs,          &
      &                        gam     = cheb,               &
      &                        params  = fpt%legToChebParams )

    do iDof = 1, nIndeps*n, n
      ! Normalize the coefficients of the Chebyshev polynomials due
      ! to the unnormalized version of DCT in the FFTW.
      cheb(iDof+1:n+iDof-2) = 0.5_rk * cheb(iDof+1:n+iDof-2)
    end do

    call fftw_execute_r2r( fpt%planChebToPnt, &
      &                    cheb,              &
      &                    pntVal             )

  end subroutine ply_legToPnt_lobatto_vec
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to transform point values at Chebyshev nodes to a Legendre
  !! expansion.
  subroutine ply_pntToLeg_single( fpt, pntVal, legCoeffs, nIndeps )
    ! -------------------------------------------------------------------- !
    class(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    real(kind=rk), intent(inout) :: legCoeffs(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n)
    real(kind=rk) :: normFactor
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    normFactor = 1.0_rk / real(n,kind=rk)
    do iDof = 1, nIndeps*n, n
      call fftw_execute_r2r( fpt%planPntToCheb,     &
        &                    pntVal(iDof:iDof+n-1), &
        &                    cheb                   )
      ! Normalize the coefficients of the Chebyshev polynomials due
      ! to the unnormalized version of DCT in the FFTW.
      cheb(1) = cheb(1) * 0.5_rk * normfactor
      cheb(2:n:2) = -normFactor * cheb(2:n:2)
      cheb(3:n:2) =  normFactor * cheb(3:n:2)

      call ply_fpt_single( gam    = legCoeffs(iDof:iDof+n-1), &
        &                  alph   = cheb,                     &
        &                  params = fpt%ChebToLegParams       )
    end do

  end subroutine ply_pntToLeg_single
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Vectorizing subroutine to transform point values at Chebyshev nodes to a
  !! Legendre expansion.
  subroutine ply_pntToLeg_vec( fpt, pntVal, legCoeffs, nIndeps )
    ! -------------------------------------------------------------------- !
    class(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    real(kind=rk), intent(inout) :: legCoeffs(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n*nIndeps)
    real(kind=rk) :: normFactor
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    call fftw_execute_r2r( fpt%planPntToCheb, &
      &                    pntVal,            &
      &                    cheb               )

    normFactor = 1.0_rk / real(n,kind=rk)
    do iDof = 1, nIndeps*n, n
      ! Normalize the coefficients of the Chebyshev polynomials due
      ! to the unnormalized version of DCT in the FFTW.
      cheb(iDof) = cheb(iDof) * 0.5_rk * normfactor
      cheb(iDof+1:iDof+n-1:2) = -normFactor * cheb(iDof+1:iDof+n-1:2)
      cheb(iDof+2:iDof+n-1:2) =  normFactor * cheb(iDof+2:iDof+n-1:2)
    end do

    call ply_fpt_exec_striped( nIndeps = nIndeps,            &
      &                        alph    = cheb,               &
      &                        gam     = legCoeffs,          &
      &                        params  = fpt%ChebToLegParams )

  end subroutine ply_pntToLeg_vec
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine to transform point values at Chebyshev-Lobatto nodes to a
  !! Legendre expansion.
  subroutine ply_pntToLeg_lobatto_single( fpt, pntVal, legCoeffs, nIndeps )
    ! -------------------------------------------------------------------- !
    class(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    real(kind=rk), intent(inout) :: legCoeffs(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n)
    real(kind=rk) :: normFactor
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    normFactor = 0.5_rk / real(n-1,kind=rk)
    do iDof = 1, nIndeps*n, n
      call fftw_execute_r2r( fpt%planPntToCheb,     &
        &                    pntVal(iDof:iDof+n-1), &
        &                    cheb                   )
      ! Normalize the coefficients of the Chebyshev polynomials due
      ! to the unnormalized version of DCT in the FFTW.
      cheb(1) = cheb(1) * normFactor
      cheb(2:n-1) = 2.0_rk * normFactor * cheb(2:n-1)
      cheb(n) = cheb(n) * normFactor

      call ply_fpt_single( gam    = legCoeffs(iDof:iDof+n-1), &
        &                  alph   = cheb,                     &
        &                  params = fpt%ChebToLegParams       )
    end do

  end subroutine ply_pntToLeg_lobatto_single
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Vectorizing subroutine to transform point values at Chebyshev-Lobatto
  !! nodes to a Legendre expansion.
  subroutine ply_pntToLeg_lobatto_vec( fpt, pntVal, legCoeffs, nIndeps )
    ! -------------------------------------------------------------------- !
    class(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    real(kind=rk), intent(inout) :: legCoeffs(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n*nIndeps)
    real(kind=rk) :: normFactor
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    call fftw_execute_r2r( fpt%planPntToCheb, &
      &                    pntVal,            &
      &                    cheb               )

    normFactor = 0.5_rk / real(n-1,kind=rk)
    do iDof = 1, nIndeps*n, n
      ! Normalize the coefficients of the Chebyshev polynomials due
      ! to the unnormalized version of DCT in the FFTW.
      cheb(iDof) = cheb(iDof) * normFactor
      cheb(iDof+1:iDof+n-2) = 2.0_rk * normFactor * cheb(iDof+1:iDof+n-2)
      cheb(iDof+n-1) = cheb(iDof+n-1) * normFactor
    end do

    call ply_fpt_exec_striped( nIndeps = nIndeps,            &
      &                        alph    = cheb,               &
      &                        gam     = legCoeffs,          &
      &                        params  = fpt%ChebToLegParams )

  end subroutine ply_pntToLeg_lobatto_vec
  ! ------------------------------------------------------------------------ !

end module ply_legFpt_module
