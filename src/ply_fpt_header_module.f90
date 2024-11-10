! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2013-2014, 2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
!
! Parts of this file were written by Verena Krupp, Harald Klimach, Peter Vitt,
! Nikhil Anand, Daniel Petró and Neda Ebrahimi Pour for University of Siegen.
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
!> Parameters for the FPT method
!!
!! The Fast Polynomial Transformation implements the approach described in
!! B. K. Alpert und V. Rokhlin, „A Fast Algorithm for the Evaluation of
!! Legendre Expansions“, SIAM Journal on Scientific and Statistical Computing,
!! Vol. 12, Nr. 1, pp. 158–179, Jan. 1991, doi: 10.1137/0912009.
!!
!! It utilizes the Fast Fourier Transformation by first converting the
!! Legendre modes into Chebyshev modes.
!! The conversion between Legendre and Chebyshev modes is done approximately
!! by approximating increasingly larger blocks away from the diagonal.
!! This method is only available if the executable is linked against the
!! [FFTW](http://fftw.org/).
!!
!! As with the other projection methods, a `factor` can be specified to
!! use more points in the nodal representation and achieve some de-aliasing.
!! Because the FFT works especially well for powers of two, it is possible to
!! choose the oversampling such, that the oversampled modes have count of the
!! next larger power of 2.
!! To achieve this, the option `adapt_factor_pow2` has to be set to `true`,
!! by default it is assumed to be `false`.
!! Further it is possible to make use of `lobattoPoints` with this
!! transformation method.
!! If `lobattoPoints = true`, the points on the interval boundary will be
!! included in the set of points in the nodal representation. This is for
!! example necessary when positivity is to be preserved for the numerical
!! fluxes.
!! By default this setting is `false`.
!!
!! All other options tune the transformation algorithm:
!!
!! * `blocksize` defines the minimal block size that is to be approximated in
!!   the transformation matrix. It defaults to 64, which is the recommendation
!!   for double precision computations, but requires high polynomial degrees to
!!   attain any approximation at all. The (oversampled) number of nodes needs to
!!   be larger than two times the blocksize, to have at least one approximated
!!   block. As long as the number modes is below this threshold the method is
!!   not "fast" and a computational complexity of number of modes squared is
!!   required for the operation. Smaller values push the approximation closer
!!   to the diagonal. This can speed up the computation for smaller number of
!!   modes but also detoriates the accuracy of the transformation.
!! * `approx_terms` number of terms to use in the approximation of blocks,
!!   defaults to 18, which is recommended for double precision computations.
!!   In each block only `approx_terms` will be used to represent rows in the
!!   block. Smaller values make the transformation faster, but less accurate
!!   (if blocks are actually approximated).
!! * `implementation` selects the implementation for the transformation. There
!!   are two variants of the implementation: `'scalar'` this is the default
!!   and treats the transformations with an outer loop over the independent
!!   operations. The `'vector'` implementation on the other hand gathers
!!   multiple independent operations together and performs them all at once
!!   with an inner loop over blocks length `striplen`. While the `'vector'`
!!   variant may exploit vector instructions, it utilizes a larger amount of
!!   temporary memory. The default setting for this option is `'scalar'`.
!! * `striplen` determines the length for vectorized loops to be used in the
!!   matrix operation. It defaults to the `vlen` setting defined in
!!   [[tem_compileconf_module]] during compilation.
!!   Depending on the computing architecture, different values may provide
!!   more efficient computations.
!! * `subblockingWidth` defines striding in the multiplication of the diagonal
!!   elements in the transformation matrix. The default for this setting is 8.
!!
!! The configuration table for the FPT table may, for example, look as follows:
!!
!!```lua
!!  projection = {
!!    kind              = 'fpt',
!!    factor            = 1.5,
!!    adapt_factor_pow2 = true,
!!    lobattoPoints     = false,
!!    blocksize         = 16,
!!    approx_terms      = 12,
!!    implementation    = 'scalar',
!!    striplen          = 256,
!!    subblockingWidth  = 8
!!  }
!!```
module ply_fpt_header_module

  use env_module,              only: rk, labelLen
  use aotus_module,            only: flu_State, aot_get_val
  use aot_out_module,          only: aot_out_type, aot_out_val

  use tem_aux_module,          only: tem_abort
  use tem_tools_module,        only: upper_to_lower
  use tem_logging_module,      only: logUnit
  use tem_compileconf_module,  only: vlen
  use tem_float_module
  use ply_nodes_header_module

  implicit none

  private

  !> The recommended minimal blocksize for double precision.
  integer, public, parameter :: ply_fpt_default_blocksize = 64

  !> The default width of the subblocking of the diagonal calculation of the
  !! fpt projection
  integer, public, parameter :: ply_fpt_default_subblockingWidth = 8

  !> Default number of terms to use in FPT blocks. 18 is recommended for
  !! double precision.
  integer, public, parameter :: ply_fpt_default_approx_terms = 18

  !> Value to signify the use of the scalar FPT implementation.
  integer, public, parameter :: ply_fpt_scalar = 1

  !> Value to signify the use of the vector FPT implementation.
  integer, public, parameter :: ply_fpt_vector = 2

  !> Type for the fpt header, stores all information needed to initialize the
  !! fpt method later on
  type ply_fpt_header_type
    type(ply_nodes_header_type) :: nodes_header
    !> In case of nonlinear equations, aliasing occurs if the projections
    !! of the nonlinear terms on the testfunctions are not calculated
    !! accurately enough. To avoid these errors it is possible to
    !! extend the transformation vectors of the FPT with zeros. This
    !! factor determines by how many zeros the modal vector is extended
    !! before transformation. This factor has to be chosen properly with
    !! respect of the type of nonlinearity of your equation.
    real(kind=rk) :: factor = 1.0_rk

    !> The blockisze of the fast bases exchange algorithm from
    !! Legendre to Chebyshev polynomials.
    !! A negative number indicates to use the default blocksize of the
    !! algorithm.
    integer :: blocksize = ply_fpt_default_blocksize

    !> The number of approximation terms to use for blocks apart from the
    !! diagonal.
    !!
    !! This defaults to 18, which is recommended for double precision.
    integer :: approx_terms = ply_fpt_default_approx_terms

    !> The implementation variant to use for the transformation computation.
    !!
    !! The computation can be done either by a `'vector'` implementation or
    !! by a `'scalar'` implementation.
    !! We indicate the respective implementations by the integers
    !! [[ply_fpt_header_module:ply_fpt_scalar]] or
    !! [[ply_fpt_header_module:ply_fpt_vector]].
    integer :: implementation

    !> The striplen, that should be used for vectorized simultaneous
    !! computations of the matrix operation.
    !!
    !! This defaults to the vlen from the TEM_compileconf_module, it might
    !! be set differently here, as we are dealing with a twodimensional
    !! problem here, and the optimal setting might be different from the code
    !! parts.
    integer :: striplen = vlen

    !> The width of the subblocks used during the unrolled base exchange to
    !! ensure a better cache usage.
    !!
    !! The default is a subblocking width of 8.
    integer :: subblockingWidth = ply_fpt_default_subblockingWidth

    !> Should the oversampling factor be adapted to ensure a power of 2
    !! in the oversampled polynomial?
    !!
    !! If this is true, the factor will be increased to ensure
    !! an oversampled representation with a power of 2.
    !! Default is false.
    logical :: adapt_factor_pow2 = .false.
  end type ply_fpt_header_type

  interface assignment(=)
    module procedure Copy_fpt_header
  end interface

  interface operator(==)
    module procedure isEqual
  end interface

  interface operator(/=)
    module procedure isUnequal
  end interface

  interface operator(<)
    module procedure isSmaller
  end interface

  interface operator(<=)
    module procedure isSmallerOrEqual
  end interface

  interface operator(>)
    module procedure isGreater
  end interface

  interface operator(>=)
    module procedure isGreaterOrEqual
  end interface

  public :: operator(==), operator(/=), operator(<), operator(<=)
  public :: operator(>), operator(>=)

  public :: ply_fpt_header_load, ply_fpt_header_display
  public :: ply_fpt_header_define
  public :: ply_fpt_header_type
  public :: ply_fpt_header_out

  public :: assignment(=)


contains


  ! ------------------------------------------------------------------------ !
  !> Copy the FPT header information.
  pure subroutine Copy_fpt_header(left,right)
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_Fpt_header_type), intent(out) :: left
    !> fpt to copy from
    type(ply_Fpt_header_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%nodes_header = right%nodes_header
    left%factor = right%factor
    left%blocksize = right%blocksize
    left%adapt_factor_pow2 = right%adapt_factor_pow2
    left%approx_terms = right%approx_terms
    left%implementation = right%implementation
    left%striplen = right%striplen
    left%subblockingWidth = right%subblockingWidth

  end subroutine Copy_fpt_header
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Read the FPT configuration options from the provided Lua script in
  !! `conf`.
  subroutine ply_fpt_header_load(me, conf, thandle)
    ! -------------------------------------------------------------------- !
    type(ply_fpt_header_type), intent(out) :: me
    type(flu_State), intent(inout) :: conf
    integer, intent(in) :: thandle
    ! -------------------------------------------------------------------- !
    integer :: iError
    character(len=labelLen) :: impl_variant
    ! -------------------------------------------------------------------- !
    ! check for fpt lib

    ! for fpt chebyshev nodes are used
    me%nodes_header%nodes_kind = 'chebyshev'

    ! fill up the fpt_header
    call aot_get_val( L       = conf,                      &
      &               thandle = thandle,                   &
      &               key     = 'blocksize',               &
      &               val     = me%blocksize ,             &
      &               default = ply_fpt_default_blocksize, &
      &               ErrCode = iError                     )

    if (me%blocksize <= 0) then
      write(logUnit(1),*) 'ERROR in loading projection: blocksize for FPT has' &
        & // ' to be larger than 0!'
      write(logUnit(1),*) 'But it is set to ', me%blocksize
      call tem_abort()
    end if

    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'factor',  &
      &               val     = me%factor, &
      &               default = 1.0_rk,    &
      &               ErrCode = iError     )

    if (me%factor <= 0) then
      write(logUnit(1),*) 'ERROR in loading projection: factor for projection' &
        & // ' has to be larger than 0!'
      write(logUnit(1),*) 'But it is set to ', me%factor
      call tem_abort()
    end if

    call aot_get_val( L       = conf,                        &
      &               thandle = thandle,                     &
      &               key     = 'approx_terms',              &
      &               val     = me%approx_terms,             &
      &               ErrCode = iError,                      &
      &               default = ply_fpt_default_approx_terms )

    call aot_get_val( L       = conf,             &
      &               thandle = thandle,          &
      &               key     = 'implementation', &
      &               val     = impl_variant,     &
      &               ErrCode = iError,           &
      &               default = 'scalar'          )
    impl_variant = upper_to_lower(impl_variant)
    impl_variant = adjustl(impl_variant)
    select case(trim(impl_variant))
    case('scalar')
      me%implementation = ply_fpt_scalar
    case('vector')
      me%implementation = ply_fpt_vector
    case default
      write(logUnit(1),*) 'ERROR in loading projection: implementation' &
        & // ' has to be either "scalar" or "vector"!'
      write(logUnit(1),*) 'But it is set to ', trim(impl_variant)
      call tem_abort()
    end select

    call aot_get_val( L       = conf,        &
      &               thandle = thandle,     &
      &               key     = 'striplen',  &
      &               val     = me%striplen, &
      &               ErrCode = iError,      &
      &               default = vlen         )

    call aot_get_val( L       = conf,                            &
      &               thandle = thandle,                         &
      &               key     = 'subblockingWidth',              &
      &               val     = me%subblockingWidth,             &
      &               ErrCode = iError,                          &
      &               default = ply_fpt_default_subblockingWidth )

    write(logUnit(1), *) 'subblockingWidth = ', me%subblockingWidth

    call aot_get_val( L       = conf,                 &
      &               thandle = thandle,              &
      &               key     = 'adapt_factor_pow2',  &
      &               val     = me%adapt_factor_pow2, &
      &               ErrCode = iError,               &
      &               default = .false.               )

    ! check for lobatto Points
    call aot_get_val( L       = conf,                          &
      &               thandle = thandle,                       &
      &               key     = 'lobattoPoints',               &
      &               val     = me%nodes_header%lobattoPoints, &
      &               ErrCode = iError,                        &
      &               default = .false.                        )

  end subroutine ply_fpt_header_load
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Define settings for the Fast Polynomial Transformation.
  subroutine ply_fpt_header_define(                                           &
    &          me, blocksize, factor, approx_terms, implementation, striplen, &
    &          subBlockingWidth, adapt_factor_pow2, lobattoPoints             )
    ! -------------------------------------------------------------------- !
    !> FPT header to hold the defined settings.
    type(ply_fpt_header_type), intent(out) :: me

    !> Blocksize to use in approximation algorithm. Defaults to
    !! [[ply_fpt_header_module:ply_fpt_default_blocksize]].
    integer, optional, intent(in) :: blocksize

    !> Oversampling factor to use.
    !!
    !! This can be used to reduce aliasing when projecting functions that
    !! will be truncated in the polynomial series expansion.
    !! Default is a factor of 1, so no oversampling.
    real(kind=rk), optional, intent(in) :: factor

    !> Number of approximation terms to use for the representation of the
    !! blocks in the Legendre to Chebyshev transformation algorithm.
    !! Defaults to [[ply_fpt_header_module:ply_fpt_default_approx_terms]].
    integer, optional, intent(in) :: approx_terms

    !> Implementation to use in the computation.
    !!
    !! Select the implementation variant to use. Either scalar
    !! ([[ply_fpt_header_module:ply_fpt_scalar]]) or vectorized
    !! ([[ply_fpt_header_module:ply_fpt_vector]]).
    !! Default is [[ply_fpt_header_module:ply_fpt_scalar]].
    integer, optional, intent(in) :: implementation

    !> Length of strips to use in the transformation implementation.
    !! Defaults to [[tem_compileconf_module:vlen]].
    integer, optional, intent(in) :: striplen

    !> Width for subblocks in unrolling the approximate Legendre to
    !! Chebyshev transformation. Defaults to
    !! [[ply_fpt_header_module:ply_fpt_default_subblockingWidth]].
    integer, optional, intent(in) :: subBlockingWidth

    !> Adapt the oversampling factor such, that oversampled space has a
    !! number of degrees of freedoms in one direction that is a power of 2.
    !! Default is .false..
    logical, optional, intent(in) :: adapt_factor_pow2

    !> Wether to use Chebyshev-Lobatto points (include boundary points) or
    !! not. Defaults to .false..
    logical, optional, intent(in) :: lobattoPoints
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ! for fpt chebyshev nodes are used
    me%nodes_header%nodes_kind = 'chebyshev'

    ! Defaults
    me%blocksize         = ply_fpt_default_blocksize
    me%factor            = 1.0_rk
    me%approx_terms      = ply_fpt_default_approx_terms
    me%implementation    = ply_fpt_scalar
    me%striplen          = vlen
    me%subBlockingWidth  = ply_fpt_default_subblockingWidth
    me%adapt_factor_pow2 = .false.

    me%nodes_header%lobattoPoints = .false.

    ! Overwrite defaults if set by caller.
    if (present(blocksize)) me%blocksize = blocksize
    if (present(factor)) me%factor = factor
    if (present(approx_terms)) me%approx_terms = approx_terms
    if (present(implementation)) me%implementation = implementation
    if (present(striplen)) me%striplen = striplen
    if (present(subBlockingWidth)) me%subBlockingWidth = subBlockingWidth
    if (present(adapt_factor_pow2)) me%adapt_factor_pow2 = adapt_factor_pow2

    if (present(lobattoPoints)) me%nodes_header%lobattoPoints = lobattoPoints

    if (me%factor <= 0) then
      write(logUnit(1),*) 'ERROR in defining projection: factor for projection' &
        & // ' has to be larger than 0!'
      write(logUnit(1),*) 'But it is set to ', me%factor
      call tem_abort()
    end if

    if ( me%implementation /=  ply_fpt_scalar &
      &  .and. me%implementation /= ply_fpt_vector ) then
      write(logUnit(1),*) 'ERROR in defining projection: implementation' &
        & // ' has to be either "scalar" or "vector"!'
      write(logUnit(1),*) 'But it is set to unknown value ', me%implementation
      call tem_abort()
    end if

  end subroutine ply_fpt_header_define
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Write FPT settings into a Lua table.
  subroutine ply_fpt_header_out(me, conf)
    ! -------------------------------------------------------------------- !
    type(ply_fpt_header_type), intent(in) :: me
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ! fill up the fpt_header
    call aot_out_val( put_conf = conf,        &
       &              vname    = 'blocksize', &
       &              val      = me%blocksize )

    call aot_out_val( put_conf = conf,     &
       &              vname    = 'factor', &
       &              val      = me%factor )

    call aot_out_val( put_conf = conf,           &
       &              vname    = 'approx_terms', &
       &              val      = me%approx_terms )

    select case(me%implementation)
    case(ply_fpt_scalar)
      call aot_out_val( put_conf = conf,             &
         &              vname    = 'implementation', &
         &              val      = 'scalar'          )
    case(ply_fpt_vector)
      call aot_out_val( put_conf = conf,             &
         &              vname    = 'implementation', &
         &              val      = 'vector'          )
    end select

    call aot_out_val( put_conf = conf,       &
       &              vname    = 'striplen', &
       &              val      = me%striplen )

    call aot_out_val( put_conf = conf,               &
       &              vname    = 'subblockingWidth', &
       &              val      = me%subblockingWidth )

    call aot_out_val( put_conf = conf,                &
       &              vname    = 'adapt_factor_pow2', &
       &              val      = me%adapt_factor_pow2 )

    call aot_out_val( put_conf = conf,                         &
       &              vname    = 'lobattoPoints',              &
       &              val      = me%nodes_header%lobattoPoints )

  end subroutine ply_fpt_header_out
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Print the FPT settings to the log output.
  subroutine ply_fpt_header_display (me)
    ! -------------------------------------------------------------------- !
    type(ply_fpt_header_type), intent(in) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) ' Using fast polynomial transforms for projection.'
    write(logUnit(1),*)
    write(logUnit(1),*) ' * Kind of projection method = fpt'
    write(logUnit(1),*) ' * Dealising factor to use in projection = ', &
      &                 me%factor
    write(logUnit(1),*) ' * Adapt factor to ensure power of 2 order = ', &
      &                 me%adapt_factor_pow2
    write(logUnit(3),*) '     This setting is only relevant for' &
      &                 //' polynomial degrees < 2*blocksize.'
    write(logUnit(1),*) ' * using LobattoPoints =', &
      &                 me%nodes_header%lobattoPoints
    if (me%implementation == ply_fpt_scalar) then
      write(logUnit(1),*) ' * Using scalar implementation'
    else
      write(logUnit(1),*) ' * Using VECTOR implementation'
    end if
    write(logUnit(1),*) ' * Block approximation:'
    write(logUnit(1),*) '   * Blocksize for FPT =', me%blocksize
    write(logUnit(1),*) '   * Number of approximation terms = ', me%approx_terms
    write(logUnit(1),*) '   * Strip length = ', me%striplen
    write(logUnit(1),*) '   * Subblocking width = ', me%subblockingWidth
    write(logUnit(1),*) ''
  end subroutine ply_fpt_header_display
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function provides the test for equality of two projections.
  !!
  !! Two fpt header are considered to be equal, if their  node_header,
  !! fpt_blocksize or the factor are equal.
  pure function isEqual(left, right) result(equality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is equal??
    logical :: equality
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    equality = ( left%nodes_header == right%nodes_header )           &
      & .and. ( left%factor .feq. right%factor )                     &
      & .and. ( left%blocksize == right%blocksize )                  &
      & .and. ( left%approx_terms == right%approx_terms )            &
      & .and. ( left%striplen == right%striplen )                    &
      & .and. ( left%subblockingWidth == right%subblockingWidth )    &
      & .and. ( left%adapt_factor_pow2 .eqv. right%adapt_factor_pow2 )

  end function isEqual
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function provides the test for unequality of two projections.
  !!
  !! Two fpt header are considered to be unequal, if their  node_header,
  !! fpt_blocksize or the factor are not equal.
  pure function isUnequal(left, right) result(unequality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    unequality = ( left%nodes_header /= right%nodes_header )         &
      & .or. ( left%factor .fne. right%factor )                      &
      & .or. ( left%blocksize /= right%blocksize )                   &
      & .or. ( left%approx_terms /= right%approx_terms )             &
      & .or. ( left%striplen /= right%striplen )                     &
      & .or. ( left%subblockingWidth /= right%subblockingWidth )     &
      & .or. ( left%adapt_factor_pow2 .neqv. right%adapt_factor_pow2 )

  end function isUnequal
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function provides a < comparison of two projections.
  !!
  !! Sorting of fpt header is given by node_header, fpt_blocksize and
  !! last by factor.
  pure function isSmaller(left, right) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !

    small = .false.

    if (left%adapt_factor_pow2) then
      left_log = 1
    else
      left_log = 0
    end if
    if (right%adapt_factor_pow2) then
      right_log = 1
    else
      right_log = 0
    end if

    small = left%nodes_header < right%nodes_header
    if (left%nodes_header == right%nodes_header) then
      small = left%factor < right%factor
      if (left%factor .feq. right%factor) then
        small = left%blocksize < right%blocksize
        if (left%blocksize == right%blocksize) then
          small = left%approx_terms < right%approx_terms
          if ( left%approx_terms ==right%approx_terms) then
            small = left%striplen < right%striplen
            if (left%striplen == right%striplen) then
              small = left%subblockingWidth < right%subblockingWidth
              if (left%subblockingWidth == right%subblockingWidth) then
                small = (left_log < right_log)
              end if
            end if
          end if
        end if
      end if
    end if

  end function isSmaller
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function provides a <= comparison of two projections.
  !!
  !! Sorting of fpt header is given by node_header, fpt_blocksize and
  !! last by factor.
  pure function isSmallerOrEqual(left, right) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !

    small = .false.

    if (left%adapt_factor_pow2) then
      left_log=1
    else
      left_log=0
    end if
    if (right%adapt_factor_pow2) then
      right_log=1
    else
      right_log=0
    end if

    small = left%nodes_header < right%nodes_header
    if (left%nodes_header == right%nodes_header) then
      small = left%factor < right%factor
      if (left%factor .feq. right%factor) then
        small = left%blocksize < right%blocksize
        if (left%blocksize == right%blocksize) then
          small = left%approx_terms < right%approx_terms
          if ( left%approx_terms ==right%approx_terms) then
            small = left%striplen < right%striplen
            if (left%striplen == right%striplen) then
              small = left%subblockingWidth < right%subblockingWidth
              if (left%subblockingWidth == right%subblockingWidth) then
                small = (left_log <= right_log)
              end if
            end if
          end if
        end if
      end if
    end if

  end function isSmallerOrEqual
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function provides a > comparison of two projections.
  !!
  !! Sorting of fpt header is given by node_header, fpt_blocksize and
  !! last by factor.
  pure function isGreater(left, right) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !

    great = .false.

    if (left%adapt_factor_pow2) then
      left_log=1
    else
      left_log=0
    end if
    if (right%adapt_factor_pow2) then
      right_log=1
    else
      right_log=0
    end if

    great = left%nodes_header > right%nodes_header
    if (left%nodes_header == right%nodes_header) then
      great = left%factor > right%factor
      if (left%factor .feq. right%factor) then
        great = left%blocksize > right%blocksize
        if (left%blocksize == right%blocksize) then
          great = left%approx_terms > right%approx_terms
          if ( left%approx_terms ==right%approx_terms) then
            great = left%striplen > right%striplen
            if (left%striplen == right%striplen) then
              great = left%subblockingWidth > right%subblockingWidth
              if (left%subblockingWidth == right%subblockingWidth) then
                great = (left_log > right_log)
              end if
            end if
          end if
        end if
      end if
    end if

  end function isGreater
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function provides a >= comparison of two projections.
  !!
  !! Sorting of fpt header is given by node_header, fpt_blocksize and
  !! last by factor.
  pure function isGreaterOrEqual(left, right) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !

    great = .false.

    if (left%adapt_factor_pow2) then
      left_log=1
    else
      left_log=0
    end if
    if (right%adapt_factor_pow2) then
      right_log=1
    else
      right_log=0
    end if

    great = left%nodes_header > right%nodes_header
    if (left%nodes_header == right%nodes_header) then
      great = left%factor > right%factor
      if (left%factor .feq. right%factor) then
        great = left%blocksize > right%blocksize
        if (left%blocksize == right%blocksize) then
          great = left%approx_terms > right%approx_terms
          if ( left%approx_terms ==right%approx_terms) then
            great = left%striplen > right%striplen
            if (left%striplen == right%striplen) then
              great = left%subblockingWidth > right%subblockingWidth
              if (left%subblockingWidth == right%subblockingWidth) then
                great = (left_log >= right_log)
              end if
            end if
          end if
        end if
      end if
    end if
  end function isGreaterOrEqual
  ! ------------------------------------------------------------------------ !

end module ply_fpt_header_module
