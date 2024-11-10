! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2018, 2020 Harald Klimach <harald.klimach@uni-siegen.de.de>
! Copyright (c) 2013-2014 Verena Krupp <v.krupp@grs-sim.de>
! Copyright (c) 2015 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
! Copyright (c) 2016-2017, 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop for German Research School
! for Simulation Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Verena Krupp, Peter Vitt,
! Tobias Girresser, Nikhil Anand, Kay Langhammer, Kannan Masilamani and Neda
! Ebrahimi Pour for University of Siegen.
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
!> The polynomial projection describes the change from modal to nodal
!! representation and vice-versa.
!!
!! The method to be used for this transformation is configured in a table that
!! contains the `kind` of transformation and then further options to the
!! transformation.
!! There are three kinds of transformation available:
!!
!! * 'fpt' Fast polynomial transformation, only available if linked against
!!   FFTW.
!! * 'l2p' Standard projection, evaluate each polynomial for each point.
!! * 'fxt' Fast transformation based on multipole method, enabled by the
!!   (included) FXTPACK library.
!!
!! By default the nodal representation uses as many nodes as there are modes
!! in the modal representation. This can be changed with the `factor` setting.
!! There will `factor` times number of modes nodes used in the transformation.
!! This oversampling allows for dealiasing of the nodal representation.
!! A different set of points will be used depending on the kind of
!! transformation to be used:
!!
!! * FPT uses Chebyshev integration nodes
!! * L2P may make use of either Chebyshev or Gauss-Legendre nodes
!! * FXT Gauss-Legendre integration nodes
!!
!! If the Chebyshev nodes are used for the nodal representation, the interval
!! boundaries can be included by setting `lobattoPoints=true`. In this case
!! there will also be at least one point more be used in the nodal
!! representation than there are modes in the Legendre series.
!!
!! If FPT is configured but the executable is not linked against the FFTW, there
!! will be a warning, and the simulation uses the L2P method instead.
!! The L2P method is also the default when no kind is provided at all.
!!
!! For further options of the individual projection kinds, see their respective
!! descriptions:
!!
!! * L2P: [[ply_l2p_header_module]]
!! * FPT: [[ply_fpt_header_module]]
!! * FXT: [[ply_fxt_header_module]]
!!
!! A minimalistic configuration for the projection method looks as follows:
!!
!!```lua
!!  projection = {
!!    kind = 'l2p',
!!    factor = 1.5
!!  }
!!```
!!
!! The name of the table (here `projection`) is arbitrary and defined by the
!! application.
!! An application may load multiple projection definitions with differing names.
!!
module ply_poly_project_module
  use env_module,                only: rk, labelLen

  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logUnit
  use tem_tools_module,          only: tem_horizontalSpacer
  use ply_dof_module,            only: Q_space, &
   &                                   P_space
  use ply_prj_header_module
  use ply_dynArray_project_module, only: dyn_ProjectionArray_type, &
    &                                    ply_fill_dynProjectArray, &
    &                                    ply_prj_init_type

  use ply_LegFpt_module,           only: ply_legFpt_type, &
    &                                    ply_init_legFpt, &
    &                                    assignment(=)
  use ply_legFpt_2D_module,        only: ply_pntToLeg_2D, &
    &                                    ply_legToPnt_2D
  use ply_legFpt_3D_module,        only: ply_pntToLeg_3D, &
    &                                    ply_legToPnt_3D
  use ply_l2p_module,              only: ply_l2p_type,     &
    &                                    ply_init_l2p,     &
    &                                    assignment(=),    &
    &                                    ply_l2p_trafo_3d, &
    &                                    ply_l2p_trafo_2d, &
    &                                    ply_l2p_trafo_1d
  use ply_nodes_header_module,     only: ply_nodes_header_type, &
    &                                    assignment(=)
  use ply_nodes_module,            only: ply_nodes_create, &
    &                                    ply_facenodes_type
  use ply_fxt_module,              only: ply_fxt_type,   &
    &                                    ply_init_fxt,   &
    &                                    ply_fxt_m2n_1D, &
    &                                    ply_fxt_m2n_3D, &
    &                                    ply_fxt_m2n_2D, &
    &                                    ply_fxt_n2m_1D, &
    &                                    ply_fxt_n2m_3D, &
    &                                    ply_fxt_n2m_2D

  implicit none

  private

  !> Additional data, required for the projection.
  type ply_prj_body_type
    !> The fast polynomial transformation which will be used in case
    !! of nonlinear equations. It is used if fpt is choses as projection
    !! method in the lua file
    type(ply_legFpt_type)  :: fpt
    !> The Legendre Polynomial type for the Fast Orthogonal Function
    !! Transform via fxtpack. It is used if 'fxt' is chosen as projection
    !! method in the lua file
    type(ply_fxt_type) :: fxt
    !> Projection method which cam be used for transfoamtion from modal to
    !! nodal space and vice versa. It is used if 'l2p' is chosen as projection
    !! method in the lua file
    type(ply_l2p_type) :: l2p
    !> Volume quadrature points in the reference element
    real(kind=rk), allocatable :: nodes(:,:)
    !> Facial quadrature nodes (reference element) for all three spatial
    !! direction and left and right face.
    !! These points are necessary to transfer boundary conditions given
    !! in physical space to modal space by projection (l2p or fpt)
    type(ply_faceNodes_type), allocatable :: faces(:,:)
    !> quadrature points including oversampling factor
    integer                    :: nQuadPoints
    !> degree of freedom of the scheme depending on maxPolyDegree
    integer                    :: ndofs
    ! the oversamp_dofs are the degrees of freedom for this 'oversampling
    ! degree and equal to (s*(m+1))**d
    integer                    :: oversamp_dofs
    ! minimal number of dofs, used to enable the flexibilty to
    ! set the oversampling factor < 1
    integer                    :: min_dofs

  end type ply_prj_body_type


  !> Projection definition.
  type ply_poly_project_type
    !> Polynomial basis type.
    !!
    !! 3D Monomials have the form x^i * y^j * z^k
    !! - Q_space: quadratic polynomial space (i,j,k) <= maxPolyDegree
    !! - P_space: polynomial space i+j+k <= maxPolyDegree
    integer                       :: basisType

    !> Kind of projection. Currently available:
    !! - 'l2p', L2-Projection
    !! - 'fpt', Fast Polynomial Transformation. Requires the FFTW.
    !! - 'fxt', Fast Polynomial Transformation. uses FXTPACK
    character(len=labelLen)       :: kind

    !> The maximal polynomial degree per spatial direction.
    integer                       :: maxPolyDegree

    !> Using oversampling, the modal space need to be extended according
    ! to the oversampling factor, thus the oversampling degree is (s*m+1)-1
    integer                       :: oversamp_degree

    ! minimal number of dofs, used to enable the flexibilty to
    ! set the oversampling factor < 1
    integer                       :: min_degree

    !> quadrature points including oversampling factor
    ! per spatial direction
    integer                       :: nQuadPointsPerDir

    !> Logical to indicate whether Chebyshev-Lobatto points or simple
    !! Chebyshev points are used
    logical                       :: lobattoPoints = .false.

    !> projection header consits of general information like which kind
    !! of projection is used

    !> In the body datatype, there is for each dimension the main data
    !! for the projection method stored
    type(ply_prj_body_type)   :: body_1d
    type(ply_prj_body_type)   :: body_2d
    type(ply_prj_body_type)   :: body_3d

  end type ply_poly_project_type


  interface assignment(=)
    module procedure Copy_poly_project
    module procedure Copy_poly_project_body
  end interface

  interface ply_poly_project_m2n
    module procedure ply_poly_project_m2n_multiVar
  end interface

  interface ply_poly_project_n2m
    module procedure ply_poly_project_n2m_multiVar
  end interface

  public :: assignment(=)
  public :: ply_fill_project_list
  public :: ply_poly_project_fillbody
  public :: ply_poly_project_m2n
  public :: ply_poly_project_n2m
  public :: ply_poly_project_type
  public :: ply_faceNodes_type
  public :: ply_get_quadpoints_faces
  public :: ply_get_quadpoints_faces_2d
  public :: ply_get_quadpoints_faces_1d
  public :: ply_prj_body_type


contains


  ! ------------------------------------------------------------------------ !
  subroutine Copy_poly_project(left,right)
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_poly_project_type), intent(out) :: left
    !> fpt to copy from
    type(ply_poly_project_type), intent(in) :: right
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    left%body_1d = right%body_1d
    left%body_2d = right%body_2d
    left%body_3d = right%body_3d

    left%kind = right%kind
    left%maxPolyDegree = right%maxPolyDegree
    left%basisType = right%basisType
    left%oversamp_degree = right%oversamp_degree
    left%min_degree = right%min_degree
    left%nquadpointsPerDir = right%nquadpointsPerDir
    left%lobattoPoints = right%lobattoPoints

  end subroutine copy_poly_project
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine Copy_poly_project_body(left,right)
    ! -------------------------------------------------------------------- !
    ! fpt to copy to
    type(ply_prj_body_type), intent(out) :: left
    ! fpt to copy from
    type(ply_prj_body_type), intent(in) :: right
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    left%fpt = right%fpt
    left%l2p = right%l2p
    left%fxt = right%fxt
    left%nodes = right%nodes
    left%faces = right%faces
    left%nquadpoints = right%nquadpoints
    left%ndofs = right%ndofs
    left%oversamp_dofs = right%oversamp_dofs

  end subroutine copy_poly_project_body
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Fill ups the bodys accroding to the DA.
  subroutine ply_fill_project_list( proj_list, dyn_projectArray, scheme_dim )
    ! -------------------------------------------------------------------- !
    type(ply_poly_project_type), intent(out), allocatable :: proj_list(:)
    type(dyn_ProjectionArray_type), intent(in) :: dyn_projectArray
    integer, intent(in) :: scheme_dim
    ! -------------------------------------------------------------------- !
    integer :: ipos
    ! -------------------------------------------------------------------- !
    call tem_horizontalSpacer(fUnit=logUnit(2))
    write(logUnit(2),*) 'Loading list of projection methods ... '

    ! allocate the poly_proj_list, the maximum value of all positions will give
    ! the number of elements in the DA
    allocate (proj_list(dyn_projectArray%nVals) )
    write(logUnit(5),*) 'the number of elements in proj_pos is =', &
      &                 dyn_projectArray%nVals
    do ipos=1, dyn_projectArray%nVals
      write(logUnit(5),*) 'for pos=', ipos, 'dyn array is',         &
        &                 dyn_projectArray%val(ipos)%basisType,     &
        &                 dyn_projectArray%val(ipos)%maxpolydegree, &
        &                 dyn_projectArray%val(ipos)%header%kind

      call ply_poly_project_fillbody(                 &
        &    me         = proj_list(ipos),            &
        &    proj_init  = dyn_projectArray%val(ipos), &
        &    scheme_dim = scheme_dim                  )

      write(logUnit(5),*) ' for position', ipos, &
        &                 'of projection list, projection type is'
      write(logUnit(5),*) '  kind = ', proj_list(ipos)%kind
      write(logUnit(5),*) '  maxPolyDegree =', proj_list(ipos)%maxPolyDegree
      write(logUnit(5),*) '  basisType = ', proj_list(ipos)%basisType
      write(logUnit(5),*) '  oversamp_degree = ', &
        &                 proj_list(ipos)%oversamp_degree
      write(logUnit(5),*) '  min_degree = ', proj_list(ipos)%min_degree
      write(logUnit(5),*) '  nquadpointsPerDir = ', &
        &                 proj_list(ipos)%nquadpointsPerDir
      write(logUnit(5),*) '  lobattoPoints = ', proj_list(ipos)%lobattoPoints

    end do

  end subroutine ply_fill_project_list
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Fill the body of the projection with all required data,
  !! ply_poly_project_define has to be used beforehand to set necessary header
  !! information.
  subroutine ply_poly_project_fillbody(me, proj_init, scheme_dim)
    ! -------------------------------------------------------------------- !
    type(ply_poly_project_type), intent(inout)  :: me
    type(ply_prj_init_type), intent (in)    :: proj_init
    integer, intent(in) :: scheme_dim
    ! -------------------------------------------------------------------- !
    ! the oversampling order, need to get the number of modes in when
    ! oversampling is used
    integer :: oversampling_order
    ! the number of volume quadrature points per spatial direction
    integer :: numQuadPointsPerDir
    real(kind=rk) :: log_order, rem_log
    real(kind=rk) :: over_factor
    integer :: lb_log
    type(ply_nodes_header_type) :: nodes_header
    ! -------------------------------------------------------------------- !
    ! set the kind in the final projection type
    me%kind = trim(proj_init%header%kind)
    me%maxPolyDegree = proj_init%maxPolyDegree
    me%basisType = proj_init%basisType

    select case(trim(proj_init%header%kind))
    case('fpt')
      over_factor = proj_init%header%fpt_header%factor
      nodes_header = proj_init%header%fpt_header%nodes_header
      me%lobattopoints = nodes_header%lobattopoints
    case('l2p')
      over_factor = proj_init%header%l2p_header%factor
      nodes_header = proj_init%header%l2p_header%nodes_header
      me%lobattopoints = nodes_header%lobattopoints
    case('fxt')
      over_factor = proj_init%header%fxt_header%factor
      nodes_header = proj_init%header%fxt_header%nodes_header
      me%lobattopoints = .false.
    end select

    ! Find the oversampling order
    oversampling_order = ceiling( over_factor*(me%maxpolyDegree+1) )
    if (trim(proj_init%header%kind) == 'fpt') then
      ! Ensure an appropriate order in the oversampled polynomial
      ! representation.
      if ( proj_init%header%fpt_header%adapt_factor_pow2               &
        &  .and. (iand(oversampling_order, oversampling_order-1) /= 0) &
        &  ) then

        write(logUnit(1),*) '*** NOTE: oversampling order is increased to' &
          &                 //' next power of 2! ***'
        write(logUnit(2),*) '          original oversampling order would' &
          &                 // ' have been: ', oversampling_order
        ! Oversampling_order is not a power of 2, find the next power of 2
        ! and use that one instead...
        log_order = log(real(oversampling_order, kind=rk))/log(2.0_rk)
        lb_log = max(floor(log_order), 0)
        rem_log = log_order - lb_log
        if (rem_log > epsilon(log_order)*lb_log) then
          oversampling_order = 2**(lb_log+1)
        else
          oversampling_order = 2**lb_log
        end if

      end if
    end if

    if (me%lobattoPoints) then
      ! If lobatto Points are to be used, use at least one point more
      ! in the nodal representation.
      oversampling_order = max(oversampling_order, me%maxPolyDegree+2)
    end if

    if (trim(proj_init%header%kind) == 'fxt') then
      ! Ensure an even order for the oversampled polynomial representation.
      ! There is a bug in FXTPACK resulting faulty conversions for odd
      ! numbers of evaluation points.
      ! To avoid this, we increase the oversampling order by 1, if it is odd.
      oversampling_order = oversampling_order + mod(oversampling_order,2)
    end if

    write(logUnit(1),*) 'Using an oversampled order of: ', oversampling_order
    write(logUnit(1),*) 'maxpolydegree is: ', me%maxPolyDegree
    write(logUnit(2),*) '(Actual oversampling factor: ',                   &
      &                 real(oversampling_order)/real(me%maxpolydegree+1), &
      &                 ')'

    numQuadPointsPerDir = oversampling_order

    me%basisType = proj_init%basisType
    me%maxPolyDegree = proj_init%maxPolyDegree
    me%nquadpointsPerDir = numQuadPointsPerDir
    me%oversamp_degree = oversampling_order-1
    me%min_degree = min(me%maxPolyDegree, me%oversamp_degree)

    ! number of dof depending on q_space or p_space
    if (me%basisType == Q_space) then
       me%body_3d%ndofs = (me%maxPolyDegree+1)**3
       me%body_2d%ndofs = (me%maxPolyDegree+1)**2
       me%body_1d%ndofs = me%maxPolyDegree+1
       me%body_3d%min_dofs = (me%min_degree+1)**3
       me%body_2d%min_dofs = (me%min_degree+1)**2
       me%body_1d%min_dofs = me%min_degree+1
    else !p_space
  me%body_3d%ndofs = (((me%maxpolydegree) + 1) &
    &   * ((me%maxpolydegree) + 2) &
    &   * ((me%maxpolydegree) + 3)) &
    & / 6
  me%body_2d%ndofs = ((me%maxpolydegree)+1)*((me%maxpolydegree)+2)/2
  me%body_1d%ndofs = ((me%maxpolydegree)+1)
  me%body_3d%min_dofs = (((me%min_degree) + 1) &
    &   * ((me%min_degree) + 2) &
    &   * ((me%min_degree) + 3)) &
    & / 6
  me%body_2d%min_dofs = ((me%min_degree)+1)*((me%min_degree)+2)/2
  me%body_1d%min_dofs = ((me%min_degree)+1)
    end if

    me%body_3d%nquadpoints = numQuadPointsPerDir**3
    me%body_2d%nquadpoints = numQuadPointsPerDir**2
    me%body_1d%nquadpoints = numQuadPointsPerDir

    me%body_3d%oversamp_dofs = (oversampling_order)**3
    me%body_2d%oversamp_dofs = (oversampling_order)**2
    me%body_1d%oversamp_dofs = oversampling_order

    select case (trim(proj_init%header%kind))
    case('fpt')
      ! Fill fpt datatype
      call ply_init_legfpt(                                 &
        &    maxPolyDegree    = me%oversamp_degree,         &
        &    nIndeps          = 1,                          &
        &    fpt              = me%body_1d%fpt,             &
        &    header           = proj_init%header%fpt_header )

      if (scheme_dim >= 2) then
        call ply_init_legfpt(                                 &
          &    maxPolyDegree    = me%oversamp_degree,         &
          &    nIndeps          = me%oversamp_degree+1,       &
          &    fpt              = me%body_2d%fpt,             &
          &    header           = proj_init%header%fpt_header )
      end if

      if (scheme_dim >= 3) then
        call ply_init_legfpt(                                 &
          &    maxPolyDegree    = me%oversamp_degree,         &
          &    nIndeps          = (me%oversamp_degree+1)**2,  &
          &    fpt              = me%body_3D%fpt,             &
          &    header           = proj_init%header%fpt_header )
      end if

    case('l2p')
      ! Fill the L2 projection datatype
      if (scheme_dim >= 3) then
        call ply_init_l2p( l2p    = me%body_3d%l2p,              &
          &                header = proj_init%header%l2p_header, &
          &                degree = me%oversamp_degree           )
      end if

      if (scheme_dim >= 2) then
        call ply_init_l2p( l2p    = me%body_2d%l2p,              &
          &                header = proj_init%header%l2p_header, &
          &                degree = me%oversamp_degree           )
      end if

      call ply_init_l2p( l2p    = me%body_1d%l2p,              &
        &                header = proj_init%header%l2p_header, &
        &                degree = me%oversamp_degree           )

    case ('fxt')
      ! Fill the fxt Legendre Polynomial datatype
      if (scheme_dim >= 3) then
        call ply_init_fxt( fxt    = me%body_3d%fxt,              &
          &                header = proj_init%header%fxt_header, &
          &                degree = me%oversamp_degree           )
      end if

      if (scheme_dim >= 2) then
        call ply_init_fxt( fxt    = me%body_2d%fxt,              &
          &                header = proj_init%header%fxt_header, &
          &                degree = me%oversamp_degree           )
      end if
        call ply_init_fxt( fxt    = me%body_1d%fxt,              &
          &                header = proj_init%header%fxt_header, &
          &                degree = me%oversamp_degree           )

    case default
      write(logUnit(1),*) 'ERROR in initializing projection:'
      write(logUnit(1),*) 'Unknown projection method <', &
        &                 trim(proj_init%header%kind), &
        &                 '>'
      write(logUnit(1),*) 'Stopping....'
      call tem_abort()

    end select

    !> Initialization/Create of the volume quadrature nodes and the
    !! quadrature points on the face
    call ply_nodes_create(                           &
      &    me                = nodes_header,         &
      &    nodes             = me%body_1d%nodes,     &
      &    faces             = me%body_1d%faces,     &
      &    nQuadPointsPerDir = me%nQuadPointsPerDir, &
      &    nDims             = 1                     )

    if (scheme_dim >= 2) then
      call ply_nodes_create(                           &
        &    me                = nodes_header,         &
        &    nodes             = me%body_2d%nodes,     &
        &    faces             = me%body_2d%faces,     &
        &    nQuadPointsPerDir = me%nQuadPointsPerDir, &
        &    nDims             = 2                     )
    end if

    if (scheme_dim >= 3) then
      call ply_nodes_create(                           &
        &    me                = nodes_header,         &
        &    nodes             = me%body_3d%nodes,     &
        &    faces             = me%body_3d%faces,     &
        &    nQuadPointsPerDir = me%nQuadPointsPerDir, &
        &    nDims             = 3                     )
    end if

  end subroutine ply_poly_project_fillbody
  ! ------------------------------------------------------------------------ !


  ! ***************** MODAL to NODAL ***************************************** !
  ! ------------------------------------------------------------------------ !
  !> Convert nDoF modes to nodal values.
  subroutine ply_poly_project_m2n_multiVar(me, dim, nVars, modal_data, &
    &                                      nodal_data)
    ! -------------------------------------------------------------------- !
    type(ply_poly_project_type), intent(inout) :: me
    integer, intent(in) :: dim
    !> The number of variables to project. If a variable consists of more than
    !! one component, the number of components has to be passed. If there are
    !! more than one variable, the sum of all components has to be passed (e.g.
    !! 6 when there are two three-dimensional vectors).
    integer, intent(in) :: nVars
    real(kind=rk), intent(inout) :: modal_data(:,:)
    real(kind=rk), intent(inout) :: nodal_data(:,:)
    ! -------------------------------------------------------------------- !
    integer :: iVar
    ! -------------------------------------------------------------------- !

    select case(trim(me%kind))
    case ('l2p')
      ! for the projection modal to nodal, we do not need to distingusih
      ! between the spaces since the modal values results from the computation
      ! and the projection is on evluation of the nodes at the points with
      ! additional summation
      select case(dim)
      case (1)
        do iVar = 1, nVars
          call ply_l2p_trafo_1D( trafo = me%body_1D%l2p%leg2node, &
            &                    projected = nodal_data(:,iVar),  &
            &                    original  = modal_data(:,iVar)   )
        end do
      case (2)
        do iVar = 1, nVars
          call ply_l2p_trafo_2D( trafo = me%body_2D%l2p%leg2node, &
            &                    projected = nodal_data(:,iVar),  &
            &                    original  = modal_data(:,iVar)   )
        end do
      case (3)
        do iVar = 1, nVars
          call ply_l2p_trafo_3D( trafo = me%body_3D%l2p%leg2node, &
            &                    projected = nodal_data(:,iVar),  &
            &                    original  = modal_data(:,iVar)   )
        end do
      end select

    case ('fpt')

      select case (dim)
      case (3)
        call ply_LegToPnt_3D( fpt       = me%body_3d%fpt, &
          &                   pntVal    = nodal_data,     &
          &                   legCoeffs = modal_data,     &
          &                   nVars     = nVars           )
      case (2)
        call ply_LegToPnt_2D( fpt       = me%body_2d%fpt, &
          &                   pntVal    = nodal_data,     &
          &                   legCoeffs = modal_data,     &
          &                   nVars     = nVars           )
      case (1)
        do iVar = 1,nVars
          call me%body_1d%fpt%LegToPnt(          &
            &    pntVal    = nodal_data(:,iVar), &
            &    legCoeffs = modal_data(:,iVar), &
            &    nIndeps   = 1                   )
        end do
      end select

    case ('fxt')
      select case (dim)
      case (3)
        do iVar = 1,nVars
          call ply_fxt_m2n_3D( fxt = me%body_3d%fxt,            &
            &               modal_data = modal_data(:,iVar),    &
            &               nodal_data = nodal_data(:,iVar),    &
            &              oversamp_degree = me%oversamp_degree )
        end do
      case (2)
        do iVar = 1,nVars
          call ply_fxt_m2n_2D( fxt = me%body_2d%fxt,            &
            &               modal_data = modal_data(:,iVar),    &
            &               nodal_data = nodal_data(:,iVar),    &
            &              oversamp_degree = me%oversamp_degree )
        end do

      case (1)
        do iVar = 1,nVars
          call ply_fxt_m2n_1D( fxt = me%body_1d%fxt,         &
            &               modal_data = modal_data(:,iVar), &
            &               nodal_data = nodal_data(:,iVar)  )
        end do
      end select
    end select

  end subroutine ply_poly_project_m2n_multivar
  ! ------------------------------------------------------------------------ !



  ! **************** NODAL to MODAL ****************************************** !
  ! ------------------------------------------------------------------------ !
  !> Convert nodal values to nDoFs modes.
  subroutine ply_poly_project_n2m_multiVar(me, dim, nVars, nodal_data, &
    &                                      modal_data)
    ! -------------------------------------------------------------------- !
    type(ply_poly_project_type), intent(inout) :: me
    integer, intent(in) :: dim
    integer, intent(in) :: nVars
    real(kind=rk), intent(inout) :: nodal_data(:,:)
    real(kind=rk), intent(inout) :: modal_data(:,:)
    ! -------------------------------------------------------------------- !
    integer :: iVar
    ! -------------------------------------------------------------------- !

    select case(trim(me%kind))
    case ('l2p')
      select case (dim)
      case (1)
        do iVar = 1, nVars
          call ply_l2p_trafo_1D( trafo = me%body_1D%l2p%node2leg, &
            &                    projected = modal_data(:,iVar),  &
            &                    original  = nodal_data(:,iVar)   )
        end do
      case (2)
        do iVar = 1, nVars
          call ply_l2p_trafo_2D( trafo = me%body_2D%l2p%node2leg, &
            &                    projected = modal_data(:,iVar),  &
            &                    original  = nodal_data(:,iVar)   )
        end do
      case (3)
        do iVar = 1, nVars
          call ply_l2p_trafo_3D( trafo = me%body_3D%l2p%node2leg, &
            &                    projected = modal_data(:,iVar),  &
            &                    original  = nodal_data(:,iVar)   )
        end do
      end select

     case ('fpt')
      !projection via fpt
      select case (dim)
      case (3)
        call ply_pntToLeg_3D( fpt       = me%body_3d%fpt, &
          &                   nVars     = nVars,          &
          &                   pntVal    = nodal_data,     &
          &                   legCoeffs = modal_data      )
      case (2)
        call ply_pntToLeg_2D( fpt       = me%body_2d%fpt, &
          &                   nVars     = nVars,          &
          &                   pntVal    = nodal_data,     &
          &                   legCoeffs = modal_data      )
      case (1)
        do iVar = 1,nVars
          call me%body_1d%fpt%pntToLeg(          &
            &    nIndeps   = 1,                  &
            &    pntVal    = nodal_data(:,iVar), &
            &    legCoeffs = modal_data(:,iVar)  )
        end do
      end select

    case ('fxt')
      select case (dim)
      case (3)
        do iVar = 1, nVars
          call ply_fxt_n2m_3D(                          &
            &    fxt              = me%body_3d%fxt,     &
            &    nodal_data       = nodal_data(:,iVar), &
            &    modal_data       = modal_data(:,iVar), &
            &    oversamp_degree  = me%oversamp_degree  )
        end do

      case (2)
        do iVar = 1, nVars
          call ply_fxt_n2m_2D(                          &
            &    fxt              = me%body_2d%fxt,     &
            &    nodal_data       = nodal_data(:,iVar), &
            &    modal_data       = modal_data(:,iVar), &
            &    oversamp_degree  = me%oversamp_degree  )
        end do

      case (1)
        do iVar = 1, nVars
          call ply_fxt_n2m_1D(                          &
            &    fxt              = me%body_1d%fxt,     &
            &    nodal_data       = nodal_data(:,iVar), &
            &    modal_data       = modal_data(:,iVar)  )
        end do
      end select

    case default
       write(logUnit(1),*) 'ERROR in projection nodal to modal'
    end select

  end subroutine ply_poly_project_n2m_multivar
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> function to provide the coordinates from the quadrature points on the faces
  ! idir and ialign are inputs which identify which face is needed
  ! faces is allocated as face(dir,align)
  subroutine ply_get_quadpoints_faces(poly_proj, idir, ialign, points)
    ! -------------------------------------------------------------------- !
    type(ply_poly_project_type), intent(in) :: poly_proj
    integer, intent(in)                     :: idir
    integer, intent(in)                     :: ialign
    real(kind=rk), allocatable, intent(inout) :: points (:,:)
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

     allocate (points(poly_proj%body_3d%faces(idir,iAlign)%nquadpoints,3))
     points = poly_proj%body_3d%faces(idir,iAlign)%points

  end subroutine ply_get_quadpoints_faces
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine ply_get_quadpoints_faces_2d(poly_proj, idir, ialign, points)
    ! -------------------------------------------------------------------- !
    type(ply_poly_project_type), intent(in) :: poly_proj
    integer, intent(in)                     :: idir
    integer, intent(in)                     :: ialign
    real(kind=rk), allocatable, intent(out) :: points (:,:)
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

     allocate (points(poly_proj%body_2d%faces(idir,iAlign)%nquadpoints,3))
     points = poly_proj%body_2d%faces(idir,iAlign)%points

  end subroutine ply_get_quadpoints_faces_2d
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine ply_get_quadpoints_faces_1d(poly_proj, idir, ialign, points)
    ! -------------------------------------------------------------------- !
    type(ply_poly_project_type), intent(in) :: poly_proj
    integer, intent(in)                     :: idir
    integer, intent(in)                     :: ialign
    real(kind=rk), allocatable, intent(out) :: points (:,:)
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
     allocate (points(poly_proj%body_1d%faces(idir,iAlign)%nquadpoints,3))
     points = poly_proj%body_1d%faces(idir,iAlign)%points
  end subroutine ply_get_quadpoints_faces_1d
  ! ------------------------------------------------------------------------ !

end module ply_poly_project_module
