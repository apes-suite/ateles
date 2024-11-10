! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
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
!> author: Jens Zudrop
!! Module provides routines to set boundary values for
!! the co-volume stabilization.
module atl_covolume_boundary_module

  ! Treelm modules
  use env_module,                    only: rk
  use tem_time_module,               only: tem_time_type
  use tem_aux_module,                only: tem_abort
  use tem_faceData_module,           only: tem_invFace_map, &
                                         & tem_left, &
                                         & tem_right
  use tem_logging_module,            only: logUnit
  use treelmesh_module,              only: treelmesh_type

  ! Ateles modules
  use atl_cube_elem_module,          only: atl_cube_elem_type
  use atl_boundary_module,           only: atl_level_boundary_type
  use atl_bc_header_module,          only: atl_boundary_type
  use atl_equation_module,           only: atl_equations_type
  use atl_scheme_module,             only: atl_scheme_type, &
                                         & atl_modg_scheme_prp, &
                                         & atl_modg_2d_scheme_prp, &
                                         & atl_modg_1d_scheme_prp
  use ply_poly_project_module,       only: ply_poly_project_type
  use ply_poly_project_module,       only: ply_poly_project_type, assignment(=)
  use atl_modg_bnd_module,           only: atl_modg_bnd
  use atl_modg_2d_bnd_module,        only: atl_modg_2d_bnd
  use atl_modg_1d_bnd_module,        only: atl_modg_1d_bnd
  use ply_modg_basis_module,         only: ply_faceValLeftBndAns
  use atl_kerneldata_module,         only: atl_statedata_type

  implicit none
  private

  public :: atl_set_covolume_bnd

contains

  !> Routine to set boundary values for the covolume filter.
  subroutine atl_set_covolume_bnd( bc, boundary, state, scheme, equation, &
    & tree, time, mesh, poly_proj, nodalBnd, iDir, iLevel                 )
    ! ---------------------------------------------------------------------------
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary
    !> The face data on the current level
    type(atl_statedata_type), intent(inout) :: state
    !> The parameters of th the modg scheme.
    type(atl_scheme_type), intent(inout) :: scheme
    !> The underlying equation system
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> The description of the mesh on the current level.
    type(atl_cube_elem_type),  intent(in) :: mesh
    !> Data for the projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Set boundaries in nodal fashion by default? If set to false,
    !! the boundaries may still be set in nodal way whenever necessary (e.g.
    !! boundaries which have space-time dependence, etc.)
    logical, intent(in) :: nodalBnd
    !> The spatial direction to set
    integer, intent(in) :: iDir
    !> current Level working on
    integer, intent(in) :: iLevel
    ! ---------------------------------------------------------------------------
    integer :: nBCs, facePos, neighPos, neighAlign, nScalars
    integer :: iBC, iFace, iAlign
    real(kind=rk), allocatable :: faceOp(:,:), stateOp(:,:), faceRep_bnd(:,:)
    real(kind=rk) :: elemLength
    real(kind=rk) :: bndBaryCoord(1:3)
    integer :: nquadpoints, ndofs, ndofsface
    integer :: oversamp_dofs, oversamp_degree, degree
    ! ---------------------------------------------------------------------------

    select case(scheme%scheme)
    case(atl_modg_scheme_prp, atl_modg_2d_scheme_prp, atl_modg_1d_scheme_prp)
    case default
      write(logUnit(1),*) 'ERROR: Not able to set co-volume boundary'
      write(logUnit(1),*) 'values for this scheme, stopping ...'
      call tem_abort()
    end select

    nScalars = equation%varSys%nScalars

    nBCs = boundary%nBCs

    ! initialized just to silence a compiler warning, will be overwritten later
    ! on
    ndofs = 0

    ! The length of an element
    elemLength = mesh%length

    ! get correct amount of quadrature points, dofs, degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    select case(equation%eq_kind)
    case( 'maxwell','maxwelldivcorrection','loclineuler','euler','acoustic', &
      &   'heat', 'navier_stokes', 'filtered_navier_stokes' )
      ndofs = poly_proj%body_3D%ndofs
      nquadpoints = poly_proj%body_2D%nquadpoints
      ndofsface = poly_proj%body_2D%ndofs
      oversamp_dofs = poly_proj%body_2D%oversamp_dofs
      oversamp_degree = poly_proj%oversamp_degree
      degree = scheme%modg%maxPolyDegree
    case('maxwell_2d','euler_2d', 'navier_stokes_2d', 'acoustic_2d', 'heat_2d', &
      &  'filtered_navier_stokes_2d' )
      ndofs = poly_proj%body_2D%ndofs
      nquadpoints = poly_proj%body_1D%nquadpoints
      ndofsface = poly_proj%body_1D%ndofs
      oversamp_dofs = poly_proj%body_1D%oversamp_dofs
      oversamp_degree = poly_proj%oversamp_degree
      degree = scheme%modg_2d%maxPolyDegree
    case('euler_1d','advection_1d','heat_1d')
      ndofs = poly_proj%body_1D%ndofs
      nquadpoints = 1
      ndofsface = 1
      oversamp_dofs = 1
      oversamp_degree = 1
      degree = scheme%modg_1d%maxPolyDegree
    case default
      write(logUnit(1),*) 'ERROR: Unknwon equation for setting '
      write(logUnit(1),*) 'co-volume boundary states: '//trim(equation%eq_kind)
      write(logUnit(1),*) ' stopping ...'
      call tem_abort()
    end select

    ! @todo add other variables to private if necessary
    !!!!OMP PARALLEL &
    !!!!OMP PRIVATE(iBC, iAlign, iFace) &
    !!!!OMP DEFAULT(shared)

    ! Iterate over all the boundaries and set the right face values for
    ! the boundaries on all relevant faces.
    do iBC = 1, nBCs

      ! Now, we iterate over all the faces with this boundary conditions and
      ! set the corresponding face values
      allocate( faceOp(ndofsface, equation%varSys%nScalars) )
      allocate( faceRep_bnd(ndofsface, equation%varSys%nScalars) )
      allocate( stateOp(ndofs, equation%varSys%nScalars) )

      do iAlign = 1,2
        do iFace = 1, boundary%bnd(iBC)%faces(iDir, iAlign)%facePos%nVals

          ! The position of the face with the current boundary condition
          ! inside the face representation.
          facePos = boundary%bnd(iBC)%faces(iDir, iAlign)%facePos%val(iFace)

          ! Create the modal representation on the face for the current
          ! face. We need the modal representation of the neighboring fluid
          ! element for that.
          neighPos = boundary%bnd(iBC)%faces(iDir, iAlign)%neighPos%val(iFace)
          neighAlign = tem_invFace_map(iAlign)
          stateOp = state%state(neighPos,:,:)

          ! get the barycentric coordinate of the (virtual) boundary element
          bndBaryCoord(1:3) = mesh%bary_coord(neighPos,1:3)
          bndBaryCoord(iDir) = bndBaryCoord(iDir) + ((-1.0_rk)**neighAlign)*elemLength

          ! Evaluate state at the face
          select case(equation%nDimensions)
          !! ----------------------------------------------------------------------------
          !! MODG
          !! ----------------------------------------------------------------------------
          case(3)
            ! Evaluate neighboring element on the face.
            faceOp = atl_volToFace( state = stateOp, iDir = iDir, iAlign = neighAlign, &
                                  & nDofsFace = ndofsface, maxPolyDeg = degree, &
                                  & nScalars = equation%varSys%nScalars )
            ! Set boundary values on the face.
            call atl_modg_bnd( bc              = bc(iBC),                    &
              &                faceOp          = faceOp,                     &
              &                poly_proj       = poly_proj,                  &
              &                normalRot       = equation%varRotation(iDir), &
              &                nDerivatives    = 0,                          &
              &                equation        = equation,                   &
              &                isNodalScheme   = nodalBnd ,                  &
              &                time            = time,                       &
              &                currentFace     = iFace,                      &
              &                currentLevel    = iLevel,                     &
              &                nquadpoints     = nquadpoints,                &
              &                ndofs           = ndofsface ,                 &
              &                oversamp_dofs   = oversamp_dofs,              &
              &                modalFace       = faceRep_bnd(:,:)            )

            ! Lift the face information to the volume.
            state%state(facePos,:,:) = &
                & atl_extend_covol_face( faceState = faceRep_bnd, &
                                       & iDir = iDir, &
                                       & maxPolyDeg = degree, &
                                       & nScalars = equation%varSys%nScalars )

          !! ----------------------------------------------------------------------------
          !! MODG_2D
          !! ----------------------------------------------------------------------------
          case(2)
            ! Evaluate neighboring element on the face.
            faceOp = atl_volToFace_2d(state = stateOp, iDir = iDir, iAlign = neighAlign, &
                                 & nDofsFace = ndofsface, maxPolyDeg = degree, &
                                 & nScalars = equation%varSys%nScalars )
            ! Set boundary values on the face.

            call atl_modg_2d_bnd( bc            = bc(iBC),                    &
              &                   faceOp        = faceOp,                     &
              &                   poly_proj     = poly_proj,                  &
              &                   normalRot     = equation%varRotation(iDir), &
              &                   equation      = equation,                   &
              &                   nDerivatives  = 0,                          &
              &                   isNodalScheme = nodalBnd ,                  &
              &                   time          = time,                       &
              &                   currentFace   = iFace,                      &
              &                   currentLevel  = iLevel,                     &
              &                   nquadpoints   = nquadpoints,                &
              &                   ndofs         = ndofsface ,                 &
              &                   oversamp_dofs = oversamp_dofs,              &
              &                   modalFace     = faceRep_bnd(:,:)            )

            ! Lift the face information to the volume.
            state%state(facePos,:,:) = &
                   & atl_extend_covol_face_2d(faceState = faceRep_bnd, &
                                          & iDir = iDir, &
                                          & maxPolyDeg = degree, &
                                          & nScalars = equation%varSys%nScalars )

           !! ----------------------------------------------------------------------------
           !! MODG_1D
           !! ----------------------------------------------------------------------------
           case(1)
             ! Evaluate neighboring element on the face.
             faceOp = atl_volToFace_1d(state = stateOp, iAlign = neighAlign, &
               &        nDofsFace = ndofsface, maxPolyDeg = degree,          &
               &        nScalars = equation%varSys%nScalars                  )
             ! Set boundary values on the face.
             faceRep_bnd(:,:) = atl_modg_1d_bnd(             &
               & bc            = bc(iBC),                    &
               & faceOp        = faceOp,                     &
               & poly_proj     = poly_proj,                  &
               & normalRot     = equation%varRotation(iDir), &
               & equation      = equation,                   &
               & tree          = tree,                       &
               & isNodalScheme = nodalBnd ,                  &
               & time          = time,                       &
               & faceDir       = iDir,                       &
               & leftOrRight   = neighAlign,                 &
               & bndBaryCoord  = bndBaryCoord,               &
               & elemLength    = elemLength                  )
            ! Lift the face information to the volume.
            state%state(facePos,:,:) = &
                   & atl_extend_covol_face_1d(faceState = faceRep_bnd, &
                                          & maxPolyDeg = degree, &
                                          & nScalars = equation%varSys%nScalars )
          case default
            write(logUnit(1), *) 'ERROR: Unknown spatial dimension!'
            write(logUnit(1), *) 'Not able to set co-volume boundaries, stopping ...'
            call tem_abort()
          end select



        end do ! iFace
      end do ! iAlign
      deallocate(faceOp)
      deallocate(stateOp)
      deallocate(faceRep_bnd)
    end do ! iBC

    !!!!OMP END PARALLEL

  end subroutine atl_set_covolume_bnd

  !> Lift face data to volume data.
  function atl_extend_covol_face(faceState, iDir, maxPolyDeg, nScalars) &
                                & result( volState )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: faceState(:,:)
    !> The spatial direction:
    !! 1 -> x direction.
    !! 2 -> y direction.
    !! 3 -> z direction.
    integer, intent(in) :: iDir
    integer, intent(in) :: maxPolyDeg, nScalars
    real(kind=rk) :: volState((maxPolyDeg+1)**3, nScalars)
    ! --------------------------------------------------------------------------
    integer :: iAnsX, iAnsY, iAnsZ, pos, pos_face
    ! --------------------------------------------------------------------------

    ! We inititalize the lifted volumetric data with 0.
    volState(:,:) = 0.0_rk

    select case(iDir)
    ! Extend x-direction (i.e. 2D facewise modal representation in y-z-direction
    ! is available).
    case(1)
      do iAnsZ = 1, maxPolyDeg+1
        do iAnsY = 1, maxPolyDeg+1
  pos = 1                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansy                                      &
    &      + ( ( iansz-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
          volState(pos,:) = faceState(pos_face,:)
        end do
      end do

    ! Extend y-direction (i.e. 2D facewise modal representation in x-z-direction
    ! is available).
    case(2)
      do iAnsX = 1, maxPolyDeg+1
        do iAnsZ = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansx                                      &
    &      + ( ( iansz-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
          volState(pos,:) = faceState(pos_face,:)
        end do
      end do

    ! Extend z-direction (i.e. 2D facewise modal representation in x-y-direction
    ! is available).
    case(3)
      do iAnsX = 1, maxPolyDeg+1
        do iAnsY = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
          volState(pos,:) = faceState(pos_face,:)
        end do
      end do

    case default
      write(logUnit(1),*) 'ERROR in atl_extend_covol_face: '
      write(logUnit(1),*) 'Unknown direction, stopping ...'
      stop
    end select

  end function atl_extend_covol_face


  !> Lift face data to volume data.
  function atl_extend_covol_face_2d(faceState, iDir, maxPolyDeg, nScalars) &
                                & result( volState )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: faceState(:,:)
    !> The spatial direction:
    !! 1 -> x direction.
    !! 2 -> y direction.
    !! 3 -> z direction.
    integer, intent(in) :: iDir
    integer, intent(in) :: maxPolyDeg, nScalars
    real(kind=rk) :: volState((maxPolyDeg+1)**2, nScalars)
    ! --------------------------------------------------------------------------
    integer :: iAnsX, iAnsY, pos
    ! --------------------------------------------------------------------------

    ! We inititalize the lifted volumetric data with 0.
    volState(:,:) = 0.0_rk

    select case(iDir)
    ! Extend x-direction (i.e. 1D facewise modal representation in y-direction
    ! is available).
    case(1)
      do iAnsY = 1, maxPolyDeg+1
  pos = 1                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
        volState(pos,:) = faceState(iAnsY,:)
      end do


    ! Extend y-direction (i.e. 1D facewise modal representation in x-direction
    ! is available).
    case(2)
      do iAnsX = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
        volState(pos,:) = faceState(iAnsX,:)
      end do
    case default
      write(logUnit(1),*) 'ERROR in atl_extend_covol_face_2d: '
      write(logUnit(1),*) 'Unknown direction, stopping ...'
      stop
    end select

  end function atl_extend_covol_face_2d


  !> Lift face data to volume data.
  function atl_extend_covol_face_1d(faceState, maxPolyDeg, nScalars) &
                                & result( volState )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: faceState(:,:)
    integer, intent(in) :: maxPolyDeg, nScalars
    real(kind=rk) :: volState((maxPolyDeg+1), nScalars)
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    volState(:,:) = 0.0_rk
    volState(1,:) = faceState(1,:)

  end function atl_extend_covol_face_1d


  !> Project elemental state to a particular face.
  function atl_volToFace(state, iDir, iAlign, maxPolyDeg, nDofsFace, nScalars) &
                        & result(face)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: state(:,:)
    integer, intent(in) :: iDir, iAlign, maxPolyDeg, nDofsFace, nScalars
    real(kind=rk) :: face(nDofsFace,nScalars)
    ! --------------------------------------------------------------------------
    integer :: pos, pos_face, iAnsX, iAnsY, iAnsZ
    real(kind=rk) :: faceVal
    ! --------------------------------------------------------------------------

    face(:,:) = 0.0_rk

    select case(iDir)
    ! Evaluate in x direction
    case(1)
      if(iAlign == tem_left) then ! Eval on ref. elem. at x=-1
        do iAnsZ = 1, maxPolyDeg+1
          do iAnsY = 1, maxPolyDeg+1
            do iAnsX = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansy                                      &
    &      + ( ( iansz-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
              faceVal = ply_faceValLeftBndAns(iAnsX)
              face(pos_face,:) = face(pos_face,:) + faceVal * state(pos,:)
            end do
          end do
        end do
      elseif(iAlign == tem_right) then ! Eval on ref. elem. at x=+1
        do iAnsZ = 1, maxPolyDeg+1
          do iAnsY = 1, maxPolyDeg+1
            do iAnsX = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansy                                      &
    &      + ( ( iansz-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
              face(pos_face,:) = face(pos_face,:) + state(pos,:)
            end do
          end do
        end do
      else
        write(logUnit(1),*) 'ERROR in atl_volToFace: unknown face'
        write(logUnit(1),*) 'alignment, stopping ...'
        stop
      end if
     case(2)
       if(iAlign == tem_left) then ! Eval on ref. elem. at y=-1
         do iAnsZ = 1, maxPolyDeg+1
           do iAnsX = 1, maxPolyDeg+1
             do iAnsY = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansx                                      &
    &      + ( ( iansz-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
               faceVal = ply_faceValLeftBndAns(iAnsY)
               face(pos_face,:) = face(pos_face,:) + faceVal * state(pos,:)
             end do
           end do
         end do
       elseif(iAlign == tem_right) then ! Eval on ref. elem. at y=+1
         do iAnsZ = 1, maxPolyDeg+1
           do iAnsX = 1, maxPolyDeg+1
             do iAnsY = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansx                                      &
    &      + ( ( iansz-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
               face(pos_face,:) = face(pos_face,:) + state(pos,:)
             end do
           end do
         end do
       else
         write(logUnit(1),*) 'ERROR in atl_volToFace: unknown face'
         write(logUnit(1),*) 'alignment, stopping ...'
         stop
       end if
     case(3)
       if(iAlign == tem_left) then ! Eval on ref. elem. at z=-1
         do iAnsY = 1, maxPolyDeg+1
           do iAnsX = 1, maxPolyDeg+1
             do iAnsZ = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
               faceVal = ply_faceValLeftBndAns(iAnsZ)
               face(pos_face,:) = face(pos_face,:) + faceVal * state(pos,:)
             end do
           end do
         end do
        elseif(iAlign == tem_right) then ! Eval on ref. elem. at z=+1
          do iAnsY = 1, maxPolyDeg+1
            do iAnsX = 1, maxPolyDeg+1
              do iAnsZ = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (iansz-1)*(maxpolydeg+1))*(maxpolydeg+1)
  pos_face = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
                face(pos_face,:) = face(pos_face,:) + state(pos,:)
              end do
            end do
          end do
        else
          write(logUnit(1),*) 'ERROR in atl_volToFace: unknown face'
          write(logUnit(1),*) 'alignment, stopping ...'
          stop
        end if
    case default
      write(logUnit(1),*) 'ERROR in atl_volToFace: unknown direction, stopping ...'
      stop
    end select


  end function atl_volToFace


  !> Project elemental state to a particular face.
  function atl_volToFace_2d(state, iDir, iAlign, maxPolyDeg, nDofsFace, nScalars) &
                        & result(face)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: state(:,:)
    integer, intent(in) :: iDir, iAlign, maxPolyDeg, nDofsFace, nScalars
    real(kind=rk) :: face(nDofsFace,nScalars)
    ! --------------------------------------------------------------------------
    integer :: pos, iAnsX, iAnsY
    real(kind=rk) :: faceVal
    ! --------------------------------------------------------------------------

    face(:,:) = 0.0_rk

    select case(iDir)
    ! Evaluate in x direction
    case(1)
      if(iAlign == tem_left) then ! Eval on ref. elem. at x=-1
        do iAnsY = 1, maxPolyDeg+1
          do iAnsX = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
            faceVal = ply_faceValLeftBndAns(iAnsX)
            face(iAnsY,:) = face(iAnsY,:) + faceVal * state(pos,:)
          end do
        end do
      elseif(iAlign == tem_right) then ! Eval on ref. elem. at x=+1
        do iAnsY = 1, maxPolyDeg+1
          do iAnsX = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
            face(iAnsY,:) = face(iAnsY,:) + state(pos,:)
          end do
        end do
      else
        write(logUnit(1),*) 'ERROR in atl_volToFace_2d: unknown face'
        write(logUnit(1),*) 'alignment, stopping ...'
        stop
      end if
    case(2)
      if(iAlign == tem_left) then ! Eval on ref. elem. at y=-1
        do iAnsX = 1, maxPolyDeg+1
          do iAnsY = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
            faceVal = ply_faceValLeftBndAns(iAnsY)
            face(iAnsX,:) = face(iAnsX,:) + faceVal * state(pos,:)
          end do
        end do
      elseif(iAlign == tem_right) then ! Eval on ref. elem. at y=+1
        do iAnsX = 1, maxPolyDeg+1
          do iAnsY = 1, maxPolyDeg+1
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydeg+1))*(maxpolydeg+1)
            face(iAnsX,:) = face(iAnsX,:) + state(pos,:)
          end do
        end do
      else
        write(logUnit(1),*) 'ERROR in atl_volToFace_2d: unknown face'
        write(logUnit(1),*) 'alignment, stopping ...'
        stop
      end if
    case default
      write(logUnit(1),*) 'ERROR in atl_volToFace_2d: unknown direction, stopping ...'
      stop
    end select

  end function atl_volToFace_2d


  !> Project elemental state to a particular face.
  function atl_volToFace_1d(state, iAlign, maxPolyDeg, nDofsFace, nScalars) &
    &        result(face)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: state(:,:)
    !> The face alignment.
    !! 1 -> left face (ref. element at -1)
    !! 2 -> right face (ref. element at +1)
    integer, intent(in) :: iAlign
    integer, intent(in) :: maxPolyDeg, nDofsFace, nScalars
    real(kind=rk) :: face(nDofsFace,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: faceVal
    integer :: iAnsX
    ! --------------------------------------------------------------------------

    face(:,:) = 0.0_rk

    if(iAlign == tem_left) then ! Eval on ref. elem. at x=-1
      do iAnsX = 1, maxPolyDeg+1
        faceVal = ply_faceValLeftBndAns(iAnsX)
        face(1,:) = face(1,:) + faceVal * state(iAnsX,:)
      end do
    elseif(iAlign == tem_right) then ! Eval on ref. elem. at x=+1
      do iAnsX = 1, maxPolyDeg+1
        face(1,:) = face(1,:) + state(iAnsX,:)
      end do
    else
      write(logUnit(1),*) 'ERROR in atl_volToFace_1d: unknown face'
      write(logUnit(1),*) 'alignment, stopping ...'
      stop
    end if


  end function atl_volToFace_1d

end module atl_covolume_boundary_module
