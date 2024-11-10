! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> author: Jens Zudrop
!! Module collects all routines, datatypes, etc. to set boundary
!! contitions for the MOdal Discontinuous Galerkin scheme.
module atl_modg_bnd_module
  use, intrinsic :: iso_c_binding,    only: c_f_pointer

  ! Treelm modules
  use env_module,                     only: rk
  use tem_aux_module,                 only: tem_abort
  use tem_faceData_module,            only: tem_invFace_map
  use tem_coordinate_module,          only: coordRotation_type
  use tem_bc_module,                  only: tem_bc_state_type
  use tem_time_module,                only: tem_time_type
  use tem_logging_module,             only: logUnit
  use tem_timer_module,               only: tem_startTimer, &
    &                                       tem_stopTimer
  use tem_varSys_module,              only: tem_varSys_type
  use tem_spacetime_fun_module,       only: tem_st_fun_listElem_type

  ! Ateles modules
  use atl_kerneldata_module,          only: atl_statedata_type
  use atl_bc_header_module,           only: atl_boundary_type
  use atl_boundary_module,            only: atl_level_boundary_type
  use atl_facedata_module,            only: atl_facedata_type
  use atl_equation_module,            only: atl_equations_type
  use atl_cube_elem_module,           only: atl_cube_elem_type
  use ply_poly_project_module,        only: ply_poly_project_type, &
    &                                       assignment(=),         &
    &                                       ply_poly_project_m2n,  &
    &                                       ply_poly_project_n2m
  use atl_materialPrp_module,         only: atl_faceMaterialData_type
  use atl_timer_module,               only: atl_timerHandles
  use ply_oversample_module,          only: ply_convert2oversample, &
    &                                       ply_convertFromOversample

  implicit none

  private

  public :: atl_modg_set_bnd, atl_modg_bnd

contains

  ! ***************************************************************************
  !> Subroutine to set face values to impose boundary conditions
  !! at a certain point of the domain. The subroutine is operating levelwise.
  !! @todo: For zero gradient BC, we need stateData passed into this routine
  subroutine atl_modg_set_bnd( bc, boundary, facedata, equation,  time, mesh, &
    &                          poly_proj, nodalBnd, material, currentLevel,   &
    &                          statedata                                      )
    ! ---------------------------------------------------------------------------
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary
    !> The face data on the current level
    type(atl_facedata_type), intent(inout) :: facedata
    !> The state data on the current level
    type(atl_statedata_type), intent(inout) :: statedata
    !> The underlying equation system
    type(atl_equations_type), intent(in) :: equation
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
    !> Material description of the faces contained in boundary. One
    !! entry for each spatial direction, i.e. x,y,z.
    type(atl_faceMaterialData_type), intent(in), optional :: material(3)
    !> the level to compute on
    integer, intent(in) :: currentLevel
    ! ---------------------------------------------------------------------------
    integer :: facePos, neighPos, neighAlign, nscalars
    integer :: iBC, iDir, iFace, iAlign, iglobFace
    real(kind=rk), allocatable :: faceOp(:,:)
    real(kind=rk) :: elemLength
    real(kind=rk) :: bndBaryCoord(1:3)
    real(kind=rk), allocatable :: faceMaterial(:,:)
    integer :: nQuadPoints, nDofs
    integer :: oversamp_dofs, oversamp_degree
    ! varaibles for zero gradient BC
    integer :: signFact, iMode, jMode, nModes, modePos, offset_face, offset_vol
    integer :: kMode, k, iStride, mfacepos
    real(kind=rk) :: corrector(poly_proj%body_2D%nDofs,size(statedata%state,3))
    real(kind=rk) :: sidefact
    ! ---------------------------------------------------------------------------

    if (present(material)) then
      allocate( faceMaterial( size(material(1)%leftElemMaterialDat,2), &
        &                     size(material(1)%leftElemMaterialDat,3)) )
    end if

    nScalars = equation%varSys%nScalars

    ! The length of an element
    elemLength = mesh%length

    ! get correct amount of quadrature points, dofs, degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nQuadPoints = poly_proj%body_2D%nQuadPoints
    nDofs = poly_proj%body_1D%nDofs
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs
    oversamp_degree = poly_proj%oversamp_degree

    ! Iterate over all the boundaries and set the right face values for
    ! the boundaries on all relevant faces.

    do iBC = 1, boundary%nBCs

      ! Now, we iterate over all the faces with this boundary conditions and
      ! set the corresponding face values
      iglobFace = 0

      ! calculate the highest mode which will be modified for zero grd BC
       nModes = ceiling(nDofs * bc(iBC)%neumann_mode_fraction)
      do iDir = 1,3
        allocate( faceOp(size(facedata%faceRep(iDir)%dat,2), &
          &              size(facedata%faceRep(iDir)%dat,3)) )
        do iAlign = 1,2

          neighAlign = tem_invFace_map(iAlign)
          ! for zero gradient BC calculation
           sidefact = (-1)**neighAlign

          do iFace = 1, boundary%bnd(iBC)%faces(iDir, iAlign)%facePos%nVals
            ! count global iterator
            iglobFace = iglobFace +1

            ! The position of the face with the current boundary condition
            ! inside the face representation.
            facePos = boundary%bnd(iBC)%faces(iDir, iAlign)%facePos%val(iFace)

            ! Create the modal representation on the face for the current
            ! face. We need the modal representation of the neighboring fluid
            ! element for that.
            neighPos = boundary%bnd(iBC)%faces(iDir, iAlign)%neighPos%val(iFace)

            if (bc(iBC)%enforce_zero_grad) then
            !if ( .false. ) then

               ! faceOp is the state that we tell the boundary condition is
               ! present at the face.
               ! Initially set this to the first mode as this will be used if only
               ! one mode is to be used for the extrapolation.
               if ( iDir == 1 ) then
                 ! mode coefficients: a_{1,j,k}, j,k = 1,...,nDoFs
                 faceOp = statedata%state(neighPos,1::nDoFs,:)
               else if ( iDir == 2 ) then
                 ! mode coefficients: a_{i,1,k}, i,k = 1,...,nDoFs
                do k = 1, nDoFs
                   offset_face = (k-1)*nDoFs
                   offset_vol  = (k-1)*nDoFs*nDoFs
                   faceOp(offset_face+1:offset_face+nDofs,:) = &
                     & statedata%state(neighPos,offset_vol+1:offset_vol+nDoFs,:)
                 end do
               else ! iDir == 3
                 ! mode coefficients: a_{i,j,1}, i,j = 1,...,nDoFs
                 faceOp(:,:) = statedata%state(neighPos,1:nDoFs*nDoFs,:)
               end if

               if (nModes > 1) then
                 ! Only need to actually compute something if we are to use more
                 ! then just the first mode for the extrapolation
                 ! signfact is used to keep track of the alternating sign on the
                 ! left end of the element.
                 ! (sidefact decides whether we are left (=-1) or right(=1))

                 ! To enforce a zero gradient we compute the last mode
                 ! (at nModes) to match all previous modes contributions in the
                 ! derivative. This "correcting" mode is stored in corrector.
                 corrector = 0.0_rk

                 ! For x direction, jMode
                 ! For y direction, iMode
                 ! For z dircetion, kMode
                 do kMode = 1, nModes
                   do jMode = 1, nModes
                       mfacepos = jMode + (kMode-1)*nModes
                     ! By taking the derivative all modes are shifted down by 1.
                     ! Thus, the relevant modes are those up to nModes-1. However,
                     ! the one we want to set to enforce the 0 gradient, is one in
                     ! nModes-1. Thus, we need to sum up to nModes-2...
                     ! We set signfact before the iMode loop starts, to have the
                     ! right sign for the computation, which depends on the iMode
                     signfact = int(sidefact)
                     do iMode=1,nModes-2

                       ! mode position in state array
                       if ( iDir == 1 ) then
                         modePos = (iMode+1) + ((jMode-1) + (kMode-1)*nModes)*nModes
                         iStride =1
                       else if ( iDir == 2 ) then
                         modePos = jMode + ((iMode) + (kMode-1)*nModes)*nModes
                         iStride = nModes
                       else !iDir==3
                         modePos = jMode + ((kMode-1) + (iMode)*nModes)*nModes
                         iStride = nModes**2
                       end if

                       ! At the same time we can compute the value at the face by
                       ! summing the modes.
                       faceOp(mfacepos,:nscalars) = faceOp(mfacepos,:nscalars)    &
                         &         + signfact * statedata%state(neighpos,modePos, &
                         &                                               :nscalars)
                       ! The derivative contribution of iMode+1, we subtract it from
                       ! the corrector to obtain a corrector, which will balance out
                       ! all other modes.
                       ! do it for kMode and jMode over all x

                       corrector(mfacepos,:) = corrector(mfacepos,:)   &
                         &       + (signfact*sidefact*(iMode*(iMode+1))/2)      &
                         &       * statedata%state(neighpos,modePos,:)

                       signfact = int(signfact * sidefact)
                     end do ! iMode


                   ! The last mode is not used, instead we use the corrector as
                   ! last mode. As that was computed in the derivative, it needs to
                   ! "integrated" (divison by the factor for the Legendre derivative)
                   faceOp(mfacepos,:nscalars) = faceOp(mfacepos,:nscalars) &
                     &                         - (sidefact*2*corrector(mfacepos,:)) &
                     &                                 / ((nModes)*(nModes-1))

                 end do ! jMode
               end do !kMode
               end if ! nModes > 1

            else ! zero gradient BC is not available yet, so it always goes to here
              faceOp = facedata%faceRep(iDir)%dat(neighPos,:,:,neighAlign)
            end if

            ! get the barycentric coordinate of the (virtual) boundary element
            bndBaryCoord(1:3) = mesh%bary_coord(neighPos,1:3)
            bndBaryCoord(iDir) = bndBaryCoord(iDir) &
              & + ( (-1.0_rk)**neighAlign ) * elemLength

            ! If material information is passed to this function,
            ! we read out the material of the boundary face and pass it
            ! to the atl_modg_bnd function. If not, we call atl_modg_bnd without
            ! this information.
            if( present(material) ) then
              ! We read out the left face material (at a boundary face left
              ! and right face material should be equal).
              faceMaterial(:,:) = material(iDir)%leftElemMaterialDat(iFace,:,:)
              ! @todo JZ: replace by subroutine call, due to OpenMP
              call atl_modg_bnd(                                 &
                & bc              = bc(iBC),                     &
                & faceOp          = faceOp,                      &
                & poly_proj       = poly_proj,                   &
                & normalRot       = equation%varRotation(iDir),  &
                & nDerivatives    = equation%nDerivatives,       &
                & equation        = equation,                    &
                & isNodalScheme   = nodalBnd,                    &
                & time            = time,                        &
                & currentFace     = iglobFace,                   &
                & currentLevel    = currentLevel,                &
                & nQuadPoints     = nQuadPoints,                 &
                & nDofs           = poly_proj%body_2D%nDofs,     &
                & oversamp_dofs   = oversamp_dofs,               &
                & modalFace       = facedata%faceRep(iDir)       &
                  &                  %dat(facePos,:,:,iAlign),   &
                & faceMaterial    = faceMaterial                 )

            else

              call  atl_modg_bnd(                                      &
                &  bc              = bc(iBC),                          &
                &  faceOp          = faceOp,                           &
                &  poly_proj       = poly_proj,                        &
                &  normalRot       = equation%varRotation(iDir),       &
                &  nDerivatives    = equation%nDerivatives,            &
                &  equation        = equation,                         &
                &  isNodalScheme   = nodalBnd,                         &
                &  time            = time,                             &
                &  currentFace     = iglobFace,                        &
                &  currentLevel    = currentLevel,                     &
                &  nquadpoints     = nquadpoints,                      &
                &  ndofs           = poly_proj%body_2D%nDofs,          &
                &  oversamp_dofs   = oversamp_dofs,                    &
                &  modalFace       =  facedata%faceRep(iDir)           &
                  &                           %dat(facePos,:,:,iAlign) )
            end if! material
          end do ! iFace
        end do ! iAlign
        deallocate(faceOp)
      end do ! iDir
    end do ! iBC

  end subroutine atl_modg_set_bnd
  ! ***************************************************************************


  ! ***************************************************************************
  !> Subroutine to create the modal representation for a ceratin boundary face.
  subroutine atl_modg_bnd( bc, faceOp, poly_proj, nDerivatives, equation,  &
    &                      normalRot, isNodalScheme, time, currentFace,    &
    &                      currentLevel, nQuadPoints, nDofs,oversamp_dofs, &
    &                      modalFace,faceMaterial                          )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> Data for the projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> Rotation indices to rotate global coordinate system into face normal
    !! coordinate system.
    type(coordRotation_type), intent(in) :: normalRot
    !> The number of derivative boundaries to be set
    integer, intent(in) :: nDerivatives
    !> Does the solver require isNodalScheme information anyway?
    logical, intent(in) :: isNodalScheme
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> current face to compute on, used for index array
    integer, intent(in) :: currentFace
    !> the level to compute on
    integer, intent(in) :: currentLevel
    !> Number of quadurature points on the face
    integer, intent(in) :: nQuadPoints, ndofs, oversamp_dofs
    !> The modal representation on the boundary face.
    real(kind=rk), intent(inout) :: modalFace(:,:)
    !> The material of the boundary face.
    !! First dimension is the number of points on the face.
    !! Second dimension is the number of material parameters.
    real(kind=rk), intent(in), optional :: faceMaterial(:,:)
    ! ---------------------------------------------------------------------------

    if(isNodalScheme) then
      ! For nonlinear, nodal schemes the boundary conditions will be imposed in a
      ! pointwise way.
      call  modg_nodal_bnd( bc              = bc,              &
        &                   faceOp          = faceOp,          &
        &                   poly_proj       = poly_proj,       &
        &                   equation        = equation,        &
        &                   normalRot       = normalRot,       &
        &                   time            = time,            &
        &                   currentFace     = currentFace,     &
        &                   currentLevel    = currentLevel,    &
        &                   nDerivatives    = nDerivatives,    &
        &                   nQuadPoints     = nQuadPoints,     &
        &                   oversamp_dofs   = oversamp_dofs,   &
        &                   modalFace       = modalFace,       &
        &                   faceMaterial    = faceMaterial     )

    else
      ! For linear equations it could be possible to impose boundary conditions
      ! in modal space directly.
      call modg_modal_bnd( bc            = bc,            &
        &                  faceOp        = faceOp,        &
        &                  poly_proj     = poly_proj,     &
        &                  equation      = equation,      &
        &                  normalRot     = normalRot,     &
        &                  time          = time,          &
        &                  currentFace   = currentFace,   &
        &                  currentLevel  = currentLevel,  &
        &                  nQuadPoints   = nQuadPoints,   &
        &                  nDofs         = nDofs,         &
        &                  oversamp_dofs = oversamp_dofs, &
        &                  modalFace     = modalFace,     &
        &                  faceMaterial  = faceMaterial   )

      if( equation%nDerivatives > 0 ) then
        write(logUnit(1),*) 'ERROR in atl_modg_bnd: not able to ' // &
          & 'boundary conditions for higher order equations in modal way,' // &
          & 'stopping ...'
        call tem_abort()
      end if

    end if

  end subroutine atl_modg_bnd
  ! ***************************************************************************


  ! ***************************************************************************
  !> Set boundary values in a nodal way
  subroutine modg_nodal_bnd( bc, faceOp, poly_proj, equation, normalRot,    &
    &                        time, currentFace, currentLevel, nDerivatives, &
    &                        nQuadPoints, oversamp_dofs, modalFace,         &
    &                        faceMaterial                                   )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> Data for the projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> Rotation indices to rotate global coordinate system into face normal
    !! coordinate system.
    type(coordRotation_type), intent(in) :: normalRot
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> current face to compute on, used for index
    integer, intent(in) :: currentFace
    !> the level to compute on
    integer, intent(in) :: currentLevel
    !> The number of derivative boundaries to be set
    integer, intent(in) :: nDerivatives
    !> Number of quadurature points on the face andi Number of Dofs for the face
    integer, intent(in) :: nQuadPoints, oversamp_dofs
    !> result of the bnd routine, modal coefficent on the boundary faces
    real(kind=rk), intent(out) :: modalFace(:,:)
    !> The material of the boundary face
    !! First dimension is the number of points on the face.
    !! Second dimension is the number of material parameters.
    real(kind=rk), intent(in), optional :: faceMaterial(:,:)
    ! ---------------------------------------------------------------------------
    integer :: iBcVar, bcIndex, nVars
    ! for the oversampling
    ! tmp arrays
    real(kind=rk), allocatable :: pointValOp(:,:), pointFace(:,:)
    real(kind=rk), allocatable :: oversamp_Face(:,:)
    real(kind=rk), allocatable :: pointValOp_derX(:,:), pointValOp_derY(:,:), &
      &                           pointFace_derX(:,:), pointFace_derY(:,:)
    ! ---------------------------------------------------------------------------

    nVars = equation%varSys%nScalars+3*nDerivatives*equation%varSys%nScalars

    allocate(pointFace(nQuadPoints,  nVars))
    allocate(pointValOp(nQuadPoints, nVars))
    allocate(oversamp_face(oversamp_dofs, nVars))
    pointFace = 0.0_rk
    pointValOp = 0.0_rk
    oversamp_face = 0.0_rk

    ! --> modal space
    call ply_convert2oversample(                           &
      & state       = faceOp( :poly_proj%body_2d%min_dofs, &
      &                       :equation%varSys%nScalars ), &
      & poly_proj   = poly_proj,                           &
      & nDim        = 2,                                   &
      & modalcoeffs = oversamp_face                        )
    ! --> oversamp modal space
    ! Now, we transform the modal representation of this element to nodal
    call ply_poly_project_m2n(me         = poly_proj,     &
      &                       dim        = 2 ,            &
      &                       nVars      = nVars,         &
      &                       nodal_data = pointValOp,    &
      &                       modal_data = oversamp_face  )
    ! --> oversamp nodal space


    ! take care of bc_trafo, i.e. transform from conservative to primitive.
    if(.not.bc%bc_trafo%identity) then
      call bc%bc_trafo%to( equation = equation,    &
        &                  instate  = pointValOp,  &
        &                  material = faceMaterial )
    end if

    ! Loop over all the STATE quantitites of the equation system and generate the
    ! modal representation for it.
    do iBcVar = 1, equation%varSys%nScalars
      ! Is a rotation to face normal representation necessary?
      ! So, we get the index for the normal direction on the face.
      if(bc%bc_normal_vec) then
        bcIndex = normalRot%varTransformIndices(iBcVar)
      else
        bcIndex = iBcVar
      end if

      ! Create the modal representation on the face.
      select case(bc%state(iBcVar)%style)
      case('neumann')
        ! In case of dirichlet boundary condition we would like
        ! to prescribe the derivative of the boundary value, so:
        ! Extrapolate the values
        ! @todo replace by subroutine call, due to OpenMP
        pointFace(:,bcIndex) = modg_bnd_extrapolate(            &
          &                      nVals  = nQuadPoints,          &
          &                      faceOp = pointValOp(:,bcIndex) )
      case('dirichlet')
        ! In case of dirichlet boundary condition we would like
        ! to prescribe the boundary value, so:
        ! Mirror the values around the prescribed one.
        ! @todo replace by subroutine call, due to OpenMP

        ! If boundary variable is refered to zero_const then
        ! set modalFace value to zero
        if (bc%varDict%val(iBcVar)%value == 'zero_const') then
          pointFace(:, bcIndex) = 0.0_rk
        else
          pointFace(:,bcIndex) = modg_bnd_mirrorPoint(                   &
            &                      bc           = bc%state(iBcVar),      &
            &                      varSys       = equation%varSys,       &
            &                      time         = time,                  &
            &                      currentFace  = currentFace,           &
            &                      currentLevel = currentLevel,          &
            &                      nPoints      = nQuadPoints            )
        end if
      case default
        write(logUnit(1),*) 'ERROR in modg_nodalBndRep:'
        write(logUnit(1),*)  ' BC style for ', iBcVar, ' is ', &
          &                    bc%state(iBcVar)%style
        write(logUnit(1),*)  ' unknown bnd style, stopping ...'
        call tem_abort()
      end select
    end do !iBcVar
    ! --> pointFace is filled from 1:nScalar with nodal state values on face

    ! check if there are gradients to be set
    if (nDerivatives > 0) then
      ! For gradients temp array
      allocate(pointFace_derX(nQuadPoints, equation%varSys%nScalars))
      allocate(pointFace_derY(nQuadPoints, equation%varSys%nScalars))
      allocate(pointValOp_derX(nQuadPoints,equation%varSys%nScalars))
      allocate(pointValOp_derY(nQuadPoints,equation%varSys%nScalars))

      ! Loop over all the quantitites of the equation system and generate the
      ! modal representation for it.
      do iBcVar = 1, equation%varSys%nScalars
        ! get the nodal gradient from the already transformed array
        pointValOp_derX(:poly_proj%body_2d%min_dofs,iBcVar) = &
          & pointValOp( :poly_proj%body_2d%min_dofs,          &
          &             iBcVar+equation%varSys%nScalars       )
        pointValOp_derY(:poly_proj%body_2d%min_dofs,iBcVar) = &
          & pointValOp( :poly_proj%body_2d%min_dofs,          &
          &             iBcVar+3*equation%varSys%nScalars     )

        ! Is a rotation to face normal representation necessary?
        ! So, we get the index for the normal direction on the face.
        if(bc%bc_normal_vec_gradient) then
          !@todo bcIndex = normalRot%varTransformIndices(iBcVar)
          !@todo ! consider rotation of the derivatives
          !@todo if(equation%eq_kind .eq. 'navier_stokes_2d') then
          !@todo bcIndex_grad = normalRot%derTransformIndices(2:3)
          !@todo -normalRot%derTransformIndices(1)
          !@todo else
            write(logUnit(1),*) 'ERROR in modg_nodal_bnd_gradient: ' // &
              &                 'rotation of normal gradients ' //         &
              &                 'supported for navier_stokes only, stopping ...'
            call tem_abort()
        end if
        !@todo else
          bcIndex = iBcVar
          !@todo bcIndex_grad(1) = 1
          !@todo bcIndex_grad(2) = 2
        !@todo end if

        ! Create the modal representation on the face.
        select case(bc%state_gradient(iBcVar)%style)
        case('neumann')
          ! @todo replace by subroutine call, due to OpenMP
          ! Extrapolate the values
          ! ... derivatives in the first spatial direction
          pointFace_derX(:,bcIndex) = modg_bnd_extrapolate(                 &
            &                           nVals  = nQuadPoints,               &
            &                           faceOp = pointValOp_derX(:,bcIndex) )
          ! ... derivatives in the second spatial direction
          pointFace_derY(:,bcIndex) = modg_bnd_extrapolate(                 &
            &                           nVals  = nQuadPoints,               &
            &                           faceOp = pointValOp_derY(:,bcIndex) )
        case('dirichlet')
            pointFace_derX(:,bcIndex) = modg_bnd_mirrorPoint( &
              &  bc           = bc%state_gradient(iBcVar),    &
              &  nPoints      = nQuadPoints,                  &
              &  time         = time,                         &
              &  varSys       = equation%varSys,              &
              &  currentFace  = currentFace,                  &
              &  currentLevel = currentLevel                  )

            pointFace_derY(:,bcIndex) = modg_bnd_mirrorPoint( &
              &  bc           = bc%state_gradient(iBcVar),    &
              &  nPoints      = nQuadPoints,                  &
              &  time         = time,                         &
              &  varSys       = equation%varSys,              &
              &  currentFace  = currentFace,                  &
              &  currentLevel = currentLevel                  )

        case default
          write(logUnit(1),*) 'ERROR in modg_nodal_bnd for gradient: &
            & unknown bnd style, stopping ...'
          call tem_abort()
        end select
        ! copy the derivatives to the large array to convert all
        ! variables n2m at once
        pointFace(:poly_proj%body_2d%min_dofs,iBcVar+equation%varSys%nScalars) &
          &  =  pointFace_derX(:poly_proj%body_2d%min_dofs,iBcVar)
        pointFace(:poly_proj%body_2d%min_dofs,iBcVar+3*equation%varSys%nScalars)&
          &  =  pointFace_derY(:poly_proj%body_2d%min_dofs,iBcVar)
      end do !iBCVar
    end if ! derivative

    ! take care of bc_trafo, i.e. transform from primitive to conservative.
    if(.not.bc%bc_trafo%identity) then
      call bc%bc_trafo%from( equation = equation,    &
        &                    instate  = pointFace,   &
        &                    material = faceMaterial )
    end if

    ! transform back to modal space if necessary
    ! --> oversamp nodal space
    call ply_poly_project_n2m(me         = poly_proj,    &
      &                       dim        = 2 ,           &
      &                       nVars      = nVars,        &
      &                       nodal_data = pointFace,    &
      &                       modal_data = oversamp_Face )
    ! --> oversamp modal space
    call ply_convertFromOversample(modalCoeffs = oversamp_face, &
      &                            poly_proj   = poly_proj,     &
      &                            nDim        = 2,             &
      &                            state       = modalFace      )
    ! --> modal space


  end subroutine modg_nodal_bnd
  ! ***************************************************************************


  ! ***************************************************************************
  !> Set boundary values in a modal way
  subroutine modg_modal_bnd( bc, faceOp, poly_proj, equation, normalRot,   &
    &                        time, currentFace, currentLevel, nQuadPoints, &
    &                        nDofs, oversamp_dofs, modalFace, faceMaterial )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> Data for the projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> Rotation indices to rotate global coordinate system into face normal
    !! coordinate system.
    type(coordRotation_type), intent(in) :: normalRot
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    integer, intent(in) :: currentFace
    !> the level to compute on
    integer, intent(in) :: currentLevel
    !> Number of quadurature points and  Number of Dofs for the face
    integer, intent(in) :: nQuadPoints, nDofs, oversamp_dofs
    !> The modal representation on the boundary face.
    real(kind=rk), intent(inout):: modalFace(:,:)
    !> current face to compute on, used for getting index
    !> The material of the boundary face.
    !! First dimension is the number of points on the face.
    !! Second dimension is the number of material parameters.
    real(kind=rk), intent(in), optional :: faceMaterial(:,:)
    ! ---------------------------------------------------------------------------
    integer :: iBcVar, bcIndex
    real(kind=rk), allocatable :: pointFace(:,:)
    real(kind=rk), allocatable :: tmpModal(:,:)
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! ---------------------------------------------------------------------------

    ! take care of bc_trafo, i.e. transform from conservative to primitive.
    if(.not.bc%bc_trafo%identity) then
      call bc%bc_trafo%to( equation = equation,    &
        &                  instate  = faceOp,      &
        &                  material = faceMaterial )
    end if

    ! Loop over all the quantitites of the equation system and generate the
    ! modal representation for it.
    do iBcVar = 1, equation%varSys%nScalars

      ! Is a rotation to face normal representation necessary?
      ! So, we get the index for the normal direction on the face.
      if(bc%bc_normal_vec) then
        bcIndex = normalRot%varTransformIndices(iBcVar)
      else
        bcIndex = iBcVar
      end if

      ! @todo PV 20160131 Boundary variable need not to be spacetime
      ! function so it is better to introduce vartype in tem_varSys_op_type to
      ! check for st_fun.
      ! convert c_ptr stfunList to fortran pointer
      call c_f_pointer( equation%varSys%method%val( bc%state(iBcVar)%varPos ) &
        &                              %method_data, fPtr )

      ! Create the modal representation on the face.
      select case(bc%state(iBcVar)%style)
      case('neumann')
        ! In case of dirichlet boundary condition we would like
        ! to prescribe the derivative of the boundary value, so:
        ! Extrapolate the values
        ! @todo replace by subroutine call, due to OpenMP
        modalFace(:,bcIndex) = modg_bnd_extrapolate(        &
          &                      nVals  = nDofs,            &
          &                      faceOp = faceOp(:,bcIndex) )
                           !    & (modg%maxPolyDegree+1)**2, faceOp(:,bcIndex) )
      case('dirichlet')
        ! In case of dirichlet boundary condition we would like
        ! to prescribe the boundary value, so:
        ! Mirror the values around the prescribed one.
        ! @todo replace by subroutine call, due to OpenMP

        ! If boundary variable is refered to zero_const then
        ! set modalFace value to zero
        if (bc%varDict%val(iBcVar)%value == 'zero_const') then
          modalFace(:, bcIndex) = 0.0_rk

        else

          ! the function is not constant, we transfer to point values,
          ! mirror pointwise and tranfer back to a modal representation.
          allocate(pointFace(nQuadPoints,1) )
          allocate(tmpModal(oversamp_dofs,1) )

          ! Mirror pointwise
          pointFace(:,1) = modg_bnd_mirrorPoint(              &
            &                bc           = bc%state(iBcVar), &
            &                varSys       = equation%varSys,  &
            &                time         = time,             &
            &                currentFace  = currentFace,      &
            &                currentLevel = currentLevel,     &
            &                nPoints      = nQuadPoints       )

          call ply_poly_project_n2m(me         = poly_proj, &
            &                       dim        = 2 ,        &
            &                       nVars      = 1,         &
            &                       nodal_data = pointFace, &
            &                       modal_data = tmpModal   )
          ! --> oversamp modal space

          call ply_convertFromOversample(                     &
            &      modalCoeffs = tmpModal,                    &
            &      poly_proj   = poly_proj,                   &
            &      nDim        = 2,                           &
            &      state       = modalFace(:,bcIndex:bcIndex) )
          ! --> modal space

          deallocate(pointFace)
          deallocate(tmpModal)
        end if
      case default
        write(logUnit(1),*) 'ERROR in modg_modalBndRep: unknown bnd style, stopping ...'
        call tem_abort()
      end select

    end do !BCIndex

    ! take care of bc_trafo, i.e. transform from primitive to conservative.
    if(.not.bc%bc_trafo%identity) then
      call bc%bc_trafo%from( equation = equation,    &
        &                    instate  = modalFace,   &
        &                    material = faceMaterial )
    end if

  end subroutine modg_modal_bnd
  ! ***************************************************************************


  ! ***************************************************************************
  !> Function to extrapolate face values for a given boundary condition in
  !! physical or modal space.
  function modg_bnd_extrapolate( nVals, faceOp ) result(modalFace)
    ! ---------------------------------------------------------------------------
    !> The number of coefficients to extrapolate
    integer, intent(in) :: nVals
    !> Modal representation on the face of the neighboring element.
    real(kind=rk), intent(in) :: faceOp(:)
    !> The extrapolated modal representation.
    real(kind=rk), allocatable :: modalFace(:)
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    allocate( modalFace(nVals) )

    ! Just extrapolate the values on the face, which means: copy the
    ! modal representation.

    modalFace(:) = faceOp(:)

  end function modg_bnd_extrapolate
  ! ***************************************************************************


  ! ***************************************************************************
  !> Function to mirror pointvalues for a given boundary conditions.
  function modg_bnd_mirrorPoint( bc,  varSys, time,  currentFace,        &
      &                          currentLevel, nPoints ) result(pointFace)
    ! ---------------------------------------------------------------------------
    !> The boundary state.
    type(tem_bc_state_type), intent(in) :: bc
    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> The current absolute time.
    type(tem_time_type), intent(in) :: time
    !> current face used to compute correct index in indices array
    integer, intent(in) :: currentFace
    !> current level
    integer, intent(in) :: currentLevel
    !> The number of point values to be mirrored.
    integer, intent(in) :: nPoints
    !> The mirrored isNodalScheme representation.
    real(kind=rk) :: pointFace(npoints)
    ! ---------------------------------------------------------------------------
    integer :: idx_start, idx_end
    ! ---------------------------------------------------------------------------

    ! the boundary values at the quadrature nodes which are stored fore each
    ! variable in the pointData and are  accessed via Indices
    idx_start = (currentFace-1)*nPoints+1
    idx_end = currentFace*nPoints
    call tem_startTimer( timerHandle = atl_timerHandles%readBC )
    call varSys%method%val(bc%varPos)%get_valOfIndex(                        &
      & varSys  = varSys,                                                    &
      & time    = time,                                                      &
      & iLevel  = currentLevel,                                              &
      & idx     = bc%pntIndex%indexLvl(currentLevel)%val(idx_start:idx_end), &
      & nVals   = nPoints,                                                   &
      & res     = pointFace                                                  )
    call tem_stopTimer( timerHandle = atl_timerHandles%readBC )

  end function modg_bnd_mirrorPoint
  ! ***************************************************************************


end module atl_modg_bnd_module
