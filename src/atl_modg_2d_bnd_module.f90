! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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
module atl_modg_2d_bnd_module
  use, intrinsic :: iso_c_binding,    only: c_f_pointer

  ! Treelm modules
  use env_module,                     only: rk
  use tem_aux_module,                 only: tem_abort
  use tem_faceData_module,            only: tem_invFace_map
  use tem_coordinate_module,          only: coordRotation_type
  use tem_bc_module,                  only: tem_bc_state_type
  use tem_time_module,                only: tem_time_type
  use tem_spacetime_fun_module,       only: tem_spacetime_for
  use tem_logging_module,             only: logUnit
  use tem_timer_module,               only: tem_startTimer, &
    &                                       tem_stopTimer
  use tem_varSys_module,              only: tem_varSys_type
  use tem_spacetime_fun_module,       only: tem_st_fun_listElem_type


  use atl_bc_header_module,           only: atl_boundary_type
  use atl_boundary_module,            only: atl_level_boundary_type
  use atl_facedata_module,            only: atl_facedata_type
  use atl_kerneldata_module,          only: atl_statedata_type
  use atl_modg_2d_scheme_module,      only: atl_modg_2d_scheme_type
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

  public :: atl_modg_2d_set_bnd, atl_modg_2d_bnd


contains


  ! ****************************************************************************
  !> Subroutine to set face values to impose boundary conditions
  !! at a certain point of the domain. The subroutine is operating levelwise.
  subroutine atl_modg_2d_set_bnd( bc, boundary, facedata, statedata, modg,   &
    &                             equation, time, mesh, poly_proj, nodalBnd, &
    &                             material, currentLevel                     )
    ! ---------------------------------------------------------------------------
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary
    !> The face data on the current level
    type(atl_facedata_type), intent(inout) :: facedata
    !> The state data on the current level
    type(atl_statedata_type), intent(inout) :: statedata
    !> The parameters of th the modg scheme.
    type(atl_modg_2d_scheme_type), intent(inout) :: modg
    !> Data for the projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The underlying equation system
    type(atl_equations_type), intent(in) :: equation
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> The description of the mesh on the current level.
    type(atl_cube_elem_type),  intent(in) :: mesh
    !> Set boundaries in nodal fashion by default? If set to false,
    !! the boundaries may still be set in nodal way whenever necessary (e.g.
    !! boundaries which have space-time dependence, etc.)
    logical, intent(in) :: nodalBnd
    !> Material description of the faces contained in boundary. One
    !! entry for each spatial direction, i.e. x,y.
    type(atl_faceMaterialData_type), intent(in), optional :: material(2)
    !> the level to compute on
    integer, intent(in) :: currentLevel
    ! ---------------------------------------------------------------------------
    integer :: facePos, neighPos, neighAlign, nScalars
    integer :: iBC, iDir, iFace, iAlign, iglobFace
    integer :: nquadpoints, ndofs, oversamp_dofs
    real(kind=rk), allocatable :: faceOp(:,:)
    real(kind=rk) :: elemLength
    real(kind=rk) :: bndBaryCoord(1:3)
    real(kind=rk), allocatable :: faceMaterial(:,:)
    integer :: signFact, iMode, jMode, nModes, modePos
    real(kind=rk) :: corrector(poly_proj%body_1D%ndofs,size(statedata%state,3))
    real(kind=rk) :: sidefact
    ! ---------------------------------------------------------------------------

    if(present(material)) then
      allocate( faceMaterial( size(material(1)%leftElemMaterialDat,2), &
                            & size(material(1)%leftElemMaterialDat,3)) )
    end if

    nScalars = equation%varSys%nScalars

    ! The length of an element
    elemLength = mesh%length

    ! get correct amount of quadrature points on the face and degree
    ! due to projection  method. oversamp_dof are
    ! used for the oversampling loop
    nquadpoints   = poly_proj%body_1D%nQuadPoints
    ndofs         = poly_proj%body_1D%ndofs
    oversamp_dofs = poly_proj%body_1D%overSamp_dofs

! @todo add other variables to private if necessary
    !!!!OMP PARALLEL &
    !!!!OMP PRIVATE(iBC, iDir, iAlign, iFace) &
    !!!!OMP DEFAULT(shared)

    ! Iterate over all the boundaries and set the right face values for
    ! the boundaries on all relevant faces.
    do iBC = 1, boundary%nBCs

      ! Now, we iterate over all the faces with this boundary conditions and
      ! set the corresponding face values
      iglobFace = 0

      ! calculate the highest mode which will be modified for zero grd BC
      nModes = ceiling(nDofs * bc(iBC)%neumann_mode_fraction)

      do iDir = 1,2
        allocate( faceOp(size(facedata%faceRep(iDir)%dat,2), &
          &              size(facedata%faceRep(iDir)%dat,3)) )

        do iAlign = 1,2

          neighAlign = tem_invFace_map(iAlign)
          sidefact = (-1)**neighAlign

          do iFace = 1, boundary%bnd(iBC)%faces(iDir, iAlign)%facePos%nVals
            ! count global iterator
            iglobFace = iglobFace +1

            ! Create the modal representation on the face for the current
            ! face. We need the modal representation of the neighboring fluid
            ! element for that
            neighPos = boundary%bnd(iBC)%faces(iDir, iAlign)%neighPos%val(iFace)

            faceOp = facedata%faceRep(iDir)%dat(neighPos,:,:,neighAlign)

            if (bc(iBC)%enforce_zero_grad) then

              ! faceOp is the state that we tell the boundary condition is
              ! present at the face.
              ! Initially set this to the first mode as this will be used if only
              ! one mode is to be used for the extrapolation.
              if ( iDir == 1 ) then
                ! mode coefficients: a_{1,j}, j = 1,...,nDoFs
                faceOp(:,:nScalars) = statedata%state(neighPos,1::nDoFs,:)
              else if ( iDir == 2 ) then
                ! mode coefficients: a_{i,1}, i = 1,...,nDoFs
                faceOp(:,:nScalars) = statedata%state(neighPos,1:nDoFs,:)
              end if

              if (nModes > 1) then
                ! Only need to actually compute something if we use more
                ! then just the first mode for the extrapolation

                ! For x direction, jMode is the y mode
                ! For y direction, jMode is the x mode

                ! signfact is used to keep track of the alternating sign on the
                ! left end of the element.
                ! (sidefact decides whether we are left (=-1) or right(=1)),
                ! while we start with signfact=sidefact to have the right augury
                signfact = int(sidefact)

                ! To enforce a zero gradient we compute the last mode
                ! (at nModes) to match all previous modes contributions in the
                ! derivative. This "correcting" mode is stored in corrector.
                corrector = 0.0_rk

                do iMode=1,nModes-2

                  ! By taking the derivative all modes are shifted down by 1.
                  ! Thus, the relevant modes are those up to nModes-1. However,
                  ! the one we want to set to enforce the 0 gradient, is one in
                  ! nModes-1. Thus, we need to sum up to nModes-2...
                  ! For x direction, sum iMode i.e. the x mode
                  ! For y direction, sum iMode i.e. the y mode
                  !!signfact = int(signfact * sidefact)

                  do jMode = 1, ndofs
                    ! mode position in state array
                    if ( iDir == 1 ) then
                      modePos = (iMode+1) + (jMode-1)*(nModes)
                    else if ( iDir == 2 ) then
                      modePos = jMode + (iMode)*(nModes)
                    end if
                    ! At the same time we can compute the value at the face by
                    ! summing the modes, while considering the right signfact for that.
                    faceOp(jMode,:nscalars) = faceOp(jMode,:nscalars)       &
                      &                       + signfact                    &
                      &                         * statedata%state(neighpos, &
                      &                                           modePos,  &
                      &                                           :nscalars )

                    ! The derivative contribution of iMode+1, we subtract it from
                    ! the corrector to obtain a corrector, which will balance out
                    ! all other modes.
                    corrector(jMode,:) = corrector(jMode,:)            &
                      &       +(signfact*sidefact*(iMode*(iMode+1))/2) &
                      &         * statedata%state(neighpos,modePos,:)

                  end do ! jMode
                    signfact = int(signfact * sidefact)
                end do ! iMode

                ! The last mode is not used, instead we use the corrector as
                ! last mode. As that was computed in the derivative, it needs to
                ! "integrated" (divison by the factor for the Legendre derivative)
                do jMode = 1, nDofs
                  faceOp(jMode,:nScalars) = faceOp(jMode,:nScalars) &
                    &                        - (sidefact*2*corrector(jMode,:)) &
                    &                       / ((nModes)*(nModes-1))

                end do ! jMode

              end if ! nModes > 1

            end if

            ! The position of the face with the current boundary condition
            ! inside the face representation.
            facePos = boundary%bnd(iBC)%faces(iDir, iAlign)%facePos%val(iFace)

            ! get the barycentric coordinate of the (virtual) boundary element
            bndBaryCoord(1:3) = mesh%bary_coord(neighPos,1:3)
            bndBaryCoord(iDir) = bndBaryCoord(iDir) &
              & + ( (-1.0_rk)**neighAlign ) * elemLength

            ! If material information is passed to this function,
            ! we read out the material of the boundary face and pass it
            ! to the modg_bnd function. If not, we call modg_bnd without
            ! this information.
            if( present(material) ) then

              ! We read out the left face material (at a boundary face left
              ! and right face material should be equal).
              faceMaterial(:,:) = material(iDir)%leftElemMaterialDat(iFace,:,:)
              ! @todo JZ: replace by subroutine call, due to OpenMP
              call atl_modg_2d_bnd(                             &
                & bc              = bc(iBC),                    &
                & faceOp          = faceOp,                     &
                & poly_proj       = poly_proj,                  &
                & normalRot       = equation%varRotation(iDir), &
                & nDerivatives    = equation%nDerivatives,      &
                & equation        = equation,                   &
                & isNodalScheme   = nodalBnd,                   &
                & time            = time,                       &
                & currentFace     = iglobFace,                  &
                & currentLevel    = currentLevel,               &
                & nQuadPoints     = nQuadPoints,                &
                & nDofs           = nDofs,                      &
                & oversamp_dofs   = oversamp_dofs,              &
                & modalFace       = facedata%faceRep(iDir)      &
                  &                  %dat(facePos,:,:,iAlign),  &
                & faceMaterial    = faceMaterial                )
            else
              ! @todo JZ: replace by subroutine call, due to OpenMP
              call atl_modg_2d_bnd(                              &
                &  bc              = bc(iBC),                    &
                &  faceOp          = faceOp,                     &
                &  poly_proj       = poly_proj,                  &
                &  normalRot       = equation%varRotation(iDir), &
                &  nDerivatives    = equation%nDerivatives,      &
                &  equation        = equation,                   &
                &  isNodalScheme   = nodalBnd,                   &
                &  time            = time,                       &
                &  currentFace     = iglobFace,                  &
                &  currentLevel    = currentLevel,               &
                &  nquadpoints     = nquadpoints,                &
                &  ndofs           = ndofs,                      &
                &  oversamp_dofs   = oversamp_dofs,              &
                &  modalFace       = facedata%faceRep(iDir)      &
                  &                  %dat(facePos,:,:,iAlign)    )
            end if !material
          end do ! iFace
        end do ! iAlign
        deallocate(faceOp)
      end do ! iDir
    end do ! iBC

    !!!!OMP END PARALLEL

  end subroutine atl_modg_2d_set_bnd
  ! ****************************************************************************


  ! ****************************************************************************
  !> Subroutine to create the modal representation for a ceratin boundary face.
  subroutine atl_modg_2d_bnd( bc, faceOp, poly_proj, nDerivatives, equation,  &
    &                         normalRot, isNodalScheme, time, currentFace,    &
    &                         currentLevel,nquadpoints, ndofs, oversamp_dofs, &
    &                         modalFace, faceMaterial                         )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> The parameters for projection method.
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
    !> Current face to compute, needed to get correct position in index array
    integer, intent(in) :: currentFace
    !> the level to compute on
    integer, intent(in) :: currentLevel
    !> integers for allocation of temp arrays, depend on number of quadrature
    !! points and for modal values number of dofs
    integer, intent(in) :: nquadpoints, ndofs, oversamp_dofs
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
      call modg_2d_nodal_bnd(bc              = bc,              &
        &                    faceOp          = faceOp,          &
        &                    poly_proj       = poly_proj,       &
        &                    equation        = equation,        &
        &                    normalRot       = normalRot,       &
        &                    time            = time,            &
        &                    currentFace     = currentFace,     &
        &                    currentLevel    = currentLevel,    &
        &                    nDerivatives    = nDerivatives,    &
        &                    nQuadPoints     = nQuadPoints,     &
        &                    oversamp_dofs   = oversamp_dofs,   &
        &                    modalFace       = modalFace,       &
        &                    faceMaterial    = faceMaterial     )

    else

      ! For linear equations it could be possible to impose boundary conditions
      ! in modal space directly.
      call modg_2d_modal_bnd( bc            = bc,             &
        &                     faceOp        = faceOp,         &
        &                     poly_proj     = poly_proj,      &
        &                     equation      = equation,       &
        &                     normalRot     = normalRot,      &
        &                     time          = time,           &
        &                     currentFace   = currentFace,    &
        &                     currentLevel  = currentLevel,   &
        &                     nQuadPoints   = nQuadPoints,    &
        &                     nDofs         = nDofs,          &
        &                     oversamp_dofs = oversamp_dofs,  &
        &                     modalFace     = modalFace,      &
        &                     faceMaterial  = faceMaterial    )

      if( equation%nDerivatives > 0 ) then
        write(logUnit(1),*) 'ERROR in atl_modg_2d_bnd: not able to ' // &
          & 'boundary conditions for higher order equations in modal way,' // &
          & 'stopping ...'
        call tem_abort()
      end if

    end if

  end subroutine atl_modg_2d_bnd
  ! ****************************************************************************


  ! ****************************************************************************
  !> Set boundary values in a nodal way
  subroutine modg_2d_nodal_bnd( bc, faceOp, poly_proj, equation, normalRot, &
    &                           time, currentFace, currentLevel,            &
    &                           nDerivatives, nquadpoints, oversamp_dofs,   &
    &                           modalFace, faceMaterial                     )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> The parameters for projection method.
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> Rotation indices to rotate global coordinate system into face normal
    !! coordinate system.
    type(coordRotation_type), intent(in) :: normalRot
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> current face to cpmoute on, needed to access correct index
    integer, intent(in) :: currentFace
    !> the level to compute on
    integer, intent(in) :: currentLevel
    !> The number of derivative boundaries to be set
    integer, intent(in) :: nDerivatives
    !> integers for allocation of temp arrays, depend on number of quadrature
    !! points and for modal values number of dofs
    integer, intent(in) :: nquadpoints, oversamp_dofs
    !> result of the bnd routine, modal coefficent on the boundary faces
    real(kind=rk), intent(inout) :: modalFace(:,:)
    !> The material of the boundary face.
    !! First dimension is the number of points on the face.
    !! Second dimension is the number of material parameters.
    real(kind=rk), intent(in), optional :: faceMaterial(:,:)
    ! ---------------------------------------------------------------------------
    integer :: iBcVar, bcIndex, iter, nVars
    real(kind=rk), allocatable :: pointValOp(:,:), pointFace(:,:), tmpFace(:,:)
    real(kind=rk), allocatable :: pointValOp_derX(:,:), pointValOp_derY(:,:), &
      &                           pointFace_derX(:,:), pointFace_derY(:,:)
  ! ---------------------------------------------------------------------------

    nVars = equation%varSys%nScalars+2*nDerivatives*equation%varSys%nScalars

    ! temporary arrays
    allocate(pointFace(nQuadPoints, nVars))
    allocate(pointValOp(nQuadPoints,nVars))
    allocate(tmpFace(oversamp_dofs, nVars))

    do iter=lbound(modalFace,2),ubound(modalFace,2)
      modalFace(:,iter) = 0.0_rk
    end do

    ! --> modal space
    call ply_convert2oversample(                             &
      & state       = faceOp(:poly_proj%body_1d%nDofs,:), &
      & poly_proj   = poly_proj,                             &
      & nDim        = 1,                                     &
      & modalCoeffs = tmpface,                               &
      & nScalars    = nVars  )
    ! --> oversamp modal space
    ! Now, we transform the modal representation of this element to nodal
    call ply_poly_project_m2n( me         = poly_proj,  &
      &                        dim        = 1 ,         &
      &                        nVars      = nVars,      &
      &                        nodal_data = pointValOp, &
      &                        modal_data = tmpFace     )
    ! --> oversamp nodal space


    ! take care of bc_trafo, i.e. transform from conservative to primitive.
    if(.not.bc%bc_trafo%identity) then
      if( present(faceMaterial) ) then
        call bc%bc_trafo%to( equation = equation,    &
          &                  instate  = pointValOp,  &
          &                  material = faceMaterial )
      else
        call bc%bc_trafo%to( equation = equation, instate = pointValOp )
      end if
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
        pointFace(:,bcIndex) = modg_2d_bnd_extrapolate(         &
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
          call modg_2d_bnd_mirrorPoint( bc           = bc%state(iBcVar),      &
            &                           nPoints      = nQuadPoints,           &
            &                           time         = time,                  &
            &                           varSys       = equation%varSys,       &
            &                           currentFace  = currentFace,           &
            &                           pointFace    = pointFace(:,bcIndex),  &
            &                           currentLevel = currentLevel           )
        end if
      case default
        write(logUnit(1),*) 'ERROR in modg_2d_nodalBndRep:'
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
        pointValOp_derX(:nQuadPoints,iBcVar) = &
          & pointValOp( :nQuadPoints,          &
          &             iBcVar+equation%varSys%nScalars       )
        pointValOp_derY(:nQuadPoints,iBcVar) = &
          & pointValOp( :nQuadPoints,          &
          &             iBcVar+2*equation%varSys%nScalars     )

        ! Is a rotation to face normal representation necessary?
        ! So, we get the index for the normal direction on the face.
        if(bc%bc_normal_vec_gradient) then
          !@todo bcIndex = normalRot%varTransformIndices(iBcVar)
          !@todo ! consider rotation of the derivatives
          !@todo if(equation%eq_kind .eq. 'navier_stokes_2d') then
          !@todo bcIndex_grad = normalRot%derTransformIndices(2:3)
          !@todo -normalRot%derTransformIndices(1)
          !@todo else
            write(logUnit(1),*) 'ERROR in modg_2d_nodal_bnd_gradient: ' // &
              &                 'rotation of normal gradients ' //         &
              &                 ' supported for navier_stokes_2d only, stopping ...'
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
          pointFace_derX(:,bcIndex) = modg_2d_bnd_extrapolate(              &
            &                           nVals  = nQuadPoints,               &
            &                           faceOp = pointValOp_derX(:,bcIndex) )
          ! ... derivatives in the second spatial direction
          pointFace_derY(:,bcIndex) = modg_2d_bnd_extrapolate(              &
            &                           nVals  = nQuadPoints,               &
            &                           faceOp = pointValOp_derY(:,bcIndex) )
        case('dirichlet')
            call modg_2d_bnd_mirrorPoint(                  &
              &  bc           = bc%state_gradient(iBcVar), &
              &  nPoints      = nQuadPoints,               &
              &  time         = time,                      &
              &  varSys       = equation%varSys,           &
              &  currentFace  = currentFace,               &
              &  pointFace    = pointFace_derX(:,bcIndex), &
              &  currentLevel = currentLevel               )

            call modg_2d_bnd_mirrorPoint(                  &
              &  bc           = bc%state_gradient(iBcVar), &
              &  nPoints      = nQuadPoints,               &
              &  time         = time,                      &
              &  varSys       = equation%varSys,           &
              &  currentFace  = currentFace,               &
              &  pointFace    = pointFace_derY(:,bcIndex), &
              &  currentLevel = currentLevel               )

        case default
          write(logUnit(1),*) 'ERROR in modg_2d_nodal_bnd for gradient: &
            & unknown bnd style, stopping ...'
          call tem_abort()
        end select
        ! copy the derivatives to the large array to convert all
        ! variables n2m at once
        pointFace(:nQuadPoints,iBcVar+equation%varSys%nScalars) &
          & = pointFace_derX(:,iBcVar)
        pointFace(:nQuadPoints,iBcVar+2*equation%varSys%nScalars) &
          & = pointFace_derY(:,iBcVar)
      end do !iBCVar
    end if ! derivative

    ! take care of bc_trafo, i.e. transform from primitive to conservative.
    if(.not.bc%bc_trafo%identity) then
      if(present(faceMaterial)) then
        call bc%bc_trafo%from( equation = equation,    &
          &                    instate  = pointFace,   &
          &                    material = faceMaterial )
      else
        call bc%bc_trafo%from( equation = equation, instate = pointFace  )
      end if
    end if

    ! Convert back to modal values (1D FPT has to reside in an OMP parallel region)
    call ply_poly_project_n2m(me         = poly_proj,  &
      &                       dim        = 1 ,         &
      &                       nVars      = nVars,      &
      &                       nodal_data = pointFace,  &
      &                       modal_data = tmpFace     )

   ! --> oversamp modal space
    call ply_convertFromOversample(modalCoeffs = tmpFace,   &
      &                            poly_proj   = poly_proj, &
      &                            nDim        = 1,         &
      &                            state       = modalFace  )
    ! --> modal space

  end subroutine modg_2d_nodal_bnd
  ! ****************************************************************************


  ! ****************************************************************************
  !> Set boundary values in a modal way
  subroutine modg_2d_modal_bnd( bc, faceOp, poly_proj, equation, normalRot, &
      &                       time, currentFace, currentLevel, nquadpoints, &
      &                       ndofs, oversamp_dofs, modalFace, faceMaterial )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> The parameters for projection method.
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> Rotation indices to rotate global coordinate system into face normal
    !! coordinate system.
    type(coordRotation_type), intent(in) :: normalRot
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> Current Face to compute
    integer, intent(in) :: currentFace
    !> the level to compute on
    integer, intent(in) :: currentLevel
    !> integers for allocation of temp arrays, depend on number of quadrature
    !! points and for modal values number of dofs
    integer, intent(in) :: nquadpoints, ndofs, oversamp_dofs
    !> result of the bnd routine, modal coefficent on the boundary faces
    real(kind=rk), intent(inout) :: modalFace(:,:)
    !> The material of the boundary face.
    !! First dimension is the number of points on the face.
    !! Second dimension is the number of material parameters.
    real(kind=rk), intent(in), optional :: faceMaterial(:,:)
    ! ---------------------------------------------------------------------------
    integer :: iBcVar, bcIndex, iDeg, iter
    real(kind=rk), allocatable :: pointFace(:,:)
    real(kind=rk), allocatable :: tmpModal(:,:)
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! ---------------------------------------------------------------------------

    do iter=lbound(modalFace,2),ubound(modalFace,2)
    modalFace(:,iter) = 0.0_rk
       end do


    ! take care of bc_trafo, i.e. transform from conservative to primitive.
    if(.not.bc%bc_trafo%identity) then
      if( present(faceMaterial) ) then
        call bc%bc_trafo%to( equation = equation,    &
          &                  instate  = faceOp,      &
          &                  material = faceMaterial )
      else
        call bc%bc_trafo%to( equation = equation, instate = faceOp )
      end if
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
        modalFace(:,bcIndex) = modg_2d_bnd_extrapolate(     &
          &                      nVals  = ndofs,            &
          &                      faceOp = faceOp(:,bcIndex) )
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
          ! @todo JZ: Mirroring is not thread safe yet!
          call modg_2d_bnd_mirrorPoint( bc           = bc%state(iBcVar), &
            &                           nPoints      = nQuadPoints,      &
            &                           time         = time,             &
            &                           varSys       = equation%varSys,  &
            &                           currentFace  = currentFace,      &
            &                           pointFace    = pointFace(:,1),   &
            &                           currentLevel = currentLevel      )

          ! transform back to modal values
          call ply_poly_project_n2m(me         = poly_proj, &
            &                       dim        = 1 ,        &
            &                       nVars      = 1,         &
            &                       nodal_data = pointFace, &
            &                       modal_data = tmpModal   )

          !! usually generic call to convertFromOversample
          do iDeg = 1, poly_proj%body_1d%min_dofs
            modalFace(iDeg,bcIndex) = tmpModal(iDeg,1)
          end do


          deallocate(pointFace)
          deallocate(tmpModal)

        end if
      case default
        write(logUnit(1),*) 'ERROR in modg_2d_modalBndRep: unknown bnd style, stopping ...'
        call tem_abort()
      end select

    end do

    ! take care of bc_trafo, i.e. transform from primitive to conservative.
    if(.not.bc%bc_trafo%identity) then
      if( present(faceMaterial) ) then
        call bc%bc_trafo%from( equation = equation,    &
          &                    instate  = modalFace,   &
          &                    material = faceMaterial )
      else
        call bc%bc_trafo%from( equation = equation, instate = modalFace )
      end if
    end if

  end subroutine modg_2d_modal_bnd
  ! ****************************************************************************


  ! ****************************************************************************
  !> Function to extrapolate face values for a given boundary condition in
  !! physical or modal space.
  function modg_2d_bnd_extrapolate( nVals, faceOp ) result(modalFace)
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

  end function modg_2d_bnd_extrapolate
  ! ****************************************************************************


  ! ****************************************************************************
  !> Function to mirror pointvalues for a given boundary conditions.
  subroutine modg_2d_bnd_mirrorPoint( bc, nPoints, time,  varSys,          &
      &                               currentFace, pointFace, currentLevel )
    ! ---------------------------------------------------------------------------
    !> The boundary state.
    type(tem_bc_state_type), intent(in) :: bc
    !> The number of point values to be mirrored.
    integer, intent(in) :: nPoints
    !> The current absolute time.
    type(tem_time_type), intent(in) :: time
    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> The mirrored isNodalScheme representation.
    real(kind=rk), intent(inout):: pointFace(:)
    !> Current Face to get the correct position in index array
    integer,intent(in) :: currentFace
    !> current level
    integer, intent(in) :: currentLevel
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

  end subroutine modg_2d_bnd_mirrorPoint
  ! ****************************************************************************

end module atl_modg_2d_bnd_module


