! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
module atl_modg_1d_bnd_module
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
  use tem_varSys_module,              only: tem_varSys_type
  use tem_spacetime_fun_module,       only: tem_spacetime_fun_type, &
    &                                       tem_st_fun_listElem_type
  use treelmesh_module,               only: treelmesh_type

  use atl_bc_header_module,           only: atl_boundary_type
  use atl_boundary_module,            only: atl_level_boundary_type
  use atl_facedata_module,            only: atl_facedata_type
  use atl_kerneldata_module,          only: atl_statedata_type
  use atl_equation_module,            only: atl_equations_type
  use atl_reference_element_module,   only: atl_refToPhysCoord
  use atl_cube_elem_module,           only: atl_cube_elem_type
  use atl_materialPrp_module,         only: atl_faceMaterialData_type
  use ply_poly_project_module,        only: ply_poly_project_type,      &
    &                                       assignment(=),              &
    &                                       ply_get_quadpoints_faces_1d

  implicit none

  private

  public :: atl_modg_1d_set_bnd, atl_modg_1d_bnd

contains

  !> Subroutine to set face values to impose boundary conditions.
  !!
  !! We set the "outer" state according to the configure boundary condition.
  !!
  !! For all variables, where a Dirichlet condition is imposed, this value
  !! is fixed and simply set. It might be that a variable transformation is
  !! necessary, or we need to perform a transformation to nodal space,
  !! all of this is taken care of in [[atl_modg_1d_bnd]], which is called
  !! in this routine.
  !!
  !! For Neumann boundaries, the default approach is to just use the value
  !! at the boundary also for the "outer" face value, such that in the flux
  !! left and right state of the variable is always the same.
  !! However, this is sensitive to oscillations and may easily cause stability
  !! issues. Alternatively we can also compute a new value by modifying the
  !! polynomial in the element to have a zero gradient on the boundary
  !! enforced and then use this value for the extrapolated "outer" state
  !! instead. For details see the
  !! [Neumann boundary conditions](page/neumann_boundaries.md].
  !!
  !! The subroutine is operating levelwise.
  subroutine atl_modg_1d_set_bnd( bc, boundary, facedata, statedata,         &
    &                             poly_proj, material, equation, tree, time, &
    &                             mesh                                       )
    ! ---------------------------------------------------------------------------
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary
    !> The face data on the current level
    type(atl_facedata_type), intent(inout) :: facedata
    !> The state data on the current level
    type(atl_statedata_type), intent(inout) :: statedata
    !> Data for the projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The underlying equation system
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> The description of the mesh on the current level.
    type(atl_cube_elem_type),  intent(in) :: mesh
    !> Material description of the faces contained in boundary. One
    !! entry for each spatial direction, i.e. x,y.
    type(atl_faceMaterialData_type), intent(in), optional :: material
    ! ---------------------------------------------------------------------------
    integer :: nBCs, facePos, neighPos, neighAlign, nScalars
    integer :: nModes
    integer :: iMode
    integer :: iBC, iDir, iFace, iAlign
    real(kind=rk), allocatable :: faceOp(:,:)
    real(kind=rk) :: dx
    real(kind=rk) :: sigmod(size(statedata%state,3))
    real(kind=rk) :: bndBaryCoord(1:3)
    real(kind=rk) :: corrector(size(statedata%state,3))
    real(kind=rk), allocatable :: faceMaterial(:,:)
    integer :: sidefact, signfact
    ! ---------------------------------------------------------------------------

    if(present(material)) then
      allocate( faceMaterial( size(material%leftElemMaterialDat,2), &
                            & size(material%leftElemMaterialDat,3)) )
    end if

    nScalars = equation%varSys%nScalars

    nBCs = boundary%nBCs

    ! The length of an element
    dx = mesh%length

    ! There is just one direction in 1D
    iDir = 1

    !>@todo add other variables to private if necessary
    !!!!OMP PARALLEL &
    !!!!OMP PRIVATE(iBC, iDir, iAlign, iFace) &
    !!!!OMP DEFAULT(shared)

    !  facedata%faceRep(iDir)%dat( nFaces, nFaceDoFs, nScalars, 2 )
    allocate( faceOp(size(facedata%faceRep(iDir)%dat,2), &
      &              size(facedata%faceRep(iDir)%dat,3)) )

    ! Iterate over all the boundaries and set the right face values for
    ! the boundaries on all relevant faces.
    do iBC = 1, nBCs
      ! Compute the number of modes to use for the boundary extrapolation of
      ! Neumann boundaries.
      nModes = ceiling(size(statedata%state,1) * bc(iBC)%neumann_mode_fraction)

      ! Now, we iterate over all the faces with this boundary condition and
      ! set the corresponding face values
      do iAlign = 1,2
        neighAlign = tem_invFace_map(iAlign)
        sidefact = (-1)**neighAlign

        ! iDir is x, y, z direction
        ! iAlign is left or right face
        do iFace = 1,boundary%bnd(iBC)%faces(iDir, iAlign)%facePos%nVals


          ! Create the modal representation on the face for the current
          ! face. We need the modal representation of the neighboring (to the
          ! face) fluid element for that.
          neighPos = boundary%bnd(iBC)%faces(iDir, iAlign)%neighPos%val(iFace)

          if (bc(iBC)%enforce_zero_grad) then

            ! faceOp is the state that we tell the boundary condition is
            ! present at the face.
            ! Initially set this to the first mode as this will be used if only
            ! one mode is to be used for the extrapolation.
            faceOp = statedata%state(neighPos,1:1,:)

            if (nModes > 1) then
              ! Only need to actually compute something if we are to use more
              ! then just the first mode for the extrapolation

              ! signfact is used to keep track of the alternating sign on the
              ! left end of the element.
              ! (sidefact decides whether we are left (=-1) or right(=1))
              signfact = sidefact

              ! To enforce a zero gradient we compute the last mode
              ! (at nModes) to match all previous modes contributions in the
              ! derivative. This "correcting" mode is stored in corrector.
              corrector = 0.0_rk

              ! By taking the derivative all modes are shifted down by 1.
              ! Thus, the relevant modes are those up to nModes-1. However,
              ! the one we want to set to enforce the 0 gradient, is one in
              ! nModes-1. Thus, we need to sum up to nModes-2...
              do iMode=1,nModes-2
                sigmod = signfact*statedata%state(neighpos,iMode+1,:)
                ! The derivative contribution of iMode+1, we add it to
                ! the corrector to obtain a corrector, which will balance out
                ! all other modes contributions to the derivative.
                corrector = corrector                           &
                  &       + ( sidefact*iMode*(iMode+1) ) * sigmod

                ! At the same time we can compute the value at the face by
                ! summing the modes.
                faceOp(1,:) = faceOp(1,:) + sigmod
                signfact = signfact * sidefact
              end do

              ! The last mode is not used, instead we use the corrector as
              ! last mode. As that was computed in the derivative, it needs to
              ! "integrated" (divison by the factor for the Legendre derivative)
              faceOp(1,:) = faceOp(1,:) - ( (sidefact * corrector) &
                &                           / (nModes*(nModes-1))  )
            end if

          else

            ! If we do not enforce a zero gradient in any way, we just copy
            ! the face representation which was already projected beforehand
            ! into faceOp.
            faceOp = facedata%faceRep(iDir)%dat(neighPos,:,:,neighAlign)

          end if

          ! get the barycentric coordinate of the (virtual) boundary element
          bndBaryCoord(iDir) = mesh%bary_coord(neighPos,iDir)
          bndBaryCoord(iDir) = bndBaryCoord(iDir) + sidefact*dx

          ! The position of the face with the current boundary condition
          ! inside the face representation.
          facePos = boundary%bnd(iBC)%faces(iDir, iAlign)%facePos%val(iFace)
          if (present(material)) then
            faceMaterial = material%leftElemMaterialDat(iFace,:,:)
            !>@todo replace by subroutine call for OpenMP
            facedata%faceRep(iDir)%dat(facePos,:,:,iAlign)                     &
              & = atl_modg_1d_bnd( bc            = bc(iBC),                    &
              &                    faceOp        = faceOp,                     &
              &                    poly_proj     = poly_proj,                  &
              &                    normalRot     = equation%varRotation(iDir), &
              &                    equation      = equation,                   &
              &                    tree          = tree,                       &
              &                    isNodalScheme = equation%isNonlinear,       &
              &                    time          = time,                       &
              &                    faceDir       = iDir,                       &
              &                    leftOrRight   = neighAlign,                 &
              &                    bndBaryCoord  = bndBaryCoord,               &
              &                    elemLength    = dx,                         &
              &                    faceMaterial  = faceMaterial                )
          else
            !>@todo replace by subroutine call for OpenMP
            facedata%faceRep(iDir)%dat(facePos,:,:,iAlign)                     &
              & = atl_modg_1d_bnd( bc            = bc(iBC),                    &
              &                    faceOp        = faceOp,                     &
              &                    poly_proj     = poly_proj,                  &
              &                    normalRot     = equation%varRotation(iDir), &
              &                    equation      = equation,                   &
              &                    tree          = tree,                       &
              &                    isNodalScheme = equation%isNonlinear,       &
              &                    time          = time,                       &
              &                    faceDir       = iDir,                       &
              &                    leftOrRight   = neighAlign,                 &
              &                    bndBaryCoord  = bndBaryCoord,               &
              &                    elemLength    = dx                          )
          end if
        end do ! iFace
      end do ! iAlign
    end do ! iBC

    !!!!OMP END PARALLEL

    deallocate(faceOp)

  end subroutine atl_modg_1d_set_bnd

  !> Subroutine to create the modal representation for a ceratin boundary face.
  function atl_modg_1d_bnd( bc, faceOp, poly_proj, equation, tree, normalRot, &
    &                       isNodalScheme, time, faceDir, leftOrRight,        &
    &                       bndBaryCoord, elemLength, facematerial )          &
    & result( modalFace )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> The parameters for projection method.
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> Rotation indices to rotate global coordinate system into face normal
    !! coordinate system.
    type(coordRotation_type), intent(in) :: normalRot
    !> The modal representation on the boundary face.
    real(kind=rk) :: modalFace(1, equation%varSys%nScalars)
    !> Does the solver require isNodalScheme information anyway?
    logical, intent(in) :: isNodalScheme
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> The spatial direction of the boundary face, i.e.: \n
    !! 1 -> x direction
    integer, intent(in) :: faceDir
    !> Is left or right of the fluid element a boundary face.
    integer, intent(in) :: leftOrRight
    !> The barycentric boundary element coordinates
    real(kind=rk), intent(in) :: bndBaryCoord(1:3)
    !> The length of an element on the current level
    real(kind=rk), intent(in) :: elemLength
    !> The material of the boundary face.
    !! First dimension is the number of points on the face.
    !! Second dimension is the number of material parameters.
    real(kind=rk), intent(in), optional :: faceMaterial(:,:)
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    if (isNodalScheme) then
      ! For nonlinear, nodal schemes the boundary conditions will be imposed in a
      ! pointwise way.
      modalFace = modg_1d_nodal_bnd( bc           = bc,           &
        &                            faceOp       = faceOp,       &
        &                            poly_proj    = poly_proj,    &
        &                            equation     = equation,     &
        &                            tree         = tree,         &
        &                            normalRot    = normalRot,    &
        &                            time         = time,         &
        &                            faceDir      = faceDir,      &
        &                            leftOrRight  = leftOrRight,  &
        &                            bndBaryCoord = bndBaryCoord, &
        &                            elemLength   = elemLength,   &
        &                            facematerial = faceMaterial  )
    else
      ! For linear equations it could be possible to impose boundary conditions
      ! in modal space directly.
      modalFace = modg_1d_modal_bnd( bc           = bc,           &
        &                            faceOp       = faceOp,       &
        &                            poly_proj    = poly_proj,    &
        &                            equation     = equation,     &
        &                            tree         = tree,         &
        &                            normalRot    = normalRot,    &
        &                            time         = time,         &
        &                            faceDir      = faceDir,      &
        &                            leftOrRight  = leftOrRight,  &
        &                            bndBaryCoord = bndBaryCoord, &
        &                            elemLength   = elemLength    )
    end if

  end function atl_modg_1d_bnd


  !> Set boundary values in a nodal way
  function modg_1d_nodal_bnd( bc, faceOp, poly_proj, equation, tree, &
    &                         normalRot, time, faceDir, leftOrRight, &
    &                         bndBaryCoord, elemLength, facematerial )     &
    &        result( modalFace )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> The parameters for projection method.
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> Rotation indices to rotate global coordinate system into face normal
    !! coordinate system.
    type(coordRotation_type), intent(in) :: normalRot
    !> The modal representation on the boundary face.
    real(kind=rk) :: modalFace(1, equation%varSys%nScalars)
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> The spatial direction of the boundary face, i.e.: \n
    !! 1 -> x direction
    integer, intent(in) :: faceDir
    !> Is left or right of the fluid element a boundary face.
    integer, intent(in) :: leftOrRight
    !> The barycentric boundary element coordinates
    real(kind=rk), intent(in) :: bndBaryCoord(1:3)
    !> The length of an element on the current level
    real(kind=rk), intent(in) :: elemLength
    !> The material of the boundary face.
    !! First dimension is the number of points on the face.
    !! Second dimension is the number of material parameters.
    real(kind=rk), intent(in), optional :: faceMaterial(:,:)
    ! ---------------------------------------------------------------------------
    integer :: iBcVar, bcIndex, iVar
    real(kind=rk), allocatable :: pointValOp(:,:), pointFace(:,:), tmpFace(:)
    ! ---------------------------------------------------------------------------

    allocate(pointFace(1, equation%varSys%nScalars))
    allocate(pointValOp(1,equation%varSys%nScalars))
    allocate(tmpFace(1))

    modalFace(:,:) = 0.0_rk

    do iVar = 1, equation%varSys%nScalars
      pointValOp(1,iVar) = faceOp(1, iVar)
    end do

    ! take care of bc_trafo, i.e. transform from conservative to primitive.
    if(.not.bc%bc_trafo%identity) then
      if (present(faceMaterial)) then
        call bc%bc_trafo%to( equation = equation,    &
          &                  instate  = pointValOp,  &
          &                  material = faceMaterial )
      else
        call bc%bc_trafo%to( equation = equation, instate = pointValOp )
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

      ! Create the modal representation on the face.
      select case(bc%state(iBcVar)%style)
      case('neumann')
        ! In case of dirichlet boundary condition we would like
        ! to prescribe the derivative of the boundary value, so:
        ! Extrapolate the values
        ! @todo replace by subroutine call, due to OpenMP
        pointFace(:,bcIndex) = modg_1d_bnd_extrapolate(          &
          &                      faceOp =  pointValOp(:,bcIndex) )
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
          pointFace(:,bcIndex) = modg_1d_bnd_mirrorPoint(                &
            &                      bc           = bc%state(iBcVar),      &
            &                      poly_proj    = poly_proj,             &
            &                      time         = time,                  &
            &                      varSys       = equation%varSys,       &
            &                      tree         = tree,                  &
            &                      faceDir      = faceDir,               &
            &                      leftOrRight  = leftOrRight,           &
            &                      bndBaryCoord = bndBaryCoord,          &
            &                      elemLength   = elemLength             )
        end if
      case default
        write(logUnit(1),*) 'ERROR in modg_1d_modalBndRep: unknown bnd style, stopping ...'
        call tem_abort()
      end select


    end do

    ! take care of bc_trafo, i.e. transform from primitive to conservative.
    if (.not.bc%bc_trafo%identity) then
      if (present(faceMaterial)) then
        call bc%bc_trafo%from( equation = equation,    &
          &                    instate  = pointFace,   &
          &                    material = faceMaterial )
      else
        call bc%bc_trafo%from( equation = equation, instate = pointFace )
      end if
    end if

    ! transform back to modal space if necessary
    do iVar = 1, equation%varSys%nScalars
      modalFace(1, iVar) = pointFace(1,iVar)
    end do

  end function modg_1d_nodal_bnd


  !> Set boundary values in a modal way
  function modg_1d_modal_bnd( bc, faceOp, poly_proj, equation, tree, &
                            & normalRot, time, faceDir, leftOrRight, &
                            & bndBaryCoord, elemLength )             &
                            & result( modalFace )
    ! ---------------------------------------------------------------------------
    !> The boundary condition to generate the modal representation for.
    type(atl_boundary_type), intent(in) :: bc
    !> The modal representation on the face of the neighboring element.
    real(kind=rk), intent(inout) :: faceOp(:,:)
    !> The parameters for projection method.
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> Rotation indices to rotate global coordinate system into face normal
    !! coordinate system.
    type(coordRotation_type), intent(in) :: normalRot
    !> The modal representation on the boundary face.
    real(kind=rk) :: modalFace(1, equation%varSys%nScalars)
    !> The absolute time point.
    type(tem_time_type), intent(in) :: time
    !> The spatial direction of the boundary face, i.e.: \n
    !! 1 -> x direction
    integer, intent(in) :: faceDir
    !> Is left or right of the fluid element a boundary face.
    integer, intent(in) :: leftOrRight
    !> The barycentric boundary element coordinates
    real(kind=rk), intent(in) :: bndBaryCoord(1:3)
    !> The length of an element on the current level
    real(kind=rk), intent(in) :: elemLength
    ! ---------------------------------------------------------------------------
    integer :: iBcVar, bcIndex
    type(tem_st_fun_listElem_type), pointer :: fPtr
    ! ---------------------------------------------------------------------------

    modalFace(:,:) = 0.0_rk

    ! take care of bc_trafo, i.e. transform from conservative to primitive.
    if(.not.bc%bc_trafo%identity) then
      write(logUnit(1),*) 'ERROR in modg_1d_modalBndRep: not able to convert modal ' // &
        & 'face representation by to desired face values, stopping ...'
      call tem_abort()
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
        modalFace(:,bcIndex) = modg_1d_bnd_extrapolate(     &
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

        else if ( fPtr%val(1)%fun_kind == 'const' ) then
          ! if spacetime function of the user variable
          ! is constant then set modalFace(1) to constant
          modalFace(:,bcIndex) = modg_1d_bnd_mirrorModalConst(       &
            &                      st_Fun        = fPtr%val(1)       )

        else

          ! The function is NOT constant, we transfer to point values,
          ! mirror pointwise and transfer back to a modal representation.
          ! (trivial in 1D)

          ! Mirror pointwise
          modalFace(:,bcIndex) = modg_1d_bnd_mirrorPoint( &
            &              bc           = bc%state(iBcVar), &
            &              poly_proj    = poly_proj,        &
            &              time         = time,             &
            &              varSys       = equation%varSys,  &
            &              tree         = tree,             &
            &              faceDir      = faceDir,          &
            &              leftOrRight  = leftOrRight,      &
            &              bndBaryCoord = bndBaryCoord,     &
            &              elemLength   = elemLength        )

        end if

      case default
        write(logUnit(1),*) 'ERROR in modg_1d_modalBndRep: unknown bnd style, stopping ...'
        call tem_abort()
      end select


    end do

    ! take care of bc_trafo, i.e. transform from primitive to conservative.
    if(.not.bc%bc_trafo%identity) then
      write(logUnit(1),*) 'ERROR in modg_1d_modalBndRep: not able to convert back modal ' // &
        & 'face representation by to desired face values, stopping ...'
      call tem_abort()
    end if

  end function modg_1d_modal_bnd


  !> Function to extrapolate face values for a given boundary condition in
  !! physical or modal space.
  function modg_1d_bnd_extrapolate( faceOp ) result(modalFace)
    ! ---------------------------------------------------------------------------
    !> Modal representation on the face of the neighboring element.
    real(kind=rk), intent(in) :: faceOp(:)
    !> The extrapolated modal representation.
    real(kind=rk) :: modalFace(1)
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    ! Just extrapolate the values on the face, which means: copy the
    ! modal representation.
    modalFace = faceOp

  end function modg_1d_bnd_extrapolate



  !> Function to mirror a modal representation around a given boundary condition
  !! in modal space.
  function modg_1d_bnd_mirrorModalConst( st_fun ) result(modalFace)
    ! ---------------------------------------------------------------------------
    !> List of constant spacetime functions
    type(tem_spacetime_fun_type), intent(in) :: st_fun
    !> The mirrored modal representation.
    real(kind=rk) :: modalFace(1)
    ! ---------------------------------------------------------------------------

    modalFace(1) = st_fun%const(1)

  end function modg_1d_bnd_mirrorModalConst


  !> Function to mirror pointvalues for a given boundary conditions.
  function modg_1d_bnd_mirrorPoint( bc, poly_proj, time, varSys, tree,  &
    &                               faceDir, leftOrRight, bndBaryCoord, &
    &                               elemLength ) result(pointFace)
    ! ---------------------------------------------------------------------------
    !> The boundary state.
    type(tem_bc_state_type), intent(in) :: bc
    !> Data for the projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The current absolute time.
    type(tem_time_type), intent(in) :: time
    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The mirrored isNodalScheme representation.
    real(kind=rk) :: pointFace(1)
    !> The spatial direction of the boundary face, i.e.: \n
    !! 1 -> x direction
    integer, intent(in) :: faceDir
    !> Is left or right of the fluid element a boundary face.
    integer, intent(in) :: leftOrRight
    !> The barycentric boundary element coordinates
    real(kind=rk) :: bndBaryCoord(1:3)
    !> The element length on the current refinement level
    real(kind=rk) :: elemLength
    ! ---------------------------------------------------------------------------
    real(kind=rk), allocatable :: bndVal(:), bndPhysCoord(:,:), bndCoords(:,:)
    ! ---------------------------------------------------------------------------

    allocate( bndVal(1) )
    allocate( bndPhysCoord(1,1:3) )

    ! get the quadrature points on the boundary faces
    call ply_get_quadpoints_faces_1d( poly_proj = poly_proj,   &
      &                               iDir      = faceDir,     &
      &                               ialign    = LeftOrRight, &
      &                               points    = bndCoords    )

    ! Move the Chebyshev coordinates from the reference element to the
    ! physical element.
    call atl_refToPhysCoord( refpoints  = bndCoords,    &
      &                      nPoints    = 1,            &
      &                      baryCoord  = bndBaryCoord, &
      &                      elemLength = elemLength,   &
      &                      physPoints = bndPhysCoord  )

    ! the boundary values at the Chebyshev nodes
    !>@todo This can now be used within a OpenMP parallel block, using
    ! solver%conf(iThread).
    call varSys%method%val(bc%varPos)%get_point( varSys = varSys,       &
      &                                          point  = bndPhysCoord, &
      &                                          time   = time,         &
      &                                          tree   = tree,         &
      &                                          nPnts  = 1,            &
      &                                          res    = bndVal        )

    pointFace(:) = bndVal(:)

    deallocate( bndVal )
    deallocate( bndPhysCoord )
    deallocate( bndCoords )

  end function modg_1d_bnd_mirrorPoint

end module atl_modg_1d_bnd_module


