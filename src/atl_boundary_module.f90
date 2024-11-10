! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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

!> Module to collect all informations about boundary conditions.
!! author: Jens Zudrop
!!
!! We virtually create a "boundary element" outside each boundary face.
!! They get indices larger than nElems on the respective level.
!! These are then used in the face descriptions to refer to the
!! "outside" element (we put the "boundary element" index into the
!! corresponding leftpos or rightpos of the face with a boundary side).
module atl_boundary_module
  use env_module,                   only: rk

  use tem_aux_module,               only: tem_abort
  use tem_faceData_module,          only: tem_face_type, tem_invFace_map
  use tem_element_module,           only: eT_fluid
  use tem_bc_prop_module,           only: tem_bc_prop_type
  use tem_grow_array_module,        only: grw_intArray_type,     &
    &                                     init,                  &
    &                                     append
  use tem_bc_header_module,         only: tem_bc_header_type
  use tem_logging_module,           only: logUnit
  use tem_element_module,           only: eT_fluid
  use treelmesh_module,             only: treelmesh_type
  use tem_varSys_module,            only: tem_varSys_type

  use ply_poly_project_module,      only: ply_poly_project_type

  use aotus_module,                 only: flu_State

  use atl_equation_module,          only: atl_equations_type
  use atl_bc_header_module,         only: atl_boundary_type, atl_load_bc
  use atl_scheme_module,            only: atl_scheme_type,        &
    &                                     atl_modg_scheme_prp,    &
    &                                     atl_modg_2d_scheme_prp, &
    &                                     atl_modg_1d_scheme_prp
  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_reference_element_module, only: atl_refToPhysCoord

  implicit none

  private

  !> Facewise description of the boundaries.
  type atl_faceBnd_type
    !> Position of the faces for which the boundary condition has to be
    !! applied.
    type(grw_intArray_type) :: facePos
    !> Position of the neighboring fluid element.
    type(grw_intArray_type) :: neighPos
  end type atl_faceBnd_type

  !> Description of a certain boundary condition.
  type atl_bndDesc_type
    !> Facewise description of the boundaries.
    !! First dimension is 3 for the three spatial directions.
    !! Second dimension is 2 for left and right faces.
    type(atl_faceBnd_type) :: faces(3,2)
  end type atl_bndDesc_type

  !> type to represent the different boundary conditions on
  !! the same refinement level
  type atl_level_boundary_type
    !> The number of boundary conditions in this description.
    integer :: nBCs
    !> The boundary description. One for each different kind of
    !! boundary. To access this array, use the boundary id as index.
    type(atl_bndDesc_type), allocatable :: bnd(:)
    !> Postition of individual projection method in the projection list
    integer :: poly_proj_pos
  end type atl_level_boundary_type

  public :: atl_level_boundary_type, atl_bndDesc_type
  public :: atl_init_bndList, atl_init_elem_bnd, atl_get_numBndElems
  public :: atl_fill_BCIndex

contains

  ! ************************************************************************ !
  !> Get the number of (virtual) boundary elements for each level and each
  !! direction.
  function atl_get_numBndElems( minLevel, maxLevel, boundary_list ) &
      & result(nBndElems)
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: minLevel !< minimal level in the tree
    integer, intent(in) :: maxLevel !< maximal level in the tree

    !> Boundary description for all levels.
    type(atl_level_boundary_type), intent(in) :: &
      & boundary_list(minLevel:maxLevel)

    !> The number of (virtual) boundary elements for each level and each
    !! spatial direction
    integer :: nBndElems(minLevel:maxLevel,1:3)
    ! -------------------------------------------------------------------- !
    integer :: iLevel, iBc, iDir, lb_bc, ub_bc
    ! -------------------------------------------------------------------- !

    nBndElems(:,:) = 0
    do iLevel = minLevel, maxLevel

      lb_bc = lbound(boundary_list(iLevel)%bnd,1)
      ub_bc = ubound(boundary_list(iLevel)%bnd,1)

      do iBc = lb_bc, ub_bc
        do iDir = 1,3
          nBndElems(iLevel, iDir) = nBndElems(iLevel, iDir)                &
            & + boundary_list(iLevel)%bnd(iBc)%faces(iDir,1)%facePos%nVals &
            & + boundary_list(iLevel)%bnd(iBc)%faces(iDir,2)%facePos%nVals
        end do
      end do

    end do

  end function atl_get_numBndElems
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Creates boundary informations for the faces.
  !!
  !! Used for covolumes.
  subroutine atl_init_elem_bnd( minLevel, maxLevel, bc_header, face_list, &
    &                           boundary_list, scheme_list                )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    !> The boundary condition header data as given in the mesh.
    type(tem_bc_header_type), intent(in) :: bc_header
    !> Description of the faces
    type(tem_face_type), intent(inout) :: face_list(minLevel:maxLevel)
    !> The created boundary information.
    type(atl_level_boundary_type), intent(inout) :: &
      & boundary_list(minLevel:maxLevel)
    !> Scheme information
    type(atl_scheme_type), intent(in) :: scheme_list(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer, allocatable :: nBndElems(:,:)
    integer :: iLevel, iDir, iElem, iBc, iAlign
    integer :: nElems, neighPos, bndElemPos, bc_id, dimen, opSide
    ! -------------------------------------------------------------------- !


    ! Iterate over all the levels and correct the boundary faces in
    ! the description of the faces. Additionally, the boundary elements have to
    ! be created.
    allocate(nBndElems(minLevel:maxLevel,3))
    nBndElems(:,:) = 0
    do iLevel = minLevel, maxLevel

      ! Check the number of polynmial degrees of freedom per spatial direction.
      select case(scheme_list(iLevel)%scheme)
      case(atl_modg_scheme_prp)
        dimen = 3
      case(atl_modg_2d_scheme_prp)
        dimen = 2
      case(atl_modg_1d_scheme_prp)
        dimen = 1
      case default
        call tem_abort( 'ERROR in atl_init_elem_bnd: not able to generate bnd' &
          & // 'info for this scheme, stopping ...'                            )
      end select

      ! set the number of boundary conditions we have.
      boundary_list(iLevel)%nBCs = bc_header%nBCs
      allocate( boundary_list(iLevel)%bnd(boundary_list(iLevel)%nBCs) )

      do iDir = 1, dimen

        ! Init the growing array for the face positions.
        do iBc = 1, boundary_list(iLevel)%nBCs
          call init( me     = boundary_list(iLevel)%bnd(iBc)      &
            &                                      %faces(iDir,1) &
            &                                      %facePos,      &
            &        length = 16                                  )
          call init( me     = boundary_list(iLevel)%bnd(iBc)      &
            &                                      %faces(iDir,2) &
            &                                      %facePos,      &
            &        length = 16                                  )
          call init( me     = boundary_list(iLevel)%bnd(iBc)      &
            &                                      %faces(iDir,1) &
            &                                      %neighPos,     &
            &        length = 16                                  )
          call init( me     = boundary_list(iLevel)%bnd(iBc)      &
            &                                      %faces(iDir,2) &
            &                                      %neighPos,     &
            &        length = 16                                  )
        end do

        ! Get the number of elements on the current level in the
        ! dimension-by-dimension element description of the faces.
        nElems = face_list(iLevel)%dimByDimDesc(iDir)%nElems

        do iElem = 1,face_list(iLevel)%dimByDimDesc(iDir)%elem%nElems(eT_fluid)

          do iAlign = 1, 2

            neighPos = face_list(iLevel)%dimByDimDesc(iDir)   &
              &                         %neigh(1)             &
              &                         %nghElems(iAlign,iElem)

            ! Is element on this side (iAlign) a boundary?
            ! If neighbor is a boundary, then neighPos is negative and -BC_ID
            if ( neighPos < 0 ) then
              ! This element has a boundary, increment the counter for
              ! boundary elements.
              nBndElems(iLevel, iDir) = nBndElems(iLevel, iDir) + 1

              ! Set the correct boundary element position in the neighbor array.
              bndElemPos = nBndElems(iLevel, iDir) + nElems
              face_list(iLevel)%dimByDimDesc(iDir)     &
                &              %neigh(1)               &
                &              %nghElems(iAlign,iElem) &
                & = bndElemPos

              ! Get BC ID
              bc_id = (-1) * neighPos
              ! the outside face is in the opposite side of current side
              opSide = tem_invFace_map( iAlign )
              ! append the face position to the boundary description.
              call append( me  = boundary_list(iLevel)%bnd(bc_id)         &
                &                                     %faces(iDir,opSide) &
                &                                     %facePos,           &
                &          val = bndElemPos                               )
              ! Set neighPos to myself
              call append( me  = boundary_list(iLevel)%bnd(bc_id)         &
                &                                     %faces(iDir,opSide) &
                &                                     %neighPos,          &
                &          val = iElem                                    )
            end if
          end do ! iAlign = 1, 2

        end do
      end do
    end do

  end subroutine atl_init_elem_bnd
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to create the levelwise list of boundaries.
  subroutine atl_init_bndList( conf, tree, equation, bc_prop, face_list, &
    &                          boundary_list, bc, bc_header, scheme_list )
    ! -------------------------------------------------------------------- !
    !> Lua script to obtain the configuration data from.
    type(flu_State) :: conf
    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree
    !> The equation to initialize the boundary info for.
    type(atl_equations_type), intent(inout) :: equation
    !> The boundary properties.
    type(tem_bc_prop_type), intent(in) :: bc_prop
    !> The description of the faces as provided by the [[tem_face_module]]
    type(tem_face_type), intent(inout) ::                  &
      & face_list(tree%global%minLevel:tree%global%maxLevel)
    !> The created boundary type.
    type(atl_level_boundary_type), intent(inout) ::                  &
      & boundary_list(tree%global%minLevel:tree%global%maxLevel)
    !> Global description of all boundaries
    type(atl_boundary_type), allocatable, intent(out) :: bc(:)
    !> The boundary condition header data as given in the mesh.
    type(tem_bc_header_type), intent(out) :: bc_header
    !> The scheme you are using.
    type(atl_scheme_type), intent(in) ::                &
      & scheme_list(tree%global%minLevel:tree%global%maxLevel)
    ! -------------------------------------------------------------------- !

    ! Load the boundary description from the lua configuration file
    call atl_load_bc( bc        = bc,       &
      &               bc_header = bc_header,&
      &               bc_prop   = bc_prop,  &
      &               equation  = equation, &
      &               conf      = conf      )

    ! create the boundary informations for the faces
    call atl_init_face_bnd( minLevel      = tree%global%minLevel, &
      &                     maxLevel      = tree%global%maxLevel, &
      &                     bc_header     = bc_header,            &
      &                     face_list     = face_list,            &
      &                     boundary_list = boundary_list,        &
      &                     scheme_list   = scheme_list           )

  end subroutine atl_init_bndList
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Creates boundary informations for the faces.
  !!
  !! For each boundary face we create a virtual boundary element outside
  !! the fluid domain.
  !! Those virtual elements are indexed with an offset of the number of
  !! elements in the level.
  !! And these indices are than used for the faces to refer to this virtual
  !! element as a neighbor on the corresponding side.
  subroutine atl_init_face_bnd( minLevel, maxLevel, bc_header, face_list, &
    &                           boundary_list, scheme_list                )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    !> The boundary condition header data as given in the mesh.
    type(tem_bc_header_type), intent(in) :: bc_header
    !> Description of the faces
    type(tem_face_type), intent(inout) :: face_list(minLevel:maxLevel)
    !> The created boundary information.
    type(atl_level_boundary_type), intent(inout) :: &
      & boundary_list(minLevel:maxLevel)
    !> Scheme information
    type(atl_scheme_type), intent(in) :: scheme_list(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer, allocatable :: nBndElems(:,:)
    integer :: iLevel, iDir, iFace, iBc
    integer :: nElems, leftPos, rightPos, bc_id, dimen
    integer :: bndElemPos
    ! -------------------------------------------------------------------- !


    ! Iterate over all the levels and correct the boundary faces in
    ! the description of the faces. Additionally, the boundary elements have to
    ! be created.
    allocate(nBndElems(minLevel:maxLevel,3))
    nBndElems(:,:) = 0
    do iLevel = minLevel, maxLevel

      ! Check the number of polynmial degrees of freedom per spatial direction.
      select case(scheme_list(iLevel)%scheme)
      case(atl_modg_scheme_prp)
        dimen = 3
      case(atl_modg_2d_scheme_prp)
        dimen = 2
      case(atl_modg_1d_scheme_prp)
        dimen = 1
      case default
        call tem_abort( 'ERROR in atl_init_face_bnd: not able to generate bnd' &
          & // ' info for this scheme, stopping ...'                           )
      end select

      ! set the number of boundary conditions we have.
      boundary_list(iLevel)%nBCs = bc_header%nBCs
      allocate( boundary_list(iLevel)%bnd(boundary_list(iLevel)%nBCs) )

      do iDir = 1, dimen

        ! Init the growing array for the face positions.
        do iBc = 1, boundary_list(iLevel)%nBCs
          call init(                                                         &
            & me     = boundary_list(iLevel)%bnd(iBc)%faces(iDir,1)%facePos, &
            & length = 16                                                    )
          call init(                                                         &
            & me     = boundary_list(iLevel)%bnd(iBc)%faces(iDir,2)%facePos, &
            & length = 16                                                    )
          call init(                                                          &
            & me     = boundary_list(iLevel)%bnd(iBc)%faces(iDir,1)%neighPos, &
            & length = 16                                                     )
          call init(                                                          &
            & me     = boundary_list(iLevel)%bnd(iBc)%faces(iDir,2)%neighPos, &
            & length = 16                                                     )
        end do

        ! Get the number of elements on the current level in the
        ! dimension-by-dimension element description of the faces.
        nElems = face_list(iLevel)%dimByDimDesc(iDir)%nElems

        ! Iterate over all computed faces.
        do iFace = 1, size(face_list(iLevel)%faces(iDir)%computeFace%leftPos)

          leftPos = face_list(iLevel)%faces(iDir)%computeFace%leftPos(iFace)
          rightPos = face_list(iLevel)%faces(iDir)%computeFace%rightPos(iFace)

          ! Is the left side of this face a boundary?
          if (leftPos < 0) then
            ! We count the number of boundary elements we have to create and
            ! offset the index by the total number of elements on the level.
            nBndElems(iLevel, iDir) = nBndElems(iLevel, iDir) + 1
            bndElemPos = nBndElems(iLevel, iDir) + nElems
            ! There is a boundary left to this face, point to the (virtual)
            ! boundary element as left neighboring element.
            face_list(iLevel)%faces(iDir)%computeFace%leftPos(iFace) &
              & = bndElemPos
            ! append the boundary element position to the boundary description.
            bc_id = (-1) * leftPos
            ! The boundary face is right of the virtual boundary element.
            call append(                                                      &
              & me  = boundary_list(iLevel)%bnd(bc_id)%faces(iDir,2)%facePos, &
              & val = bndElemPos                                              )
            ! The element right of the face is the neighboring element for this
            ! boundary element
            call append(                                                       &
              & me  = boundary_list(iLevel)%bnd(bc_id)%faces(iDir,2)%neighPos, &
              & val = rightPos                                                 )
          end if

          ! Is the right side of this face a boundary?
          if (rightPos < 0) then
            ! We count the number of boundary elements we have to create and
            ! offset the index by the total number of elements on the level.
            nBndElems(iLevel, iDir) = nBndElems(iLevel, iDir) + 1
            bndElemPos = nBndElems(iLevel, iDir) + nElems
            ! There is a boundary right to this face, point to the (virtual)
            ! boundary element as right neighboring element.
            face_list(iLevel)%faces(iDir)%computeFace%rightPos(iFace) = &
              & bndElemPos
            ! append the boundary element position to the boundary description.
            bc_id = (-1) * rightPos
            ! The boundary face is left of the virtual boundary element.
            call append(                                                      &
              & me  = boundary_list(iLevel)%bnd(bc_id)%faces(iDir,1)%facePos, &
              & val = bndElemPos                                              )
            ! The element left of the face is the neighboring element for this
            ! boundary element
            call append(                                                       &
              & me  = boundary_list(iLevel)%bnd(bc_id)%faces(iDir,1)%neighPos, &
              & val = leftPos                                                  )
          end if

        end do
      end do
    end do

  end subroutine atl_init_face_bnd
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_fill_BCIndex( tree, bc, boundary_list, nDim, varSys, &
    &                          poly_proj_list, mesh_list              )
    ! -------------------------------------------------------------------- !
    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree
    !> Global description of all boundaries
    type(atl_boundary_type),  intent(inout) :: bc(:)
    !> The created boundary information.
    type(atl_level_boundary_type), intent(in) ::               &
      & boundary_list(tree%global%minLevel:tree%global%maxLevel)
    integer, intent(in) :: nDim
    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys
    !> unique list for projection methods
    type(ply_poly_project_type), target, intent(inout)  :: poly_proj_list(:)
    !> The mesh  you are using.
    type(atl_cube_elem_type), intent(in) :: mesh_list(tree%global%minLevel &
      &                                           :tree%global%maxLevel)
    ! -------------------------------------------------------------------- !
    integer:: iLevel, iVar, iBc
    integer:: nFaces, nQuadPnts
    integer, allocatable :: idx(:)
    real(kind=rk), allocatable :: BCPnts(:,:)
    character, allocatable :: BCoffsetBit(:)
    type(ply_poly_project_type), pointer :: poly_proj
    integer :: varPos
    ! -------------------------------------------------------------------- !
    write(logUnit(3),*) ' Set params in bc variable'
    do iBC = 1, size(bc)
      do iVar = 1, bc(iBC)%varDict%nVals
        varPos = bc(iBC)%state(iVar)%varPos
        call varSys%method%val(varPos)%set_params( &
          & varSys   = varSys,                     &
          & instring = 'isSurface = true'          )
      end do

      do iVar = 1, bc(iBC)%varDict_gradient%nVals
        varPos = bc(iBC)%state_gradient(iVar)%varPos
        call varSys%method%val(varPos)%set_params( &
          & varSys   = varSys,                     &
          & instring = 'isSurface = true'          )
      end do

    end do

    write(logUnit(3),*) ' Setup indices for boundary condition'
    ! run over all levels
    do iLevel = tree%global%minLevel, tree%global%maxLevel

      ! set the correct projection type
      poly_proj => poly_proj_list(boundary_list(iLevel)%poly_proj_pos)

      ! First get the amount of points and faces to allocate arrays accordingly
      select case(nDim)
      case (1)
        nQuadPnts = poly_proj%body_1d%faces(1,1)%nquadPoints
      case (2)
        nQuadPnts = poly_proj%body_2d%faces(1,1)%nquadPoints
      case (3)
        nQuadPnts = poly_proj%body_3d%faces(1,1)%nquadPoints
      end select

      write(logUnit(9),*) 'Nr. Quad points at faces: ', nQuadPnts

      ! run over all boundary conditions
      do iBc = 1,boundary_list(iLevel)%nBCs

        write(LogUnit(10),*) 'for boundary', iBC, &
          & 'number of variables are ', Bc(iBC)%varDict%nVals

        ! number of faces with this boundary
        nFaces = sum(boundary_list(iLevel)%bnd(iBc)%faces(:,:)%facepos%nVals)

        ! if there are no faces, we can jump to the next bc
        if (nFaces == 0 ) cycle

        call atl_get_points_from_BC(                     &
          & bnd        = boundary_list(iLevel)%bnd(iBC), &
          & points     = BCPnts,                         &
          & offset_bit = BCoffsetBit,                    &
          & poly_proj  = poly_proj,                      &
          & mesh       = mesh_list(ilevel),              &
          & nDim       = nDim,                           &
          & nQuadPnts  = nQuadPnts,                      &
          & nFaces     = nFaces                          )
        allocate( idx(nFaces*nQuadPnts) )
        ! run over all variables for this boundary condition
        do iVar = 1, Bc(iBC)%varDict%nVals
          idx = 0
          call varSys%method%val(bc(iBc)%state(iVar)%varPos)%setup_indices( &
            & varSys     = varSys,                                          &
            & point      = BCPnts,                                          &
            & offset_bit = BCoffsetBit,                                     &
            & iLevel     = iLevel,                                          &
            & tree       = tree,                                            &
            & nPnts      = nQuadPnts*nFaces,                                &
            & idx        = idx                                              )

          call append( me  = bc(iBc)%state(iVar)%pntIndex%indexLvl(iLevel), &
            &          val = idx                                            )
        end do !iVar

        ! run over all gradient variables for this boundary condition
        do iVar = 1, Bc(iBC)%varDict_gradient%nVals
          idx = 0
          call varSys%method%val(bc(iBc)%state_gradient(iVar)%varPos) &
            &                           %setup_indices(               &
            & varSys     = varSys,                                    &
            & point      = BCPnts,                                    &
            & offset_bit = BCoffsetBit,                               &
            & iLevel     = iLevel,                                    &
            & tree       = tree,                                      &
            & nPnts      = nQuadPnts*nFaces,                          &
            & idx        = idx                                        )
          call append( me  = bc(iBc)%state_gradient(iVar)%pntIndex &
            &                       %indexLvl(iLevel),             &
            &          val = idx                                   )
        end do !iVar

        deallocate( BCPnts )
        deallocate( BCoffsetBit )
        deallocate( idx )
      end do !iBc

    end do !iLevel

  end subroutine atl_fill_BCIndex
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get all the surface points for a specific boundary.
  subroutine atl_get_points_from_BC( bnd, points, offset_bit, poly_proj, &
    &                                mesh, nDim, nQuadPnts, nFaces )
    ! -------------------------------------------------------------------- !
    type(atl_bndDesc_type), intent(in) :: bnd
    !> array of points on the specific boundray
    real(kind=rk), allocatable, intent(out) :: points(:,:)
    !> offset vector for all points on the boundary, requiered for coupling
    ! transform into character
    character, allocatable, intent(out) :: offset_bit(:)
    !> projection list
    type(ply_poly_project_type), intent(in) :: poly_proj
    !> The mesh  you are using.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> Equation nDimensions
    integer, intent(in) :: nDim
    !> Number of quadrature points on faces
    integer, intent(in) :: nQuadPnts
    !> Number if faces on this boundary
    integer, intent(in) :: nFaces
    ! -------------------------------------------------------------------- !
    integer :: idir, iAlign, iFace
    integer :: neighPos, neighAlign
    integer :: coord(3)
    real(kind=rk) :: bndBaryCoord(1:3)
    real(kind=rk), allocatable :: facePnts(:,:)
    integer :: iglobFace
    ! -------------------------------------------------------------------- !

    ! allocate array for all exchange points and the offset array
    allocate( points (nFaces*nQuadPnts,3) )
    allocate( offset_bit (nFaces*nQuadPnts) )
    ! arrays for the points of this specific iDir,iAlign
    allocate( facePnts (nQuadPnts,3) )

    iglobFace=0
    ! gather the points of that boundary
    ! run over the faces
    do idir = 1, nDim
      coord = 0
      ! run over left and right face
      do iAlign = 1,2

        ! run over the faces of this kind
        do iFace = 1, bnd%faces(iDir,iAlign)%facePos%nVals
          iglobFace = iglobFace + 1

          ! barycenter coordinates of the boundary element.
          ! since elemPos of local boundary element is not known directly,
          ! local barycenter is computed from neighbor barycenter
          neighPos = bnd%faces(iDir,iAlign)%neighPos%val(iFace)
          neighAlign = tem_invFace_map(iAlign)
          bndBaryCoord(1:3) = mesh%bary_coord(neighPos,1:3)
          bndBaryCoord(iDir) = bndBaryCoord(iDir)                      &
            &                    + ((-1.0_rk)**neighAlign) * mesh%length

          ! convert quad points to physical space
          select case (nDim)
          case (3)
            call atl_refToPhysCoord(                                      &
              & refPoints  = poly_proj%body_3d%faces(iDir,iAlign)%points, &
              & nPoints    = nQuadPnts,                                   &
              & baryCoord  = bndBaryCoord,                                &
              & elemLength = mesh%length,                                 &
              & physPoints = facePnts                                     )
          case (2)
            call atl_refToPhysCoord(                                      &
              & refPoints  = poly_proj%body_2d%faces(iDir,iAlign)%points, &
              & nPoints    = nQuadPnts,                                   &
              & baryCoord  = bndBaryCoord,                                &
              & elemLength = mesh%length,                                 &
              & physPoints = facePnts                                     )
          case (1)
            call atl_refToPhysCoord(                                      &
              & refPoints  = poly_proj%body_1d%faces(iDir,iAlign)%points, &
              & nPoints    = nQuadPnts,                                   &
              & baryCoord  = bndBaryCoord,                                &
              & elemLength = mesh%length,                                 &
              & physPoints = facePnts                                     )
          end select
          ! save the points onto the overall points array
          points( (iglobFace-1)*nQuadPnts+1:iglobFace*nQuadPnts,:) &
            & = facePnts(:,:)

          ! the offset bit is the same for all points on this iFace
          ! since musubi is linkwise and this is more general
          ! we transform iDir,iAlign into coord and use that as offset bit
          coord(iDir) = neighAlign*2 - 3
          offset_bit( (iglobFace-1)*nQuadPnts+1:iglobFace*nQuadPnts ) = &
            & achar((coord(1)+1) + (coord(2)+1)*4 + (coord(3)+1)*16 )

        end do !iFaces
      end do !iAlign
    end do !iDir

  end subroutine atl_get_points_from_BC
  ! ************************************************************************ !

end module atl_boundary_module

