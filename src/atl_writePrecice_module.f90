! Copyright (c) 2013-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017-2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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


!> This module includes subroutines, to use different point distribution on the
!! coupling interface, when using preCICE for the coupling. The point
!! distribution can either be equidistant or non-equidistant. Beside that the
!! Nearest-Projection (an interpolation method) can be used with both point
!! distribution.
module atl_writePrecice_module
  use, intrinsic :: iso_c_binding,    only: c_loc, C_F_POINTER

  use env_module,                     only: rk, long_k

  use tem_logging_module,             only: logUnit
  use tem_grow_array_module,          only: append
  use tem_time_module,                only: tem_time_type
  use treelmesh_module,               only: treelmesh_type
  use tem_precice_module,             only: tem_precice_write,            &
    &                                       tem_precice_coupling_type,    &
    &                                       tem_precice_set_edge,         &
    &                                       tem_precice_set_triangle,     &
    &                                       tem_precice_set_vertices_pos
  use tem_precice_interpolation_module,            &
    & only: tem_create_surface_equidistant_points, &
    &       tem_create_edges_triangles
  use tem_geometry_module,            only: tem_CoordOfReal, &
    &                                       tem_BaryOfID
  use tem_topology_module,            only: tem_IdOfCoord
  use tem_spacetime_fun_module,       only: tem_st_fun_linkedList_type, &
   &                                        tem_st_fun_listElem_type
  use tem_varSys_module,              only: tem_varSys_type

  use ply_poly_project_module,        only: ply_poly_project_type, &
    &                                       assignment(=)

  use atl_cube_elem_module,           only: atl_cube_elem_type
  use atl_varSys_module,              only: atl_varSys_data_type
  use atl_reference_element_module,   only: atl_refToPhysCoord

  implicit none

  private

  public :: atl_write_precice
  public :: atl_write_precice_getPoint
  public :: atl_read_precice

contains


  ! ************************************************************************ !
  subroutine atl_nearest_projection(nPointsPerDir, nPnts, nPointsPerFace, &
    &                               scheme_dim, iLevel, preciceCpl        )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: nPointsPerDir
    integer, intent(in) :: nPnts
    integer, intent(in) :: nPointsPerFace
    integer, intent(in) :: scheme_dim
    integer, intent(in) :: iLevel
    integer :: nFaces
    integer :: nTris
    integer :: nEdges
    integer :: iEdge, iFace, iEdge_loc, offset, iTri
    integer, allocatable :: edgeIDs (:)
    type(tem_precice_coupling_type), pointer :: preciceCpl
    ! -------------------------------------------------------------------- !
    call tem_create_edges_triangles(                          &
      & me                = preciceCpl%interpolation          &
      &                               %interpol_data(iLevel), &
      & scheme_dim        = scheme_dim,                       &
      & nQuadPointsPerDir = nPointsPerDir                     )
      write(logUnit(5),*) 'scheme_dim', scheme_dim
      write(logUnit(5),*) 'Number of Points per direction', nPointsPerDir
      nFaces = nPnts / nPointsPerFace
      write(logUnit(5),*) 'nFaces and nPointsPerFace ', nFaces, &
         &                nPointsPerFace
      nEdges = preciceCpl%interpolation%interpol_data(iLevel)%nEdges
      write(logUnit(5),*) 'nEdges', nEdges
      !> For a 2D testcase just the edges are of importance, in 3D also
      !! triangles have to be provided
      associate(                                                             &
        & posIDs    => preciceCpl%writeVar%posIDLvl(iLevel)%posIDs,          &
        & edges     => preciceCpl%interpolation%interpol_data(iLevel)%edges, &
        & triangles => preciceCpl%interpolation%interpol_data(iLevel)        &
        &                                      %triangles                    )
        if (scheme_dim > 1) then
          allocate(edgeIDs(nEdges*nFaces))
          write(logUnit(5),*) 'Setting edges for precice'
         !> loop over all faces and edges. Set an offset, to get from one face
         !! to the other
          iEdge = 0
          do iFace = 1, nFaces
            do iEdge_loc = 1, nEdges
              iEdge = iEdge + 1
              offset = (iFace-1)*nPointsPerFace
              call tem_precice_set_edge(                                   &
                & meshID         = preciceCpl%writeVar%meshID,             &
                & firstVertexID  = posIDs( offset + edges(iEdge_loc, 1) ), &
                & secondVertexID = posIDs( offset + edges(iEdge_loc, 2) ), &
                & edgeID         = edgeIDs( iEdge )                        )
            end do
          end do
          write(logUnit(5),*) 'Done setting edges for precice'
        end if

        !> For a 3D testcase in addition to the edges triangles have to be
        !! created
        nTris = preciceCpl%interpolation%interpol_data(iLevel)%nTriangles
        if (scheme_dim > 2) then
          write(logUnit(5),*) 'Starting setting the EdgeIDs for the triangles'
          do iFace = 1, nFaces
            do iTri = 1, nTris
              offset = (iFace-1)*nEdges
              call tem_precice_set_triangle(                             &
                & meshID       = preciceCpl%writeVar%meshID,             &
                & firstEdgeID  = edgeIDs( offset + triangles(iTri, 1) ), &
                & secondEdgeID = edgeIDs( offset + triangles(iTri, 2) ), &
                & thirdEdgeID  = edgeIDs( offset + triangles(iTri, 3) )  )
            end do
          end do
          write(logUnit(5),*) 'Done setting triangles for precice'
        end if
      end associate
      deallocate( edgeIDs )

  end subroutine atl_nearest_projection
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> ceate equidistant points for write to precice and use the number of these
  !! points for the setup_indices routine. If its set by the user, the
  !! Nearest-Pojection can be used for these points
  subroutine atl_write_equiPoints( stFun, tree, mesh_list, varSys, preciceCpl )
    ! -------------------------------------------------------------------- !
    !> This type contains all information, which are needed for the coupling
    !! with precice
    type(tem_precice_coupling_type), pointer :: preciceCpl
    type(tem_varSys_type), intent(in) :: varSys
    !> This list contains space-time-function
    type(tem_st_fun_ListElem_type), intent(inout), target :: stFun
    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(in)            :: tree
    type(atl_cube_elem_type), intent(in) ::                &
      & mesh_list(tree%global%minLevel:tree%global%maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iLevel, iSerPnt, iCoord, iPnt, iElem, iVar
    integer :: nPnts, nEqPnts, nEqPnts_face, nPnts_face, nElems, varPos
    integer :: iEqPnt
    integer :: pos
    integer :: iAlign, iDira(1), iDir, iEqu
    integer :: offsetX, offsetY, offsetZ
    integer(kind=long_k) :: treeID
    !> list of indices for all coupling points
    integer, allocatable :: idx(:)
    !> Array with all equidistant points
    real(kind=rk), allocatable :: eQuePoints_all(:,:)
    !> list of equidistant points for creating the indices
    real(kind=rk), allocatable :: serializedPnts(:)
    real(kind=rk) :: point(3), baryCoord(3), diff(3), shifted(3)
    real(kind=rk), allocatable :: eQuePoint(:,:)
    type(atl_varSys_data_type), pointer :: fPtr => NULL()
    type(ply_poly_project_type), pointer :: ply_poly => NULL()
    ! -------------------------------------------------------------------- !
    write(logUnit(5),*) 'Using equidistant point distribution for the' &
      & // ' interpolation'
    call C_F_POINTER( stFun%solver_bundle, fPtr )

    do iLevel = tree%global%minLevel, tree%global%maxLevel
      pos = fPtr%solverData%poly_proj_posPtr(iLevel)
      ply_poly => fPtr%solverData%polyProjectPtr(pos)

      !> get the number of unquie points, since in all dircetions the number of
      !! points is equal, just choose one direction (x, y or z)
      nPnts = stFun%pntData%pntLvl(iLevel)%grwPnt%coordX%nVals
      write(logUnit(5),*) 'Total number of Points', nPnts
      !create the equdistant points on the reference element
      call tem_create_surface_equidistant_points(                    &
        & me       = preciceCpl%interpolation%interpol_data(iLevel), &
        & nfactor  = preciceCpl%interpolation,                       &
        & nDir     = fPtr%solverData%equationPtr%nDimensions,        &
        & nPoly    = ply_poly%oversamp_degree                        )

      !> when there are no points, you dont write them to precice
      !! can happen in parallel
      if ( nPnts == 0) then
      !!VK  write(*,*) 'number of Points to writ to precice is 0, &
      !!VK                  cycle this spacetimefunction'
        cycle
      end if
      ! initialize nPnts_face just to silence compiler warning
      nPnts_face = 0
      ! select dimension in to get the right number of Points
      select case (fPtr%solverData%equationPtr%nDimensions)
      case (1)
        nEqPnts_face = preciceCpl%interpolation%interpol_data(iLevel)%nPoints
        nPnts_face = ply_poly%body_1d%faces(1, 1)%nQuadPoints
      case (2)
        nEqPnts_face = preciceCpl%interpolation%interpol_data(iLevel)%nPoints
        nPnts_face = ply_poly%body_2d%faces(1, 1)%nQuadPoints
      case (3)
        nEqPnts_face = preciceCpl%interpolation%interpol_data(iLevel)%nPoints
        nPnts_face = ply_poly%body_3d%faces(1, 1)%nQuadPoints
      end select

      nElems = nPnts/nPnts_face
      nEqPnts = nElems*nEqPnts_face
      allocate( eQuePoint(nEqPnts_face,3))
      allocate( eQuePoints_all(nEqPnts,3))
      allocate( serializedPnts(3*nEqPnts) )
      iSerPnt = 0
      iEqPnt = 0
      do iElem = 1, nElems
        iPnt = (iElem-1)*nPnts_face + 1

        point(1) = stFun%pntData%pntLvl(iLevel)%grwPnt%coordX%val(iPnt)
        point(2) = stFun%pntData%pntLvl(iLevel)%grwPnt%coordY%val(iPnt)
        point(3) = stFun%pntData%pntLvl(iLevel)%grwPnt%coordZ%val(iPnt)
        !> transform offset bit back
        offsetX = mod(ichar(stFun%pntData%pntLvl(iLevel) &
          &                               %offset_bit%val(iPnt)),4) - 1
        offsetY = mod(ichar(stFun%pntData%pntLvl(iLevel) &
          &                              %offset_bit%val(iPnt)),16)/4 - 1
        offsetZ = ichar(stFun%pntData%pntLvl(iLevel)     &
          &                          %offset_bit%val(iPnt))/16 - 1
        !> shift the points according to the offset bit * small fraction
        !! of the boundingCubeLength hence it is mesh independent
        shifted(1) = point(1) &
          & - offsetX*spacing(tree%global%BoundingCubeLength)
        shifted(2) = point(2) &
          & - offsetY*spacing(tree%global%BoundingCubeLength)
        shifted(3) = point(3) &
          & - offsetZ*spacing(tree%global%BoundingCubeLength)


        treeID = tem_IdOfCoord(                               &
          &          tem_CoordOfReal(tree, shifted(:),iLevel) )

        baryCoord = tem_BaryOfId( tree, treeID )

        diff = shifted - baryCoord

        iDira = maxloc(abs(diff))
        iDir = iDira(1)
        if (diff(iDir) > 0.0_rk) then
          iAlign = 2
        else
          iAlign = 1
        end if
        !> We assume here, that all points for a specific face are grouped
        !! together in order to make use of the posIDs

        select case (fPtr%solverData%equationPtr%nDimensions)
        case (1)
          call atl_refToPhysCoord(                                        &
            & refPoints  = preciceCpl%interpolation%interpol_data(iLevel) &
            &                                      %faces(iDir,iAlign)    &
            &                                      %Points,               &
            & nPoints    = nEqPnts_face,                                  &
            & baryCoord  = baryCoord,                                     &
            & elemLength = mesh_list(iLevel)%length,                      &
            & physPoints = eQuePoint                                      )
        case (2)
          call atl_refToPhysCoord(                                        &
            & refPoints  = preciceCpl%interpolation%interpol_data(iLevel) &
            &                                      %faces(iDir,iAlign)    &
            &                                      %Points,               &
            & nPoints    = nEqPnts_face,                                  &
            & baryCoord  = baryCoord,                                     &
            & elemLength = mesh_list(iLevel)%length,                      &
            & physPoints = eQuePoint                                      )
        case (3)
          call atl_refToPhysCoord(                                        &
            & refPoints  = preciceCpl%interpolation%interpol_data(iLevel) &
            &                                      %faces(iDir,iAlign)    &
            &                                      %Points,               &
            & nPoints    = nEqPnts_face,                                  &
            & baryCoord  = baryCoord,                                     &
            & elemLength = mesh_list(iLevel)%length,                      &
            & physPoints = eQuePoint                                      )
        end select

        !> Give the correct exchange points to precice (store it there)
        !! return indices which will be used for all further calls
        do iEqu = 1, nEqPnts_face
          do iCoord = 1,3
            iSerPnt = iSerPnt + 1
            serializedPnts(iSerPnt) = eQuePoint(iEqu,iCoord)
          end do
          iEqPnt = iEqPnt + 1
          eQuePoints_all(iEqPnt,1) = eQuePoint(iEqu,1) &
            & - offsetX*spacing(tree%global%BoundingCubeLength)
          eQuePoints_all(iEqPnt,2) = eQuePoint(iEqu,2) &
            & - offsetY*spacing(tree%global%BoundingCubeLength)
          eQuePoints_all(iEqPnt,3) = eQuePoint(iEqu,3) &
            & - offsetZ*spacing(tree%global%BoundingCubeLength)
        end do !Equ
      end do !iElem
      write(logUnit(5),*) 'Done creating EquiPoints'
      allocate (idx(nEqPnts))

      ! loop over the variables to write to precice
      do iVar  = 1, preciceCpl%writeVar%nVars
        idx = 0
        varPos = preciceCpl%writeVar%varPos(iVar)
        !> setup_indices for boundary variables are read_var in precice
        !! space-time function we are not able to use 'idx' created here
        !! as input for get_valOfIndex for write_vars.Therefore we have
        !! to do setup_indices for write_vars once during initialization
        !! and then store the idx in precice_coupling type.
        !! get the setup_indices for the shifted points
        call varSys%method%val(varPos)%setup_indices( &
          & varSys = varSys,                          &
          & point  = eQuePoints_all,                  &
          & iLevel = iLevel,                          &
          & tree   = tree,                            &
          & nPnts  = nEqPnts,                         &
          & idx    = idx                              )

        !> Since we are sending 5 times the same idx, for the unique
        !! Points and all of them have the same idx, we just need to
        !! store it once
        if (iVar == 1) then
          call append(                                    &
            & me  = preciceCpl%pntIndex%indexLvl(iLevel), &
            & val = idx                                   )
        end if
      end do !iVar

      allocate( preciceCpl%writeVar%posIDLvl(iLevel)%posIDs(nEqPnts) )
      call tem_precice_set_vertices_pos(               &
        & meshID      = preciceCpl%writeVar%meshID,    &
        & nPoints     = nEqPnts,                       &
        & points      = serializedPnts,                &
        & positionID  = preciceCpl%writeVar            &
        &                         %posIDLvl(iLevel)    &
        &                         %posIDs(1:nEqPnts) )
      write(logUnit(5),*) 'Done with precice set_vertices for equidistant &
        &                  Points'
      ! ****************************************************************** !
      !> Create edges and triangles for the interpolation with equidistant
      !! points, if the user wants to make use of it, by setting that in
      !! the ateles config file
      ! ****************************************************************** !
      if ( preciceCpl%interpolation%use_NP) then
        write(logUnit(5),*) 'Using Nearest-Projection with equidistant Points'
        call atl_nearest_projection(                                        &
          & nPointsPerDir  = preciceCpl%interpolation%interpol_data(iLevel) &
          &                                          %nPointsPerDir,        &
          & nPointsPerFace = nEqPnts_face,                                  &
          & nPnts          = nEqPnts,                                       &
          & preciceCpl     = preciceCpl,                                    &
          & scheme_dim     = fPtr%solverData%equationPtr%nDimensions,       &
          & iLevel         = iLevel                                         )
      end if
      deallocate( serializedPnts )
      deallocate( eQuePoint )
      deallocate( eQuePoints_all )
      deallocate( idx )
    end do !iLevel
  end subroutine atl_write_equiPoints
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> setup_indices for boundary variables are read_var in precice
  !! space-time function we are not able to use 'idx' created here as input for
  !! get_valOfIndex for write_vars.Therefore we have to do setup_indices for
  !! write_vars once during initialization and then store the idx in
  !! precice_coupling type. This is done for the non-equidistant points. We also
  !! set the vertices and the edges as well as the triangles for the
  !! Interpolation between the domains.
  subroutine atl_write_nonequiPoints( varSys, tree, stFun, preciceCpl )
    ! -------------------------------------------------------------------- !
    !> information about the mesh
    type(treelmesh_type), intent(in) :: tree
    !> contains the elements in spacetime-function
    type(tem_st_fun_listElem_type), intent(inout) :: stFun
    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    type(tem_precice_coupling_type), pointer :: preciceCpl
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: point(:,:)
    real(kind=rk), allocatable :: shiftpoints(:,:)
    real(kind=rk), allocatable :: serializedPnts(:)
    integer, allocatable :: idx(:)
    type(atl_varSys_data_type), pointer :: fPtr => NULL()
    type(ply_poly_project_type), pointer :: ply_poly => NULL()
    integer :: iLevel, iVar, iPnt, i, iCoord
    integer :: nPnts, varPos
    integer :: pos
    integer :: offsetX, offsetY, offsetZ
    ! -------------------------------------------------------------------- !
    write(logUnit(7),*) 'Write to preCICE Non-equidistant Points'
    call C_F_POINTER( stFun%solver_bundle, fPtr )

    do iLevel = tree%global%minLevel, tree%global%maxLevel
      pos = fPtr%solverData%poly_proj_posPtr(iLevel)
      !> Obtain the point coordinates to get values for
      ply_poly => fPtr%solverData%polyProjectPtr(pos)
      !> get the number of unquie points ( xyz should be same number)
      nPnts = stFun%pntData%pntLvl(iLevel)%grwPnt%coordX%nVals
      if ( nPnts == 0) then
!!VK        write(*,*) 'number of Points to write to precice is 0, cycle &
!!VK          & this spacetimefunction'
        cycle
      end if

      !> shift the points via the offset bit for the get_point routine
      allocate( shiftpoints(nPnts,3) )
      allocate( point(nPnts,3) )
      allocate( serializedPnts(3*nPnts) )
      do iPnt = 1, nPnts
        !> get the unquie list of points
        point(iPnt,1) = stFun%pntData%pntLvl(iLevel)%grwPnt%coordX%val(iPnt)
        point(iPnt,2) = stFun%pntData%pntLvl(iLevel)%grwPnt%coordY%val(iPnt)
        point(iPnt,3) = stFun%pntData%pntLvl(iLevel)%grwPnt%coordZ%val(iPnt)
        !> get the unquie list of points serialize them since precice needs
        !! a 1d array transform offet bit back
        offsetX = mod(ichar(stFun%pntData%pntLvl(iLevel)%offset_bit      &
          &                                             %val(iPnt)),4) - 1
        offsetY = mod(ichar(stFun%pntData%pntLvl(iLevel)%offset_bit         &
          &                                             %val(iPnt)),16)/4 - 1
        offsetZ = ichar(stFun%pntData%pntLvl(iLevel)%offset_bit      &
          &                                         %val(iPnt))/16 - 1
        !> shift the points according to the offset bit * small fraction
        !! of the boundingCubeLength hence it is mesh independent
        shiftPoints(iPnt,1) = point(iPnt,1) &
          & - offsetX*spacing(tree%global%BoundingCubeLength)
        shiftPoints(iPnt,2) = point(iPnt,2) &
          & - offsetY*spacing(tree%global%BoundingCubeLength)
        shiftPoints(iPnt,3) = point(iPnt,3) &
          & - offsetZ*spacing(tree%global%BoundingCubeLength)
      end do

      i=0
      do iPnt= 1, nPnts
        do iCoord = 1,3
          i=i+1
          serializedPnts(i)=point(iPnt,iCoord)
        end do
      end do

      allocate(idx(nPnts))
      !> loop over the variables to write to precice
      do iVar  = 1, preciceCpl%writeVar%nVars
        idx = 0
        !> get_point routine for unquie array of points of that variable
        varPos = preciceCpl%writeVar%varPos(iVar)
        !> get the setup_indices for the shifted points
        call varSys%method%val(varPos)%setup_indices( &
          & varSys = varSys,                          &
          & point  = shiftPoints,                     &
          & iLevel = iLevel,                          &
          & tree   = tree,                            &
          & nPnts  = nPnts,                           &
          & idx    = idx                              )

         !> Since we are sending 5 times the same idx, for the unique
         !! Points and all of them have the same idx, we just need to
         !! store it once
        if (iVar == 1) then
          call append(                                    &
            & me  = preciceCpl%pntIndex%indexLvl(iLevel), &
            & val = idx                                   )
        end if
      end do !iVar

      !> send the points to precice and get back the positionID
      !! allocate posID array
      allocate( preciceCpl%writeVar%posIDLvl(iLevel)%posIDs(nPnts) )
      write(logUnit(5),*) 'Send unique points to precice to get position' &
        & // ' IDs for non-equidistant Points'
      write(logUnit(5),*) 'nPnts ', nPnts
      call tem_precice_set_vertices_pos(                                    &
        & meshID     = preciceCpl%writeVar%meshID,                          &
        & nPoints    = nPnts,                                               &
        & points     = serializedPnts,                                      &
        & positionID = preciceCpl%writeVar%posIDLvl(iLevel)%posIDs(1:nPnts) )
      write(logUnit(1),*) 'Done set vertices with non-equidistant Points'

      !*****************************************************************
      !> From here on the Nearest-Projection is introduced, which has to
      !! be activated by the user in the ateles-config file.
      !*****************************************************************
      if (preciceCpl%interpolation%use_NP) then
        write(logUnit(5),*) 'Using Nearest-Projection, with non-equidistant' &
          & // 'Points'
        !> We assume here, that all points for a specific face are grouped
        !! together in order to make use of the posIDs
        !! Call to create the edges and triangles acoording to the
        !! dimension
        call atl_nearest_projection(                                        &
          & nPointsPerDir  = ply_poly%nQuadPointsPerDir,                    &
          & nPointsPerFace = ply_poly%nQuadPointsPerDir                     &
          &                   **(fPtr%solverData%equationPtr%nDimensions-1),&
          & nPnts          = nPnts,                                         &
          & preciceCpl     = preciceCpl,                                    &
          & scheme_dim     = fPtr%solverData%equationPtr%nDimensions,       &
          & iLevel         = iLevel                                         )
      end if !Nearest-Projection
      deallocate( point )
      deallocate( serializedPnts )
      deallocate( idx )
      deallocate( shiftpoints )
    end do !iLevel
  end subroutine atl_write_nonequiPoints
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_write_precice( stFunList, varSys, time, tree )
    ! -------------------------------------------------------------------- !
    !> This list contains all space-time-function
    type(tem_st_fun_linkedList_type), intent(inout), target :: stFunList
    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> The current absolute time.
    type(tem_time_type), intent(in) :: time
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: res(:)
    integer :: iLevel, iVar
    integer :: varPos
    integer :: iStFun, nPnts
    type(tem_precice_coupling_type), pointer :: preciceCpl => NULL()
    type(tem_st_fun_listElem_type), pointer :: current
    ! -------------------------------------------------------------------- !
    write(logUnit(7),*) 'Write to preCICE the values for the unique points '
    !> loop over linked list of spacetimefunction
    !> set the current to head
    current => stFunList%head
    do
      if( .not. associated(current) ) Exit

      !> loop over all shapes of that spacetimefunction
      do iStFun = 1, current%nVals
        !> check if there is a precice kind
        if ( any(current%val(:)%fun_kind == 'precice') ) then

          preciceCpl => current%val(iStFun)%precice_coupling
          !> check if this one has write flag
          if (preciceCpl%write) then
            do iLevel = tree%global%minLevel, tree%global%maxLevel
              !> get the number of unquie points ( xyz should be same number)
              nPnts = preciceCpl%pntIndex%indexLvl(iLevel)%nVals
              write(logUnit(7),*) 'nPnts = ', nPnts
              if ( nPnts == 0) then
                cycle
              end if
              allocate ( res(nPnts) )
              !> loop over the variables to write to precice
              do iVar  = 1, preciceCpl%writeVar%nVars
                !> Get_point routine for unquie array of points of that
                !! variable.
                !! Since the unique list of points is stored before, through the
                !! setup indices, we know depending on which kind of points
                !! we use, how many points we have. Thus we always access the
                !! right number of points, therefore we just have to call
                !! getvalofIndex just once independent from what kind of points
                !! we take for write to precice.
                varPos = preciceCpl%writeVar%varPos(iVar)
                call varSys%method%val(varPos)%get_valOfIndex(                 &
                  & varSys = varSys,                                           &
                  & nVals  = nPnts,                                            &
                  & iLevel = iLevel,                                           &
                  & time   = time,                                             &
                  & idx    = preciceCpl%pntIndex%indexLvl(iLevel)%val(:nPnts), &
                  & res    = res                                               )

                !> write the scalar value to precice
                call tem_precice_write(                              &
                  & dataID  = preciceCpl%writeVar%IDs(iVar),         &
                  & posID   = preciceCpl%writeVar%posIDLvl(iLevel)   &
                  &                              %posIDs(1:nPnts),   &
                  & npoints = nPnts,                                 &
                  & val     = res                                    )

              end do!iVar
              deallocate ( res )
            end do !iLevel
          end if !write flag
        end if ! kind=precice
      end do !iStFun
      current => current%next
    end do !linklist
  end subroutine atl_write_precice
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_write_precice_getPoint( stFunList, varSys, tree, mesh_list )
    ! -------------------------------------------------------------------- !
    !> This list contains all space-time-function
    type(tem_st_fun_linkedList_type), intent(inout), target :: stFunList
    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    type(atl_cube_elem_type), intent(in) ::               &
     & mesh_list(tree%global%minLevel:tree%global%maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iStFun
    type(tem_precice_coupling_type), pointer :: preciceCpl => NULL()
    type(tem_st_fun_listElem_type), pointer :: current
    ! -------------------------------------------------------------------- !
    write(logUnit(7),*) 'Write to preCICE the setup indices'
    !> loop over linked list of spacetimefunction
    !> set the current to head
    current => stFunList%head
    do
      if( .not. associated(current) ) Exit

      !>check if there is a precice kind
      if ( any(current%val(:)%fun_kind == 'precice') ) then
        !>loop over all shapes of that spacetimefunction
        do iStFun = 1, current%nVals

          preciceCpl => current%val(iStFun)%precice_coupling
          !> check if this one has write flag
          if (preciceCpl%write) then
            if (preciceCpl%interpolation%use_EQ_points) then
              call atl_write_equiPoints(          &
                & stFun      = current,           &
                & varSys     = varSys,            &
                & preciceCpl = preciceCpl,        &
                & mesh_list  = mesh_list,         &
                & tree       = tree               )
            else
              call atl_write_nonequiPoints(       &
                & stFun      = current,           &
                & varSys     = varSys,            &
                & preciceCpl = preciceCpl,        &
                & tree   = tree                   )
            end if !Equidistant Points or non-equidistant
          end if !write flag
        end do !iStFun
      end if ! kind=precice
      current => current%next
    end do !linklist
  end subroutine atl_write_precice_getPoint
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine atl_read_precice( stFunList, tree )
    ! -------------------------------------------------------------------- !
    !> This list contains all space-time-function
    type(tem_st_fun_linkedList_type), intent(inout), target :: stFunList
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: points(:,:)
    real(kind=rk), allocatable :: serializedPnts(:)
    integer :: iLevel, iCoord, iPnt, i
    integer :: iStFun, nPnts
    type(tem_precice_coupling_type), pointer :: preciceCpl => NULL()
    type(tem_st_fun_listElem_type), pointer :: current
    ! -------------------------------------------------------------------- !
    write(logUnit(7),*) 'Reading from preCICE'
    !> loop over linked list of spacetimefunction
    !> set the current to head
    current => stFunList%head
    do
      if( .not. associated(current) ) Exit

      !> loop over all shapes of that spacetimefunction
      do iStFun = 1, current%nVals
        !> check if there is a precice kind
        if ( any(current%val(:)%fun_kind == 'precice') ) then
          preciceCpl => current%val(iStFun)%precice_coupling
          do iLevel = tree%global%minLevel, tree%global%maxLevel
            !> get the number of unquie points ( xyz should be same number)
            nPnts = current%pntData%pntLvl(iLevel)      &
            &                                %grwPnt%coordX%nVals
            write(logUnit(7),*) 'nPnts for read = ', nPnts
            if ( nPnts == 0) then
              cycle
            end if
            allocate ( points(nPnts,3) )
            ! get the unquie list of points
            points(:,1) = current%pntData%pntLvl(iLevel)   &
              &                            %grwPnt%coordX%val(1:nPnts)
            points(:,2) = current%pntData%pntLvl(iLevel)   &
              &                            %grwPnt%coordY%val(1:nPnts)
            points(:,3) = current%pntData%pntLvl(iLevel)   &
              &                            %grwPnt%coordZ%val(1:nPnts)
            ! serialize them since precice needs a 1d array
            allocate( serializedPnts(3*nPnts) )
            i=0
            do iPnt= 1, nPnts
              do iCoord = 1,3
                i=i+1
                serializedPnts(i)=points(iPnt,iCoord)
              end do
            end do
            allocate( preciceCpl%readVar%posIDLvl(iLevel)%posIDs(nPnts) )
            call tem_precice_set_vertices_pos(                    &
              & meshID      = preciceCpl%readVar%meshID,          &
              & nPoints     = nPnts,                              &
              & points      = serializedPnts,                     &
              & positionID  = preciceCpl%readVar%posIDLvl(iLevel) &
              &                                 %posIDs(1:nPnts)  )
            deallocate( points )
            deallocate( serializedPnts )
          end do !iLevel
        end if ! kind=precice
      end do !iStFun
      current => current%next
    end do !linklist
  end subroutine atl_read_precice
  ! ************************************************************************ !

  end module atl_writePrecice_module
