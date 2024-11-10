! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2016,2018,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

!> Handling of restart data.
!! Restart data is used to provide the state fields at given points in time.
!! They can be used to resume the simulation again in a further simulation run
!! or for post-processing to create visualizations.
!!
!! The restart information is split into two parts, a header file describing
!! the contained data and a binary file holding the degrees of freedom from the
!! polynomial representation of the state fields for each element in the mesh.
!! Please note, that the maximal polynomial degree of all elements will be used
!! for each element when writing the field state to disk to allow for efficient
!! uniform file access.
!! It is also possible to configure a different polynomial degree in the
!! simulation run than what is available in the restart data. This allows for
!! higher resolved simulations to follow on coarser simulations.
!! The degrees of freedoms are sorted correctly into the different representations
!! by the routines in the [[ply_transfer_module]].
!!
!! When to write restarts is controlled by a time_control object. Details on the
!! configuration of restart files is found in the [[tem_restart_module]].
module atl_restart_module
  use env_module,                only: rk, io_buffer_size, newUnit

  use aotus_module,              only: flu_state, aot_get_val
  use aot_table_module,          only: aot_table_open, &
    &                                  aot_table_close
  use aot_out_module,            only: aot_out_type,        &
    &                                  aot_out_open,        &
    &                                  aot_out_open_table,  &
    &                                  aot_out_close_table, &
    &                                  aot_out_val,         &
    &                                  aot_out_close

  use tem_restart_module,        only: tem_restart_type,       &
    &                                  tem_restart_writeData,  &
    &                                  tem_restart_openRead,   &
    &                                  tem_restart_closeRead,  &
    &                                  tem_restart_closeWrite, &
    &                                  tem_restart_openWrite,  &
    &                                  tem_restart_readData,   &
    &                                  tem_init_restart
  use treelmesh_module,          only: treelmesh_type
  use tem_topology_module,       only: tem_levelof
  use tem_aux_module,            only: tem_abort, tem_open_distconf
  use tem_tools_module,          only: upper_to_lower, &
    &                                  tem_getOptValOrDef
  use tem_logging_module,        only: logUnit
  use tem_comm_env_module,       only: tem_comm_env_type
  use tem_time_module,           only: tem_time_type
  use tem_timer_module,          only: tem_startTimer, tem_stopTimer
  use tem_timeControl_module,    only: tem_timeControl_check
  use tem_simControl_module,     only: tem_simControl_type, &
    &                                  tem_simControl_dump_now
  use tem_solveHead_module,      only: tem_solveHead_type
  use tem_status_module,         only: tem_stat_interval
  use tem_varMap_module,         only: tem_varMap_type, &
    &                                  tem_create_varMap
  use tem_varSys_module,         only: tem_varSys_type

  use atl_cube_container_module, only: atl_cube_container_type
  use atl_equation_module,       only: atl_equations_type
  use atl_scheme_module,         only: atl_modg_scheme_prp,    &
    &                                  atl_modg_2d_scheme_prp, &
    &                                  atl_modg_1d_scheme_prp
  use atl_timer_module,          only: atl_timerHandles

  use ply_dof_module,            only: P_space, Q_space
  use ply_transfer_module,       only: ply_transfer_dofs

  implicit none

  private

  public :: atl_initRestart
  public :: atl_writeRestart
  public :: atl_writeRestartIfNecessary
  public :: atl_readRestart
  public :: atl_writeSolverSpecInfo


contains


  ! ------------------------------------------------------------------------ !
  !> Initializes writing restart files, if activated.
  !!
  !! Prepares the solver specific information for the restart header and calls
  !! initialization of underlying restart functionality in treelm.
  subroutine atl_initRestart( restart, equation, solver, mesh, tree, nDofsIO, &
    &                         scheme_dim )
    ! -------------------------------------------------------------------- !
    !> Description of the restart structure to initialize.
    type(tem_restart_type), intent(inout) :: restart
    !> The variable system of the current equation system.
    type(atl_equations_type), intent(in) :: equation
    !> General description of the deployed solver.
    type(tem_solveHead_type),  intent(in) :: solver
    !> Levelwise representation of the mesh.
    type(atl_cube_container_type), intent(in) :: mesh
    !> The octree of the mesh.
    type(treelmesh_type),      intent(in) :: tree
    !> The number of degrees of freedoms to write to/read from a restart file
    integer, intent(in) :: nDofsIO
    !> Dimensionality to be used for the state fields.
    integer, intent(in) :: scheme_dim
    ! -------------------------------------------------------------------- !
    integer :: solSpec_unit
    type(tem_varMap_type) :: varMap
    ! -------------------------------------------------------------------- !

    if (restart%controller%writeRestart    &
      & .or. restart%controller%readRestart) then
      write(logUnit(1),*) 'initializing restart'
      solSpec_unit = -1
      ! Create stateVarmap to dump the state variables in restart
      call tem_create_varMap(                                                  &
        &    varName = equation%varsys                                         &
        &                      %varName                                        &
        &                      %val(equation%stateVar(1:equation%varSys        &
        &                                                       %nStateVars)), &
        &    varSys  = equation%varSys,                                        &
        &    varMap  = varMap                                                  )

      if (tree%global%myPart == 0 ) then
        solSpec_Unit = newUnit()
        open(unit=solSpec_Unit, status='scratch')
        call atl_writeSolverSpecInfo( mesh       = mesh,                 &
          &                           scheme_dim = scheme_dim,           &
          &                           minlevel   = tree%global%minlevel, &
          &                           outUnit    = solSpec_unit          )
      end if

      call tem_init_restart( me           = restart,     &
        &                    solver       = solver,      &
        &                    varMap       = varMap,      &
        &                    tree         = tree,        &
        &                    nDofs_write  = nDofsIO,     &
        &                    solSpec_unit = solSpec_unit )

      if (restart%comm%rank == 0 ) then
        close(solSpec_Unit)
      end if

    end if

  end subroutine atl_initRestart
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Writes a restart file, if necessary.
  !!
  !! This routine checks whether writing restart files is activated at all. If
  !! so, it further checks whether a restart file has to be written at current
  !! point in time. If this also applies, a third check ensures that there
  !! wasn't already a restart file written for the current point in time.
  subroutine atl_writeRestartIfNecessary( restart, simControl, equation, tree, &
      &                                   mesh, force                          )
    ! -------------------------------------------------------------------- !
    !> Informations and states of the cubic mesh.
    type(atl_cube_container_type),intent(inout) :: mesh
    type(tem_simControl_type), intent(inout) :: simControl
    !> the restart infotmation
    type(tem_restart_type) :: restart
    !> mesh, provided in treelm format
    type(treelmesh_type) :: tree
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> Set this to true to override the restart interval check.
    logical, optional :: force
    ! -------------------------------------------------------------------- !
    logical :: doRestart, lForce
    ! -------------------------------------------------------------------- !

    lForce = tem_getOptValOrDef(force, .false.)

    if (restart%controller%writeRestart) then
      ! Check if we have to write a restart file
      ! (only if writeRestart is active)
      call tem_timeControl_check(                        &
        &    me        = restart%controller%timeControl, &
        &    now       = simControl%now,                 &
        &    comm      = restart%comm%comm,              &
        &    triggered = doRestart                       )

      if (doRestart .or. lForce) then
        if ( .not. simControl%status%bits(tem_stat_interval) ) then
          ! If this is not an interval (where the current time is anyway
          ! printed), print the current time.
          call tem_simControl_dump_now(simControl, logUnit(1))
        end if
        ! Check whether there was already a restart file written for the current
        ! timestep. If so, return immediately.
        if ( restart%lastWritten%iter == simControl%now%iter ) then
          write(logUnit(1),*) 'Restart file already written for current' &
            &                 // ' point in time!'
        else
          write(logUnit(1),*) 'Writing restart in current iteration'
          call atl_writeRestart(                         &
            &    mesh        = mesh,                     &
            &    restart     = restart,                  &
            &    tree        = tree,                     &
            &    equation    = equation,                 &
            &    timing      = simControl%now,           &
            &    varSys      = equation%varSys,          &
            &    timerHandle = atl_timerHandles%wRestart )
        end if
      end if
    end if

  end subroutine atl_writeRestartIfNecessary
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Writes a restart file for the current point in time.
  subroutine atl_writeRestart( mesh, restart, tree, equation, timing, &
    &                          varSys, timerHandle                    )
    ! -------------------------------------------------------------------- !
    !> Informations and states of the cubic mesh.
    type(atl_cube_container_type),intent(inout) :: mesh
    !> the restart infotmation
    type(tem_restart_type) :: restart
    !> mesh, provided in treelm format
    type(treelmesh_type) :: tree
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    !> current simulation time information
    type(tem_time_type), intent(in) :: timing
    type(tem_varSys_type), intent(in) :: varSys
    !> Timer handle to measure write restart time
    integer, intent(in) :: timerHandle
    ! -------------------------------------------------------------------- !
    ! local variables
    real(kind=rk), allocatable :: chunk(:)
    integer :: iChunk
    integer :: elemOff
    integer :: nDims
    integer :: polyspace
    integer :: maxpolydegree
    ! -------------------------------------------------------------------- !
    call tem_startTimer(timerHandle = timerHandle )

    allocate( chunk(io_buffer_size) )

    !>@todo HK: Not so sure about this output, is probably also already done
    !!          by the treelm routines.
    write(logUnit(1),'(A23,E14.6,A15,I10)') 'Writing restart at t = ', &
      &                                     timing%sim, ', iterations = ', &
      &                                     timing%iter

    ! open the file, write header and prepare buffering
    call tem_restart_openWrite( me     = restart, &
      &                         tree   = tree,    &
      &                         timing = timing,  &
      &                         varsys = varsys   )


    select case(mesh%scheme_list(tree%global%minlevel)%scheme)
    case(atl_modg_scheme_prp)
      ndims = 3
      polyspace = mesh%scheme_list(tree%global%minlevel)%modg%basisType
      maxPolyDegree = maxval(mesh%scheme_list(:)%modg%maxPolyDegree)
    case(atl_modg_2d_scheme_prp)
      ndims = 2
      polyspace = mesh%scheme_list(tree%global%minlevel)%modg_2d%basisType
      maxPolyDegree = maxval(mesh%scheme_list(:)%modg_2d%maxPolyDegree)
    case(atl_modg_1d_scheme_prp)
      ndims = 1
      polyspace = mesh%scheme_list(tree%global%minlevel)%modg_1d%basisType
      maxPolyDegree = maxval(mesh%scheme_list(:)%modg_1d%maxPolyDegree)
    end select

    ! loop over chunks
    elemOff = 0
    do iChunk = 1, restart%write_file%nChunks

      ! set the number of elements that are on the stack
      restart%nChunkElems = min( restart%write_file%chunkSize, &
        &                        tree%nElems-elemOff           )

      call serializeData( mesh        = mesh,          &
        &                 restart     = restart,       &
        &                 tree        = tree,          &
        &                 equation    = equation,      &
        &                 chunk       = chunk,         &
        &                 iChunk      = iChunk,        &
        &                 chunkdegree = maxpolydegree, &
        &                 chunkspace  = polyspace,     &
        &                 statespace  = polyspace,     &
        &                 ndims       = ndims,         &
        &                 reverse     = .false.        )

      ! write the data to file
      call tem_restart_writeData( restart, chunk )

      elemOff = elemOff + restart%nChunkElems

    end do

    ! close the file and write last header
    call tem_restart_closeWrite( me     = restart, &
      &                          tree   = tree,    &
      &                          timing = timing,  &
      &                          varsys = varsys   )

    deallocate( chunk )

    write(logUnit(1),*) '... finished writing restart data'

    call tem_stopTimer( timerHandle = timerHandle )

  end subroutine atl_writeRestart
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Read the serialized restart file into the state vectors
  !!
  subroutine atl_readRestart( mesh, restart, tree, equation, proc )
    ! -------------------------------------------------------------------- !
    !> Informations and states of the cubic mesh.
    type(atl_cube_container_type),intent(inout) :: mesh
    type(tem_restart_type) :: restart !< the restart infotmation
    !> mesh, provided in treelm format
    type(treelmesh_type) :: tree
    type(atl_Equations_type), intent(in) :: equation
    !> The schemes on the levels
    type(tem_comm_env_type), intent(in) :: proc
    ! -------------------------------------------------------------------- !
    ! local variables
    type(flu_State) :: conf
    real(kind=rk), allocatable :: chunk(:)
    integer :: maxPolyDegree
    integer :: iChunk
    integer :: elemOff
    integer :: ndims
    integer :: polyspace
    integer :: read_space
    integer :: read_dims
    integer :: read_degree
    integer :: ply_handle
    character :: space_char
    character :: defspace
    integer :: iError
    ! -------------------------------------------------------------------- !

    allocate( chunk( io_buffer_size ))

    write(logUnit(1),*) 'Reading restart...'

    ! open the file, read header and prepare buffering
    call tem_restart_openRead( me = restart )

    select case(mesh%scheme_list(tree%global%minlevel)%scheme)
    case(atl_modg_scheme_prp)
      ndims = 3
      polyspace = mesh%scheme_list(tree%global%minlevel)%modg%basisType
      maxPolyDegree = maxval(mesh%scheme_list(:)%modg%maxPolyDegree)
    case(atl_modg_2d_scheme_prp)
      ndims = 2
      polyspace = mesh%scheme_list(tree%global%minlevel)%modg_2d%basisType
      maxPolyDegree = maxval(mesh%scheme_list(:)%modg_2d%maxPolyDegree)
    case(atl_modg_1d_scheme_prp)
      ndims = 1
      polyspace = mesh%scheme_list(tree%global%minlevel)%modg_1d%basisType
      maxPolyDegree = maxval(mesh%scheme_list(:)%modg_1d%maxPolyDegree)
    end select

    if (polyspace == P_Space) then
      defspace = 'P'
    else
      defspace = 'Q'
    end if

    call tem_open_distconf( L        = conf,                                  &
      &                     filename = trim(restart%controller%readFileName), &
      &                     proc     = proc                                   )
    call aot_table_open( L       = conf,        &
      &                  thandle = ply_handle,  &
      &                  key     = 'polynomial' )
    call aot_get_val( L       = conf,             &
      &               thandle = ply_handle,       &
      &               key     = 'dimensionality', &
      &               val     = read_dims,        &
      &               ErrCode = iError,           &
      &               default = ndims             )
    call aot_get_val( L       = conf,       &
      &               thandle = ply_handle, &
      &               key     = 'space',    &
      &               val     = space_char, &
      &               ErrCode = iError,     &
      &               default = defspace    )
    call aot_get_val( L       = conf,         &
      &               thandle = ply_handle,   &
      &               key     = 'degree',     &
      &               val     = read_degree,  &
      &               ErrCode = iError,       &
      &               default = maxPolyDegree )
    call aot_table_close( L       = conf,      &
      &                   thandle = ply_handle )

    if (upper_to_lower(space_char) == 'p') then
      read_space = P_Space
    else
      read_space = Q_Space
    end if

    ! loop over chunks
    elemOff = 0
    do iChunk = 1, restart%read_file%nChunks
      ! set the number of elements that are on the stack
      restart%nChunkElems = min( restart%read_file%chunkSize, &
        &                        tree%nElems-elemOff          )

      call tem_restart_readData( restart, chunk )

      call serializeData( mesh        = mesh,        &
        &                 restart     = restart,     &
        &                 tree        = tree,        &
        &                 equation    = equation,    &
        &                 chunk       = chunk,       &
        &                 iChunk      = iChunk,      &
        &                 chunkdegree = read_degree, &
        &                 chunkspace  = read_space,  &
        &                 statespace  = polyspace,   &
        &                 ndims       = ndims,       &
        &                 reverse     = .true.       )

      elemOff = elemOff + restart%nChunkElems

    end do

    call tem_restart_closeRead( me = restart )

    deallocate( chunk )

  end subroutine atl_readRestart
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This subroutine serializes the given data to perform a restart.
  !!
  subroutine serializeData( mesh, restart, tree, equation, chunk, iChunk, &
    &                       chunkdegree, chunkspace, statespace, ndims,   &
    &                       reverse                                       )
    ! -------------------------------------------------------------------- !
    !> Informations and states of the cubic mesh.
    type(atl_cube_container_type),intent(inout) :: mesh
    type(tem_restart_type), intent(inout) :: restart
    type(treelmesh_type), intent(in) :: tree
    !> The equation system you use.
    type(atl_equations_type), intent(in) :: equation
    real(kind=rk), intent(inout) :: chunk(:)
    integer, intent(in) :: iChunk !< the current chunk

    !> Polynomial degree of the data in the IO chunk
    integer, intent(in) :: chunkdegree

    !> Polynomial space of the serialized data
    integer, intent(in) :: chunkspace

    !> Polynomial space of the state data
    integer, intent(in) :: statespace

    !> Dimensionality of the polynomials
    integer, intent(in) :: nDims

    !> deserialize the data instead of serializing it?
    logical, intent(in) :: reverse
    ! -------------------------------------------------------------------- !
    integer :: elemOff
    integer :: elemPos
    integer :: nDofs
    integer :: level_deg
    integer :: iDoF
    integer :: iLevel
    integer :: iElem
    integer :: iIndex
    integer :: iScalar
    integer :: nScalars
    integer :: dof_off
    real(kind=rk), allocatable :: modal_background(:)
    real(kind=rk), allocatable :: singleVar(:)
    ! -------------------------------------------------------------------- !

    nScalars = equation%varSys%nScalars

    ! Global counter for all words put into the buffer.
    iIndex = 0

    if (reverse) then
      ! Number of elements written by previous chunks.
      elemOff = ((iChunk-1)*restart%read_file%chunkSize)
      allocate(singleVar(restart%read_file%nDofs))
      ndofs = restart%read_file%ndofs

    else
      ! Number of elements written by previous chunks.
      elemOff = ((iChunk-1)*restart%write_file%chunkSize)

      ! We fill the chunk with zeros in advanced to make sure: if we have
      ! less number of polynomial degrees of freedom on the current level,
      ! then we append zeros for the missing higher order modes
      chunk = 0.0_rk
      allocate(singleVar(restart%write_file%nDofs))
      ndofs = restart%write_file%ndofs
    end if

    !!@todo The modal background should not be applied here, but rather later
    !!      on.
    allocate(modal_background(nScalars))

    select case(equation%eq_kind)
    case('acoustic_2d')
      modal_background(1) = equation%acoustic%density_0
      modal_background(2) = equation%acoustic%velocity_0(1)
      modal_background(3) = equation%acoustic%velocity_0(2)
    case('lineareuler_2d')
      modal_background(1) = equation%lineareuler%density_0
      modal_background(2) = equation%lineareuler%velocity_0(1)
      modal_background(3) = equation%lineareuler%velocity_0(2)
      modal_background(4) = equation%lineareuler%pressure_0
    case('acoustic')
      modal_background(1) = equation%acoustic%density_0
      modal_background(2) = equation%acoustic%velocity_0(1)
      modal_background(3) = equation%acoustic%velocity_0(2)
      modal_background(4) = equation%acoustic%velocity_0(3)
    case('lineareuler')
      modal_background(1) = equation%lineareuler%density_0
      modal_background(2) = equation%lineareuler%velocity_0(1)
      modal_background(3) = equation%lineareuler%velocity_0(2)
      modal_background(4) = equation%lineareuler%velocity_0(3)
      modal_background(5) = equation%lineareuler%pressure_0
    case default
      modal_background = 0.0_rk
    end select

    ! For the current chunk (given from the caller), treat the corresponding
    ! elements from the tree ID list
    ! Run over all variable systems (should be only one here!)
    elems: do iElem =  elemOff+1, elemOff+restart%nChunkElems
      ! On which level lies the current element
      iLevel = tem_LevelOf( tree%treeID( iElem ))
      elemPos = mesh%levelPointer(iElem)
      select case(mesh%scheme_list(tree%global%minlevel)%scheme)
      case(atl_modg_1d_scheme_prp)
        level_deg = mesh%scheme_list(iLevel)%modg_1d%maxPolyDegree
      case(atl_modg_2d_scheme_prp)
        level_deg = mesh%scheme_list(iLevel)%modg_2d%maxPolyDegree
      case(atl_modg_scheme_prp)
        level_deg = mesh%scheme_list(iLevel)%modg%maxPolyDegree
      case default
        write(logUnit(1),*) 'ERROR in serializeData: ' &
          &                 // 'unknown spatial scheme, stopping ...'
        call tem_abort()
      end select

      if (reverse) then
        ! De-Serialize
        do iScalar = 1,nScalars
          dof_off = iScalar + (iElem-(elemOff+1)) * nScalars * nDofs
          do iDoF=1,nDofs
            iIndex = dof_off + (iDof-1) * nScalars
            singleVar(iDoF) = chunk(iIndex)
          end do
          singleVar(1) = singleVar(1) - modal_background(iScalar)
          call ply_transfer_dofs( indat     = singleVar,                       &
            &                     inspace   = chunkspace,                      &
            &                     indegree  = chunkDegree,                     &
            &                     outdat    = mesh%statedata_list(iLevel)      &
            &                                     %state(elemPos, :, iScalar), &
            &                     outspace  = statespace,                      &
            &                     outdegree = level_deg,                       &
            &                     ndims     = ndims                            )
        end do
      else
        ! Serialize
        do iScalar = 1,nScalars
          call ply_transfer_dofs( outdat    = singleVar,                       &
            &                     outspace  = chunkspace,                      &
            &                     outdegree = chunkDegree,                     &
            &                     indat     = mesh%statedata_list(iLevel)      &
            &                                     %state(elemPos, :, iScalar), &
            &                     inspace   = statespace,                      &
            &                     indegree  = level_deg,                       &
            &                     ndims     = ndims                            )
          singleVar(1) = singleVar(1) + modal_background(iScalar)
          dof_off = iScalar + (iElem-(elemOff+1)) * nScalars * nDofs
          do iDoF=1,nDofs
            iIndex = dof_off + (iDof-1) * nScalars
            chunk(iIndex) = singleVar(iDof)
          end do
        end do
      end if

    end do elems

    deallocate(modal_background)
    deallocate(singleVar)

  end subroutine serializeData
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Write solver specific info to scratch file
  subroutine atl_writeSolverSpecInfo( mesh, scheme_dim, minlevel, outUnit )
    ! -------------------------------------------------------------------- !
    !> The mesh type
    type(atl_cube_container_type), intent(in) :: mesh
    !> The scheme dimension
    integer, intent(in) :: scheme_dim
    !> Minimal level in the domain to pick the polynomial basis type.
    integer, intent(in) :: minlevel
    !> unit to output solver info in lua format
    integer, intent(inout) :: outUnit
    !> The lua chunk to create the solver specific header parts with.
    type(aot_out_type) :: out_conf
    ! -------------------------------------------------------------------- !
    integer :: maxPolyDegree
    integer :: polyspace
    character :: polychar
    ! -------------------------------------------------------------------- !
    ! Write scratch file only from root process and
    ! outUnit is not initialized before
    call aot_out_open( put_conf = out_conf, outUnit = outUnit )

    polyspace = 0
    select case(scheme_dim)
    case(3)
      polyspace = mesh%scheme_list(minlevel)%modg%basisType
      maxPolyDegree = maxval( mesh%scheme_list(:)%modg%maxPolyDegree )
    case(2)
      polyspace = mesh%scheme_list(minlevel)%modg_2d%basisType
      maxPolyDegree = maxval( mesh%scheme_list(:)%modg_2d%maxPolyDegree )
    case(1)
      polyspace = mesh%scheme_list(minlevel)%modg_1d%basisType
      maxPolyDegree = maxval( mesh%scheme_list(:)%modg_1d%maxPolyDegree )
    case default
      call tem_abort('Wrong scheme dimension! Select 1, 2 or 3.')
    end select

    if (polyspace == P_Space) then
      polychar = 'P'
    else
      polychar = 'Q'
    end if

    call aot_out_open_table( put_conf = out_conf,    &
      &                      tname    = 'polynomial' )

    call aot_out_val( put_conf = out_conf,         &
      &               vname    = 'dimensionality', &
      &               val      = scheme_dim        )

    call aot_out_val( put_conf = out_conf, &
      &               vname    = 'space',  &
      &               val      = polychar  )

    call aot_out_val( put_conf = out_conf,     &
      &               vname    = 'degree',     &
      &               val      = maxpolydegree )

    call aot_out_close_table( put_conf = out_conf )

    call aot_out_close( put_conf = out_conf )

  end subroutine atl_writeSolverSpecInfo
  ! ------------------------------------------------------------------------ !

end module atl_restart_module
