! Copyright (c) 2011-2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011-2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014, 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

!> author: Jens Zudrop
!! This module provides central Ateles data type, containing the various data
!! of the simulation.
!!
!! The idea is to support different element types in the mesh, and collect
!! them in the atl_container_module.atl_element_container_type.
module atl_container_module
  use env_module,              only: rk, globalMaxLevels
  use treelmesh_module,        only: treelmesh_type, load_tem
  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit
  use tem_bc_prop_module,      only: tem_bc_prop_type, init_tem_bc_prop
  use tem_global_module,       only: tem_global_mesh_read
  use tem_time_module,         only: tem_time_type
  use tem_comm_module,         only: tem_commPattern_type
  use tem_comm_env_module,     only: tem_comm_env_type

  use flu_binding,     only : flu_state

  use atl_equation_module,                only: atl_equations_type
  use atl_facedata_module,                only: atl_init_facedata
  use atl_scheme_module,                  only: atl_init_scheme,        &
    &                                           atl_modg_scheme_prp,    &
    &                                           atl_modg_2d_scheme_prp, &
    &                                           atl_modg_1d_scheme_prp
  use atl_time_integration_module,        only: atl_global_timestep_type
  use atl_global_time_integration_module,       &
    &   only: atl_init_global_time_integration, &
    &         atl_global_time_integration_load
  use atl_modg_kernel_module,             only: atl_init_modg_kernel
  use atl_modg_2d_kernel_module,          only: atl_init_modg_2d_kernel
  use atl_modg_1d_kernel_module,          only: atl_init_modg_1d_kernel
  use atl_physCheck_module,               only: atl_physCheck_type, &
    &                                           atl_init_physCheck
  use atl_cube_container_module,          only: atl_init_cube_container, &
    &                                           atl_cube_container_type

  implicit none
  private

  !> container to collect all elements within our simulations domain
  type atl_element_container_type
    !> The cubes inside the domain
    type(atl_cube_container_type) :: cubes

    !> Global timediscretization type
    type(atl_global_timestep_type) :: time

    !> Physical checks information
    type(atl_physCheck_type) :: physCheck
  end type atl_element_container_type

  public :: atl_element_container_type, atl_init_elem_container


contains


  ! ****************************************************************************
  !> Initialize the container module.
  !!
  !! This routine builds up the necessary information as requested by the
  !! configuration in the Lua script <tt>conf</tt>.
  subroutine atl_init_elem_container( me, equation, conf, tree, time,          &
    &                                 readRestart, proc, commPattern, boundary )
    ! --------------------------------------------------------------------------

    !> complete domain
    type(atl_element_container_type), intent(inout) :: me

    !> Handle for the Lua config file
    type(flu_State) :: conf

    !> Description of the equation system
    type(atl_equations_type), intent(inout) :: equation

    !> Representation of the current time
    type(tem_time_type) :: time

    !> Should a restart be read?
    logical, intent(in) :: readRestart

    !> The mesh in treelm format.
    type(treelmesh_type), intent(inout) :: tree

    !> mpi communication environment with mpi communicator
    type(tem_comm_env_type) :: proc

    !> mpi communication pattern type
    type(tem_commPattern_type) :: commPattern

    type(tem_bc_prop_type), intent(inout) :: boundary
    ! --------------------------------------------------------------------------
    integer :: iLevel, dimen, nScalarFace, nScalarFlux
    logical :: use_levelWeights
    logical :: need_element_deviations
    real(kind=rk) :: modg_levelWeight(globalMaxLevels)
    ! --------------------------------------------------------------------------

    ! Depending on the equation system we may need to compute the
    ! maximal polynomial deviation in each element.
    need_element_deviations = equation%requiresDeviation

    write(logUnit(1),*) 'Setting up container parameters ...'

    if (.not. readRestart) then
      ! Not doing a restart, load the global mesh found in the configuration.
      call tem_global_mesh_read( me     = tree%global,   &
        &                        conf   = conf,          &
        &                        comm   = proc%comm,     &
        &                        myPart = proc%rank,     &
        &                        nParts = proc%comm_size )
    end if

    ! we init the scheme (same on all the levels)
    allocate(me%cubes%scheme_list(tree%global%minLevel:tree%global%maxLevel))
    write(logUnit(2),*) 'Inititalizing the scheme ...'
    call atl_init_scheme( me%cubes%scheme_list, &
      &                   conf, &
      &                   tree%global%minLevel, tree%global%maxLevel )

    ! init the general projection method into the list of projection methods
    allocate(me%cubes%poly_proj_pos(tree%global%minLevel:tree%global%maxLevel))

    use_levelWeights = .false.
    modg_levelWeight = 1.0_rk
    if ( (me%cubes%scheme_list(tree%global%minLevel)%scheme &
      &   == atl_modg_scheme_prp) .and. (proc%comm_size > 1) ) then
      use_levelWeights = (minval(me%cubes%scheme_list(:)%modg%maxPolyDegree) &
        &                 /= maxval(me%cubes%scheme_list(:)%modg%maxPolyDegree))
    end if

    if (use_levelWeights) then
      do iLevel=tree%global%minLevel,tree%global%maxLevel
        modg_levelWeight(iLevel) = real( me%cubes%scheme_list(iLevel)%nDofs, &
          &                              kind=rk )
      end do
    end if

    if (.not. readRestart) then
      if (associated(tree%global%Property)) deallocate(tree%global%property)
      if (associated(tree%Property)) deallocate(tree%property)
      ! Not doing a restart, load the mesh found in the configuration.
      if (use_levelWeights) then
        call load_tem( me          = tree,            &
          &            conf        = conf,            &
          &            comm        = proc%comm,       &
          &            myPart      = proc%rank,       &
          &            nParts      = proc%comm_size,  &
          &            levelWeight = modg_levelWeight )
      else
        call load_tem( me     = tree,          &
          &            conf   = conf,          &
          &            comm   = proc%comm,     &
          &            myPart = proc%rank,     &
          &            nParts = proc%comm_size )
      end if
    end if

    ! Load boundary conditions by treelm routines
    call init_tem_bc_prop( tree   = tree,      &
      &                    mypart = proc%rank, &
      &                    comm   = proc%comm, &
      &                    bc     = boundary   )

    call atl_global_time_integration_load(           &
      &    me             = me%time,                 &
      &    conf           = conf,                    &
      &    equation       = equation                 )

    if (me%time%control%use_modal_estimate) then
      ! If the timestep control should be done with a modal estimation,
      ! we need to compute the maximal polynomial deviation in each
      ! element.
      need_element_deviations = .true.
    end if

    ! Right now, we only support cubical elements, initialize them
    ! now, with the mesh read so far.
    call atl_init_cube_container(                            &
      &    tree                    = tree,                   &
      &    boundary                = boundary,               &
      &    cube_container          = me%cubes,               &
      &    conf                    = conf,                   &
      &    equation                = equation,               &
      &    proc                    = proc,                   &
      &    commPattern             = commPattern,            &
      &    need_element_deviations = need_element_deviations )

    ! Initialize the kernels states (including the scheme)
    call init_kernel( time           = time,       &
      &               cube_container = me%cubes,   &
      &               equation       = equation,   &
      &               tree           = tree,       &
      &               commPattern    = commPattern )

    ! Initialize the face states
    select case( me%cubes%scheme_list(tree%global%minLevel)%scheme )
    case(atl_modg_scheme_prp)
      dimen = 3
    case(atl_modg_2d_scheme_prp)
      dimen = 2
    case(atl_modg_1d_scheme_prp)
      dimen = 1
    case default
      write(logUnit(1),*) 'ERROR in atl_init_elem_container: '
      write(logUnit(1),*) 'Not able to determine number of scalars on the' &
        & // ' face, stopping ... '
      stop
    end select
    ! For the state on the face, we need the state itself plus
    ! all the derivatives of all state variables in all spatial directions
    nScalarFace = equation%varSys%nScalars * (1+equation%nDerivatives*dimen)
    ! For the flux, we need the flux plus the stabilization flux
    nScalarFlux = equation%varSys%nScalars * (1+equation%nDerivatives)
    call atl_init_facedata( faces         = me%cubes%mesh_list(:)%faces,       &
      &                     facedata      = me%cubes%facedata_list,            &
      &                     minLevel      = tree%global%minLevel,              &
      &                     maxLevel      = tree%global%maxLevel,              &
      &                     nDim          = dimen,                             &
      &                     nScalarsState = nScalarFace,                       &
      &                     nScalarsFlux  = nScalarFlux,                       &
      &                     nFaceDofs     = me%cubes%scheme_list(:)%nFaceDofs, &
      &                     boundary      = me%cubes%boundary_list             )

    ! Initialize the global timestepping scheme
    call atl_init_global_time_integration(           &
      &    me             = me%time,                 &
      &    minLevel       = tree%global%minLevel,    &
      &    maxLevel       = tree%global%maxLevel,    &
      &    statedata_list = me%cubes%statedata_list, &
      &    mesh_list      = me%cubes%mesh_list,      &
      &    equation       = equation                 )

    ! Init the physical checks
    call atl_init_physCheck( me%physCheck, conf, equation )

  end subroutine atl_init_elem_container
  ! ****************************************************************************


  ! ****************************************************************************
  !> Initialize the kernel states for all parts of the mesh.
  !!
  !! Intializes the kernel for the cubic elements of the mesh.
  subroutine init_kernel(time, cube_container, equation, tree, commPattern)
    ! -------------------------------------------------------------------------!
    !> current time
    type(tem_time_type), intent(in)             :: time
    !> container of cubic elements.
    type(atl_cube_container_type), intent(inout) :: cube_container
    !> The equation we try to solve.
    type(atl_equations_type), intent(in) :: equation
    !> The tree representation of the mesh.
    type(treelmesh_type)                        :: tree
    !> mpi communication pattern type
    type(tem_commPattern_type)                  :: commPattern
    ! -------------------------------------------------------------------------!
    integer :: iList
    integer, allocatable :: nTotal(:)
    ! -------------------------------------------------------------------------!

    write(logUnit(1),*) '...init kernel states...'

    ! get the number of total cells per level
    allocate(nTotal(tree%global%minLevel:tree%global%maxLevel))
    do iList = tree%global%minLevel, tree%global%maxLevel
      nTotal(iList) = size(cube_container%mesh_list(iList)%descriptor%total, 1)
    end do

    ! now we continue with the kernel itself, since the informations
    ! depend on the data stored in scheme. This routine also inits the
    ! kerneldata fields.
    select case(cube_container%scheme_list(tree%global%minLevel)%scheme)
    case(atl_modg_scheme_prp)
      ! modg scheme
      call atl_init_modg_kernel(                                       &
        & kerneldata_list     = cube_container%kerneldata_list,        &
        & statedata_list      = cube_container%statedata_list,         &
        & statedata_stab_list = cube_container%statedata_stab_list,    &
        & mesh_list           = cube_container%mesh_list,              &
        & scheme_list         = cube_container%scheme_list,            &
        & boundary_list       = cube_container%boundary_list,          &
        & boundary_stab_list  = cube_container%boundary_stab_list,     &
        & equation            = equation,                              &
        & tree                = tree,                                  &
        & time                = time,                                  &
        & commPattern         = commPattern,                           &
        & need_deviation      = cube_container%need_element_deviations )

    case(atl_modg_2d_scheme_prp)
      ! modg scheme
      call atl_init_modg_2d_kernel(                                    &
        & kerneldata_list     = cube_container%kerneldata_list,        &
        & statedata_list      = cube_container%statedata_list,         &
        & statedata_stab_list = cube_container%statedata_stab_list,    &
        & mesh_list           = cube_container%mesh_list,              &
        & scheme_list         = cube_container%scheme_list,            &
        & boundary_list       = cube_container%boundary_list,          &
        & boundary_stab_list  = cube_container%boundary_stab_list,     &
        & equation            = equation,                              &
        & tree                = tree,                                  &
        & time                = time,                                  &
        & commPattern         = commPattern,                           &
        & need_deviation      = cube_container%need_element_deviations )

    case(atl_modg_1d_scheme_prp)
      ! modg scheme
      call atl_init_modg_1d_kernel(                                    &
        & kerneldata_list     = cube_container%kerneldata_list,        &
        & statedata_list      = cube_container%statedata_list,         &
        & statedata_stab_list = cube_container%statedata_stab_list,    &
        & mesh_list           = cube_container%mesh_list,              &
        & scheme_list         = cube_container%scheme_list,            &
        & boundary_list       = cube_container%boundary_list,          &
        & boundary_stab_list  = cube_container%boundary_stab_list,     &
        & equation            = equation,                              &
        & tree                = tree,                                  &
        & time                = time,                                  &
        & commPattern         = commPattern,                           &
        & need_deviation      = cube_container%need_element_deviations )

    case default
      write(logUnit(1),*) 'ERROR in init_kernel: ' // &
        & ' initialization of this kernel is not supported (yet?), stopping....'
      call tem_abort()

    end select

  end subroutine init_kernel
  ! ****************************************************************************

end module atl_container_module
