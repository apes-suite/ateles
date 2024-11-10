! Copyright (c) 2011-2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
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

!> Module to provide general solver parameters.
!!
!! This module provides a solver module variable, to collect general data on
!! the simulation and the solver in a central place.
module atl_solver_param_module
  use env_module,             only: PathLen
  use tem_general_module,     only: tem_general_type, tem_load_general
  use tem_restart_module,     only: tem_load_restart
  use ply_sampled_tracking_module, only: ply_sampled_tracking_type, &
    &                                    ply_sampled_tracking_load
  use treelmesh_module,       only: treelmesh_type
  use tem_bc_prop_module,     only: tem_bc_prop_type
  use tem_timeControl_module, only: tem_timeControl_start_at_sim
  use tem_logging_module,     only: tem_logging_load_primary
  use tem_debug_module,       only: tem_debug_load_main

  implicit none

  private

  type atl_solver_param_type
    !> general data coming from treelem
    type(tem_general_type) :: general

    !> Tracking objects to capture subsets of the overall simulation,
    !! or derived quantities.
    type(ply_sampled_tracking_type) :: plySampleTrack

    !> Boundary properties of elements in the mesh.
    type(tem_bc_prop_type) :: boundary

    !> Polynomial degree for each variable in the variable system.
    !!
    !! Needed for subsampled tracking, currently all have the same degree.
    integer, allocatable :: var_degree(:)

    !> Polynomial degree for each level in the mesh.
    integer, allocatable :: lvl_degree(:)

    !> Polynomial space for each variable in the variable system.
    !!
    !! Needed for subsampled tracking, currently all have the same space.
    integer, allocatable :: var_space(:)

    !> Type for balanciing weights
    character(len=PathLen) :: weights_file
  end type

  public :: atl_solver_param_type, atl_load_solver_parameters

contains

  !> Routine to initialize the global parameters, sets the solver
  !! module variable.
  !!
  !! Subroutine to initialize the global parameters like simulation name, etc.
  !! as specified in the given lua configuration file. The configuration file
  !! is passed as a handle to this subroutine.
  subroutine atl_load_solver_parameters(me, tree)
    ! --------------------------------------------------------------------------!
    type(atl_solver_param_type), intent(inout) :: me
    type(treelmesh_type), intent(inout) :: tree
    ! --------------------------------------------------------------------------!

    ! load and initialize logUnit
    call tem_logging_load_primary(conf = me%general%solver%conf(1), &
      &                           rank = me%general%proc%rank       )

    ! load and initialize debug unit
    call tem_debug_load_main(conf = me%general%solver%conf(1), &
      &                      rank = me%general%proc%rank       )

    ! Loading general settings.
    !>@todo HK: this is strange, the array of conf, should probably be moved
    !!          somewhere else.
    call tem_load_general( me   = me%general,               &
      &                    conf = me%general%solver%conf(1) )

    ! Init restart, this will fill also the tree, if restart is read!
    ! Now will be updated to refer to the time provided in the restart file.
    call tem_load_restart( me       = me%general%restart,        &
      &                    conf     = me%general%solver%conf(1), &
      &                    tree     = tree,                      &
      &                    timing   = me%general%simControl%now, &
      &                    globProc = me%general%proc            )

    if (me%general%restart%controller%readRestart) then
      ! If we are to restart the simulation, let the simulation control start
      ! at that the given time.
      call tem_timeControl_start_at_sim(              &
        &    me  = me%general%simControl%timeControl, &
        &    now = me%general%simControl%now          )
    end if

    ! Loading the tracker for the variable at a point or line
    call ply_sampled_tracking_load( me   = me%plySampleTrack,        &
      &                             conf = me%general%solver%conf(1) )

  end subroutine atl_load_solver_parameters

end module atl_solver_param_module
