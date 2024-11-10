! Copyright (c) 2018 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> The precice module contains basic routines to interact with the precice
!! coupling tool.
module atl_precice_module

  use env_module,                    only: rk, labelLen
  use treelmesh_module,              only: treelmesh_type
  use tem_tools_module,              only: tem_horizontalSpacer
  use tem_logging_module,            only: logUnit
  use tem_timer_module,              only: tem_startTimer, tem_stopTimer
  use tem_precice_module,            only: tem_precice_initialize_data,       &
    &                                      tem_precice_fulfilled_Action,      &
    &                                      tem_precice_init,                  &
    &                                      tem_precice_action_write_initData, &
    &                                      tem_precice_action_required

  use atl_solver_param_module,       only: atl_solver_param_type
  use atl_equation_module,           only: atl_equations_type
  use atl_cube_elem_module,          only: atl_cube_elem_type
  use atl_timer_module,              only: atl_timerHandles
  use atl_writePrecice_module,       only: atl_write_precice_getPoint, &
     &                                     atl_write_precice,          &
     &                                     atl_read_precice

  implicit none

  private

  public :: atl_init_precice

contains

  ! ************************************************************************ !
  subroutine atl_init_precice( params, equation, tree, mesh_list, &
    &                          timestepLimit                      )
    ! -------------------------------------------------------------------- !
    type(atl_solver_param_type), intent(inout) :: params

    !> Description on the equation system to solve.
    type(atl_Equations_type), intent(inout), target :: equation

    !> Mesh data in treelmesh format.
    type(treelmesh_type), intent(in) :: tree

    type(atl_cube_elem_type), intent(in) ::               &
     & mesh_list(tree%global%minLevel:tree%global%maxLevel)

    !> timesteplimitation due to the precice timestep
    real(kind=rk), intent(inout) :: timestepLimit
    ! -------------------------------------------------------------------- !
    character(len=labelLen) :: nameaction
    integer :: action_isrequired
    ! -------------------------------------------------------------------- !
    call tem_startTimer( timerHandle = atl_timerHandles%preciceInit )
    write(logUnit(1),*) 'Initialize preCICE'

    call atl_write_precice_getPoint(    &
      & stFunlist = equation%stFunlist, &
      & tree      = tree,               &
      & mesh_list = mesh_list,          &
      & varSys    = equation%varSys     )

    !> Read from preCICE the requested points
    call atl_read_precice(              &
      & stFunList = equation%stFunList, &
      & tree      = tree                )

    !> If coupling with precice is used, precice need to be initilized
    !! to calcualte e.g. the mapping, and communication and the
    !! precice_timestep loaded from the config file as input
    !! need to be called after precice_create and set?read/write/mesh positions
    !! to precice
    !! and before initial condition is writen to precice
    call tem_horizontalSpacer(fUnit = logUnit(1))

    call tem_precice_init(timesteplimit = timestepLimit)
    call tem_horizontalSpacer(fUnit = logUnit(1))

    !> get the internal precice name for the action writining initial data
    call tem_precice_action_write_initData(nameAction = Nameaction )
    !> with this name we check is this action is required
    call tem_precice_action_required( Nameaction = Nameaction,       &
      &                               isRequired = action_isRequired )

    !> check if action is required and if so, write initial condition to precice
    if (action_isRequired .eq. 1) then
      write(logUnit(1),*) 'Initialize data for preCICE'
      call atl_write_precice( stFunList = equation%stFunList,            &
        &                     varSys    = equation%varSys,               &
        &                     time      = params%general%simControl%now, &
        &                     tree      = tree                           )

       ! call to precice -->  fulfill the action to signal that these data was
       ! the initial data for the simulation
       call tem_precice_fulfilled_action(nameAction = nameAction)
    end if

    !> call to precice to exchange the initialize Data
    call tem_precice_initialize_data()

    write(logUnit(1),*) 'Done initializing preCICE'
    call tem_stopTimer( timerHandle = atl_timerHandles%preciceInit )

  end subroutine atl_init_precice
  ! ************************************************************************ !

end module atl_precice_module
