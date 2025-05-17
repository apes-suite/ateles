! Copyright (c) 2011-2016, 2018-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2011-2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011-2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> ATELES
!!
!! This is the main program of the Adaptive Tree-based Efficient and Lithe
!! Equation Solver.
!! It provides the Modal Discontinuous Galerkin solver in the APES framework.
!! (c) 2011 German Research School for Simulation Sciences GmbH
!! (c) 2013-2016 University of Siegen
!!
program ateles

  use mpi

  use env_module,              only: rk
  use tem_general_module,      only: tem_start, tem_finalize
  use tem_precice_module,      only: precice_available
  use treelmesh_module,        only: treelmesh_type

  use atl_timer_module,        only: atl_addTimers
  use atl_aux_module,          only: atl_banner
  use atl_solver_param_module, only: atl_solver_param_type
  use atl_container_module,    only: atl_element_container_type
  use atl_equation_module,     only: atl_equations_type
  use atl_varSys_module,       only: atl_varSys_solverData_type
  use atl_program_module,      only: atl_initialize_program, &
    &                                atl_solve_program,      &
    &                                atl_finalize_program,   &
    &                                atl_load_config

  use ply_poly_project_module, only: ply_poly_project_type

  implicit none

  ! ------------------------------------------------------------------------ !
  ! The structure that holds the solver parameter
  type(atl_solver_param_type) :: params

  ! Description of the equation system to solve
  type(atl_equations_type) :: equation

  ! The treelmesh data structure
  type(treelmesh_type) :: tree

  ! Data Infomation of the variable System
  type(atl_varSys_solverData_type) :: varSys_data

  ! Number of cells on each levels
  integer, allocatable :: nCellsNoBnd(:)

  ! Main data structure of Ateles describing the mesh elements
  type(atl_element_container_type) :: element_container

  ! Desribe the projetion methods for the polynomials
  type(ply_poly_project_type), allocatable :: poly_proj_list(:)

  ! Timestep specified from precice
  real(kind=rk) :: precice_dt

  integer :: iError
  ! ------------------------------------------------------------------------ !
  ! Initialize the treelm environment
  call tem_start( codeName   = 'ATELES',                 &
    &             general    = params%general,           &
    &             simcontrol = params%general%simControl )

  if (params%general%proc%rank == 0) then
    call atl_banner(trim(params%general%solver%version))
  end if

  ! init the timer
  call atl_addTimers()

  ! load the config file
  call atl_load_config(params = params, tree = tree)

  if (.not. precice_available) then
    ! call initialize of the program without optional precice variable
    call atl_initialize_program(               &
      & params            = params,            &
      & equation          = equation,          &
      & tree              = tree,              &
      & varSys_data       = varSys_data,       &
      & nCellsNoBnd       = nCellsNoBnd,       &
      & element_container = element_container, &
      & poly_proj_list    = poly_proj_list     )
    call MPI_Barrier(MPI_COMM_WORLD, iError)

    ! this is function for solve the program including calc_timestep, meshstep.
    call atl_solve_program( params            = params,            &
      &                     equation          = equation,          &
      &                     tree              = tree,              &
      &                     nCellsNoBnd       = nCellsNoBnd,       &
      &                     element_container = element_container, &
      &                     poly_proj_list    = poly_proj_list     )
  else
    ! call initialize of the program with precice variable
    call atl_initialize_program(               &
      & params            = params,            &
      & equation          = equation,          &
      & tree              = tree,              &
      & varSys_data       = varSys_data,       &
      & nCellsNoBnd       = nCellsNoBnd,       &
      & element_container = element_container, &
      & poly_proj_list    = poly_proj_list,    &
      & precice_dt        = precice_dt         )
    call MPI_Barrier(MPI_COMM_WORLD, iError)

    ! this is function for solve the program including calc_timestep, meshstep.
    call atl_solve_program(                    &
      & params            = params,            &
      & equation          = equation,          &
      & tree              = tree,              &
      & nCellsNoBnd       = nCellsNoBnd,       &
      & element_container = element_container, &
      & poly_proj_list    = poly_proj_list,    &
      & precice_dt        = precice_dt         )
  end if

  ! finialize ateles with veything output related, finialize treelm
  call atl_finalize_program( params            = params,           &
      &                      equation          = equation,         &
      &                      tree              = tree,             &
      &                      nCellsNoBnd       = nCellsNoBnd,      &
      &                      element_container = element_container )

  call tem_finalize(params%general)

end program ateles
