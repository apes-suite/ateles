! Copyright (c) 2011-2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2017,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
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
!! This module provides the definition and methods for time step computations.
module atl_calc_time_module
  use mpi
  use env_module,                            only: rk, rk_mpi
  use tem_aux_module,                        only: tem_abort
  use tem_float_module,                      only: operator(.fgt.), &
    &                                              operator(.flt.)
  use tem_logging_module,                    only: logUnit
  use tem_comm_env_module,                   only: tem_comm_env_type
  use tem_general_module,                    only: tem_general_type
  use treelmesh_module,                      only: treelmesh_type
  use tem_isNaN_module,                      only: tem_isNaN
  use tem_time_module,                       only: tem_time_type
  use tem_precice_module,                    only: precice_handle

  use atl_equation_module,                   only: atl_equations_type
  use atl_kerneldata_module,                 only: atl_statedata_type
  use atl_elemental_time_integration_module, only: atl_timestep_type
  use atl_scheme_module,                     only: atl_scheme_type,        &
    &                                              atl_modg_scheme_prp,    &
    &                                              atl_modg_2d_scheme_prp, &
    &                                              atl_modg_1d_scheme_prp
  use atl_materialPrp_module,                only: atl_material_type
  use atl_eqn_linearEuler_module,            only: atl_eqn_update_background
  use atl_cube_elem_module,                  only: atl_cube_elem_type
  use atl_time_integration_module,           only: atl_global_timestep_type

  implicit none

  private

  public :: atl_get_timestep

contains


  ! ************************************************************************ !
  !> Subrountine which gather all calls to get the timestep for the current
  !! iteration
  subroutine atl_get_timestep( tree, mesh_list, scheme_list, material_list,    &
    &                          equation, time, statedata_list, nCellsNoBnd,    &
    &                          general, adaptive_timestep, initial, precice_dt )
    ! -------------------------------------------------------------------- !
    !> The treelmesh data structure
    type(treelmesh_type), intent(in) :: tree

    !> List of meshes for different kernels
    type(atl_cube_elem_type), intent(in) :: mesh_list(tree%global%minLevel: &
      &                                               tree%global%maxLevel  )

    !> scheme desription
    type(atl_scheme_type), intent(inout) :: scheme_list(tree%global%minLevel: &
      &                                                 tree%global%maxLevel  )

    !> List of material parameter information for the mesh. One entry for
    !! level, running from minlevel to maxlevel.
    type(atl_material_type), intent(in) :: material_list(tree%global%minLevel: &
      &                                                  tree%global%maxLevel  )

    !> Global timediscretization type
    type(atl_global_timestep_type), intent(in) :: time

    !> Local time
    type(atl_statedata_type), intent(in)                             &
      & :: statedata_list(tree%global%minLevel:tree%global%maxLevel  )

    ! Description of the equation system to solve
    type(atl_equations_type), intent(inout) :: equation

    ! Number of cells on each levels
    integer, intent(in) :: nCellsNoBnd(:)

    !> general data coming from treelem
    type(tem_general_type), intent(in) :: general

    !> Flag for adaptive timestep calculation
    logical, intent(in) :: adaptive_timestep

    !> Flag for timestep calculation because it is the init step
    logical, intent(in) :: initial

    !> Timestep specified from precice. If this value is present, precice is
    !! considered active.
    real(kind=rk), intent(out), optional :: precice_dt
    ! -------------------------------------------------------------------- !
    !> Iterating over the various levels
    integer :: iList
    !> Global timestep introduced for precice
    real(kind=rk) :: glob_dt
    !> variable needed for truncation of last time step
    real(kind=rk) :: time2end
    !>parameter needed to fix the subcycling for the couplinging in precice
    real(kind=rk) ::para_phi
    ! -------------------------------------------------------------------- !

    ! for some equation systems the timestep is not changing over simulation
    ! time, hence there is no need to calculate the timestep and get the global
    ! timestep in every iteration also for the initialize step the calcualtion
    ! need to be done

    if (time%control%fixed_dt > 0.0_rk) then
      !> There is an fixed_dt requested, use that instead of adaptive
      !! timestepping.
      scheme_list(:)%time%dt = time%control%fixed_dt
      write(logUnit(2),*) "Writing fixed time step: ", time%control%fixed_dt
    else

      if ( (adaptive_timestep) .or. (initial) )then
        ! Calc timestep length for this iteration on each level.
        ! This could also give a local timestep for each level
        do iList = tree%global%minLevel, tree%global%maxLevel
          call calculate_cfl_timestep(                     &
            & length    = mesh_list(iList)%length,         &
            & cfl       = time%control%cfl,                &
            & cfl_visc  = time%control%cfl_visc,           &
            & equation  = equation,                        &
            & dt        = scheme_list(ilist)%time%dt,      &
            & timestep  = time%elementSteps(iList),        &
            & material  = material_list(iList),            &
            & scheme    = scheme_list(iList),              &
            & localtime = statedata_list(iList)%local_time )
        end do

        ! No local timestep so far, so we have to use the smallest one.
        ! Find the smallest time step across all levels and partitions.
        ! (involves MPI communication)
        call create_global_timestep(           &
          & dt       = scheme_list(:)%time%dt, &
          & nCells   = nCellsNoBnd,            &
          & minlevel = tree%global%minLevel,   &
          & maxlevel = tree%global%maxLevel,   &
          & proc     = general%proc            )

        ! Check for nan timestep size, and abort, if one is found.
        if ( tem_isnan(scheme_list(tree%global%minlevel)%time%dt) ) then
          call tem_abort( 'ERROR: detected NAN timestep, stopping ...' )
        end if
      end if

    end if


    ! checking for precice timestep as well as trunctaction of the timestep due
    ! to overall simulation time need to be done in every iteration
    glob_dt = minval(scheme_list(:)%time%dt)
    ! when coupling via precice is used, there are two time steps:
    ! one specified by the cfl condition and one by the delta t to the next sync
    ! timestep need to define which timestep to use
    if ( present(precice_dt) ) then
      if (precice_handle%use_RK2_inter) then
        write(*,*) "Making use of intermediate timestep"
        write(logUnit(2),*) 'Ateles timestep controlled by cfl number :', &
          &                 glob_dt
        write(logUnit(2),*) 'preCICE timestep precice_dt', precice_dt
        glob_dt = glob_dt
        write(logUnit(2),*) 'global timestep', glob_dt
      else
        write(logUnit(2),*) 'Ateles timestep controlled by cfl number :', &
          &                 glob_dt
        write(logUnit(2),*) 'preCICE timestep precice_dt', precice_dt
        glob_dt = min (glob_dt, precice_dt)
        write(logUnit(2),*) 'Actually used timestep in solver glob_dt=', glob_dt
      end if
    end if
    ! Tthe parameter para_phi helps us to aviod small time steps at the end of
    ! the simulation, which leads to subcycling. Thus, we check wether time2end
    ! is only slightly bigger than our dt. If this is the case, than our
    ! global_dt is equal to time2end otherwise 2 half time steps are done.
    para_phi = 1.000001
    ! Truncate the last time step if the final time is reached.
    time2end = general%simControl%timeControl%max%sim &
      &          - general%simControl%now%sim
    if ( (time2end .flt. (1.5_rk*glob_dt)) &
      &  .and. (time2end .fgt. para_phi*glob_dt)  ) then
      glob_dt = 0.5_rk * time2end
    end if
    if (time2end .flt. para_phi*glob_dt) then
      glob_dt = time2end
    end if
    scheme_list(:)%time%dt = glob_dt


  end subroutine atl_get_timestep
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate the timestep for a whole part of a cubic mesh by a CFL condition.
  !!
  !! This routine takes the primitive variables of the equation system and
  !! calculates the next timestep.
  !! The calculation is based on the cfl condition and the next restart
  !! timepoints.
  subroutine calculate_cfl_timestep( length, cfl, cfl_visc, equation, dt,  &
    &                                timestep, scheme, material, localtime )
    ! -------------------------------------------------------------------- !

    !> The length of the cubes you are calculating the cfl condition for.
    real(kind=rk), intent(in) :: length

    !> The CFL factor to apply (for the convective part of the equation).
    real(kind=rk), intent(in) :: cfl

    !> The CFL factor to apply (for the viscous part of the equation).
    real(kind=rk), intent(in) :: cfl_visc

    !> The equation system to be used in the simulation.
    type(atl_equations_type), intent(inout) :: equation

    !> Resulting timestep.
    real(kind=rk), intent(out) :: dt

    !> Timestep information
    type(atl_timestep_type), intent(in) :: timestep

    !> Info about the scheme.
    type(atl_scheme_type), intent(in) :: scheme

    !> Material information for all elements on the current level
    type(atl_material_type), intent(in) :: material

    !> Local time, required for update of temporal background in linear Euler
    type(tem_time_type), intent(in) :: localtime
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: speedOfLight
    integer :: nPoly
    ! -------------------------------------------------------------------- !

    nPoly = 0 ! default value
    ! Get the polynomial nPoly
    select case(scheme%scheme)
    case(atl_modg_scheme_prp)
      nPoly = scheme%modg%maxPolyDegree + 1
    case(atl_modg_2d_scheme_prp)
      nPoly = scheme%modg_2d%maxPolyDegree + 1
    case(atl_modg_1d_scheme_prp)
      nPoly = scheme%modg_1d%maxPolyDegree + 1
    case default
      call tem_abort(                                                 &
        & 'ERROR in calc_time: Unknown spatial scheme, stopping ... ' )
    end select

    ! We calculate the new timestep based on the cfl condition with
    ! a given cfl-coefficient limit. We include fluid, ghost, halo
    ! and boundary cells here.
    !
    ! Select the timestep computation based upon the equation system to solve.
    select case(trim(equation%eq_kind))
    case('navier_stokes', 'filtered_navier_stokes')
      call calc_timestep_viscflow_cube(                 &
        & cfl_conv   = cfl,                             & ! in
        & cfl_visc   = cfl_visc,                        & ! in
        & length     = length,                          & ! in
        & dt         = dt,                              & ! out
        & timestep   = timestep,                        & ! in
        & nPoly      = nPoly,                           & ! in
        & mu         = equation%NavierStokes%mu,        & ! in
        & therm_cond = equation%NavierStokes%therm_cond ) ! in

    case('euler')
      call calc_timestep_flow_cube(cfl      = cfl,      & ! in
        &                          length   = length,   & ! in
        &                          dt       = dt,       & ! out
        &                          timestep = timestep, & ! in
        &                          nPoly    = nPoly     )

    case('loclineuler','loclineuler_1d')
      call calc_timestep_flow_cube_mod(cfl      = cfl,        & !in
        &                              length   = length,     & !in
        &                              dt       = dt,         & !out
        &                              timestep = timestep,   & !in
        &                              nPoly    = nPoly       ) !in

    case('lineareuler')
      ! update background
      call atl_eqn_update_background(         &
        & me          = equation%linearEuler, &
        & time        = localtime,            &
        & nDimensions = 3                     )

      call calc_timestep_linearEuler_cube(                  &
        & cfl          = cfl,                               & ! in
        & length       = length,                            & ! in
        & vel          = equation%linearEuler%velocity_0,   & ! in
        & SpeedofSound = equation%linearEuler%SpeedOfSound, & ! in
        & dt           = dt,                                & ! out
        & nPoly        = nPoly                              )

    case('acoustic')
      call calc_timestep_acoustic_cube(                  &
        & cfl          = cfl,                            & ! in
        & length       = length,                         & ! in
        & vel          = equation%acoustic%velocity_0,   & ! in
        & SpeedofSound = equation%acoustic%SpeedOfSound, & ! in
        & dt           = dt,                             & ! out
        & nPoly        = nPoly                           )

    case('acoustic_2d')
      call calc_timestep_acoustic_2d_cube(               &
        & cfl          = cfl,                            & ! in
        & length       = length,                         & ! in
        & vel          = equation%acoustic%velocity_0,   & ! in
        & SpeedofSound = equation%acoustic%SpeedOfSound, & ! in
        & dt           = dt,                             & ! out
        & nPoly        = nPoly                           )


    case('lineareuler_2d')
      ! update background
      call atl_eqn_update_background(         &
        & me          = equation%linearEuler, &
        & time        = localtime,            &
        & nDimensions = 2                     )

      call calc_timestep_linearEuler_2d_cube(               &
        & cfl          = cfl,                               & ! in
        & length       = length,                            & ! in
        & vel          = equation%linearEuler%velocity_0,   & ! in
        & SpeedofSound = equation%linearEuler%SpeedOfSound, & ! in
        & dt           = dt,                                & ! out
        & nPoly        = nPoly                              )

    case('euler_2d')
      call calc_timestep_flow_cube_2d( cfl      = cfl,      & ! in
        &                              length   = length,   & ! in
        &                              dt       = dt,       & ! out
        &                              timestep = timestep, & ! in
        &                              nPoly    = nPoly     )

    case('navier_stokes_2d', 'filtered_navier_stokes_2d')
      call calc_timestep_viscflow_cube_2d(              &
        & cfl_conv   = cfl,                             & ! in
        & cfl_visc   = cfl_visc,                        & ! in
        & length     = length,                          & ! in
        & dt         = dt,                              & ! out
        & timestep   = timestep,                        & ! in
        & nPoly      = nPoly,                           & ! in
        & mu         = equation%NavierStokes%mu,        & ! in
        & therm_cond = equation%NavierStokes%therm_cond ) ! in

    case('euler_1d')
      call calc_timestep_flow_cube_1d( cfl      = cfl,      & ! in
        &                              length   = length,   & ! in
        &                              dt       = dt,       & ! out
        &                              timestep = timestep, & ! in
        &                              nPoly    = nPoly     )

    case('maxwell','maxwell_2d')
      ! Read out the fastest wave speed as provided by the material description
      speedOfLight = material%maxPropSpeed
      ! Use a CFL condition to determine the next timestep
      call calc_timestep_ed_cube( cfl          = cfl,          & ! in
        &                         length       = length,       & ! in
        &                         speedOfLight = speedOfLight, & ! in
        &                         dt           = dt,           &
        &                         nPoly        = nPoly         )

    case('maxwelldivcorrection')
      ! Read out the fastest wave speed as provided by the material description
      speedOfLight = material%maxPropSpeed
      ! Use a CFL condition to determine the next timestep
      call calc_timestep_ed_cube( cfl          = cfl,          & ! in
        &                         length       = length,       & ! in
        &                         speedOfLight = speedOfLight, & ! in
        &                         dt           = dt,           & ! out
        &                         nPoly        = nPoly         )

    case('nernstplanck')
      call calc_timestep_nerplanck_cube( cfl    = cfl,    & ! in
        &                                length = length, & ! in
        &                                dt     = dt,     & ! out
        &                                nPoly  = nPoly   )

    case('advection_1d')
      call calc_timestep_adv_cube( cfl    = cfl,                         & ! in
        &                          length = length,                      & ! in
        &                          vel    = equation%advection%velocity, & ! in
        &                          dt     = dt,                          & ! out
        &                          nPoly  = nPoly                        )

    case('heat_1d','heat_2d','heat')
      call calc_timestep_heat_cube_1d( cfl      = cfl_visc, & ! in
        &                              length   = length,   & ! in
        &                              dt       = dt,       & ! out
        &                              equation = equation, & ! in
        &                              nPoly    = nPoly     )

    case default
      write(logUnit(1),*) 'Equation kind: ' // trim(equation%eq_kind)
      write(logUnit(1),*) 'Timestep calculation does not support this' &
        & // ' equation, stopping...'
      call tem_abort()
    end select

  end subroutine calculate_cfl_timestep
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a flow
  !! simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  subroutine calc_timestep_heat_cube_1d(cfl, length, dt, equation, nPoly )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in)  :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in)  :: length
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> The number of polynomials per spatial direction The equation system to
    !! be used in the simulation.
    integer, intent(in) :: nPoly
    type(atl_equations_type), intent(in) :: equation
    integer, parameter :: ik = selected_int_kind(16)
    integer(kind=ik)  :: npow
    ! -------------------------------------------------------------------- !

    npow = nPoly
    npow = npow**4
    ! Calculate timestep from CFL-condition
    dt = 0.07* (cfl * length * length) / (equation%heat%k*nPow)

  end subroutine calc_timestep_heat_cube_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL number of cubes in a
  !! advection simulation.
  !!
  !! This subroutine calculates the time step with respect to a given cfl
  !! number for simulations of advection processes.
  !! Please notice that this routine assumes constant advection speed.
  subroutine calc_timestep_adv_cube(cfl, length, vel, dt, nPoly)
    ! ---------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Upper bound of the advection velocity
    real(kind=rk), intent(in) :: vel
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! ---------------------------------------------------------------------- !

    dt = cfl * ( length / (vel * nPoly**2) )

  end subroutine calc_timestep_adv_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL number of cubes in a
  !! electrodynamic simulation.
  !!
  !! This subroutine calculates the time step with respect to a given cfl
  !! number for simulations of electrodynamic processes.
  !! Please notice that this routine assumes constant material parameters.
  subroutine calc_timestep_ed_cube(cfl, length, speedOfLight, dt, nPoly)
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Speed of light
    real(kind=rk), intent(in) :: speedOfLight
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !

    dt = cfl * ( length / (speedOfLight * nPoly**2) )

  end subroutine calc_timestep_ed_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL number of cubes in a
  !! Nernst-Planck simulation.
  !!
  !! This subroutine calculates the time step with respect to a given cfl
  !! number for simulations of Nernst-Planck equations.
  !! Please notice that this routine assumes constant material parameters.
  subroutine calc_timestep_nerplanck_cube(cfl, length, dt, nPoly)
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !

    dt = cfl * length / (nPoly**2)

  end subroutine calc_timestep_nerplanck_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a linear
  !! Euler simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for linear Euler
  !! simulations only.
  subroutine calc_timestep_linearEuler_cube( cfl, length, vel, SpeedOfSound, &
    &                                        dt, nPoly                       )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Background velocity in the domain (x,y,z)
    real(kind=rk), intent(in) :: vel(3)
    !> Speed of sound, based on background density and pressure
    real(kind=rk), intent(in) :: SpeedofSound
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !
    real(kind=rk)  :: max_vel
    ! -------------------------------------------------------------------- !

    max_vel = sqrt(sum(vel**2))
    ! Calculate timestep from CFL-condition
    dt =  cfl * length / abs(max_vel + SpeedOfSound) / 2._rk / (nPoly**2)

  end subroutine calc_timestep_linearEuler_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a linear
  !! Euler 2d simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for linear Euler 2d
  !! simulations only.
  subroutine calc_timestep_linearEuler_2d_cube( cfl, length, vel,        &
    &                                           SpeedOfSound, dt,  nPoly )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Background velocity in the domain (x,y,z)
    real(kind=rk), intent(in) :: vel(2)
    !> Speed of sound, based on background density and pressure
    real(kind=rk), intent(in) :: SpeedofSound
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !
    real(kind=rk)  :: max_vel
    ! -------------------------------------------------------------------- !

    max_vel = sqrt(sum(vel**2))
    ! Calculate timestep from CFL-condition
    dt =  cfl * length / abs(max_vel + SpeedOfSound) / 2._rk / (nPoly**2)

  end subroutine calc_timestep_linearEuler_2d_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a acoustic
  !! simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for acoustic simulations
  !! only.
  subroutine calc_timestep_acoustic_2d_cube( cfl, length, vel, SpeedOfSound, &
    &                                        dt,  nPoly )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Background velocity in the domain (x,y,z)
    real(kind=rk), intent(in) :: vel(2)
    !> Speed of sound, based on background density and pressure
    real(kind=rk), intent(in) :: SpeedofSound
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !
    real(kind=rk)  :: max_vel
    ! -------------------------------------------------------------------- !

    max_vel = sqrt(sum(vel**2))
    ! Calculate timestep from CFL-condition
    dt =  cfl * length / abs(max_vel + SpeedOfSound) / 2._rk / (nPoly**2)

  end subroutine calc_timestep_acoustic_2d_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a acoustic
  !! simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for acoustic simulations
  !! only.
  subroutine calc_timestep_acoustic_cube( cfl, length, vel, SpeedOfSound, &
    &                                     dt,  nPoly )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Background velocity in the domain (x,y,z)
    real(kind=rk), intent(in) :: vel(3)
    !> Speed of sound, based on background density and pressure
    real(kind=rk), intent(in) :: SpeedofSound
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !
    real(kind=rk)  :: max_vel
    ! -------------------------------------------------------------------- !

    max_vel = sqrt(sum(vel**2))
    ! Calculate timestep from CFL-condition
    dt =  cfl * length / abs(max_vel + SpeedOfSound) / 2._rk / (nPoly**2)

  end subroutine calc_timestep_acoustic_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a flow
  !! simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for flow simulations only.
  subroutine calc_timestep_flow_cube( cfl, length, dt, timestep, nPoly )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> Info about the timestep type
    type(atl_timestep_type), intent(in) :: timestep
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !
    real(kind=rk)  :: max_velocity
    ! -------------------------------------------------------------------- !

    max_velocity = maxval(                                                  &
      & abs(timestep%euler%maxVel(:)) + abs(timestep%euler%speedOfSound(:)) )

    ! Calculate timestep from CFL-condition
    dt =  cfl * length / max_velocity / 2._rk / (nPoly**2)

  end subroutine calc_timestep_flow_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine calc_timestep_flow_cube_mod( cfl, length, dt, timestep, nPoly )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> Info about the timestep type
    type(atl_timestep_type), intent(in) :: timestep
    ! The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !
    real(kind=rk)  :: max_velocity
    ! -------------------------------------------------------------------- !

    max_velocity =  maxval( abs(timestep%LoclinEuler%meanVel(:))          &
      &                       + abs(timestep%LoclinEuler%speedOfSound(:)) )
    ! Calculate timestep from CFL-condition
    dt =  cfl * length / max_velocity / 2._rk / (nPoly**2)

  end subroutine calc_timestep_flow_cube_mod
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a flow
  !! simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for flow simulations only.
  subroutine calc_timestep_flow_cube_2d( cfl, length, dt, timestep, nPoly )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> Info about the timestep type
    type(atl_timestep_type), intent(in) :: timestep
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !
    real(kind=rk)  :: max_velocity
    ! -------------------------------------------------------------------- !

    max_velocity =  maxval( abs(timestep%euler_2d%maxVel(:))           &
      &                       + abs(timestep%euler_2d%speedOfSound(:)) )
    ! Calculate timestep from CFL-condition
    dt =  cfl * length / max_velocity / 2._rk / (nPoly**2)

  end subroutine calc_timestep_flow_cube_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a viscous
  !! flow simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for flow simulations only.
  subroutine calc_timestep_viscflow_cube( cfl_conv, cfl_visc, length, dt, &
    &                                     timestep, nPoly, mu, therm_cond )
    ! -------------------------------------------------------------------- !
    !> CFL number for the convective part
    real(kind=rk), intent(in) :: cfl_conv
    !> CFL number for the viscous part
    real(kind=rk), intent(in) :: cfl_visc
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> Info about the timestep type
    type(atl_timestep_type), intent(in) :: timestep
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    !> The dynamic viscosity
    real(kind=rk), intent(in) :: mu
    !> The thermal conductivity
    real(kind=rk), intent(in) :: therm_cond
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: max_velocity, dt_conv, dt_visc
    real(kind=rk) :: disc_fact
    ! -------------------------------------------------------------------- !

    max_velocity =  maxval( abs(timestep%euler%maxVel(:))           &
      &                       + abs(timestep%euler%speedOfSound(:)) )

    disc_fact = 0.5_rk * length / (nPoly**2)

    ! Calculate timestep from CFL-condition for the convective part
    ! of the Navier-Stokes equations
    dt_conv =  cfl_conv * disc_fact / max_velocity


    dt_visc = cfl_visc * disc_fact**2   &
      &       / ( max( mu, therm_cond ) )

    dt = min(dt_conv, dt_visc)

  end subroutine calc_timestep_viscflow_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a viscous
  !! flow simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for flow simulations only.
  subroutine calc_timestep_viscflow_cube_2d( cfl_conv, cfl_visc, length, dt, &
    &                                        timestep, nPoly, mu, therm_cond )
    ! -------------------------------------------------------------------- !
    !> CFL number for the convective part
    real(kind=rk), intent(in) :: cfl_conv
    !> CFL number for the viscous part
    real(kind=rk), intent(in) :: cfl_visc
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> Info about the timestep type
    type(atl_timestep_type), intent(in) :: timestep
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    !> The dynamic viscosity
    real(kind=rk), intent(in) :: mu
    ! The thermal conductivity
    real(kind=rk), intent(in) :: therm_cond
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: max_velocity, dt_conv, dt_visc
    real(kind=rk) :: disc_fact
    ! -------------------------------------------------------------------- !

    max_velocity =  maxval( abs(timestep%euler_2d%maxVel(:))           &
      &                       + abs(timestep%euler_2d%speedOfSound(:)) )

    disc_fact = 0.5_rk * length / (nPoly**2)

    ! Calculate timestep from CFL-condition for the convective part
    ! of the Navier-Stokes equations
    dt_conv =  cfl_conv * disc_fact / max_velocity


    dt_visc = cfl_visc * disc_fact**2   &
      &       / ( max( mu, therm_cond ) )

    dt = min(dt_conv, dt_visc)

  end subroutine calc_timestep_viscflow_cube_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate time step based on a given CFL condition for a cube in a flow
  !! simulation.
  !!
  !! This subroutine calculates the timstep according to the CFL condition.
  !! Please notice, that this routine can be applied for flow simulations only.
  subroutine calc_timestep_flow_cube_1d(cfl, length, dt, timestep, nPoly )
    ! -------------------------------------------------------------------- !
    !> CFL number
    real(kind=rk), intent(in) :: cfl
    !> Reference length of all elements
    real(kind=rk), intent(in) :: length
    !> Resulting time step width
    real(kind=rk), intent(out) :: dt
    !> Info about the timestep type
    type(atl_timestep_type), intent(in) :: timestep
    !> The number of polynomials per spatial direction
    integer, intent(in) :: nPoly
    ! -------------------------------------------------------------------- !
    real(kind=rk)  :: max_velocity
    ! -------------------------------------------------------------------- !

    max_velocity =  maxval( abs(timestep%euler_1d%maxVel(:))           &
      &                       + abs(timestep%euler_1d%speedOfSound(:)) )

    ! Calculate timestep from CFL-condition
    dt =  cfl * length / max_velocity / 2._rk / (nPoly**2)

  end subroutine calc_timestep_flow_cube_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> \brief subroutine to create a single global timestep.
  subroutine create_global_timestep( dt, nCells, minLevel, maxLevel, proc )
    ! -------------------------------------------------------------------- !
    !> The minimum level of your mesh.
    integer, intent(in) :: minLevel

    !> The maximum level of your mesh.
    integer, intent(in) :: maxLevel

    !> The delta t on the levels of your mesh. This will be update with
    !! the new global, local timestep.
    real(kind=rk), intent(inout) :: dt(minLevel:maxLevel)

    !> The number of cells on the different levels.
    integer,intent(in) :: nCells(minLevel:maxLevel)
    !> mpi communication enviroment with mpi communicator
    type(tem_comm_env_type), intent(in) :: proc
    ! -------------------------------------------------------------------- !
    integer :: iList
    real(kind=rk), allocatable :: local_dt(:)
    ! -------------------------------------------------------------------- !

    allocate(local_dt(minLevel:maxLevel))
    do iList = minLevel, maxLevel
      if(nCells(iList).gt.0.0_rk) then
        local_dt(iList) = dt(iList)
      else
        local_dt(iList) = huge(dt(iList))
      end if
    end do
    local_dt(minLevel) = calc_common_global_timestep(local_dt, proc)
    do iList = minLevel, maxLevel
      dt(iList) = local_dt(minLevel)
    end do

    deallocate(local_dt)

  end subroutine create_global_timestep
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Function to find a single global time step for all levels and processes.
  function calc_common_global_timestep(local_dt, proc) result(dt)
    ! -------------------------------------------------------------------- !
    !> Process local time steps on each level.
    real(kind=rk), allocatable, intent(inout) :: local_dt(:)
    !> mpi communication enviroment with mpi communicator
    type(tem_comm_env_type), intent(in) :: proc

    !> Resulting global time step for all processes and levels.
    real(kind=rk) :: dt
    ! -------------------------------------------------------------------- !
    integer :: ierr
    real(kind=rk) :: local_min
    ! -------------------------------------------------------------------- !

    local_min = minval(local_dt)

    call mpi_allreduce(local_min, dt, 1, rk_mpi, mpi_min, proc%comm, ierr)

  end function calc_common_global_timestep
  ! ************************************************************************ !

end module atl_calc_time_module
