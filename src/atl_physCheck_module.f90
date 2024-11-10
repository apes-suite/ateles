! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Parid Ndreka
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
!! Module provides routines and datatypes to check/correct the physical
!! values of a given state.
module atl_physCheck_module
  ! Treelm modules
  use env_module,                   only: rk, long_k
  use tem_tools_module,             only: tem_horizontalSpacer
  use tem_aux_module,               only: tem_abort
  use tem_status_module,            only: tem_stat_nonPhysical, tem_status_type
  use tem_isNaN_module,             only: tem_isNaN
  use tem_logging_module,           only: logUnit
  use tem_element_module,           only: eT_fluid

  ! Aotus module
  use aotus_module,                 only: flu_State, aot_get_val
  use aot_table_module,             only: aot_table_open, aot_table_close

  ! Ateles module
  use atl_kerneldata_module,        only: atl_statedata_type
  use atl_scheme_module,            only: atl_scheme_type,        &
    &                                     atl_modg_2d_scheme_prp, &
    &                                     atl_modg_scheme_prp,    &
    &                                     atl_modg_1d_scheme_prp
  use atl_equation_module,          only: atl_equations_type
  use atl_eqn_euler_module,         only: atl_Euler_type
  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_time_integration_module,  only: atl_global_timestep_type

  use ply_poly_project_module,      only: ply_poly_project_type, &
    &                                     assignment(=),         &
    &                                     ply_poly_project_m2n
  use ply_oversample_module,        only: ply_convert2oversample


  implicit none
  private

  !> Datatype to describe the physical checks.
  type atl_physCheck_type
    !> Is the check active
    logical :: isActive = .false.
    !> The interval for the physical checks.
    integer :: iterationInterval = huge(int(1))
    !> The tolerance below which a state will be marked as unphysical
    real(kind=rk) :: tolerance = 1.0e-13_rk
  end type atl_physCheck_type

  public :: atl_check_val, atl_physCheck_type, atl_init_physCheck


contains

  !> summary: Read the info for the physical checks from the configuration file
  subroutine atl_init_physCheck( check, conf, equation )
    ! ---------------------------------------------------------------------------
    !> The physical check to init.
    type(atl_physCheck_type), intent(out) :: check
    !> Handle for the Lua config file
    type(flu_State), intent(in) :: conf
    !> The equation you are simulating.
    type(atl_equations_type),intent(in) :: equation
    ! ---------------------------------------------------------------------------
    integer :: check_table, iError
    ! ---------------------------------------------------------------------------

    check%isActive = .false.

    call aot_table_open(L=conf, thandle=check_table, key='check')
    if(check_table.eq.0) then
      return
    end if

    call tem_horizontalSpacer(fUnit=logUnit(1))
    write(logUnit(1),*) 'Check for physical values is active:'

    ! Read the check interval
    call aot_get_val(L = conf, thandle = check_table, &
      &              key = 'interval', &
      &              val = check%iterationInterval, &
      &              ErrCode = iError)

    if(iError.ne.0) then
      call tem_abort( 'ERROR in atl_init_physCheck: not able to find' &
        & // ' interval in check table, stopping ...'                 )
    end if
    write(logUnit(1),*) 'Check interval is: ', check%iterationInterval


    select case(equation%eq_kind)
    case('maxwell', 'maxwelldivcorrection','maxwell_2d')
      ! No need for a tolrance, we set it to zero
      check%tolerance = 0.0_rk
    case default
      ! Read the check tolerance
      call aot_get_val(L = conf, thandle = check_table, &
        &              key = 'tolerance', &
        &              val = check%tolerance, &
        &              ErrCode = iError, &
        &              default = 1.0e-13_rk )
      write(logUnit(1),*) 'Check tolerance is: ', check%tolerance
    end select

    check%isActive = .true.

    call aot_table_close(L = conf, thandle = check_table)
    call tem_horizontalSpacer(fUnit=logUnit(1))

  end subroutine atl_init_physCheck


  !> Routine to check if the physical values of a state are physically
  !! meaningful or not.
  subroutine atl_check_val( minlevel, maxlevel, statedata_list, mesh_list,     &
                          & stat, equation, scheme_list, poly_proj_pos,        &
                          & poly_proj_list,check, iteration, time              )
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: minlevel
    integer, intent(in) :: maxlevel
    type(atl_statedata_type), intent(in) :: statedata_list(minlevel:maxlevel)
    type(atl_cube_elem_type), intent(in) :: mesh_list(minlevel:maxlevel)
    !> The current status bits from the treelm general configuration
    type(tem_status_type), intent(inout) :: stat
    type(atl_equations_type), intent(in) :: equation
    type(atl_scheme_type), intent(in) :: scheme_list(minlevel:maxlevel)
    integer, intent(in) :: poly_proj_pos(minlevel:maxlevel)
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    type(atl_physCheck_type), intent(in) :: check
    integer, intent(in) :: iteration
    type(atl_global_timestep_type), intent(inout) :: time
    ! ---------------------------------------------------------------------------
    !> max velocity to calculate the cfl in euler and linearEuler
    real(kind=rk) :: maxVel
    integer :: iLevel
    integer :: nLevelFluids
    logical :: isPhysical
    ! ---------------------------------------------------------------------------

    ! Check if the checking for this iteration has to be done!
    if (check%isActive .and. mod(iteration,check%iterationInterval).eq.0) then
      isPhysical = .true.
      ! Check which equation we have
      select case(trim(equation%eq_kind))
      case('euler_1d')
        do iLevel = minlevel, maxlevel
          nLevelFluids = mesh_list(iLevel)%descriptor     &
            &                             %elem           &
            &                             %nElems(eT_fluid)
          if (nLevelFluids > 0) then
            isPhysical = isPhysical                                     &
              & .and. atl_physCheck_euler1d(                            &
              &   statedata    = statedata_list(iLevel),                &
              &   euler1d      = equation%Euler,                        &
              &   scheme       = scheme_list(iLevel),                   &
              &   tolerance    = check%tolerance,                       &
              &   poly_proj    = poly_proj_list(poly_proj_pos(iLevel)), &
              &   nElems_fluid = nLevelFluids                           )

            if (time%control%fixed_dt > 0.0_rk) then

              !calculate the max velocity for the cfl
              maxVel = maxval(                                                &
                & abs(time%elementSteps(iLevel)%euler_1d%maxVel(:))           &
                &   + abs(time%elementSteps(iLevel)%euler_1d%speedOfSound(:)) )
              isPhysical = isPhysical                                        &
                & .and. atl_cflCheck_euler(                                  &
                &   fixed_dt = time%control%fixed_dt,                        &
                &   maxVel   = maxVel,                                       &
                &   length   = mesh_list(iLevel)%length,                     &
                &   nPoly    = scheme_list(iLevel)%modg_1d%maxPolyDegree + 1 )
            end if

          end if

        end do

      case('euler_2d', 'navier_stokes_2d')
        do iLevel = minlevel, maxlevel
          nLevelFluids = mesh_list(iLevel)%descriptor     &
            &                             %elem           &
            &                             %nElems(eT_fluid)
          if (nLevelFluids > 0) then
            isPhysical = isPhysical                                       &
              & .and. atl_physCheck_euler2d(                              &
              &   statedata    = statedata_list(iLevel),                  &
              &   euler2d      = equation%Euler,                          &
              &   scheme       = scheme_list(iLevel),                     &
              &   tolerance    = check%tolerance,                         &
              &   poly_proj    = poly_proj_list(poly_proj_pos(iLevel)),   &
              &   total = mesh_list(iLevel)%descriptor%total,             &
              &   nElems_fluid = nLevelFluids                             )

            if (time%control%fixed_dt > 0.0_rk) then
              !calculate the max velocity for the cfl
              maxVel = maxval(                                               &
                & abs(time%elementSteps(iLevel)%euler_2d%maxVel(:))          &
                &   + abs(time%elementSteps(iLevel)%euler_2d%speedOfSound(:)))
              if ( trim(equation%eq_kind) == 'euler_2d') then
                isPhysical = isPhysical                                        &
                  & .and. atl_cflCheck_euler(                                  &
                  &   fixed_dt = time%control%fixed_dt,                        &
                  &   maxVel   = maxVel,                                       &
                  &   length   = mesh_list(iLevel)%length,                     &
                  &   nPoly    = scheme_list(iLevel)%modg_2d%maxPolyDegree + 1 )

              else ! is NS
                isPhysical = isPhysical                                          &
                  & .and. atl_cflCheck_navier(                                   &
                  &   fixed_dt   = time%control%fixed_dt,                        &
                  &   maxVel     = maxVel,                                       &
                  &   length     = mesh_list(iLevel)%length,                     &
                  &   mu         = equation%NavierStokes%mu,                     &
                  &   therm_cond = equation%NavierStokes%therm_cond,             &
                  &   nPoly      = scheme_list(iLevel)%modg_2d%maxPolyDegree + 1 )
              end if !  euler or NS

            end if ! fixed_dt > 0

          end if

        end do

      case('euler', 'loclineuler', 'navier_stokes')
        do iLevel = minlevel, maxlevel
          nLevelFluids = mesh_list(iLevel)%descriptor     &
            &                             %elem           &
            &                             %nElems(eT_fluid)
          if (nLevelFluids > 0) then
            isPhysical = isPhysical                                       &
              & .and. atl_physCheck_euler(                                &
              &   statedata    = statedata_list(iLevel),                  &
              &   euler        = equation%Euler,                          &
              &   scheme       = scheme_list(iLevel),                     &
              &   tolerance    = check%tolerance,                         &
              &   poly_proj    = poly_proj_list(poly_proj_pos(iLevel)),   &
              &   nElems_fluid = nLevelFluids                             )

            if (time%control%fixed_dt > 0.0_rk) then
              !calculate the max velocity for the cfl
              maxVel = maxval(                                             &
                & abs(time%elementSteps(iLevel)%euler%maxVel(:))           &
                &   + abs(time%elementSteps(iLevel)%euler%speedOfSound(:)) )
             if ( trim(equation%eq_kind) == 'euler') then
              isPhysical = isPhysical                                       &
                  & .and. atl_cflCheck_euler(                               &
                  &   fixed_dt = time%control%fixed_dt,                     &
                  &   maxVel   = maxVel,                                    &
                  &   length   = mesh_list(iLevel)%length,                  &
                  &   nPoly    = scheme_list(iLevel)%modg%maxPolyDegree + 1 )
              else

                isPhysical = isPhysical                                        &
                  & .and. atl_cflCheck_navier(                                 &
                  &   fixed_dt   = time%control%fixed_dt,                      &
                  &   maxVel     = maxVel,                                     &
                  &   length     = mesh_list(iLevel)%length,                   &
                  &   mu         = equation%NavierStokes%mu,                   &
                  &   therm_cond = equation%NavierStokes%therm_cond,           &
                  &   nPoly      = scheme_list(iLevel)%modg%maxPolyDegree + 1  )
              end if !euler or NS

            end if !fixed_dt > 0

          end if

        end do

      case('filtered_navier_stokes')
        select case(trim(equation%FiltNavierStokes%model_type))
          case('rans')
            do iLevel = minlevel, maxlevel
              isPhysical = isPhysical                                     &
                & .and. atl_physCheck_Rans(                               &
                &   statedata    = statedata_list(iLevel),                &
                &   euler        = equation%Euler,                        &
                &   scheme       = scheme_list(iLevel),                   &
                &   tolerance    = check%tolerance,                       &
                &   poly_proj    = poly_proj_list(poly_proj_pos(iLevel)), &
                &   nElems_fluid = mesh_list(iLevel)%descriptor           &
                &                                   %elem                 &
                &                                   %nElems(eT_fluid)     )
            end do
          case default
            write(logUnit(1),*) 'Turbulence model not defined ... stopping'
            call tem_abort()
          end select
      case('filtered_navier_stokes_2d')
        select case(trim(equation%FiltNavierStokes%model_type))
          case('rans_2d')
            do iLevel = minlevel, maxlevel
              isPhysical = isPhysical                                 &
                & .and. atl_physCheck_Rans_2d(                        &
                &   statedata    = statedata_list(iLevel),            &
                &   euler        = equation%Euler,                    &
                &   scheme       = scheme_list(iLevel),               &
                &   tolerance    = check%tolerance,                   &
                &   poly_proj    = poly_proj_list                     &
                &                  (poly_proj_pos(iLevel)),           &
                &   nElems_fluid = mesh_list(iLevel)%descriptor       &
                &                                   %elem             &
                &                                   %nElems(eT_fluid) )
            end do
          case default
            write(logUnit(1),*) 'Turbulence model not defined ... stopping'
            call tem_abort()
          end select
      case('acoustic_2d')
        do iLevel = minlevel, maxlevel
          isPhysical = isPhysical                                     &
            & .and. atl_physCheck_acoustic_2d(                        &
            &   statedata    = statedata_list(iLevel),                &
            &   scheme       = scheme_list(iLevel),                   &
            &   poly_proj    = poly_proj_list(poly_proj_pos(iLevel)), &
            &   nElems_fluid = mesh_list(iLevel)%descriptor           &
            &                                   %elem                 &
            &                                   %nElems(eT_fluid)     )
        end do
      case('acoustic' )
        do iLevel = minlevel, maxlevel
          isPhysical = isPhysical                                     &
            & .and. atl_physCheck_acoustic(                           &
            &   statedata    = statedata_list(iLevel),                &
            &   scheme       = scheme_list(iLevel),                   &
            &   poly_proj    = poly_proj_list(poly_proj_pos(iLevel)), &
            &   nElems_fluid = mesh_list(iLevel)%descriptor           &
            &                                   %elem                 &
            &                                   %nElems(eT_fluid)     )
        end do
      case('lineareuler' )
        do iLevel = minlevel, maxlevel
          nLevelFluids = mesh_list(iLevel)%descriptor     &
            &                             %elem           &
            &                             %nElems(eT_fluid)
          if (nLevelFluids > 0) then
            isPhysical = isPhysical                                     &
              & .and. atl_physCheck_lineareuler(                        &
              &   statedata    = statedata_list(iLevel),                &
              &   scheme       = scheme_list(iLevel),                   &
              &   poly_proj    = poly_proj_list(poly_proj_pos(iLevel)), &
              &   nElems_fluid = nLevelFluids                           )

             if (time%control%fixed_dt > 0.0_rk) then
              !calculate the max velocity for the cfl
              maxVel = abs(sqrt(sum((equation%linearEuler%velocity_0)**2)) &
                &    + equation%linearEuler%SpeedOfSound)
              isPhysical = isPhysical                                     &
                & .and. atl_cflCheck_euler(                               &
                &   fixed_dt = time%control%fixed_dt,                     &
                &   maxVel   = maxVel,                                    &
                &   length   = mesh_list(iLevel)%length,                  &
                &   nPoly    = scheme_list(iLevel)%modg%maxPolyDegree + 1 )
            end if

          end if

        end do

      case('lineareuler_2d' )
        do iLevel = minlevel, maxlevel
          nLevelFluids = mesh_list(iLevel)%descriptor     &
            &                             %elem           &
            &                             %nElems(eT_fluid)
          if (nLevelFluids > 0) then
            isPhysical = isPhysical  &
              & .and. atl_physCheck_lineareuler_2d(                     &
              &   statedata    = statedata_list(iLevel),                &
              &   scheme       = scheme_list(iLevel),                   &
              &   poly_proj    = poly_proj_list(poly_proj_pos(iLevel)), &
              &   nElems_fluid = nLevelFluids                           )

            if (time%control%fixed_dt > 0.0_rk) then
              !calculate the max velocity for the cfl
              maxVel = abs(sqrt(sum((equation%linearEuler%velocity_0)**2))   &
                &    + equation%linearEuler%SpeedOfSound)
              isPhysical = isPhysical                                        &
                & .and. atl_cflCheck_euler(                                  &
                &   fixed_dt = time%control%fixed_dt,                        &
                &   maxVel   = maxVel,                                       &
                &   length   = mesh_list(iLevel)%length,                     &
                &   nPoly    = scheme_list(iLevel)%modg_2d%maxPolyDegree + 1 )
            end if
          end if
        end do

      case('maxwell','pec_maxwell','pec_maxwell_scatter', 'pec_maxwell_pml', &
        & 'maxwelldivcorrection', 'maxwell_2d','pec_maxwell_2d'              )
        do iLevel = minlevel, maxlevel
          isPhysical = isPhysical                                 &
            & .and. atl_physCheck_maxwell(                        &
            &   statedata    = statedata_list(iLevel),            &
            &   nDofs        = scheme_list(iLevel)%nDofs,         &
            &   nElems_fluid = mesh_list(iLevel)%descriptor       &
            &                                   %elem             &
            &                                   %nElems(eT_fluid) )
        end do

      case default
        write(logUnit(1),*) 'ERROR in atl_check_val: no physical check for ' &
          &                 // ' this equation type defined, stopping ...'
        call tem_abort()
      end select

      stat%bits(tem_stat_nonPhysical) = .not. isPhysical

    end if

  end subroutine atl_check_val

  !> Routine to check if the physicle values of the state are physically
  !! meaningful or not for the Euler and Linear Euler equation, checking
  !! the cfl for a fixed timestep
  function atl_cflCheck_euler(fixed_dt, nPoly, length, maxVel) &
    &                            result(isPhysical)
    !--------------------------------------------------------------------------
    !> Fixed timestep prescribed in ateles
    real(kind=rk), intent(in) :: fixed_dt
    !> Length of the element
    real(kind=rk), intent(in) :: length
    !> Order per spatial direction
    integer, intent(in) :: nPoly
    !> max Velocity for calculating the cfl
    real(kind=rk), intent(in) :: maxVel
    logical :: isPhysical
    !--------------------------------------------------------------------------
    real(kind=rk) :: cfl
   !---------------------------------------------------------------------------

    isPhysical = .true.

     cfl = (fixed_dt * maxVel * 2._rk * (nPoly**2)) / length

     ! Check wether cfl number is greater or equal to 1
     if (cfl <= 0.0_rk .or. cfl >= 1.0_rk) then
       write(logUnit(1),*)'WARNING! Please check cfl number'
       write(logUnit(1),*) 'cfl =', cfl
       isPhysical = .false.
    end if

  end function atl_cflCheck_euler

  !> Routine to check if the physicle values of the state are physically
  !! meaningful or not for the Navier-Stokes equation, checking the cfl
  !! for a fixed timestep
  function atl_cflCheck_navier(fixed_dt, nPoly, length, maxVel, mu, therm_cond) &
    &                             result(isPhysical)
    !--------------------------------------------------------------------------
    !> Fixed timestep prescribed in ateles
    real(kind=rk), intent(in) :: fixed_dt
    !> Length of the element
    real(kind=rk), intent(in) :: length
    !> Order per spatial direction
    integer, intent(in) :: nPoly
    !> max Velocity for calculating the cfl
    real(kind=rk), intent(in) :: maxVel
    !> vicosity variable
    real (kind=rk), intent(in) :: mu
    !> thermal conductivity
    real (kind=rk), intent(in) :: therm_cond
    logical :: isPhysical
    !--------------------------------------------------------------------------
    !> cfl calculated from the vicosity part
    real(kind=rk) :: cfl_visc
    !> cfl considering the convective part
    real(kind=rk) :: cfl_conv
   !---------------------------------------------------------------------------

    isPhysical = .true.

     cfl_conv = (fixed_dt * maxVel * 2._rk * (nPoly**2)) / length

     cfl_visc = (fixed_dt * (max(mu, therm_cond)) * real(nPoly**4,rk)) / (length**2)

     ! Check wether cfl number is greater or equal to 1
     if (cfl_conv <= 0.0_rk .or. cfl_conv >= 1.0_rk) then
       write(logUnit(1),*)'WARNING! Please check cfl_conv number '
       write(logUnit(1),*)'cfl_conv =', cfl_conv
       isPhysical = .false.
     elseif (cfl_visc <= 0.0_rk .or. cfl_visc >= 1.0_rk) then
       write(logUnit(1),*)'WARNING! Please check cfl_visc number'
       write(logUnit(1),*)'cfl_visc =', cfl_visc
       isPhysical = .false.
    end if
  end function atl_cflCheck_navier

  !> Routine to check if the physical values of a state are physically
  !! meaningful or not for the Euler equation.
  function atl_physCheck_euler1d( statedata, euler1d, scheme, poly_proj,      &
    &                             nElems_fluid, tolerance                   ) &
    &                             result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_euler_type), intent(in) :: euler1d
    type(atl_scheme_type), intent(in) :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in) :: nElems_fluid
    real(kind=rk), intent(in) :: tolerance
    logical :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:), pressure(:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints
    integer :: oversamp_dofs
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_1D%nquadpoints
    oversamp_dofs= poly_proj%body_1D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_1d_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,3) )
      allocate( pointVal(nQuadPoints,3) )
      allocate( pressure(nQuadPoints) )

      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! Transform to point values
        ! --> modal space
        call ply_convert2oversample(state       = statedata%state(iElem,:,:), &
          &                         poly_proj   = poly_proj,                  &
          &                         nDim        = 1,                          &
          &                         modalCoeffs = modalCoeffs                 )
        ! --> oversamp modal space
        call ply_poly_project_m2n(me = poly_proj,        &
          &                       dim = 1 ,              &
          &                       nVars = 3,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! --> oversamp nodal space

        !  Calculate the pressure
        do iPoint = 1, nquadpoints
          pressure(iPoint) = (euler1d%isen_coef - 1.0_rk) * &
          &            ( pointVal(iPoint,3) - (0.5_rk / pointVal(iPoint,1) ) &
          &                               * (pointVal(iPoint,2)**2) )
        end do

        ! Check for density, energy and pressure
        do iPoint = 1, nquadpoints
          if( any(pointVal(iPoint,[1,3]).le.tolerance)    &
            & .or. pressure(iPoint).le.tolerance          &
            & .or. any(tem_isNan(pointVal(iPoint,[1,3]))) &
            & .or. tem_isNan(pressure(iPoint))            ) then
            isPhysical = .false.
          end if
        end do

      end do
    case default
      call tem_abort( 'ERROR in atl_physCheck_euler1d: not able to check' &
        & // ' physical values for this scheme, stopping ...'             )
    end select

  end function atl_physCheck_euler1d

  !> Routine to check if the physical values of a state are physically
  !! meaningful or not for the Euler equation.
  function atl_physCheck_euler2d( statedata, euler2d, scheme, poly_proj,   &
    &                             nElems_fluid, tolerance, total)result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_euler_type), intent(in) :: euler2d
    type(atl_scheme_type), intent(in) :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in) :: nElems_fluid
    real(kind=rk), intent(in) :: tolerance
    integer(kind=long_k), intent(in) :: total(:)
    logical :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:), pressure(:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints
    integer :: oversamp_dofs
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_2d_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,4) )
      allocate( pointVal(nQuadPoints,4) )
      allocate( pressure(nQuadPoints) )



      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! --> modal space
        call ply_convert2oversample(state       = statedata%state(iElem,:,:), &
          &                         poly_proj   = poly_proj,                  &
          &                         nDim        = 2,                          &
          &                         modalCoeffs = modalCoeffs                 )
        ! --> oversamp modal space

        ! Transform to point values
        call ply_poly_project_m2n(me = poly_proj ,       &
          &                       dim = 2 ,              &
          &                       nVars = 4,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! --> oversample nodal space

        !  Calculate the pressure

        do iPoint = 1,nQuadPoints
          pressure(iPoint) = (euler2d%isen_coef - 1.0_rk) * &
          &            (pointVal(iPoint,4) - (0.5_rk / pointVal(iPoint,1) ) &
          &                               * sum(pointVal(iPoint,2:3)**2))
        end do


        ! Check for density, energy and pressure and cfl number

        do iPoint = 1, nQuadPoints
          if( pointVal(iPoint,1).le.tolerance ) then
            write(logunit(1),*) 'Density value for treeid ', total(iElem), &
              & ' at point ', iPoint, ' is below tolerance. '
            isPhysical = .false.
          end if
          if( pointVal(iPoint,4).le.tolerance ) then
            write(logunit(1),*) 'Energy value for treeid ',  total(iElem), &
              & ' at point ' , iPoint, ' is below tolerance. '
            isPhysical = .false.
          end if
          if( pressure(iPoint).le.tolerance ) then
            write(logunit(1),*) 'Pressure value for treeid ',  total(iElem), &
              & ' at point ' , iPoint, ' is below tolerance. '
            isPhysical = .false.
          end if
          if( any(tem_isNan(pointVal(iPoint,[1,4]))) ) then
            write(logunit(1),*) 'Density or energy for treeid ',  &
              & total(iElem), ' at point ' , iPoint, ' is NAN. '
            isPhysical = .false.
          end if
          if( tem_isNan(pressure(iPoint)) ) then
            write(logunit(1),*) 'Pressure for treeid ',  total(iElem), &
              & ' at point ', iPoint, ' is NAN. '
            isPhysical = .false.
          end if
        end do

      end do



    case default
      call tem_abort( 'ERROR in atl_physCheck_euler2d: not able to check' &
        & // ' physical values for this scheme, stopping ...'             )
    end select


  end function atl_physCheck_euler2d


  !> Routine to check if the physical values of a state are physically
  !! meaningful or not for the Euler equation.
  function atl_physCheck_euler( statedata, euler, scheme, poly_proj,     &
    &                           nElems_fluid, tolerance)result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_euler_type), intent(in) :: euler
    type(atl_scheme_type), intent(in) :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in) :: nElems_fluid
    real(kind=rk), intent(in) :: tolerance
    logical :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:), pressure(:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints
    integer :: oversamp_dofs
    integer :: mpd1, mpd1_square
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_3D%nquadpoints
    oversamp_dofs = poly_proj%body_3D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,5) )
      allocate( pointVal(nquadpoints,5) )
      allocate( pressure(nquadpoints) )

      ! used in oversampling loop
      mpd1 = poly_proj%min_degree+1
      mpd1_square = mpd1**2

      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! Transform to point values
        ! --> modal space
        call ply_convert2oversample(state       = statedata%state(iElem,:,:), &
          &                         poly_proj   = poly_proj,                  &
          &                         nDim        = 3,                          &
          &                         modalCoeffs = modalCoeffs                 )
        ! --> oversamp modal space
        call ply_poly_project_m2n(me = poly_proj ,       &
          &                       dim = 3 ,              &
          &                       nVars = 5,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! -->oversamp nodal space

        !  Calculate the pressure

        do iPoint = 1, nquadpoints
          pressure(iPoint) = (euler%isen_coef - 1.0_rk) * &
          &            (pointVal(iPoint,5) - (0.5_rk / pointVal(iPoint,1) ) &
          &                               * sum(pointVal(iPoint,2:4)**2))
        end do


        ! Check for density, energy and pressure

        do iPoint = 1, nquadpoints
          if( any(pointVal(iPoint,[1,5]).le.tolerance)    &
            & .or. pressure(iPoint).le.tolerance          &
            & .or. any(tem_isNan(pointVal(iPoint,[1,5]))) &
            & .or. tem_isNan(pressure(iPoint))            ) then
            isPhysical = .false.
          end if
        end do



      end do



    case default
      call tem_abort( 'ERROR in atl_physCheck_euler: not able to check' &
        & // ' physical values for this scheme, stopping ...'           )
    end select


  end function atl_physCheck_euler


  !> summary: Routine to check if the physical values of a state are
  !! physically meaningful or not for the Filtered Navier Stokes equation.
  function atl_physCheck_Rans( statedata, euler, scheme, poly_proj, &
    &                     nElems_fluid, tolerance )result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_euler_type), intent(in) :: euler
    type(atl_scheme_type), intent(in) :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in) :: nElems_fluid
    real(kind=rk), intent(in) :: tolerance
    logical :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:), pressure(:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints
    integer :: oversamp_dofs
    integer :: mpd1, mpd1_square
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_3D%nquadpoints
    oversamp_dofs = poly_proj%body_3D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,7) )
      allocate( pointVal(nquadpoints,7) )
      allocate( pressure(nquadpoints) )

      ! used in oversampling loop
      mpd1 = poly_proj%min_degree+1
      mpd1_square = mpd1**2

      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! Transform to point values
        ! --> modal space
        call ply_convert2oversample(state       = statedata%state(iElem,:,:), &
          &                         poly_proj   = poly_proj,                  &
          &                         nDim        = 3,                          &
          &                         modalCoeffs = modalCoeffs                 )
        ! --> oversamp modal space
        call ply_poly_project_m2n(me = poly_proj ,       &
          &                       dim = 3 ,              &
          &                       nVars = 7,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! -->oversamp nodal space

        !  Calculate the pressure
        do iPoint = 1, nquadpoints
          pressure(iPoint) = (euler%isen_coef - 1.0_rk) * &
          &            (pointVal(iPoint,5) - (0.5_rk / pointVal(iPoint,1) ) &
          &                               * sum(pointVal(iPoint,2:4)**2)    &
          &            - pointVal(iPoint,6))
        end do

        ! Check for density, energy and pressure
        do iPoint = 1, nquadpoints
          if( any(pointVal(iPoint,[1,7]).le.tolerance)    &
            & .or. pressure(iPoint).le.tolerance          &
            & .or. any(tem_isNan(pointVal(iPoint,[1,7]))) &
            & .or. tem_isNan(pressure(iPoint))            ) then
            isPhysical = .false.
          end if
        end do


      end do

    case default
      call tem_abort( 'ERROR in atl_physCheck_Rans: not able to check' &
        & // ' physical values for this scheme, stopping ...'          )
    end select


  end function atl_physCheck_Rans


  function atl_physCheck_Rans_2d( statedata, euler, scheme, poly_proj, &
    &                     nElems_fluid, tolerance )result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_euler_type), intent(in) :: euler
    type(atl_scheme_type), intent(in) :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in) :: nElems_fluid
    real(kind=rk), intent(in) :: tolerance
    logical :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:), pressure(:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: nquadpoints
    integer :: oversamp_dofs, iElem, iPoint
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_2d_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,6) )
      allocate( pointVal(nquadpoints,6) )
      allocate( pressure(nquadpoints) )



      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! Transform to point values
        ! --> modal space
        call ply_convert2oversample(state       = statedata%state(iElem,:,:), &
          &                         poly_proj   = poly_proj,                  &
          &                         nDim        = 2,                          &
          &                         modalCoeffs = modalCoeffs                 )
        ! --> oversamp modal space
        call ply_poly_project_m2n(me = poly_proj ,       &
          &                       dim = 2 ,              &
          &                       nVars = 6,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! -->oversamp nodal space

        !  Calculate the pressure

        do iPoint = 1, nquadpoints
          pressure(iPoint) = (euler%isen_coef - 1.0_rk) * &
          &            (pointVal(iPoint,4) - (0.5_rk / pointVal(iPoint,1) ) &
          &                               * sum(pointVal(iPoint,2:3)**2)    &
          &            - pointVal(iPoint,5))
        end do


        ! Check for density, energy and pressure

        do iPoint = 1, nquadpoints
          if( any(pointVal(iPoint,[1,4]).le.tolerance) .or. pressure(iPoint).le.tolerance &
            & .or. any(tem_isNan(pointVal(iPoint,[1,4]))) .or. tem_isNan(pressure(iPoint)) ) then
            isPhysical = .false.
          end if
        end do



      end do



    case default
      write(logUnit(1),*) 'ERROR in atl_physCheck_Rans: not able to check physical values' // &
                    & ' for this scheme, stopping ...'
      call tem_abort()
    end select


  end function atl_physCheck_Rans_2d


  !> Routine to check if the physical values of a state are physically
  !! meaningful or not for the acoustic 2d equation.
  function atl_physCheck_acoustic_2d( statedata, scheme, poly_proj, &
    &                                 nElems_fluid                ) &
    &                                 result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_scheme_type), intent(in)          :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in)                        :: nElems_fluid
    logical                                    :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints
    integer :: oversamp_dofs
    integer :: mpd1, mpd1_square
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_2d_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,3) )
      allocate( pointVal(nquadpoints,3) )

      ! used in oversampling loop
      mpd1 = poly_proj%min_degree+1
      mpd1_square = mpd1**2

      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! --> modal space
        call ply_convert2oversample(state       = statedata%state(iElem,:,:), &
          &                         poly_proj   = poly_proj,                  &
          &                         nDim        = 2,                          &
          &                         modalCoeffs = modalCoeffs                 )
        ! --> oversamp modal space

        ! Transform to point values
        call ply_poly_project_m2n(me = poly_proj,        &
          &                       dim = 2 ,              &
          &                       nVars = 3,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! -->oversamp nodal space

        ! Check for density is negative or any NAN value

        do iPoint = 1, nquadpoints
          if( (any(tem_isNan(pointVal(iPoint,[1,2,3]))) ) ) then
            isPhysical = .false.
          end if
        end do

      end do

    case default
      call tem_abort( 'ERROR in atl_physCheck_acoustic_2d: not able to check' &
        & // 'physical values for this scheme, stopping ...'                  )
    end select

  end function atl_physCheck_acoustic_2d

  !> Routine to check if the physical values of a state are physically
  !! meaningful or not for the acoustic equation.
  function atl_physCheck_acoustic( statedata, scheme, poly_proj, &
    &                              nElems_fluid                ) &
    &                              result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_scheme_type), intent(in)          :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in)                        :: nElems_fluid
    logical                                    :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints
    integer :: oversamp_dofs
    integer :: mpd1, mpd1_square
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_3D%nquadpoints
    oversamp_dofs = poly_proj%body_3D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,4) )
      allocate( pointVal(nquadpoints,4) )

      ! used in oversampling loop
      mpd1 = poly_proj%min_degree+1
      mpd1_square = mpd1**2

      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! --> modal space
        call ply_convert2oversample(state       = statedata%state(iElem,:,:), &
          &                         poly_proj   = poly_proj,                  &
          &                         nDim        = 2,                          &
          &                         modalCoeffs = modalCoeffs                 )
        ! --> oversamp modal space

        ! Transform to point values
        call ply_poly_project_m2n(me = poly_proj,        &
          &                       dim = 3 ,              &
          &                       nVars = 4,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! -->oversamp nodal space

        ! Check for density is negative or any NAN value
        do iPoint = 1, nquadpoints
          if( (any(tem_isNan(pointVal(iPoint,[1,2,3,4]))) ) ) then
            isPhysical = .false.
          end if
        end do


      end do

    case default
      call tem_abort( 'ERROR in atl_physCheck_acoustic: not able to check' &
        & // 'physical values for this scheme, stopping ...'               )
    end select

  end function atl_physCheck_acoustic

  !> Routine to check if the physical values of a state are physically
  !! meaningful or not for the linear euler equation.
  function atl_physCheck_lineareuler( statedata, scheme, poly_proj, &
    &                                 nElems_fluid                ) &
    &                                 result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_scheme_type), intent(in)          :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in)                        :: nElems_fluid
    logical                                    :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints
    integer :: oversamp_dofs
    integer :: mpd1, mpd1_square
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_3D%nquadpoints
    oversamp_dofs = poly_proj%body_3D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,5) )
      allocate( pointVal(nquadpoints,5) )

      ! used in oversampling loop
      mpd1 = poly_proj%min_degree+1
      mpd1_square = mpd1**2

      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! --> modal space
        call ply_convert2oversample( state       = statedata%state(iElem,:,:), &
          &                          poly_proj   = poly_proj,                  &
          &                          nDim        = 3,                          &
          &                          modalCoeffs = modalCoeffs                 )
        ! --> oversamp modal space

        ! Transform to point values
        call ply_poly_project_m2n(me = poly_proj,        &
          &                       dim = 3 ,              &
          &                       nVars = 5,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! -->oversamp nodal space

        ! Check for density is negative or any NAN value

        do iPoint = 1, nquadpoints
          if( any(tem_isNan(pointVal(iPoint,[1,5])))) then
            isPhysical = .false.
          end if
        end do


      end do


    case default
      call tem_abort( 'ERROR in atl_physCheck_linearEuler: not able to check' &
        & // 'physical values for this scheme, stopping ...'               )
    end select



  end function atl_physCheck_lineareuler


  !> summary: Routine to check if the physical values of a state are physically meaningful
  !! or not for the linear euler equation.
  function atl_physCheck_lineareuler_2d( statedata, scheme, poly_proj, &
    &                                    nElems_fluid                ) &
    &                                    result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    type(atl_scheme_type), intent(in)          :: scheme
    type(ply_poly_project_type), intent(inout) :: poly_proj
    integer, intent(in)                        :: nElems_fluid
    logical                                    :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:)
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints
    integer :: oversamp_dofs
    integer :: mpd1, mpd1_square
    ! ---------------------------------------------------------------------------

    isPhysical = .true.

    nquadpoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    select case(scheme%scheme)
    case(atl_modg_2d_scheme_prp)

      allocate( modalCoeffs(oversamp_dofs,4) )
      allocate( pointVal(nquadpoints,4) )

      ! used in oversampling loop
      mpd1 = poly_proj%min_degree+1
      mpd1_square = mpd1**2

      ! For each element we transform to point values and check
      ! whether density, pressure and energy are still positive.
      do iElem = 1, nElems_fluid

        ! --> modal space
        call ply_convert2oversample(state= statedata%state(iElem,:,:), &
          &                         ndim = 2,                          &
          &                         poly_proj= poly_proj,              &
          &                         modalCoeffs = modalCoeffs          )
        ! --> oversamp modal space

        ! Transform to point values
        call ply_poly_project_m2n(me = poly_proj,        &
          &                       dim = 2 ,              &
          &                       nVars = 4,             &
          &                       nodal_data=pointVal,   &
          &                       modal_data=modalCoeffs )
        ! -->oversamp nodal space

        ! Check for density is negative or any NAN value

        do iPoint = 1, nquadpoints
          if( any(tem_isNan(pointVal(iPoint,[1,4]))))  then
            isPhysical = .false.
          end if
        end do


      end do


    case default
      write(logUnit(1),*)'ERROR in atl_physCheck_linearEuler_2d: not able to check physical values' // &
                    & ' for this scheme, stopping ...'
      call tem_abort()
    end select


  end function atl_physCheck_lineareuler_2d

  !> Routine to check if the physical values of a state are physically
  !! meaningful or not for the Maxwell equation.
  function atl_physCheck_maxwell( statedata, nElems_fluid, nDofs) &
    &                             result(isPhysical)
    ! ---------------------------------------------------------------------------
    type(atl_statedata_type), intent(in) :: statedata
    integer, intent(in) :: nElems_fluid
    integer, intent(in) :: nDofs
    logical :: isPhysical
    ! ---------------------------------------------------------------------------
    ! Loop vars
    integer :: iDof
    ! ---------------------------------------------------------------------------

    isPhysical = .true.


    ! For each element we transform to point values and check
    ! whether density, pressure and energy are still positive.


    ! Check for density, energy and pressure

    do iDof = 1, nDofs
      if( any(tem_isNan( statedata%state(:nElems_fluid,iDof,:) )) ) then
        isPhysical = .false.
      end if
    end do


  end function atl_physCheck_maxwell

end module atl_physCheck_module
