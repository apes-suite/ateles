! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
!! Routines, functions and datatypes to steer the
!! time stepping mechanism in ATELES.
module atl_global_time_integration_module
  use env_module,               only: rk
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit
  use tem_element_module,       only: eT_fluid

  use aotus_module,             only: flu_State, aot_get_val
  use aot_table_module,         only: aot_table_open, aot_table_close

  use atl_time_integration_module, only:       &
    &   atl_global_timestep_type,              &
    &   atl_timestep_control_type,             &
    &   explicitEuler,                         &
    &   explicitRungeKutta,                    &
    &   explicitSSPRungeKutta,                 &
    &   explicitLocalPredictorGlobalCorrector, &
    &   explicitRungeKuttaTaylor,              &
    &   imexRungeKutta
  use atl_kerneldata_module,    only: atl_statedata_type, &
    &                                 atl_kerneldata_type
  use atl_fwdEuler_module,      only: atl_init_explicitEuler
  use atl_rk4_module,           only: atl_init_explicitRungeKutta
  use atl_rktaylor_module,      only: atl_init_explicitRungeKuttaTaylor
  use atl_imexrk_module,        only: atl_init_imexRungeKutta
  use atl_predcor_cerk4_module, only: &
    &   atl_init_explicitLocalPredictorGlobalCorrector
  use atl_ssprk2_module,        only: atl_init_explicit_ssprk
  use atl_equation_module,      only: atl_equations_type
  use atl_cube_elem_module,     only: atl_cube_elem_type
  use atl_elemental_time_integration_module, only: atl_timestep_type
  use atl_eqn_euler_derive_module,           only: atl_eqn_euler_cons2prim
  use atl_eqn_euler_2d_derive_module,        only: atl_eqn_euler_2d_cons2prim
  use atl_eqn_euler_1d_derive_module,        only: atl_eqn_euler_1d_cons2prim
  use ply_poly_project_module,   only: ply_poly_project_type, ply_poly_project_m2n
  use ply_oversample_module,     only: ply_convert2oversample

  implicit none
  private

  public :: atl_global_time_integration_load
  public :: atl_init_global_time_integration, atl_initTimeStepInfo


contains


  ! ****************************************************************************
  !> Routine to load the timestepping scheme from the config.
  subroutine atl_global_time_integration_load( me, conf, equation )
    ! --------------------------------------------------------------------------
    !> the scheme you want to initialize.
    type(atl_global_timestep_type), intent(inout) :: me
    !> flu binding to lua configuration file.
    type(flu_State), intent(in) :: conf
    !> The equaton you are using.
    type(atl_equations_type), intent(in) :: equation
    ! --------------------------------------------------------------------------
    character(len=100) :: timestepName, controlName
    integer :: scheme_table, temporal_table, control_table
    integer :: iError, steps
    ! --------------------------------------------------------------------------

    call aot_table_open(L=conf, thandle=scheme_table, key='scheme')
    if (scheme_table == 0) then
      write(logUnit(1),*) 'ERROR in init_time_integarion: no scheme table in ' &
        &                 // 'Lua configuration file found,stopping...'
      call tem_abort()
    end if

    ! open the temporal subtable
    call aot_table_open(L=conf, parent=scheme_table, thandle=temporal_table, &
      &                 key='temporal')

    ! get the name of timestepping scheme
    call aot_get_val(L = conf, thandle = temporal_table, &
      &              key = 'name', &
      &              val = timestepName, &
      &              ErrCode = iError)

    select case(timestepName)
    case('explicitEuler')
      me%timestepType = explicitEuler
      me%nSteps = 1

    case('imexRungeKutta')
      ! get the number of steps of the runge kutta method
      call aot_get_val( L       = conf,           &
        &               thandle = temporal_table, &
        &               key     = 'steps',        &
        &               val     = steps,          &
        &               ErrCode = iError          )
      me%nSteps = steps
      me%timestepType = imexRungeKutta

    case('explicitRungeKutta')
      ! get the number of steps of the runge kutta method
      call aot_get_val(L = conf, thandle = temporal_table, &
      &                key = 'steps', &
      &                val = steps, &
      &                ErrCode = iError)
      me%nSteps = steps
      me%timestepType = explicitRungeKutta

    case('explicitRungeKuttaTaylor')
      ! get the number of steps of the runge kutta method
      call aot_get_val(L = conf, thandle = temporal_table, &
      &                key = 'steps', &
      &                val = steps, &
      &                ErrCode = iError)
      me%nSteps = steps
      me%timestepType = explicitRungeKuttaTaylor

    case('explicitSSPRungeKutta')
      ! get the number of steps of the runge kutta method
      call aot_get_val(L = conf, thandle = temporal_table, &
      &                key = 'steps', &
      &                val = steps, &
      &                ErrCode = iError)
      me%nSteps = steps
      me%timestepType = explicitSSPRungeKutta

    case('explicitLocalPredictorGlobalCorrector')
      ! get the number of steps of the predictor runge kutta method
      call aot_get_val(L = conf, thandle = temporal_table, &
      &                key = 'steps', &
      &                val = steps, &
      &                ErrCode = iError)
      me%nSteps = steps
      me%timestepType = explicitLocalPredictorGlobalCorrector

    case default
      write(logUnit(1),*) 'Unknown timestepping algorithm!'
      write(logUnit(1),*) 'Please use one of the following timestepping' &
        &                 //' algorithms:'
      write(logUnit(1),*) ' - explicitEuler '
      write(logUnit(1),*) ' - imexRungeKutta '
      write(logUnit(1),*) ' - explicitRungeKutta '
      write(logUnit(1),*) ' - explicitRungeKuttaTaylor '
      write(logUnit(1),*) ' - explicitSSPRungeKutta '
      write(logUnit(1),*) ' - explicitLocalPredictorGlobalCorrector '
      write(logUnit(1),*) ', stopping ...'
      call tem_abort()
    end select

    write(logunit(1),*) 'Time integration scheme:', timestepName
    write(logunit(1),*) '    with steps=', me%nSteps

    ! open the control subtable
    call aot_table_open(L=conf, parent=temporal_table, thandle=control_table, &
      &                 key='control')

    ! Set fixed_dt to a negative value to deactivate it.
    me%control%fixed_dt = -1.0_rk

    controlName = 'UNKNOWN'

    ! get control mechanism of the timestep
    call aot_get_val(L = conf, thandle = control_table, &
      &              key = 'name', &
      &              val = controlName, &
      &              ErrCode = iError)
    if (iError.ne.0) then
      write(logUnit(1),*) 'ERROR in atl_init_global_time_integration:'
      write(logUnit(1),*) 'Could not get a name for the timestep control ...'
      call tem_abort()
    end if

    select case(trim(controlName))
    case('cfl')
      ! get the cfl number (for the convective part of the equation)
      call aot_get_val(L = conf, thandle = control_table, &
        &              key = 'cfl', &
        &              val = me%control%cfl, &
        &              ErrCode = iError)
      if (iError.ne.0) then
        write(logUnit(1),*) 'ERROR in atl_init_global_time_integration:'
        write(logUnit(1),*) 'No CFL number specified, stopping ...'
        call tem_abort()
      end if
      write(logunit(1),*) '    dynamically computed timestep,' &
        &                 // ' with CFL number=', me%control%cfl

      call aot_get_val( L       = conf,                          &
        &               thandle = control_table,                 &
        &               key     = 'use_modal_estimate',          &
        &               val     = me%control%use_modal_estimate, &
        &               default = .false.,                       &
        &               ErrCode = iError                         )

      if (equation%nDerivatives > 0) then
        ! get the cfl number (for the viscous part of the equation)
        call aot_get_val(L = conf, thandle = control_table, &
          &              key = 'cfl_visc', &
          &              val = me%control%cfl_visc, &
          &              ErrCode = iError )
        if (iError.ne.0) then
          write(logUnit(1),*) 'ERROR in atl_global_time_integration_load:'
          write(logUnit(1),*) 'No cfl_visc number specified, stopping ...'
          call tem_abort()
        end if
        write(logunit(1),*) '                 and for viscous terms' &
          &                 //' CFL number =',                       &
          &                 me%control%cfl_visc
      end if
      if (me%control%use_modal_estimate) then
        write(logunit(1),*) '    using a Modal estimate of the state for' &
          &                 //' CFL condition'
      else
        write(logunit(1),*) '    using a Nodal values to limit timestep'
      end if

    case ('fixed')
      call aot_get_val( L       = conf,                &
        &               thandle = control_table,       &
        &               key     = 'dt',                &
        &               val     = me%control%fixed_dt, &
        &               ErrCode = iError               )
        if ( (iError /= 0) .or. (me%control%fixed_dt <= 0.0_rk) ) then
          write(logUnit(1),*) 'ERROR in atl_global_time_integration_load:'
          write(logUnit(1),*) 'Something wrong with the fixed timestep!'
          write(logUnit(1),*) 'Please provide a dt > 0.'
          write(logUnit(1),*) 'Stopping ...'
          call tem_abort()
        end if
      write(logunit(1),*) '    fixed timestep of dt=', me%control%fixed_dt
    case default
      write(logUnit(1),*) 'Unknwon timestep control mechanism'
      write(logUnit(1),*) 'Should be one of:'
      write(logUnit(1),*) '* cfl'
      write(logUnit(1),*) '* fixed'
      write(logUnit(1),*) 'But is: ' // trim(controlName)
      write(logUnit(1),*) ' stopping...'
      call tem_abort()
    end select

    call aot_table_close(L = conf, thandle = control_table)
    call aot_table_close(L = conf, thandle = temporal_table)
    call aot_table_close(L = conf, thandle = scheme_table)


  end subroutine atl_global_time_integration_load
  ! ****************************************************************************


  ! ****************************************************************************
  !> Routine to init the timestepping scheme.
  !!
  !! This allocates datafields as needed. Be aware that the configuration
  !! has to be loaded beforehand!
  subroutine atl_init_global_time_integration( me, minLevel, maxLevel,    &
    &                                          statedata_list, mesh_list, &
    &                                          equation                   )
    ! --------------------------------------------------------------------------
    !> the scheme you want to initialize.
    type(atl_global_timestep_type), intent(inout) :: me
    !> The minimum level of the mesh
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The state list used in your solver, for each level one entry.
    type(atl_statedata_type) , intent(in) :: statedata_list(minLevel:maxLevel)
    !> Mesh description
    type(atl_cube_elem_type), intent(in) :: mesh_list(minLevel:maxLevel)
    !> The equaton you are using.
    type(atl_equations_type), intent(in) :: equation
    ! --------------------------------------------------------------------------
    integer :: iLevel, n_fluids
    ! --------------------------------------------------------------------------

    select case(me%timestepType)
    case(explicitEuler)
      call atl_init_explicitEuler(me,minLevel, maxLevel)

    case(imexRungeKutta)
      call atl_init_imexRungeKutta( me, minLevel, maxLevel, me%nsteps, &
        &                           statedata_list                     )

    case(explicitRungeKutta)
      call atl_init_explicitRungeKutta(me, minLevel, maxLevel, me%nsteps, &
        &                              statedata_list                     )

    case(explicitRungeKuttaTaylor)
      call atl_init_explicitRungeKuttaTaylor(me, minLevel, maxLevel,   &
        &                                    me%nsteps, statedata_list )

    case(explicitSSPRungeKutta)
      call atl_init_explicit_ssprk(me, minLevel, maxLevel, me%nsteps, &
        &                          statedata_list                     )

    case(explicitLocalPredictorGlobalCorrector)
      call atl_init_explicitLocalPredictorGlobalCorrector(me,     &
        &                            minLevel, maxLevel,          &
        &                            me%nsteps, statedata_list    )

    end select

    select case(trim(equation%eq_kind))
    case( 'maxwell', 'maxwelldivcorrection','maxwell_2d','advection_1d',     &
      & 'acoustic', 'acoustic_2d','heat_1d','heat_2d','heat', 'lineareuler', &
      & 'lineareuler_2d'                                                     )
      ! nothing to be done here: the timestep depends only on material
      ! parameters
    case('euler', 'navier_stokes', 'filtered_navier_stokes')
      do iLevel = minlevel, maxlevel
        ! get the number of fluid elements in the mesh
        n_fluids = mesh_list(iLevel)%descriptor%elem%nElems(eT_fluid)
        allocate(me%elementSteps(iLevel)%euler%maxVel(n_fluids))
        allocate(me%elementSteps(iLevel)%euler%speedOfSound(n_fluids))
      end do
    case('loclineuler','loclineuler_1d')
      do iLevel = minlevel, maxlevel
        ! get the number of fluid elements in the mesh
        n_fluids = mesh_list(iLevel)%descriptor%elem%nElems(eT_fluid)
        allocate(me%elementSteps(iLevel)%LoclinEuler%meanVel(n_fluids))
        allocate(me%elementSteps(iLevel)%LoclinEuler%speedOfSound(n_fluids))
      end do
    case('euler_2d','navier_stokes_2d', 'filtered_navier_stokes_2d')
      do iLevel = minlevel, maxlevel
        ! get the number of fluid elements in the mesh
        n_fluids = mesh_list(iLevel)%descriptor%elem%nElems(eT_fluid)
        allocate(me%elementSteps(iLevel)%euler_2d%maxVel(n_fluids))
        allocate(me%elementSteps(iLevel)%euler_2d%speedOfSound(n_fluids))
      end do
    case('euler_1d')
      do iLevel = minlevel, maxlevel
        ! get the number of fluid elements in the mesh
        n_fluids = mesh_list(iLevel)%descriptor%elem%nElems(eT_fluid)
        allocate(me%elementSteps(iLevel)%euler_1d%maxVel(n_fluids))
        allocate(me%elementSteps(iLevel)%euler_1d%speedOfSound(n_fluids))
      end do
    case default
      write(logUnit(1),*) 'ERROR in atl_init_global_time_integration:'
      write(logUnit(1),*) 'Unknown equation ', trim(equation%eq_kind), &
        &                 ', stopping ...'
      call tem_abort()
    end select

  end subroutine atl_init_global_time_integration
  ! ****************************************************************************


  ! ****************************************************************************
  !> Subroutine to initialize the timestep information for the first iteration
  subroutine atl_initTimeStepInfo( equation, mesh, poly_proj, timestep, &
    &                              control, statedata, kerneldata )
    ! --------------------------------------------------------------------------
    !> The equaton you are using.
    type(atl_equations_type), intent(in) :: equation
    !> Mesh description
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The projection methos of the simulation.
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Timestep information. This routine will update information about the
    !! the maximum velocity and the maximum speed of sound in the system.
    type(atl_timestep_type), intent(inout) :: timestep
    !> Description of the timestep control.
    type(atl_timestep_control_type), intent(in) :: control
    !> The current state for all elements on the current level
    type(atl_statedata_type), intent(in) :: statedata
    !> Additional information on the state.
    type(atl_kerneldata_type), intent(in) :: kerneldata
    ! --------------------------------------------------------------------------
    integer :: iElem, n_fluids
    real(kind=rk), allocatable :: pointVal(:,:), velAbs(:), modalCoeffs(:,:)
    real(kind=rk) :: rho_min, rho_max
    real(kind=rk) :: m_max(3), m_min(3)
    real(kind=rk) :: m_limited
    real(kind=rk) :: e_max, gamprod
    integer :: nquadpoints
    integer :: oversamp_dofs
    ! --------------------------------------------------------------------------

    select case(trim(equation%eq_kind))
    case( 'maxwell', 'maxwelldivcorrection','maxwell_2d','advection_1d',        &
      & 'acoustic', 'acoustic_2d', 'heat_1d', 'heat_2d', 'heat', 'lineareuler', &
      & 'lineareuler_2d'                                                        )
      ! nothing to be done here: the timestep depends only on material
      ! parameters
    case('euler','navier_stokes', 'filtered_navier_stokes')

      ! get the number of fluid elements in the mesh
      n_fluids = mesh%descriptor%elem%nElems(eT_fluid)
      gamprod = (equation%euler%isen_coef-1)*equation%euler%isen_coef

      m_or_n: if (control%use_modal_estimate) then
        ! Use modal estimation for timestep limitation.
        elemModal: do iElem = 1, n_fluids
          rho_min = max( statedata%state(iElem,1,1)                  &
            &             - kerneldata%deviation(iElem,1),           &
            &            statedata%state(iElem,1,1)*epsilon(rho_min) )
          m_max(1) = statedata%state(iElem,1,2) &
            &        + kerneldata%deviation(iElem,2)
          m_max(2) = statedata%state(iElem,1,3) &
            &        + kerneldata%deviation(iElem,3)
          m_max(3) = statedata%state(iElem,1,4) &
            &        + kerneldata%deviation(iElem,4)
          m_min(1) = max( statedata%state(iElem,1,2)       &
            &             - kerneldata%deviation(iElem,2), &
            &             0.0_rk                           )
          m_min(2) = max( statedata%state(iElem,1,3)       &
            &             - kerneldata%deviation(iElem,3), &
            &             0.0_rk                           )
          m_min(3) = max( statedata%state(iElem,1,4)       &
            &             - kerneldata%deviation(iElem,4), &
            &             0.0_rk                           )
          rho_max = statedata%state(iElem,1,1) &
            &       + kerneldata%deviation(iElem,1)
          e_max = statedata%state(iElem,1,5) &
            &     + kerneldata%deviation(iElem,5)
          timestep%euler%maxVel(iElem) = sqrt(m_max(1)**2    &
            &                                 + m_max(2)**2  &
            &                                 + m_max(3)**2) &
            &                            / rho_min
          m_limited = min(e_max, m_min(1)**2 + m_min(2)**2 + m_min(3)**2)
          timestep%euler%speedOfSound(iElem)         &
            &  = sqrt( gamprod*(e_max/rho_min        &
            &                   - m_limited/rho_max) )
        end do elemModal
      else m_or_n
        ! Need to compute timestep limitation with nodal values.

        ! detect correct amount of quadpoints and modes due to projection
        ! method
        nquadpoints = poly_proj%body_3D%nquadpoints
        oversamp_dofs= poly_proj%body_3D%oversamp_dofs

        allocate(modalCoeffs(oversamp_dofs,equation%varSys%nScalars))
        allocate(pointVal(nquadpoints,equation%varSys%nScalars))
        allocate(velAbs(nquadpoints))
        elemNodal: do iElem = 1, n_fluids

          ! get the modal coefficients of the current element.
          ! --> modal space
          if (equation%euler%ensure_positivity) then
            call ply_convert2oversample(                           &
              &    state             = statedata%state(iElem,:,:), &
              &    poly_proj         = poly_proj,                  &
              &    nDim              = 3,                          &
              &    modalCoeffs       = modalCoeffs,                &
              &    ensure_positivity = [.true.,                    &
              &                         .false., .false., .false., &
              &                         .true. ]                   )
          else
            call ply_convert2oversample(                     &
              &    state       = statedata%state(iElem,:,:), &
              &    poly_proj   = poly_proj,                  &
              &    nDim        = 3,                          &
              &    modalCoeffs = modalCoeffs                 )
          end if
          ! --> oversamp modal space

          ! transform to point values
          !do iVar=1,equation%varSys%nScalars
          call ply_poly_project_m2n(me = poly_proj,                   &
            &                       dim = 3 ,                         &
            &                       nVars = equation%varSys%nScalars, &
            &                       nodal_data= pointVal,             &
            &                       modal_data= modalCoeffs           )
          ! --> oversampling nodal space

          ! transform to primitive variables (in place)
          call atl_eqn_euler_cons2prim(equation, pointVal)

          ! calculate the max velocity and the speed of sound
          ! ... build the absolute velocity
          velAbs(:) = sqrt( (pointVal(:,2))**2 + (pointVal(:,3))**2 &
                              & + (pointVal(:,4))**2 )
          ! assign the values to the timestepping
          timestep%euler%maxVel(iElem) = maxval( velAbs(:)  )
          timestep%euler%speedOfSound(iElem) =                               &
            & sqrt(                                                          &
            &   maxval(                                                      &
            &     equation%euler%isen_coef * pointVal(:,5) / pointVal(:,1) ) )
        end do elemNodal

        deallocate(modalCoeffs)
        deallocate(pointVal)
        deallocate(velAbs)

      end if m_or_n

    case('loclineuler')
      ! subroutine that computes the absolute magnitude of mean velocity
      ! and speed of sound for every element

      ! get the number of fluid elements in the mesh
      n_fluids = mesh%descriptor%elem%nElems(eT_fluid)

      oversamp_dofs= poly_proj%body_3D%oversamp_dofs

      allocate( velAbs(n_fluids) )

      !do iElem = 1, n_fluids

        ! division by mass density because conservative variables are used
        ! velAbs is used as the square of the velocity vector. This avoids
        ! redundant operations

        velAbs  = (1._rk/statedata%state(:,1,1)**2)                       &
          &        * ( (statedata%state(:,1,2)**2)                        &
          &        +   (statedata%state(:,1,3)**2)                        &
          &        +   (statedata%state(:,1,4)**2)                        &
          &          )

        timestep%LoclinEuler%meanVel = sqrt(velAbs)

        timestep%LoclinEuler%speedOfSound     =                           &
          & sqrt(                                                         &
          &        equation%euler%isen_coef                               &
          &      *(equation%euler%isen_coef - 1)                          &
          &      * ( 1._rk/statedata%state(:,1,1) )                       &
          &      * ( statedata%state(:,1,5)                               &
          &      -   0.5_rk*statedata%state(:,1,1)*velAbs                 &
          &        )                                                      &
          &     )

      !end do

    case('loclineuler_1d')
      ! subroutine that computes the absolute magnitude of mean velocity
      ! and speed of sound for every element

      ! get the number of fluid elements in the mesh
      n_fluids = mesh%descriptor%elem%nElems(eT_fluid)

      oversamp_dofs= poly_proj%body_1D%oversamp_dofs

      allocate( velAbs(n_fluids) )

        ! division by mass density because conservative variables are used

        velAbs  = abs((1._rk/statedata%state(:,1,1))                       &
          &     *       statedata%state(:,1,2)     )

        timestep%LoclinEuler%meanVel = velAbs

        timestep%LoclinEuler%speedOfSound     =                           &
          & sqrt(                                                         &
          &        equation%euler%isen_coef                               &
          &      *(equation%euler%isen_coef - 1)                          &
          &      * ( 1._rk/statedata%state(:,1,1) )                       &
          &      * ( statedata%state(:,1,3)                               &
          &      -   0.5_rk*statedata%state(:,1,1)*(velAbs**2)            &
          &        )                                                      &
          &     )


    case('euler_2d', 'navier_stokes_2d', 'filtered_navier_stokes_2d')

      ! get the number of fluid elements in the mesh
      n_fluids = mesh%descriptor%elem%nElems(eT_fluid)
      gamprod = (equation%euler%isen_coef-1)*equation%euler%isen_coef

      m_or_n_2D: if (control%use_modal_estimate) then
        ! Use modal estimation for timestep limitation.
        elemModal_2D: do iElem = 1, n_fluids
          rho_min = max( statedata%state(iElem,1,1)                  &
            &             - kerneldata%deviation(iElem,1),           &
            &            statedata%state(iElem,1,1)*epsilon(rho_min) )
          m_max(1) = statedata%state(iElem,1,2) &
            &        + kerneldata%deviation(iElem,2)
          m_max(2) = statedata%state(iElem,1,3) &
            &        + kerneldata%deviation(iElem,3)
          m_min(1) = max( statedata%state(iElem,1,2)       &
            &             - kerneldata%deviation(iElem,2), &
            &             0.0_rk                          )
          m_min(2) = max( statedata%state(iElem,1,3)       &
            &             - kerneldata%deviation(iElem,3), &
            &             0.0_rk                          )
          rho_max = statedata%state(iElem,1,1) &
            &       + kerneldata%deviation(iElem,1)
          e_max = statedata%state(iElem,1,4) &
            &     + kerneldata%deviation(iElem,4)
          timestep%euler_2d%maxVel(iElem) = sqrt(m_max(1)**2 + m_max(2)**2) &
            &                            / rho_min
          m_limited = min(e_max, m_min(1)**2 + m_min(2)**2)
          timestep%euler_2d%speedOfSound(iElem)      &
            &  = sqrt( gamprod*(e_max/rho_min        &
            &                   - m_limited/rho_max) )
        end do elemModal_2D
      else m_or_n_2D
        ! Need to compute timestep limitation with nodal values.

        ! detect correct amount of quadpoints and modes due to projection
        ! method
        nquadpoints = poly_proj%body_2D%nquadpoints
        oversamp_dofs= poly_proj%body_2D%oversamp_dofs

        allocate(modalCoeffs(oversamp_dofs,equation%varSys%nScalars))
        allocate(pointVal(nquadpoints,equation%varSys%nScalars))
        allocate(velAbs(nquadpoints))
        elemNodal_2D: do iElem = 1, n_fluids

          ! get the modal coefficients of the current element.
          ! --> modal space
          if (equation%euler%ensure_positivity) then
            call ply_convert2oversample(                           &
              &    state             = statedata%state(iElem,:,:), &
              &    poly_proj         = poly_proj,                  &
              &    nDim              = 2,                          &
              &    modalCoeffs       = modalCoeffs,                &
              &    ensure_positivity = [.true.,                    &
              &                         .false., .false.,          &
              &                         .true. ]                   )
          else
            call ply_convert2oversample(                     &
              &    state       = statedata%state(iElem,:,:), &
              &    poly_proj   = poly_proj,                  &
              &    nDim        = 2,                          &
              &    modalCoeffs = modalCoeffs                 )
          end if
          ! --> oversamp modal space

          ! transform to point values
          !do iVar=1,equation%varSys%nScalars
          call ply_poly_project_m2n(me = poly_proj,                   &
            &                       dim = 2 ,                         &
            &                       nVars = equation%varSys%nScalars, &
            &                       nodal_data= pointVal,             &
            &                       modal_data= modalCoeffs           )
          ! --> oversampling nodal space

          ! transform to primitive variables (in place)
          call atl_eqn_euler_2d_cons2prim(equation, pointVal)

          ! calculate the max velocity and the speed of sound
          ! ... build the absolute velocity
          velAbs(:) = sqrt( (pointVal(:,2))**2 + (pointVal(:,3))**2  )
          ! assign the values to the timestepping
          timestep%euler_2d%maxVel(iElem) = maxval( velAbs(:)  )
          timestep%euler_2d%speedOfSound(iElem) =                            &
            & sqrt(                                                          &
            &   maxval(                                                      &
            &     equation%euler%isen_coef * pointVal(:,4) / pointVal(:,1) ) )
        end do elemNodal_2D

        deallocate(modalCoeffs)
        deallocate(pointVal)
        deallocate(velAbs)

      end if m_or_n_2D


    case('euler_1d')

      ! get the number of fluid elements in the mesh
      n_fluids = mesh%descriptor%elem%nElems(eT_fluid)
      gamprod = (equation%euler%isen_coef-1)*equation%euler%isen_coef

      m_or_n_1D: if (control%use_modal_estimate) then
        ! Use modal estimation for timestep limitation.
        elemModal_1D: do iElem = 1, n_fluids
          rho_min = max( statedata%state(iElem,1,1)                  &
            &             - kerneldata%deviation(iElem,1),           &
            &            statedata%state(iElem,1,1)*epsilon(rho_min) )
          m_max(1) = statedata%state(iElem,1,2) &
            &        + kerneldata%deviation(iElem,2)
          m_min(1) = max( statedata%state(iElem,1,2)       &
            &             - kerneldata%deviation(iElem,2), &
            &             0.0_rk                           )
          rho_max = statedata%state(iElem,1,1) &
            &       + kerneldata%deviation(iElem,1)
          e_max = statedata%state(iElem,1,3) &
            &     + kerneldata%deviation(iElem,3)
          timestep%euler_1d%maxVel(iElem) = m_max(1) / rho_min
          m_limited = min(e_max, m_min(1)**2)
          timestep%euler_1d%speedOfSound(iElem)      &
            &  = sqrt( gamprod*(e_max/rho_min        &
            &                   - m_limited/rho_max) )
        end do elemModal_1D
      else m_or_n_1D
        ! Need to compute timestep limitation with nodal values.

        ! detect correct amount of quadpoints and modes due to projection
        ! method
        nquadpoints = poly_proj%body_1D%nquadpoints
        oversamp_dofs= poly_proj%body_1D%oversamp_dofs

        allocate(modalCoeffs(oversamp_dofs,equation%varSys%nScalars))
        allocate(pointVal(nquadpoints,equation%varSys%nScalars))
        allocate(velAbs(nquadpoints))
        elemNodal_1D: do iElem = 1, n_fluids

          ! get the modal coefficients of the current element.
          ! --> modal space
          if (equation%euler%ensure_positivity) then
            call ply_convert2oversample(                           &
              &    state             = statedata%state(iElem,:,:), &
              &    poly_proj         = poly_proj,                  &
              &    nDim              = 1,                          &
              &    modalCoeffs       = modalCoeffs,                &
              &    ensure_positivity = [.true.,                    &
              &                         .false.,                   &
              &                         .true. ]                   )
          else
            call ply_convert2oversample(                     &
              &    state       = statedata%state(iElem,:,:), &
              &    poly_proj   = poly_proj,                  &
              &    nDim        = 1,                          &
              &    modalCoeffs = modalCoeffs                 )
          end if
          ! --> oversamp modal space

          ! transform to point values
          !do iVar=1,equation%varSys%nScalars
          call ply_poly_project_m2n(me = poly_proj,                   &
            &                       dim = 1 ,                         &
            &                       nVars = equation%varSys%nScalars, &
            &                       nodal_data= pointVal,             &
            &                       modal_data= modalCoeffs           )
          ! --> oversampling nodal space

          ! transform to primitive variables (in place)
          call atl_eqn_euler_1d_cons2prim(equation, pointVal)

          ! calculate the max velocity and the speed of sound
          ! ... build the absolute velocity
          velAbs(:) = sqrt( (pointVal(:,2))**2 )
          ! assign the values to the timestepping
          timestep%euler_1d%maxVel(iElem) = maxval( velAbs(:)  )
          timestep%euler_1d%speedOfSound(iElem) =                            &
            & sqrt(                                                          &
            &   maxval(                                                      &
            &     equation%euler%isen_coef * pointVal(:,3) / pointVal(:,1) ) )
        end do elemNodal_1D

        deallocate(modalCoeffs)
        deallocate(pointVal)
        deallocate(velAbs)

      end if m_or_n_1D

    case default
      write(logUnit(1),*) 'ERROR in atl_initTimeStepInfo:'
      write(logUnit(1),*) 'unknown equation to solve, stopping ...'
      call tem_abort()
    end select

  end subroutine atl_initTimeStepInfo
  ! ****************************************************************************

end module atl_global_time_integration_module
