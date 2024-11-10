! Copyright (c) 2014-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016, 2018, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

!> The linearized Euler equations of compressible inviscid flows.
!!
!! These are the Euler equations but linearized around a global background
!! state.
!! In contrast to the nonlinear Euler equations in [[atl_eqn_euler_module]],
!! The state is represented in terms of the primitive variables density,
!! velocity and pressure.
!! To define the properties of the fluid the following parameters need to be
!! defined:
!!
!! * The isentropic expansion coefficient `isen_coef`
!! * The `background` state around the Euler equations are linearized each
!!   state variable needs to be given, and may be a function of time, see
!!   [[tem_temporal_module:tem_load_temporal]]
!!   * `density` the background density
!!   * `velocityX` the background velocity in X direction
!!   * `velocityY` the background velocity in Y direction
!!   * If 3D: `velocityZ` for the background velocity in Z direction
!!   * `pressure` the background pressure
!!
!! Hence, the definition for the linearized euler equation takes a form like
!!
!!```lua
!!  equation = {
!!    name   = 'linearEuler',
!!    isen_coef = 1.4,
!!    background = {
!!      density = 1.225,
!!      velocityX = 100.0,
!!      velocityY = 0.0,
!!      velocityZ = 0.0,
!!      pressure = 100000.0
!!    }
!!  }
!!```
!!
!! The following equation names are implementing the linearized Euler equations:
!!
!! * `lineareuler`
!! * `lineareuler_2d`
!!
module atl_eqn_LinearEuler_module
  ! Treelm modules
  use env_module,             only: rk
  use tem_aux_module,         only: tem_abort
  use tem_logging_module,     only: logUnit
  use tem_tools_module,       only: tem_horizontalSpacer
  use tem_temporal_module,    only: tem_temporal_type, &
    &                               tem_load_temporal, &
    &                               tem_temporal_for
  use tem_time_module,        only: tem_time_type

  ! Aotus modules
  use aotus_module,           only: flu_State, aot_get_val
  use aot_out_module,         only: aot_out_type, aot_out_val
  use aot_table_module,       only: aot_table_open, aot_table_close

  ! Ateles modules
  use atl_materialFun_module, only: atl_materialFun_type

  implicit none

  private

  public :: atl_LinearEuler_type
  public :: atl_load_LinearEuler
  public :: atl_save_LinearEuler
  public :: atl_eqn_update_background
  public :: atl_lineuler_numflux

  !> Type to store the temporal function for each background state
  type temporal_background_type
    type(tem_temporal_type) :: density
    type(tem_temporal_type) :: velocityX
    type(tem_temporal_type) :: velocityY
    type(tem_temporal_type) :: velocityZ
    type(tem_temporal_type) :: pressure
  end type temporal_background_type

  type dir_proc
    !> Procedure to compute the numerical flux for the equation at hand.
    !!
    !! What kind of fluxes are available depends on the equation that is
    !! being solved.
    procedure(atl_lineuler_numflux), pointer, nopass :: numflux => NULL()
  end type dir_proc

  !> The Euler equation properties are stored here
  type atl_LinearEuler_type
    !> isentropic coefficient
    real(kind=rk) :: isen_coef
    !> background density
    real(kind=rk) :: density_0
    !> background velocity (x,y,z) direction
    real(kind=rk),allocatable :: velocity_0(:)
    !> background pressure
    real(kind=rk) :: pressure_0
    !> speedofSound, depends on temporal background
    ! need to be updated every timestep
    real(kind=rk) :: speedOfSound
    !> type for the temporal function of background, used to update background
    ! in every timestep
    type(temporal_background_type) :: temporal_background
    !> The functions for the penalizations
    type(atl_materialFun_type) :: penalization
    !> type for direction specific procedure like the numerical flux
    type(dir_proc) :: dir_proc(3)
  end type atl_LinearEuler_type

  ! ------------------------------------------------------------------------ !
  abstract interface
    !> Interface definition for numerical fluxes.
    !! @todo HK: should be vectorized to reduce overheads.
    subroutine atl_lineuler_numflux( nSides, nFaceDofs,faceRep, faceFlux,      &
      &                              leftPos, rightPos, var, LinearEuler, iDir )
      import :: atl_LinearEuler_type, rk
      !> Datatype for LinearEuler equation include all background data
      type(atl_LinearEuler_type), intent(in) :: LinearEuler
      integer, intent(in) :: nFaceDofs, nSides
      real(kind=rk), intent(in) :: faceRep(:, :, :, :)
      real(kind=rk), intent(inout) :: faceFlux(:, :, :, :)
      integer, intent(in) :: leftPos(nSides), rightPos(nsides)
      integer, intent(in) :: var(:)
      !> Direction of the flow, used for background velocity
      integer, intent(in) :: idir
    end subroutine atl_lineuler_numflux
  end interface
  ! ------------------------------------------------------------------------ !


contains


  ! ------------------------------------------------------------------------ !
  !> subroutine to initialize an equation of type linear euler equation

  !! as defined in the configuration file
  subroutine atl_load_LinearEuler( linearEuler, conf, eq_table, spatial_dim )
    ! -------------------------------------------------------------------- !
    !> Resulting description of the Euler equation parameters.
    type(atl_LinearEuler_type), intent(out) :: LinearEuler

    !> Handle to the configuration script, to load the parameters from.
    type(flu_State) :: conf

    !> Handle to the table containing the description for the equation
    !! system.
    integer, intent(in) :: eq_table

    !> The spatial dimension of the Euler equation
    integer :: spatial_dim
    ! -------------------------------------------------------------------- !
    integer :: iError, LinearEuler_table
    ! -------------------------------------------------------------------- !

    ! allocate the dimension of background velocity array according to dimension
    allocate(LinearEuler%velocity_0(spatial_dim))

    !read the data from the equation table of the lua file
    call aot_get_val(L = conf, thandle = eq_table, key = 'isen_coef', &
      &              val = LinearEuler%isen_coef, &
      &              ErrCode = iError)
    if(iError /= 0) then
      call tem_abort( 'ERROR: not able to find isentropic coefficient in ' &
        & // 'equation table, stopping ...'                                )
    end if

    call tem_horizontalSpacer(funit=logUnit(5))
    write(logUnit(5),*) 'Loading background parameters for linearized euler ' &
      & // 'equation from config file. '

    ! Open subtable for background properties
    call aot_table_open( L       = conf,              &
      &                  parent  = eq_table,          &
      &                  tHandle = LinearEuler_table, &
      &                  key     = 'background'       )

    if(LinearEuler_table.eq.0) then
      write(logUnit(1),*) 'ERROR in init_LinearEuler: no background ' &
        & // 'properties defined, stopping ...'
    end if

    ! load the data from the equation table of the lua file
    call tem_load_temporal (                               &
      &  me     = LinearEuler%temporal_background%density, &
      &  conf   = conf,                                    &
      &  parent = LinearEuler_table,                       &
      &  key    = 'density'                                )
    call tem_load_temporal (                                 &
      &  me     = LinearEuler%temporal_background%velocityX, &
      &  conf   = conf,                                      &
      &  parent = LinearEuler_table,                         &
      &  key    = 'velocityX'                                )
    call tem_load_temporal (                                 &
      &  me     = LinearEuler%temporal_background%velocityY, &
      &  conf   = conf,                                      &
      &  parent = LinearEuler_table,                         &
      &  key    = 'velocityY'                                )
    if (spatial_dim == 3) then
      call tem_load_temporal (                                 &
        &  me     = LinearEuler%temporal_background%velocityZ, &
        &  conf   = conf,                                      &
        &  parent = LinearEuler_table,                         &
        &  key    = 'velocityZ'                                )
    end if
    call tem_load_temporal ( &
      &  me     = LinearEuler%temporal_background%pressure, &
      &  conf   = conf,                                     &
      &  parent = LinearEuler_table,                        &
      &  key    = 'pressure'                                )

    ! Close the Lua table with the background information
    call aot_table_close( L = conf, tHandle = LinearEuler_table )


  end subroutine atl_load_LinearEuler
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! dump the equation variables into the lua file
  subroutine atl_save_LinearEuler(me, eqn_name, nDimensions, conf)
    ! -------------------------------------------------------------------- !
    type(atl_LinearEuler_type), intent(in) :: me
    character(len=*), intent(in) :: eqn_name
    integer, intent(in) :: nDimensions
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    call aot_out_val( put_conf = conf,          &
      &               vname    = 'name',        &
      &               val      = trim(eqn_name) )

    ! Dump equation Properties
    call aot_out_val( put_conf = conf,        &
      &               vname    = 'isen_coef', &
      &               val      = me%isen_coef )
    call aot_out_val( put_conf = conf,        &
      &               vname    = 'density',   &
      &               val      = me%density_0 )
    call aot_out_val( put_conf = conf,            &
      &               vname    = 'velocityX',     &
      &               val      = me%velocity_0(1) )
    call aot_out_val( put_conf = conf,            &
      &               vname    = 'velocityY',     &
      &               val      = me%velocity_0(2) )
    if (nDimensions .eq. 3) then
      call aot_out_val( put_conf = conf,            &
        &               vname    = 'velocityZ',     &
        &               val      = me%velocity_0(3) )
    end if
    call aot_out_val( put_conf = conf,         &
      &               vname    = 'pressure',   &
      &               val      = me%pressure_0 )

  end subroutine atl_save_LinearEuler
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Routine which updates the background since it is a temporal function and
  ! vary in time
  subroutine atl_eqn_update_background( me, time, nDimensions)
    ! -------------------------------------------------------------------- !
    !> linearEuler type including background
    type(atl_LinearEuler_type), intent(inout) :: me
    !> timer object incl. the current time information
    type(tem_time_type), intent( in )  :: time
    !> spatial dimension
    integer, intent(in) :: nDimensions
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    write(logUnit(5),*) "Update background for linear euler equation"

    ! background density
    me%density_0 = tem_temporal_for(               &
      & temporal = me%temporal_background%density, &
      & time     = time                            )
    ! background velocity X
    me%velocity_0(1) = tem_temporal_for(                 &
      & temporal     = me%temporal_background%velocityX, &
      & time         = time                              )
    ! background velocity Y
    me%velocity_0(2) = tem_temporal_for(                 &
      & temporal     = me%temporal_background%velocityY, &
      & time         = time                              )
    if (nDimensions == 3) then
      ! background velocity Z
      me%velocity_0(3) = tem_temporal_for(                 &
        & temporal     = me%temporal_background%velocityZ, &
        & time         = time                              )
    end if
    ! background pressure
    me%pressure_0 = tem_temporal_for(               &
      & temporal = me%temporal_background%pressure, &
      & time     = time                             )

    ! calculate speed of sound
    me%SpeedOfSound = sqrt( me%isen_coef * me%pressure_0/ me%density_0 )

    write(logUnit(10),*) 'Updated Background of linearized Euler equation'
    write(logUnit(10),*) ' * isen_coef: ', me%isen_coef
    write(logUnit(10),*) ' * background density: ', me%density_0
    write(logUnit(10),*) ' * background velocityX: ', me%velocity_0(1)
    write(logUnit(10),*) ' * background velocityY: ', me%velocity_0(2)
    if ( nDimensions == 3) then
      write(logUnit(10),*) ' * background velocity_Z: ', me%velocity_0(3)
    end if
    write(logUnit(10),*) ' * background pressure: ', me%pressure_0
    write(logUnit(10),*) ' * speed of sound: ', me%speedofsound

  end subroutine atl_eqn_update_background
  ! ------------------------------------------------------------------------ !

end module atl_eqn_LinearEuler_module
