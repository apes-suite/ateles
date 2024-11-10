! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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
!! Module collecting all informations regarding stabilzation.
module atl_stabilization_module

  use aotus_module,             only: flu_State, aot_get_val
  use aot_table_module,         only: aot_table_open, aot_table_close, aot_table_length

  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit

  use atl_spectral_viscosity_module, only: atl_ini_spectral_visc, atl_spectral_visc_type
  use atl_covolume_module, only: atl_ini_covolume, atl_covolume_type
  use atl_positivity_preserv_module, only: atl_ini_positivity_preserv, atl_positivity_preserv_type
  use atl_cons_positivity_preserv_module, only: atl_ini_cons_positivity_preserv, &
                                              & atl_cons_positivity_preserv_type

  implicit none
  private

  !> Property for no stabilzation
  integer, parameter :: atl_no_stab_prp = 0

  !> Property for spectral viscosity stabilization
  integer, parameter :: atl_spectral_visc_prp = 1

  !> Property for the pointwise positivity preserving limiter
  integer, parameter :: atl_positivity_preserv_prp = 2

  !> Property for the conservative, positivity preserving limiter
  integer, parameter :: atl_cons_positivity_preserv_prp = 3

  !> Property for spectral viscosity stabilization in Chebyshev basis
  integer, parameter :: atl_cheb_spectral_visc_prp = 4

  !> Property for colvolume filter stabilization
  integer, parameter :: atl_covolume_prp = 5

  !> Datatype representing the stabilization procedure
  !! of a numerical scheme.
  type atl_stabilization_type

    !> The stabilization type
    integer :: stab_kind = atl_no_stab_prp

    !> Logical to indicate if neighbor information is required by the
    !! stabilization.
    logical :: reqNeigh = .false.

    !> Parameters of the spectral viscosity propery (if active).
    type(atl_spectral_visc_type) :: spectral_visc

    !> Parameters of the pointwise positivity preserving limiter (if active).
    type(atl_positivity_preserv_type) :: positivity_preserv

    !> Parameters of the conservative positivity preserving limiter (if active).
    type(atl_cons_positivity_preserv_type) :: cons_positivity_preserv

    !> Parameters of the covolume filter (if active).
    type(atl_covolume_type) :: covolume

  end type atl_stabilization_type


  public :: atl_stabilization_type, atl_no_stab_prp, atl_spectral_visc_prp, &
          & atl_ini_stabilization, atl_positivity_preserv_prp, &
          & atl_cons_positivity_preserv_prp, atl_cheb_spectral_visc_prp, &
          & atl_covolume_prp

contains

  subroutine atl_ini_stabilization(conf, parent_table, filter)
    ! --------------------------------------------------------------------------
    !> flu binding to lua configuration file.
    type(flu_State), intent(in) :: conf
    !> The parent table in the config file
    integer, intent(in) :: parent_table
    !> The stabilization type to be initialized
    type(atl_stabilization_type), allocatable, intent(out) :: filter(:)
    ! --------------------------------------------------------------------------
    character(len=128) :: stab_name
    integer :: iError, stab_table, filter_table, nSubtables, iStab
    ! --------------------------------------------------------------------------

    ! open the stabilization subtable
    call aot_table_open(L=conf, parent=parent_table, &
      &                 thandle=stab_table, key='stabilization')

    ! try to read the field name
    call aot_get_val(L = conf, thandle = stab_table, &
      &              key = 'name', &
      &              val = stab_name, &
      &              default = 'none', &
      &              ErrCode = iError)

    ! Check if we were able to read the name entry, if not we proceed
    ! with subtables
    if(iError.eq.0) then

      allocate(filter(1))
      iStab = 1

      select case(stab_name)
      case('none')
        filter(iStab)%stab_kind = atl_no_stab_prp
      case('spectral_viscosity')
        filter(iStab)%stab_kind = atl_spectral_visc_prp
        call atl_ini_spectral_visc( conf = conf, parent_table = stab_table, &
                                  & filter = filter(iStab)%spectral_visc )
      case('cheb_spectral_viscosity')
        filter(iStab)%stab_kind = atl_cheb_spectral_visc_prp
        call atl_ini_spectral_visc( conf = conf, parent_table = stab_table, &
                                  & filter = filter(iStab)%spectral_visc )
      case('positivity_preserv')
        filter(iStab)%stab_kind = atl_positivity_preserv_prp
        call atl_ini_positivity_preserv( conf = conf, parent_table = stab_table, &
                                  & filter = filter(iStab)%positivity_preserv )
      case('cons_positivity_preserv')
        filter(iStab)%stab_kind = atl_cons_positivity_preserv_prp
        call atl_ini_cons_positivity_preserv( conf = conf, parent_table = stab_table, &
                                  & filter = filter(iStab)%cons_positivity_preserv )
      case('covolume')
        filter(iStab)%stab_kind = atl_covolume_prp
        filter(iStab)%reqNeigh = .true.
        call atl_ini_covolume( conf = conf, parent_table = stab_table, &
                                  & filter = filter(iStab)%covolume )
      case default
        write(logUnit(1),*) 'ERROR in atl_ini_stabilization: Unknown stabilization ' // &
          & 'procedure, stopping ... '
        call tem_abort()
      end select

    else

      ! Get the number of subtables in stabilization
      nSubtables = aot_table_length(L = conf , thandle = stab_table )
      if(nSubtables>0) then
        allocate(filter(nSubtables))
      else
        allocate(filter(0))
      end if

      do iStab = 1, nSubtables

        call aot_table_open(L=conf, parent=stab_table, &
         &                 thandle=filter_table, pos= iStab )

        ! Read the name of the stabilization method
        ! get the name of the scheme
        call aot_get_val(L = conf, thandle = filter_table, &
          &              key = 'name', &
          &              val = stab_name, &
          &              default = 'none', &
          &              ErrCode = iError)
        select case(stab_name)
        case('none')
          filter(iStab)%stab_kind = atl_no_stab_prp
        case('spectral_viscosity')
          filter(iStab)%stab_kind = atl_spectral_visc_prp
          call atl_ini_spectral_visc( conf = conf, parent_table = filter_table, &
                                    & filter = filter(iStab)%spectral_visc )
        case('positivity_preserv')
          filter(iStab)%stab_kind = atl_positivity_preserv_prp
          call atl_ini_positivity_preserv( conf = conf, parent_table = filter_table, &
                                    & filter = filter(iStab)%positivity_preserv )
        case('cons_positivity_preserv')
          filter(iStab)%stab_kind = atl_cons_positivity_preserv_prp
          call atl_ini_cons_positivity_preserv( conf = conf, parent_table = filter_table, &
                                    & filter = filter(iStab)%cons_positivity_preserv )
        case('covolume')
          filter(iStab)%stab_kind = atl_covolume_prp
          filter(iStab)%reqNeigh = .true.
          call atl_ini_covolume( conf = conf, parent_table = filter_table, &
                                    & filter = filter(iStab)%covolume )
        case default
          write(logUnit(1),*) 'ERROR in atl_ini_stabilization: Unknown stabilization ' // &
            & 'procedure, stopping ... '
          call tem_abort()
        end select

        call aot_table_close(L = conf, thandle = filter_table)
      end do

    end if

    call aot_table_close(L = conf, thandle = stab_table)

  end subroutine atl_ini_stabilization

end module atl_stabilization_module
