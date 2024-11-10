! Copyright (c) 2012-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013, 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

!> Boundary conditions.
!!
!! For each boundary label found in the mesh, there needs to be a definition of
!! the boundary condition to apply at the corresponding faces.
!! What kind of boundary conditions are available depends on the equation
!! system, but each boundary condition needs to be identified by a label, that
!! matches the one defined in the mesh.
!!
!! If there are values that are to be extrapolated (Neumann boundary
!! condition), you can set `enforce_zero_grad` to true, to use an extrapolation
!! of a polynomial with zero gradient at the boundary.
!! This is achieved by computing the last mode to fulfill this condition.
!! If you set `neumann_mode_fraction` to a smaller value than 1, then only
!! this fraction of lower modes will be used in the `enforce_zero_grad`
!! procedure and higher modes will be set to 0.
!!
!! Simple example with two different boundary conditions, one without further
!! options (slipwall) and one with further settings (inflow):
!!
!!```lua
!! boundary condition = {
!!    { -- SLIPWALL
!!      --   Velocity in normal direction 0, other values extrapolated.
!!      label = 'cylinder',
!!      kind  = 'slipwall', -- or 'wall'
!!      enforce_zero_grad = true,
!!      neumann_mode_fraction = 1.0
!!    },
!!    { -- INFLOW
!!      --   Prescribe density and velocity, extrapolate pressure.
!!      label = 'outside',
!!      kind = 'conservatives',
!!      density = 1.23,
!!      velocityX = 0.2,
!!      velocityY = 0.3,
!!      velocityZ = 0.4
!!      enforce_zero_grad = true,
!!      neumann_mode_fraction = 1.0
!!    }
!!  }
!!```
!!
!! If there is a boundary in the mesh, for which no boundary condition is
!! defined the application will stop.
!! Boundary conditions that are defined, but for which there is no label in
!! the mesh, will be ignored.
module atl_bc_header_module
  use, intrinsic :: iso_c_binding,    only: c_f_pointer

  use env_module,                     only: LabelLen, rk
  use tem_bc_module,                  only: tem_bc_state_type
  use tem_bc_header_module,           only: tem_bc_header_type, &
    &                                       tem_load_bc_header
  use tem_bc_prop_module,             only: tem_bc_prop_type
  use tem_aux_module,                 only: tem_abort
  use tem_logging_module,             only: logUnit
  use tem_stringKeyValuePair_module,  only: grw_stringKeyValuePairArray_type
  use tem_varSys_module,              only: tem_varSys_type
  use tem_dyn_array_module,           only: PositionOfVal

  use aotus_module,                   only: flu_State, &
    &                                       aoterr_Fatal
  use aot_table_module,               only: aot_table_open,  &
    &                                       aot_table_close, &
    &                                       aot_get_val

  use atl_equation_module,            only: atl_equations_type, &
    &                                       atl_eqn_var_trafo_type


  implicit none

  private

  public :: atl_boundary_type, atl_load_bc
  public :: atl_store_bcVarPos


  !> This type describes a single boundary condition, which is described in
  !! the configuration files and attached to elements in the mesh file.
  type atl_boundary_type

    !> A label identifying this boundary condition
    character(len=LabelLen) :: label

    !> The kind of this boundary condition, mainly used to describe predefined
    !! boundary conditions with some default settings for some of the variables.
    character(len=LabelLen) :: BC_kind

    !> Method to use for the extrapolation of Neumann boundaries.
    !!
    !! Enforce_zero_grad indicates, that the value to set for the extrapolated
    !! value of Neumann boundary conditions will be obtained by correcting the
    !! polynomial in face normal direction to have a zero gradient at the
    !! face. Otherwise the value at the face of the unmodified polynomial will
    !! be used.
    !! If this is true, the Neumann_mode_fraction will indicate how many modes
    !! to use for the extrapolation.
    logical :: enforce_zero_grad

    !> Fraction of modes to use for the extrapolation to obtain Neumann BCs.
    !!
    !! Has to be a value between 0 and 1, where 0 means, just the integral mean
    !! (first mode) will be used for the extrapolation, and 1 means all modes
    !! except the last one are used. The last mode will be computed such, that
    !! the gradient at the boundary is null. Then the resulting value on the
    !! face is used as extrapolated value.
    real(kind=rk) :: neumann_mode_fraction

    !> A flag to indicate if vectorial quantities are defined in the boundary
    !! normal system.
    logical :: bc_normal_vec

    !> Boundary condition description for each of the state variables.
    !! The size of this array depends on the equation system and covers all
    !! required variables.
    !! The variables are expected to occur in the same order as in the equation
    !! system.
    !! Furthermore, in case of face normal boundary conditions, we assume that
    !! the variable normal to the face is occurring first in this array.
    type(tem_bc_state_type), allocatable :: state(:)

    !> Dictionary of boundary state variable with
    !! varDict%val()%key is the name of boundary variable and
    !! varDict%val()%value is the name of spacetime function variable
    type(grw_stringKeyValuePairArray_type) :: varDict

    !> Pointer to function for the necessary state variable transformation.
    type(atl_eqn_var_trafo_type) :: bc_trafo

    !> A flag to indicate if derivatives of vectorial quantities are defined in
    !! the boundary normal system.
    logical :: bc_normal_vec_gradient

    !> Pointer to function for the necessary state gradient variable
    !! transformation.
    !!
    !! Undefined if equations does not involve higher order derivatives
    type(atl_eqn_var_trafo_type) :: bc_trafo_gradient

    !> Boundary condition description for each of higher order terms of the
    !! equations.
    !!
    !! The size of this array depends on the equation system and covers all
    !! required variables and its derivatives.
    !! The variables are expected to occur in the same order as in the equation
    !! system.
    !! Furthermore, in case of face normal boundary conditions, we assume that
    !! the variable normal to the face is occurring first in this array.
    type(tem_bc_state_type), allocatable :: state_gradient(:)

    !> Dictionary of boundary state gradient variable with
    !! varDict%val()%key is the name of boundary variable and
    !! varDict%val()%value is the name of spacetime function variable.
    !!
    !! Odd count refer to state_gradient(:,1) and
    !! Even count refer to state_gradient(:,2)
    type(grw_stringKeyValuePairArray_type) :: varDict_gradient

  end type atl_boundary_type


contains


  ! ************************************************************************ !
  !> Get the boundary configuration.
  subroutine atl_load_bc(bc, bc_header, bc_prop, equation, conf, parent)
    ! -------------------------------------------------------------------- !
    !> Array of boundary conditions, will be allocated in this routine and
    !! has a length according to the number of boundary conditions in the mesh.
    type(atl_boundary_type), allocatable, intent(out) :: bc(:)

    !> The boundary condition header data as given in the mesh.
    type(tem_bc_header_type), intent(out) :: bc_header

    !> The boundary property object, describing the given boundaries in the
    !! mesh, this has to be provided to allow the matching of boundary settings
    !! in the configuration to the boundary set in the mesh.
    type(tem_bc_prop_type), intent(in) :: bc_prop

    !> Description of the equation system, to read the boundary conditions in
    !! dependency on the equation system to solve.
    type(atl_equations_type), intent(inout) :: equation

    !> Lua script to obtain the configuration data from.
    type(flu_State) :: conf

    !> A parent Lua table, in which the boundary conditions are to be found.
    integer, optional, intent(in) :: parent
    ! -------------------------------------------------------------------- !
    ! Local Variables
    integer :: bc_handle
    integer :: sub_thandle
    integer :: iBC
    integer :: bid
    integer :: iError
    ! -------------------------------------------------------------------- !


    write(logUnit(2),*) ' Loading boundary information'

    if (present(parent)) then
      call tem_load_bc_header(me           = bc_header, &
        &                     conf         = conf,      &
        &                     parentHandle = parent,    &
        &                     BC_prop      = bc_prop    )
      call aot_table_open(L       = conf,                &
        &                 parent  = parent,              &
        &                 thandle = bc_handle,           &
        &                 key     = 'boundary_condition' )
    else
      call tem_load_bc_header(me      = bc_header, &
        &                     conf    = conf,      &
        &                     BC_prop = bc_prop    )

      call aot_table_open(L       = conf,                &
        &                 thandle = bc_handle,           &
        &                 key     = 'boundary_condition' )
    end if

    allocate(bc(bc_header%nBCs))

    do iBC=1,bc_header%nBCs

      call aot_table_open(L=conf, parent=bc_handle, thandle=sub_thandle, &
        &                 pos=iBC)

      ! Skip BC definitions that are not present in the mesh.
      if (sub_thandle /= 0 .and. bc_header%BC_ID(iBC)>0) then

        bid = bc_header%BC_ID(iBC)

        ! bc_header%BC_ID(iBC) yields the position of boundary label iBC in the
        ! mesh. (identified by tem_load_bc_header)
        bc(bid)%label   = bc_header%label(iBC)
        bc(bid)%BC_kind = bc_header%BC_kind(iBC)

        ! Use an equation specific function to set the actual boundary.
        ! This load_bc routine allocates the state description and fills it
        ! according to the configuration in the table provided with sub_handle.
        ! Remember, that equation itself is passed into the subroutine
        ! automatically, as load_bc is a component of equation.
        call equation%load_bc(                                         &
          & bc_state               = bc(bid)%state,                    &
          & bc_state_gradient      = bc(bid)%state_gradient,           &
          & bc_varDict             = bc(bc_header%BC_ID(iBC))%varDict, &
          & bc_varDict_gradient    = bc(bc_header%BC_ID(iBC))          &
          &                                      %varDict_gradient,    &
          & bc_normal_vec          = bc(bid)%bc_normal_vec,            &
          & bc_normal_vec_gradient = bc(bid)%bc_normal_vec_gradient,   &
          & bc_trafo               = bc(bid)%bc_trafo,                 &
          & bc_trafo_gradient      = bc(bid)%bc_trafo_gradient,        &
          & bc_label               = bc_header%label(iBC),             &
          & bc_kind                = bc_header%BC_kind(iBC),           &
          & thandle                = sub_thandle,                      &
          & conf                   = conf                              )

        call aot_get_val( L       = conf,                      &
          &               thandle = sub_thandle,               &
          &               val     = bc(bid)%enforce_zero_grad, &
          &               ErrCode = iError,                    &
          &               key     = 'enforce_zero_grad',       &
          &               default = .false.                    )

        if (btest(iError, aoterr_Fatal)) then
          write(logunit(1), *) 'ERROR in reading boundary condition:'
          write(logunit(1), *) 'Unable to get enforce_zero_grad.'
          write(logunit(1), *) 'Note, it has to be a logical!'
          call tem_abort()
        end if

        if (bc(bid)%enforce_zero_grad) then
          write(logUnit(2),*) '     Enforcing 0 gradient for Neumann conditions!'
          call aot_get_val( L       = conf,                          &
            &               thandle = sub_thandle,                   &
            &               val     = bc(bid)%neumann_mode_fraction, &
            &               ErrCode = iError,                        &
            &               key     = 'neumann_mode_fraction',       &
            &               default = 1.0_rk                         )
          write(logUnit(2),*) '     With ', bc(bid)%neumann_mode_fraction, &
            &                 ' of the modes.'
        else
          bc(bid)%neumann_mode_fraction = 0.0_rk
        end if

      end if

      call aot_table_close(L=conf, thandle=sub_thandle)

    end do

    call aot_table_close(L=conf, thandle=bc_handle)


  end subroutine atl_load_bc
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to store position of user variable defined state and
  !! state_gradient boundary variable in bc(iBC)%state(iVar)%varPos
  !! and bc(iBC)%state_gradient
  subroutine atl_store_bcVarPos( bc, varSys )
    ! -------------------------------------------------------------------- !
    !> Array of boundary conditions contains varDict with same size as
    !! state variables
    type(atl_boundary_type), intent(inout) :: bc(:)

    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    ! -------------------------------------------------------------------- !
    integer :: iBC, iVar, nBCs, user_varPos
    character(len=labelLen) :: user_varName
    ! -------------------------------------------------------------------- !
    nBCs = size(bc)

    ! Store position of state variable
    do iBC = 1, nBCs
      do iVar = 1, bc(iBC)%varDict%nVals
        user_varName = bc(iBC)%varDict%val(iVar)%value
        user_varPos = PositionOfVal( me = varSys%varName,     &
          &                          val = trim(user_varName) )

        ! continue only if this variable exist in varSys
        if (user_varPos>0) then
          bc(iBC)%state(iVar)%varPos = user_varPos
        else
          call tem_abort( 'Error: User defined space-time function variable "' &
            & // trim(user_varName)                                            &
            & // '" not found in varSys.'                                      )
        end if

      end do !iVar
    end do !iBC

    ! Store position of state gradient variables
    do iBC = 1, nBCs
      do iVar = 1, bc(iBC)%varDict_gradient%nVals
        user_varName = bc(iBC)%varDict_gradient%val(iVar)%value
        user_varPos = PositionOfVal( me = varSys%varName,     &
          &                          val = trim(user_varName) )

        ! continue only if this variable exist in varSys
        if (user_varPos>0) then
          ! store odd iVar in (iVar,1) and even in (iVar,2)
!!VK          if (mod(iVar,2)==1) then
!!VK            bc(iBC)%state_gradient(iVar,1)%varPos = user_varPos
!!VK          else
            bc(iBC)%state_gradient(iVar)%varPos = user_varPos
!!VK          end if
        else
          call tem_abort( 'Error: User defined space-time function variable "' &
            & // trim(user_varName)                                            &
            & // '" not found in varSys.'                                      )
        end if

      end do !iVar
    end do !iBC

  end subroutine atl_store_bcVarPos
  ! ************************************************************************ !

end module atl_bc_header_module
