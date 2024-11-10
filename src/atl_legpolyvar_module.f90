! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> Definition of a spatial function based on a polynomial description.
!!
!! This module provides a variable that can be used to describe a function with
!! a Legendre series, as found in single elements of the solver.
!! The Legendre modes are read from a file in restart format, but only a single
!! polynomial can be considered (not a mesh).
module atl_legpolyvar_module
  use iso_c_binding, only: c_ptr, c_null_ptr, c_f_pointer, c_loc
  use mpi
  use aotus_module, only: flu_State, aot_get_val

  use env_module, only: rk, rk_mpi, pathLen, labelLen, newunit
  use tem_aux_module, only: tem_abort
  use tem_logging_module, only: logUnit
  use tem_operation_var_module, only: tem_get_new_varSys_data_ptr, &
    &                                 tem_free_varSys_data_ptr
  use tem_variable_module, only: tem_variable_type
  use tem_varSys_module, only: tem_varSys_type,                     &
    &                          tem_varSys_solverData_evalElem_type, &
    &                          tem_varSys_proc_point,               &
    &                          tem_varSys_proc_element,             &
    &                          tem_varSys_proc_setParams,           &
    &                          tem_varSys_proc_getParams,           &
    &                          tem_varSys_proc_setupIndices,        &
    &                          tem_varSys_proc_getValOfIndex,       &
    &                          tem_varSys_append_derVar

  use ply_dof_module, only: Q_space, P_space, ply_degree_2dof

  implicit none

  private

  !> Configuration for a legpolyvar.
  !!
  !! The legpolyvar lets you describe a spatial function with the help of a
  !! multidimensional Legendre polynomial series, as found as basis functions
  !! in the solver elements.
  !!
  !! The Legendre modes are read from a file in restart format (all modes
  !! written unformatted consecutively into a file).
  type, public :: atl_legpolyvar_type
    !> Polynomial space of the multidimensional polynomial (Q or P)
    integer :: poly_space

    !> Dimensionality of the polynomial (needs to be between 1 and 3).
    !!
    !! onedimensional polynomial series only vary in X,
    !! twodimensional in X and Y.
    integer :: nDims

    !> Maximal polynomial degree in the polynomial series.
    integer :: degree

    !> State component in the given file for the data to read.
    !!
    !! The polynomial definition in the file may contain multiple quantities.
    !! iComp specifies, which of them to use for this function.
    integer :: iComp

    !> Position of the origin corner for the box, the Legendre polynomials
    !! are to live in.
    real(kind=rk) :: origin(3)

    !> Extent of the box the Legendre polynomial series is defined in.
    real(kind=rk) :: extent

    !> Name of the file to read the polynomial coefficients from.
    !!
    !! This file needs to contain the polynomial data in the same format as
    !! in the Ateles restart files, but only a single element is considered.
    !! Layout of the data is given by poly_space, nDims and degree.
    character(len=pathLen) :: filename
  end type atl_legpolyvar_type

  public :: atl_legpolyvar_load
  public :: atl_legpolyvar_append


  ! Private parts

  type polydata_type
    !> Maximal polynomial degree in the polynomial series.
    integer :: degree

    !> Position of the origin corner for the box, the Legendre polynomials
    !! are to live in.
    real(kind=rk) :: origin(3)

    !> Extent of the box the Legendre polynomial series is defined in.
    real(kind=rk) :: extent

    !> Legendre modes.
    real(kind=rk), allocatable :: modes(:)
  end type polydata_type


contains


  ! ************************************************************************* !
  !> Load the definition of a Legendre polynomial variable from a Lua script.
  subroutine atl_legpolyvar_load(L, parent, legpolyvar)
    !> Lua script to load the polyvar definition from.
    type(flu_State) :: L

    !> Parent table in the Lua script to read the variable data from.
    integer, intent(in) :: parent

    !> Resulting Legendre polynomial description to fill.
    type(atl_legpolyvar_type), intent(out) :: legpolyvar
    ! --------------------------------------------------------------------- !
    integer :: iError
    integer :: vError(3)
    character :: poly_space_char
    ! --------------------------------------------------------------------- !

    call aot_get_val( L       = L,                &
      &               thandle = parent,           &
      &               key     = 'ndims',          &
      &               val     = legpolyvar%nDims, &
      &               ErrCode = iError            )

    if (iError /= 0) then
      write(logUnit(1),*) 'Could not read ndims for legpolyvar.'
      call tem_abort('ERROR in reading legpolyvar definition.')
    end if

    call aot_get_val( L       = L,                 &
      &               thandle = parent,            &
      &               key     = 'degree',          &
      &               val     = legpolyvar%degree, &
      &               ErrCode = iError             )

    if (iError /= 0) then
      write(logUnit(1),*) 'Could not read degree for legpolyvar.'
      call tem_abort('ERROR in reading legpolyvar definition.')
    end if

    call aot_get_val( L       = L,               &
      &               thandle = parent,          &
      &               key     = 'poly_space',    &
      &               val     = poly_space_char, &
      &               ErrCode = iError           )

    if (iError /= 0) then
      write(logUnit(1),*) 'Could not read poly_space for legpolyvar.'
      call tem_abort('ERROR in reading legpolyvar definition.')
    end if

    select case(poly_space_char)
    case ('q','Q')
      legpolyvar%poly_space = Q_space
    case ('p','P')
      legpolyvar%poly_space = P_space
    case default
      write(logUnit(1),*) 'Unknown polynomial space: ', poly_space_char
      write(logUnit(1),*) 'Should be either Q or P!'
      call tem_abort('ERROR in reading legpolyvar definition.')
    end select

    call aot_get_val( L       = L,                &
      &               thandle = parent,           &
      &               key     = 'component',      &
      &               val     = legpolyvar%iComp, &
      &               ErrCode = iError            )

    if (iError /= 0) then
      write(logUnit(1),*) 'Could not read component for legpolyvar.'
      call tem_abort('ERROR in reading legpolyvar definition.')
    end if

    call aot_get_val( L       = L,                        &
      &               thandle = parent,                   &
      &               key     = 'origin',                 &
      &               val     = legpolyvar%origin,        &
      &               default = [0.0_rk, 0.0_rk, 0.0_rk], &
      &               ErrCode = vError                    )

    if (any(vError(:legpolyvar%nDims) /= 0)) then
      write(logUnit(1),*) 'Could not read origin for legpolyvar.'
      call tem_abort('ERROR in reading legpolyvar definition.')
    end if

    call aot_get_val( L       = L,                 &
      &               thandle = parent,            &
      &               key     = 'extent',          &
      &               val     = legpolyvar%extent, &
      &               ErrCode = iError             )

    if (iError /= 0) then
      write(logUnit(1),*) 'Could not read extent for legpolyvar.'
      call tem_abort('ERROR in reading legpolyvar definition.')
    end if

    call aot_get_val( L       = L,                   &
      &               thandle = parent,              &
      &               key     = 'filename',          &
      &               val     = legpolyvar%filename, &
      &               ErrCode = iError               )

    if (iError /= 0) then
      write(logUnit(1),*) 'Could not read filename for legpolyvar.'
      call tem_abort('ERROR in reading legpolyvar definition.')
    end if


  end subroutine atl_legpolyvar_load


  subroutine atl_legpolyvar_append(var, varSys, pos, solverData_evalElem)
    !> Data describing the variable to append (needs to be a legpolyvar).
    class(tem_variable_type), intent(in) :: var

    !> Variable system to append the variable to.
    type(tem_varSys_type), intent(inout) :: varsys

    !> Position of the appended variable in the system.
    integer, optional, intent(out) :: pos

    !> A setter to allow the caller to define a routine for the construction
    !! of an element representation.
    type(tem_varSys_solverData_evalElem_type), &
      &  optional, intent(in) :: solverData_evalElem
    ! -------------------------------------------------------------------- !
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => NULL()
    type(atl_legpolyvar_type), pointer :: legpolyvar
    type(polydata_type), pointer :: polydat
    type(c_ptr) ::  method_data = c_null_ptr
    type(c_ptr) ::  c_polydat = c_null_ptr
    character(len=labelLen) :: input_varname(0)
    integer :: input_varIndex(0)
    integer :: addedPos
    logical :: wasAdded
    ! -------------------------------------------------------------------- !

    call c_f_pointer(var%solver_specifics, legpolyvar)

    allocate(polydat)
    polydat%degree = legpolyvar%degree
    polydat%origin = legpolyvar%origin
    polydat%extent = legpolyvar%extent

    c_polydat = c_loc(polydat)

    method_data = tem_get_new_varSys_data_ptr(c_polydat)

    call tem_varSys_append_derVar(          &
      &    me             = varSys,         &
      &    varname        = var%label,      &
      &    opertype       = var%opertype,   &
      &    nComponents    = 1,              &
      &    input_varname  = input_varname,  &
      &    input_varIndex = input_varIndex, &
      &    method_data    = method_data,    &
      &    get_point      = get_point,      &
      &    get_element    = get_element,    &
      &    set_params     = set_params,     &
      &    get_params     = get_params,     &
      &    setup_indices  = setup_indices,  &
      &    get_valOfIndex = get_valOfIndex, &
      &    pos            = addedPos,       &
      &    wasAdded       = wasAdded        )

    if (wasAdded) then
      call read_legpolyvar_modes( legpolyvar = legpolyvar,   &
        &                         modes      = polydat%modes )
      if (present(solverData_evalElem)) then
        call solverData_evalElem%opVar_setter(varSys%method%val(addedPos))
      end if
      write(logUnit(8),*) 'Successfully appended variable "' &
        &                 // trim(var%label) // '" to the variable system'
    else
      deallocate(polydat)
      call tem_free_varSys_data_ptr(method_data)
      if (addedPos < 1) then
        write(logUnit(7),*) 'WARNING: variable '//trim(var%label)// &
          &                 ' is not added to the variable system!'
      end if
    end if

    if (present(pos)) pos = addedPos
  end subroutine atl_legpolyvar_append


  ! Private parts

  subroutine read_legpolyvar_modes( legpolyvar, modes )
    type(atl_legpolyvar_type), intent(in) :: legpolyvar
    real(kind=rk), allocatable, intent(out) :: modes(:)
    ! -------------------------------------------------------------------- !
    integer :: iError
    integer :: myrank
    integer :: funit
    integer :: polyrecl
    integer :: ndofs
    ! -------------------------------------------------------------------- !

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, iError)

    ndofs = ply_degree_2dof( deg   = legpolyvar%degree,     &
      &                      space = legpolyvar%poly_space, &
      &                      ndims = legpolyvar%ndims       )
    allocate(modes(ndofs))

    if (myrank == 0) then
      inquire(iolength=polyrecl) modes
      funit = newunit()
      open( file   = trim(legpolyvar%filename), &
        &   unit   = funit,                     &
        &   access = 'DIRECT',                  &
        &   action = 'READ',                    &
        &   form   = 'UNFORMATTED',             &
        &   recl   = polyrecl,                  &
        &   iostat = iError                     )

      if (iError == 0) then
        read(unit=funit, rec=legpolyvar%iComp, iostat=iError) modes
        if (iError /= 0) then
          write(*,*) 'ERROR while trying to read ', ndofs, ' modes from ' &
            &        // trim(legpolyvar%filename)
          write(*,*) 'Attempted to read record ', legpolyvar%iComp
        end if
      else
        write(*,*) 'ERROR in opening file '//trim(legpolyvar%filename)
      end if

      if (iError /= 0) then
        call tem_abort('ERROR while reading data for legpolyvar from file.')
      end if

      close(funit)
    end if

    call MPI_Bcast(modes, ndofs, rk_mpi, 0, MPI_COMM_WORLD, iError)

  end subroutine read_legpolyvar_modes

end module atl_legpolyvar_module
