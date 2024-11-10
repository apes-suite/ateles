! Copyright (c) 2012-2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
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
!! Routines and datatypes realted to specification and definition
!! of the MODG scheme.
module atl_modg_scheme_module
  use env_module,         only: rk, labelLen

  use aotus_module,       only: flu_State, &
    &                           aot_get_val
  use aot_top_module,     only: aoterr_Fatal, &
    &                           aoterr_NonExistent
  use aot_table_module,   only: aot_table_open,  &
    &                           aot_table_close, &
    &                           aot_table_length

  use tem_aux_module,     only: tem_abort
  use tem_tools_module,   only: upper_to_lower
  use tem_logging_module, only: logUnit

  use ply_dof_module,     only: Q_space, P_space

  implicit none

  private

  !> Datatype for the modg scheme
  type atl_modg_scheme_type
    !> The maximal polynomial degree per spatial direction.
    integer :: maxPolyDegree
    !> Polynomial basis type.
    !!
    !! 3D Monomials have the form x^i * y^j * z^k
    !! - Q_space: quadratic polynomial space (i,j,k) <= maxPolyDegree
    !! - P_space: polynomial space i+j+k <= maxPolyDegree
    integer :: basisType
    !> Facial Chebyshev nodes (reference element) for all three spatial
    !! direction and left and right face.
    !! These points are necessary to transfer boundary conditions given
    !! in physical space to modal space by
    !! means of an FPT.
  !!  type(atl_faceChebNodes_type), allocatable :: chebNodesFace(:,:)
  end type atl_modg_scheme_type

  public :: atl_modg_scheme_type
  public :: atl_modg_scheme_init

contains

  subroutine atl_modg_scheme_init(me, nDofs, nFaceDofs, &
    &                             conf, thandle, currentLevel, maxLevel )
    type(atl_modg_scheme_type), intent(out) :: me
    integer, intent(out) :: nDofs
    integer, intent(out) :: nFaceDofs
    type(flu_State) :: conf
    integer, intent(in) :: thandle
    !> The current level of the mesh
    integer, intent(in) :: currentLevel
    !> The maximal level of the mesh
    integer, intent(in) :: maxLevel

    character(len=labelLen) :: modg_space
    character(len=labelLen) :: polyspace
    integer :: degree
    integer :: iError

    ! Find the maximal polynomial degree for the current level.
    call atl_modg_scheme_load_polyDegree(degree, conf, thandle, &
      &                                  currentlevel, maxlevel)

    ! Read the chosen polynomial space
    call aot_get_val(L = conf, thandle = thandle, &
      &              key = 'modg_space', &
      &              val = modg_space, &
      &              ErrCode = iError, &
      &              default = 'Q' )

    polyspace = upper_to_lower(modg_space)
    polyspace = adjustl(polyspace)

    select case(trim(polyspace))
    case('q')
      me%basisType = Q_space

      nDofs = (degree+1)**3
      nFaceDofs = (degree+1)**2
      me%maxPolyDegree = degree

    case('p')
      me%basisType = P_space

      nDofs = ((degree+1)*(degree+2)*(degree+3))/6
      nFaceDofs = ((degree+1)*(degree+2))/2
      me%maxPolyDegree = degree

    case default
      write(logUnit(1),*) 'ERROR in init_kernel_state: unknown modg_space: ' &
        &            // trim(modg_space) // ' !'
      write(logUnit(1),*) 'Supported are:'
      write(logUnit(1),*) '* Q (quadratic with i,j,k <= maxDegree'
      write(logUnit(1),*) '* P (with i+j+k <= maxDegree'
      write(logUnit(1),*) 'Stopping....'
      call tem_abort()
    end select

    write(logUnit(1),*) 'MODG scheme with '//trim(modg_space)//' basis on level ', &
      &            currentLevel
    write(logUnit(1),*) '... maximal polynomial degree: ', degree

  end subroutine atl_modg_scheme_init



  subroutine atl_modg_scheme_load_polyDegree(degree, conf, thandle, &
    &                                        currentlevel, maxlevel)
    integer, intent(out) :: degree
    type(flu_state) :: conf
    integer, intent(in) :: thandle
    integer, intent(in) :: currentLevel
    integer, intent(in) :: maxLevel

    character(len=labelLen) :: predefDegrees
    character(len=labelLen) :: predeg
    integer :: mtable
    integer :: degTable
    integer :: iError
    integer :: nDegTables
    integer :: iTable
    integer :: tableLevel
    logical :: foundLevel
    real(kind=rk) :: levelFact, base_order

    ! Get the degree of the scheme
    call aot_get_val(L = conf, thandle = thandle, &
      &              key = 'm', &
      &              val = degree, &
      &              ErrCode = iError)

    ! if no error occurs, we have a global polynomial degree
    if (iError.ne.0) then
      ! open the spatial subtable
      call aot_table_open(L=conf, parent=thandle, &
        &                 thandle = mtable, key = 'm')

      call aot_get_val(L = conf, thandle = mtable, key = 'predefined', &
        &              val = predefDegrees, ErrCode = iError)

      predef: if (btest(iError, aoterr_nonExistent)) then
        ! No predefined set of polynomial degrees, try to get them
        ! individually for each level.
        ! The number of tables in the m table
        nDegTables = aot_table_length(L = conf, thandle = mTable)

        ! Iterate over all the sub tables and check if one has
        ! the current mesh level
        foundLevel = .false.
        do iTable = 1, nDegTables
          call aot_table_open(L = conf, parent = mTable, thandle = degTable, &
            &                 pos = iTable)

          ! get the level of the current table
          call aot_get_val( L = conf, thandle = degTable, key = 'level', &
            &               val = tableLevel, ErrCode = iError )
          if (iError.ne.0) then
            write(logUnit(1),*) 'ERROR in atl_modg_scheme_init: not able to read' &
              &            //' level for '
            write(logUnit(1),*) '      the level specific polynomial degree!'
            write(logUnit(1),*) 'STOPPING ...'
            call tem_abort()
          end if

          if (tableLevel.eq.currentLevel) then
            ! read the degree for the level
            call aot_get_val( L = conf, thandle = degTable, key = 'm', &
                          & val = degree, ErrCode = iError )
            if (iError.ne.0) then
              write(logUnit(1),*) 'ERROR in atl_modg_scheme_init: not able to ' &
                &            //' read degree for '
              write(logUnit(1),*) '      the level specific polynomial degree!'
              write(logUnit(1),*) '      Failing level: ', currentLevel
              write(logUnit(1),*) 'STOPPING ...'
              call tem_abort()
            end if
            foundLevel = .true.
          end if

          call aot_table_close(L = conf, thandle = degtable)

        end do

        ! Check if we found a degree for the current level of the mesh
        if (.not.foundlevel) then
          write(logUnit(1),*) 'ERROR in atl_modg_scheme_init:'
          write(logUnit(1),*) '      level not found for the level specific' &
            &         // '      polynomial degree!'
          write(logUnit(1),*) '      Failing level: ', currentLevel
          write(logUnit(1),*) 'STOPPING ...'
          call tem_abort()
        end if

      else predef

        if (btest(iError, aoterr_Fatal)) then
          write(logUnit(1),*) 'ERROR in atl_modg_scheme_init: ' &
            &            // 'could not read the predefined polynomial'
          write(logUnit(1),*) '      degrees for different levels!'
          write(logUnit(1),*) 'ABORTING...'
          call tem_abort()
        end if

        predeg = upper_to_lower(predefDegrees)
        predeg = adjustl(predeg)

        select case(trim(predeg))
        case ('fixedfact')
          ! This set of degrees increases the polynomial degrees by a fixed
          ! factor from level to level according to the following formula:
          ! m(iLevel) = nint(baseOrder*(factor**(maxLevel-iLevel)) - 1)

          ! Get the factor to increase the order level by level.
          ! Defaults to sqrt(2), as this maintains approximately the same
          ! time step restriction throughout the domain.
          call aot_get_val(L = conf, thandle = mtable, key = 'factor', &
            &              val = levelFact, default = sqrt(2.0_rk), &
            &              ErrCode = iError)
          if (btest(iError, aoterr_Fatal)) then
            write(logUnit(1),*) 'ERROR in atl_modg_scheme_init: ' &
              &            // 'could not read the factor for the'
            write(logUnit(1),*) '      predefined fixedfact set of degrees!'
            write(logUnit(1),*) 'ABORTING...'
            call tem_abort()
          end if

          ! Get the base order to use on the maximal level.
          ! Defaults to 1.
          call aot_get_val(L = conf, thandle = mtable, key = 'base_order', &
            &              val = base_order, default = 1.0_rk, &
            &              ErrCode = iError)
          if (btest(iError, aoterr_Fatal)) then
            write(logUnit(1),*) 'ERROR in atl_modg_scheme_init: ' &
              &            // 'could not read the base_order for the'
            write(logUnit(1),*) '      predefined fixedfact set of degrees!'
            write(logUnit(1),*) 'ABORTING...'
            call tem_abort()
          end if

          degree = nint( base_order*(levelFact**(maxLevel-currentLevel)) &
            &            - 1.0_rk )
        case default
          write(logUnit(1),*) 'ERROR in atl_modg_scheme_init: ' &
            &            // 'unknoen predefined polynomial degrees:'
          write(logUnit(1),*) trim(predefDegrees)
          write(logUnit(1),*) 'Currently supported predefined sets are:'
          write(logUnit(1),*) ' * fixedFact (use a fixed factor between levels)'
          write(logUnit(1),*) 'ABORTING...'
          call tem_abort()
        end select

      end if predef

      ! Close the m table
      call aot_table_close(L = conf, thandle = mtable)
    end if

  end subroutine atl_modg_scheme_load_polyDegree

end module atl_modg_scheme_module
