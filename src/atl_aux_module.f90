! Copyright (c) 2012-2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013, 2015-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Some auxilary routines.
module atl_aux_module
  use env_module,           only: rk, labelLen

  ! Aotus modules
  use aotus_module,         only: flu_State
  use aot_table_module,     only: aot_table_open,  &
    &                             aot_table_close, &
    &                             aot_get_val
  use aot_table_ops_module, only: aot_table_length

  ! Treelm modules
  use tem_logging_module,   only: logUnit, llerror, lldebug
  use tem_tools_module,     only: upper_to_lower
  use tem_aux_module,       only: tem_print_execInfo, utc_date_string
  use tem_logging_module,   only: logUnit

  implicit none

  private

  public :: atl_banner
  public :: atl_bubbleSortArray

contains

  ! ************************************************************************ !
  !> Prominently let the user now, what he actually is running right now.
  subroutine atl_banner(version_string)
    ! -------------------------------------------------------------------- !
    character(len=*), intent(in) :: version_string
    ! -------------------------------------------------------------------- !
    character(len=26) :: dat_string
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) "                                               "
    write(logUnit(1),*) "          ____________     ______              "
    write(logUnit(1),*) "          ___    |_  /________  /____________  "
    write(logUnit(1),*) "          __  /| |  __/  _ \_  /_  _ \_  ___/  "
    write(logUnit(1),*) "          _  ___ / /_ /  __/  / /  __/(__  )   "
    write(logUnit(1),*) "          /_/  |_\__/ \___//_/  \___//___"      &
      &            // trim(version_string)
    write(logUnit(1),*) "                                             "
    write(logUnit(1),*) " (C) 2012 German Research School for Simulation" &
      &                 // " Sciences"
    write(logUnit(1),*) " (C) 2013,2014-2020 University of Siegen"
    write(logUnit(1),*) "                                             "
    ! Write the information about the executable, gathered at build time to
    ! the screen.
    call tem_print_execInfo()
    write(logUnit(1),*) "                                             "
    dat_string = utc_date_string()
    write(logUnit(1),*) "Run at: " // dat_string
    write(logUnit(1),*) "                                             "
    write(logUnit(1),*) "                                             "

  end subroutine atl_banner
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Bubble sorting of array of real numbers of size n
  ! returns array A sorted in Ascending order
  subroutine atl_bubbleSortArray( A, n )
    ! -------------------------------------------------------------------- !
    !> The number of elements in array A
    integer,       intent(in)    :: n
    !> The array to sort
    real(kind=rk), intent(inout) :: A(1:n)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: temp
    integer i, j
    ! -------------------------------------------------------------------- !

    do i = 1, n
      do j = n, i + 1, -1
        if( A(j-1) > A(j) ) then
          temp = A(j-1)
          A(j-1) = A(j)
          A(j) = temp
        end if
      end do
    end do

  end subroutine atl_bubbleSortArray
  ! ************************************************************************ !


end module atl_aux_module
