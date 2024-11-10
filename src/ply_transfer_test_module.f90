! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach for University of Siegen.
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

!> Test module for the routines in [[ply_transfer_module]]
module ply_transfer_test_module
  use env_module, only: rk

  use ply_dof_module, only: P_Space, Q_space
  use ply_transfer_module

  implicit none


contains


  ! ************************************************************************ !
  !> Test for ply_transfer_dofs_1d.
  subroutine ply_test_transfer_1d(success, lu)
    ! -------------------------------------------------------------------- !
    !> Indicator whether the check was successful.
    !!
    !! True if all tests pass successfully, otherwise false.
    logical, intent(out) :: success

    !> A logunit to write messages to.
    integer, intent(in) :: lu
    ! -------------------------------------------------------------------- !

    real(kind=rk) :: smallpoly(3)
    real(kind=rk) :: largepoly(5)
    real(kind=rk) :: tmp(22)
    real(kind=rk) :: eps
    logical :: test_ok

    success = .true.

    eps = epsilon(smallpoly(1))

    ! Test small to large
    write(lu,*) '1D check transfer small to large:'
    smallpoly = [16.0_rk, 8.0_rk, 4.0_rk]
    call ply_transfer_dofs_1d( indat     = smallpoly, &
      &                        indegree  = 2,         &
      &                        outdat    = largepoly, &
      &                        outdegree = 4          )

    test_ok = ( all(abs(largepoly(:3) - smallpoly) < eps) &
      &         .and. all(abs(largepoly(4:)) < eps)       )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallpoly
      write(lu,*) '   large:', largepoly
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test section (equally sized)
    tmp = -5.0_rk
    write(lu,*) '1D check transfer to equal sized but section:'
    call ply_transfer_dofs_1d( indat     = smallpoly, &
      &                        indegree  = 2,         &
      &                        outdat    = tmp(6:8),  &
      &                        outdegree = 2          )
    test_ok = ( all(abs(tmp(6:8) - smallpoly) < eps)  &
      &         .and. all(abs(tmp(:5)+5.0_rk) < eps)  &
      &         .and. all(abs(tmp(10:)+5.0_rk) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallpoly
      write(lu,*) '     tmp:', tmp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small
    write(lu,*) '1D check transfer large to small:'
    largepoly = [16.0_rk, -8.0_rk, 4.0_rk, -2.0_rk, 1.0_rk]
    call ply_transfer_dofs_1d( indat     = largepoly, &
      &                        indegree  = 4,         &
      &                        outdat    = smallpoly, &
      &                        outdegree = 2          )

    test_ok = ( all(abs(largepoly(:3) - smallpoly) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   large:', largepoly
      write(lu,*) '   small:', smallpoly
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    if (success) then
      write(lu,*) 'Successfully passed checks for ply_transfer_dofs_1d.'
    else
      write(lu,*) 'FAILED checks for ply_transfer_dofs_1d!'
    end if

  end subroutine ply_test_transfer_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Test for ply_transfer_dofs_2d.
  subroutine ply_test_transfer_2d(success, lu)
    ! -------------------------------------------------------------------- !
    !> Indicator whether the check was successful.
    !!
    !! True if all tests pass successfully, otherwise false.
    logical, intent(out) :: success

    !> A logunit to write messages to.
    integer, intent(in) :: lu
    ! -------------------------------------------------------------------- !

    real(kind=rk) :: smallq(9)
    real(kind=rk) :: largeq(25)
    real(kind=rk) :: smallp(6)
    real(kind=rk) :: largep(15)
    real(kind=rk) :: tmp(22)
    real(kind=rk) :: eps
    logical :: test_ok

    ! Note: 2D P-Poly is numbered as:
    !
    ! y-mode
    ! | x->  1  2  3  4  5
    ! v
    ! 1      1  2  4  7 11
    ! 2      3  5  8 12
    ! 3      6  9 13
    ! 4     10 14
    ! 5     15

    success = .true.

    eps = epsilon(smallq(1))

    ! Test small to large Q-Q
    write(lu,*) '2D Q check transfer small to large:'
    smallq = [ 16.0_rk,  8.0_rk, 4.0_rk, &
      &        27.0_rk,  9.0_rk, 3.0_rk, &
      &       125.0_rk, 25.0_rk, 5.0_rk  ]
    call ply_transfer_dofs_2d( indat     = smallq,  &
      &                        indegree  = 2,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = largeq,  &
      &                        outdegree = 4        )

    test_ok = ( all(abs(largeq(:3) - smallq(:3)) < eps) &
      &         .and. all(abs(largeq(4:5)) < eps)       )
    test_ok = test_ok                                           &
      &       .and. ( all(abs(largeq(6:8) - smallq(4:6)) < eps) &
      &               .and. all(abs(largeq(9:10)) < eps)        )
    test_ok = test_ok                                             &
      &       .and. ( all(abs(largeq(11:13) - smallq(7:9)) < eps) &
      &               .and. all(abs(largeq(14:)) < eps)           )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallq
      write(lu,*) '   large:', largeq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test section (equally sized) Q-Q
    tmp = -5.0_rk
    write(lu,*) '2D Q check transfer to equal sized but section:'
    call ply_transfer_dofs_2d( indat     = smallq,    &
      &                        indegree  = 2,         &
      &                        inspace   = Q_Space,   &
      &                        outspace  = Q_Space,   &
      &                        outdat    = tmp(6:14), &
      &                        outdegree = 2          )
    test_ok = ( all(abs(tmp(6:14) - smallq) < eps)    &
      &         .and. all(abs(tmp(:5)+5.0_rk) < eps)  &
      &         .and. all(abs(tmp(15:)+5.0_rk) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallq
      write(lu,*) '     tmp:', tmp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small Q-Q
    write(lu,*) '2D Q check transfer large to small:'
    largeq = [ 16.0_rk, -8.0_rk,  4.0_rk, -2.0_rk, 1.0_rk, &
      &        -8.0_rk,  4.0_rk, -2.0_rk,  1.0_rk, 0.5_rk, &
      &         4.0_rk, -2.0_rk,  1.0_rk, -0.5_rk, 0.3_rk, &
      &        -2.0_rk,  1.0_rk, -0.5_rk,  0.3_rk, 0.1_rk, &
      &         1.0_rk, -0.5_rk,  0.3_rk, -0.1_rk, 0.0_rk  ]
    call ply_transfer_dofs_2d( indat     = largeq, &
      &                        indegree  = 4,      &
      &                        inspace   = Q_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = smallq, &
      &                        outdegree = 2       )

    test_ok = ( all(abs(largeq(:3) - smallq(:3)) < eps) )
    test_ok = test_ok                                           &
      &       .and. ( all(abs(largeq(6:8) - smallq(4:6)) < eps) )
    test_ok = test_ok                                             &
      &       .and. ( all(abs(largeq(11:13) - smallq(7:9)) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   large:', largeq
      write(lu,*) '   small:', smallq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test small to large P-P
    write(lu,*) '2D P check transfer small to large:'
    smallp = [ 6.0_rk, 5.0_rk, 4.0_rk, &
      &        3.0_rk, 2.0_rk,         &
      &        1.0_rk                  ]
    call ply_transfer_dofs_2d( indat     = smallp,  &
      &                        indegree  = 2,       &
      &                        inspace   = P_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = largep,  &
      &                        outdegree = 4        )

    test_ok = ( all(abs(largep(:6) - smallp) < eps) &
      &         .and. all(abs(largep(7:)) < eps)    )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallp
      write(lu,*) '   large:', largep
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test section (equally sized) P-P
    tmp = -5.0_rk
    write(lu,*) '2D P check transfer to equal sized but section:'
    call ply_transfer_dofs_2d( indat     = smallp,    &
      &                        indegree  = 2,         &
      &                        inspace   = P_Space,   &
      &                        outspace  = P_Space,   &
      &                        outdat    = tmp(6:11), &
      &                        outdegree = 2          )
    test_ok = ( all(abs(tmp(6:11) - smallp) < eps)    &
      &         .and. all(abs(tmp(:5)+5.0_rk) < eps)  &
      &         .and. all(abs(tmp(12:)+5.0_rk) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallp
      write(lu,*) '     tmp:', tmp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small P-P
    write(lu,*) '2D P check transfer large to small:'
    largep = [ 15.0_rk,                                  &
      &        14.0_rk, 13.0_rk,                         &
      &        12.0_rk, 11.0_rk, 10.0_rk,                &
      &         9.0_rk,  8.0_rk,  7.0_rk, 6.0_rk,        &
      &         5.0_rk,  4.0_rk,  3.0_rk, 2.0_rk, 1.0_rk ]
    call ply_transfer_dofs_2d( indat     = largep,  &
      &                        indegree  = 4,       &
      &                        inspace   = P_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = smallp,  &
      &                        outdegree = 2        )

    test_ok = ( all(abs(largep(:6) - smallp(:6)) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   large:', largep
      write(lu,*) '   small:', smallp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test small to large P2Q
    write(lu,*) '2D P2Q check transfer small to large:'
    largeq = -1.0_rk
    smallp = [ 16.0_rk,                   & ! Diagonal starting x=0
      &         8.0_rk, 27.0_rk,          & ! Diagonal starting x=1
      &         4.0_rk,  9.0_rk, 125.0_rk ] ! Diagonal starting x=2
    call ply_transfer_dofs_2d( indat     = smallp,  &
      &                        indegree  = 2,       &
      &                        inspace   = P_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = largeq,  &
      &                        outdegree = 4        )

    test_ok = ( all(abs(largeq(:3) - smallp([1,2,4])) < eps) &
      &         .and. all(abs(largeq(4:5)) < eps)            )
    test_ok = test_ok                                             &
      &       .and. ( all(abs(largeq(6:7) - smallp([3,5])) < eps) &
      &               .and. all(abs(largeq(8:10)) < eps)          )
    test_ok = test_ok                                   &
      &       .and. ( abs(largeq(11) - smallp(6)) < eps &
      &               .and. all(abs(largeq(12:)) < eps) )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallp
      write(lu,*) '   large:', largeq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test equal max degree P2Q
    write(lu,*) '2D P2Q check transfer to equal sized:'
    call ply_transfer_dofs_2d( indat     = smallp,  &
      &                        indegree  = 2,       &
      &                        inspace   = P_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = smallq,  &
      &                        outdegree = 2        )
    test_ok = ( all(abs(smallq(:3) - smallp([1,2,4])) < eps) )
    test_ok = test_ok                                             &
      &       .and. ( all(abs(smallq(4:5) - smallp([3,5])) < eps) &
      &               .and. (abs(smallq(6)) < eps)                )
    test_ok = test_ok                                    &
      &       .and. ( (abs(smallq(7) - smallp(6)) < eps) &
      &               .and. all(abs(smallq(8:)) < eps)   )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   Ppoly:', smallp
      write(lu,*) '   Qpoly:', smallq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small P2Q
    write(lu,*) '2D P2Q check transfer large to small:'
    largep = [ 16.0_rk,                                   & ! diag start x=0
      &        -8.0_rk, -8.0_rk,                          & ! diag start x=1
      &         4.0_rk,  4.0_rk, 4.0_rk,                  & ! diag start x=2
      &        -2.0_rk, -2.0_rk, -2.0_rk, -2.0_rk,        & ! diag start x=3
      &         1.0_rk,  1.0_rk,  1.0_rk,  1.0_rk, 1.0_rk ] ! diag start x=4
    call ply_transfer_dofs_2d( indat     = largep,  &
      &                        indegree  = 4,       &
      &                        inspace   = P_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = smallq,  &
      &                        outdegree = 2        )

    test_ok = ( all(abs(smallq(:3) - largep([1,2,4])) < eps) )
    test_ok = test_ok                                               &
      &       .and. ( all(abs(smallq(4:6) - largep([3,5,8])) < eps) )
    test_ok = test_ok                                                &
      &       .and. ( all(abs(smallq(7:9) - largep([6,9,13])) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   largeP:', largep
      write(lu,*) '   smallQ:', smallq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test small to large Q2P
    write(lu,*) '2D Q2P check transfer small to large:'
    largep = -1.0_rk
    smallq = [ 16.0_rk,  8.0_rk, 4.0_rk, &
      &        27.0_rk,  9.0_rk, 3.0_rk, &
      &       125.0_rk, 25.0_rk, 5.0_rk  ]
    call ply_transfer_dofs_2d( indat     = smallq,  &
      &                        indegree  = 2,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = largep,  &
      &                        outdegree = 4        )

    test_ok = ( all(abs(smallq(:3) - largep([1,2,4])) < eps) &
      &         .and. all(abs(largep([7,11])) < eps)         )
    test_ok = test_ok                                               &
      &       .and. ( all(abs(smallq(4:6) - largep([3,5,8])) < eps) &
      &               .and. (abs(largep(12)) < eps)                 )
    test_ok = test_ok                                           &
      &       .and. ( all(abs(smallq(7:9) - largep([6,9,13])) < eps) &
      &               .and. all(abs(largep([10,14,15])) < eps)       )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallq
      write(lu,*) '   large:', largep
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test equal max degree Q2P
    tmp = -5.0_rk
    write(lu,*) '2D Q2P check transfer to equal sized:'
    call ply_transfer_dofs_2d( indat     = smallq,  &
      &                        indegree  = 2,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = smallp,  &
      &                        outdegree = 2        )
    test_ok = ( all(abs(smallq(:3) - smallp([1,2,4])) < eps) )
    test_ok = test_ok                                         &
      &       .and. all(abs(smallq(4:5) - smallp([3,5])) < eps)
    test_ok = test_ok                                &
      &       .and. (abs(smallq(7) - smallp(6)) < eps)
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   Qpoly:', smallq
      write(lu,*) '   Ppoly:', smallp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small Q2P
    write(lu,*) '2D Q2P check transfer large to small:'
    largeq = [ 16.0_rk, -8.0_rk,  4.0_rk, -2.0_rk, 1.0_rk, &
      &        -8.0_rk,  4.0_rk, -2.0_rk,  1.0_rk, 0.5_rk, &
      &         4.0_rk, -2.0_rk,  1.0_rk, -0.5_rk, 0.3_rk, &
      &        -2.0_rk,  1.0_rk, -0.5_rk,  0.3_rk, 0.1_rk, &
      &         1.0_rk, -0.5_rk,  0.3_rk, -0.1_rk, 0.0_rk  ]
    call ply_transfer_dofs_2d( indat     = largeq,  &
      &                        indegree  = 4,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = smallp,  &
      &                        outdegree = 2        )

    test_ok = ( all(abs(largeq(:3) - smallp([1,2,4])) < eps) )
    test_ok = test_ok                                             &
      &       .and. ( all(abs(largeq(6:7) - smallp([3,5])) < eps) )
    test_ok = test_ok                                 &
      &       .and. (abs(largeq(11) - smallp(6)) < eps)

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   largeQ:', largeq
      write(lu,*) '   smallP:', smallp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    if (success) then
      write(lu,*) 'Successfully passed checks for ply_transfer_dofs_2d.'
    else
      write(lu,*) 'FAILED checks for ply_transfer_dofs_2d!'
    end if


  end subroutine ply_test_transfer_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Test for ply_transfer_dofs_3d.
  subroutine ply_test_transfer_3d(success, lu)
    ! -------------------------------------------------------------------- !
    !> Indicator whether the check was successful.
    !!
    !! True if all tests pass successfully, otherwise false.
    logical, intent(out) :: success

    !> A logunit to write messages to.
    integer, intent(in) :: lu
    ! -------------------------------------------------------------------- !

    real(kind=rk) :: smallq(27)
    real(kind=rk) :: largeq(125)
    real(kind=rk) :: smallp(10)
    real(kind=rk) :: largep(35)
    real(kind=rk) :: tmp(200)
    real(kind=rk) :: eps
    real(kind=rk) :: maxmode
    integer :: iMode
    logical :: test_ok

    ! Note: 3D P-Poly is numbered as:
    !
    ! Z = 1                 !   Z = 2
    !                       !
    ! y-mode                !   y-mode
    ! | x->  1  2  3  4  5  !   | x->  1  2  3  4
    ! v                     !   v
    ! 1      1  2  5 11 21  !   1      4  8 15 26
    ! 2      3  6 12 22     !   2      9 16 27
    ! 3      7 13 23        !   3     17 28
    ! 4     14 24           !   4     29
    ! 5     25              !
    !                       !
    ! ------------------------------------------------
    !                 !                !
    ! Z = 3           !   Z = 4        !   Z = 5
    !                 !                !
    ! y-mode          !   y-mode       !   y-mode
    ! | x->  1  2  3  !   | x->  1  2  !   | x->  1
    ! v               !   v            !   v
    ! 1     10 18 30  !   1     20 33  !   1     35
    ! 2     19 31     !   2     34     !
    ! 3     32        !                !
    !                 !                !
    ! ------------------------------------------------

    success = .true.

    eps = epsilon(smallq(1))

    ! Test small to large Q-Q
    largeq = 3.14_rk
    write(lu,*) '3D Q check transfer small to large:'
    maxmode = real(size(smallq),kind=rk)
    do iMode=1,size(smallq)
      smallq(iMode) = maxmode - real(iMode-1, kind=rk)
    end do
    call ply_transfer_dofs_3d( indat     = smallq,  &
      &                        indegree  = 2,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = largeq,  &
      &                        outdegree = 4        )

    ! Z=1
    test_ok = ( all(abs(largeq(:3) - smallq(:3)) < eps) &
      &         .and. all(abs(largeq(4:5)) < eps)       )
    test_ok = test_ok                                           &
      &       .and. ( all(abs(largeq(6:8) - smallq(4:6)) < eps) &
      &               .and. all(abs(largeq(9:10)) < eps)        )
    test_ok = test_ok                                             &
      &       .and. ( all(abs(largeq(11:13) - smallq(7:9)) < eps) &
      &               .and. all(abs(largeq(14:25)) < eps)         )

    ! Z=2
    test_ok = test_ok                                               &
      &       .and. ( all(abs(largeq(26:28) - smallq(10:12)) < eps) &
      &               .and. all(abs(largeq(29:30)) < eps)           )
    test_ok = test_ok                                               &
      &       .and. ( all(abs(largeq(31:33) - smallq(13:15)) < eps) &
      &               .and. all(abs(largeq(34:35)) < eps)           )
    test_ok = test_ok                                               &
      &       .and. ( all(abs(largeq(36:38) - smallq(16:18)) < eps) &
      &               .and. all(abs(largeq(39:50)) < eps)           )

    ! Z=3
    test_ok = test_ok                                               &
      &       .and. ( all(abs(largeq(51:53) - smallq(19:21)) < eps) &
      &               .and. all(abs(largeq(54:55)) < eps)           )
    test_ok = test_ok                                               &
      &       .and. ( all(abs(largeq(56:58) - smallq(22:24)) < eps) &
      &               .and. all(abs(largeq(59:60)) < eps)           )
    test_ok = test_ok                                               &
      &       .and. ( all(abs(largeq(61:63) - smallq(25:27)) < eps) &
      &               .and. all(abs(largeq(64:)) < eps)             )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallq
      write(lu,*) '   large:', largeq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test section (equally sized) Q-Q
    tmp = -5.0_rk
    write(lu,*) '3D Q check transfer to equal sized but section:'
    call ply_transfer_dofs_3d( indat     = smallq,     &
      &                        indegree  = 2,          &
      &                        inspace   = Q_Space,    &
      &                        outspace  = Q_Space,    &
      &                        outdat    = tmp(11:37), &
      &                        outdegree = 2          )
    test_ok = ( all(abs(tmp(11:37) - smallq) < eps)   &
      &         .and. all(abs(tmp(:10)+5.0_rk) < eps) &
      &         .and. all(abs(tmp(38:)+5.0_rk) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallq
      write(lu,*) '     tmp:', tmp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small Q-Q
    write(lu,*) '3D Q check transfer large to small:'
    maxmode = real(size(largeq),kind=rk)
    do iMode=1,size(largeq)
      largeq(iMode) = maxmode - real(iMode-1, kind=rk)
    end do
    call ply_transfer_dofs_3d( indat     = largeq,  &
      &                        indegree  = 4,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = smallq,  &
      &                        outdegree = 2        )

    ! Z=1
    test_ok = all(abs(largeq(:3) - smallq(:3)) < eps)
    test_ok = test_ok .and. all(abs(largeq(6:8) - smallq(4:6)) < eps)
    test_ok = test_ok .and. all(abs(largeq(11:13) - smallq(7:9)) < eps)

    ! Z=2
    test_ok = test_ok .and. all(abs(largeq(26:28) - smallq(10:12)) < eps)
    test_ok = test_ok .and. all(abs(largeq(31:33) - smallq(13:15)) < eps)
    test_ok = test_ok .and. all(abs(largeq(36:38) - smallq(16:18)) < eps)

    ! Z=3
    test_ok = test_ok .and. all(abs(largeq(51:53) - smallq(19:21)) < eps)
    test_ok = test_ok .and. all(abs(largeq(56:58) - smallq(22:24)) < eps)
    test_ok = test_ok .and. all(abs(largeq(61:63) - smallq(25:27)) < eps)

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   large:', largeq
      write(lu,*) '   small:', smallq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test small to large P-P
    write(lu,*) '3D P check transfer small to large:'
    maxmode = real(size(smallp),kind=rk)
    do iMode=1,size(smallp)
      smallp(iMode) = maxmode - real(iMode-1, kind=rk)
    end do
    call ply_transfer_dofs_3d( indat     = smallp,  &
      &                        indegree  = 2,       &
      &                        inspace   = P_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = largep,  &
      &                        outdegree = 4        )

    test_ok = ( all(abs(largep(:10) - smallp) < eps) &
      &         .and. all(abs(largep(11:)) < eps)    )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallp
      write(lu,*) '   large:', largep
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test section (equally sized) P-P
    tmp = -5.0_rk
    write(lu,*) '3D P check transfer to equal sized but section:'
    call ply_transfer_dofs_3d( indat     = smallp,    &
      &                        indegree  = 2,         &
      &                        inspace   = P_Space,   &
      &                        outspace  = P_Space,   &
      &                        outdat    = tmp(6:15), &
      &                        outdegree = 2          )
    test_ok = ( all(abs(tmp(6:15) - smallp) < eps)    &
      &         .and. all(abs(tmp(:5)+5.0_rk) < eps)  &
      &         .and. all(abs(tmp(16:)+5.0_rk) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallp
      write(lu,*) '     tmp:', tmp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small P-P
    write(lu,*) '3D P check transfer large to small:'
    maxmode = real(size(largep),kind=rk)
    do iMode=1,size(largep)
      largep(iMode) = maxmode - real(iMode-1, kind=rk)
    end do
    call ply_transfer_dofs_3d( indat     = largep,  &
      &                        indegree  = 4,       &
      &                        inspace   = P_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = smallp,  &
      &                        outdegree = 2        )

    test_ok = ( all(abs(largep(:10) - smallp) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   large:', largep
      write(lu,*) '   small:', smallp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test small to large P2Q
    write(lu,*) '3D P2Q check transfer small to large:'
    largeq = -1.0_rk
    maxmode = real(size(smallp),kind=rk)
    do iMode=1,size(smallp)
      smallp(iMode) = maxmode - real(iMode-1, kind=rk)
    end do
    call ply_transfer_dofs_3d( indat     = smallp,  &
      &                        indegree  = 2,       &
      &                        inspace   = P_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = largeq,  &
      &                        outdegree = 4        )

    ! Z=1
    test_ok = ( all(abs(largeq(:3) - smallp([1,2,5])) < eps) &
      &         .and. all(abs(largeq(4:5)) < eps)            )
    test_ok = test_ok                                             &
      &       .and. ( all(abs(largeq(6:7) - smallp([3,6])) < eps) &
      &               .and. all(abs(largeq(8:10)) < eps)          )
    test_ok = test_ok                                     &
      &       .and. ( abs(largeq(11) - smallp(7)) < eps   &
      &               .and. all(abs(largeq(12:25)) < eps) )

    ! Z=2
    test_ok = test_ok                                               &
      &       .and. ( all(abs(largeq(26:27) - smallp([4,8])) < eps) &
      &               .and. all(abs(largeq(28:30)) < eps)          )
    test_ok = test_ok                                     &
      &       .and. ( abs(largeq(31) - smallp(9)) < eps   &
      &               .and. all(abs(largeq(32:50)) < eps) )

    ! Z=3
    test_ok = test_ok                                      &
      &       .and. ( (abs(largeq(51) - smallp(10)) < eps) &
      &               .and. all(abs(largeq(52:)) < eps)    )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallp
      write(lu,*) '   large:', largeq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test equal max degree P2Q
    write(lu,*) '3D P2Q check transfer to equal sized:'
    call ply_transfer_dofs_3d( indat     = smallp,  &
      &                        indegree  = 2,       &
      &                        inspace   = P_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = smallq,  &
      &                        outdegree = 2        )

    ! Z=1
    test_ok = ( all(abs(smallq(:3) - smallp([1,2,5])) < eps) )
    test_ok = test_ok                                             &
      &       .and. ( all(abs(smallq(4:5) - smallp([3,6])) < eps) &
      &               .and. (abs(smallq(6)) < eps)                )
    test_ok = test_ok                                    &
      &       .and. ( (abs(smallq(7) - smallp(7)) < eps) &
      &               .and. all(abs(smallq(8:9)) < eps)  )

    ! Z=2
    test_ok = test_ok                                               &
      &       .and. ( all(abs(smallq(10:11) - smallp([4,8])) < eps) &
      &               .and. (abs(smallq(12)) < eps)                 )
    test_ok = test_ok                                     &
      &       .and. ( (abs(smallq(13) - smallp(9)) < eps) &
      &               .and. all(abs(smallq(14:18)) < eps) )

    ! Z=3
    test_ok = test_ok                                      &
      &       .and. ( (abs(smallq(19) - smallp(10)) < eps) &
      &               .and. all(abs(smallq(20:)) < eps)    )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   Ppoly:', smallp
      write(lu,*) '   Qpoly:', smallq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small P2Q
    write(lu,*) '3D P2Q check transfer large to small:'
    maxmode = real(size(largep),kind=rk)
    do iMode=1,size(largep)
      largep(iMode) = maxmode - real(iMode-1, kind=rk)
    end do
    call ply_transfer_dofs_3d( indat     = largep,  &
      &                        indegree  = 4,       &
      &                        inspace   = P_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = smallq,  &
      &                        outdegree = 2        )

    ! Z=1
    test_ok = ( all(abs(smallq(:3) - largep([1,2,5])) < eps) )
    test_ok = test_ok                                                &
      &       .and. ( all(abs(smallq(4:6) - largep([3,6,12])) < eps) )
    test_ok = test_ok                                                 &
      &       .and. ( all(abs(smallq(7:9) - largep([7,13,23])) < eps) )

    ! Z=2
    test_ok = test_ok                                                  &
      &       .and. ( all(abs(smallq(10:12) - largep([4,8,15])) < eps) )
    test_ok = test_ok                                                  &
      &       .and. ( all(abs(smallq(13:15) - largep([9,16,27])) < eps) )
    test_ok = test_ok                                                 &
      &       .and. ( all(abs(smallq(16:17) - largep([17,28])) < eps) &
      &               .and. (abs(smallq(18)) < eps)                   )

    ! Z=3
    test_ok = test_ok                                                    &
      &       .and. ( all(abs(smallq(19:21) - largep([10,18,30])) < eps) )
    test_ok = test_ok                                                 &
      &       .and. ( all(abs(smallq(22:23) - largep([19,31])) < eps) &
      &               .and. (abs(smallq(24)) < eps)                   )
    test_ok = test_ok                                      &
      &       .and. ( (abs(smallq(25) - largep(32)) < eps) &
      &               .and. all(abs(smallq(26:)) < eps)    )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   largeP:', largep
      write(lu,*) '   smallQ:', smallq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test small to large Q2P
    write(lu,*) '3D Q2P check transfer small to large:'
    largep = -1.0_rk
    maxmode = real(size(smallq),kind=rk)
    do iMode=1,size(smallq)
      smallq(iMode) = maxmode - real(iMode-1, kind=rk)
    end do
    call ply_transfer_dofs_3d( indat     = smallq,  &
      &                        indegree  = 2,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = largep,  &
      &                        outdegree = 4        )

    ! Z=1
    test_ok = ( all(abs(smallq(:3) - largep([1,2,5])) < eps) &
      &         .and. all(abs(largep([11,21])) < eps)        )
    test_ok = test_ok                                                &
      &       .and. ( all(abs(smallq(4:6) - largep([3,6,12])) < eps) &
      &               .and. (abs(largep(22)) < eps)                 )
    test_ok = test_ok                                           &
      &       .and. ( all(abs(smallq(7:9) - largep([7,13,23])) < eps) &
      &               .and. all(abs(largep([14,24,25])) < eps)        )

    ! Z=2
    test_ok = test_ok                                                  &
      &       .and. ( all(abs(smallq(10:12) - largep([4,8,15])) < eps) &
      &               .and. (abs(largep(26)) < eps)                    )
    test_ok = test_ok .and. all(abs(smallq(13:15) - largep([9,16,27])) < eps)
    test_ok = test_ok .and. all(abs(smallq(16:17) - largep([17,28])) < eps)
    test_ok = test_ok .and. (abs(largep(29)) < eps)

    ! Z=3
    test_ok = test_ok .and. all(abs(smallq(19:21) - largep([10,18,30])) < eps)
    test_ok = test_ok .and. all(abs(smallq(22:23) - largep([19,31])) < eps)
    test_ok = test_ok .and. (abs(smallq(25) - largep(32)) < eps)
    test_ok = test_ok .and. all(abs(largep([20,33,34,35])) < eps)

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallq
      write(lu,*) '   large:', largep
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test equal max degree Q2P
    tmp = -5.0_rk
    write(lu,*) '3D Q2P check transfer to equal sized:'
    call ply_transfer_dofs_3d( indat     = smallq,  &
      &                        indegree  = 2,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = smallp,  &
      &                        outdegree = 2        )

    ! Z=1
    test_ok = ( all(abs(smallq(:3) - smallp([1,2,5])) < eps) )
    test_ok = test_ok .and. all(abs(smallq(4:5) - smallp([3,6])) < eps)
    test_ok = test_ok .and. (abs(smallq(7) - smallp(7)) < eps)

    ! Z=2
    test_ok = test_ok .and. all(abs(smallq(10:11) - smallp([4,8])) < eps)
    test_ok = test_ok .and. (abs(smallq(13) - smallp(9)) < eps)

    ! Z=3
    test_ok = test_ok .and. (abs(smallq(19) - smallp(10)) < eps)

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   Qpoly:', smallq
      write(lu,*) '   Ppoly:', smallp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small Q2P
    write(lu,*) '3D Q2P check transfer large to small:'
    maxmode = real(size(largeq),kind=rk)
    do iMode=1,size(largeq)
      largeq(iMode) = maxmode - real(iMode-1, kind=rk)
    end do
    call ply_transfer_dofs_3d( indat     = largeq,  &
      &                        indegree  = 4,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = P_Space, &
      &                        outdat    = smallp,  &
      &                        outdegree = 2        )

    ! Z=1
    test_ok = ( all(abs(largeq(:3) - smallp([1,2,5])) < eps) )
    test_ok = test_ok .and. all(abs(largeq(6:7) - smallp([3,6])) < eps)
    test_ok = test_ok .and. (abs(largeq(11) - smallp(7)) < eps)

    ! Z=2
    test_ok = test_ok .and. all(abs(largeq(26:27) - smallp([4,8])) < eps)
    test_ok = test_ok .and. (abs(largeq(31) - smallp(9)) < eps)

    ! Z=3
    test_ok = test_ok .and. (abs(largeq(51) - smallp(10)) < eps)

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   largeQ:', largeq
      write(lu,*) '   smallP:', smallp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    if (success) then
      write(lu,*) 'Successfully passed checks for ply_transfer_dofs_3d.'
    else
      write(lu,*) 'FAILED checks for ply_transfer_dofs_3d!'
    end if


  end subroutine ply_test_transfer_3d
  ! ************************************************************************ !


end module ply_transfer_test_module
