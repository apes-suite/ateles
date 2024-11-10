! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Module to compute the exact solution for the Riemann Problem for
!! the Euler equations.
module atl_exact_riemann_euler_module
  use env_module, only: rk

  implicit none

  private

  type atl_ere_solState1D_type
    !> Auxilary values describing the gas (given by the isentropic coefficient).
    real(kind=rk) :: G(9)

    !> Isentropic coefficient
    real(kind=rk) :: gam

    !> Abort criterion for the Newton-Raphson scheme.
    real(kind=rk) :: tolerance

    !> Maximum number of iterations for the Newton-Raphson scheme.
    integer :: nMaxIter

    !> Number of actually made Newton-Raphson iterations.
    integer :: nIter

    !> Flag if this is a critical state combination, that would result in a
    !! vacuum.
    logical :: critical

    !> Primitive variables left of discontinuity
    !!
    !! - 1: density
    !! - 2: velocity
    !! - 3: pressure
    real(kind=rk) :: prim_left(3)

    !> Primitive variables right of discontinuity
    !!
    !! - 1: density
    !! - 2: velocity
    !! - 3: pressure
    real(kind=rk) :: prim_right(3)

    real(kind=rk) :: prim_rare2cont(3)
    real(kind=rk) :: prim_cont2shock(3)
    real(kind=rk) :: prim_rare(3)

    real(kind=rk) :: p_mean
    real(kind=rk) :: u_mean

    real(kind=rk) :: rare_c
    real(kind=rk) :: rare_c_fac

    !> Speed of sound left
    real(kind=rk) :: cl

    !> Speed of sound right
    real(kind=rk) :: cr

    !> Characteristic velocities
    !!
    !! Ordered from left most to right most.
    !! Using 4 characteristic velocities to describe the solution of the
    !! Riemann problem:
    !! - 2 for the rarefication
    !! - 1 contact discontinuity
    !! - 1 shock
    real(kind=rk) :: charvel(4)

    !> Kind of state between the characteristic velocities.
    !!
    !! Ordered from left to right.
    !! There are three fields between the characteristics, the states outside
    !! the characteristic velocities do not depend on the characteristics.
    !! For each state the kind is one of the following:
    !! - rarefication
    !! - rare2cont
    !! - cont2shock
    integer :: charfield(3)

  end type atl_ere_solState1D_type


  public :: atl_ere_solState1D_type
  public :: atl_ere_init
  public :: atl_ere_init_consts
  public :: atl_ere_sample
  public :: atl_ere_eval_onEdge
  public :: atl_ere_dump_solState


  integer, parameter :: density = 1
  integer, parameter :: velocity = 2
  integer, parameter :: pressure = 3

  integer, parameter :: rarefication = 1
  integer, parameter :: rare2cont = 2
  integer, parameter :: cont2shock = 3


contains


  !> Compute the exact solution of the Riemann problem for the Euler equation.
  !!
  !! For the given state in primitive variables left and right, the solution is
  !! computed and stored in me.
  !! A Newton-Raphson iteration is used to find the solution in this
  !! non-linear system.
  !! Necessary parameters besides left and right state are the isentropic
  !! coefficient gam, and optionally a target tolerance to reach in the
  !! iterative scheme and the maximal number of iterations to perform.
  subroutine atl_ere_init( me, prim_left, prim_right, gam, tolerance, nMaxIter )
    !> Description of the Riemann solution.
    !!
    !! If the input results in a critical state, me%critical will be set to
    !! true. nIter will contain the actual number of Newton-Raphson iterations,
    !! that were performed.
    type(atl_ere_solState1D_type), intent(out) :: me

    !> Primitive variables left of discontinuity
    !!
    !! - 1: density
    !! - 2: velocity
    !! - 3: pressure
    real(kind=rk), intent(in) :: prim_left(3)

    !> Primitive variables right of discontinuity
    real(kind=rk), intent(in) :: prim_right(3)

    !> Isentropic expansion coefficient of the gas
    real(kind=rk), intent(in) :: gam

    !> Tolerated variation in the Newton-Raphson iterations, default: 8 epsilon
    real(kind=rk), intent(in), optional :: tolerance

    !> Maximal number of iterations to do, default: 100000
    integer, intent(in), optional :: nMaxIter

    real(kind=rk) :: du
    real(kind=rk) :: u
    real(kind=rk) :: cl, cr
    real(kind=rk) :: p, p0
    real(kind=rk) :: cha
    real(kind=rk) :: fl, fld
    real(kind=rk) :: fr, frd
    real(kind=rk) :: cm
    real(kind=rk) :: pm
    integer :: kk

    call atl_ere_init_consts( me        = me,        &
      &                       gam       = gam,       &
      &                       tolerance = tolerance, &
      &                       nMaxIter  = nMaxIter   )

    ! Velocity difference over discontinuity.
    du = prim_right(velocity) - prim_left(velocity)

    ! Speed of Sound left.
    cl = sqrt(gam * prim_left(pressure)/prim_left(density))

    ! Speed of Sound right.
    cr = sqrt(gam * prim_right(pressure)/prim_right(density))

    me%critical = (me%G(4) * (cl+cr) - du <= 0.0_rk)

    if (.not. me%critical) then
      ! Only start the Newton-Raphson scheme, if the input does not result
      ! in a critical state.
      call nr_start(p, prim_left(density),  prim_right(density),  &
        &              prim_left(velocity), prim_right(velocity), &
        &              prim_left(pressure), prim_right(pressure), &
        &              cl, cr, me)

      do kk=1,me%nMaxIter

        p0 = p

        call nr_1side(fl, fld, p, prim_left(density), &
          &                       prim_left(pressure), &
          &           cl, me)
        call nr_1side(fr, frd, p, prim_right(density), &
          &                       prim_right(pressure), &
          &           cr, me)

        p = max(p - (fl+fr+du) / (fld+frd), me%tolerance)

        ! Compute relative change in pressure
        cha = 2.0_rk * abs( (p-p0) / (p+p0) )

        ! Leave iterations, when relative change is sufficiently small
        if (cha <= me%tolerance) EXIT

      end do

      u = 0.5_rk*(prim_left(velocity)+prim_right(velocity) &
        &         + fr-fl)

      me%nIter = kk
      me%prim_left = prim_left
      me%prim_right = prim_right
      me%cl = cl
      me%cr = cr
      me%p_mean = p
      me%u_mean = u

      if (me%p_mean <= me%prim_left(pressure)) then
        ! Shock moves to the right.
        me%charfield(1) = rarefication
        me%charfield(2) = rare2cont
        me%charfield(3) = cont2shock

        cm = me%cl*(me%p_mean/me%prim_left(pressure))**me%G(1)
        pm = me%p_mean/me%prim_right(pressure)

        ! Left end of rarefication
        me%charvel(1) = me%prim_left(velocity) - me%cl

        ! Right end of rarefication
        me%charvel(2) = me%u_mean - cm

        ! Contact discontinuity
        me%charvel(3) = me%u_mean

        ! Shock
        me%charvel(4) = me%prim_right(velocity) &
          &             + me%cr*sqrt(me%G(2)*pm + me%G(1))

        ! Compute the state in the different areas of the solution
        me%prim_rare2cont(density) = me%prim_left(density) &
          &                          * (me%p_mean &
          &                             / me%prim_left(pressure))**me%G(8)
        me%prim_rare2cont(velocity) = me%u_mean
        me%prim_rare2cont(pressure) = me%p_mean

        me%prim_cont2shock(density) = me%prim_right(density) * (pm+me%G(6)) &
          &                           / (pm*me%G(6) + 1.0_rk)
        me%prim_cont2shock(velocity) = me%u_mean
        me%prim_cont2shock(pressure) = me%p_mean

        me%prim_rare(velocity) = me%cl + me%G(7)*me%prim_left(velocity)
        me%rare_c = me%cl
        me%rare_c_fac = me%G(5)
        me%prim_rare(density) = me%prim_left(density)
        me%prim_rare(pressure) = me%prim_left(pressure)

      else
        ! Shock moves to the left.
        me%charfield(1) = cont2shock
        me%charfield(2) = rare2cont
        me%charfield(3) = rarefication

        cm = me%cr*(me%p_mean/me%prim_right(pressure))**me%G(1)
        pm = me%p_mean/me%prim_left(pressure)

        ! Shock
        me%charvel(1) = me%prim_left(velocity) &
          &             - me%cl*sqrt(me%G(2)*pm + me%G(1))

        ! Contact discontinuity
        me%charvel(2) = me%u_mean

        ! Left end of rarefication
        me%charvel(3) = me%u_mean + cm

        ! Right end of rarefication
        me%charvel(4) = me%prim_right(velocity) + me%cr

        ! Compute the state in the different areas of the solution
        me%prim_rare2cont(density) = me%prim_right(density) &
          &                          * (me%p_mean &
          &                             / me%prim_right(pressure))**me%G(8)
        me%prim_rare2cont(velocity) = me%u_mean
        me%prim_rare2cont(pressure) = me%p_mean

        me%prim_cont2shock(density) = me%prim_left(density) * (pm+me%G(6)) &
          &                           / (pm*me%G(6) + 1.0_rk)
        me%prim_cont2shock(velocity) = me%u_mean
        me%prim_cont2shock(pressure) = me%p_mean

        me%prim_rare(velocity) = me%G(7)*me%prim_right(velocity) - me%cr
        me%rare_c = me%cr
        me%rare_c_fac = -me%G(5)
        me%prim_rare(density) = me%prim_right(density)
        me%prim_rare(pressure) = me%prim_right(pressure)
      end if
    end if

  end subroutine atl_ere_init
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Set the general constants for the exact riemann solver.
  subroutine atl_ere_init_consts( me, gam, tolerance, nMaxIter )
    !> Description of the Riemann solution to set the constants in
    type(atl_ere_solState1D_type), intent(out) :: me

    !> Isentropic expansion coefficient of the gas
    real(kind=rk), intent(in) :: gam

    !> Tolerated variation in the Newton-Raphson iterations, default: 8 epsilon
    real(kind=rk), intent(in), optional :: tolerance

    !> Maximal number of iterations to do, default: 100000
    integer, intent(in), optional :: nMaxIter

    if (present(tolerance)) then
      me%tolerance = tolerance
    else
      me%tolerance = 8.0_rk * epsilon(me%tolerance)
    end if

    if (present(nMaxIter)) then
      me%nMaxIter = nMaxIter
    else
      me%nMaxIter = 100000
    end if

    me%gam = gam
    me%G = [ (gam - 1.0_rk) / (2.0_rk * gam), &
             (gam + 1.0_rk) / (2.0_rk * gam), &
             2.0_rk * gam / (gam - 1.0_rk),   &
             2.0_rk / (gam - 1.0_rk),         &
             2.0_rk / (gam + 1.0_rk),         &
             (gam - 1.0_rk) / (gam + 1.0_rk), &
             0.5_rk * (gam-1.0_rk),           &
             1.0_rk / gam,                    &
             gam - 1.0_rk                     ]

  end subroutine atl_ere_init_consts
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Evaluate the state on the edge using the exact riemann solver
  elemental subroutine atl_ere_eval_onEdge( me,                               &
    &                       rho_left,   vn_left,  v1_left,  v2_left,  p_left, &
    &                       rho_right, vn_right, v1_right, v2_right, p_right, &
    &                             rho,       vn,       v1,       v2,       p  )
    !> Description of the general constants for the Riemann solver.
    type(atl_ere_solState1D_type), intent(in) :: me
    real(kind=rk),intent(in) :: rho_left, vn_left, v1_left, v2_left, p_left
    real(kind=rk),intent(in) :: rho_right, vn_right, v1_right, v2_right, p_right
    real(kind=rk),intent(out) :: rho, vn, v1, v2, p

    real(kind=rk) :: du
    real(kind=rk) :: u
    real(kind=rk) :: cl, cr
    real(kind=rk) :: p0
    real(kind=rk) :: cha
    real(kind=rk) :: fl, fld
    real(kind=rk) :: fr, frd
    real(kind=rk) :: pratio, cratio
    integer :: kk

    ! Set a bad state to return, if there is a critical state.
    rho = 0.0_rk
    p   = 0.0_rk

    du = vn_right - vn_left
    cl = sqrt(me%gam * p_left / rho_left)
    cr = sqrt(me%gam * p_right / rho_right)

    if (me%g(4) * (cl+cr) - du > 0.0_rk) then
      ! Only proceed, if we do not run into a critical state.
      call nr_start(p, rho_left,  rho_right, &
        &               vn_left,   vn_right, &
        &                p_left,    p_right, &
        &                    cl,         cr, me)

      do kk=1,me%nMaxIter

        p0 = p

        call nr_1side(fl, fld, p, rho_left, &
          &                         p_left, &
          &           cl, me)
        call nr_1side(fr, frd, p, rho_right, &
          &                         p_right, &
          &           cr, me)

        p = max(p - (fl+fr+du) / (fld+frd), me%tolerance)

        ! Compute relative change in pressure
        cha = 2.0_rk * abs( (p-p0) / (p+p0) )

        ! Leave iterations, when relative change is sufficiently small
        if (cha <= me%tolerance) EXIT

      end do

      u = 0.5_rk*(vn_left + vn_right + fr - fl)

      if (u >= 0.0_rk) then

        ! Contact discontinuity is right of the edge
        if (p <= p_left) then
          ! Shock moves to the right
          if ((vn_left - cl) >= 0.0_rk) then
            ! All characteristics move to the right, simply use the
            ! left state on the edge.
            rho = rho_left
            vn  = vn_left
            v1  = v1_left
            v2  = v2_left
            p   = p_left
          else
            pratio = p/p_left
            if ((u - cl*pratio**me%G(1)) < 0.0_rk) then
              ! The edge is between rarefication and contact.
              rho = rho_left * pratio**me%G(8)
              vn  = u
              !>@todo: clarify what to do with v1 and v2
              v1 = v1_left
              v2 = v2_left
            else
              ! The edge is within the rarefication.
              vn  = me%G(5)*(cl + me%G(7)*vn_left)
              cratio = vn/cl
              rho = rho_left * cratio**me%G(4)
              p   = p_left * cratio**me%G(3)
              !>@todo: clarify what to do with v1 and v2
              v1 = v1_left
              v2 = v2_left
            end if
          end if

        else
          ! Shock moves to the left
          pratio = p / p_left
          if ((vn_left - cl*sqrt(me%G(2)*pratio + me%G(1))) >= 0.0_rk) then
            ! Edge is left of the shock, just use the left state.
            rho = rho_left
            vn  = vn_left
            v1  = v1_left
            v2  = v2_left
            p   = p_left
          else
            ! Edge is between shock and contact
            vn = u
            rho = rho_left*(pratio+me%G(6)) / (pratio*me%G(6) + 1.0_rk)
            !>@todo: clarify what to do with v1 and v2
            v1 = v1_right
            v2 = v2_right
          end if
        end if

      else

        ! Contact discontinuity is left of the edge
        if (p <= p_right) then
          ! Shock moves to the left
          if ((vn_right + cr) <= 0.0_rk) then
            ! All characteristics move to the left, simply use the right
            ! state on the edge.
            rho = rho_right
            vn  = vn_right
            v1  = v1_right
            v2  = v2_right
            p   = p_right
          else
            pratio = p / p_right
            if ((u + cr*pratio**me%G(1)) > 0.0_rk) then
              ! The edge is between rarefication and contact.
              rho = rho_right * pratio**me%G(8)
              vn = u
              !>@todo: clarify what to do with v1 and v2
              v1 = v1_right
              v2 = v2_right
            else
              ! The edge is within the rarefication.
              vn = me%G(5)*(me%G(7)*vn_right - cr)
              cratio = -vn/cr
              rho = rho_right * cratio**me%G(4)
              p   = p_right * cratio**me%G(3)
              !>@todo: clarify what to do with v1 and v2
              v1 = v1_right
              v2 = v2_right
            end if
          end if
        else
          ! Shock moves to the right
          pratio = p / p_right
          if ((vn_right + cr*sqrt(me%G(2)*pratio + me%G(1))) <= 0.0_rk) then
            ! Edge is right of the shock, use right state
            rho = rho_right
            vn  = vn_right
            v1  = v1_right
            v2  = v2_right
            p   = p_right
          else
            ! Edge is between shock and contact
            vn = u
            rho = rho_right*(pratio+me%G(6)) / (pratio*me%G(6) + 1.0_rk)
            !>@todo: clarify what to do with v1 and v2
            v1 = v1_left
            v2 = v2_left
          end if
        end if

      end if

    end if
  end subroutine atl_ere_eval_onEdge
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Evaluate the solution to the 1D Riemann problem for a given sample point s.
  !!
  !! The sample s has to be a velocity, which is compared to the various
  !! characteristics in the solution of the Riemann problem.
  !! Depending on the location of s between those characteristics, the
  !! appropriate state prim_state is returned.
  elemental subroutine atl_ere_sample( me, s, rho, v, p )
    !> Description of the solution.
    type(atl_ere_solState1D_type), intent(in) :: me

    !> Velocity along which to find the state.
    real(kind=rk), intent(in) :: s

    !> Density at s.
    real(kind=rk), intent(out) :: rho

    !> Velocity at s.
    real(kind=rk), intent(out) :: v

    !> Pressure at s.
    real(kind=rk), intent(out) :: p

    integer :: state_field
    real(kind=rk) :: c

    ! Decide which area of the solution s is located in.
    if (s <= me%charvel(1)) then
      ! Left of leftmost characteristic, just return the left state.
      rho = me%prim_left(density)
      v   = me%prim_left(velocity)
      p   = me%prim_left(pressure)
      state_field = 0
    else if (s <= me%charvel(2)) then
      state_field = me%charfield(1)
    else if (s <= me%charvel(3)) then
      state_field = me%charfield(2)
    else if (s < me%charvel(4)) then
      state_field = me%charfield(3)
    else
      ! Right of rightmost characteristic, just return the right state.
      rho = me%prim_right(density)
      v   = me%prim_right(velocity)
      p   = me%prim_right(pressure)
      state_field = 0
    end if

    ! Set the solution according to its location in the fields of the solution,
    ! if s is located in between the characteristics.
    select case(state_field)
    case(rare2cont)
      rho = me%prim_rare2cont(density)
      v   = me%prim_rare2cont(velocity)
      p   = me%prim_rare2cont(pressure)
    case(cont2shock)
      rho = me%prim_cont2shock(density)
      v   = me%prim_cont2shock(velocity)
      p   = me%prim_cont2shock(pressure)
    case(rarefication)
      c = me%rare_c_fac * (me%prim_rare(velocity) - me%G(7)*s)
      v   = me%G(5) * (me%prim_rare(velocity) + s)
      rho = me%prim_rare(density) * (c/me%rare_c)**me%G(4)
      p   = me%prim_rare(pressure) * (c/me%rare_c)**me%G(3)
    end select

  end subroutine atl_ere_sample
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  subroutine atl_ere_dump_solState( me, outunit )
    type(atl_ere_solState1D_type), intent(in) :: me
    integer, intent(in), optional :: outunit

    if (present(outunit)) then
      write(outunit,*) 'Riemann solution:'
      write(outunit,*) 'tolerance  = ', me%tolerance
      write(outunit,*) 'nMaxIter   = ', me%nMaxIter
      write(outunit,*) 'nIter      = ', me%nIter
      write(outunit,*) 'critical   = ', me%critical
      write(outunit,*) 'left       = ', me%prim_left
      write(outunit,*) 'right      = ', me%prim_right
      write(outunit,*) 'rare2cont  = ', me%prim_rare2cont
      write(outunit,*) 'cont2shock = ', me%prim_cont2shock
      write(outunit,*) 'rare       = ', me%prim_rare
      write(outunit,*) 'p_mean     = ', me%p_mean
      write(outunit,*) 'u_mean     = ', me%u_mean
      write(outunit,*) 'rare_c     = ', me%rare_c
      write(outunit,*) 'rare_c_fac = ', me%rare_c_fac
      write(outunit,*) 'cl         = ', me%cl
      write(outunit,*) 'cr         = ', me%cr
      write(outunit,*) 'charvel    = ', me%charvel
      write(outunit,*) 'charfield  = ', me%charfield
    else
      write(*,*) 'Riemann solution:'
      write(*,*) 'tolerance  = ', me%tolerance
      write(*,*) 'nMaxIter   = ', me%nMaxIter
      write(*,*) 'nIter      = ', me%nIter
      write(*,*) 'critical   = ', me%critical
      write(*,*) 'left       = ', me%prim_left
      write(*,*) 'right      = ', me%prim_right
      write(*,*) 'rare2cont  = ', me%prim_rare2cont
      write(*,*) 'cont2shock = ', me%prim_cont2shock
      write(*,*) 'rare       = ', me%prim_rare
      write(*,*) 'p_mean     = ', me%p_mean
      write(*,*) 'u_mean     = ', me%u_mean
      write(*,*) 'rare_c     = ', me%rare_c
      write(*,*) 'rare_c_fac = ', me%rare_c_fac
      write(*,*) 'cl         = ', me%cl
      write(*,*) 'cr         = ', me%cr
      write(*,*) 'charvel    = ', me%charvel
      write(*,*) 'charfield  = ', me%charfield
    end if

  end subroutine atl_ere_dump_solState
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Initial setup for the iterative computation of the solution to the
  !! Riemann problem.
  elemental subroutine nr_start(p, rhol, rhor, &
    &                              ul,   ur,   &
    &                              pl,   pr,   &
    &                              al,   ar,   &
    &                           me)
    real(kind=rk), intent(out) :: p
    real(kind=rk), intent(in) :: rhol, rhor
    real(kind=rk), intent(in) :: ul,   ur
    real(kind=rk), intent(in) :: pl,   pr
    real(kind=rk), intent(in) :: al,   ar
    type(atl_ere_solState1D_type), intent(in) :: me

    real(kind=rk) :: pv, pmin,pmax
    real(kind=rk) :: pnu, pde, qmax
    real(kind=rk) :: qrat, gel, ger

    qmax = 2.0_rk

    pv   = 0.5_rk*(pl+pr) - 0.125_rk*(ur-ul)*(rhol+rhor)*(al+ar)
    pmin = min(pl,pr)
    pmax = max(pl,pr)
    qrat = pmax/pmin

    if ( (qrat <= qmax) .and. ((pmin <= pv) .and. (pv <= pmax)) ) then
       p = max(me%tolerance,pv)
    else
      if (pv <= pmin) then
         pnu = al + ar - me%G(7)*(ur-ul)
         pde = al/pl**me%G(1) + ar/pr**me%G(1)
         p   = (pnu/pde)**me%G(3)
      else
         gel = sqrt((me%G(5)/rhol) / (me%G(6)*pl + max(me%tolerance,pv)))
         ger = sqrt((me%G(5)/rhor) / (me%G(6)*pr + max(me%tolerance,pv)))
         p   = (gel*pl + ger*pr - (ur-ul)) / (gel+ger)
         p   = max(me%tolerance,p)
      end if
    end if
  end subroutine nr_start
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!


  ! ----------------------------------------------------------------------------!
  !> Compute one-sided iteration
  elemental subroutine nr_1side(f, fd, p, rhok, pk, ck, me)
    real(kind=rk), intent(out) :: f, fd
    real(kind=rk), intent(in) :: p
    real(kind=rk), intent(in) :: rhok, pk, ck
    type(atl_ere_solState1D_type), intent(in) :: me

    real(kind=rk) :: ak, bk
    real(kind=rk) :: qrt, prat

    if (p <= pk) then
      prat = p/pk
      f    = me%G(4)*ck*(prat**me%G(1) - 1.0_rk)
      fd   = (1.0_rk/(rhok*ck))*prat**(-me%G(2))
    else
      ak  = me%G(5)/rhok
      bk  = me%G(6)*pk
      qrt = sqrt(ak/(bk+p))
      f   = (p-pk) * qrt
      fd  = (1.0_rk - 0.5_rk*(p-pk)/(bk+p)) * qrt
    END IF

  end subroutine nr_1side
  ! ----------------------------------------------------------------------------!
  ! ----------------------------------------------------------------------------!

end module atl_exact_riemann_euler_module
