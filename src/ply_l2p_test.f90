! Copyright (c) 2020 Harald Klimach <harald.klimach@uni-siegen.de>
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
!
!> Unit test to check the L2P implementation.
program ply_l2p_test
  use env_module, only: rk
  use ply_nodeset_module, only: ply_nodeset_chebyshev, &
    &                           ply_nodeset_chebyloba, &
    &                           ply_nodeset_legendre
  use ply_l2p_header_module, only: ply_l2p_header_type, &
    &                              ply_l2p_header_define
  use ply_l2p_module, only: ply_l2p_type,  &
    &                       ply_init_l2p,  &
    &                       ply_l2p_trafo_1d

  implicit none

  type(ply_l2p_header_type) :: header
  type(ply_l2p_type) :: trafo
  real(kind=rk), allocatable :: modaldat(:)
  real(kind=rk), allocatable :: recovered(:)
  real(kind=rk), allocatable :: nodevals(:)
  real(kind=rk), allocatable :: coords(:)
  real(kind=rk) :: maxdiff
  real(kind=rk) :: diff
  integer :: degree

  degree = 4
  allocate(modaldat(degree+1))
  allocate(recovered(degree+1))
  allocate(nodevals(degree+1))
  allocate(coords(degree+1))

  ! ---- Checking trafo with Chebyshev-Lobatto
  ! Defining a polynomial of the form f(x) = 1.5 - 1.5*x*x
  modaldat = 0.0_rk
  modaldat(1) = 1.0_rk
  modaldat(3) = -1.0_rk

  call ply_l2p_header_define( me            = header,      &
    &                         nodes_kind    = 'chebyshev', &
    &                         lobattoPoints = .true.       )
  coords = ply_nodeset_chebyloba(degree+1)

  call ply_init_l2p( l2p    = trafo,  &
    &                header = header, &
    &                degree = degree  )

  call ply_l2p_trafo_1d( trafo     = trafo%leg2node, &
    &                    original  = modaldat,       &
    &                    projected = nodevals        )

  maxdiff = maxval(abs(nodevals - (1.5_rk - 1.5_rk*coords*coords)))
  write(*,*) 'Maximal difference in nodal data for Chebyshev-Lobatto:', &
    &        maxdiff

  ! ---- Checking trafo and back again
  degree = 24
  deallocate(modaldat, recovered, nodevals, coords)
  allocate(modaldat(degree+1))
  allocate(recovered(degree+1))
  allocate(nodevals(degree+1))

  call random_number(modaldat)

  ! -- Chebyshev-Lobatto (header from above)
  call ply_init_l2p( l2p    = trafo,  &
    &                header = header, &
    &                degree = degree  )

  call ply_l2p_trafo_1d( trafo     = trafo%leg2node, &
    &                    original  = modaldat,       &
    &                    projected = nodevals        )
  call ply_l2p_trafo_1d( trafo     = trafo%node2leg, &
    &                    original  = nodevals,       &
    &                    projected = recovered       )
  diff = maxval(abs(recovered - modaldat))
  write(*,*) 'Maximal difference after trafo and back for Chebyshev-Lobatto:', &
    &        diff
  maxdiff = max(diff, maxdiff)

  call ply_l2p_header_define( me            = header,      &
    &                         nodes_kind    = 'chebyshev', &
    &                         lobattoPoints = .false.      )
  call ply_init_l2p( l2p    = trafo,  &
    &                header = header, &
    &                degree = degree  )
  call ply_l2p_trafo_1d( trafo     = trafo%leg2node, &
    &                    original  = modaldat,       &
    &                    projected = nodevals        )
  call ply_l2p_trafo_1d( trafo     = trafo%node2leg, &
    &                    original  = nodevals,       &
    &                    projected = recovered       )
  diff = maxval(abs(recovered - modaldat))
  write(*,*) 'Maximal difference after trafo and back for Chebyshev:', &
    &        diff
  maxdiff = max(diff, maxdiff)

  call ply_l2p_header_define( me            = header,           &
    &                         nodes_kind    = 'gauss-legendre', &
    &                         lobattoPoints = .false.           )
  call ply_init_l2p( l2p    = trafo,  &
    &                header = header, &
    &                degree = degree  )
  call ply_l2p_trafo_1d( trafo     = trafo%leg2node, &
    &                    original  = modaldat,       &
    &                    projected = nodevals        )
  call ply_l2p_trafo_1d( trafo     = trafo%node2leg, &
    &                    original  = nodevals,       &
    &                    projected = recovered       )
  diff = maxval(abs(recovered - modaldat))
  write(*,*) 'Maximal difference after trafo and back for Gauss-Legendre:', &
    &        diff
  maxdiff = max(diff, maxdiff)

  if (maxdiff < sqrt(epsilon(maxdiff))) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED!'
  end if

end program ply_l2p_test
