! Copyright (c) 2015, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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

! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt and Harald Klimach for
! University of Siegen.
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
! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the number of degrees of freedom for Q polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * A variable to store the number of degrees of freedom for a P tensor product
!   polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.

! The x ansatz degree is turned into the degree of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x ansatz degree is turned into the degree of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
module atl_project_physflux_module
  use env_module,                only: rk
  ! ateles
  use atl_cube_elem_module,      only: atl_cube_elem_type
  use atl_equation_module,       only: atl_equations_type
  use atl_kerneldata_module,     only: atl_kerneldata_type
  use atl_scheme_module,         only: atl_modg_scheme_type
  use atl_modg_2d_scheme_module, only: atl_modg_2d_scheme_type
  ! polynomials
  use ply_dof_module,            only: Q_space, P_space
  use ply_modg_basis_module,     only: ply_scalProdDualLegDiff

  implicit none

  private

  public :: atl_modg_project_physflux_testfunc, &
    &       atl_modg_2d_project_physFlux_testfunc


contains


  !> Subroutine to project modal representations of physical flux, numerical flux
  !! and source terms onto test functions.
  subroutine atl_modg_project_PhysFlux_testFunc( mesh, equation, kerneldata, &
    &                                           scheme, iDir, dl_prod,       &
    &                                           dirVec, iElem, state_der     )
    ! --------------------------------------------------------------------------
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The direction
    integer, intent(in) :: iDir
    !> The parameters of the MODG scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    !> stored scalar products of the testfunction and ansatz function
    real(kind=rk), intent(in) :: dl_prod(2, scheme%maxPolyDegree+1)
    !> vector for direction indicators
    integer, intent(in) :: dirVec(3,3)
    integer, intent(in) :: iElem
    real(kind=rk), intent(in)  :: state_der(kerneldata%nDofs,equation%varSys%nScalars)
    ! --------------------------------------------------------------------------

    ! Projection of the physical flux
      select case(scheme%basisType)
      case(Q_space)

        select case(equation%varSys%nScalars)
        case(6)
          select case(iDir)
          case(1)
            call modg_prj_pFlux1_Q_6( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          case(2)
            call modg_prj_pFlux2_Q_6( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          case(3)
            call modg_prj_pFlux3_Q_6( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          end select
        case(5)
          select case(iDir)
          case(1)
            call modg_prj_pFlux1_Q_5( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          case(2)
            call modg_prj_pFlux2_Q_5( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          case(3)
            call modg_prj_pFlux3_Q_5( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          end select
        case default
          call modg_project_physFlux_Q( nScalars = equation%varSys%nScalars, &
            &                           maxPolyDegree = scheme%maxPolyDegree, &
            &                           length = mesh%length, &
            &                           dl_prod = dl_prod, &
            &                           state = kerneldata%state_der, &
            &                           dirVec = dirVec(:,iDir),     &
            &                           iElem  = iElem,              &
            &                           state_der = state_der        )
        end select

      case(P_space)
        call modg_project_physFlux_P( nScalars = equation%varSys%nScalars, &
          &                           maxPolyDegree = scheme%maxPolyDegree, &
          &                           length = mesh%length, &
          &                           dl_prod = dl_prod, &
          &                           iElem = iElem, &
          &                           state = kerneldata%state_der, &
          &                           dirVec = dirVec(:,iDir), &
          &                           state_der = state_der                  )
      end select


  end subroutine atl_modg_project_PhysFlux_testFunc


  !> Projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_project_physFlux_Q( nScalars, maxPolyDegree, length, state, &
    &                                 dl_prod, dirVec, iElem, state_der )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> ordering of xyz for current direction
    integer, intent(in) :: dirVec(3)
    integer, intent(in) :: iElem
    !> The state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj, scalProd1(maxPolyDegree)
    integer :: testPos
    integer :: iTest2, iTest1, iTest3, iTestVec(3), min2mpd
    integer :: iAnsVec(3)
    integer :: iVar
    integer :: ansPos(4)
    real(kind=rk) :: scalProd(4)
    integer :: jk
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**2

    do iTest1=2,maxPolyDegree+1
      scalProd1(iTest1-1) = ply_scalProdDualLegDiff(iTest1-1, iTest1) &
        &                   * jacobiDetStiffProj
    end do

    min2mpd = min(2, maxPolyDegree+1)

    ! unrolled loop

    do iTest3 = 1, min2mpd

      do jk=1,min2mpd*maxPolyDegree
        iTest2 = mod(jk-1,min2mpd) + 1
        iTest1 = (jk-1)/min2mpd + 2

        ! one entry

        iTestVec = (/iTest1, iTest2, iTest3/)
  testpos = itestvec(dirvec(1))                                      &
    &      + ( ( itestvec(dirvec(2))-1)                             &
    &      + (itestvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)

        iAnsVec = (/iTest1-1, iTest2, iTest3/)
  anspos(1) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(1) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(2,iTest3)

        do iVar=1,nScalars
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(ansPos(1),iVar) * scalProd(1)
        end do

      end do

      do jk=1,(maxPolyDegree-1)*maxPolyDegree
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 3

        ! two entries

        iTestVec = (/iTest1, iTest2, iTest3/)
  testpos = itestvec(dirvec(1))                                      &
    &      + ( ( itestvec(dirvec(2))-1)                             &
    &      + (itestvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAnsVec = (/iTest1-1, iTest2-2, iTest3/)
  anspos(1) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(1) =   scalProd1(iTest1-1) &
          &           * dl_prod(1,iTest2) &
          &           * dl_prod(2,iTest3)

        iAnsVec = (/iTest1-1, iTest2, iTest3/)
  anspos(2) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(2) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(2,iTest3)

        do iVar=1,nScalars
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(ansPos(1),iVar) * scalProd(1) &
            & + state_der(ansPos(2),iVar) * scalProd(2)
        end do

      end do

    end do

    do iTest3 = 3, maxPolyDegree+1

      do jk=1,min2mpd*maxPolyDegree
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 1

        ! two entries

        iTestVec = (/iTest1, iTest2, iTest3/)
  testpos = itestvec(dirvec(1))                                      &
    &      + ( ( itestvec(dirvec(2))-1)                             &
    &      + (itestvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAnsVec = (/iTest1-1, iTest2, iTest3-2/)
  anspos(1) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(1) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(1,iTest3)

        iAnsVec = (/iTest1-1, iTest2, iTest3/)
  anspos(2) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(2) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(2,iTest3)

        do iVar=1,nScalars
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar) &
            & + state_der(ansPos(1),iVar) * scalProd(1) &
            & + state_der(ansPos(2),iVar) * scalProd(2)
        end do
      end do

      do jk=1,maxPolyDegree*(maxPolyDegree-1)
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 3

        ! four entries

        iTestVec = (/iTest1, iTest2, iTest3/)
  testpos = itestvec(dirvec(1))                                      &
    &      + ( ( itestvec(dirvec(2))-1)                             &
    &      + (itestvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAnsVec = (/iTest1-1, iTest2-2, iTest3-2/)
  anspos(1) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(1) =   scalProd1(iTest1-1) &
          &           * dl_prod(1,iTest2) &
          &           * dl_prod(1,iTest3)

        iAnsVec = (/iTest1-1, iTest2, iTest3-2/)
  anspos(2) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(2) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(1,iTest3)

        iAnsVec = (/iTest1-1, iTest2-2, iTest3/)
  anspos(3) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(3) =   scalProd1(iTest1-1) &
          &           * dl_prod(1,iTest2) &
          &           * dl_prod(2,iTest3)

        iAnsVec = (/iTest1-1, iTest2, iTest3/)
  anspos(4) = iansvec(dirvec(1))                                      &
    &      + ( ( iansvec(dirvec(2))-1)                             &
    &      + (iansvec(dirvec(3))-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalProd(4) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(2,iTest3)

        do iVar=1,nScalars
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(ansPos(1),iVar) * scalProd(1) &
            & + state_der(ansPos(2),iVar) * scalProd(2) &
            & + state_der(ansPos(3),iVar) * scalProd(3) &
            & + state_der(ansPos(4),iVar) * scalProd(4)
        end do

      end do

    end do


  end subroutine modg_project_physFlux_Q


  !> Projection of the physical flux onto the testfunctions, with unrolled loops
  !! => fewer loop-overhead/instructions, but more "random" memory accesses
  !! MZ: perhaps this version is faster for low order (or always, depending on the machine?)
  subroutine modg_project_physFlux_P( nScalars, iElem, maxPolyDegree, length, &
    &                                 state, dl_prod, dirVec, state_der       )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> ordering of xyz for current direction
    integer, intent(in) :: dirVec(3)
    !> The state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    !> The element index
    integer, intent(in) :: iElem
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj, scalProd1(maxPolyDegree)
    integer :: testPos
    integer :: iTest2, iTest1, iTest3, iTestVec(3), min2mpd
    integer :: iAnsVec(3)
    integer :: iVar
    integer :: ansPos(4)
    real(kind=rk) :: scalProd(4)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! This is the stiffness term!
    !
    ! We have cubiic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**2

    do iTest1=2,maxPolyDegree+1
      scalProd1(iTest1-1) = ply_scalProdDualLegDiff(iTest1-1, iTest1) &
        &                   * jacobiDetStiffProj
    end do

    min2mpd = min(2, maxPolyDegree+1)

    ! unrolled loop

    do iTest3 = 1, min2mpd
      do iTest2 = 1, min(2, maxPolyDegree+1 - (iTest3-1))
        do iTest1 = 2, maxPolyDegree+1 - (iTest3-1) - (iTest2-1)

          ! one entry

          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)

          iAnsVec = (/iTest1-1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(1) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(2,iTest3)

          do iVar=1,nScalars
            state(iElem,testPos,iVar) &
              & = state(iElem,testPos,iVar) &
              & + state_der(ansPos(1),iVar) * scalProd(1)
          end do

        end do
      end do

      do iTest2 = 3, maxPolyDegree+1 - (iTest3-1)
        do iTest1 = 2, maxPolyDegree+1 - (iTest3-1) - (iTest2-1)

          ! two entries

          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)


          iAnsVec = (/iTest1-1, iTest2-2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(1) =   scalProd1(iTest1-1) &
            &           * dl_prod(1,iTest2) &
            &           * dl_prod(2,iTest3)

          iAnsVec = (/iTest1-1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(2) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(2) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(2,iTest3)

          do iVar=1,nScalars
            state(iElem,testPos,iVar) &
              & = state(iElem,testPos,iVar) &
              & + state_der(ansPos(1),iVar) * scalProd(1) &
              & + state_der(ansPos(2),iVar) * scalProd(2)
          end do

        end do
      end do
    end do

    do iTest3 = 3, maxPolyDegree+1
      do iTest2 = 1, min(2, maxPolyDegree+1 - (iTest3-1))
        do iTest1 = 2, maxPolyDegree+1 - (iTest3-1) - (iTest2-1)

          ! two entries

          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)


          iAnsVec = (/iTest1-1, iTest2, iTest3-2/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(1) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(1,iTest3)

          iAnsVec = (/iTest1-1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(2) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(2) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(2,iTest3)

          do iVar=1,nScalars
            state(iElem,testPos,iVar) &
              & = state(iElem,testPos,iVar) &
              & + state_der(ansPos(1),iVar) * scalProd(1) &
              & + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      end do

      do iTest2 = 3, maxPolyDegree+1 - (iTest3-1)
        do iTest1 = 2, maxPolyDegree+1 - (iTest3-1) - (iTest2-1)

          ! four entries

          iTestVec = (/iTest1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  testpos = (((itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 3) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 2) &
    &     * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((itestvec(dirvec(3))-1) * (itestvec(dirvec(1)) + itestvec(dirvec(2)) + itestvec(dirvec(3)) -2) &
    &   - ((itestvec(dirvec(3))-2) * (itestvec(dirvec(3))-1)) / 2) &
    & + (itestvec(dirvec(2))-1)


          iAnsVec = (/iTest1-1, iTest2-2, iTest3-2/)
  ! integer divisions are no mistake here.
  anspos(1) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(1) =   scalProd1(iTest1-1) &
            &           * dl_prod(1,iTest2) &
            &           * dl_prod(1,iTest3)

          iAnsVec = (/iTest1-1, iTest2, iTest3-2/)
  ! integer divisions are no mistake here.
  anspos(2) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(2) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(1,iTest3)

          iAnsVec = (/iTest1-1, iTest2-2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(3) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(3) =   scalProd1(iTest1-1) &
            &           * dl_prod(1,iTest2) &
            &           * dl_prod(2,iTest3)

          iAnsVec = (/iTest1-1, iTest2, iTest3/)
  ! integer divisions are no mistake here.
  anspos(4) = (((iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 3) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 2) &
    &     * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) - 1)) &
    &   / 6 + 1)             &
    & + ((iansvec(dirvec(3))-1) * (iansvec(dirvec(1)) + iansvec(dirvec(2)) + iansvec(dirvec(3)) -2) &
    &   - ((iansvec(dirvec(3))-2) * (iansvec(dirvec(3))-1)) / 2) &
    & + (iansvec(dirvec(2))-1)
          scalProd(4) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(2,iTest3)

          do iVar=1,nScalars
            state(iElem,testPos,iVar) &
              & = state(iElem,testPos,iVar) &
              & + state_der(ansPos(1),iVar) * scalProd(1) &
              & + state_der(ansPos(2),iVar) * scalProd(2) &
              & + state_der(ansPos(3),iVar) * scalProd(3) &
              & + state_der(ansPos(4),iVar) * scalProd(4)
          end do

        end do
      end do
    end do

  end subroutine modg_project_physFlux_P


  !> Subroutine to project modal representations of physical flux, numerical
  !! flux and source terms onto test functions.
  subroutine atl_modg_2d_project_physFlux_testFunc( mesh, equation,        &
    &                                               kerneldata, iElem,     &
    &                                               dl_prod, iDir, scheme, &
    &                                               state_data, ndofs      )
    ! --------------------------------------------------------------------------
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The parameters of the MODG scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    !> The total degrees of freedom
    integer, intent(in) :: nDofs
    !> The direction
    integer, intent(in) :: iDir
    !> The element index
    integer, intent(in) :: iElem
    real(kind=rk) , intent(in) :: dl_prod(2, scheme%maxPolyDegree+1)
    !> The physical fluxes that needs to be projected
    real(kind=rk), intent(in)  :: state_data(nDofs,equation%varSys%nScalars)
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Iterate over all the elements and do the following:
    ! 1. Project physical fluxes (3 directions) to test functions (stiffness
    !    terms).
    !    Attention: physical x flux: state_der(:,:,2,:)
    !    Attention: physical y flux: state_der(:,:,3,:)
    ! 2. Project numerical fluxes (4 faces) to test functions.
    ! 3. Project source terms to test functions.
    ! Attention: write the projections to state_der(:,:,1,:) because inverse
    ! of mass matrix will be applied to these entries.

    select case(scheme%basisType)
    case(Q_space)
      if (iDir == 1) then
        ! Projection of the physical flux
        ! ... x direction
        call modg_2d_project_physFluxX_Q(             &
          & nScalars      = equation%varSys%nScalars, &
          & maxPolyDegree = scheme%maxPolyDegree,     &
          & length        = mesh%length,              &
          & dl_prod       = dl_prod,                  &
          & iElem         = iElem,                    &
          & state         = kerneldata%state_der,     &
          & state_der     = state_data                )
      else
        ! ... y direction
        call modg_2d_project_physFluxY_Q(             &
          & nScalars      = equation%varSys%nScalars, &
          & maxPolyDegree = scheme%maxPolyDegree,     &
          & length        = mesh%length,              &
          & dl_prod       = dl_prod,                  &
          & iElem         = iElem,                    &
          & state         = kerneldata%state_der,     &
          & state_der     = state_data                )
      endif

    case(P_space)
      if (iDir == 1) then
        ! Projection of the physical flux
        ! ... x direction
        call modg_2d_project_physFluxX_P(             &
          & nScalars      = equation%varSys%nScalars, &
          & maxPolyDegree = scheme%maxPolyDegree,     &
          & length        = mesh%length,              &
          & dl_prod       = dl_prod,                  &
          & iElem         = iElem,                    &
          & state         = kerneldata%state_der,     &
          & nDofs         = nDofs,                    &
          & state_der     = state_data                )
      else
        ! ... y direction
        call modg_2d_project_physFluxY_P(             &
          & nScalars      = equation%varSys%nScalars, &
          & maxPolyDegree = scheme%maxPolyDegree,     &
          & length        = mesh%length,              &
          & dl_prod       = dl_prod,                  &
          & iElem         = iElem,                    &
          & state         = kerneldata%state_der,     &
          & nDofs         = nDofs,                    &
          & state_der     = state_data                )
      endif
    end select

  end subroutine atl_modg_2d_project_physFlux_testfunc


  !> Projection of the physical flux in x direction onto the testfunctions.
  subroutine modg_2d_project_physFluxX_Q( nScalars, maxPolyDegree, length, &
    &                                     dl_prod, state, iElem, state_der )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> The element index
    integer, intent(in) :: iElem
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der((maxPolyDegree+1)**2,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj
    real(kind=rk) :: scalProdX(maxPolyDegree)
    integer :: ansPos(2), testPos
    integer :: iTestY
    integer :: iAnsX
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(2)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test
    ! functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**1

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    ! for x direction (x test function differentiated)
    ! get the relevant indices for the ansatz function
    do iAnsX=1,maxPolyDegree
      scalProdX(iAnsX) = ply_scalProdDualLegDiff(iAnsX, iAnsX+1) &
        &                * jacobiDetStiffProj
    end do

    ! now, project onto all test functions
    YLoop: do iTestY = 1, maxPolyDegree+1

  testpos = 1                                      &
    &      + ( ( itesty-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
      do iVar=var_lb,var_ub
        state(iElem,testPos,iVar) = 0.0_rk
      end do

      if (iTestY > 2) then
        ! Need to add two terms
        do iAnsX = 1, maxPolyDegree
          ! the position of the current test functions
  testpos = iansx+1                                      &
    &      + ( ( itesty-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          scalProd(1) = dl_prod(1, iTestY) * scalProdX(iAnsX)
          scalProd(2) = dl_prod(2, iTestY) * scalProdX(iAnsX)
  anspos(1) = iansx                                      &
    &      + ( ( itesty-2-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(2) = iansx                                      &
    &      + ( ( itesty-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                      &
              &  = state(iElem,testPos,iVar)               &
              &  + state_der(ansPos(1),iVar) * scalProd(1) &
              &  + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      else
        ! Need to add one term
        do iAnsX = 1, maxPolyDegree
          ! the position of the current test functions
  testpos = iansx+1                                      &
    &      + ( ( itesty-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          scalProd(1) = dl_prod(2, iTestY) * scalProdX(iAnsX)
  anspos(1) = iansx                                      &
    &      + ( ( itesty-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                     &
              & = state(iElem,testPos,iVar)               &
              &   + state_der(ansPos(1),iVar) * scalProd(1)
          end do

        end do
      end if
    end do YLoop

  end subroutine modg_2d_project_physFluxX_Q


  !> Projection of the physical flux in y direction onto the testfunctions.
  subroutine modg_2d_project_physFluxY_Q( nScalars, maxPolyDegree, length,  &
    &                                     dl_prod, state , ielem, state_der )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> The element index
    integer, intent(in) :: iElem
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der((maxPolyDegree+1)**2,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj
    real(kind=rk) :: scalProdY(maxPolyDegree)
    integer :: ansPos(2), testPos
    integer :: iTestX
    integer :: iAnsY
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(2)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test
    ! functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**1

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    !  for y direction (y test function differentiated)
    do iAnsY=1,maxPolyDegree
      scalProdY(iAnsY) = ply_scalProdDualLegDiff(iAnsY, iAnsY+1) &
        &                * jacobiDetStiffProj
    end do

    ! now, project onto all test functions
    XLoop: do iTestX = 1, maxPolyDegree+1

      if (iTestX > 2) then
        ! Need to add two terms
        do iAnsY = 1, maxPolyDegree
  testpos = itestx                                      &
    &      + ( ( iansy+1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          scalProd(1) = dl_prod(1, iTestX) * scalProdY(iAnsY)
          scalProd(2) = dl_prod(2, iTestX) * scalProdY(iAnsY)
  anspos(1) = itestx-2                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
  anspos(2) = itestx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                      &
              &  = state(iElem,testPos,iVar)               &
              &  + state_der(ansPos(1),iVar) * scalProd(1) &
              &  + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      else
        ! Need to add one term
        do iAnsY = 1, maxPolyDegree
  testpos = itestx                                      &
    &      + ( ( iansy+1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          scalProd(1) = dl_prod(2, iTestX) * scalProdY(iAnsY)
  anspos(1) = itestx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iVar = var_lb,var_ub
            state(iElem,testPos,iVar)                    &
              &  = state(iElem,testPos,iVar)             &
              &  + state_der(ansPos(1),iVar) * scalProd(1)
          end do
        end do
      end if
    end do XLoop

  end subroutine modg_2d_project_physFluxY_Q


  !> Projection of the physical flux in x direction onto the testfunctions.
  subroutine modg_2d_project_physFluxX_P( nScalars, maxPolyDegree, length, &
    &                                     dl_prod, state, iElem, nDofs,    &
    &                                     state_der                        )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> The element index
    integer, intent(in) :: iElem
    !> Number of degrees of freedom
    integer, intent(in) :: nDofs
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der(nDofs,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj
    real(kind=rk) :: scalProdX(maxPolyDegree)
    integer :: ansPos(2), testPos
    integer :: iTestY
    integer :: iAnsX
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(2)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! This is the stiffness term!
    !
    ! We have cubiic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.

    jacobiDetStiffProj = (0.5_rk*length)**1

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    ! for x direction (x test function differentiated)
    ! get the relevant indices for the ansatz function
    do iAnsX=1,maxPolyDegree
      scalProdX(iAnsX) = ply_scalProdDualLegDiff(iAnsX, iAnsX+1) &
        &                * jacobiDetStiffProj
    end do

    ! now, project onto all test functions
    YLoop: do iTestY = 1, maxPolyDegree+1

  ! integer divisions are no mistake here.
  testpos = ((((1 - 1) + (itesty - 1))            &
    &   * (((1 - 1) + (itesty - 1)) + 1)) / 2 + 1) &
    & + (itesty - 1)
      do iVar=var_lb,var_ub
        state(iElem,testPos,iVar) = 0.0_rk
      end do

      if (iTestY > 2) then
        ! Need to add two terms
        do iAnsX = 1, maxPolyDegree-(iTestY-1)
          ! the position of the current test functions
  ! integer divisions are no mistake here.
  testpos = ((((iansx+1 - 1) + (itesty - 1))            &
    &   * (((iansx+1 - 1) + (itesty - 1)) + 1)) / 2 + 1) &
    & + (itesty - 1)
          scalProd(1) = dl_prod(1, iTestY) * scalProdX(iAnsX)
          scalProd(2) = dl_prod(2, iTestY) * scalProdX(iAnsX)
  ! integer divisions are no mistake here.
  anspos(1) = ((((iansx - 1) + (itesty-2 - 1))            &
    &   * (((iansx - 1) + (itesty-2 - 1)) + 1)) / 2 + 1) &
    & + (itesty-2 - 1)
  ! integer divisions are no mistake here.
  anspos(2) = ((((iansx - 1) + (itesty - 1))            &
    &   * (((iansx - 1) + (itesty - 1)) + 1)) / 2 + 1) &
    & + (itesty - 1)
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                      &
              &  = state(iElem,testPos,iVar)               &
              &  + state_der(ansPos(1),iVar) * scalProd(1) &
              &  + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      else
        ! Need to add one term
        do iAnsX = 1, maxPolyDegree - (iTestY-1)
          ! the position of the current test functions
  ! integer divisions are no mistake here.
  testpos = ((((iansx+1 - 1) + (itesty - 1))            &
    &   * (((iansx+1 - 1) + (itesty - 1)) + 1)) / 2 + 1) &
    & + (itesty - 1)
          scalProd(1) = dl_prod(2, iTestY) * scalProdX(iAnsX)
  ! integer divisions are no mistake here.
  anspos(1) = ((((iansx - 1) + (itesty - 1))            &
    &   * (((iansx - 1) + (itesty - 1)) + 1)) / 2 + 1) &
    & + (itesty - 1)
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                     &
              & = state(iElem,testPos,iVar)               &
              &   + state_der(ansPos(1),iVar) * scalProd(1)
          end do

        end do
      end if
    end do YLoop

  end subroutine modg_2d_project_physFluxX_P


  !> Projection of the physical flux in y direction onto the testfunctions.
  subroutine modg_2d_project_physFluxY_P( nScalars, maxPolyDegree, length, &
    &                                     dl_prod, state , iElem, nDofs,   &
    &                                     state_der                        )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> The element index
    integer, intent(in) :: iElem
    !> Number of degrees of freedom
    integer, intent(in) :: nDofs
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der(nDofs,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj
    real(kind=rk) :: scalProdY(maxPolyDegree)
    integer :: ansPos(2), testPos
    integer :: iTestX
    integer :: iAnsY
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(2)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test
    ! functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**1

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    !  for y direction (y test function differentiated)
    do iAnsY=1,maxPolyDegree
      scalProdY(iAnsY) = ply_scalProdDualLegDiff(iAnsY, iAnsY+1) &
        &                * jacobiDetStiffProj
    end do

    ! now, project onto all test functions
    XLoop: do iTestX = 1, maxPolyDegree+1

      if (iTestX > 2) then
        ! Need to add two terms
        do iAnsY = 1, maxPolyDegree - (iTestX-1)
  ! integer divisions are no mistake here.
  testpos = ((((itestx - 1) + (iansy+1 - 1))            &
    &   * (((itestx - 1) + (iansy+1 - 1)) + 1)) / 2 + 1) &
    & + (iansy+1 - 1)
          scalProd(1) = dl_prod(1, iTestX) * scalProdY(iAnsY)
          scalProd(2) = dl_prod(2, iTestX) * scalProdY(iAnsY)
  ! integer divisions are no mistake here.
  anspos(1) = ((((itestx-2 - 1) + (iansy - 1))            &
    &   * (((itestx-2 - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
  ! integer divisions are no mistake here.
  anspos(2) = ((((itestx - 1) + (iansy - 1))            &
    &   * (((itestx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                      &
              &  = state(iElem,testPos,iVar)               &
              &  + state_der(ansPos(1),iVar) * scalProd(1) &
              &  + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      else
        ! Need to add one term
        do iAnsY = 1, maxPolyDegree - (iTestX-1)
  ! integer divisions are no mistake here.
  testpos = ((((itestx - 1) + (iansy+1 - 1))            &
    &   * (((itestx - 1) + (iansy+1 - 1)) + 1)) / 2 + 1) &
    & + (iansy+1 - 1)
          scalProd(1) = dl_prod(2, iTestX) * scalProdY(iAnsY)
  ! integer divisions are no mistake here.
  anspos(1) = ((((itestx - 1) + (iansy - 1))            &
    &   * (((itestx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)
          do iVar = var_lb,var_ub
            state(iElem,testPos,iVar)                    &
              &  = state(iElem,testPos,iVar)             &
              &  + state_der(ansPos(1),iVar) * scalProd(1)
          end do
        end do
      end if
    end do XLoop

  end subroutine modg_2d_project_physFluxY_P



  !!dirVec(:,1) = [ 1,2,3 ]
  !!dirVec(:,2) = [ 2,1,3 ]
  !!dirVec(:,3) = [ 2,3,1 ]

  !> X direction for 6 scalars
  !> projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_prj_pflux1_q_6( maxpolydegree, length, state, dl_prod, &
    &                                 ielem, state_der )
    ! --------------------------------------------------------------------------
    !> the maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxpolydegree
    !> the length of the cubes.
    real(kind=rk), intent(in) :: length
    !> the state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> precomputed dual legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxpolydegree+1)
    integer, intent(in) :: ielem
    !> the state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobidetstiffproj, prescal !!(maxpolydegree)
    integer :: testpos
    integer :: itest2, itest1, itest3, min2mpd
    integer :: ians1, ians2, ians3
    integer :: ivar
    integer :: anspos1, anspos2, anspos3, anspos4
    real(kind=rk) :: scalprod1, scalprod2, scalprod3, scalprod4
    integer :: jk
    ! --------------------------------------------------------------------------

    ! jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! this is the stiffness term!
    !
    ! we have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobidetstiffproj = (0.5_rk*length)**2

    !!do itest1=2,maxpolydegree+1
    !!  prescal(itest1-1) = ply_scalprodduallegdiff(itest1-1, itest1) &
    !!    &                 * jacobidetstiffproj
    !!end do
    ! we only consider the non-zeroes here, and these are all 2.
    prescal = 2*jacobidetstiffproj

    min2mpd = min(2, maxpolydegree+1)

    ! unrolled loop

    do itest3 = 1, min2mpd

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest2 = mod(jk-1,min2mpd) + 1
        itest1 = (jk-1)/min2mpd + 2

        ! one entry

        testpos = itest1                                      &
          &      + ( ( itest2-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos1 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1
        end do

      end do

      !$nec ivdep
      do jk=1,(maxpolydegree-1)*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! two entries

        testpos = itest1                                      &
          &      + ( ( itest2-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos1 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do

      end do

    end do

    do itest3 = 3, maxpolydegree+1

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 1

        ! two entries

        testpos = itest1                                      &
          &      + ( ( itest2-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos1 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar) &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do
      end do

      !$nec ivdep
      do jk=1,maxpolydegree*(maxpolydegree-1)
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! four entries

        testpos = itest1                                      &
          &      + ( ( itest2-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3-2
        anspos1 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos2 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos3 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod3 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos4 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod4 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2 &
            & + state_der(anspos3,ivar) * scalprod3 &
            & + state_der(anspos4,ivar) * scalprod4
        end do

      end do

    end do


  end subroutine modg_prj_pflux1_q_6

  !> Y direction for 6 scalars
  !> projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_prj_pflux2_q_6( maxpolydegree, length, state, dl_prod, &
    &                                 ielem, state_der )
    ! --------------------------------------------------------------------------
    !> the maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxpolydegree
    !> the length of the cubes.
    real(kind=rk), intent(in) :: length
    !> the state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> precomputed dual legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxpolydegree+1)
    integer, intent(in) :: ielem
    !> the state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobidetstiffproj, prescal !!(maxpolydegree)
    integer :: testpos
    integer :: itest2, itest1, itest3, min2mpd
    integer :: ians1, ians2, ians3
    integer :: ivar
    integer :: anspos1, anspos2, anspos3, anspos4
    real(kind=rk) :: scalprod1, scalprod2, scalprod3, scalprod4
    integer :: jk
    ! --------------------------------------------------------------------------

    ! jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! this is the stiffness term!
    !
    ! we have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobidetstiffproj = (0.5_rk*length)**2

    !!do itest1=2,maxpolydegree+1
    !!  prescal(itest1-1) = ply_scalprodduallegdiff(itest1-1, itest1) &
    !!    &                 * jacobidetstiffproj
    !!end do
    ! we only consider the non-zeroes here, and these are all 2.
    prescal = 2*jacobidetstiffproj

    min2mpd = min(2, maxpolydegree+1)

    ! unrolled loop

    do itest3 = 1, min2mpd

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest2 = mod(jk-1,min2mpd) + 1
        itest1 = (jk-1)/min2mpd + 2

        ! one entry

        testpos = itest2                                      &
          &      + ( ( itest1-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos1 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1
        end do

      end do

      !$nec ivdep
      do jk=1,(maxpolydegree-1)*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! two entries

        testpos = itest2                                      &
          &      + ( ( itest1-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos1 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do

      end do

    end do

    do itest3 = 3, maxpolydegree+1

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 1

        ! two entries

        testpos = itest2                                      &
          &      + ( ( itest1-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos1 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar) &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do
      end do

      !$nec ivdep
      do jk=1,maxpolydegree*(maxpolydegree-1)
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! four entries

        testpos = itest2                                      &
          &      + ( ( itest1-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3-2
        anspos1 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos2 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos3 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod3 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos4 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod4 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2 &
            & + state_der(anspos3,ivar) * scalprod3 &
            & + state_der(anspos4,ivar) * scalprod4
        end do

      end do

    end do


  end subroutine modg_prj_pflux2_q_6

  !> Z direction for 6 scalars
  !> projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_prj_pflux3_q_6( maxpolydegree, length, state, dl_prod, &
    &                                 ielem, state_der )
    ! --------------------------------------------------------------------------
    !> the maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxpolydegree
    !> the length of the cubes.
    real(kind=rk), intent(in) :: length
    !> the state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> precomputed dual legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxpolydegree+1)
    integer, intent(in) :: ielem
    !> the state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobidetstiffproj, prescal !!(maxpolydegree)
    integer :: testpos
    integer :: itest2, itest1, itest3, min2mpd
    integer :: ians1, ians2, ians3
    integer :: ivar
    integer :: anspos1, anspos2, anspos3, anspos4
    real(kind=rk) :: scalprod1, scalprod2, scalprod3, scalprod4
    integer :: jk
    ! --------------------------------------------------------------------------

    ! jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! this is the stiffness term!
    !
    ! we have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobidetstiffproj = (0.5_rk*length)**2

    !!do itest1=2,maxpolydegree+1
    !!  prescal(itest1-1) = ply_scalprodduallegdiff(itest1-1, itest1) &
    !!    &                 * jacobidetstiffproj
    !!end do
    ! we only consider the non-zeroes here, and these are all 2.
    prescal = 2*jacobidetstiffproj

    min2mpd = min(2, maxpolydegree+1)

    ! unrolled loop

    do itest3 = 1, min2mpd

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest2 = mod(jk-1,min2mpd) + 1
        itest1 = (jk-1)/min2mpd + 2

        ! one entry

        testpos = itest2                                      &
          &      + ( ( itest3-1)                             &
          &      + (itest1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos1 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1
        end do

      end do

      !$nec ivdep
      do jk=1,(maxpolydegree-1)*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! two entries

        testpos = itest2                                      &
          &      + ( ( itest3-1)                             &
          &      + (itest1-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos1 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do

      end do

    end do

    do itest3 = 3, maxpolydegree+1

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 1

        ! two entries

        testpos = itest2                                      &
          &      + ( ( itest3-1)                             &
          &      + (itest1-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos1 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar) &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do
      end do

      !$nec ivdep
      do jk=1,maxpolydegree*(maxpolydegree-1)
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! four entries

        testpos = itest2                                      &
          &      + ( ( itest3-1)                             &
          &      + (itest1-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3-2
        anspos1 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos2 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos3 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod3 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos4 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod4 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(6)
        do ivar=1,6
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2 &
            & + state_der(anspos3,ivar) * scalprod3 &
            & + state_der(anspos4,ivar) * scalprod4
        end do

      end do

    end do


  end subroutine modg_prj_pflux3_q_6

  !> X direction for 5 scalars
  !> projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_prj_pflux1_q_5( maxpolydegree, length, state, dl_prod, &
    &                                 ielem, state_der )
    ! --------------------------------------------------------------------------
    !> the maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxpolydegree
    !> the length of the cubes.
    real(kind=rk), intent(in) :: length
    !> the state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> precomputed dual legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxpolydegree+1)
    integer, intent(in) :: ielem
    !> the state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobidetstiffproj, prescal !!(maxpolydegree)
    integer :: testpos
    integer :: itest2, itest1, itest3, min2mpd
    integer :: ians1, ians2, ians3
    integer :: ivar
    integer :: anspos1, anspos2, anspos3, anspos4
    real(kind=rk) :: scalprod1, scalprod2, scalprod3, scalprod4
    integer :: jk
    ! --------------------------------------------------------------------------

    ! jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! this is the stiffness term!
    !
    ! we have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobidetstiffproj = (0.5_rk*length)**2

    !!do itest1=2,maxpolydegree+1
    !!  prescal(itest1-1) = ply_scalprodduallegdiff(itest1-1, itest1) &
    !!    &                 * jacobidetstiffproj
    !!end do
    ! we only consider the non-zeroes here, and these are all 2.
    prescal = 2*jacobidetstiffproj

    min2mpd = min(2, maxpolydegree+1)

    ! unrolled loop

    do itest3 = 1, min2mpd

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest2 = mod(jk-1,min2mpd) + 1
        itest1 = (jk-1)/min2mpd + 2

        ! one entry

        testpos = itest1                                      &
          &      + ( ( itest2-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos1 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1
        end do

      end do

      !$nec ivdep
      do jk=1,(maxpolydegree-1)*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! two entries

        testpos = itest1                                      &
          &      + ( ( itest2-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos1 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do

      end do

    end do

    do itest3 = 3, maxpolydegree+1

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 1

        ! two entries

        testpos = itest1                                      &
          &      + ( ( itest2-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos1 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar) &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do
      end do

      !$nec ivdep
      do jk=1,maxpolydegree*(maxpolydegree-1)
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! four entries

        testpos = itest1                                      &
          &      + ( ( itest2-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3-2
        anspos1 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos2 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos3 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod3 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos4 = ians1                                      &
          &      + ( ( ians2-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod4 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2 &
            & + state_der(anspos3,ivar) * scalprod3 &
            & + state_der(anspos4,ivar) * scalprod4
        end do

      end do

    end do


  end subroutine modg_prj_pflux1_q_5

  !> Y direction for 5 scalars
  !> projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_prj_pflux2_q_5( maxpolydegree, length, state, dl_prod, &
    &                                 ielem, state_der )
    ! --------------------------------------------------------------------------
    !> the maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxpolydegree
    !> the length of the cubes.
    real(kind=rk), intent(in) :: length
    !> the state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> precomputed dual legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxpolydegree+1)
    integer, intent(in) :: ielem
    !> the state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobidetstiffproj, prescal !!(maxpolydegree)
    integer :: testpos
    integer :: itest2, itest1, itest3, min2mpd
    integer :: ians1, ians2, ians3
    integer :: ivar
    integer :: anspos1, anspos2, anspos3, anspos4
    real(kind=rk) :: scalprod1, scalprod2, scalprod3, scalprod4
    integer :: jk
    ! --------------------------------------------------------------------------

    ! jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! this is the stiffness term!
    !
    ! we have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobidetstiffproj = (0.5_rk*length)**2

    !!do itest1=2,maxpolydegree+1
    !!  prescal(itest1-1) = ply_scalprodduallegdiff(itest1-1, itest1) &
    !!    &                 * jacobidetstiffproj
    !!end do
    ! we only consider the non-zeroes here, and these are all 2.
    prescal = 2*jacobidetstiffproj

    min2mpd = min(2, maxpolydegree+1)

    ! unrolled loop

    do itest3 = 1, min2mpd

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest2 = mod(jk-1,min2mpd) + 1
        itest1 = (jk-1)/min2mpd + 2

        ! one entry

        testpos = itest2                                      &
          &      + ( ( itest1-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos1 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1
        end do

      end do

      !$nec ivdep
      do jk=1,(maxpolydegree-1)*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! two entries

        testpos = itest2                                      &
          &      + ( ( itest1-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos1 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do

      end do

    end do

    do itest3 = 3, maxpolydegree+1

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 1

        ! two entries

        testpos = itest2                                      &
          &      + ( ( itest1-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos1 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar) &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do
      end do

      !$nec ivdep
      do jk=1,maxpolydegree*(maxpolydegree-1)
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! four entries

        testpos = itest2                                      &
          &      + ( ( itest1-1)                             &
          &      + (itest3-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3-2
        anspos1 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos2 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos3 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod3 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos4 = ians2                                      &
          &      + ( ( ians1-1)                             &
          &      + (ians3-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod4 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2 &
            & + state_der(anspos3,ivar) * scalprod3 &
            & + state_der(anspos4,ivar) * scalprod4
        end do

      end do

    end do


  end subroutine modg_prj_pflux2_q_5

  !> Z direction for 5 scalars
  !> projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_prj_pflux3_q_5( maxpolydegree, length, state, dl_prod, &
    &                                 ielem, state_der )
    ! --------------------------------------------------------------------------
    !> the maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxpolydegree
    !> the length of the cubes.
    real(kind=rk), intent(in) :: length
    !> the state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> precomputed dual legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxpolydegree+1)
    integer, intent(in) :: ielem
    !> the state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobidetstiffproj, prescal !!(maxpolydegree)
    integer :: testpos
    integer :: itest2, itest1, itest3, min2mpd
    integer :: ians1, ians2, ians3
    integer :: ivar
    integer :: anspos1, anspos2, anspos3, anspos4
    real(kind=rk) :: scalprod1, scalprod2, scalprod3, scalprod4
    integer :: jk
    ! --------------------------------------------------------------------------

    ! jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! this is the stiffness term!
    !
    ! we have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobidetstiffproj = (0.5_rk*length)**2

    !!do itest1=2,maxpolydegree+1
    !!  prescal(itest1-1) = ply_scalprodduallegdiff(itest1-1, itest1) &
    !!    &                 * jacobidetstiffproj
    !!end do
    ! we only consider the non-zeroes here, and these are all 2.
    prescal = 2*jacobidetstiffproj

    min2mpd = min(2, maxpolydegree+1)

    ! unrolled loop

    do itest3 = 1, min2mpd

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest2 = mod(jk-1,min2mpd) + 1
        itest1 = (jk-1)/min2mpd + 2

        ! one entry

        testpos = itest2                                      &
          &      + ( ( itest3-1)                             &
          &      + (itest1-1)*(maxpolydegree+1))*(maxpolydegree+1)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos1 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1
        end do

      end do

      !$nec ivdep
      do jk=1,(maxpolydegree-1)*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! two entries

        testpos = itest2                                      &
          &      + ( ( itest3-1)                             &
          &      + (itest1-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos1 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do

      end do

    end do

    do itest3 = 3, maxpolydegree+1

      !$nec ivdep
      do jk=1,min2mpd*maxpolydegree
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 1

        ! two entries

        testpos = itest2                                      &
          &      + ( ( itest3-1)                             &
          &      + (itest1-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos1 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos2 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar) &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2
        end do
      end do

      !$nec ivdep
      do jk=1,maxpolydegree*(maxpolydegree-1)
        itest1 = mod(jk-1,maxpolydegree) + 2
        itest2 = (jk-1)/maxpolydegree + 3

        ! four entries

        testpos = itest2                                      &
          &      + ( ( itest3-1)                             &
          &      + (itest1-1)*(maxpolydegree+1))*(maxpolydegree+1)


        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3-2
        anspos1 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3-2
        anspos2 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,itest2) * dl_prod(1,itest3)

        ians1 = itest1-1
        ians2 = itest2-2
        ians3 = itest3
        anspos3 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod3 = prescal * dl_prod(1,itest2) * dl_prod(2,itest3)

        ians1 = itest1-1
        ians2 = itest2
        ians3 = itest3
        anspos4 = ians2                                      &
          &      + ( ( ians3-1)                             &
          &      + (ians1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod4 = prescal * dl_prod(2,itest2) * dl_prod(2,itest3)

        !$nec unroll(5)
        do ivar=1,5
          state(ielem,testpos,ivar) &
            & = state(ielem,testpos,ivar)  &
            & + state_der(anspos1,ivar) * scalprod1 &
            & + state_der(anspos2,ivar) * scalprod2 &
            & + state_der(anspos3,ivar) * scalprod3 &
            & + state_der(anspos4,ivar) * scalprod4
        end do

      end do

    end do


  end subroutine modg_prj_pflux3_q_5

end module atl_project_physflux_module
