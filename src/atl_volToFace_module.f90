! Copyright (c) 2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> module routines to convert the data/gradient of data  from volume to face
!! and all related routines
module atl_volToFace_module
  use env_module,               only: rk
  use tem_faceData_module,      only: tem_dirToFace_map
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit
  use tem_param_module,         only: q__E, q__W, q__N, q__S, q__T, q__B, &
    &                                 qAxis
  use tem_element_module,       only: eT_fluid

  use ply_leg_diff_module,      only: ply_calcDiff_leg,    &
    &                                 ply_calcDiff_leg_2d, &
    &                                 ply_calcDiff_leg_1d
  use ply_modg_basis_module,    only: ply_faceValLeftBndAns
  use ply_dof_module,           only: Q_space

  use atl_cube_elem_module,     only: atl_cube_elem_type
  use atl_facedata_module,      only: atl_facedata_type
  use atl_kerneldata_module,    only: atl_statedata_type
  use atl_equation_module,      only: atl_equations_type

  private

  public ::                            &
    & atl_modg_volToFace_Q,            &
    & atl_modg_volToFace_grad_Q,       &
    & atl_modg_2d_volToFace_Q,         &
    & atl_modg_2d_volToFace_grad_Q,    &
    & atl_modg_1d_volToFace_Q,         &
    & atl_modg_1d_volToFace_grad_Q,    &
    & atl_modg_1d_modalVolToModalFace, &
    & atl_modg_volToFace_P,            &
    & atl_modg_2d_volToFace_P

  contains

  !> Project modal representation of an element to one of its faces for Q space.
  !!
  !! Project modal representation of an element onto one of its faces.
  !! Therefore, this function returns the modal representation of the solution
  !! on the face. This function can project onto an arbitrary face direction.
  subroutine atl_modg_volToFace_Q( nScalars, volState, maxPolyDegree, faceDir, &
    &                              nELems, faceState                           )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: nScalars
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink,
    !! e.g. \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: pos, facePos, iAnsX, iAnsY, iAnsZ, leftOrRight
    integer :: iVar
    real(kind=rk) :: faceVal
    integer :: mpd1, mpd1_square
    ! --------------------------------------------------------------------------

    leftOrRight = tem_dirToFace_map(faceDir)
    mpd1 = maxpolydegree+1
    mpd1_square = mpd1**2

    ! now, project to a certain face.
    if( faceDir.eq.q__E ) then

      do facepos = 1,mpd1_square
        do iAnsX = 1, maxPolyDegree+1
          ! get position of the current ansatz function
          pos = (facepos-1)*mpd1 + iAnsX

          do iVar = 1, nScalars
            faceState(:nElems,facePos,iVar,leftOrRight) &
              &  = faceState(:nElems,facePos,iVar,leftOrRight) &
              &  + volState(:nElems,pos,iVar)
          end do

        end do

      end do

    elseif(faceDir.eq.q__W) then

      do facepos = 1,mpd1_square

        do iAnsX = 1, maxPolyDegree+1
          ! get position of the current ansatz function
          pos = (facepos-1)*mpd1 + iAnsX
          ! get the face value of the ansatz function with fixed coordinate
          faceVal = ply_faceValLeftBndAns(iAnsX)

          do iVar = 1, nScalars
            faceState(:nElems,facePos,iVar,leftOrRight) &
              &  = faceState(:nElems,facePos,iVar,leftOrRight) &
              &  + faceVal*volState(:nElems,pos,iVar)
          end do

        end do

      end do

    elseif(faceDir.eq.q__S) then

      do facepos = 1,mpd1_square
        iAnsX = mod(facepos-1, mpd1) + 1
        iAnsZ = (facepos-1)/mpd1 + 1

        do iAnsY = 1, maxPolyDegree+1
          ! get position of the current ansatz function
          pos = iAnsX + (iAnsY-1)*mpd1 + (iAnsZ-1)*mpd1_square
          ! get the face value of the ansatz function with fixed coordinate
          faceVal = ply_faceValLeftBndAns(iAnsY)

          do iVar = 1, nScalars
            faceState(:nElems,facePos,iVar,leftOrRight) &
              &  = faceState(:nElems,facePos,iVar,leftOrRight) &
              &  + faceVal*volState(:nElems,pos,iVar)
          end do

        end do

      end do

    elseif(faceDir.eq.q__N) then

      do facepos = 1,mpd1_square
        iAnsX = mod(facepos-1, mpd1) + 1
        iAnsZ = (facepos-1)/mpd1 + 1

        do iAnsY = 1, maxPolyDegree+1
          ! get position of the current ansatz function
          pos = iAnsX + (iAnsY-1)*mpd1 + (iAnsZ-1)*mpd1_square

          do iVar = 1, nScalars
            faceState(:nElems,facePos,iVar,leftOrRight) &
              &  = faceState(:nElems,facePos,iVar,leftOrRight) &
              &  + volState(:nElems,pos,iVar)
          end do

        end do

      end do

    elseif(faceDir.eq.q__T) then

      do facepos = 1,mpd1_square

        do iAnsZ = 1, maxPolyDegree+1
          ! get position of the current ansatz function
          pos = facepos + (iAnsZ-1)*mpd1_square

          do iVar = 1, nScalars
            faceState(:nElems,facePos,iVar,leftOrRight) &
              &  = faceState(:nElems,facePos,iVar,leftOrRight) &
              &  + volState(:nElems,pos,iVar)
          end do

        end do

      end do

    elseif(faceDir.eq.q__B) then

      do facepos = 1,mpd1_square

        do iAnsZ = 1, maxPolyDegree+1
          ! get position of the current ansatz function
          pos = facepos + (iAnsZ-1)*mpd1_square
          ! get the face value of the ansatz function with fixed coordinate
          faceVal = ply_faceValLeftBndAns(iAnsZ)

          do iVar = 1, nScalars
            faceState(:nElems,facePos,iVar,leftOrRight) &
              &  = faceState(:nElems,facePos,iVar,leftOrRight) &
              &  + faceVal*volState(:nElems,pos,iVar)
          end do

        end do

      end do

    else
      write(logUnit(1),*) 'ERROR in atl_modg_volToFace_Q:'
      write(logUnit(1),*) 'Unsupported face direction, stopping...'
      call tem_abort()
    end if
  end subroutine atl_modg_volToFace_Q
  ! ****************************************************************************


  ! ****************************************************************************
  !> Project modal representation of gradients of an element to one of its
  !! faces for Q space.
  !!
  !! Project modal representation of gradient of an element onto one of
  !! its faces. Therefore, this function returns the modal representation
  !! of the solution on the face. This function can project onto an arbitrary
  !! face direction.
  subroutine atl_modg_volToFace_grad_Q( nScalars, volState, maxPolyDegree,     &
    &                                   faceDir, nElems, elemLength, faceState )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: nScalars
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink,
    !! e.g. \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of elements
    integer, intent(in) :: nElems
    !> The lenght of an element
    real(kind=rk), intent(in) :: elemLength
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: pos, facePos, iAnsX, iAnsY, iAnsZ, leftOrRight
    integer :: iElem
    integer :: lb_x, lb_y, lb_z
    integer :: ub_x, ub_y, ub_z
    real(kind=rk) :: faceVal
    real(kind=rk), allocatable :: state_gradient(:,:,:), modalCoeffs(:,:)
    integer :: mpd1, mpd1_square
    ! --------------------------------------------------------------------------

    leftOrRight = tem_dirToFace_map(faceDir)


    ! The lower and upper bounds for the indices of the variables on the face.
    ! We append the gradients to the original variables on the face
    lb_x = nScalars + 1
    ub_x = lb_x + nScalars-1
    lb_y = ub_x + 1
    ub_y = lb_y + nScalars-1
    lb_z = ub_y + 1
    ub_z = lb_z + nScalars-1

    ! Allocate memory for the gradient buffer
    allocate( state_gradient((maxPolyDegree+1)**3,nScalars,3) )
    allocate( modalCoeffs((maxPolyDegree+1)**3,nScalars) )

    mpd1 = maxpolydegree+1
    mpd1_square = mpd1**2


    do iElem = 1, nElems

      ! --> modal space
      modalCoeffs(:,:) = volState(iElem,:,:)

      ! Build gradient inside the current element
      call ply_calcDiff_leg( legCoeffs     = modalCoeffs,    &
        &                    legCoeffsDiff = state_gradient, &
        &                    maxPolyDegree = maxPolyDegree,  &
        &                    nVars         = nScalars,       &
        &                    elemLength    = elemLength      )

      ! now, project to a certain face.
      if( faceDir.eq.q__E ) then

        do facepos = 1,mpd1_square
          do iAnsX = 1, maxPolyDegree+1
            ! get position of the current ansatz function
            pos = (facepos-1)*mpd1 + iAnsX

            ! ... x derivatives
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,1)
            ! ... y derivatives
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,2)
            ! ... z derivatives
            faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  = faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,3)

          end do

        end do

      elseif(faceDir.eq.q__W) then

        do facepos = 1,mpd1_square

          do iAnsX = 1, maxPolyDegree+1
            ! get position of the current ansatz function
            pos = (facepos-1)*mpd1 + iAnsX
            ! get the face value of the ansatz function with fixed coordinate
            faceVal = ply_faceValLeftBndAns(iAnsX)

            ! ... x derivatives
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,1)
            ! ... y derivatives
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,2)
            ! ... z derivatives
            faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  = faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,3)

          end do

        end do

      elseif(faceDir.eq.q__S) then

        do facepos = 1,mpd1_square
          iAnsX = mod(facepos-1, mpd1) + 1
          iAnsZ = (facepos-1)/mpd1 + 1

          do iAnsY = 1, maxPolyDegree+1
            ! get position of the current ansatz function
            pos = iAnsX + (iAnsY-1)*mpd1 + (iAnsZ-1)*mpd1_square
            ! get the face value of the ansatz function with fixed coordinate
            faceVal = ply_faceValLeftBndAns(iAnsY)

            ! ... x derivatives
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
                &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
                &  + faceVal*state_gradient(pos,1:nScalars,1)
            ! ... y derivatives
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
                &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
                &  + faceVal*state_gradient(pos,1:nScalars,2)
            ! ... z derivatives
            faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
                &  = faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
                &  + faceVal*state_gradient(pos,1:nScalars,3)

          end do

        end do

      elseif(faceDir.eq.q__N) then

        do facepos = 1,mpd1_square
          iAnsX = mod(facepos-1, mpd1) + 1
          iAnsZ = (facepos-1)/mpd1 + 1

          do iAnsY = 1, maxPolyDegree+1
            ! get position of the current ansatz function
            pos = iAnsX + (iAnsY-1)*mpd1 + (iAnsZ-1)*mpd1_square

            ! ... x derivatives
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,1)
            ! ... y derivatives
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,2)
            ! ... z derivatives
            faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  = faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,3)

          end do

        end do

      elseif(faceDir.eq.q__T) then

        do facepos = 1,mpd1_square

          do iAnsZ = 1, maxPolyDegree+1
            ! get position of the current ansatz function
            pos = facepos + (iAnsZ-1)*mpd1_square

            ! ... x derivatives
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,1)
            ! ... y derivatives
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,2)
            ! ... z derivatives
            faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  = faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,3)

          end do

        end do

      elseif(faceDir.eq.q__B) then

        do facepos = 1,mpd1_square

          do iAnsZ = 1, maxPolyDegree+1
            ! get position of the current ansatz function
            pos = facepos + (iAnsZ-1)*mpd1_square
            ! get the face value of the ansatz function with fixed coordinate
            faceVal = ply_faceValLeftBndAns(iAnsZ)

            ! ... x derivatives
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,1)
            ! ... y derivatives
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,2)
            ! ... z derivatives
            faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  = faceState(iElem,facePos,lb_z:ub_z,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,3)

          end do

        end do

      else
        write(logUnit(1),*) 'ERROR in atl_modg_volToFace_grad_Q:'
        write(logUnit(1),*) 'Unsupported face direction, stopping...'
        call tem_abort()
      end if

    end do

  end subroutine atl_modg_volToFace_grad_Q
  ! ****************************************************************************


  ! ****************************************************************************
  !> Project modal representation of an element to one of its faces for Q space.
  !!
  !! Project modal representation of an element onto one of its faces.
  !! Therefore, this function returns the modal representation of the solution
  !! on the face. This function can project onto an arbitrary face direction.
  subroutine atl_modg_2d_volToFace_Q( volState, maxPolyDegree, faceDir, &
    &                                 nScalars, nELems, faceState       )
    ! --------------------------------------------------------------------------
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink,
    !! e.g. \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: iElem, lb, ub
    integer :: pos, facePos, iAnsX, iAnsY, leftOrRight
    real(kind=rk) :: faceVal
    ! --------------------------------------------------------------------------

    leftOrRight = tem_dirToFace_map(faceDir)

    lb = 1
    ub = lb + nScalars-1

    ! now, project to a certain face.
    if( faceDir.eq.q__E ) then
      do iAnsY = 1, maxPolyDegree+1
        ! get the position of the modal dof in 2D
  facepos = iansy                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        do iAnsX = 1, maxPolyDegree+1
          ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iElem=1,nElems
            faceState(iElem,facePos,lb:ub,leftOrRight) &
              &  = faceState(iElem,facePos,lb:ub,leftOrRight) &
              &  + volState(iElem,pos,1:nScalars)
          end do
        end do
      end do
    elseif(faceDir.eq.q__W) then
      do iAnsY = 1, maxPolyDegree+1
        ! get the position of the modal dof in 2D
  facepos = iansy                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        do iAnsX = 1, maxPolyDegree+1
          ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          ! get the face value of the ansatz function with fixed coordinate
          faceVal = ply_faceValLeftBndAns(iAnsX)
          do iElem=1,nElems
            faceState(iElem,facePos,lb:ub,leftOrRight) &
              &  = faceState(iElem,facePos,lb:ub,leftOrRight) &
              &  + faceVal*volState(iElem,pos,1:nScalars)
          end do
        end do
      end do
    elseif(faceDir.eq.q__S) then
      do iAnsX = 1, maxPolyDegree+1
        ! get the position of the modal dof in 2D
  facepos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        do iAnsY = 1, maxPolyDegree+1
          ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          ! get the face value of the ansatz function with fixed coordinate
          faceVal = ply_faceValLeftBndAns(iAnsY)
          do iElem=1,nElems
            faceState(iElem,facePos,lb:ub,leftOrRight) &
              &  = faceState(iElem,facePos,lb:ub,leftOrRight) &
              &  + faceVal*volState(iElem,pos,1:nScalars)
          end do
        end do
      end do
    elseif(faceDir.eq.q__N) then
      do iAnsX = 1, maxPolyDegree+1
        ! get the position of the modal dof in 2D
  facepos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        do iAnsY = 1, maxPolyDegree+1
          ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iElem=1,nElems
            faceState(iElem,facePos,lb:ub,leftOrRight) &
              &  = faceState(iElem,facePos,lb:ub,leftOrRight) &
              &  + volState(iElem,pos,1:nScalars)
          end do
        end do
      end do
    else
      write(logUnit(1),*) 'ERROR in atl_modg_2d_volToFace_Q:'
      write(logUnit(1),*) 'Unsupported face direction, stopping...'
      call tem_abort()
    end if

  end subroutine atl_modg_2d_volToFace_Q
  ! ****************************************************************************


  ! ****************************************************************************
  !> Project modal representation of gradients of an element to one of its
  !! faces for Q space.
  !!
  !! Project modal representation of gradients of an element onto one of its
  !! faces. Therefore, this function returns the modal representation of the
  !! solution's gradient on the face. This function can project onto an
  !! arbitrary face direction.
  subroutine atl_modg_2d_volToFace_grad_Q( volState, maxPolyDegree, faceDir, &
    &                                      nScalars, nELems, elemLength,     &
    &                                      faceState                         )
    ! --------------------------------------------------------------------------
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink,
    !! e.g. \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The number of elements
    integer, intent(in) :: nElems
    !> The lenght of an element
    real(kind=rk), intent(in) :: elemLength
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: iElem, lb_x, ub_x, lb_y, ub_y
    integer :: pos, facePos, iAnsX, iAnsY, leftOrRight
    real(kind=rk) :: faceVal
    real(kind=rk), allocatable :: state_gradient(:,:,:), modalCoeffs(:,:)
    ! --------------------------------------------------------------------------

    leftOrRight = tem_dirToFace_map(faceDir)


    ! The lower and upper bounds for the indices of the variables on the face.
    ! We append the gradients to the original variables on the face
    lb_x = nScalars + 1
    ub_x = lb_x + nScalars-1
    lb_y = ub_x + 1
    ub_y = lb_y + nScalars-1

    ! Allocate memory for the gradient buffer
    allocate( state_gradient((maxPolyDegree+1)**2,nScalars,2) )
    allocate( modalCoeffs((maxPolyDegree+1)**2,nScalars) )


    elemLoop: do iElem=1,nElems

      ! --> modal space
      modalCoeffs(:,:) = volState(iElem,:,:)

      ! Build gradient inside the current element
      call ply_calcDiff_leg_2d( legCoeffs     = modalCoeffs,    &
        &                       legCoeffsDiff = state_gradient, &
        &                       maxPolyDegree = maxPolyDegree,  &
        &                       nVars         = nScalars,       &
        &                       elemLength    = elemLength      )

      ! now, project to a certain face.
      if( faceDir.eq.q__E ) then
        do iAnsY = 1, maxPolyDegree+1
          ! get the position of the modal dof in 2D
  facepos = iansy                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iAnsX = 1, maxPolyDegree+1
            ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

            ! Derivatives in x direction
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,1)
            ! Derivatives in y direction
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,2)

          end do
        end do
      elseif(faceDir.eq.q__W) then
        do iAnsY = 1, maxPolyDegree+1
          ! get the position of the modal dof in 2D
  facepos = iansy                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iAnsX = 1, maxPolyDegree+1
            ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
            ! get the face value of the ansatz function with fixed coordinate
            faceVal = ply_faceValLeftBndAns(iAnsX)

            ! Derivatives in x direction
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,1)
            ! Derivatives in y direction
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,2)
          end do
        end do
      elseif(faceDir.eq.q__S) then
        do iAnsX = 1, maxPolyDegree+1
          ! get the position of the modal dof in 2D
  facepos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iAnsY = 1, maxPolyDegree+1
            ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
            ! get the face value of the ansatz function with fixed coordinate
            faceVal = ply_faceValLeftBndAns(iAnsY)

            ! Derivatives in x direction
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,1)
            ! Derivatives in y direction
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + faceVal*state_gradient(pos,1:nScalars,2)
          end do
        end do
      elseif(faceDir.eq.q__N) then
        do iAnsX = 1, maxPolyDegree+1
          ! get the position of the modal dof in 2D
  facepos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          do iAnsY = 1, maxPolyDegree+1
            ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( iansy-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)

            ! Derivatives in x direction
            faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  = faceState(iElem,facePos,lb_x:ub_x,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,1)
            ! Derivatives in y direction
            faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  = faceState(iElem,facePos,lb_y:ub_y,leftOrRight) &
              &  + state_gradient(pos,1:nScalars,2)
          end do
        end do
      else
        write(logUnit(1),*) 'ERROR in atl_modg_2d_volToFace_grad_Q:'
        write(logUnit(1),*) 'Unsupported face direction, stopping...'
        call tem_abort()
      end if

    end do elemLoop


  end subroutine atl_modg_2d_volToFace_grad_Q
  ! ****************************************************************************

  ! ****************************************************************************
  !> Project modal representation of an element to one of its faces for Q space.
  !!
  !! Project modal representation of an element onto one of its faces.
  !! Therefore, this function returns the modal representation of the solution
  !! on the face. This function can project onto an arbitrary face direction.
  subroutine atl_modg_1d_volToFace_Q( volState, maxPolyDegree, faceDir, &
    &                                 nScalars, nELems, faceState       )
    ! --------------------------------------------------------------------------
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink,
    !! e.g. \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iVar
    integer :: pos, iAnsX,leftOrRight
    real(kind=rk) :: faceVal
    ! --------------------------------------------------------------------------

    leftOrRight = tem_dirToFace_map(faceDir)

    ! now, project to a certain face.
    if( faceDir.eq.q__E ) then
      do iAnsX = 1, maxPolyDegree+1
        ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        do iElem=1,nElems
          do iVar = 1, nScalars
            faceState(iElem,1,iVar,leftOrRight) &
              &  = faceState(iElem,1,iVar,leftOrRight) &
              &  + volState(iElem,pos,iVar)
          end do
        end do
      end do
    elseif(faceDir.eq.q__W) then
      do iAnsX = 1, maxPolyDegree+1
        ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
        ! get the face value of the ansatz function with fixed coordinate
        faceVal = ply_faceValLeftBndAns(iAnsX)
        do iElem=1,nElems
          do iVar = 1, nScalars
            faceState(iElem,1,iVar,leftOrRight) &
              &  = faceState(iElem,1,iVar,leftOrRight) &
              &  + faceVal*volState(iElem,pos,iVar)
          end do
        end do
      end do
    else
      write(logUnit(1),*) 'ERROR in atl_modg_1d_volToFace_Q:'
      write(logUnit(1),*) 'Unsupported face direction, stopping...'
      call tem_abort()
    end if

  end subroutine atl_modg_1d_volToFace_Q
  ! ****************************************************************************


  ! ****************************************************************************
  !> Projects modal representation of each cell to its faces, i.e.
  !! this subroutine creates a modal representation on the faces.
  subroutine atl_modg_1d_modalVolToModalFace( mesh, statedata, facedata, &
    &                                         nScalars, maxPolyDegree,   &
    &                                         basisType, equation        )
    ! --------------------------------------------------------------------------
    !> The elements we apply the projection for.
    type(atl_cube_elem_type),  intent(in) :: mesh
    !> Volumetric, modal states for each element.
    type(atl_statedata_type), intent(in) :: statedata
    !> Modal representation on the face (will be updated by this routine for all
    !! fluid elements in mesh).
    type(atl_facedata_type), intent(inout) :: facedata
    !> The number of scalars varaibales in your equation system.
    integer, intent(in) :: nScalars
    !> The parameters of your modg scheme.
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: basisType
    !> The equation you solve.
    type(atl_equations_type)          :: equation
    ! --------------------------------------------------------------------------
    integer :: spaceDir, iDir
    ! --------------------------------------------------------------------------

    ! Iterate over all the fluid elements and project its modal representations
    ! to the faces of the element.

    select case(basisType)
      case(Q_space) ! Q tensor product ansatz functions

        do iDir = 1, equation%nDimensions
          facedata%faceRep(iDir)                                             &
            &     %dat(1:mesh%descriptor%elem%nElems(eT_fluid),:,:,:) = 0.0_rk
        end do

        ! Now, iterate over all the faces and project to this face.
        ! ... x faces
        ! Project to the face in the current direction.
        spaceDir = qAxis(1)
        call atl_modg_1d_volToFace_Q(                              &
          & volState      = statedata%state,                       &
          & maxPolyDegree = maxPolyDegree,                         &
          & faceDir       = 1,                                     &
          & nScalars      = nScalars,                              &
          & nElems        = mesh%descriptor%elem%nElems(eT_fluid), &
          & faceState     = facedata%faceRep(spaceDir)%dat         )
        ! ... x
        ! Project to the face in the current direction.
        spaceDir = qAxis(4)
        call atl_modg_1d_volToFace_Q(                              &
          & volState      = statedata%state,                       &
          & maxPolyDegree = maxPolyDegree,                         &
          & faceDir       = 4,                                     &
          & nScalars      = nScalars,                              &
          & nElems        = mesh%descriptor%elem%nElems(eT_fluid), &
          & faceState     = facedata%faceRep(spaceDir)%dat         )


        ! Now, iterate over all the faces and project the gradients
        ! to this face.
        if(equation%nDerivatives == 1 ) then
          spaceDir = qAxis(1)

          call atl_modg_1d_VolToFace_grad_Q(                         &
            & volState      = statedata%state,                       &
            & maxPolyDegree = maxPolyDegree,                         &
            & faceDir       = 1,                                     &
            & nScalars      = nScalars,                              &
            & nElems        = mesh%descriptor%elem%nElems(eT_fluid), &
            & elemLength    = mesh%length,                           &
            & faceState     = facedata%faceRep(spaceDir)%dat         )

          spaceDir = qAxis(4)

          call atl_modg_1d_VolToFace_grad_Q(                         &
            & volState      = statedata%state,                       &
            & maxPolyDegree = maxPolyDegree,                         &
            & faceDir       = 4,                                     &
            & nScalars      = nScalars,                              &
            & nElems        = mesh%descriptor%elem%nElems(eT_fluid), &
            & elemLength    = mesh%length,                           &
            & faceState     = facedata%faceRep(spaceDir)%dat         )
          end if
      case default
        call tem_abort( 'ERROR in atl_modg_1d_modalVolToModalFace: ' &
          & // 'Unknown tensor product'                              )
    end select
  end subroutine atl_modg_1d_modalVolToModalFace
  ! ****************************************************************************


  ! ****************************************************************************
  subroutine atl_modg_1d_VolToFace_grad_Q( volstate, maxPolyDegree, faceDir, &
    &                                      nScalars, nElems,elemLength,      &
    &                                      faceState                         )
    ! --------------------------------------------------------------------------
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink,
    !! e.g. \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The number of elements
    integer, intent(in) :: nElems
    !> Length of elements
    real(kind=rk), intent(in) :: elemLength
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: iElem, lb, ub
    integer :: pos, iAnsX,leftOrRight
    real(kind=rk) :: faceVal
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    real(kind=rk), allocatable :: state_gradient(:,:)
    ! --------------------------------------------------------------------------

    ! Allocate memory for the gradient buffer
    allocate( state_gradient(maxPolyDegree+1,nScalars) )
    allocate( modalCoeffs(maxPolyDegree+1,nScalars) )

    lb = nScalars + 1
    ub = lb + nScalars -1

    leftOrRight = tem_dirToFace_map(faceDir)

    do iElem = 1, nElems

      modalCoeffs(:,:) = volState(iElem,:,:)

      ! Now we calculate the gradient of the modal representation required
      call ply_calcDiff_leg_1d( legCoeffs     = modalCoeffs,    &
        &                       legcoeffsDiff = state_gradient, &
        &                       maxPolyDegree = maxPolyDegree,  &
        &                       elemLength    = elemLength      )

        ! now, project to a certain face.
      if( faceDir.eq.q__E ) then
        do iAnsX = 1, maxPolyDegree+1
          ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          faceState(iElem,1,lb:ub,leftOrRight) &
            &  = faceState(iElem,1,lb:ub,leftOrRight) &
            &  + state_gradient(pos,1:nScalars)
        end do
      elseif(faceDir.eq.q__W) then
        do iAnsX = 1, maxPolyDegree+1
          ! get position of the current ansatz function
  pos = iansx                                      &
    &      + ( ( 1-1)                             &
    &      + (1-1)*(maxpolydegree+1))*(maxpolydegree+1)
          ! get the face value of the ansatz function with fixed coordinate
          faceVal = ply_faceValLeftBndAns(iAnsX)
          faceState(iElem,1,lb:ub,leftOrRight) &
            &  = faceState(iElem,1,lb:ub,leftOrRight) &
            &  + faceVal*state_gradient(pos,1:nScalars)
        end do
      else
        write(logUnit(1),*) 'ERROR in atl_modg_1d_volToFace_grad_Q:'
        write(logUnit(1),*) 'Unsupported face direction, stopping...'
        call tem_abort()
      end if

    end do  ! elems loop

  end subroutine atl_modg_1d_VolToFace_grad_Q
  ! ****************************************************************************


  !> Project modal representation of an element to one of its faces for P space.
  !!
  !! Project modal representation of an element onto one of its faces.
  !! Therefore, this function returns the modal representation of the solution
  !! on the face. This function can project onto an arbitrary face direction.
  subroutine atl_modg_volToFace_P( nTotalElems, nTotalFaces, nDofs, nFaceDofs, &
      &                            nScalars, volState, maxPolyDegree, faceDir, &
      &                            nELems, faceState                           )
    ! --------------------------------------------------------------------------
    !> dimensions
    integer, intent(in) :: nTotalElems, nDofs, nScalars, nTotalFaces, nFaceDofs
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(nTotalElems,nDofs,nScalars)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink,
    !! e.g. \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(nTotalFaces,nFaceDofs,nScalars,2)
    ! --------------------------------------------------------------------------
    integer :: pos, posMax, facePos, iAnsX, iAnsY, iAnsZ, leftOrRight
    real(kind=rk) :: faceVal
    integer :: iVar
    ! --------------------------------------------------------------------------

    leftOrRight = tem_dirToFace_map(faceDir)

    ! now, project to a certain face.
    if( faceDir.eq.q__E ) then

      iAnsX = 1
      iAnsY = 1
      iAnsZ = 1
  posmax = (((maxpolydegree) + 1) &
    &   * ((maxpolydegree) + 2) &
    &   * ((maxpolydegree) + 3)) &
    & / 6
      do pos = 1, posMax
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = ((((iansy - 1) + (iansz - 1))            &
    &   * (((iansy - 1) + (iansz - 1)) + 1)) / 2 + 1) &
    & + (iansz - 1)

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  elseif (iansy .ne. 1) then
    ! next block
    iansx = iansy - 1
    iansy = 1
    iansz = iansz + 1
  else
    ! next layer
    iansx = iansz + 1
    iansy = 1
    iansz = 1
  end if
      end do

    elseif(faceDir.eq.q__W) then

      iAnsX = 1
      iAnsY = 1
      iAnsZ = 1
  posmax = (((maxpolydegree) + 1) &
    &   * ((maxpolydegree) + 2) &
    &   * ((maxpolydegree) + 3)) &
    & / 6
      do pos = 1, posMax
        ! get the face value of the ansatz function with fixed coordinate
        faceVal = ply_faceValLeftBndAns(iAnsX)
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = ((((iansy - 1) + (iansz - 1))            &
    &   * (((iansy - 1) + (iansz - 1)) + 1)) / 2 + 1) &
    & + (iansz - 1)

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + faceVal*volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  elseif (iansy .ne. 1) then
    ! next block
    iansx = iansy - 1
    iansy = 1
    iansz = iansz + 1
  else
    ! next layer
    iansx = iansz + 1
    iansy = 1
    iansz = 1
  end if
      end do

    elseif(faceDir.eq.q__S) then

      iAnsX = 1
      iAnsY = 1
      iAnsZ = 1
  posmax = (((maxpolydegree) + 1) &
    &   * ((maxpolydegree) + 2) &
    &   * ((maxpolydegree) + 3)) &
    & / 6
      do pos = 1, posMax
        ! get the face value of the ansatz function with fixed coordinate
        faceVal = ply_faceValLeftBndAns(iAnsY)
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = ((((iansx - 1) + (iansz - 1))            &
    &   * (((iansx - 1) + (iansz - 1)) + 1)) / 2 + 1) &
    & + (iansz - 1)

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + faceVal*volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  elseif (iansy .ne. 1) then
    ! next block
    iansx = iansy - 1
    iansy = 1
    iansz = iansz + 1
  else
    ! next layer
    iansx = iansz + 1
    iansy = 1
    iansz = 1
  end if
      end do

    elseif(faceDir.eq.q__N) then

      iAnsX = 1
      iAnsY = 1
      iAnsZ = 1
  posmax = (((maxpolydegree) + 1) &
    &   * ((maxpolydegree) + 2) &
    &   * ((maxpolydegree) + 3)) &
    & / 6
      do pos = 1, posMax
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = ((((iansx - 1) + (iansz - 1))            &
    &   * (((iansx - 1) + (iansz - 1)) + 1)) / 2 + 1) &
    & + (iansz - 1)

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  elseif (iansy .ne. 1) then
    ! next block
    iansx = iansy - 1
    iansy = 1
    iansz = iansz + 1
  else
    ! next layer
    iansx = iansz + 1
    iansy = 1
    iansz = 1
  end if
      end do

    elseif(faceDir.eq.q__T) then

      iAnsX = 1
      iAnsY = 1
      iAnsZ = 1
  posmax = (((maxpolydegree) + 1) &
    &   * ((maxpolydegree) + 2) &
    &   * ((maxpolydegree) + 3)) &
    & / 6
      do pos = 1, posMax
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = ((((iansx - 1) + (iansy - 1))            &
    &   * (((iansx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  elseif (iansy .ne. 1) then
    ! next block
    iansx = iansy - 1
    iansy = 1
    iansz = iansz + 1
  else
    ! next layer
    iansx = iansz + 1
    iansy = 1
    iansz = 1
  end if
      end do

    elseif(faceDir.eq.q__B) then

      iAnsX = 1
      iAnsY = 1
      iAnsZ = 1
  posmax = (((maxpolydegree) + 1) &
    &   * ((maxpolydegree) + 2) &
    &   * ((maxpolydegree) + 3)) &
    & / 6
      do pos = 1, posMax
        ! get the face value of the ansatz function with fixed coordinate
        faceVal = ply_faceValLeftBndAns(iAnsZ)
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = ((((iansx - 1) + (iansy - 1))            &
    &   * (((iansx - 1) + (iansy - 1)) + 1)) / 2 + 1) &
    & + (iansy - 1)

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + faceVal*volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  elseif (iansy .ne. 1) then
    ! next block
    iansx = iansy - 1
    iansy = 1
    iansz = iansz + 1
  else
    ! next layer
    iansx = iansz + 1
    iansy = 1
    iansz = 1
  end if
      end do
    else
      write(logUnit(1),*) 'ERROR in atl_modg_volToFace_P:'
      write(logUnit(1),*) 'Unsupported face direction, stopping...'
      call tem_abort()
    end if

  end subroutine atl_modg_volToFace_P


  ! ****************************************************************************
  !> Project modal representation of an element to one of its faces for P space.
  !!
  !! Project modal representation of an element onto one of its faces.
  !! Therefore, this function returns the modal representation of the solution
  !! on the face. This function can project onto an arbitrary face direction.
  subroutine atl_modg_2d_volToFace_P( volState, maxPolyDegree, faceDir, &
    &                                 nScalars, nELems, faceState       )
    ! --------------------------------------------------------------------------
    !> The modal representation in the volume. First dimension is the number of
    !! voluemtrix numbers of degrees of freedom and second dimension is the
    !! number of scalar variables in the equation system.
    real(kind=rk), intent(in) :: volState(:,:,:)
    !> The maximal polynomial degree per spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The face to project the modal representation to.
    !! Use one of the first six directions of \link tem_param_module \endlink,
    !! e.g. \link tem_param_module::q__e \endlink
    integer, intent(in) :: faceDir
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The number of elements
    integer, intent(in) :: nElems
    !> The modal representation on the face
    real(kind=rk), intent(inout) :: faceState(:,:,:,:)
    ! --------------------------------------------------------------------------
    integer :: pos, posMax, facePos
    integer :: iAnsX, iAnsY, iVar
    integer :: leftOrRight
    real(kind=rk) :: faceVal
    ! --------------------------------------------------------------------------

    leftOrRight = tem_dirToFace_map(faceDir)

    ! now, project to a certain face.
    if( faceDir.eq.q__E ) then

      iAnsX = 1
      iAnsY = 1
  posmax = ((maxpolydegree)+1)*((maxpolydegree)+2)/2
      do pos = 1, posMax
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = iansy

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  else
    ! next layer
    iansx = iansy + 1
    iansy = 1
  end if
      end do

    elseif(faceDir.eq.q__W) then

      iAnsX = 1
      iAnsY = 1
  posmax = ((maxpolydegree)+1)*((maxpolydegree)+2)/2
      do pos = 1, posMax
        ! get the face value of the ansatz function with fixed coordinate
        faceVal = ply_faceValLeftBndAns(iAnsX)
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = iansy

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + faceVal*volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  else
    ! next layer
    iansx = iansy + 1
    iansy = 1
  end if
      end do

    elseif(faceDir.eq.q__S) then

      iAnsX = 1
      iAnsY = 1
  posmax = ((maxpolydegree)+1)*((maxpolydegree)+2)/2
      do pos = 1, posMax
        ! get the face value of the ansatz function with fixed coordinate
        faceVal = ply_faceValLeftBndAns(iAnsY)
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = iansx

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + faceVal*volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  else
    ! next layer
    iansx = iansy + 1
    iansy = 1
  end if
      end do

    elseif(faceDir.eq.q__N) then

      iAnsX = 1
      iAnsY = 1
  posmax = ((maxpolydegree)+1)*((maxpolydegree)+2)/2
      do pos = 1, posMax
        ! get the position of the modal dof in 2D
  ! integer divisions are no mistake here.
  facepos = iansx

        do iVar=1,nScalars
          faceState(:nElems,facePos,iVar,leftOrRight) &
            &  = faceState(:nElems,facePos,iVar,leftOrRight) &
            &  + volState(:nElems,pos,iVar)
        end do

  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (iansx .ne. 1) then
    ! next item
    iansx = iansx - 1
    iansy = iansy + 1
  else
    ! next layer
    iansx = iansy + 1
    iansy = 1
  end if
      end do

    else
      write(logUnit(1),*) 'ERROR in atl_modg_2d_volToFace_P:'
      write(logUnit(1),*) 'Unsupported face direction, stopping...'
      call tem_abort()
    end if

  end subroutine atl_modg_2d_volToFace_P
  ! ****************************************************************************


end module atl_volToFace_module
