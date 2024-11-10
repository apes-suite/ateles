! Copyright (c) 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
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
module atl_subresolution_module
  use env_module, only: pathLen, rk, isLittleEndian, long_k, newunit

  use aotus_module, only: flu_State, aot_get_val, aoterr_Fatal, close_config

  use treelmesh_module,       only: treelmesh_type
  use tem_aux_module,         only: tem_open_distconf, tem_abort
  use tem_color_prop_module,  only: tem_color_prop_type
  use tem_comm_env_module,    only: tem_comm_env_type
  use tem_logging_module,     only: logunit
  use tem_subres_prop_module, only: tem_subres_prop_type, tem_subres_prop_load
  use tem_tools_module,       only: upper_to_lower

  use ply_dof_module, only: P_Space, Q_space

  implicit none

  private

  type atl_subresolution_type
    integer :: polydegree
    integer :: basisType
    type(tem_subres_prop_type) :: subres_prop
  end type atl_subresolution_type

  public :: atl_subresolution_type, atl_subresolution_load
  public :: atl_subres_import_color


contains


  ! ****************************************************************************!
  !> Subroutine to load subresolution information for a given tree.
  subroutine atl_subresolution_load(me, tree, proc, coloring)
    ! --------------------------------------------------------------------------!
    type(atl_subresolution_type), intent(out) :: me
    type(treelmesh_type), intent(in) :: tree
    type(tem_comm_env_type), intent(in) :: proc
    type(tem_color_prop_type), intent(in) :: coloring
    ! --------------------------------------------------------------------------!
    character(len=pathLen) :: configfile
    character :: polyspace
    type(flu_State) :: conf
    integer :: iError
    ! --------------------------------------------------------------------------!

    configfile = trim(tree%global%dirname)//'subresolution.lua'

    call tem_subres_prop_load( me       = me%subres_prop, &
      &                        tree     = tree,           &
      &                        coloring = coloring        )

    call tem_open_distconf( L        = conf,             &
      &                     filename = trim(configfile), &
      &                     proc     = proc              )

    call aot_get_val( L       = conf,          &
      &               key     = 'polydegree',  &
      &               val     = me%polydegree, &
      &               ErrCode = iError         )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*) &
        &  'FATAL Error occured, while retrieving subresolution polydegree'
      call tem_abort()
    end if

    call aot_get_val( L       = conf,        &
      &               key     = 'polyspace', &
      &               val     = polyspace,   &
      &               default = 'q',         &
      &               ErrCode = iError       )

    polyspace = upper_to_lower(polyspace)
    select case(polyspace)
    case('q')
      me%basisType = Q_space

    case('p')
      me%basisType = P_space

    case default
      write(logunit(1),*) 'ERROR in subresolution loading!'
      write(logunit(1),*) 'Unknown polyspace ', trim(polyspace)
      write(logUnit(1),*) 'Supported are:'
      write(logUnit(1),*) '* Q (quadratic with i,j,k <= maxDegree)'
      write(logUnit(1),*) '* P (with i+j+k <= maxDegree)'
      write(logUnit(1),*) 'Stopping....'
      call tem_abort()

    end select

    call close_config(conf)

  end subroutine atl_subresolution_load
  ! ****************************************************************************!


  ! ****************************************************************************!
  !> Get the subresolution data for all elements for a given color and in the
  !! requested format.
  subroutine atl_subres_import_color( me, tree, coloring, iColor,              &
    &                                 target_degree, target_space, target_dim, &
    &                                 subresdat )
    ! --------------------------------------------------------------------------!
    type(atl_subresolution_type), intent(in) :: me
    type(treelmesh_type), intent(in) :: tree
    type(tem_color_prop_type), intent(in) :: coloring
    integer, intent(in) :: iColor
    integer, intent(in) :: target_degree
    integer, intent(in) :: target_space
    integer, intent(in) :: target_dim
    real(kind=rk), allocatable, intent(out) :: subresdat(:,:)
    ! --------------------------------------------------------------------------!
    character(len=pathLen) :: datfile
    integer :: target_Dofs
    integer :: in_dofs
    integer :: min_dofs
    integer :: pdim
    real(kind=rk), allocatable :: indat(:)
    character(len=4) :: datext
    integer :: rl
    integer :: fUnit
    integer :: nElems
    integer(kind=long_k) :: offset
    integer :: iElem
    integer :: iStep, iDof
    integer :: tt_X, tt_Y, tt_Z
    integer :: in_X, in_Y, in_Z
    integer :: tt_pos, in_pos
    integer :: tt_off, tt_zoff
    integer :: in_off, in_zoff
    integer :: minOrd
    ! --------------------------------------------------------------------------!


    ! Only need to do anything, if there is actually a subresolution property.
    ! Checking this via the association status of the property header.
    if (associated(me%subres_prop%header)) then

      if (isLittleEndian) then
        datext = '.lsb'
      else
        datext = '.msb'
      end if

      nElems = me%subres_prop%nElems(iColor)
      offset = me%subres_prop%offset(iColor)

      ! Figure out the target polynomial representation.
      select case(target_space)
      case (Q_Space)
        target_dofs = (target_degree+1)**target_dim

      case (P_Space)
        target_dofs = target_degree+1
        do pdim=2,target_dim
          target_dofs = (target_dofs * (target_degree+pdim) ) / pdim
        end do

      end select

      ! Figure out input polynomial representation.
      select case(me%basistype)
      case (Q_Space)
        in_dofs = (me%polydegree+1)**3

      case (P_Space)
        in_dofs = me%polydegree+1
        do pdim=2,3
          in_dofs = (in_dofs * (me%polydegree+pdim) ) / pdim
        end do

      end select

      ! Allocate arrays accordingly.
      allocate(indat(in_dofs))
      allocate(subresdat(target_dofs, nElems))

      inquire(iolength=rl) indat

      datfile = trim(tree%global%dirname) // 'subresdata_' &
        &       // trim(coloring%color_label(iColor)) // datext

      fUnit = newunit()

      open( unit = fUnit, file = datfile, action = 'read', &
        &   access = 'direct', form = 'unformatted',       &
        &   recl = rl, status = 'old'                      )

      minord = min(target_degree+1, me%polydegree+1)

      subresdat = 0.0_rk

      do iElem=1,nElems

        read(fUnit, rec=offset+iElem) indat

        select case(target_space)
        case (Q_Space)
          if (me%basistype == Q_Space) then
            ! Both, target and input are Q Polynomials
            select case(target_dim)
            case (1)
              subresdat(:minord, iElem) = indat(:minord)

            case (2)
              do tt_Y=0,minord-1
                tt_off = tt_Y*(target_degree+1)
                in_off = tt_Y*(me%polydegree+1)
                subresdat(tt_off+1:tt_off+minord, iElem) &
                  &  = indat(in_off+1:in_off+minord)
              end do

            case (3)
              do tt_Z=0,minord-1
                tt_zoff = tt_Z*(target_degree+1)**2
                in_zoff = tt_Z*(me%polydegree+1)**2
                do tt_Y=0,minord-1
                  tt_off = tt_Y*(target_degree+1) + tt_zoff
                  in_off = tt_Y*(me%polydegree+1) + in_zoff
                  subresdat(tt_off+1:tt_off+minord, iElem) &
                    &  = indat(in_off+1:in_off+minord)
                end do
              end do
            end select

          else

            ! Target is Q, but input is P
            select case(target_dim)
            case (1)
              subresdat(:minord, iElem) = indat(:minord)

            case (2)
              in_X = 1
              in_Y = 1
              do iDof=1,in_dofs
  ! integer divisions are no mistake here.
  in_pos = ((((in_x - 1) + (in_y - 1))            &
    &   * (((in_x - 1) + (in_y - 1)) + 1)) / 2 + 1) &
    & + (in_y - 1)
                tt_pos = in_X + (target_degree+1)*(in_Y-1)
                subresdat(tt_pos, iElem) = indat(in_pos)
                ! Ensure, that next iteration is in the target range
                do istep=iDof,in_dofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (in_x .ne. 1) then
    ! next item
    in_x = in_x - 1
    in_y = in_y + 1
  else
    ! next layer
    in_x = in_y + 1
    in_y = 1
  end if
                  if ((in_X <= minord) .and. (in_Y <= minord)) EXIT
                end do
                if ((in_X > minord) .or. (in_Y > minord)) EXIT
              end do

            case (3)
              in_X = 1
              in_Y = 1
              in_Z = 1
              do iDof=1,in_dofs
  ! integer divisions are no mistake here.
  in_pos = (((in_x + in_y + in_z - 3) &
    &     * (in_x + in_y + in_z - 2) &
    &     * (in_x + in_y + in_z - 1)) &
    &   / 6 + 1)             &
    & + ((in_z-1) * (in_x + in_y + in_z -2) &
    &   - ((in_z-2) * (in_z-1)) / 2) &
    & + (in_y-1)
                tt_pos = in_X + (target_degree+1)*( (in_Y-1) &
                  &                                + (target_degree+1)*(in_Z-1))
                subresdat(tt_pos, iElem) = indat(in_pos)
                ! Ensure, that next iteration is in the target range
                do istep=iDof,in_dofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (in_x .ne. 1) then
    ! next item
    in_x = in_x - 1
    in_y = in_y + 1
  elseif (in_y .ne. 1) then
    ! next block
    in_x = in_y - 1
    in_y = 1
    in_z = in_z + 1
  else
    ! next layer
    in_x = in_z + 1
    in_y = 1
    in_z = 1
  end if
                  if ( (in_X <= minord) .and. (in_Y <= minord) &
                    &                   .and. (in_Z <= minord) ) EXIT
                end do
                if ( (in_X > minord) .or. (in_Y > minord) &
                  &                  .or. (in_Z > minord) ) EXIT
              end do
            end select
          end if

        case (P_Space)
          if (me%basistype == Q_Space) then
            ! Target is P, input is Q
            select case(target_dim)
            case (1)
              subresdat(:minord, iElem) = indat(:minord)

            case (2)
              tt_X = 1
              tt_Y = 1
              do iDof=1,target_dofs
  ! integer divisions are no mistake here.
  tt_pos = ((((tt_x - 1) + (tt_y - 1))            &
    &   * (((tt_x - 1) + (tt_y - 1)) + 1)) / 2 + 1) &
    & + (tt_y - 1)
                in_pos = tt_X + (me%polydegree+1)*(tt_Y-1)
                subresdat(tt_pos, iElem) = indat(in_pos)
                ! Ensure, that next iteration is in the target range
                do istep=iDof,target_dofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (tt_x .ne. 1) then
    ! next item
    tt_x = tt_x - 1
    tt_y = tt_y + 1
  else
    ! next layer
    tt_x = tt_y + 1
    tt_y = 1
  end if
                  if ((tt_X <= minord) .and. (tt_Y <= minord)) EXIT
                end do
                if ((tt_X > minord) .or. (tt_Y > minord)) EXIT
              end do

            case (3)
              tt_X = 1
              tt_Y = 1
              tt_Z = 1
              do iDof=1,target_dofs
  ! integer divisions are no mistake here.
  tt_pos = (((tt_x + tt_y + tt_z - 3) &
    &     * (tt_x + tt_y + tt_z - 2) &
    &     * (tt_x + tt_y + tt_z - 1)) &
    &   / 6 + 1)             &
    & + ((tt_z-1) * (tt_x + tt_y + tt_z -2) &
    &   - ((tt_z-2) * (tt_z-1)) / 2) &
    & + (tt_y-1)
                in_pos = tt_X + (me%polydegree+1)*( (tt_Y-1) &
                  &                                + (me%polydegree+1)*(tt_Z-1))
                subresdat(tt_pos, iElem) = indat(in_pos)
                ! Ensure, that next iteration is in the target range
                do istep=iDof,target_dofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (tt_x .ne. 1) then
    ! next item
    tt_x = tt_x - 1
    tt_y = tt_y + 1
  elseif (tt_y .ne. 1) then
    ! next block
    tt_x = tt_y - 1
    tt_y = 1
    tt_z = tt_z + 1
  else
    ! next layer
    tt_x = tt_z + 1
    tt_y = 1
    tt_z = 1
  end if
                  if ((tt_X <= minord) .and. (tt_Y <= minord) &
                    &                  .and. (tt_Z <= minord) ) EXIT
                end do
                if ((tt_X > minord) .or. (tt_Y > minord) &
                  &                 .or. (tt_Z > minord) ) EXIT
              end do

            end select

          else
            ! Both input and target ar P polynomials
            select case(target_dim)
            case (1)
              subresdat(:minord, iElem) = indat(:minord)

            case (2)
              min_dofs = (minord*(minord+1))/2
              tt_X = 1
              tt_Y = 1
              do iDof=1,min_dofs
  ! integer divisions are no mistake here.
  tt_pos = ((((tt_x - 1) + (tt_y - 1))            &
    &   * (((tt_x - 1) + (tt_y - 1)) + 1)) / 2 + 1) &
    & + (tt_y - 1)
                in_pos = tt_pos
                subresdat(tt_pos, iElem) = indat(in_pos)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (tt_x .ne. 1) then
    ! next item
    tt_x = tt_x - 1
    tt_y = tt_y + 1
  else
    ! next layer
    tt_x = tt_y + 1
    tt_y = 1
  end if
              end do

            case (3)
              min_dofs = ( (minord+2)*((minord*(minord+1))/2) ) / 3
              tt_X = 1
              tt_Y = 1
              tt_Z = 1
              do iDof=1,min_dofs
  ! integer divisions are no mistake here.
  tt_pos = (((tt_x + tt_y + tt_z - 3) &
    &     * (tt_x + tt_y + tt_z - 2) &
    &     * (tt_x + tt_y + tt_z - 1)) &
    &   / 6 + 1)             &
    & + ((tt_z-1) * (tt_x + tt_y + tt_z -2) &
    &   - ((tt_z-2) * (tt_z-1)) / 2) &
    & + (tt_y-1)
  ! integer divisions are no mistake here.
  in_pos = (((tt_x + tt_y + tt_z - 3) &
    &     * (tt_x + tt_y + tt_z - 2) &
    &     * (tt_x + tt_y + tt_z - 1)) &
    &   / 6 + 1)             &
    & + ((tt_z-1) * (tt_x + tt_y + tt_z -2) &
    &   - ((tt_z-2) * (tt_z-1)) / 2) &
    & + (tt_y-1)
                subresdat(tt_pos, iElem) = indat(in_pos)
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (tt_x .ne. 1) then
    ! next item
    tt_x = tt_x - 1
    tt_y = tt_y + 1
  elseif (tt_y .ne. 1) then
    ! next block
    tt_x = tt_y - 1
    tt_y = 1
    tt_z = tt_z + 1
  else
    ! next layer
    tt_x = tt_z + 1
    tt_y = 1
    tt_z = 1
  end if
              end do
            end select
          end if

        end select

      end do

      close(fUnit)

    else

      ! No subresolution data at all, allocate the array with size 0.
      allocate(subresdat(0,0))

    end if

  end subroutine atl_subres_import_color


end module atl_subresolution_module
