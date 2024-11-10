! Copyright (c) 2015-2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Peter Vitt, Tobias
! Girresser and Daniel Petró for University of Siegen.
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

module ply_subresolution_module
  use, intrinsic :: iso_c_binding, only: c_f_pointer

  use env_module,                  only: pathLen,        &
    &                                    rk,             &
    &                                    isLittleEndian, &
    &                                    long_k,         &
    &                                    newunit

  use aotus_module,                only: flu_State,    &
    &                                    aot_get_val,  &
    &                                    aoterr_Fatal, &
    &                                    close_config

  use treelmesh_module,            only: treelmesh_type
  use tem_aux_module,              only: tem_open_distconf, &
    &                                    tem_abort
  use tem_color_prop_module,       only: tem_color_prop_type, &
    &                                    colors_per_char
  use tem_comm_env_module,         only: tem_comm_env_type
  use tem_logging_module,          only: logunit
  use tem_subres_prop_module,      only: tem_subres_prop_type, &
    &                                    tem_subres_prop_load
  use tem_tools_module,            only: upper_to_lower
  use tem_time_module,             only: tem_time_type
  use tem_varSys_module,           only: tem_varSys_type, &
    &                                    tem_varSys_op_type

  use ply_dof_module,              only: P_Space, &
    &                                    Q_space
  use ply_transfer_module,         only: ply_transfer_P_dim, &
    &                                    ply_transfer_dofs

  implicit none

  private

  type ply_subresolution_type
    integer :: polydegree = 0
    integer :: basisType
    type(tem_subres_prop_type) :: subres_prop
  end type ply_subresolution_type

  !> Self contained description of color data to be used for method data.
  type ply_subres_colvar_type
    !> Pointer to the overall color property description.
    type(tem_color_prop_type), pointer :: color => NULL()

    !> Pointer to the overall subresolution property description.
    type(ply_subresolution_type), pointer :: subres => NULL()

    !> Degrees of freedom for subresolved elements.
    real(kind=rk), allocatable :: subresdat(:,:)

    !> Number of degrees of freedom of the subresolution data.
    integer :: nsubdofs

    !> Position of this color in the list of colors, needed to make use of the
    !! pointers above.
    integer :: colpos
  end type ply_subres_colvar_type

  public :: ply_subresolution_type
  public :: ply_subresolution_load
  public :: ply_subres_import_color
  public :: ply_subres_get_elemcolor
  public :: ply_subres_colvar_type


contains


  ! ************************************************************************ !
  !> Subroutine to load subresolution information for a given tree.
  subroutine ply_subresolution_load( me, tree, proc, coloring )
    ! --------------------------------------------------------------------- !
    type(ply_subresolution_type), intent(out) :: me
    type(treelmesh_type), intent(in) :: tree
    type(tem_comm_env_type), intent(in) :: proc
    type(tem_color_prop_type), intent(in) :: coloring
    ! --------------------------------------------------------------------- !
    character(len=pathLen) :: configfile
    character :: polyspace
    type(flu_State) :: conf
    integer :: iError
    ! --------------------------------------------------------------------- !

    configfile = trim(tree%global%dirname) // 'subresolution.lua'

    call tem_subres_prop_load( me       = me%subres_prop, &
      &                        tree     = tree,           &
      &                        coloring = coloring        )

    ! Set the polydegree initially to 0 to ensure a proper setting even if
    ! no subresolution property is present.
    me%polydegree = 0

    ! Only need to do anything, if there is actually a subresolution property.
    ! Checking this via the association status of the property header.
    if (associated(me%subres_prop%header)) then

      call tem_open_distconf( L        = conf,             &
        &                     filename = trim(configfile), &
        &                     proc     = proc              )

      call aot_get_val( L       = conf,          &
        &               key     = 'polydegree',  &
        &               val     = me%polydegree, &
        &               ErrCode = iError         )

      if (btest(iError, aoterr_Fatal)) then
        call tem_abort(                                                      &
          & 'FATAL Error occured, while retrieving subresolution polydegree' )
      end if

      call aot_get_val( L       = conf,        &
        &               key     = 'polyspace', &
        &               val     = polyspace,   &
        &               default = 'q',         &
        &               ErrCode = iError       )

      select case(upper_to_lower(trim(polyspace)))
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
        call tem_abort()

      end select

      call close_config(conf)

    end if

  end subroutine ply_subresolution_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get the subresolution data for all elements for a given color and in the
  !! requested format.
  subroutine ply_subres_import_color( me, tree, coloring, iColor,              &
    &                                 target_degree, target_space, target_dim, &
    &                                 subresdat )
    ! --------------------------------------------------------------------- !
    type(ply_subresolution_type), intent(in) :: me
    type(treelmesh_type), intent(in) :: tree
    type(tem_color_prop_type), intent(in) :: coloring
    integer, intent(in) :: iColor
    integer, intent(in) :: target_degree
    integer, intent(in) :: target_space
    integer, intent(in) :: target_dim
    real(kind=rk), allocatable, intent(out) :: subresdat(:,:)
    ! --------------------------------------------------------------------- !
    character(len=pathLen) :: datfile
    integer :: target_Dofs
    integer :: in_dofs
    integer :: read_dofs
    integer :: pre_dofs
    integer :: pdim
    real(kind=rk), allocatable :: indat(:)
    real(kind=rk), allocatable :: predat(:)
    character(len=4) :: datext
    integer :: rl
    integer :: fUnit
    integer :: nElems
    integer(kind=long_k) :: offset
    integer :: iElem
    integer :: minOrd
    integer :: in_dim
    integer :: recs_per_elem
    ! --------------------------------------------------------------------- !

    ! Seeder always writes three-dimensional data.
    ! Assume an input dimension of 3:
    in_dim = 3

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
      !
      ! We moved this allocation in front of the loop to silence a compiler
      ! warning about a potentially uninitialized variable.
      target_dofs = target_degree+1
      select case(target_space)
        case (Q_Space)
          target_dofs = target_dofs**target_dim

        case (P_Space)
          do pdim=2,target_dim
            target_dofs = (target_dofs * (target_degree+pdim) ) / pdim
          end do

        case default
          call tem_abort( 'Wrong target_space! Select Q_Space or P_Space.' )

      end select

      ! We read the complete data per element in most cases.
      recs_per_elem = 1

      ! Figure out input polynomial representation.
      select case(me%basistype)
        case (Q_Space)
          in_dofs = (me%polydegree+1)**3

          if (target_dim < in_dim) then
            read_dofs = (me%polydegree+1)**target_dim
            recs_per_elem = (me%polydegree+1)**(in_dim-target_dim)
          else
            read_dofs = in_dofs
          end if

        case (P_Space)
          in_dofs = me%polydegree+1
          do pdim=2,in_dim
            in_dofs = (in_dofs * (me%polydegree+pdim) ) / pdim
          end do
          read_dofs = in_dofs

        case default
          call tem_abort( 'Wrong basistype! Select Q_Space or P_Space.' )

      end select

      ! Allocate arrays accordingly.
      allocate(indat(read_dofs))
      allocate(subresdat(target_dofs, nElems))

      inquire(iolength=rl) indat

      datfile = trim(tree%global%dirname) // 'subresdata_' &
        &       // trim(coloring%color_label(iColor)) // datext

      fUnit = newunit()

      open( unit   = fUnit,         &
        &   file   = datfile,       &
        &   action = 'read',        &
        &   access = 'direct',      &
        &   form   = 'unformatted', &
        &   recl   = rl,            &
        &   status = 'old'          )

      minord = min(target_degree+1, me%polydegree+1)

      subresdat = 0.0_rk

      if  ( (me%basistype == P_Space) .and. (target_dim < in_dim) ) then

        pre_dofs = me%polydegree+1
        do pdim=2,target_dim
          pre_dofs = (pre_dofs * (me%polydegree+pdim) ) / pdim
        end do
        allocate(predat(pre_dofs))

        do iElem=1,nElems

          read(fUnit, rec=offset+iElem) indat

          call ply_transfer_P_dim( indat  = indat,        &
            &                      indim  = in_dim,       &
            &                      outdat = predat,       &
            &                      outdim = target_dim,   &
            &                      degree = me%polydegree )
          call ply_transfer_dofs( indat     = predat,              &
            &                     inspace   = me%basistype,        &
            &                     indegree  = me%polydegree,       &
            &                     outdat    = subresdat(:, iElem), &
            &                     outspace  = target_space,        &
            &                     outdegree = target_degree,       &
            &                     ndims     = target_dim           )
        end do

      else

        do iElem=1,nElems

          read(fUnit, rec=(offset+iElem-1)*recs_per_elem+1) indat

          call ply_transfer_dofs( indat     = indat,               &
            &                     inspace   = me%basistype,        &
            &                     indegree  = me%polydegree,       &
            &                     outdat    = subresdat(:, iElem), &
            &                     outspace  = target_space,        &
            &                     outdegree = target_degree,       &
            &                     ndims     = target_dim           )

        end do

      end if

      close(fUnit)

    else

      ! No subresolution data at all, allocate the array with size 0.
      allocate(subresdat(0,0))

    end if

  end subroutine ply_subres_import_color
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Get the color of an element.
  !!
  !! This routine provides the get_element for the variable definition.
  !! It returns the coloring value for the elements in elempos with the given
  !! number of degrees of freedom.
  !!
  !! The header of this subroutine must be same as tem_varSys_proc_element
  subroutine ply_subres_get_elemcolor( fun, varSys, elempos, time, tree, &
    &                                  nElems, nDofs, res                )
    ! --------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of elements to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (nComponents of resulting variable) x (nDegrees of freedom) x (nElems)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! --------------------------------------------------------------------- !
    type(ply_subres_colvar_type), pointer :: p
    integer :: ipos, iElem
    integer :: icol_pos
    integer :: icol_elem
    integer :: isub_pos
    integer :: isub_elem
    integer :: colchar, colbit
    integer :: cpos
    integer :: first
    integer :: minsubdofs
    logical :: has_subres
    real(kind=rk) :: void, fill
    ! --------------------------------------------------------------------- !
    !! Dummy code to touch unused arguments
    type(tem_time_type) :: now
    if (varSys%nStateVars > 0) then
      if (tree%nElems == nElems) then
        now = time
      end if
    end if
    !! end dummy code to touch unused arguments

    call c_f_pointer(fun%method_data, p)

    cpos = p%colpos

    icol_pos = 1
    icol_elem = p%color%property%elemID(icol_pos)

    has_subres = (size(p%subresdat) > 0)

    ! Byte and bit to probe for this color
    colchar = 1 + (cpos-1) / colors_per_char
    colbit = mod((cpos-1), colors_per_char)

    ! Initialize all degrees of freedom with 0.
    res = 0.0_rk

    void = p%color%color_void(cpos)
    fill = p%color%color_fill(cpos)

    do iPos = 1, nElems

      iElem = elempos(ipos)
      ! iElem is smaller than the first colored element's ID icol_elem and thus
      ! can't be colored.
      if ( iElem < icol_elem ) then
        ! Not colored element, set first degree of freedom to the void value.
        res(nDofs*(iPos-1)+1) = void
      else
        ! The current element is larger or equal than the current colored one.
        ! Keep on increasing the colored counter, until we reach the current
        ! element or the end of the list of colored elements
        !>\todo HK: we could do a binary search on the colpos, to find the next
        !!          colored element equal or greater to the current element.
        do
          if ( (iElem <= icol_elem) &
            &  .or. (icol_pos == p%color%property%nElems) ) EXIT
          icol_pos = icol_pos + 1
          icol_elem = p%color%property%elemID(icol_pos)
        end do
        ! Now the current colored element is either the same as the current
        ! element, further ahead or the last element in the colored elements.
        if (iElem == icol_elem) then
          ! For the bit checking, we need to convert the character into an
          ! integer via ichar.
          if ( btest(ichar(p%color%colored_bit(colchar, icol_pos)), &
            &        colbit) ) then
            res(nDofs*(iPos-1)+1) = fill

            ! If there is subresolution data: Check wether this element is
            ! subresolved for this color.
            if (has_subres) then
              isub_pos = 1
              isub_elem = p%subres%subres_prop%elem(cpos)%ID(isub_pos)
              minsubdofs = min(nDofs, p%nSubdofs)

              do
                if ( (iElem <= isub_elem) &
                  &  .or. (isub_pos == p%subres%subres_prop%nElems(cpos)) ) EXIT
                isub_pos = isub_pos + 1
                isub_elem = p%subres%subres_prop%elem(cpos)%ID(isub_pos)
              end do
              if (iElem == isub_elem) then
                if ( btest(ichar(p%subres                                  &
                  &               %subres_prop                             &
                  &               %subresolved_colors(colchar, isub_pos)), &
                  &        colbit) ) then
                  ! This element has subresolution information for this color.
                  ! Set it from the subresdat accordingly.
                  first = nDofs*(iPos-1)
                  res(first+1:first+minsubdofs) &
                    & = p%subresdat(:minsubdofs, isub_pos)
                end if
                isub_pos = min(isub_pos+1, p%subres%subres_prop%nElems(cpos))
                isub_elem = p%subres%subres_prop%elem(cpos)%ID(isub_pos)
              end if

            end if

          else

            res(nDofs*(iPos-1)+1) = void

          end if

          ! We can push the position of the colored element one further, as
          ! we found the matching element already.
          icol_pos = min(icol_pos+1, p%color%property%nElems)
          icol_elem = p%color%property%elemID(icol_pos)
        else
          res(nDofs*(iPos-1)+1) = void
        end if
      end if

    end do

  end subroutine ply_subres_get_elemcolor
  ! ************************************************************************ !


end module ply_subresolution_module
